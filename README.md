# Atlantic Water Subduction Pathways in Fram Strait

This Repository contains the code to perform Lagrangian particle tracking on the native FESOM1.4 grid an analyze the results.
Due to their size model data, model mesh as well as the resulting trajectory files and processed files can't be included in this repository.

The code is adjusted to run on AWIs HPC albedo. Running it will require adjusting paths and possibly HPC specific settings.

Most analysis notebooks are named after the plot they are generating. The exception in `Descent_Types.ipynb` which contains multiple plots and prints.
There are three main steps that have to be performed before running any of the plotting code:
1. Running the trajectory code and creating the particle netCDF file with `trajectories_fesom.f90`
2. Extract the descent windows for each particle with `Pick_subduction.ipynb`
3. Sort in to horizontal pathways with `Extract-pathways.ipynb`

These pre-steps will be described in the following. The description is mostly taylored to my future self and albedo but it might also help others to run the code.

---

## Running the Tracking Code

### Required data

- FESOM1.4 model output
- mesh data

### Required Code

- `derivative.ipynb`
- `plot_FESOM_79_tansect_TS.ipynb`
- `ocean_helper.py`
- `Create_particle_starts.ipynb`
- Tracking code with all of the modules
    - `trajectories_fesom.F90`
    - GSW library
    - …
- `auto_submit.py` (optional)

### Preparation steps

1. **Calculate basis function** and triangle area (`derivative.ipynb`) to allow for the calculation of derivatives in the tracking code
2. If its the first time running:
    1. **Plot FESOM data** (TS-diagramm and sections) using `plot_FESOM_79_tansect_TS.ipynb` to make sure the water mass definitions in `ocean_helper.py` make sense
    2. **Adjust water mass definitions** (/thresholds) as necessary in `ocean_helper.py` 
    At this point it is only important that the ones for Atlantic Water are correct because they are the only ones that are used before tracking
3. **Write Input file**
Use `Create_particle_starts.ipynb` to generate the input files for the particle tracking.
Here you can set the starting position for each of your drifters as well as the starting times.
The starting time has to be the same for all drifters in an input file in the current version of the tracking code. 
    - (Alternatively you can use `trajectories_fesom_byParticle.F90` , which is the previous version. Here you can set start dates for individually per particle, but be aware that its much, much slower as every model day is loaded for each particle, no matter how often that day might be used.)
    - Per default input files are names after the following naming convention:
    `drifter_input_*releaseyear*_*releaseday*_lat_*releaselatitude*_particles_*particle count of the unfiltered grid*.dat`
        - We use the particle count of the unfiltered grid instead of the actual count in the file so that it is possible to relate files from different release days. Because since we filter to have only AW particles starting, the actual count will vary from day to day.
4. **Compiling**
    1. As usual for fortran, we have to compile the code for running. Make sure you are in the right folder, 
    2. `module load intel-oneapi-compilers/2022.1.0`  did run once 
    3. and compile using the `Makefile` , by entering first `make clean`  to get rid of any old stuff and then `make` in your terminal.

### Running the code

Here you have two options:

- Start the tracking manually
    1. Adjusting paths and parameters in the `namelist.nml` file
    2. Submit the job with `sbatch job.batch`  from you terminal
- Auto-start multiple jobs in a row
    1. make sure you have all the required input files following the naming convention
    2. Open `auto_submit.py`  and adjust paths and parameters as needed (and save it :D)
    3. start a tmux session using `tmux new -s particles`
    4. Navigate to your code directory and run the script with `python auto_submit.py` 
    5. Detach from tmux (keeps it running in background) with: Press `Ctrl+B`, then `D`
    6. Reattach later to check progress with `tmux attach -t particles`
    
    WARNING: 
    
    1. If your submitted job would overwrite an existing output file it will pause and ask for confirmation. So yes it starts multiple jobs in a row but it might require user input in between
    2. The current version of the script always overwrites the last `namelist.nml` and `job.batch` file

Either way you can supervise the progress of your job using `squeue -u username` 
But I can also highly recommend checking out [https://github.com/Gordi42/stama](https://github.com/Gordi42/stama)
”A terminal user interface for monitoring and managing slurm jobs”. (Thanks Silvano for making my life easier :D)

### Output

So finally if everything worked out fine you will have three new files per job

- particle data as netCDF
- the log file `.out`
- and a hopefully empty error file


## Analysis

There are two main “processing” steps that are required for most plotting codes (`Fig*_*.ipynb`).

1. Running the algorithm detecting the “main” descent and PW encounter (`Pick-subduction.ipynb` )
    - This notebook calls several functions from `Subduction_picking.py` to get the Polar Water (PW) encounter and start, end and centre of the main descent for each particle track. The resulting time indices and datetimes as well as the positions and several dynamical properties at the time stamps and over the whole descent are saved in a `.csv` file following the naming convention:
    `Subduction-idxs_*releaseyear*…` And the rest of the file name of the NC-File of the trajectories
2. And sorting each particle into a group of pathways using several gates (`Extract-pathways.ipynb`)
    - This notebook checks which 'gates' the idividual drifters pass to group them into tracks (eg. northern necirculation, southern recirculation, ...)
`

Both of these produce files that will be loaded in other scripts.

Generally all notebooks use functions from `particle_loading.py` . As the name suggests it reads the netCDF data from the tracking output, with some error guards, if the time loading doesn’t work (thats mostly the case if the file comes from the old tracking code). On top of that this python file also has a number of functions that precompute variables like vorticity, Coriolis parameter, density, … from the raw data and add them to the `xarray.Dataset`. Some variables can only be calculated if full column data was saved during the tracking (eg. Richardson Number, Brunt-Väisälä frequency, …). These variables are mostly home to the `add_stratification_to_dataset(*)` function.

   


