import os
import subprocess
from datetime import datetime
import time


# ── Simulation range ──────────────────────────────────────
years         = range(2006, 2008)  # Runs 2006, 2007, 2008, 2009, 2010
seasons       = [274] #[1, 91, 182, 274] # start day of year for each season
tracking_days = 730                # number of days to track particles
 
# ── Particle release parameters ───────────────────────────
release_lat       = 76    # Latitude where particles are released
particles_per_run = 123 #subset batch # 2640  # Number of particles released in each run (before filtering for Atlantic Water only)
 
# ── Integration parameters ────────────────────────────────
steps_per_day = 240    # Integration steps per day
n_out         = 20     # Save output every xth time step
twoD_tracking = False  # Whether to track only in 2D (setting vertical velocity w=0)
backwards     = False  # Whether to track backwards in time
put_back      = False  # What to do at the boundary (put particle back to nearest element center)
save_column   = True   # Whether to save temperature and salinity of the whole water column 
                       # Saving the full column is not really feasable for either very high particle counts nor long tracking periods (file loading during each integration day takes too long)
 
# ── Paths ─────────────────────────────────────────────────
mesh_path       = "/albedo/work/projects/oce_rio/cwekerle/mesh/Arc08_sub/"
fesom_path      = "/albedo/work/projects/oce_rio/cwekerle/result/Arc40/from_1988_erai/"
fesom_runid     = "Arc40"
input_path      = "/albedo/work/user/quraus001/FESOM_particles/input/"
output_path     = "/albedo/work/user/quraus001/FESOM_particles/result/"
bafu_file       = "/albedo/work/user/quraus001/FESOM_particles/preparation/bafux_bafuy_2d.nc"
executable_path = "/albedo/work/user/quraus001/FESOM_particles/code/trajectories_fesom.x"
 
# ─────────────────────────────────────────────────────────
# Helper: convert Python bool to Fortran .true./.false.
# ─────────────────────────────────────────────────────────

 
def fortran_bool(value):
    return ".true." if value else ".false."
 

# ─────────────────────────────────────────────────────────
# Main loop
# ─────────────────────────────────────────────────────────
today_str = datetime.now().strftime("%Y-%m-%d")

for year in years:
    for doy in seasons:
 
        # Build filenames
        infile  = f"drifter_input_{year}-{doy:03d}_lat-{release_lat}_particles-{particles_per_run}.dat"
        outfile = f"drifter_output_{year}-{doy:03d}_lat-{release_lat}_particles-{particles_per_run}_tracking-{tracking_days}d_{today_str}.nc"
        #outfile = f"testing-auto-submit_{year}-{doy:03d}.nc"
 
        print(f"Year {year}, DOY {doy:03d}")
        # Check if the input file exists
        if not os.path.isfile(os.path.join(input_path, infile)):
            print(f"  SKIPPING: Failed to find {infile}")
            continue
 
        # Check if the output file already exists
        if os.path.isfile(os.path.join(output_path, outfile)):
            print(f"  Output file already exists. Do you want to OVERWRITE? (yes to continue, anything else to skip)")
            overwriting = input()
            if overwriting != "yes":
                print(f"  Skipping.")
                continue
            else:
                # Delete old file to avoid permission denied errors in the fortran code...
                try:
                    # 5-second pause *BEFORE* deleting the file
                    print("\n  ⚠️ WARNING: File will be permanently deleted and overwritten in 5 seconds...")
                    print("  Press Ctrl+C now to abort the script!")
                    
                    for i in range(5, 0, -1):
                        print(f"  Deleting in {i}...", end="\r")
                        time.sleep(1)
                    print("  Deleting old file...                         ") # Clears the countdown line
                    
                    # Delete old file to avoid permission denied errors in the fortran code...
                    full_output_path = os.path.join(output_path, outfile)
                    os.remove(full_output_path)
                    print(f"  Old file deleted successfully. Proceeding with execution.")
                    
                except KeyboardInterrupt:
                    print("\n\n  Execution aborted by user. File was NOT deleted. Exiting safely.")
                    break # or 'sys.exit()' to completely kill the script
                except OSError as e:
                    print(f"  Error deleting file: {e}. Skipping this file.")
                    continue
                



        print(f"  in:  {infile}")
        print(f"  out: {outfile}")
 
        # Write namelist.nml
        namelist = f"""&Paths
MeshPath     = "{mesh_path}"
ResultPath   = "{fesom_path}"
InputPath    = "{input_path}"
OutputPath   = "{output_path}"
BafuFile     = "{bafu_file}"
outfile      = "{outfile}"
infile       = "{infile}"
runid        = "{fesom_runid}"
/
&Parameters
syear         = {year}
num_days      = {tracking_days}
steps_per_day = {steps_per_day}
n_out         = {n_out}
twoD_tracking = {fortran_bool(twoD_tracking)}
backwards     = {fortran_bool(backwards)}
put_back      = {fortran_bool(put_back)}
save_column   = {fortran_bool(save_column)}
/
"""
        with open("namelist.nml", "w") as f:
            f.write(namelist)
 
        #print(f"  namelist.nml written")
        print()

# Write job.batch
        job_name = f"traj_{year}_{doy:03d}"
        batch = f"""#!/bin/bash
#SBATCH --account=po_physoze.oce_rio
#SBATCH --job-name={job_name}
#SBATCH --partition=smp
#SBATCH --time=12:00:00
#SBATCH --qos=12h
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=32000
#SBATCH -o %x_%j.out
#SBATCH -e error_%x_%j.out
#SBATCH --hint=nomultithread
 
export OMP_NUM_THREADS=1
ulimit -s unlimited
 
module load netcdf-fortran/4.5.4-intel-oneapi-mpi2021.6.0-intel2021.6.0
module load intel-oneapi-compilers/2022.1.0
 
JOBID=`echo $SLURM_JOB_ID | cut -d"." -f1`
srun {executable_path}
"""
        with open("job.batch", "w") as f:
            f.write(batch)
 
        # Submit and wait for it to finish before starting the next one
        print(f"  Submitting {job_name} ...")
        result = subprocess.run(["sbatch", "--wait", "job.batch"])
 
        if result.returncode == 0:
            print(f"  Done!\n")
        else:
            print(f"  Job failed! Stopping here so you can check what went wrong.\n")
            break  # stop the inner loop
    else:
        continue  # inner loop completed normally, keep going
    break          # inner loop was broken (job failed), stop outer loop too
 
print("All done.")