# Makefile for trajectories

FC = ifort

EXE = trajectories_fesom.x

# --- GSW Library Paths ---
GSW_INC = -I./GSW-Fortran/build/gsw -I./GSW-Fortran/build/toolbox
GSW_LIB = ./GSW-Fortran/build

OBJECTS =  mod_config.o \
	gen_modules_rotate_grid.o \
	trajectories_fesom.o \
	oce_read_mesh.o \
	mod_output.o

MODULES = mod_config.o \
	gen_modules_rotate_grid.o \
	oce_read_mesh.o \
	mod_output.o

FCFLAGS = -g -traceback -O -m64 -r8 
#FCFLAGS =   -O2 -m64 -fdefault-real-8#-fdefault-real-8 
#FCFLAGS = #-fdefault-real-8 
NETCDF_INC  =  /albedo/soft/spack-workspace/spack/spack-0.18.1/var/spack/environments/esm-tools-albedo/.spack-env/._view/6f5o6tosl5pzkwdborml37vkfyx7xntd/include 
NETCDF_LIB  =  /albedo/soft/spack-workspace/spack/spack-0.18.1/var/spack/environments/esm-tools-albedo/.spack-env/._view/6f5o6tosl5pzkwdborml37vkfyx7xntd/lib/


#module load netcdf

$(EXE): $(MODULES) $(OBJECTS)
	$(FC) -o $(EXE) $(FCFLAGS) -L$(NETCDF_LIB) -L$(GSW_LIB) $(OBJECTS) -lnetcdff -lgsw

 %.o: %.F90
	$(FC) -c $(FCFLAGS) -I$(NETCDF_INC) $(GSW_INC) mod_config.F90
	$(FC) -c $(FCFLAGS) -I$(NETCDF_INC) $(GSW_INC) gen_modules_rotate_grid.F90
	$(FC) -c $(FCFLAGS) -I$(NETCDF_INC) $(GSW_INC) oce_read_mesh.F90
	$(FC) -c $(FCFLAGS) -I$(NETCDF_INC) $(GSW_INC) mod_output.F90
	$(FC) -c $(FCFLAGS) -I$(NETCDF_INC) $(GSW_INC) trajectories_fesom.F90

clean :
	rm -f *~ *.o *.mod *.x $(EXE)


