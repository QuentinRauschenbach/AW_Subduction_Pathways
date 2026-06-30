module g_mod_output



  implicit none
  save


  contains


subroutine init_output(filename, ndrifter, nsteps, ndepth_levels, depth_levels, start_year, steps_per_day, n_out, backwards, num_days, first_day, infile, save_column)
 
  implicit none

#include "netcdf.inc" 

  character(200), INTENT(IN):: filename
  integer, INTENT(IN)       :: ndrifter, nsteps, ndepth_levels
  integer, INTENT(IN)       :: start_year
  real(kind=8), dimension(ndepth_levels), INTENT(IN) :: depth_levels
  integer                   :: status, ncid, dimid_rec, dimid_drifter, dimid_depth_levels
  integer                   :: dimids(2), dimids3(3)
  integer                   :: rec_varid, blon_varid, blat_varid, bdepth_varid, depth_levels_varid
  integer                   :: btemp_varid, bsalt_varid
  integer                   :: bday_varid, byear_varid
  integer                   :: drifter_varid, T_varid
  integer                   :: btemp_surface_varid, bsalt_surface_varid
  integer                   :: btemp_column_varid, bsalt_column_varid, mld_varid
  integer                   :: du_dx_varid, dv_dx_varid, du_dy_varid, dv_dy_varid, du_dz_varid, dv_dz_varid 
  !integer                   :: bu_column_varid, bv_column_varid, bw_column_varid
  integer                   :: bu_varid, bv_varid, bw_varid
  integer                   :: bvort_varid
  character(100)            :: longname
  character(1)              :: trind
  character(40)             :: time_units
  character(10)             :: calendar_type ! For the calendar attribute value
  integer, allocatable      :: drifter_coord(:)
  integer :: i
  real(kind=8) :: fill_val = -999.0d0
  integer :: fill_vars(17)
  integer, INTENT(IN) :: steps_per_day, n_out, num_days, first_day
  logical, INTENT(IN) :: backwards, save_column
  character(200), INTENT(IN) :: infile
  integer :: actual_len


  allocate(drifter_coord(ndrifter))
  drifter_coord = 0

  print*, 'initialize new output file'

  ! number of drifters: ndrifter
  ! number of time steps: nsteps

  print*, '  number of drifters: ', ndrifter
  print*, '  number of time steps: ', nsteps
  print*, '  number of depth levels:  ', ndepth_levels


  ! create a file
  !status = nf_create(filename, nf_clobber, ncid)
  status = nf_create(filename, ior(nf_clobber, nf_netCDF4), ncid)
  if (status.ne.nf_noerr) call handle_err(status)

  ! Define the dimensions
  status = nf_def_dim(ncid, 'drifter', ndrifter, dimid_drifter)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_dim(ncid, 'time', nsteps, dimid_rec)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_dim(ncid, 'depth_levels', ndepth_levels, dimid_depth_levels)
  if (status .ne. nf_noerr) call handle_err(status)

  ! Define coordinate variable 'drifter' (1D variable of length ndrifter)
  status = nf_def_var(ncid, 'drifter', NF_INT, 1, (/dimid_drifter/), drifter_varid)
  if (status .ne. nf_noerr) call handle_err(status)

  ! Define coordinate variable 'time' (1D variable of length nsteps)
  status = nf_def_var(ncid, 'time', NF_DOUBLE, 1, (/dimid_rec/), T_varid)
  if (status .ne. nf_noerr) call handle_err(status)

  ! Define coordinate variable 'depth_levels' (1D variable of length ndepth_levels)
  status = nf_def_var(ncid, 'depth_levels', NF_DOUBLE, 1, (/dimid_depth_levels/), depth_levels_varid)
  if (status .ne. nf_noerr) call handle_err(status)

  ! Add long_name attributes to the coordinate variables.
  status = nf_put_att_text(ncid, drifter_varid, 'long_name', 12, 'Drifter ID')
  if (status .ne. nf_noerr) call handle_err(status)

  status = nf_put_att_text(ncid, T_varid, 'long_name', 4, 'Time')
  if (status .ne. nf_noerr) call handle_err(status)

  status = nf_put_att_text(ncid, depth_levels_varid, 'long_name', 13, 'Depth Levels')
  if (status .ne. nf_noerr) call handle_err(status)

  ! Construct the units string for time
  write(time_units, '(a,i4,a)') 'seconds since ', start_year, '-01-0100:00:00'

  status = nf_put_att_text(ncid, T_varid, 'units', len_trim(time_units), trim(time_units))
  if (status .ne. nf_noerr) call handle_err(status)

  calendar_type = 'noleap'
  status = nf_put_att_text(ncid, T_varid, 'calendar', len_trim(calendar_type), trim(calendar_type))
  if (status .ne. nf_noerr) call handle_err(status)

  status = nf_def_var(ncid, 'day', NF_INT, 1, (/dimid_rec/), bday_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'year', NF_INT, 1, (/dimid_rec/), byear_varid)
  if (status .ne. nf_noerr) call handle_err(status)


  
  ! Define the netCDF variables for 2D fields.
  ! In Fortran, the unlimited dimension must come
  ! last on the list of dimids.
  dimids(1) = dimid_rec
  dimids(2) = dimid_drifter

  ! For 3D variables (time, drifter, depth level)
  dimids3(1) = dimid_rec
  dimids3(2) = dimid_drifter
  dimids3(3) = dimid_depth_levels

  !status = nf_def_var(ncid, 'time_drifter', NF_DOUBLE, 2, dimids, rec_varid)
  !if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'lon', NF_DOUBLE, 2, dimids, blon_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'lat', NF_DOUBLE, 2, dimids, blat_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'depth', NF_DOUBLE, 2, dimids, bdepth_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'temp', NF_DOUBLE, 2, dimids, btemp_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'salt', NF_DOUBLE, 2, dimids, bsalt_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'u', NF_DOUBLE, 2, dimids, bu_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'v', NF_DOUBLE, 2, dimids, bv_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'w', NF_DOUBLE, 2, dimids, bw_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'du_dx', NF_DOUBLE, 2, dimids, du_dx_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'dv_dx', NF_DOUBLE, 2, dimids, dv_dx_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'du_dy', NF_DOUBLE, 2, dimids, du_dy_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'dv_dy', NF_DOUBLE, 2, dimids, dv_dy_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'du_dz', NF_DOUBLE, 2, dimids, du_dz_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'dv_dz', NF_DOUBLE, 2, dimids, dv_dz_varid)
  if (status .ne. nf_noerr) call handle_err(status)

  status = nf_def_var(ncid, 'temp_surface', NF_DOUBLE, 2, dimids, btemp_surface_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'salt_surface', NF_DOUBLE, 2, dimids, bsalt_surface_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'mld', NF_DOUBLE, 2, dimids, mld_varid)
  if (status .ne. nf_noerr) call handle_err(status)

if (save_column) then
  status = nf_def_var(ncid, 'temp_column', NF_DOUBLE, 3, dimids3, btemp_column_varid)
  if (status .ne. nf_noerr) call handle_err(status)

  status = nf_def_var(ncid, 'salt_column', NF_DOUBLE, 3, dimids3, bsalt_column_varid)
  if (status .ne. nf_noerr) call handle_err(status)
end if

!status = nf_def_var(ncid, 'bu_column', NF_DOUBLE, 3, dimids3, bu_column_varid)
!if (status .ne. nf_noerr) call handle_err(status)

!status = nf_def_var(ncid, 'bv_column', NF_DOUBLE, 3, dimids3, bv_column_varid)
!if (status .ne. nf_noerr) call handle_err(status)

!status = nf_def_var(ncid, 'bw_column', NF_DOUBLE, 3, dimids3, bw_column_varid)
!if (status .ne. nf_noerr) call handle_err(status)


  ! Assign long_name and units attributes to variables.
  !longname='model time'
  !status = nf_put_att_text(ncid, rec_varid, 'long_name', len_trim(longname), trim(longname)) 
  !if (status .ne. nf_noerr) call handle_err(status)
  !status = nf_put_att_text(ncid, rec_varid, 'units', 1, 'd')
  !if (status .ne. nf_noerr) call handle_err(status)
  ! Construct the units string for time

  !status = nf_put_att_text(ncid, rec_varid, 'units', len_trim(time_units), trim(time_units))
  !if (status .ne. nf_noerr) call handle_err(status)

  !calendar_type = 'noleap'
  !status = nf_put_att_text(ncid, rec_varid, 'calendar', len_trim(calendar_type), trim(calendar_type))
  !if (status .ne. nf_noerr) call handle_err(status)


  status = nf_put_att_text(ncid, btemp_surface_varid, 'long_name', 29, 'Surface temperature above drifter')
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, bsalt_surface_varid, 'long_name', 27, 'Surface salinity above drifter')
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, bvort_varid, 'long_name', 16, 'relative vorticity at drifter element')
  if (status .ne. nf_noerr) call handle_err(status)
  


  longname='longitude'
  status = nf_put_att_text(ncid, blon_varid, 'description', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, blon_varid, 'units', 3, 'deg')
  if (status .ne. nf_noerr) call handle_err(status)
  longname='latitude'
  status = nf_put_att_text(ncid, blat_varid, 'description', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, blat_varid, 'units', 3, 'deg')
  if (status .ne. nf_noerr) call handle_err(status)
  longname='depth'
  status = nf_put_att_text(ncid, bdepth_varid, 'description', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, bdepth_varid, 'units', 1, 'm')
  if (status .ne. nf_noerr) call handle_err(status)
  longname='potential temperature'
  status = nf_put_att_text(ncid, btemp_varid, 'description', len_trim(longname), trim(longname))
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, btemp_varid, 'units', 4, 'degC')
  if (status .ne. nf_noerr) call handle_err(status)
  longname='salinity'
  status = nf_put_att_text(ncid, bsalt_varid, 'description', len_trim(longname), longname) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, bsalt_varid, 'units', 3, 'psu')
  if (status .ne. nf_noerr) call handle_err(status)
  longname='day'
  status = nf_put_att_text(ncid, bday_varid, 'description', len_trim(longname), longname) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, bday_varid, 'units', 1, 'd')
  if (status .ne. nf_noerr) call handle_err(status)
  longname='year'
  status = nf_put_att_text(ncid, byear_varid, 'description', len_trim(longname), longname) 
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, byear_varid, 'units', 2, 'yr')
  if (status .ne. nf_noerr) call handle_err(status)
  !longname = 'relative vorticity at drifter element'
  !status = nf_put_att_text(ncid, bvort_varid, 'description', len_trim(longname), trim(longname))
  !if (status .ne. nf_noerr) call handle_err(status)
  !status = nf_put_att_text(ncid, bvort_varid, 'units', 6, '1/s')
  !if (status .ne. nf_noerr) call handle_err(status)

  ! du_dx
longname = 'weighted mean of du/dx at drifter particle position'
  status = nf_put_att_text(ncid, du_dx_varid, 'description', len_trim(longname), trim(longname))
if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, du_dx_varid, 'units', 3, '1/s')
if (status .ne. nf_noerr) call handle_err(status)

! du_dy
longname = 'weighted mean of du/dy at drifter particle position'
  status = nf_put_att_text(ncid, du_dy_varid, 'description', len_trim(longname), trim(longname))
if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, du_dy_varid, 'units', 3, '1/s')
if (status .ne. nf_noerr) call handle_err(status)

! dv_dx
longname = 'weighted mean of dv/dx at drifter particle position'
  status = nf_put_att_text(ncid, dv_dx_varid, 'description', len_trim(longname), trim(longname))
if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, dv_dx_varid, 'units', 3, '1/s')
if (status .ne. nf_noerr) call handle_err(status)

! dv_dy
longname = 'weighted mean of dv/dy at drifter particle position'
  status = nf_put_att_text(ncid, dv_dy_varid, 'description', len_trim(longname), trim(longname))
if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, dv_dy_varid, 'units', 3, '1/s')
if (status .ne. nf_noerr) call handle_err(status)

! du_dz
longname = 'weighted mean of du/dz at drifter particle position'
  status = nf_put_att_text(ncid, du_dz_varid, 'description', len_trim(longname), trim(longname))
if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, du_dz_varid, 'units', 3, '1/s')
if (status .ne. nf_noerr) call handle_err(status)

! dv_dz
longname = 'weighted mean of dv/dz at drifter particle position'
  status = nf_put_att_text(ncid, dv_dz_varid, 'description', len_trim(longname), trim(longname))
if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, dv_dz_varid, 'units', 3, '1/s')
if (status .ne. nf_noerr) call handle_err(status)

! mld
longname = 'mixed layer depth (0.03)'
  status = nf_put_att_text(ncid, mld_varid, 'description', len_trim(longname), trim(longname))
if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, mld_varid, 'units', 3, 'm')
if (status .ne. nf_noerr) call handle_err(status)

  ! Add _FillValue attribute to all variables (except time and drifter coordinate)
  fill_vars = (/ blon_varid, blat_varid, bdepth_varid, btemp_varid, bsalt_varid, &
                 btemp_surface_varid, bsalt_surface_varid, &
                 du_dx_varid, dv_dx_varid, du_dy_varid, dv_dy_varid, du_dz_varid, dv_dz_varid, &
                 bu_varid, bv_varid, bw_varid, mld_varid /)

  do i = 1, size(fill_vars)
     status = nf_put_att_double(ncid, fill_vars(i), '_FillValue', NF_DOUBLE, 1, fill_val)
     if (status .ne. nf_noerr) call handle_err(status)
  end do

  if (save_column) then
   status = nf_put_att_double(ncid, btemp_column_varid, '_FillValue', NF_DOUBLE, 1, fill_val)
   if (status .ne. nf_noerr) call handle_err(status)
   status = nf_put_att_double(ncid, bsalt_column_varid, '_FillValue', NF_DOUBLE, 1, fill_val)
   if (status .ne. nf_noerr) call handle_err(status)
  end if

  ! --- Global Attributes ---
  
  ! Integers
  status = nf_put_att_int(ncid, NF_GLOBAL, 'steps_per_day', NF_INT, 1, steps_per_day)
  status = nf_put_att_int(ncid, NF_GLOBAL, 'n_out', NF_INT, 1, n_out)
  status = nf_put_att_int(ncid, NF_GLOBAL, 'total_experiment_days', NF_INT, 1, num_days)
  ! Release Information
  status = nf_put_att_int(ncid, NF_GLOBAL, 'release_start_year', NF_INT, 1, start_year)
  status = nf_put_att_int(ncid, NF_GLOBAL, 'release_start_doy', NF_INT, 1, first_day)

  ! Add a description for clarity
  status = nf_put_att_text(ncid, NF_GLOBAL, 'output_description', 155, &
      'Particle tracking run starting on specified release date. ' // &
      'Numerical integration frequency: steps_per_day. ' // &
      'Storage frequency: every n_out integration steps.')


  ! Boolean (converted to string for readability)
  if (backwards) then
     status = nf_put_att_text(ncid, NF_GLOBAL, 'direction', 17, 'backward tracking')
  else
     status = nf_put_att_text(ncid, NF_GLOBAL, 'direction', 16, 'forward tracking')
  end if

  ! Add input file name as a global attribute for traceability
  if (actual_len <= 0) actual_len = len_trim(infile)

status = nf_put_att_text(ncid, NF_GLOBAL, 'start_position_file', actual_len, infile(1:actual_len))

  if (save_column) then
    status = nf_put_att_text(ncid, NF_GLOBAL, 'save_column', 4, 'true')
  else
    status = nf_put_att_text(ncid, NF_GLOBAL, 'save_column', 5, 'false')
  end if

  status = nf_enddef(ncid)
  if (status .ne. nf_noerr) call handle_err(status)

  ! Write drifter coordinate variable
  !allocate(drifter_coord(ndrifter))
  do i=1, ndrifter
     drifter_coord(i) = i
  end do
  status = nf_put_var_int(ncid, drifter_varid, drifter_coord)
  if (status .ne. nf_noerr) call handle_err(status)
  deallocate(drifter_coord)

  

  ! Write depth levels coordinate variable
  status = nf_put_vara_double(ncid, depth_levels_varid, (/1/), (/ndepth_levels/), depth_levels)
  if (status .ne. nf_noerr) call handle_err(status)

  ! close file
  status=nf_close(ncid)
  if (status .ne. nf_noerr) call handle_err(status)

  
end subroutine init_output
!
!--------------------------------------------------------------

subroutine write_output_byDay(filename, ndrifter, nrec, ndepth_levels, start_time_idx, time, blon, blat, bdepth, btemp, bsalt, bday, byear, btemp_surface, bsalt_surface, du_dx, dv_dx, du_dy, dv_dy, du_dz, dv_dz, bu, bv, bw, mld, save_column, btemp_column, bsalt_column)!, bu_column, bv_column, bw_column

  implicit none

#include "netcdf.inc" 

  integer                   :: status, ncid
  integer                   :: T_varid, blon_varid, blat_varid, bdepth_varid
  integer                   :: btemp_varid, bsalt_varid
  integer                   :: bday_varid, byear_varid
  integer                   :: btemp_surface_varid, bsalt_surface_varid
  !integer                   :: bvort_varid
  integer                   :: btemp_column_varid, bsalt_column_varid, mld_varid
  !integer                   :: bu_column_varid, bv_column_varid, bw_column_varid
  integer                   :: bu_varid, bv_varid, bw_varid
  integer                   :: du_dx_varid, dv_dx_varid, du_dy_varid, dv_dy_varid, du_dz_varid, dv_dz_varid

  integer                   :: start(2), count(2), count_3D(3), start_3D(3)
  integer                   :: start_T(1), count_T(1)
  real(kind=8), allocatable :: aux2(:,:)
  integer, allocatable      :: aux3(:)

  character(200), INTENT(IN)          :: filename
  integer, INTENT(IN)                 :: ndrifter, nrec, ndepth_levels, start_time_idx
  logical, intent(IN) :: save_column

  ! 1D arrays for time/day info (same for all drifters)
  real(kind=8),dimension(nrec), INTENT(IN)    :: time
  integer,dimension(nrec), INTENT(IN)         :: bday, byear

  ! 2D arrays for drifter specific info (time, drifter)
  real(kind=8),dimension(nrec, ndrifter), INTENT(IN) :: blon, blat, bdepth, btemp, bsalt
  real(kind=8), dimension(nrec, ndrifter), intent(IN) :: btemp_surface, bsalt_surface, bu, bv, bw, du_dx, dv_dx, du_dy, dv_dy, du_dz, dv_dz, mld ! bvort

  ! 3D arrays for drifter specific column info (time, drifter, depth level)
  real(kind=8), dimension(nrec, ndrifter, ndepth_levels), intent(IN), optional :: btemp_column, bsalt_column !, bu_column, bv_column, bw_column


  integer :: len_time, len_drifter, len_depth
  integer :: dimid_time, dimid_drifter_dbg, dimid_depth_dbg

  ! Allocate 2D auxiliary arrays 
  allocate(aux2(nrec, ndrifter)) 
  allocate(aux3(nrec))
 
  ! open file
  status = nf_open(filename, nf_write, ncid)
  if (status .ne. nf_noerr) call handle_err(status)

  print*, 'writing to output to file: ', trim(filename)

  status = nf_inq_dimid(ncid, 'time', dimid_time)
if (status .ne. nf_noerr) call handle_err(status)

status = nf_inq_dimlen(ncid, dimid_time, len_time)
if (status .ne. nf_noerr) call handle_err(status)

status = nf_inq_dimid(ncid, 'drifter', dimid_drifter_dbg)
if (status .ne. nf_noerr) call handle_err(status)

status = nf_inq_dimlen(ncid, dimid_drifter_dbg, len_drifter)
if (status .ne. nf_noerr) call handle_err(status)

status = nf_inq_dimid(ncid, 'depth_levels', dimid_depth_dbg)
if (status .ne. nf_noerr) call handle_err(status)

status = nf_inq_dimlen(ncid, dimid_depth_dbg, len_depth)
if (status .ne. nf_noerr) call handle_err(status)

  !print *, '  file dims:'
  !print *, '    time          = ', len_time
  !print *, '    drifter       = ', len_drifter
  !print *, '    depth_levels  = ', len_depth


  ! inquire variable id
  status=nf_inq_varid(ncid, 'time', T_varid) 
  if (status .ne. nf_noerr) call handle_err(status)
  status=nf_inq_varid(ncid, 'lon', blon_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status=nf_inq_varid(ncid, 'lat', blat_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status=nf_inq_varid(ncid, 'depth', bdepth_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status=nf_inq_varid(ncid, 'temp', btemp_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status=nf_inq_varid(ncid, 'salt', bsalt_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status=nf_inq_varid(ncid, 'day', bday_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status=nf_inq_varid(ncid, 'year', byear_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_inq_varid(ncid, 'temp_surface', btemp_surface_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_inq_varid(ncid, 'salt_surface', bsalt_surface_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  !status = nf_inq_varid(ncid, 'vort', bvort_varid)
  !if (status .ne. nf_noerr) call handle_err(status)

  ! Velocities
  status = nf_inq_varid(ncid, 'u', bu_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_inq_varid(ncid, 'v', bv_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_inq_varid(ncid, 'w', bw_varid)
  if (status .ne. nf_noerr) call handle_err(status)

  ! Derivatives
  status = nf_inq_varid(ncid, 'du_dx', du_dx_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_inq_varid(ncid, 'dv_dx', dv_dx_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_inq_varid(ncid, 'du_dy', du_dy_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_inq_varid(ncid, 'dv_dy', dv_dy_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_inq_varid(ncid, 'du_dz', du_dz_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_inq_varid(ncid, 'dv_dz', dv_dz_varid)
  if (status .ne. nf_noerr) call handle_err(status)

  status = nf_inq_varid(ncid, 'mld', mld_varid)
  if (status .ne. nf_noerr) call handle_err(status)

  if (save_column .and. present(btemp_column) .and. present(bsalt_column)) then
    status = nf_inq_varid(ncid, 'temp_column', btemp_column_varid)
    if (status .ne. nf_noerr) call handle_err(status)
    status = nf_inq_varid(ncid, 'salt_column', bsalt_column_varid)
    if (status .ne. nf_noerr) call handle_err(status)
  end if


  !status = nf_inq_varid(ncid, 'bu_column', bu_column_varid)
  !if (status .ne. nf_noerr) call handle_err(status)
  !status = nf_inq_varid(ncid, 'bv_column', bv_column_varid)
  !if (status .ne. nf_noerr) call handle_err(status)
  !status = nf_inq_varid(ncid, 'bw_column', bw_column_varid)
  !if (status .ne. nf_noerr) call handle_err(status)


  ! write variables

  start_T = (/start_time_idx/)       
  count_T = (/nrec/)         
  status = nf_put_vara_double(ncid, T_varid, start_T, count_T, time)
  if (status .ne. nf_noerr) call handle_err(status)

  ! Start writing at the first time step of the current day (start_time_idx) 
  ! and the first drifter (1).
  start = (/start_time_idx, 1/)
  ! Write all steps for all drifters
  count = (/nrec, ndrifter/)
  
  aux2 = blon
  status = nf_put_vara_double(ncid, blon_varid, start, count, aux2) 
  if (status .ne. nf_noerr) call handle_err(status)
             
  aux2 = blat
  status = nf_put_vara_double(ncid, blat_varid, start, count, aux2) 
  if (status .ne. nf_noerr) call handle_err(status)

  aux2 = bdepth
  status = nf_put_vara_double(ncid, bdepth_varid, start, count, aux2) 
  if (status .ne. nf_noerr) call handle_err(status)

  aux2 = btemp
  status = nf_put_vara_double(ncid, btemp_varid, start, count, aux2) 
  if (status .ne. nf_noerr) call handle_err(status)

  aux2 = bsalt
  status = nf_put_vara_double(ncid, bsalt_varid, start, count, aux2) 
  if (status .ne. nf_noerr) call handle_err(status)

  aux2 = bu
  status = nf_put_vara_double(ncid, bu_varid, start, count, aux2) 
  if (status .ne. nf_noerr) call handle_err(status) 
  aux2 = bv
  status = nf_put_vara_double(ncid, bv_varid, start, count, aux2) 
  if (status .ne. nf_noerr) call handle_err(status)
  aux2 = bw
  status = nf_put_vara_double(ncid, bw_varid, start, count, aux2) 
  if (status .ne. nf_noerr) call handle_err(status)

  status = nf_put_vara_double(ncid, btemp_surface_varid, start, count, btemp_surface)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_vara_double(ncid, bsalt_surface_varid, start, count, bsalt_surface)
  if (status .ne. nf_noerr) call handle_err(status)

  !status = nf_put_vara_double(ncid, bvort_varid, start, count, bvort)
  !if (status .ne. nf_noerr) call handle_err(status)

  status = nf_put_vara_double(ncid, du_dx_varid, start, count, du_dx)
  if (status .ne. nf_noerr) call handle_err(status)

  status = nf_put_vara_double(ncid, dv_dx_varid, start, count, dv_dx)
  if (status .ne. nf_noerr) call handle_err(status)

  status = nf_put_vara_double(ncid, du_dy_varid, start, count, du_dy)
  if (status .ne. nf_noerr) call handle_err(status)

  status = nf_put_vara_double(ncid, dv_dy_varid, start, count, dv_dy)
  if (status .ne. nf_noerr) call handle_err(status)

  status = nf_put_vara_double(ncid, du_dz_varid, start, count, du_dz)
  if (status .ne. nf_noerr) call handle_err(status)

  status = nf_put_vara_double(ncid, dv_dz_varid, start, count, dv_dz)
  if (status .ne. nf_noerr) call handle_err(status)

  status = nf_put_vara_double(ncid, mld_varid, start, count, mld)
  if (status .ne. nf_noerr) call handle_err(status)


  if (save_column .and. present(btemp_column) .and. present(bsalt_column)) then
  ! Write 3D column variables
    start_3D = (/start_time_idx, 1, 1/)
    count_3D = (/nrec, ndrifter, ndepth_levels/)

    !print *, '  start_3D = ', start_3D
    !print *, '  count_3D = ', count_3D


    status = nf_put_vara_double(ncid, btemp_column_varid, start_3D, count_3D, btemp_column)
    if (status .ne. nf_noerr) call handle_err(status)
    status = nf_put_vara_double(ncid, bsalt_column_varid, start_3D, count_3D, bsalt_column)
    if (status .ne. nf_noerr) call handle_err(status)
    !status = nf_put_vara_double(ncid, bu_column_varid, start_3D, count_3D, bu_column)
    !if (status .ne. nf_noerr) call handle_err(status)
    !status = nf_put_vara_double(ncid, bv_column_varid, start_3D, count_3D, bv_column)
    !if (status .ne. nf_noerr) call handle_err(status)
    !status = nf_put_vara_double(ncid, bw_column_varid, start_3D, count_3D, bw_column)
  end if


  aux3(:) = bday
  status=nf_put_vara_int(ncid, bday_varid, start, count, aux3) 
  if (status .ne. nf_noerr) call handle_err(status)

  aux3(:) = byear
  status=nf_put_vara_int(ncid, byear_varid, start, count, aux3) 
  if (status .ne. nf_noerr) call handle_err(status)

  ! close file
  status=nf_close(ncid)
  if (status .ne. nf_noerr) call handle_err(status)


  deallocate(aux2)
  deallocate(aux3)

end subroutine write_output_byDay

!
subroutine write_output(filename, ndrifter, nrec, time, blon, blat, bdepth, btemp, bsalt, bday, byear, btemp_surface, bsalt_surface, bvort)

  implicit none

#include "netcdf.inc" 

  integer                   :: status, ncid
  integer                   :: T_varid, rec_varid, blon_varid, blat_varid, bdepth_varid
  integer                   :: btemp_varid, bsalt_varid
  integer                   :: bday_varid, byear_varid
  integer                   :: btemp_surface_varid, bsalt_surface_varid
  integer                   :: bvort_varid

  integer                   :: start(2), count(2)
  integer                   :: start_T(1), count_T(1)
  real(kind=8), allocatable :: aux2(:)
  integer, allocatable      :: aux3(:)

  character(200), INTENT(IN)          :: filename
  integer, INTENT(IN)                 :: ndrifter,nrec
  real,dimension(nrec), INTENT(IN)    :: time, blat, blon, bdepth, btemp, bsalt
  real, dimension(nrec), intent(in) :: btemp_surface, bsalt_surface
  integer,dimension(nrec), INTENT(IN) :: bday, byear
  real, dimension(nrec), intent(in) :: bvort

 

  allocate(aux2(nrec)) 
  allocate(aux3(nrec)) 
 
  ! open file
  status = nf_open(filename, nf_write, ncid)
  if (status .ne. nf_noerr) call handle_err(status)

  ! inquire variable id
  status=nf_inq_varid(ncid, 'time', T_varid) 
  if (status .ne. nf_noerr) call handle_err(status)
  status=nf_inq_varid(ncid, 'time_drifter', rec_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status=nf_inq_varid(ncid, 'blon', blon_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status=nf_inq_varid(ncid, 'blat', blat_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status=nf_inq_varid(ncid, 'bdepth', bdepth_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status=nf_inq_varid(ncid, 'btemp', btemp_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status=nf_inq_varid(ncid, 'bsalt', bsalt_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status=nf_inq_varid(ncid, 'bday', bday_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status=nf_inq_varid(ncid, 'byear', byear_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_inq_varid(ncid, 'btemp_surface', btemp_surface_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_inq_varid(ncid, 'bsalt_surface', bsalt_surface_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_inq_varid(ncid, 'bvort', bvort_varid)
  if (status .ne. nf_noerr) call handle_err(status)


  ! write variables
  aux2(1:nrec)=time

  start_T = (/1/)       
  count_T = (/nrec/)         
  status = nf_put_vara_double(ncid, T_varid, start_T, count_T, time)
  if (status .ne. nf_noerr) call handle_err(status)

  start=(/1,ndrifter/)
  count=(/nrec, 1/)
  
  status=nf_put_vara_double(ncid, rec_varid, start, count, aux2) 
  if (status .ne. nf_noerr) call handle_err(status)
  
  aux2(1:nrec)=blon
  status=nf_put_vara_double(ncid, blon_varid, start, count, aux2) 
  if (status .ne. nf_noerr) call handle_err(status)
             
  aux2(1:nrec)=blat
  status=nf_put_vara_double(ncid, blat_varid, start, count, aux2) 
  if (status .ne. nf_noerr) call handle_err(status)

  aux2(1:nrec)=bdepth
  status=nf_put_vara_double(ncid, bdepth_varid, start, count, aux2) 
  if (status .ne. nf_noerr) call handle_err(status)

  aux2(1:nrec)=btemp
  status=nf_put_vara_double(ncid, btemp_varid, start, count, aux2) 
  if (status .ne. nf_noerr) call handle_err(status)

  aux2(1:nrec)=bsalt
  status=nf_put_vara_double(ncid, bsalt_varid, start, count, aux2) 
  if (status .ne. nf_noerr) call handle_err(status)

  status = nf_put_vara_double(ncid, btemp_surface_varid, start, count, btemp_surface)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_vara_double(ncid, bsalt_surface_varid, start, count, bsalt_surface)
  if (status .ne. nf_noerr) call handle_err(status)

  status = nf_put_vara_double(ncid, bvort_varid, start, count, bvort)
  if (status .ne. nf_noerr) call handle_err(status)


  aux3(1:nrec)=bday
  status=nf_put_vara_int(ncid, bday_varid, start, count, aux3) 
  if (status .ne. nf_noerr) call handle_err(status)

  aux3(1:nrec)=byear
  status=nf_put_vara_int(ncid, byear_varid, start, count, aux3) 
  if (status .ne. nf_noerr) call handle_err(status)

  ! close file
  status=nf_close(ncid)
  if (status .ne. nf_noerr) call handle_err(status)


  deallocate(aux2)
  deallocate(aux3)

end subroutine write_output


end module g_mod_output
