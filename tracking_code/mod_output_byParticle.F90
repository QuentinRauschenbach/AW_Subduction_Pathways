module g_mod_output



  implicit none
  save


  contains


subroutine init_output(filename,ndrifter,nsteps,start_year)
 
  implicit none

#include "netcdf.inc" 

  character(200), INTENT(IN):: filename
  integer, INTENT(IN)       :: ndrifter,nsteps
  integer, INTENT(IN)       :: start_year
  integer                   :: status, ncid, dimid_rec, dimid_drifter
  integer                   :: dimids(2)
  integer                   :: rec_varid, blon_varid, blat_varid, bdepth_varid
  integer                   :: btemp_varid, bsalt_varid
  integer                   :: bday_varid, byear_varid
  integer                   :: drifter_varid, T_varid
  integer                   :: btemp_surface_varid, bsalt_surface_varid
  integer                   :: bvort_varid
  character(100)            :: longname
  character(1)              :: trind
  character(40)             :: time_units
  character(10)             :: calendar_type ! For the calendar attribute value
  integer, allocatable :: drifter_coord(:)
  integer :: i

  allocate(drifter_coord(ndrifter))
  drifter_coord = 0

  print*, 'initialize new output file'

  ! number of drifters: ndrifter
  ! number of time steps: nsteps


  ! create a file
  status = nf_create(filename, nf_clobber, ncid)
  if (status.ne.nf_noerr) call handle_err(status)

  ! Define the dimensions
  status = nf_def_dim(ncid, 'drifter', ndrifter, dimid_drifter)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_dim(ncid, 'time', nsteps, dimid_rec)
  if (status .ne. nf_noerr) call handle_err(status)

  ! Define coordinate variable 'drifter' (1D variable of length ndrifter)
  status = nf_def_var(ncid, 'drifter', NF_INT, 1, (/dimid_drifter/), drifter_varid)
  if (status .ne. nf_noerr) call handle_err(status)

  ! Define coordinate variable 'time' (1D variable of length nsteps)
  status = nf_def_var(ncid, 'time', NF_DOUBLE, 1, (/dimid_rec/), T_varid)
  if (status .ne. nf_noerr) call handle_err(status)

  ! Add long_name attributes to the coordinate variables.
  status = nf_put_att_text(ncid, drifter_varid, 'long_name', 12, 'Drifter ID')
  if (status .ne. nf_noerr) call handle_err(status)

  status = nf_put_att_text(ncid, T_varid, 'long_name', 4, 'Time')
  if (status .ne. nf_noerr) call handle_err(status)

  ! Construct the units string for time
  write(time_units, '(a,i4,a)') 'seconds since ', start_year, '-01-0100:00:00'

  status = nf_put_att_text(ncid, T_varid, 'units', len_trim(time_units), trim(time_units))
  if (status .ne. nf_noerr) call handle_err(status)

  calendar_type = 'noleap'
  status = nf_put_att_text(ncid, T_varid, 'calendar', len_trim(calendar_type), trim(calendar_type))
  if (status .ne. nf_noerr) call handle_err(status)


  
  ! Define the netCDF variables for 2D fields.
  ! In Fortran, the unlimited dimension must come
  ! last on the list of dimids.
  dimids(1) = dimid_rec
  dimids(2) = dimid_drifter

  status = nf_def_var(ncid, 'time_drifter', NF_DOUBLE, 2, dimids, rec_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'blon', NF_DOUBLE, 2, dimids, blon_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'blat', NF_DOUBLE, 2, dimids, blat_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'bdepth', NF_DOUBLE, 2, dimids, bdepth_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'btemp', NF_DOUBLE, 2, dimids, btemp_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'bsalt', NF_DOUBLE, 2, dimids, bsalt_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'bday', NF_INT, 2, dimids, bday_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'byear', NF_INT, 2, dimids, byear_varid)
  if (status .ne. nf_noerr) call handle_err(status)

   status = nf_def_var(ncid, 'btemp_surface', NF_DOUBLE, 2, dimids, btemp_surface_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'bsalt_surface', NF_DOUBLE, 2, dimids, bsalt_surface_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_def_var(ncid, 'bvort', NF_DOUBLE, 2, dimids, bvort_varid)
  if (status .ne. nf_noerr) call handle_err(status)


  ! Assign long_name and units attributes to variables.
  longname='model time'
  status = nf_put_att_text(ncid, rec_varid, 'long_name', len_trim(longname), trim(longname)) 
  if (status .ne. nf_noerr) call handle_err(status)
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
  longname = 'relative vorticity at drifter element'
  status = nf_put_att_text(ncid, bvort_varid, 'description', len_trim(longname), trim(longname))
  if (status .ne. nf_noerr) call handle_err(status)
  status = nf_put_att_text(ncid, bvort_varid, 'units', 6, '1/s')
  if (status .ne. nf_noerr) call handle_err(status)

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

  status=nf_close(ncid)
  if (status .ne. nf_noerr) call handle_err(status)


end subroutine init_output
!
!--------------------------------------------------------------
!
subroutine write_output(filename,ndrifter,nrec,time,blon,blat,bdepth,btemp,bsalt,bday,byear, btemp_surface, bsalt_surface, bvort)

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
