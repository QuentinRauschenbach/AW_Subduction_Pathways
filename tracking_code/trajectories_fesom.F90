program trajectories
! written by Katja Rollenhagen / Bettina Fach / Ralph Timmermann / Claudia Wekerle (AWI)

! edited by Quentin Rauschenbach

use g_config
use g_rotate_grid
use o_param
use o_mesh
use o_elements
use g_read_mesh
use g_mod_output
use o_DATA_TYPES
use gsw_mod_toolbox

implicit none

! ---------------------------------------------------------------------
! Configuration / namelist inputs
! ---------------------------------------------------------------------
character(250) :: MeshPath, ResultPath, InputPath, OutputPath, BafuFile, runid, InputFile, OutputFile
character(200) :: outfile, infile  ! (merge dup paths if any)
integer        :: syear, num_days, steps_per_day, n_out, steps_per_day_out
logical        :: backwards, put_back, twoD_tracking, save_column
integer        :: ios

! ---------------------------------------------------------------------
! Run bookkeeping: years/days/timesteps
! ---------------------------------------------------------------------
integer        :: years, year_idx, cyear, eyear         ! current/start/end year
integer        :: day, day_prev, first_day, end_day, num_days_all, days_remaining
integer        :: cnt_time_stp, step_day, step_day_out, step ! step is either -1 or 1, depending wether tracking is backwards or forwards
real(kind=8)   :: dt, timenew
character(4)   :: cyearchar, syearchar, eyearchar, prev_yearchar
character(3)   :: first_daychar!, drifter_id_char
character(2)   :: num_dayschar

! ---------------------------------------------------------------------
! Mesh / basis functions (loaded once, read-only afterwards)
! ---------------------------------------------------------------------
integer                     :: el2d
real(kind=8), allocatable   :: bafux(:,:), bafuy(:,:), voltriangle(:)
integer, allocatable        :: nmb_el_neigh(:), el_neigh_addresses(:,:)

! ---------------------------------------------------------------------
! Ocean file and fields at 3D nodes, current day and previous day
! ---------------------------------------------------------------------
character(250)            :: filename, filename_prev
integer                   :: nrec
real(kind=8), allocatable :: u_oce(:), v_oce(:), w_oce(:), t_oce(:), s_oce(:), sinking_oce(:)
real(kind=8), allocatable :: u_oce_prev(:), v_oce_prev(:), w_oce_prev(:), t_oce_prev(:), s_oce_prev(:)
real(kind=8)              :: weight_time

! ---------------------------------------------------------------------
! Drifter state carried between timesteps ("memory")
! ---------------------------------------------------------------------
integer, allocatable :: old_elem_idx(:), old_elem_nodes(:,:), runaways(:)
real(kind=8), allocatable :: old_lat(:), old_lon(:), old_depth(:)

! ---------------------------------------------------------------------
! Drifter input file (initial positions)
! ---------------------------------------------------------------------
integer :: num_drifters
real(kind=8), allocatable :: dlon(:), dlat(:), ddepth(:), dw_oce(:), dfirst_day(:)

! ---------------------------------------------------------------------
! Drifter output buffers — one row per output step, one col per drifter
! ---------------------------------------------------------------------
real, allocatable, dimension(:,:)   :: blon, blat, bdepth, btemp, bsalt
real, allocatable, dimension(:,:)   :: btemp_surface, bsalt_surface, mld
real, allocatable, dimension(:,:)   :: bu, bv, bw
real, allocatable, dimension(:,:)   :: du_dx, dv_dx, du_dy, dv_dy, du_dz, dv_dz
real, allocatable, dimension(:,:,:) :: btemp_column, bsalt_column
real(kind=8), allocatable           :: time(:)
integer, allocatable                :: bday(:), byear(:)

! ---------------------------------------------------------------------
! Per-timestep, per-drifter scratch: position & element search
! ---------------------------------------------------------------------
real(kind=8) :: phi1, lamda1, lat1, lam1, z1, re
integer      :: el_index, elem(3), index, clayer
integer      :: found_el, found_el_nodes(3)
real(kind=8) :: deptho, depthu, alpha
logical      :: inside
real(kind=8) :: phi2, lamda2, z2, cosd ! new position


! ---------------------------------------------------------------------
! Vertical interpolation onto the drifter (triangle + layer weights)
! ---------------------------------------------------------------------
integer      :: nodes_top(3), nodes_bot(3)
real(kind=8) :: weights(3)
real(kind=8) :: u_nodes(3), v_nodes(3)
real(kind=8) :: u_top, v_top, u_bot, v_bot, dz

! ---------------------------------------------------------------------
! Tetrahedron / 3D element search
! ---------------------------------------------------------------------
logical :: find_tetra
integer :: elemtetr(4), elem_collumn(3), elems(6)
real(kind=8) :: w(4), dist(4), tolerance, distsum

! ---------------------------------------------------------------------
! Tetrahedron Variable Averages (current an previous day for interpolation)
! ---------------------------------------------------------------------
real(kind=8) :: umean, vmean, wmean, tmean, smean
real(kind=8) :: umean_prev,vmean_prev,wmean_prev, tmean_prev, smean_prev


! ---------------------------------------------------------------------
! Mixed-layer depth via potential density (TEOS-10 / gsw)
! ---------------------------------------------------------------------
integer      :: ref_level
real(kind=8) :: sigma0_ref, sigma0_k, sigma0_km1
real(kind=8) :: T_ref, S_ref, T_k, S_k
real(kind=8) :: SA_ref, CT_ref, SA_k, CT_k, SA_km1, CT_km1
real(kind=8) :: p_ref, p_k, p_km1

real(kind=8) :: mld_calc
real(kind=8) :: depth_ref

! ---------------------------------------------------------------------
! Writing output
! ---------------------------------------------------------------------
integer :: start_rec_nc

! ---------------------------------------------------------------------
! Generic loop counters
! ---------------------------------------------------------------------
integer :: i, j, k, l, m, el, e, ll, kk, drifter, cnt, counter

! ---------------------------------------------------------------------
! Parameters
! ---------------------------------------------------------------------
real(kind=8), parameter :: delta_sigma_thresh = 0.03d0  ! kg/m3
real(kind=8), parameter :: fill_val = -999.d0


interface
   subroutine interpolate_surface_ts(lon_p, lat_p, elem, glon, glat, t_oce, s_oce, temp_out, salt_out)
      real(kind=8), intent(in) :: lon_p, lat_p
      integer,     intent(in) :: elem(3)
      real,        intent(in) :: glon(:), glat(:)
      real,        intent(in) :: t_oce(:), s_oce(:)
      real(kind=8), intent(out) :: temp_out, salt_out
   end subroutine interpolate_surface_ts
end interface

interface
  subroutine load_basis_functions(filename, bafux, bafuy, el2d, voltriangle)
    character(len=*), intent(in) :: filename
    real(kind=8), allocatable, intent(out) :: bafux(:,:), bafuy(:,:), voltriangle(:)
    integer, intent(out) :: el2d
  end subroutine load_basis_functions
end interface

interface
   subroutine calculate_triangle_weights(lon_p, lat_p, elem, glon, glat, weights)
      ! Inputs
      real(kind=8), intent(in) :: lon_p, lat_p          ! Particle longitude/latitude (radians)
      integer,     intent(in) :: elem(3)                ! 3 node indices element to calculate the weights for
      real, intent(in) :: glon(:), glat(:)              ! Longitude/latitude of all mesh nodes

      ! Outputs
      real(kind=8), intent(out) :: weights(3)           ! Interpolated weights

   end subroutine calculate_triangle_weights
end interface

interface
   subroutine calculate_triangle_weighted_mean(weights, elem, field, field_out)
   ! Inputs
   real(kind=8), intent(in) :: weights(3)
   integer,      intent(in) :: elem(3)
   real(kind=8), intent(in) :: field(:)
   ! Output
   real(kind=8), intent(out) :: field_out

   end subroutine calculate_triangle_weighted_mean
end interface


!!! -------------------------------------------------------------------------------------
!!! CONFIG 

! --- Namelist blocks ---
namelist /Paths/ MeshPath, ResultPath, InputPath, OutputPath, BafuFile, outfile, infile, runid
namelist /Parameters/ syear, num_days, steps_per_day, n_out, twoD_tracking, backwards, put_back, save_column
  
! --- Open and read the namelist file ---
open(unit=10, file='namelist.nml', status='old', action='read', iostat=ios)
if (ios /= 0) then
   print *, "Error opening namelist.nml, iostat = ", ios
   stop
end if

read(10, nml=Paths, iostat=ios)
if (ios /= 0) then
   print *, "Error reading Paths namelist, iostat = ", ios
   stop
end if

read(10, nml=Parameters, iostat=ios)
if (ios /= 0) then
   print *, "Error reading Parameters namelist, iostat = ", ios
   stop
end if

close(10)

! --- Trim strings ---
MeshPath   = adjustl(MeshPath)
ResultPath = adjustl(ResultPath)
InputPath  = adjustl(InputPath)
OutputPath = adjustl(OutputPath)
BafuFile   = adjustl(BafuFile)
outfile    = adjustl(outfile)
infile     = adjustl(infile)
runid      = adjustl(runid)

! --- Print nicely ---
print *, "=== Paths ==="
print *, "MeshPath   = '", trim(MeshPath), "'"
print *, "ResultPath = '", trim(ResultPath), "'"
print *, "InputPath  = '", trim(InputPath), "'"
print *, "OutputPath = '", trim(OutputPath), "'"
print *, "BafuFile   = '", trim(BafuFile), "'"
print *, "outfile    = '", trim(outfile), "'"
print *, "infile     = '", trim(infile), "'"
print *, "runid      = '", trim(runid), "'"

print *, "=== Parameters ==="
print *, "syear         = ", syear
print *, "num_days      = ", num_days
print *, "steps_per_day = ", steps_per_day
print *, "n_out         = ", n_out
print *, "twoD_tracking = ", twoD_tracking
print *, "backwards     = ", backwards
print *, "put_back      = ", put_back
print *, "save_column   = ", save_column

! Check if steps_per_day is divisible by n_out
if (mod(steps_per_day, n_out) /= 0) then
   print *, "Error: steps_per_day must be divisible by n_out"
   stop
end if

steps_per_day_out = steps_per_day / n_out
print *, "steps_per_day_out = ", steps_per_day_out

write(syearchar, '(i4)')syear

InputFile = trim(InputPath)//trim(infile)
OutputFile = trim(OutputPath)//trim(outfile)

!!! -------------------------------------------------------------------------------------
!!! READ MESH 
print*, 'rotated_grid: ',rotated_grid
if(rotated_grid) call calculate_rotate_matrix

! call calculate_rotate_matrix
call load_basis_functions(BafuFile, bafux, bafuy, el2d, voltriangle)
print*,'bafux: ', shape(bafux), 'bafuy: ', shape(bafuy), 'voltriangle: ', shape(voltriangle)

call read_2Dmesh(MeshPath)
call read_3Dmesh(MeshPath)

print*,'elem2D_nodes',elem2D_nodes(:,:2)
print*,"elem2d",elem2d
call mesh_scaling   ! long., lat. are transf. into rad
call find_layer_elem3d
call find_tetra_in_prism
 
write(*,*) 'reading mesh: DONE'

!!! ------------------------------------------------------------------------------------
!!! Read neighbour elements

allocate(nmb_el_neigh(elem2d))
allocate(el_neigh_addresses(16, elem2d))   ! 9 = max number of neighbours in file

open(unit=20, file='elem_neighbours_rotated.out', status='old', action='read')

do el = 1, elem2d
   read(20,*) nmb_el_neigh(el), (el_neigh_addresses(m,el), m=1,16)
end do
print*,"num of elemnt neighbours",nmb_el_neigh(1)

close(20)

! --- optional: print first few elements ---
do el = 1, min(5,elem2d)
   write(*,*) 'Element', el, 'has', nmb_el_neigh(el), 'neighbours:', &
               (el_neigh_addresses(m,el), m=1,nmb_el_neigh(el))
end do

print*,"Neighbours of element 519505:", el_neigh_addresses(:,519505), "so", nmb_el_neigh(519505)

!!! -------------------------------------------------------------------------------------
!!! READ DRIFTER INPUT FILE

open(10,file=InputFile)
read(10,*)num_drifters
allocate(dlon(num_drifters),dlat(num_drifters),ddepth(num_drifters),dw_oce(num_drifters),dfirst_day(num_drifters))
do i=1,num_drifters
  read(10,*) dlon(i),dlat(i),ddepth(i),dw_oce(i),dfirst_day(i)
enddo
close(10)
print*,'finished reading input file... number of drifters:', num_drifters
first_day = dfirst_day(1)
print*,'first_day = ', first_day, 'vel: ', dw_oce(1)

!!!--------------------------------------------------------------------------------------
! NOW compute years, since first_day is known
num_days_all = num_days
days_remaining = num_days_all
eyear = syear + ceiling(real(num_days_all - (365 - first_day + 1)) / 365.0) ! syear + ceiling(real(num_days_all + first_day - 1) / 365.0)
years = eyear - syear + 1
print*, "Total years to process: ", years, " (", syear, " to ", eyear, ")"

!!! -------------------------------------------------------------------------------------
!!! INITIALIZE OUTPUT FILE
write(syearchar, '(i4)')syear
write(eyearchar, '(i4)')eyear
write(first_daychar, '(i2)')first_day
write(num_dayschar, '(i2)')num_days
call init_output(OutputFile, num_drifters, num_days_all * steps_per_day_out, z_layers, depths, syear, steps_per_day, n_out, backwards, num_days, first_day, infile, save_column)


!!! -------------------------------------------------------------------------------------
!!! ALLOCATE ARRAYS

! 3D Fields
allocate(u_oce(nod3d), v_oce(nod3d), w_oce(nod3d), sinking_oce(nod3d), t_oce(nod3d), s_oce(nod3d))
allocate(u_oce_prev(nod3d), v_oce_prev(nod3d), w_oce_prev(nod3d), t_oce_prev(nod3d), s_oce_prev(nod3d))

! Time timformation (same for all drifters)
allocate(time(num_days_all * steps_per_day))
allocate(bday(num_days_all * steps_per_day))
allocate(byear(num_days_all * steps_per_day))

! Drifter Position and other variables
allocate(blon(steps_per_day, num_drifters))
allocate(blat(steps_per_day, num_drifters))
allocate(bdepth(steps_per_day, num_drifters))
allocate(btemp(steps_per_day, num_drifters))
allocate(bsalt(steps_per_day, num_drifters))
allocate(btemp_surface(steps_per_day, num_drifters))
allocate(bsalt_surface(steps_per_day, num_drifters))
allocate(du_dx(steps_per_day, num_drifters))
allocate(dv_dx(steps_per_day, num_drifters))
allocate(du_dy(steps_per_day, num_drifters))
allocate(dv_dy(steps_per_day, num_drifters))
allocate(du_dz(steps_per_day, num_drifters))
allocate(dv_dz(steps_per_day, num_drifters))

allocate(bu(steps_per_day, num_drifters))
allocate(bv(steps_per_day, num_drifters))
allocate(bw(steps_per_day, num_drifters))

allocate(btemp_column(steps_per_day, num_drifters, z_layers))
allocate(bsalt_column(steps_per_day, num_drifters, z_layers))
allocate(mld(steps_per_day, num_drifters))
 
! Memory of the previous time step
allocate(old_elem_idx(num_drifters))
allocate(old_elem_nodes(3, num_drifters))
allocate(old_lat(num_drifters))
allocate(old_lon(num_drifters))
allocate(old_depth(num_drifters))

allocate(runaways(num_drifters))
runaways = 0  ! 0 = active, 1 = stopped

cnt_time_stp = 0

time   = 0.
bday   = 0
byear  = 0


! Find reference level index closest to 10m (or use level 1 if shallow) (for MLD calculations later)
ref_level = 1
do k = 1, z_layers
    if (depths(k) >= 10.0d0) then   ! zbar: depth of level k, negative downward
        ref_level = k
        print*,"Picked model level", k, "as reference level for potential density calculations"
        !print*,depths
        exit
    end if
end do


!!! -------------------------------------------------------------------------------------
!!! MAIN LOOP

do year_idx = 1, years
   ! determine current year
   if (backwards) then
      cyear = syear - (year_idx - 1)
   else
      cyear = syear + (year_idx - 1)
   end if

   write(cyearchar, '(i4)')cyear ! convert year to character

   ! set filename for current year
   filename = trim(ResultPath)//trim(runid)//'.'//trim(cyearchar)//trim('.oce.mean.sub.nc')

   print*, 'Processing year ', cyear, ' with file ', trim(filename)

   ! 'day' book keeping
   if (year_idx == 1) then
      first_day = dfirst_day(1)
      num_days  = min(days_remaining, 365 - first_day + 1)  ! days left in first year
   else
      if (backwards) then
         first_day = 365
         num_days  = min(days_remaining, 365)
      else
         first_day = 1
         num_days  = min(days_remaining, 365)
      end if
   end if

   days_remaining = days_remaining - num_days
   print*, "Year ", cyear, ": processing ", num_days, " days from day ", first_day, ", ", days_remaining, " days remaining after this year"

   ! Prepare day loop parameters
   if (backwards) then
      end_day   = first_day - num_days + 1
      step      = -1
   else
      end_day   = first_day + num_days - 1
      step      = 1
      print*, "Processing from day ", first_day, " to day ", end_day
   end if

   ! Loop over days
   do day = first_day, end_day, step
      print*, "calculating day", day, ". days still to process for this year:", abs(end_day - day)
      print*, "  active drifters: ", count(runaways == 0), " / ", num_drifters, " | runaways: ", count(runaways == 1)

      if (backwards) then
         day_prev = day + 1
      else
         day_prev = day - 1
      end if

      ! make sure day_prev is within valid range [1, 365] and adjust filename_prev accordingly
      if (day_prev < 1) then ! set to last day of previous year
            day_prev = 365
            ! adjust filename_prev
            write(prev_yearchar, '(i4)') cyear - 1
            filename_prev = trim(ResultPath)//trim(runid)//'.'//trim(prev_yearchar)//trim('.oce.mean.sub.nc')

         else if (day_prev > 365) then ! set to first day of next year
            day_prev = 1
            ! adjust filename_prev
            write(prev_yearchar, '(i4)') cyear + 1
            filename_prev = trim(ResultPath)//trim(runid)//'.'//trim(prev_yearchar)//trim('.oce.mean.sub.nc')
         else
            filename_prev = filename
         end if


         ! If we are on the first day, read the data for that day and the previous day
         if (day == first_day) then
            print*, 'Loading data for first day ', day, ' and previous day ', day_prev
            call oce_input_netcdf(filename, day, u_oce, v_oce, w_oce, t_oce, s_oce, nod2d, nod3d)
            call oce_input_netcdf(filename_prev, day_prev, u_oce_prev, v_oce_prev, w_oce_prev, t_oce_prev, s_oce_prev, nod2d, nod3d)

            !!! rotate horizontal velocity vector of previous day (current day will happen below)
            do ll=1,nod3d
               call vector_r2g(u_oce_prev(ll), v_oce_prev(ll), lon3d_rad(ll), lat3d_rad(ll), 1)
            end do     

            ! change sign if backwards
            if (backwards) then 
               u_oce_prev= -u_oce_prev
               v_oce_prev= -v_oce_prev
               w_oce_prev= -w_oce_prev
            end if                                                                 
         
            if (twoD_tracking) then
               w_oce_prev = 0.0
            end if

         else
            ! copy data from last current day to previous day
            u_oce_prev = u_oce
            v_oce_prev = v_oce
            w_oce_prev = w_oce
            t_oce_prev = t_oce
            s_oce_prev = s_oce

            ! load the new current day data
            print*, 'Loading data for day ', day
            call oce_input_netcdf(filename, day, u_oce, v_oce, w_oce, t_oce, s_oce, nod2d, nod3d)

         end if

         !!! rotate horizontal velocity vector
         do ll=1,nod3d
           call vector_r2g(u_oce(ll),v_oce(ll) , lon3d_rad(ll), lat3d_rad(ll), 1)
         end do     

         if (backwards) then 
            u_oce= -u_oce
            v_oce= -v_oce
            w_oce= -w_oce
         end if                                                                 
         
         if (twoD_tracking) then
            w_oce = 0.0
         end if

      blon   = fill_val
      blat   = fill_val
      bdepth = fill_val
      btemp  = fill_val
      bsalt  = fill_val
      bu     = fill_val
      bv     = fill_val
      bw     = fill_val
      btemp_surface = fill_val
      bsalt_surface = fill_val
      btemp_column  = fill_val
      bsalt_column  = fill_val
      mld    = fill_val

      ! Loop over integration steps per day
      do step_day = 1, steps_per_day
         !print*, '  step ', step_day, ' of ', steps_per_day
         
         ! calculate weight between previous day and current day
         weight_time = real(step_day-1) / real(steps_per_day-1)

         timenew = -1
         
         cnt_time_stp = cnt_time_stp + 1
         
         dt = 86400/steps_per_day

         if (backwards) then 
            time(cnt_time_stp) =  day*86400 - (step_day-1) * dt - (year_idx-1) * 365 * 86400 ! <--- not sure if this is correct
         else 
            time(cnt_time_stp) = (year_idx-1) * 365 * 86400 + (day-1)*86400 + (step_day-1) * dt
         end if

         bday(cnt_time_stp)   = day 
         byear(cnt_time_stp)  = cyear

         do drifter = 1, num_drifters
            !print*, 'Start calculating drifter ', drifter

            ! 1. CHECK: Is this drifter already dead/stuck?
            if (runaways(drifter) == 1) then
               ! The value in b*(cnt_time_stp, drifter) remains -999.0
               !print*, "Drifter", drifter, " is dead, skip to next"
               cycle  ! Skip to the next drifter immediately
            endif

            !!! -------------------------------------------------------------------
            ! Determine drifter start position for this time step
            if (cnt_time_stp .EQ. 1) then
               ! first time step for this drifter => initialize drifter position
               lamda1 = dlon(drifter)
               if(lamda1.gt.180.)lamda1 = lamda1 - 360.
               if(lamda1.lt.-180.)lamda1 = lamda1 + 360.
               phi1 = dlat(drifter)
               z1   = ddepth(drifter) 
               !print*, "Starting drifter ", drifter, " at lon=", lamda1, " lat=", phi1, " depth=", z1
            else
               ! Grab previous position 
               lamda1   = old_lon(drifter)
               phi1     = old_lat(drifter)
               z1       = old_depth(drifter)
               el_index = old_elem_idx((drifter))
               elem     = old_elem_nodes(:, drifter)
               !print*, "Continuing drifter ", drifter, " at lon=", lamda1, " lat=", phi1, " depth=", z1, " previously in element ", el_index, "with nodes", elem
            end if
            !!! -------------------------------------------------------------------

            ! radius earth at given latitude (linear combination from equator and polar radius) 
            re = phi1/90. * rpol + (90. - phi1) / 90. * req

            ! convert to radians
            lam1 = pi/180. * lamda1
            lat1 = pi/180. * phi1
	 
            !!! -------------------------------------------------------------------
	         ! Find vertical position in the water column

            !   determine upper layer (clayer)
	         do i = 1, z_layers
	            if (z1.ge.depths(i)) clayer = i
	         enddo

            ! upper and lower layer depth
	         deptho = depths(clayer)
	         depthu = depths(clayer + 1)
	             
            !!! -------------------------------------------------------------------
            ! FIND 2D POSITION OF THE PARTICLE IN THE FESOM MESH

            if (cnt_time_stp .GT. 1) then ! If we are not in the first computing step start local search for current element
               ! First check if we are still in the old element
               !print*, "Checking old element", el_index, "with nodes", elem, "for position", lamda1, phi1
               call check_old_element(pi/180.*lamda1, pi/180.*phi1, elem, el_index, glon, glat, nod2d, inside)

               if (.NOT. inside) then ! If we are not in the same element anymore we probably moved into a neighbour
                  !print*,"Particle moved outside the old element, Checking neighbors..." 
                  call check_neighbour_elements(pi/180.*lamda1, pi/180.*phi1, el_index, elem2d_nodes, nmb_el_neigh(el_index), el_neigh_addresses(1:nmb_el_neigh(el_index),el_index), glon, glat, nod2d, elem2d, found_el, found_el_nodes, inside)

                  if (.NOT. inside) then ! If we are also not inside any of the neighbours we have to fall back to global search
                     print*,'!!! Couldnt find drifter', drifter, 'in the neighbor elements'
                     !print*,"el_index",el_index
                  else ! If the particle is inside one of the neighbour elements we can update our element index and nodes
                        !print*,'Found in neighbor element' 
                        !print*, "Updating new element index from", el_index , "to", found_el
                        el_index = found_el
                        elem = found_el_nodes
                  endif

               else ! If we remain in the old element we dont need to update any element/index
                  !print*, "Particle remains in old element" 
               endif

            endif


            if ((cnt_time_stp .EQ. 1) .OR. .NOT. inside) then ! If we are in the first computing step or we didnt find the particle anywhere -> do global search
            
               !print*,"!!! Start global search for matching element for drifter ", drifter
               call nearest_point(lam1,lat1,index,nod2d,glon,glat)

               !print*,'neareast point for',lamda1,phi1,'is',index,glon(index)*180/pi,glat(index)*180./pi
            

               call search_element(pi/180.*lamda1, pi/180.*phi1, elem,el_index,index,nod2d,elem2d,glon,glat,elem2d_nodes)
               !print*,"Global search found particle in 2D element with index", el_index
            
            end if


            !!! -------------------------------------------------------------------
            !   particle outside model domain (horizontally)
            if (sum(elem).lt.6) then
               print*,"drifter ", drifter, " reached the boundary.... "
               if (put_back .eq. .false.) then
                  print*,'.... drifter stops '
	               runaways(drifter) = 1 
                  cycle ! Skip to the next drifter immediately
               elseif  (put_back .eq. .true.) then ! 
               ! set point to center of element
               ! needs to be tested some more
               ! does not work well, particle gets stuck in a corner
                  call nearest_elem_level(lam1,lat1,elem,nod2d,elem2d,z_layers,glon,glat,elem2D_nodes,nod3D_below_nod2D,clayer + 1)
                  lamda1 = (glon(elem(1)) + glon(elem(2)) + glon(elem(3)))/3. * 180/pi
                  phi1   = (glat(elem(1)) + glat(elem(2)) + glat(elem(3)))/3. * 180/pi
                  lam1   = pi/180. * lamda1
                  lat1   = pi/180. * phi1
                  !print *, 'nearest point is:',lamda1, phi1
               endif
            endif
            
            !!!! -------------------------------------------------------------------
            ! CALCULATE 2D QUANTATIES

            ! --- build the 3-D node ids for the triangle at the two bracketing layers
            do l = 1, 3
               nodes_top(l) = nod3D_below_nod2D(clayer    , elem(l))   ! 3-D node just below 2-D node at upper face
               nodes_bot(l) = nod3D_below_nod2D(clayer + 1, elem(l))   ! 3-D node just below 2-D node at lower face
            end do

            ! --- vertical weight toward the lower face (depths are 0 at surface, positive downward)
            !     alpha = 0 at top face, 1 at bottom face
            alpha = (z1 - deptho) / (depthu - deptho)
            alpha = max(0.0d0, min(1.0d0, alpha))    ! safety clamp

            ! --- get u,v at z1 at the triangle's three corner nodes
            do l = 1, 3
               u_nodes(l) = (1.0d0 - alpha) * u_oce(nodes_top(l)) + alpha * u_oce(nodes_bot(l))
               v_nodes(l) = (1.0d0 - alpha) * v_oce(nodes_top(l)) + alpha * v_oce(nodes_bot(l))
            end do

            ! --- horizontal derivatives and relative vorticity 
            du_dx(step_day, drifter) = sum( u_nodes * bafux(:, el_index) )
            dv_dx(step_day, drifter) = sum( v_nodes * bafux(:, el_index) )
            du_dy(step_day, drifter) = sum( u_nodes * bafuy(:, el_index) )
            dv_dy(step_day, drifter) = sum( v_nodes * bafuy(:, el_index) ) 

            ! calculate weights for horizontal interpolation
            call calculate_triangle_weights(lam1, lat1, elem, glon, glat, weights)

            call calculate_triangle_weighted_mean(weights, nodes_top, u_oce, u_top)
            call calculate_triangle_weighted_mean(weights, nodes_top, v_oce, v_top)
            call calculate_triangle_weighted_mean(weights, nodes_bot, u_oce, u_bot)
            call calculate_triangle_weighted_mean(weights, nodes_bot, v_oce, v_bot)

            dz = max(depthu - deptho, 1.0d0)
            du_dz(step_day, drifter) = (u_bot - u_top) / dz
            dv_dz(step_day, drifter) = (v_bot - v_top) / dz

            ! Calculate surface T, S properties
            !call interpolate_surface_ts(lam1, lat1, elem, glon, glat, t_oce, s_oce, btemp_surface(step_day, drifter), bsalt_surface(step_day, drifter))  
            call calculate_triangle_weighted_mean(weights, elem, t_oce, btemp_surface(step_day, drifter))
            call calculate_triangle_weighted_mean(weights, elem, s_oce, bsalt_surface(step_day, drifter))

            !!!!!-------------------------------------------------------------------
            ! EXTRACT WATER COLUMN
            ! grab and output the whole vertical column above and below the particle and interpolate for each level and variable
            ! loop over all vertical levels and interpolate
            ! we need this for the MLD calculation and optionally as a save_coulumn output
            do i = 1, z_layers

               do l = 1, 3 ! find the three nodes of the current element at level i
                  elem_collumn(l) = nod3d_below_nod2d(i, elem(l))
               end do

               !print*, "Interpolating at vertical level ", i, " with nodes ", elem_collumn, "under nodes", elem

               ! interpolate temperature, salinity and velocity components at this level
               call calculate_triangle_weighted_mean(weights, elem_collumn, t_oce, btemp_column(step_day, drifter, i))
               call calculate_triangle_weighted_mean(weights, elem_collumn, s_oce, bsalt_column(step_day, drifter, i))

            enddo

            !!! -------------------------------------------------------------------
            ! Calculate mixed layer depth
            !!! -------------------------------------------------------------------
            call calculate_mld(step_day, drifter, ref_level, z_layers, &
                         steps_per_day, num_drifters, &
                         btemp_column, bsalt_column, depths, &
                         phi1, lamda1, fill_val, delta_sigma_thresh, mld)

            

            !!! -------------------------------------------------------------------
            !   the two elements that define the prism
	         do l = 1, 3
	            elems(l) = nod3d_below_nod2d(clayer, elem(l))
	            elems(l + 3) = nod3d_below_nod2d(clayer + 1, elem(l))
	         end do

            !!! -------------------------------------------------------------------
            !   particle outside model domain (vertically)
	        if (elems(4).eq.-999 .or. elems(5).eq.-999 .or. elems(6).eq.-999) then
               if (put_back .eq. .false.) then
	               print*,"drifter ", drifter, 'reached the bottom.... program stops'
	               runaways(drifter) = 1 
                  cycle  ! Skip to the next drifter immediately
               elseif  (put_back .eq. .true.) then ! 
               ! set point to center of element
               ! needs to be tested some more
               ! does not work well, particle gets stuck in a corner
                  call nearest_elem_level(lam1,lat1,elem,nod2d,elem2d,z_layers,glon,glat,elem2D_nodes,nod3D_below_nod2D,clayer + 1)
                  lamda1 = ( glon(elem(1)) + glon(elem(2)) + glon(elem(3)) )/3. * 180/pi
                  phi1   = ( glat(elem(1)) + glat(elem(2)) + glat(elem(3)) )/3. * 180/pi
                  lam1   = pi/180. * lamda1
                  lat1   = pi/180. * phi1
                  !print *, 'nearest point is:', lamda1, phi1
                  do l = 1, 3
	                  elems(l) = nod3d_below_nod2d(clayer, elem(l))
	                  elems(l + 3) = nod3d_below_nod2d(clayer + 1, elem(l))
	              end do

               end if
	        end if

     
            !!! -------------------------------------------------------------------
            !   searching the tetrahedra in which 3d point lies
            call search_tetrahedron(lamda1, phi1, z1, deptho, depthu, el_index, clayer,elem3d_nodes,elem2D,z_layers,elem3d, nod3d, elemtetr, depth3d, lon3d, lat3d, tetra_in_prism,find_tetra)

            !!! -------------------------------------------------------------------
            !   in case tetrahedra not found, we repeat the procedure
            cnt = 0 
            do while (find_tetra .eq. .false. .and. cnt .le.100)  
               cnt = cnt + 1
               !print*, " ﬁ"
            !   ! 
               tolerance    = 0.0001 ! 10 m tolerance 
               lamda1 = lamda1 + tolerance ! change lamda1 a little bit
               lam1   = pi/180. * lamda1

               ! Check first if the particle is still in the old element
               call check_old_element(pi/180.*lamda1, pi/180.*phi1, elem, el_index, glon, glat, nod2d, inside)

               if (.NOT. inside) then ! If we are not in the same element anymore we probably moved into a neighbour
                  !if (cnt < 2) then
                  !   print*,"Tetrahedron tolerance moved particle to another element, Checking neighbors... (And continue searching in while loop)" 
                  ! endif
                  call check_neighbour_elements(pi/180.*lamda1, pi/180.*phi1, el_index, elem2d_nodes, nmb_el_neigh(el_index), el_neigh_addresses(1:nmb_el_neigh(el_index),el_index), glon, glat, nod2d, elem2d, found_el, found_el_nodes, inside)

                  if (.NOT. inside) then ! If we are also not inside any of the neighbours we have to fall back to global search
                     print*,'!!! Couldnt find particle in the neighbor elements, Searching globally'
                     call nearest_point(lam1, lat1, index, nod2d, glon, glat)
                     call search_element(pi/180.*lamda1, pi/180.*phi1, elem, el_index, index, nod2d, elem2d, glon, glat, elem2d_nodes)
                  else ! If the particle is inside one of the neighbour elements we can update our element index and nodes
                        !print*,'Found in neighbor element' 
                        !print*, "Updating new element index from", el_index , "to", found_el
                        el_index = found_el
                        elem = found_el_nodes
                  endif
               endif

               do l = 1, 3
	               elems(l) = nod3d_below_nod2d(clayer, elem(l))
	               elems(l + 3) = nod3d_below_nod2d(clayer + 1, elem(l))
	            enddo

               call search_tetrahedron(lamda1, phi1, z1, deptho, depthu, el_index, clayer,elem3d_nodes,elem2D,z_layers,elem3d, nod3d, elemtetr, depth3d, lon3d, lat3d, tetra_in_prism,find_tetra)

            end do

            if (find_tetra .eq. .false.) then
               print*,"!!! Tetrahedron tolerance search did not work !!!"
            endif

            !!! -------------------------------------------------------------------
            !   calculate distance of point to all four nodes that define the tetrahedra
            do kk = 1, 4
               call distance_3d(lamda1*pi/180., phi1*pi/180., lon3d(elemtetr(kk))*pi/180., lat3d(elemtetr(kk))*pi/180., z1, depth3d(elemtetr(kk)), dist(kk))
	         enddo
  
            !!! -------------------------------------------------------------------
            !   computing weights
            w = 0.
	         do kk = 1, 4
               if(dist(kk)==0.)then
                  w(kk) = 1.
                  goto 222
               endif
            end do

            distsum = dist(2)/dist(1)+dist(3)/dist(1)+dist(4)/dist(1)+dist(1)/dist(2)+dist(3)/dist(2)+dist(4)/dist(2)+&
            dist(1)/dist(3)+dist(2)/dist(3)+dist(4)/dist(3)+dist(1)/dist(4)+dist(2)/dist(4)+dist(3)/dist(4)
            !   print*,'distsum',distsum
            w(1) = (dist(2)/dist(1)+dist(3)/dist(1)+dist(4)/dist(1))/distsum
            w(2) = (dist(1)/dist(2)+dist(3)/dist(2)+dist(4)/dist(2))/distsum
            w(3) = (dist(1)/dist(3)+dist(2)/dist(3)+dist(4)/dist(3))/distsum
	         w(4) = (dist(1)/dist(4)+dist(2)/dist(4)+dist(3)/dist(4))/distsum
	 
	 
222         umean = 0.
            vmean = 0.
	         wmean = 0.
            tmean = 0.
            smean = 0.

            umean_prev = 0.
            vmean_prev = 0.
            wmean_prev = 0.
            tmean_prev = 0.
            smean_prev = 0.

	        !print*, 'w',w, sum(w)
      
            !!! -------------------------------------------------------------------
            !   computing mean velcity / temp / salt from the four nodes 
            do kk=1,4
               umean = umean + w(kk) * u_oce(elemtetr(kk))
               vmean = vmean + w(kk) * v_oce(elemtetr(kk))
               wmean = wmean + w(kk) * w_oce(elemtetr(kk))
               tmean = tmean + w(kk) * t_oce(elemtetr(kk))
               smean = smean + w(kk) * s_oce(elemtetr(kk))
               umean_prev = umean_prev + w(kk) * u_oce_prev(elemtetr(kk))
               vmean_prev = vmean_prev + w(kk) * v_oce_prev(elemtetr(kk))
               wmean_prev = wmean_prev + w(kk) * w_oce_prev(elemtetr(kk))
               tmean_prev = tmean_prev + w(kk) * t_oce_prev(elemtetr(kk))
               smean_prev = smean_prev + w(kk) * s_oce_prev(elemtetr(kk))

	           !print*, u_oce(elemtetr(kk)), v_oce(elemtetr(kk)), w_oce(elemtetr(kk))
            end do

            ! linear interpolation in time between previous day and current day
            umean = (1.0 - weight_time) * umean_prev + weight_time * umean
            vmean = (1.0 - weight_time) * vmean_prev + weight_time * vmean
            wmean = (1.0 - weight_time) * wmean_prev + weight_time * wmean
            tmean = (1.0 - weight_time) * tmean_prev + weight_time * tmean
            smean = (1.0 - weight_time) * smean_prev + weight_time * smean

            if (cnt_time_stp == 1 ) then
               blon(1, drifter)   = lamda1
               blat(1, drifter)   = phi1
	            bdepth(1, drifter) = z1
               btemp(1, drifter)  = tmean
               bsalt(1, drifter)  = smean

               old_lat(drifter) = phi1
               old_lon(drifter) = lamda1
               old_depth(drifter) = z1
               old_elem_idx(drifter) = el_index
               old_elem_nodes(:, drifter) = elem
               !bday(1, drifter)   = d ! update later because drifter independent
               !byear(1, drifter)  = syear - (year_idx - 1)
               !time(1)   = 0
               cycle
            end if

            !!! -------------------------------------------------------------------
            !   computing new position (with Euler)
            !print*,'compute new position'
            !if (j.eq.2) then
            !  call cpu_time(start_time1)
            !end if

            ! new depth: 
	         z2 = z1 - wmean * dt
	         if (z2.lt.0) z2 = 0
  
            ! new latitude: 
            phi2 = vmean * dt / re / pi * 180. + phi1
   
            ! new longitude:    
            cosd = (abs(umean) * dt / re)
            cosd = cos(cosd)
            lamda2 = cosd - sin(pi/180.*phi1) * sin(pi/180.*phi1)
            lamda2 = lamda2 / cos(pi/180.*phi1) / cos(pi/180.*phi1)

            if(lamda2.gt.1.) lamda2 = 1.
            if(umean.lt.0.)then
               lamda2 = -acos(lamda2) * 180./pi + lamda1
            else
               lamda2 = acos(lamda2) * 180./pi + lamda1
            end if
            if(lamda2.gt.180.)lamda2 = lamda2 - 360.
            if(lamda2.lt.-180.)lamda2 = lamda2 + 360.

            ! save new locations
            blon(step_day, drifter)   = lamda2
            blat(step_day, drifter)   = phi2
	         bdepth(step_day, drifter) = z2
            !print*,"Drifter", drifter, "moved to", lamda2, phi2, z2
   
            ! save properties at old position
            if (cnt_time_stp.gt.2) then
               btemp(step_day-1, drifter)=tmean
               bsalt(step_day-1, drifter)=smean
               bu(step_day-1, drifter)   = umean
               bv(step_day-1, drifter)   = vmean
               bw(step_day-1, drifter)   = wmean
            end if

            ! Update position
            z1 = z2 ! do I still need this in this loop structure
            phi1 = phi2
            lamda1 = lamda2

            ! save element index and nodes for the next time step
            old_lat(drifter) = phi1
            old_lon(drifter) = lamda1
            old_depth(drifter) = z1
            old_elem_idx(drifter) = el_index
            old_elem_nodes(:, drifter) = elem
            

            if (twoD_tracking == 0) then
               if (bdepth(step_day, drifter) == 0.) then     ! stop calculating when reaching the surface 
                  runaways(drifter) = 1
                  print*, "!!! Drifter", drifter, "jumped out of the water (it wanted to be free :D)"
                  cycle
               end if
            end if
         

         end do  ! loop over drifter
         
      end do ! step_day loop

      ! This tells NetCDF: "Start writing today's records at this position"
      print*, "Writing output for day ", day, " at record index ", start_rec_nc
      start_rec_nc = ((cnt_time_stp - steps_per_day) / n_out) + 1

      if (save_column) then
         call write_output_byDay(OutputFile, num_drifters, steps_per_day_out, z_layers, start_rec_nc, &
                                 time(cnt_time_stp-steps_per_day+1 : cnt_time_stp : n_out), & ! Slice global array
                                 blon(1:steps_per_day:n_out, :), &               ! Slice daily array ! Time is 1st, Drifter is 2nd
                                 blat(1:steps_per_day:n_out, :), &
                                 bdepth(1:steps_per_day:n_out, :), &
                                 btemp(1:steps_per_day:n_out, :), &
                                 bsalt(1:steps_per_day:n_out, :), &
                                 bday(cnt_time_stp-steps_per_day+1 : cnt_time_stp : n_out), &
                                 byear(cnt_time_stp-steps_per_day+1 : cnt_time_stp : n_out), &
                                 btemp_surface(1:steps_per_day:n_out, :), &
                                 bsalt_surface(1:steps_per_day:n_out, :), &
                                 du_dx(1:steps_per_day:n_out, :), &
                                 dv_dx(1:steps_per_day:n_out, :), &
                                 du_dy(1:steps_per_day:n_out, :), &
                                 dv_dy(1:steps_per_day:n_out, :), &
                                 du_dz(1:steps_per_day:n_out, :), &
                                 dv_dz(1:steps_per_day:n_out, :), &
                                 bu(1:steps_per_day:n_out, :), &
                                 bv(1:steps_per_day:n_out, :), &
                                 bw(1:steps_per_day:n_out, :), &
                                 mld(1:steps_per_day:n_out, :), &
                                 save_column, &
                                 btemp_column(1:steps_per_day:n_out, :, :), &     ! Slice 3D arrays ! Time is 1st, then Drifter, then Layers
                                 bsalt_column(1:steps_per_day:n_out, :, :))

      else
         call write_output_byDay(OutputFile, num_drifters, steps_per_day_out, z_layers, start_rec_nc, &
                                 time(cnt_time_stp-steps_per_day+1 : cnt_time_stp : n_out), & ! Slice global array
                                 blon(1:steps_per_day:n_out, :), &               ! Slice daily array ! Time is 1st, Drifter is 2nd
                                 blat(1:steps_per_day:n_out, :), &
                                 bdepth(1:steps_per_day:n_out, :), &
                                 btemp(1:steps_per_day:n_out, :), &
                                 bsalt(1:steps_per_day:n_out, :), &
                                 bday(cnt_time_stp-steps_per_day+1 : cnt_time_stp : n_out), &
                                 byear(cnt_time_stp-steps_per_day+1 : cnt_time_stp : n_out), &
                                 btemp_surface(1:steps_per_day:n_out, :), &
                                 bsalt_surface(1:steps_per_day:n_out, :), &
                                 du_dx(1:steps_per_day:n_out, :), &
                                 dv_dx(1:steps_per_day:n_out, :), &
                                 du_dy(1:steps_per_day:n_out, :), &
                                 dv_dy(1:steps_per_day:n_out, :), &
                                 du_dz(1:steps_per_day:n_out, :), &
                                 dv_dz(1:steps_per_day:n_out, :), &
                                 bu(1:steps_per_day:n_out, :), &
                                 bv(1:steps_per_day:n_out, :), &
                                 bw(1:steps_per_day:n_out, :), &
                                 mld(1:steps_per_day:n_out, :), &
                                 save_column)
      end if



      if (bday(cnt_time_stp) == 365) then ! read next year when reaching Jan 1st
         exit
      end if 

   end do ! day loop

end do ! year loop


end program trajectories
!
!=======================================================================
!
!
!=======================================================================
!
subroutine oce_input_netcdf(filename,nrec,u_oce,v_oce,w_oce,t_oce,s_oce,nod2d,nod3d)
  ! read fesom fields for ocean dynamics and active tracer variables
  
 
  implicit none

#include "netcdf.inc" 

  integer                   :: status, ncid, j, dimid_rec, nrec
  integer                   :: tra_varid(2)
  integer                   :: u_varid, v_varid, w_varid
  integer                   :: istart(2), icount(2), n3
  character(1)              :: trind
  real(kind=8), allocatable :: aux2(:), aux3(:) 
  integer                   :: nod2d,nod3d
  real(kind=8), DIMENSION(nod3d) :: u_oce, v_oce, w_oce, t_oce, s_oce
  character(200)            :: filename
  real(kind=8)              :: timenew


   allocate(aux2(nod2D), aux3(nod3D)) 

print*, 'nrec ',nrec
print*, 'filename ',filename

  ! open files
  status = nf_open(filename, nf_nowrite, ncid)
  if (status .ne. nf_noerr) call handle_err(status)

  ! inquire variable id
!!!  status=nf_inq_varid(ncid, 'ssh', ssh_varid)
!!!  if (status .ne. nf_noerr) call handle_err(status)
  status=nf_inq_varid(ncid, 'u', u_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status=nf_inq_varid(ncid, 'v', v_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status=nf_inq_varid(ncid, 'w', w_varid)
  if (status .ne. nf_noerr) call handle_err(status)
  status=nf_inq_varid(ncid, 'temp', tra_varid(1))
  if (status .ne. nf_noerr) call handle_err(status)
  status=nf_inq_varid(ncid, 'salt', tra_varid(2))
  if (status .ne. nf_noerr) call handle_err(status)

  ! read variables

  ! 2d fields
 !!! istart=(/1,nrec/)
 !!! icount=(/nod2d, 1/)
!!!  status=nf_get_vara_double(ncid, ssh_varid, istart, icount, aux2) 
 !!! if (status .ne. nf_noerr) call handle_err(status)
 !!! ssh=aux2(myList_nod2D)         

  ! 3d fields
  istart=(/1,nrec/)
  icount=(/nod3d, 1/)

  status=nf_get_vara_double(ncid, u_varid, istart, icount, aux3)
  if (status .ne. nf_noerr) call handle_err(status)
  u_oce(1:nod3d)=aux3   

  status=nf_get_vara_double(ncid, v_varid, istart, icount, aux3) 
  if (status .ne. nf_noerr) call handle_err(status)
  v_oce(1:nod3d)=aux3

  status=nf_get_vara_double(ncid, w_varid, istart, icount, aux3) 
  if (status .ne. nf_noerr) call handle_err(status)
  w_oce(1:nod3d)=aux3

  status=nf_get_vara_double(ncid, tra_varid(1), istart, icount, aux3) 
  if (status .ne. nf_noerr) call handle_err(status)
  t_oce(1:nod3d)=aux3

  status=nf_get_vara_double(ncid, tra_varid(2), istart, icount, aux3) 
  if (status .ne. nf_noerr) call handle_err(status)
  s_oce(1:nod3d)=aux3 

  status=nf_close(ncid)
  if (status .ne. nf_noerr) call handle_err(status)

  deallocate(aux3, aux2)   

end subroutine oce_input_netcdf
!
!=======================================================================
!
subroutine handle_err(errcode)
 !!! use g_parfe
  implicit none
  
#include "netcdf.inc" 
  
  integer errcode
  
  write(*,*) 'Error: ', nf_strerror(errcode)
  !!!call par_ex
  stop
end subroutine handle_err
!
!=======================================================================
!
subroutine point_in_element(point,a,b,c,flag)
!function arguments contain geographical coordinates (lat,lon)
  real(kind=8),dimension(2):: point,a,b,c
  logical:: flag,flaga,flagb,flagc

  call sameside(point,a,b,c,flaga)
  call sameside(point,b,a,c,flagb)
  call sameside(point,c,a,b,flagc)

  if(flaga.and.flagb.and.flagc)then 
     flag=.true.
  else
     flag=.false.
  endif
  return
end subroutine point_in_element
!
!=======================================================================
!
subroutine nearest_point(lon,lat,index,nod2d,glon,glat)

  implicit none
  integer index,j,nod2d
  real(kind=8):: lat,lon!in RAD
  real(kind=8)::dif1,phi1,phi2,lam1,lam2,dif4,dif2,pi,lon1,lat1
  real(kind=8),dimension(nod2d)::glon,glat

  lon1=lon
  lat1=lat
  pi=4.*atan(1.)
  index=0
  dif1=1000000.0
  do j=1,nod2d
     if(glon(j).lt.0.)then
        lam1=2.*pi+glon(j)
     else
        lam1=glon(j)
     endif
     if(lon1.lt.0.)then
        lam2=2.*pi+lon
     else
        lam2=lon
     endif
     phi1=glat(j)
     phi2=lat
     dif2=dsin(dble(phi1))*dsin(dble(phi2))+dcos(dble(phi1))*dcos(dble(phi2))*dcos(dble(lam2-lam1))
     dif4=dacos(dble(dif2))
     if(dif4.lt.dif1)then
        dif1=dif4
        index=j
     endif
     if(dif1.eq.0.)then
        index=j
        goto 99
     endif
  enddo
99 continue
  return
end subroutine nearest_point
!
!=======================================================================
!
subroutine nearest_elem_level(lon,lat,index,nod2d,elem2d,z_layers,glon,glat,elem2D_nodes,nod3D_below_nod2D,depthu)

  ! if the particle is outside the mesh at level depthu, we search for the nearest element 
  ! in that level, and set the particle location to the element center

  implicit none
  integer index(3),j,nod2d,elem2d,z_layers,depthu
  real(kind=8):: lat,lon !in RAD
  real(kind=8)::dif1,phi1,phi2,lam1,lam2,dif4,dif2,pi,lon1,lat1,lon_cen,lat_cen
  real(kind=8),dimension(nod2d)::glon,glat
  integer,dimension(3,elem2d)::elem2D_nodes
  integer(KIND=4), dimension(z_layers,nod2D) :: nod3D_below_nod2D  

  lon1=lon
  lat1=lat
  pi=4.*atan(1.)
  index=0
  dif1=1000000.0
  do j=1,elem2d
     if  ( (nod3D_below_nod2D(depthu,elem2D_nodes(1,j)).gt.0) .and. (nod3D_below_nod2D(depthu,elem2D_nodes(2,j)).gt.0) .and. (nod3D_below_nod2D(depthu,elem2D_nodes(3,j)).gt.0) ) then ! only check nodes in layer depthu
        
        lon_cen=(glon(elem2D_nodes(1,j)) + glon(elem2D_nodes(2,j)) + glon(elem2D_nodes(3,j)))/3.
        lat_cen=(glat(elem2D_nodes(1,j)) + glat(elem2D_nodes(2,j)) + glat(elem2D_nodes(3,j)))/3.
        if(lon_cen.lt.0.)then
           lam1=2.*pi+lon_cen
        else
           lam1=lon_cen 
        endif
        if(lon1.lt.0.)then
           lam2=2.*pi+lon
        else
           lam2=lon
        endif
        phi1=lat_cen 
        phi2=lat
        dif2=dsin(dble(phi1))*dsin(dble(phi2))+dcos(dble(phi1))*dcos(dble(phi2))*dcos(dble(lam2-lam1))
        dif4=dacos(dble(dif2))
        if(dif4.lt.dif1)then
           dif1=dif4
           index(1)=elem2D_nodes(1,j)
           index(2)=elem2D_nodes(2,j)
           index(3)=elem2D_nodes(3,j)
        endif
        if(dif1.eq.0.)then
           index(1)=elem2D_nodes(1,j)
           index(2)=elem2D_nodes(2,j)
           index(3)=elem2D_nodes(3,j)
           goto 99
        endif
     end if
  enddo
99 continue
  return
end subroutine nearest_elem_level
!
!=======================================================================
!
subroutine sameside(p1,p2,a,b,flag)
! arguments contain geographical coordinates (lat,lon)
! checks if p1 and p2 are on the same side of line of a b
  real(kind=8),dimension(2):: p1,p2,a,b,diff1,diff2,diff3
  logical:: flag
  real(kind=8),dimension(2):: cp1,cp2
  real(kind=8)::dp
  real(kind=8)::one_real

one_real=1.

  flag=.false.
  diff1(:)=b(:)-a(:)
  diff2(:)=p1(:)-a(:)
  diff3(:)=p2(:)-a(:)
  !wechsel der lon abfangen!!!
  if(abs(diff1(1)).gt.180.)diff1(1)=-1*dsign(dble(1.0),dble(diff1(1)))*360.+diff1(1)
  if(abs(diff2(1)).gt.180.)diff2(1)=-1*dsign(dble(1.0),dble(diff2(1)))*360.+diff2(1)
  if(abs(diff3(1)).gt.180.)diff3(1)=-1*dsign(dble(1.0),dble(diff3(1)))*360.+diff3(1)
!cp == cross product
  cp1(1)=diff1(2)*diff2(1)-diff1(1)*diff2(2)
  cp1(2)=diff1(1)*diff2(2)-diff1(2)*diff2(1)
  cp2(1)=diff1(2)*diff3(1)-diff1(1)*diff3(2)
  cp2(2)=diff1(1)*diff3(2)-diff1(2)*diff3(1)
!dp == dot product
  dp= cp1(1)*cp2(1)+cp1(2)*cp2(2)
  if(dp .ge.0._8) flag=.true.
  return
end subroutine sameside
!
!=======================================================================
!
subroutine check_old_element(particle_lon, particle_lat, element_nodes, element_index, grid_lon, grid_lat, nod2d, flag)

   implicit none
   integer                        :: i, element_index, element_nodes(3), nod2d    ! node-index for looping, element-index, 2D nodes of the element, number of 2D nodes in the grid
   real(kind=8)                   :: particle_lon, particle_lat, rad_inv, pi      ! longitude/latitude of the particle to check, conversion parameter, Pi
   real(kind=8), dimension(3,2)   :: nodes_geo                                    ! geographical position of the 3 element nodes
   real(kind=8), dimension(nod2d) :: grid_lon, grid_lat                           ! longitude/latitude of each 2D node
   real(kind=8), dimension(2)     :: point                                        ! converted position of the particle
   logical                        :: flag                                         ! Bool if the point is in the element (output)

   pi      = 4. * atan(1.)
   rad_inv = 180. / (4. * atan(1.))

   flag = .false. 

   ! Convert lon & lat from the particle position radians -> degrees 
   if(particle_lon.lt.0.)then
      point(1) = 2. * rad_inv * pi + rad_inv * particle_lon
   else
      point(1) = rad_inv * particle_lon
   endif
   point(2) = rad_inv * particle_lat

   ! Convert node position radians -> degrees 
   do i = 1, 3
      if (grid_lon(element_nodes(i)).lt.0.)then
         nodes_geo(i, 1) = 2. * rad_inv * pi + rad_inv * grid_lon(element_nodes(i))
      else 
         nodes_geo(i, 1) = rad_inv * grid_lon(element_nodes(i))
      endif

      nodes_geo(i, 2) = rad_inv * grid_lat(element_nodes(i))
   enddo

   call point_in_element(point, nodes_geo(1,:), nodes_geo(2,:), nodes_geo(3,:), flag)
          
  return
end subroutine check_old_element
!
!=======================================================================
!
! ---------------------------------------------------------------------
! check_neighbour_elements
! Check all neighbour elements of `current_el` for containment of the point.
! Returns the neighbour element index if found, or 0 if none matched.
! ---------------------------------------------------------------------
subroutine check_neighbour_elements(particle_lon, particle_lat, current_el, elem2d_nodes, elem_nmb, neighbours, grid_lon, grid_lat, nod2d, elem2d, found_el, found_el_nodes, is_in)
  implicit none

  ! inputs
  real(kind=8) :: particle_lon, particle_lat
  integer      :: current_el               ! element whose neighbours we test
  integer      :: elem2d
  integer      :: elem2d_nodes(3, elem2d)       ! element connectivity (3 x elem2d)
  integer      :: elem_nmb                 ! number of neighbours of the current element
  integer      :: neighbours(elem_nmb) ! neighbour ids 
  integer      :: nod2d
  integer      :: found_el_nodes(3)
  real(kind=8) :: grid_lon(nod2d), grid_lat(nod2d)

  ! output
  integer      :: found_el                 ! found element index (0 => not found)

  ! locals
  integer :: i, try_el
  logical :: is_in


  ! Default: not found
  found_el = 0
  found_el_nodes = 0
  is_in = 0

  ! Guard: if current_el is out of range or has zero neighbours, skip
  !print*, "Searching ",elem_nmb, "neighbour elements of element", current_el
  !print*, "which are...", neighbours
 
  if (current_el < 1) return
  if (elem_nmb <= 0) return
  
  ! Loop neighbours in the same style as you built them
   do i = 1, elem_nmb
      
      try_el = neighbours(i)
      !print*, "Trying element", try_el

      if (try_el <= 0) cycle        ! padded or invalid entry (e.g. -999)

      !print*, "These are the input arguments", particle_lon, particle_lat, elem2d_nodes(:,try_el), try_el, is_in

     ! call your existing routine that checks if the particle is in the element
      call check_old_element(particle_lon, particle_lat, elem2d_nodes(:,try_el), try_el, grid_lon, grid_lat, nod2d, is_in)

     if (is_in) then
        found_el = try_el
         found_el_nodes = elem2d_nodes(:,try_el)
        !print*, "Particle was found in neighbour element", found_el 
        return
     end if
  end do

  ! nothing found among neighbours => return found_el = 0#
  !print*, "NO Particle was found in neighbour element", found_el 
  return
end subroutine check_neighbour_elements


!
!=======================================================================
!
subroutine search_element(lon,lat,index,el_index,index1,nod2d,elem2d,glon,glat,elem2D_nodes )

  implicit none
  integer :: index(3),j,k,index1,nod2d,elem2d, el_index
  real(kind=8):: lat,lon,onecount,rad_inv,pi
  real(kind=8),dimension(2)::a,b,c,point
  real(kind=8),dimension(nod2d)::glon,glat
  integer,dimension(3,elem2d)::elem2D_nodes
  logical:: flag
  pi=4.*atan(1.)
  rad_inv=180./(4.*atan(1.))
  index(:)=0
  do j=1,elem2D
     do k=1,3
        if(elem2D_nodes(k,j).eq.index1)then

           if(glon(elem2D_nodes(1,j)).lt.0.)then
              a(1)=2.*rad_inv*pi+rad_inv*glon(elem2D_nodes(1,j))
           else
              a(1)=rad_inv*glon(elem2D_nodes(1,j))
           endif

           if(glon(elem2D_nodes(2,j)).lt.0.)then
              b(1)=2.*rad_inv*pi+rad_inv*glon(elem2D_nodes(2,j))
           else
              b(1)=rad_inv*glon(elem2D_nodes(2,j))
           endif

           if(glon(elem2D_nodes(3,j)).lt.0.)then
              c(1)=2.*rad_inv*pi+rad_inv*glon(elem2D_nodes(3,j))
           else
              c(1)=rad_inv*glon(elem2D_nodes(3,j))
           endif
 
           a(2)=rad_inv*glat(elem2D_nodes(1,j))
           b(2)=rad_inv*glat(elem2D_nodes(2,j))

           c(2)=rad_inv*glat(elem2D_nodes(3,j))
           if(lon.lt.0.)then
              point(1)=2.*rad_inv*pi+rad_inv*lon
           else
              point(1)=rad_inv*lon
           endif

           point(2)=rad_inv*lat
           flag=.false.

           call point_in_element(point,a,b,c,flag)

           if(flag)then
              index(1)=elem2D_nodes(1,j)
              index(2)=elem2D_nodes(2,j)
              index(3)=elem2D_nodes(3,j)
              el_index=j
              onecount=onecount+1
              return
           endif
        endif
     enddo
  enddo
  return
end subroutine search_element
!
!=======================================================================
!
subroutine search_tetrahedron(lon, lat, depth, deptho, depthu, el_index, clayer,elem3d_nodes,elem2d,z_layers,elem3d, nod3d, elemtetr, depth3d, lon3d, lat3d,tetra_in_prism,find_tetra)

  use o_DATA_TYPES

   implicit none
   

   integer:: elem2d, elem3d, nod3d, elemtetr(4),ii,j, k, l, counter, node, sort1(4), sort2(4), c1, c2, ccwt(4),clayer,z_layers, el_index
   real(kind=8):: lon, lat, depth, deptho, depthu, verh1, verh2, loncoord(4), latcoord(4)
   integer,dimension(4,elem3d):: elem3d_nodes
   real(kind = 8), dimension(nod3d):: depth3d, lon3d, lat3d
   logical     :: find_tetra
   integer       :: tetr

   type(addresstype2), dimension(z_layers,elem2d) :: tetra_in_prism

   verh1 = (depth - deptho)/(depthu - deptho)
   verh2 = (depthu - depth)/(depthu - deptho)

   find_tetra = .true.

   do j=1, tetra_in_prism(clayer,el_index)%nmb 

      tetr=tetra_in_prism(clayer,el_index)%addresses(j)

      elemtetr = elem3d_nodes(:,tetr)

      
      ! Liegt die Position des Drifters im gefundenen Tetraeder?
      !x- und y-Werte des Dreiecks bestimmen
      counter = 0
      do k = 1, 4
         if (depth3d(elemtetr(k)).eq.deptho) then
	        counter = counter + 1
	        node = k
	     endif
      enddo
      
      if (counter.eq.1) then ! blue tetrahedra
         l = 0
         do k = 1, 4
	        if (k.ne.node) then
	           l = l + 1
	           loncoord(l) = lon3d(elemtetr(node)) + ((lon3d(elemtetr(k)) - lon3d(elemtetr(node))) * verh1)
	           latcoord(l) = lat3d(elemtetr(node)) + ((lat3d(elemtetr(k)) - lat3d(elemtetr(node))) * verh1)
	        endif
         enddo
      endif
      
      if (counter.eq.2) then ! green tetrahedra
         c1 = 1
	     c2 = 3
	    do k = 1, 4
	       if (depth3d(elemtetr(k)).eq.deptho) then
	          sort1(c1) = elemtetr(k)
	          c1 = c1 + 1
	       else
	          sort1(c2) = elemtetr(k)
	          c2 = c2 + 1
	       endif
	    enddo
	 ! the two nodes that are on top of each other:
	    do k = 1, 2
	       do l = 3, 4
	          if ((lon3d(sort1(k))).eq.lon3d(sort1(l)).and.(lat3d(sort1(k))).eq.lat3d(sort1(l))) then
	             sort2(2) = sort1(k)
		         sort2(3) = sort1(l)
	             c1 = k
		         c2 = l
	          endif
	       enddo
	    enddo
	 

	    if (c1.eq.1) then
	       sort2(1) = sort1(2)
	    else
	       sort2(1) = sort1(1)
	    endif
	 
	    if (c2.eq.4) then
	       sort2(4) = sort1(3)
	    else
	       sort2(4) = sort1(4)
	    endif
	 
	    loncoord(1) = lon3d(sort2(1)) + ((lon3d(sort2(4)) - lon3d(sort2(1))) * verh1)
	    latcoord(1) = lat3d(sort2(1)) + ((lat3d(sort2(4)) - lat3d(sort2(1))) * verh1)
	    loncoord(2) = lon3d(sort2(1)) + ((lon3d(sort2(3)) - lon3d(sort2(1))) * verh1)
	    latcoord(2) = lat3d(sort2(1)) + ((lat3d(sort2(3)) - lat3d(sort2(1))) * verh1)
	    loncoord(4) = lon3d(sort2(2)) + ((lon3d(sort2(4)) - lon3d(sort2(2))) * verh1)
	    latcoord(4) = lat3d(sort2(2)) + ((lat3d(sort2(4)) - lat3d(sort2(2))) * verh1)
	    loncoord(3) = lon3d(sort2(2))
	    latcoord(3) = lat3d(sort2(2))
	 
        call ccw(loncoord(1), latcoord(1), loncoord(2), latcoord(2), lon, lat, ccwt(1))
	    call ccw(loncoord(2), latcoord(2), loncoord(3), latcoord(3), lon, lat, ccwt(2))
	    call ccw(loncoord(3), latcoord(3), loncoord(4), latcoord(4), lon, lat, ccwt(3))
	    call ccw(loncoord(4), latcoord(4), loncoord(1), latcoord(1), lon, lat, ccwt(4))
	 
	    c1 = 0
        c2 = 0
        do k = 1, 4
           if (ccwt(k).gt.0) then
	          c1 = c1 + 1
	       elseif (ccwt(k).lt.0) then
	          c2 = c2 + 1
	       endif
        enddo
      
        if ((c1.ne.0).and.(c2.ne.0)) then
        else
	    !print*, 'return'
	       return
	    endif
	 
	    !print*, loncoord
	    !print*, latcoord
     endif
      
     if (counter.eq.3) then ! red tetrahedra
        l = 0
	    do k = 1, 4
           if (depth3d(elemtetr(k)).eq.depthu) node = k
	    enddo
        do k = 1, 4
	       if (k.ne.node) then
	          l = l + 1
	          loncoord(l) = lon3d(elemtetr(node)) + ((lon3d(elemtetr(k)) - lon3d(elemtetr(node))) * verh2)
	          latcoord(l) = lat3d(elemtetr(node)) + ((lat3d(elemtetr(k)) - lat3d(elemtetr(node))) * verh2)
	       endif
        enddo
     endif
      
      !print*, 'hier5'
      !befindet sich der Drifter im gefundenen Dreieck?
     call ccw(loncoord(2), latcoord(2), loncoord(3), latcoord(3), lon, lat, ccwt(1))
     call ccw(loncoord(3), latcoord(3), loncoord(1), latcoord(1), lon, lat, ccwt(2))
     call ccw(loncoord(1), latcoord(1), loncoord(2), latcoord(2), lon, lat, ccwt(3))
     c1 = 0
     c2 = 0
     do k = 1, 3
        if (ccwt(k).gt.0) then
	       c1 = c1 + 1
	    elseif (ccwt(k).lt.0) then
	       c2 = c2 + 1
	    endif
     enddo
      
     !print*, 'hier6'
     if ((c1.ne.0).and.(c2.ne.0)) then
	 !print*, 'fast richtig'
        cycle
     endif
      
      !print*, c1, c2, counter
      
     return
      
  enddo
      !print*, elems(1), lon3d(elems(1)), lat3d(elems(1)), depth3d(elems(1))
      !print*, elems(2), lon3d(elems(2)), lat3d(elems(2)), depth3d(elems(2))
      !print*, elems(3), lon3d(elems(3)), lat3d(elems(3)), depth3d(elems(3))
      !print*, elems(4), lon3d(elems(4)), lat3d(elems(4)), depth3d(elems(4))
      !print*, elems(5), lon3d(elems(5)), lat3d(elems(5)), depth3d(elems(5))
      !print*, elems(6), lon3d(elems(6)), lat3d(elems(6)), depth3d(elems(6))
      !print*, lon, lat, depth
  !print*, 'wohl rausgefallen', c1, c2, counter

  find_tetra = .false.
     

end subroutine search_tetrahedron

!
!=======================================================================
!
subroutine distance(lon1,lat1,lon2,lat2,dist)

   real(kind=8)::lon1,lat1,lon2,lat2,dist
   real(kind = 8),parameter::mean_earth_radius=6371.0088

   dist=sin(lat1)*sin(lat2)+cos(lat1)*cos(lat2)*cos(lon2-lon1)
   !if(dist.gt.1.)dist=1.
   !if(dist.lt.-1.)dist=-1.
   dist=abs(acos(dist))*mean_earth_radius  ! distance in km
   return
end subroutine distance

!=======================================================================

subroutine distance_3d(lon1,lat1,lon2,lat2,depth1, depth2,dist)

   real(kind=8)::lon1,lat1,lon2,lat2,depth1,depth2,dist,a,b
   real(kind = 8),parameter::mean_earth_radius=6371.0088

   a=sin(lat1)*sin(lat2)+cos(lat1)*cos(lat2)*cos(lon2-lon1)
   !if(dist.gt.1.)dist=1.
   !if(dist.lt.-1.)dist=-1.
   a= abs(acos(a))*mean_earth_radius*1000 ! distance in m
   b=depth1 - depth2
   dist = sqrt( a*a + b*b)
   return
end subroutine distance_3d
!
!=======================================================================
!
subroutine ccw(p1x, p1y, p2x, p2y, x, y, ccwt)   !neu

   implicit none
   real(kind = 8):: p1x, p1y, p2x, p2y, x, y
   real(kind = 8):: dx1, dx2, dy1, dy2
   integer:: ccwt
   
   dx1 = p1x - x
   dy1 = p1y - y
   dx2 = p2x - x
   dy2 = p2y - y
   
   ccwt=0
   if (dx1*dy2.gt.dy1*dx2)ccwt=1
   if (dx1*dy2.lt.dy1*dx2)ccwt=-1
   
end subroutine ccw
!
!=======================================================================
!
subroutine distance_2(lon1, lat1, lon2, lat2, depth1, depth2, dist, pi)    !neu
!lon und lat in grad
   
   implicit none
   real(kind = 8):: lon1, lat1, lon2, lat2, depth1, depth2, dist
   real(kind = 8):: x1, x2, y1, y2, z1, z2, pi
   
   x1 = depth1 * sin((lat1 + 90)*pi/180.)*cos((lon1 + 180)*pi/180.)
   y1 = depth1 * sin((lat1 + 90)*pi/180.)*sin((lon1 + 180)*pi/180.)
   z1 = depth1 * cos((lat1 + 90)*pi/180.)
   x2 = depth2 * sin((lat2 + 90)*pi/180.)*cos((lon2 + 180)*pi/180.)
   y2 = depth2 * sin((lat2 + 90)*pi/180.)*sin((lon2 + 180)*pi/180.)
   z2 = depth2 * cos((lat2 + 90)*pi/180.)
   
   dist = sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2) + (z1 - z2)*(z1 - z2))
   
end subroutine distance_2

!=======================================================================

!
!=======================================================================
!


subroutine interpolate_surface_ts(lon_p, lat_p, elem, glon, glat, t_oce, s_oce, temp_out, salt_out)
  implicit none
  ! Inputs
  real(kind=8), intent(in) :: lon_p, lat_p          ! Particle longitude/latitude (radians)
  integer,     intent(in) :: elem(3)                ! 3 node indices of surface triangle
  real, intent(in) :: glon(:), glat(:)      ! Longitude/latitude of all mesh nodes
  real, intent(in) :: t_oce(:), s_oce(:)    ! Temperature and salinity at all nodes

  ! Outputs
  real(kind=8), intent(out) :: temp_out, salt_out   ! Interpolated temperature and salinity

  ! Locals
  real(kind=8) :: dist2d(3), wsurf(3), distsum2d
  integer :: k
  !print *, 'lon_p, lat_p:', lon_p, lat_p
  !print *, 'size(glon), size(glat), elem:', size(glon), size(glat), elem
  ! print size of t_oce and s_oce
  !print *, 'size(t_oce), size(s_oce):', size(t_oce), size(s_oce)
  ! print element nodes
  ! Compute 2D distances from particle to each triangle node
  do k = 1, 3
    call distance_2d(lon_p, lat_p, glon(elem(k)), glat(elem(k)), dist2d(k))
  end do

  !print *, 'dist2d:', dist2d, 'elem:', elem

  ! Handle zero-distance (exact node match)
  if (any(dist2d == 0.0)) then
    do k = 1, 3
      if (dist2d(k) == 0.0) then
        wsurf(k) = 1.0
      else
        wsurf(k) = 0.0
      end if
    end do
  else
    distsum2d = sum(1.0 / dist2d)
    do k = 1, 3
      wsurf(k) = (1.0 / dist2d(k)) / distsum2d
    end do
  end if

  ! Weighted average of temperature and salinity
  temp_out = 0.0
  salt_out = 0.0
  do k = 1, 3
    temp_out = temp_out + wsurf(k) * t_oce(elem(k))
    salt_out = salt_out + wsurf(k) * s_oce(elem(k))
    !print *, 'k:', k, 'wsurf(k):', wsurf(k), 't_oce(elem(k)):', t_oce(elem(k)), 'temp_out:', temp_out
  end do
end subroutine interpolate_surface_ts


subroutine distance_2d(lon1, lat1, lon2, lat2, dist)
  implicit none
  real(kind=8), intent(in)  :: lon1, lat1, lon2, lat2
  real(kind=8), intent(out) :: dist
  real(kind=8) :: dlon, dlat, a, c
  real(kind=8), parameter :: r = 6371000.0D0  ! Earth radius in meters

  dlon = lon2 - lon1
  dlat = lat2 - lat1
  a = sin(dlat/2.0D0)**2 + cos(lat1)*cos(lat2)*sin(dlon/2.0D0)**2
  c = 2.0D0 * atan2(sqrt(a), sqrt(1.0D0 - a))
  dist = r * c
end subroutine distance_2d

!=======================================================================
subroutine calculate_triangle_weights(lon_p, lat_p, elem, glon, glat, weights)
   implicit none

   ! Inputs
   real(kind=8), intent(in) :: lon_p, lat_p          ! Particle longitude/latitude (radians)
   integer,     intent(in) :: elem(3)                ! 3 node indices element to calculate the weights for
   real, intent(in) :: glon(:), glat(:)              ! Longitude/latitude of all mesh nodes

   ! Outputs
   real(kind=8), intent(out) :: weights(3)           ! Interpolated weights

   ! Locals
   integer :: k
   real(kind=8) :: dist2d(3), wsurf(3), distsum2d

   ! Compute interpolated U and V using basis functions
   do k = 1, 3
      call distance_2d(lon_p, lat_p, glon(elem(k)), glat(elem(k)), dist2d(k))
   end do

   ! Handle zero-distance (exact node match)
   if (any(dist2d == 0.0)) then
      do k = 1, 3
         if (dist2d(k) == 0.0) then
            weights(k) = 1.0
         else
            weights(k) = 0.0
         end if
      end do
   else
      distsum2d = sum(1.0 / dist2d)
      do k = 1, 3
         weights(k) = (1.0 / dist2d(k)) / distsum2d
      end do
  end if

end subroutine calculate_triangle_weights

!=======================================================================

subroutine calculate_triangle_weighted_mean(weights, elem, field, triangle_mean_out)
   implicit none
   ! Inputs
   real(kind=8), intent(in) :: weights(3)
   integer,      intent(in) :: elem(3)
   real(kind=8), intent(in) :: field(:)
   ! Output
   real(kind=8), intent(out) :: triangle_mean_out
   ! Locals
   integer :: k

   triangle_mean_out = 0.0
   do k = 1, 3
      triangle_mean_out = triangle_mean_out + weights(k) * field(elem(k))
   end do

end subroutine calculate_triangle_weighted_mean


!!! -------------------------------------------------------------------
          
subroutine calculate_mld(step_day, drifter, ref_level, z_layers, &
                         steps_per_day, num_drifters, &
                         btemp_column, bsalt_column, depths, &
                         phi1, lamda1, fill_val, delta_sigma_thresh, mld)

 
   ! ===================================================================
   ! PURPOSE:
   !  Calculates the Mixed Layer Depth (MLD) using a density threshold 
   !  method and linear interpolation between depth levels.
   !
   ! INPUTS:
   !  - step_day: Integer, index of the current time step/day
   !  - drifter: Integer, index of the current drifter/float
   !  - ref_level: Integer, the vertical layer index used as reference
   !  - z_layers: Integer, total number of vertical layers
   !  - steps_per_day: Integer, size of the time dimension of the buffers
   !  - num_drifters: Integer, size of the drifter dimension of the buffers
   !  - btemp_column: 3D real array (steps_per_day, num_drifters, z_layers)
   !  - bsalt_column: 3D real array (steps_per_day, num_drifters, z_layers)
   !  - depths: 1D real array (z_layers), depth values at vertical levels
   !  - phi1: Real, latitude of the drifter for pressure calculation
   !  - lamda1: Real, longitude of the drifter for absolute salinity
   !  - fill_val: Real, missing/fill data value indicator
   !  - delta_sigma_thresh: Real, potential density threshold for MLD
   !
   ! OUTPUTS:
   !  - mld: 2D real array (steps_per_day, num_drifters), updated at (step_day, drifter)
   ! ===================================================================
 

 
   use gsw_mod_toolbox, only: gsw_p_from_z, gsw_sa_from_sp, gsw_ct_from_pt, gsw_sigma0
   
   implicit none
 
   ! Arguments
   integer :: step_day, drifter, ref_level, z_layers
   integer :: steps_per_day, num_drifters
   real :: btemp_column(steps_per_day, num_drifters, z_layers)
   real :: bsalt_column(steps_per_day, num_drifters, z_layers)
   real(8) :: depths(z_layers)
   real(8) :: phi1, lamda1, fill_val, delta_sigma_thresh
   real :: mld(steps_per_day, num_drifters)

 
   ! Local Variables
   integer :: k
   real(8) :: T_ref, S_ref, p_ref, SA_ref, CT_ref, sigma0_ref
   real(8) :: T_k, S_k, p_k, SA_k, CT_k, sigma0_k
   real(8) :: p_km1, SA_km1, CT_km1, sigma0_km1
 
   ! Extract reference values
   T_ref = btemp_column(step_day, drifter, ref_level)
   S_ref = bsalt_column(step_day, drifter, ref_level)
 
   if (T_ref > fill_val + 1.0d0 .and. S_ref > fill_val + 1.0d0) then
 
      p_ref  = gsw_p_from_z(-depths(ref_level), phi1)
      SA_ref = gsw_sa_from_sp(S_ref, p_ref, lamda1, phi1)
      CT_ref = gsw_ct_from_pt(SA_ref, T_ref)
      sigma0_ref = gsw_sigma0(SA_ref, CT_ref)
 
      ! Set fallback to the reference depth right away
      mld(step_day, drifter) = depths(ref_level)
 
      do k = ref_level + 1, z_layers
         T_k = btemp_column(step_day, drifter, k)
         S_k = bsalt_column(step_day, drifter, k)
 
         if (T_k < fill_val + 1.0d0 .or. S_k < fill_val + 1.0d0) exit
 
         ! Guard T/S ranges before each GSW call (polar oceans adjustment)
         if (T_k < -2.5d0 .or. T_k > 25.d0 .or. S_k < 0.d0 .or. S_k > 42.d0) then
            exit  ! Treat as bottom of valid data
         end if
 
         p_k    = gsw_p_from_z(-depths(k), phi1)
         SA_k   = gsw_sa_from_sp(S_k, p_k, lamda1, phi1)
         CT_k   = gsw_ct_from_pt(SA_k, T_k)
         sigma0_k = gsw_sigma0(SA_k, CT_k)
 
         if ((sigma0_k - sigma0_ref) >= delta_sigma_thresh) then
            ! Linear interpolation between level k-1 and k
            p_km1    = gsw_p_from_z(-depths(k-1), phi1)
            SA_km1   = gsw_sa_from_sp(bsalt_column(step_day,drifter,k-1), p_km1, &
                                      lamda1, phi1)
            CT_km1   = gsw_ct_from_pt(SA_km1, btemp_column(step_day,drifter,k-1))
            sigma0_km1 = gsw_sigma0(SA_km1, CT_km1)
 
            mld(step_day, drifter) = depths(k-1) + &
               (depths(k) - depths(k-1)) * &
               (delta_sigma_thresh - (sigma0_km1 - sigma0_ref)) / &
               (sigma0_k - sigma0_km1)
            
            exit
         end if
      end do
 
   end if   
 
end subroutine calculate_mld


!!! -------------------------------------------------------------------

subroutine load_basis_functions(filename, bafux, bafuy, el2d, voltriangle)
  implicit none

#include "netcdf.inc"

  character(len=*), intent(in) :: filename
  real(kind=8), allocatable, intent(out) :: bafux(:,:), bafuy(:,:), voltriangle(:)
  integer, intent(out) :: el2d

  integer :: status, ncid
  integer :: varid_bafux, varid_bafuy, varid_voltriangle
  integer :: dimid_el2d, dimlen_el2d
  integer :: istart(2), icount(2)
  real(kind=8), allocatable :: aux(:)
  integer :: ndims, dimids(10), i
  integer :: dimlen(10)
  character(len=NF_MAX_NAME) :: dimname
  real(kind=8), allocatable :: bfx_tmp(:,:), bfy_tmp(:,:)
  integer :: xtype, natts, d0, d1
  character(len=NF_MAX_NAME) :: vname, dname0, dname1

  print *, 'load_basis_functions started'
  print *, 'filename=', filename
  !print *, 'bafux allocated? ', allocated(bafux)
  !print *, 'bafuy allocated? ', allocated(bafuy)

  ! Open NetCDF file read-only
  status = nf_open(filename, nf_nowrite, ncid)
  if (status /= nf_noerr) call handle_err(status)
  
  ! inquire dim order of bafux_2d
  !status = nf_inq_var(ncid, varid_bafux, vname, xtype, ndims, dimids, natts)
  !if (status /= nf_noerr) call handle_err(status)

  !status = nf_inq_dimlen(ncid, dimids(1), d0)
  !if (status /= nf_noerr) call handle_err(status)
  !status = nf_inq_dimlen(ncid, dimids(2), d1)
  !if (status /= nf_noerr) call handle_err(status)

  !status = nf_inq_dimname(ncid, dimids(1), dname0)
  !if (status /= nf_noerr) call handle_err(status)
  !status = nf_inq_dimname(ncid, dimids(2), dname1)
  !if (status /= nf_noerr) call handle_err(status)
  !print*, 'bafux dims = (', trim(dname0), '=', d0, ',', trim(dname1), '=', d1, ')'


  ! Get dimension ID and length for el2d
  status = nf_inq_dimid(ncid, 'el2d', dimid_el2d)
  if (status /= nf_noerr) call handle_err(status)

  status = nf_inq_dimlen(ncid, dimid_el2d, dimlen_el2d)
  if (status /= nf_noerr) call handle_err(status)

  el2d = dimlen_el2d
  print *, 'el2d=', el2d

  ! Get variable IDs
  status = nf_inq_varid(ncid, 'bafux_2d', varid_bafux)
  if (status /= nf_noerr) call handle_err(status)

  status = nf_inq_varid(ncid, 'bafuy_2d', varid_bafuy)
  if (status /= nf_noerr) call handle_err(status)

  status = nf_inq_varid(ncid, 'voltriangle', varid_voltriangle)
  if (status /= nf_noerr) call handle_err(status)

  !print *, 'bafux_2d variable has ', ndims, ' dimensions:'
  do i = 1, ndims
   status = nf_inq_dimlen(ncid, dimids(i), dimlen(i))
   status = nf_inq_dimname(ncid, dimids(i), dimname)
   !print *, 'Dim ', i, ': ', trim(dimname), ' size=', dimlen(i)
  end do

  ! Allocate output arrays now that we know size
  !allocate(bafux(3, el2d))!!!!!!!!!!!!!!!
  !allocate(bafuy(3, el2d))!!!!!!!!!!!!!!!
  allocate(voltriangle(el2d))


  !allocate(aux(3*el2d))
  allocate(bfx_tmp(el2d,3), bfy_tmp(el2d,3))   ! file likely stores (el2d,3)

  status = nf_get_var_double(ncid, varid_bafux, bfx_tmp)
  if (status /= nf_noerr) call handle_err(status)
  status = nf_get_var_double(ncid, varid_bafuy, bfy_tmp)
  if (status /= nf_noerr) call handle_err(status)

  bafux = transpose(bfx_tmp)   ! now bafux is (3, el2d)
  bafuy = transpose(bfy_tmp)

  !print*, 'bafu 11, 21, 31', bafux(1,1), bafux(2,1), bafux(3,1)
  !print*, 'bafu 11, 12, 13', bafux(1,1), bafux(1,2), bafux(1,3)
  !print*, 'bafu 11, 21, 31', bafuy(1,1), bafuy(2,1), bafuy(3,1)
  !print*, 'bafu 11, 12, 13', bafuy(1,1), bafuy(1,2), bafuy(1,3)


  deallocate(bfx_tmp, bfy_tmp)

  ! Read bafux
  !istart = (/1,1/)!!!!!!!!!!!!
  !icount = (/3, el2d/)!!!!!!!!!!!
  !status = nf_get_vara_double(ncid, varid_bafux, istart, icount, aux)
  !status = nf_get_var_double(ncid, varid_bafux, bafux)!!!!!!!!!!
  !if (status /= nf_noerr) call handle_err(status)!!!!!!!!!!!!
  !bafux = reshape(aux, shape(bafux))

  ! Read bafuy
  !status = nf_get_vara_double(ncid, varid_bafuy, istart, icount, aux)
  !status = nf_get_var_double(ncid, varid_bafuy, bafuy)!!!!!!!!
  !if (status /= nf_noerr) call handle_err(status)!!!!!!!!!!
  !bafuy = reshape(aux, shape(bafuy))

  status = nf_get_var_double(ncid, varid_voltriangle, voltriangle)
  if (status /= nf_noerr) call handle_err(status)

  ! Close NetCDF file
  status = nf_close(ncid)
  if (status /= nf_noerr) call handle_err(status)
  

  !deallocate(aux)

end subroutine load_basis_functions



