program trajectories
!von Katja Rollenhagen / Bettina Fach / Ralph Timmermann / Claudia Wekerle (AWI)


use g_config
use g_rotate_grid
use o_param
use o_mesh
use o_elements
use g_read_mesh
use g_mod_output
use o_DATA_TYPES

implicit none

 integer                                :: i, ii, j, jj, kk, ll, drifter_id, cnt, drifter
 integer                                :: id,elem(3), l, runaway, clayer, elemtetr(4)
 integer                                :: index,num, nod, years, d, d_prev, y, cyear, syear, eyear, elems(6), el_index
 real, allocatable, dimension(:)        :: time, blat, blon, bdepth, btemp, bsalt
 real, allocatable, dimension(:)        :: btemp_surface, bsalt_surface, bvort

 integer,allocatable,dimension(:)       :: bday, byear
 real(kind=8),allocatable,dimension(:)  :: dlon,dlat,ddepth,dw_oce,dfirst_day
 real(kind=8)                           :: phi1,phi2,lamda1,lamda2,z1,z2,dt,re,w(4)
 real(kind=8)                           :: umean,vmean,wmean, tmean, smean
 real(kind=8)                           :: umean_prev,vmean_prev,wmean_prev, tmean_prev, smean_prev

 real(kind=8)                           :: cosd,distsum,dist(4),lat1,lam1, distd(2), deptho, depthu
 integer                                :: yearnew, ndoyrnew
 real(kind=8), allocatable,dimension(:) :: u_oce,v_oce,w_oce, t_oce, s_oce, sinking_oce
 real(kind=8), allocatable,dimension(:) :: u_oce_prev,v_oce_prev,w_oce_prev, t_oce_prev, s_oce_prev ! data of the previous day

 real(kind=8), allocatable              :: bafux(:,:), bafuy(:,:), voltriangle(:)
 integer                                :: el2d
 real(kind=8)                           :: dv_dx, du_dy, vort, f

 integer                                :: nodes_top(3), nodes_bot(3)
 real(kind=8)                           :: u_nodes(3), v_nodes(3)
 real(kind=8)                           :: alpha
 real(kind=8)                           :: dv_dx_surf, du_dy_surf, vort_surf

 integer                                :: num_days,num_days_all,tt, nrec, steps_per_day, first_day, end_day, step

 real                                   :: start_time, stop_time, start_time1, stop_time1, timenew
 real(kind=8)                           :: tol

 character                              :: cyearchar*4, nyearchar*4, syearchar*4, eyearchar*4
 character                              :: first_daychar*3, num_dayschar*2, drifter_id_char*3
 character(250)                         :: filename, filename_prev, MeshPath, ResultPath, runid
 character(200)                         :: OutputFile, OutputPath, InputFile, InputPath, BafuFile, suffix, outfile, infile

 logical                                :: backwards, sinkingrate, find_tetra, put_back, twoD_tracking, inside
 integer                                :: e
 real(kind=8)                           :: maxsumx, maxsumy, weight

 integer :: el, counter, m
 integer, allocatable :: nmb_el_neigh(:)
 integer, allocatable :: el_neigh_addresses(:,:)
 integer :: unit
 integer :: found_el
 integer :: found_el_nodes(3)

 integer :: ios

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


!!! -------------------------------------------------------------------------------------
!!! CONFIG 

! --- Namelist blocks ---
namelist /Paths/ MeshPath, ResultPath, InputPath, OutputPath, BafuFile, outfile, infile, runid
namelist /Parameters/ syear, num_days, steps_per_day, twoD_tracking, backwards, put_back
  
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
print *, "twoD_tracking = ", twoD_tracking
print *, "backwards     = ", backwards
print *, "put_back      = ", put_back


 eyear = syear+1
 write(syearchar, '(i4)')syear
 years = eyear-syear+1
 num_days_all = num_days

 InputFile = trim(InputPath)//trim(infile)
 OutputFile = trim(OutputPath)//trim(outfile)

!!! -------------------------------------------------------------------------------------
!!! READ MESH 
print*, 'rotated_grid: ',rotated_grid
 if(rotated_grid) call calculate_rotate_matrix

 ! call calculate_rotate_matrix
 call load_basis_functions(BafuFile, bafux, bafuy, el2d, voltriangle)
 print*,'bafux: ', shape(bafux), 'bafuy: ', shape(bafuy), 'voltriangle: ', shape(voltriangle)
! quick spot checks
!e = 1
!print*, 'e=', e, sum(bafux(:,e)), sum(bafuy(:,e))

!e = 2
!print*, 'e=', e, sum(bafux(:,e)), sum(bafuy(:,e))

! optional: check the worst case over ALL elements
maxsumx = 0d0; maxsumy = 0d0
do e = 1, el2d
  maxsumx = max(maxsumx, abs(sum(bafux(:,e))))
  maxsumy = max(maxsumy, abs(sum(bafuy(:,e))))
end do
!write(*,'(A,ES12.4)') 'max |sum bafux(:,e)| = ', maxsumx
!write(*,'(A,ES12.4)') 'max |sum bafuy(:,e)| = ', maxsumy

call read_2Dmesh(MeshPath)
call read_3Dmesh(MeshPath)

print*,'elem2D_nodes',elem2D_nodes(:,:2)
print*,"elem2d",elem2d
call mesh_scaling   ! long., lat. are transf. into rad
call find_layer_elem3d
call find_tetra_in_prism
 
write(*,*) 'reading mesh: DONE'

!!! ------------------------------------------------------------------------------------

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

print*,"neighbour od element 519505:", el_neigh_addresses(:,519505), "so", nmb_el_neigh(519505)


!!! -------------------------------------------------------------------------------------
!!! READ DRIFTER INPUT FILE


 open(10,file=InputFile)
 read(10,*)num
 allocate(dlon(num),dlat(num),ddepth(num),dw_oce(num),dfirst_day(num))
 do i=1,num 
    read(10,*) dlon(i),dlat(i),ddepth(i),dw_oce(i),dfirst_day(i)
 enddo
 close(10)
 print*,'finished reading input file... number of drifters:', num
 first_day=dfirst_day(1)
 print*,'first_day = ', first_day, 'vel: ',dw_oce(1)

!!! -------------------------------------------------------------------------------------
!!! INITIALIZE OUTPUT FILE
 write(syearchar, '(i4)')syear
 write(eyearchar, '(i4)')eyear
 write(first_daychar, '(i2)')first_day
 write(num_dayschar, '(i2)')num_days
 drifter_id=999
 write(drifter_id_char, '( i3.3)') drifter_id
 call init_output(OutputFile,num,years*num_days*steps_per_day,syear)


!!! -------------------------------------------------------------------------------------
!!! ALLOCATE ARRAYS
 allocate(u_oce(nod3d),v_oce(nod3d),w_oce(nod3d),sinking_oce(nod3d),t_oce(nod3d),s_oce(nod3d))
 allocate(u_oce_prev(nod3d),v_oce_prev(nod3d),w_oce_prev(nod3d),t_oce_prev(nod3d),s_oce_prev(nod3d))

 allocate(time(years*num_days*steps_per_day))
 allocate(bday(years*num_days*steps_per_day))
 allocate(byear(years*num_days*steps_per_day))
 allocate(blon(years*num_days*steps_per_day))
 allocate(blat(years*num_days*steps_per_day))
 allocate(bdepth(years*num_days*steps_per_day))
 allocate(btemp(years*num_days*steps_per_day))
 allocate(bsalt(years*num_days*steps_per_day))
 allocate(btemp_surface(years*num_days*steps_per_day))
 allocate(bsalt_surface(years*num_days*steps_per_day))
 allocate(bvort(years*num_days*steps_per_day))




!!! -------------------------------------------------------------------------------------
!!! LOOP OVER DRIFTERS
do drifter = 1, num 

   print*, 'Start calculating drifter ', drifter

   ! read new initial drifter position 
   lamda1=dlon(drifter)
   if(lamda1.gt.180.)lamda1=lamda1-360.
   if(lamda1.lt.-180.)lamda1=lamda1+360.
   phi1=dlat(drifter)
   z1=ddepth(drifter) 
 !call cpu_time(start_time)

   time   = 0.
   bday   = 0
   byear  = 0
   blon   = 0.
   blat   = 0.
   bdepth = 0.
   btemp  = 0.
   bsalt  = 0.

   runaway = 0
   
   j=0

!!! -------------------------------------------------------------------------------------
!!! LOOP OVER YEARS
   do y = 1, years

      if (runaway == 1) then  
         print *,'runaway... exit program'
         exit
      end if 


      if ( j>0 )  then       ! second loop over years
         ! exit when it's not the first of the year, !!!or we're already at the surface
         if (backwards) then ! backwards: we start with 31 Dec
            if ( bday(j) > 1 .OR. bdepth(j) == 0.) then
            !if ( bday(j) > 1 ) then
               exit
            else  ! start the next year
               num_days=num_days-first_day
               first_day=365
            end if
         else ! forwards: we start with 1 Jan
             if ( bday(j) < 365 ) then
               exit
            else  ! start the next year
               ! recompute the number of days that still should be computed
               num_days=num_days-(365-first_day+1) 
               first_day=1
            end if
         end if

      else if ( j .EQ. 0) then
         first_day=dfirst_day(drifter)
         num_days=num_days_all
      end if


     ! NEW FESOM OUTPUT FILE
      if (backwards) then
         cyear = syear-(y-1)
      else
         cyear = syear+(y-1)
      end if
      write(cyearchar, '(i4)')cyear
      !filename=trim(ResultPath)//trim(runid)//'.'//trim(cyearchar)//trim('.oce.mean.sub.Mar_Jul_top20m.nc')
      filename=trim(ResultPath)//trim(runid)//'.'//trim(cyearchar)//trim('.oce.mean.sub.nc')
      !filename = "/albedo/work/user/quraus001/FESOM_particles/preparation/dummy_data.nc"  ! <--------------------------- FOR TESTING WITH DUMMY DATA
      !print*,  'drifter ID: ', drifter_id, 'filename: ',filename
      !print*, 'filename = [', trim(filename), ']'
      !print*, 'len_trim(filename) =', len_trim(filename)
      !print*, 'char codes:'
      !do i=1, len(filename)
      !   print*, i, iachar(filename(i:i))
      !end do

!!! -------------------------------------------------------------------------------------
!!! LOOP OVER DAYS
      if (backwards) then
         end_day   = first_day - num_days + 1
         step      = -1
      else
         end_day   = first_day + num_days - 1
         step      = 1
      end if

      !do d=first_day, first_day-num_days+1,-1 ! backwards
      !do d=first_day, first_day+num_days-1,1 ! forward time loop  
      do d = first_day, end_day, step

         print*,'Load data'
         !call cpu_time(start_time)
        
         print*, "calculating day",d

         if (backwards) then
            d_prev = d + 1
         else
            d_prev = d - 1
         end if

         ! make sure d_prev is within valid range [1, 365] and adjust filename_prev accordingly
         if (d_prev < 1) then ! set to last day of previous year
            d_prev = 365
            ! adjust filename_prev
            write(nyearchar, '(i4)') cyear - 1
            filename_prev = trim(ResultPath)//trim(runid)//'.'//trim(nyearchar)//trim('.oce.mean.sub.nc')
         else if (d_prev > 365) then ! set to first day of next year
            d_prev = 1
            ! adjust filename_prev
            write(nyearchar, '(i4)') cyear + 1
            filename_prev = trim(ResultPath)//trim(runid)//'.'//trim(nyearchar)//trim('.oce.mean.sub.nc')
         else
            filename_prev = filename
         end if


         ! If we are on the first day, read the data for that day and the next day
         if (d == first_day) then
            call oce_input_netcdf(filename,d,u_oce,v_oce,w_oce,t_oce,s_oce,nod2d,nod3d)  ! <--------------------------- CHANGED TO FIXED FIRST DAY
            call oce_input_netcdf(filename_prev,d_prev,u_oce_prev,v_oce_prev,w_oce_prev,t_oce_prev,s_oce_prev,nod2d,nod3d)
         else
            ! copy data from last current day to previous day
            u_oce_prev = u_oce
            v_oce_prev = v_oce
            w_oce_prev = w_oce
            t_oce_prev = t_oce
            s_oce_prev = s_oce

            call oce_input_netcdf(filename,d,u_oce,v_oce,w_oce,t_oce,s_oce,nod2d,nod3d)

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
         !w_oce=-dw_oce(drifter) ! constant sinking speed
         !w_oce=0.0                                                                   
         
         if (twoD_tracking) then
            w_oce = 0.0
         !else
            !w_oce=-dw_oce(drifter)   ! sinking
            ! or: w_oce already read from oce_input_netcdf
         end if

         !call cpu_time(stop_time)
         !print *, "****** Time ", &
         !   stop_time - start_time, "seconds --- call oce_input_netcdf"

!!! -------------------------------------------------------------------------------------
!!! LOOP OVER STEPS PER DAY
         do tt=1, steps_per_day 

            ! calculate weight between previous day and current day
            weight = real(tt-1) / real(steps_per_day-1)

            ndoyrnew=d-1; 
            timenew=d-1; 
            
            j=j+1 
            
            dt=86400/steps_per_day

            if (backwards) then 
               time(j)= d*86400 - (tt-1)* dt 
            else 
               time(j)= (d-1)*86400 + (tt-1)* dt
            end if

            ! radius earth (linear combination from equator and polar radius) 
            re=phi1/90.*rpol+(90.-phi1)/90.*req

            lam1=pi/180.*lamda1
            lat1=pi/180.*phi1
            !   print*,lam1,lat1, dt, re
	 
            !!! -------------------------------------------------------------------
	         !   determine upper layer (clayer)
	         do i = 1, z_layers
	            if (z1.ge.depths(i)) clayer = i
	         enddo
	 
            !!! -------------------------------------------------------------------
            ! upper and lower layer depth
	         deptho = depths(clayer)
	         depthu = depths(clayer + 1)

            !!! -------------------------------------------------------------------
            ! FIND 2D POSITION OF THE PARTICLE IN THE FESOM MESH

            if (j .GT. 1) then ! If we are not in the first computing step start local search for current element
               ! First check if we are still in the old element
               call check_old_element(pi/180.*lamda1, pi/180.*phi1, elem, el_index, glon, glat, nod2d, inside)

               if (.NOT. inside) then ! If we are not in the same element anymore we probably moved into a neighbour
                  print*,"Particle moved outside the old element, Checking neighbors..." 
                  call check_neighbour_elements(pi/180.*lamda1, pi/180.*phi1, el_index, elem2d_nodes, nmb_el_neigh(el_index), el_neigh_addresses(1:nmb_el_neigh(el_index),el_index), glon, glat, nod2d, elem2d, found_el, found_el_nodes, inside)

                  if (.NOT. inside) then ! If we are also not inside any of the neighbours we have to fall back to global search
                     print*,'!!! Couldnt find particle in the neighbor elements'
                     !print*,"el_index",el_index
                  else ! If the particle is inside one of the neighbour elements we can update our element index and nodes
                        print*,'Found in neighbor element' 
                        print*, "Updating new element index from", el_index , "to", found_el
                        el_index = found_el
                        elem = found_el_nodes
                  endif

               else ! If we remain in the old element we dont need to update any variables
                  print*, "Particle remains in old element" 
               endif

            endif


            if ((j .EQ. 1) .OR. .NOT. inside) then ! If we are in the first computing step or we didnt find the particle anywhere -> do global search
            
               print*,"!!! Start global search for matching element"
               call nearest_point(lam1,lat1,index,nod2d,glon,glat)

               print*,'neareast point for',lamda1,phi1,'is',index,glon(index)*180/pi,glat(index)*180./pi
            

               call search_element(pi/180.*lamda1, pi/180.*phi1, elem,el_index,index,nod2d,elem2d,glon,glat,elem2d_nodes)
               print*,"Global search found particle in 2D element with index", el_index
            
            end if

            !!!! -------------------------------------------------------------------
            ! calculate vorticity
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
               !t_nodes(l) = (1.0d0 - alpha) * t_oce(nodes_top(l)) + alpha * t_oce(nodes_bot(l))
               !s_nodes(l) = (1.0d0 - alpha) * s_oce(nodes_top(l)) + alpha * s_oce(nodes_bot(l))
            end do

            ! --- horizontal derivatives and relative vorticity (correct sign)
            dv_dx = sum( v_nodes * bafux(:, el_index) )
            du_dy = sum( u_nodes * bafuy(:, el_index) )
            vort  = dv_dx - du_dy
            bvort(j) = vort  
            ! --- calculate more derivatives (du/dx, dv/dy, dw/dx, dw/dy, db/dx, db/dy)
            !du_dx = sum( u_nodes * bafux(:, el_index) )
            !dv_dy = sum( v_nodes * bafuy(:, el_index) )
            !dw_dx = sum( w_nodes * bafux(:, el_index) )
            !dw_dy = sum( w_nodes * bafuy(:, el_index) )
            !db_dx = sum( b_nodes * bafux(:, el_index) )
            !db_dy = sum( b_nodes * bafuy(:, el_index) )

            call interpolate_surface_ts(lam1, lat1, elem, glon, glat, t_oce, s_oce, btemp_surface(j), bsalt_surface(j))


            !call get_surface_TS(lamda1, phi1, elem, nod3D_below_nod2D, t_oce, s_oce, lon3d, lat3d, btemp_surface(j), bsalt_surface(j))
            

            !!! -------------------------------------------------------------------
            !   particle outside model domain (horizontally)
            if (sum(elem).lt.6) then
               print*,'reached the boundary.... '
               if (put_back .eq. .false.) then
                  print*,'.... program stops '
                  do l = j, years*num_days*steps_per_day
	                 blon(l) = blon(j - 1)
	                 blat(l) = blat(j - 1)
	                 bdepth(l)=bdepth(j-1)
	                 time(l) = (j - 1) 
	              end do
	              runaway = 1
                  exit
               elseif  (put_back .eq. .true.) then ! 
               ! set point to center of element
               ! needs to be tested some more
               ! does not work well, particle gets stuck in a corner
                  call nearest_elem_level(lam1,lat1,elem,nod2d,elem2d,z_layers,glon,glat,elem2D_nodes,nod3D_below_nod2D,clayer + 1)
                  lamda1=( glon(elem(1)) + glon(elem(2))+glon(elem(3)) )/3.*180/pi
                  phi1=( glat(elem(1)) + glat(elem(2)) + glat(elem(3)) )/3.*180/pi
                  lam1=pi/180.*lamda1
                  lat1=pi/180.*phi1
                  print *, 'nearest point is:',lamda1,phi1
               endif
            endif

            !!! -------------------------------------------------------------------
            !   the two elements that define the prism
	        do l = 1, 3
	           elems(l) = nod3d_below_nod2d(clayer, elem(l))
	           elems(l + 3) = nod3d_below_nod2d(clayer + 1, elem(l))
	        enddo

            !!! -------------------------------------------------------------------
            !   particle outside model domain (vertically)
	        if (elems(4).eq.-999 .or. elems(5).eq.-999 .or. elems(6).eq.-999) then
               if (put_back .eq. .false.) then
	              print*,'reached the bottom.... program stops'
                  do l = j, years*num_days*steps_per_day
	                 blon(l) = blon(j - 1)
	                 blat(l) = blat(j - 1)
	                 bdepth(l)=bdepth(j - 1)
                     time(l) = time(j - 1)
	              end do
	              runaway = 1
                  exit 
               elseif  (put_back .eq. .true.) then ! 
               ! set point to center of element
               ! needs to be tested some more
               ! does not work well, particle gets stuck in a corner
                  call nearest_elem_level(lam1,lat1,elem,nod2d,elem2d,z_layers,glon,glat,elem2D_nodes,nod3D_below_nod2D,clayer + 1)
                  lamda1=( glon(elem(1)) + glon(elem(2))+glon(elem(3)) )/3.*180/pi
                  phi1=( glat(elem(1)) + glat(elem(2)) + glat(elem(3)) )/3.*180/pi
                  lam1=pi/180.*lamda1
                  lat1=pi/180.*phi1
                  print *, 'nearest point is:',lamda1,phi1
                  do l = 1, 3
	                 elems(l) = nod3d_below_nod2d(clayer, elem(l))
	                 elems(l + 3) = nod3d_below_nod2d(clayer + 1, elem(l))
	              enddo

               endif
	        endif

     
            !!! -------------------------------------------------------------------
            !   searching the tetrahedra in which 3d point lies
            !print *, "now computing new tetra..."
            !call cpu_time(start_time1)
            call search_tetrahedron(lamda1, phi1, z1, deptho, depthu, el_index, clayer,elem3d_nodes,elem2D,z_layers,elem3d, nod3d, elemtetr, depth3d, lon3d, lat3d, tetra_in_prism,find_tetra)
            !print*,"new tetragedra: ",elemtetr
            !call cpu_time(stop_time1)
            !print *, "****** Time ", &
            !      stop_time1 - start_time1, "seconds --- call search_tetrahedron"




            !!! -------------------------------------------------------------------
            !   in case tetrahedra not found, we repeat the procedure
            cnt=0 
            do while (find_tetra .eq. .false. .and. cnt .le.100)  
               cnt=cnt+1
               print*, " tetrahedron not found #####################"
            !   ! 
               tol = 0.0001 ! 10 m tolerance 
               lamda1=lamda1 + tol ! change lamda1 a little bit
               lam1=pi/180.*lamda1

               ! Check first if the particle is still in the old element
               call check_old_element(pi/180.*lamda1, pi/180.*phi1, elem, el_index, glon, glat, nod2d, inside)

               if (.NOT. inside) then ! If we are not in the same element anymore we probably moved into a neighbour
                  print*,"Tetrahedron tolerance moved particle to another element, Checking neighbors..." 
                  call check_neighbour_elements(pi/180.*lamda1, pi/180.*phi1, el_index, elem2d_nodes, nmb_el_neigh(el_index), el_neigh_addresses(1:nmb_el_neigh(el_index),el_index), glon, glat, nod2d, elem2d, found_el, found_el_nodes, inside)

                  if (.NOT. inside) then ! If we are also not inside any of the neighbours we have to fall back to global search
                     print*,'!!! Couldnt find particle in the neighbor elements, Searching globally'
                     call nearest_point(lam1,lat1,index,nod2d,glon,glat)
                     call search_element(pi/180.*lamda1,pi/180.*phi1,elem,el_index,index,nod2d,elem2d,glon,glat,elem2d_nodes)
                  else ! If the particle is inside one of the neighbour elements we can update our element index and nodes
                        print*,'Found in neighbor element' 
                        print*, "Updating new element index from", el_index , "to", found_el
                        el_index = found_el
                        elem = found_el_nodes
                  endif
               endif

               !else ! If we remain in the old element we dont need to update any variables
               !   print*, "Particle remains in old element" 
               !endif

               


               !if (.NOT. inside) then
               !   print*,'tetrahedron tolerance moved particle to another element'
               !   call nearest_point(lam1,lat1,index,nod2d,glon,glat)
               !   call search_element(pi/180.*lamda1,pi/180.*phi1,elem,el_index,index,nod2d,elem2d,glon,glat,elem2d_nodes)
               !endif

               do l = 1, 3
	              elems(l) = nod3d_below_nod2d(clayer, elem(l))
	              elems(l + 3) = nod3d_below_nod2d(clayer + 1, elem(l))
	            enddo
               call search_tetrahedron(lamda1, phi1, z1, deptho, depthu, el_index, clayer,elem3d_nodes,elem2D,z_layers,elem3d, nod3d, elemtetr, depth3d, lon3d, lat3d, tetra_in_prism,find_tetra)

            end do

	        !print*, elemtetr
            !do kk=1,4
	        !   print*, lon3d(elemtetr(kk)), lat3d(elemtetr(kk)), depth3d(elemtetr(kk))
            !end do

            !!! -------------------------------------------------------------------
            !   calculate distance of point to all four nodes that define the tetrahedra
            do kk = 1, 4
               call distance_3d(lamda1*pi/180., phi1*pi/180., lon3d(elemtetr(kk))*pi/180., lat3d(elemtetr(kk))*pi/180., z1, depth3d(elemtetr(kk)), dist(kk))
	         enddo
  
            !!! -------------------------------------------------------------------
            !   computing weights
            !print*,'**** dist: ',lamda1, lon3d(elemtetr(1)), phi1, lat3d(elemtetr(1)), z1, depth3d(elemtetr(1)), dist
            w=0.
	         do kk=1,4
               if(dist(kk)==0.)then
                  w(kk)=1.
                  goto 222
               endif
            end do

            distsum=dist(2)/dist(1)+dist(3)/dist(1)+dist(4)/dist(1)+dist(1)/dist(2)+dist(3)/dist(2)+dist(4)/dist(2)+&
            dist(1)/dist(3)+dist(2)/dist(3)+dist(4)/dist(3)+dist(1)/dist(4)+dist(2)/dist(4)+dist(3)/dist(4)
            !   print*,'distsum',distsum
            w(1)=(dist(2)/dist(1)+dist(3)/dist(1)+dist(4)/dist(1))/distsum
            w(2)=(dist(1)/dist(2)+dist(3)/dist(2)+dist(4)/dist(2))/distsum
            w(3)=(dist(1)/dist(3)+dist(2)/dist(3)+dist(4)/dist(3))/distsum
	         w(4)=(dist(1)/dist(4)+dist(2)/dist(4)+dist(3)/dist(4))/distsum
	 
	 
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
            umean = (1.0 - weight) * umean_prev + weight * umean
            vmean = (1.0 - weight) * vmean_prev + weight * vmean
            wmean = (1.0 - weight) * wmean_prev + weight * wmean
            tmean = (1.0 - weight) * tmean_prev + weight * tmean
            smean = (1.0 - weight) * smean_prev + weight * smean

            if (j .eq. 1 ) then
               blon(1)   = lamda1
               blat(1)   = phi1
	            bdepth(1) = z1
               btemp(1)  = tmean
               bsalt(1)  = smean
               bday(1)   = d
               byear(1)  = syear-(y-1)
               time(1)   = 0
               cycle
            end if
	 
            !print*,umean,vmean,wmean
            !print*,'old position:',lamda1,phi1,z1

            !!! -------------------------------------------------------------------
            !   computing new position (with Euler)
            print*,'compute new position'
            if (j.eq.2) then
               call cpu_time(start_time1)
            end if
            !call cpu_time(start_time1)

            ! new depth: 
	         z2=z1-wmean*dt
	         if (z2.lt.0) z2 = 0
  
            ! new latitude: 
            phi2=vmean*dt/re/pi*180.+phi1
   
            ! new longitude:    
            cosd=(abs(umean)*dt/re)
            cosd=cos(cosd)
            lamda2=cosd-sin(pi/180.*phi1)*sin(pi/180.*phi1)
            !   print*,j,lamda2
            lamda2=lamda2/cos(pi/180.*phi1)/cos(pi/180.*phi1)
            !   print*,j,lamda2!,acos(lamda2),acos(lamda2)*180./pi,lamda1
            if(lamda2.gt.1.)lamda2=1.
            if(umean.lt.0.)then
               lamda2=-acos(lamda2)*180./pi+lamda1
            else
               lamda2=acos(lamda2)*180./pi+lamda1
            end if
            if(lamda2.gt.180.)lamda2=lamda2-360.
            if(lamda2.lt.-180.)lamda2=lamda2+360.
            !   print*,j,lamda2
            !print*,'new position:',lamda2,phi2,z2
      
            blon(j)=lamda2
            blat(j)=phi2
	         bdepth(j)=z2
            bday(j)= d
	         byear(j)=syear-(y-1)
   
            !print*,lamda2,phi2,z2

            if (j.gt.2) then
               btemp(j-1)=tmean
               bsalt(j-1)=smean
            end if

            z1=z2
            phi1=phi2
            lamda1=lamda2
            
            if (j.eq.2) then
               call cpu_time(stop_time1)
               print *, "****** Time ", &
               stop_time1 - start_time1, "seconds --- compute new position"
            end if
            !call cpu_time(stop_time1)
            !print *, "****** Time ", &
            ! stop_time1 - start_time1, "seconds --- compute new position"


            if (runaway == 1) exit
         

         end do  ! loop over hours
         
         if (runaway == 1) then  ! particle outside model domain
            call write_output(OutputFile,drifter,years*num_days_all*steps_per_day,time,blon,blat,bdepth,btemp,bsalt,bday,byear, btemp_surface, bsalt_surface, bvort)
            exit
         end if   

         ! only stop for 3D tracking
         if (twoD_tracking == 0) then
            if (bdepth(j) == 0.) then     ! stop calculating when reaching the surface 
               call write_output(OutputFile,drifter,years*num_days_all*steps_per_day,time,blon,blat,bdepth,btemp,bsalt,bday,byear, btemp_surface, bsalt_surface, bvort)
               exit
            end if
         end if

         if (bday(j) ==365) then ! read next year when reaching Jan 1st
            exit
         end if  

         !call cpu_time(stop_time)
         !print *, "*** One day takes ", &
         !   stop_time - start_time, "seconds ------------------------"   

      end do  ! loop over days


   end do  ! loop over years
  

   !!! -------------------------------------------------------------------
   !   write output
   print *, "---- finished... writing output "
   call write_output(OutputFile,drifter,years*num_days_all*steps_per_day,time,blon,blat,bdepth,btemp,bsalt,bday,byear, btemp_surface, bsalt_surface, bvort)


 end do  ! loop over drifter drifter=1,num



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
  print*, "Searching ",elem_nmb, "neighbour elements of element", current_el
  !print*, "which are...", neighbours
 
  if (current_el < 1) return
  if (elem_nmb <= 0) return
  
  ! Loop neighbours in the same style as you built them
   do i = 1, elem_nmb
      
      try_el = neighbours(i)
      print*, "Trying element", try_el

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

      
      !print*, 'hier1'
      ! Liegt die Position des Drifters im gefundenen Tetraeder?
      !x- und y-Werte des Dreiecks bestimmen
      counter = 0
      do k = 1, 4
!print*, depth3d(elemtetr(k))
         if (depth3d(elemtetr(k)).eq.deptho) then
	        counter = counter + 1
	        node = k
	     endif
      enddo
      !print*, 'counter... ',counter
      
      !print*, 'hier2' 
      if (counter.eq.1) then ! blue tetrahedra
         !print*, 'counter.eq.1'
	     !print*, elemtetr
         l = 0
         do k = 1, 4
	        if (k.ne.node) then
	           l = l + 1
	           loncoord(l) = lon3d(elemtetr(node)) + ((lon3d(elemtetr(k)) - lon3d(elemtetr(node))) * verh1)
	           latcoord(l) = lat3d(elemtetr(node)) + ((lat3d(elemtetr(k)) - lat3d(elemtetr(node))) * verh1)
	        endif
         enddo
	    !print*, loncoord
	    !print*, latcoord
      endif
      
      !print*, 'hier3'
      if (counter.eq.2) then ! green tetrahedra
         !print*, 'counter.eq.2'
	     !print*, elemtetr
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
	 

        !print *, 'hier ----- c1,c2',c1,c2 
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
      
        !print*, 'hier6'
        if ((c1.ne.0).and.(c2.ne.0)) then
	       !print*, 'war wohl nicht 2', c1, c2
	      !print*, loncoord
	      !print*, latcoord
	      !print*, lon, lat
        else
	    !print*, 'return'
	       return
	    endif
	 
	    print*, loncoord
	    print*, latcoord
     endif
      
      !print*, 'hier4'
     if (counter.eq.3) then ! red tetrahedra
        !print*, 'counter.eq.3'
	    !print*, elemtetr
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
  print*, 'wohl rausgefallen', c1, c2, counter

  find_tetra = .false.
     

end subroutine search_tetrahedron
!
!=======================================================================
!
subroutine search_tetrahedron_old(lon, lat, depth, deptho, depthu, elems, elem3d, nod3d, elem3d_nodes, elemtetr, depth3d, lon3d, lat3d, find_tetra)
   implicit none
   
        ! we know the 2d element --> depth upper and lower --> we know the prism
        ! have to find tetrahedra now

   integer:: elem3d, nod3d, elems(6), elemtetr(4), j, k, l, counter, node, sort1(4), sort2(4), c1, c2, ccwt(4)
   real(kind=8):: lon, lat, depth, deptho, depthu, verh1, verh2, loncoord(4), latcoord(4)
   integer,dimension(4,elem3d):: elem3d_nodes
   real(kind = 8), dimension(nod3d):: depth3d, lon3d, lat3d
   logical     :: find_tetra
   
   verh1 = (depth - deptho)/(depthu - deptho)
   verh2 = (depthu - depth)/(depthu - deptho)
   
   !print*, 'call search_tetrahedron'
   !print*, deptho, depthu

   find_tetra = .true.


 !  print*, 'call search_tetrahedron old'
 !  print*, 'deptho, depthu',deptho, depthu


   do j=1, elem3d
      elemtetr = elem3d_nodes(:,j)
      counter = 0
      do k = 1, 4
         do l = 1, 6
	        if (elemtetr(k) == elems(l)) counter = counter + 1
         enddo
	     if (counter.lt.k) exit
      enddo
      
      if (counter.ne.4) then
         !print*, 'lag nicht drin'
         cycle  !Tetraeder liegt nicht im vorgegebenen Prisma
      endif
      
      !print*, 'hier1'
      ! Liegt die Position des Drifters im gefundenen Tetraeder?
      !x- und y-Werte des Dreiecks bestimmen
      counter = 0
      do k = 1, 4
print*, depth3d(elemtetr(k))
         if (depth3d(elemtetr(k)).eq.deptho) then
	        counter = counter + 1
	        node = k
	     endif
      enddo
      !print*, 'counter... ',counter
      
      !print*, 'hier2' 
      if (counter.eq.1) then ! blue tetrahedra
         !print*, 'counter.eq.1'
	     !print*, elemtetr
         l = 0
         do k = 1, 4
	        if (k.ne.node) then
	           l = l + 1
	           loncoord(l) = lon3d(elemtetr(node)) + ((lon3d(elemtetr(k)) - lon3d(elemtetr(node))) * verh1)
	           latcoord(l) = lat3d(elemtetr(node)) + ((lat3d(elemtetr(k)) - lat3d(elemtetr(node))) * verh1)
	        endif
         enddo
	    !print*, loncoord
	    !print*, latcoord
      endif
      
      !print*, 'hier3'
      if (counter.eq.2) then ! green tetrahedra
         !print*, 'counter.eq.2'
	     !print*, elemtetr
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
	 

        !print *, 'hier ----- c1,c2',c1,c2 
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
      
        !print*, 'hier6'
        if ((c1.ne.0).and.(c2.ne.0)) then
	       !print*, 'war wohl nicht 2', c1, c2
	      !print*, loncoord
	      !print*, latcoord
	      !print*, lon, lat
        else
	    !print*, 'return'
	       return
	    endif
	 
	    print*, loncoord
	    print*, latcoord
     endif
      
      !print*, 'hier4'
     if (counter.eq.3) then ! red tetrahedra
        !print*, 'counter.eq.3'
	    !print*, elemtetr
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
  print*, 'wohl rausgefallen', c1, c2, counter

  find_tetra = .false.
     

end subroutine search_tetrahedron_old
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



