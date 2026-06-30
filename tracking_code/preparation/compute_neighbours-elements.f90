
program main

!==============================================================================
! PROGRAM MAIN
!==============================================================================
!
! PURPOSE:
!   Reads a 2D triangular mesh, builds node-to-element and element-to-element
!   connectivity, and writes neighbour information to output files.
!
! INPUT FILES:
!   - nod2d.out   : Node coordinates and indices
!   - elem2d.out  : Triangular element definitions
!
! OUTPUT FILES:
!   - neighbors.out        : Neighbouring nodes for each node
!   - elem_neighbours.out  : Neighbouring elements for each element
!
! KEY DATA STRUCTURES:
!   - nod_in_elem2D(nod2d) : Stores elements incident to each node
!   - nghbr_nod2D(nod2d)   : Stores neighbouring nodes for each node
!   - elem_nghbr(elem2d)   : Stores neighbouring elements for each element
!
! NOTES:
!   - Uses derived types with pointer arrays for connectivity.
!   - Neighbours are written in deterministic order using maxloc sorting.
!   - Supports optional rotated grids via r2g subroutine.
! 
! Original version by Qiang, 19.05.2011
! Adapted by Quentin Rauschenbach for element neighbours
!
!==============================================================================


  implicit none

  integer                              :: i, j, m, n, k, a, tr(3),tet(4)
  integer                              :: counter, el, ml(1), cnt
  integer                              :: num_iteration
  integer                              :: nod2d, elem2d
  integer, allocatable, dimension(:,:) :: elem2d_nodes
  integer, allocatable, dimension(:)   :: ind, ind2d
  integer, dimension(100)              :: aux=0
  real(kind=8)                         :: rad
  real(kind=8), allocatable            :: xynod(:,:),lon(:),lat(:)
  logical			       :: rotated_grid=.true. 
  integer :: other_el             ! an element sharing node 'a'


  character(100)                        :: meshdir='./'

  type addresstype
     integer                                   :: nmb
     integer(KIND=4), dimension(:), pointer    :: addresses
  end type addresstype
  type(addresstype), allocatable, dimension(:) :: nod_in_elem2D 
  type(addresstype), allocatable, dimension(:) :: nghbr_nod2D

  ! --- NEW: element-to-element neighbour structure
  type(addresstype), allocatable :: elem_nghbr(:)
  integer, allocatable :: ind_el(:) ! marker array sized to elem2d
  integer, allocatable :: aux_el(:) ! temporary storage for found element ids
  !------------------------------------------------------------
  ! user specification starts here

  meshdir='/albedo/work/projects/oce_rio/cwekerle/mesh/Arc08_sub/'

  rotated_grid=.true.  

  !------------------------------------------------------------

  !------------------------------------------------------------
  ! read mesh

  open(11, file=trim(meshdir)//'nod2d.out')
  read(11,*) nod2d
  allocate(xynod(nod2d,2),ind2d(nod2d))
  do n=1,nod2d
     read(11,*) i, xynod(n,1),xynod(n,2),ind2d(n)
  enddo
  close(11)

  allocate(lon(nod2d),lat(nod2d))
  rad=3.141592653589793/180.0
  xynod=xynod*rad
  if(rotated_grid) then
  do n=1,nod2d
     call r2g(lon(n), lat(n), xynod(n,1),xynod(n,2))
  end do
  else
  	lon=xynod(:,1)
	lat=xynod(:,2)
  endif
  lon=lon/rad
  lat=lat/rad

  open(12, file=trim(meshdir)//'elem2d.out')
  read(12,*) elem2d
  allocate(elem2d_nodes(3,elem2d))
  do n=1,elem2d
     read(12,*) elem2d_nodes(1:3,n)
  end do
  close(12)


  write(*,*) 'Mesh and depth files are read'


  !----------------------------------------------------------
  ! Builds nod_in_elem2D

  allocate(ind(nod2D))

  ind=0
  do j=1,elem2D
     tr=elem2D_nodes(:,j)
     ind(tr)=ind(tr)+1
  end do
  allocate(nod_in_elem2D(nod2D))
  nod_in_elem2D%nmb=ind    
  do j=1,nod2D   
     allocate(nod_in_elem2D(j)%addresses(ind(j)))
  end do
  ind=0
  do j=1,elem2D   
     tr=elem2D_nodes(:,j)
     ind(tr)=ind(tr)+1
     do k=1,3
        nod_in_elem2D(tr(k))%addresses(ind(tr(k)))=j
     end do
  end do

  ! Builds nghbr_nod2D
  allocate(nghbr_nod2D(nod2D))
  ind=0
  do j=1, nod2D
     counter=0
     do m=1,nod_in_elem2D(j)%nmb
        el=nod_in_elem2D(j)%addresses(m)
        do k=1, 3
           a=elem2D_nodes(k,el)       
           if(a==j) cycle    ! the neighbour array of a node does not contain itself!!
           if (ind(a)==0) then  
              ind(a)=1 
              counter=counter+1         
              aux(counter)=a
           end if
        end do
     end do
     nghbr_nod2D(j)%nmb=counter
     allocate(nghbr_nod2D(j)%addresses(counter))

     ! we need to sort array aux(1:counter)
     do m=counter,1,-1
        ml=maxloc(aux(1:counter))
        a=ml(1)
        nghbr_nod2D(j)%addresses(m)=aux(a)
        ind(aux(a))=0
        aux(a)=-999
     end do
  end do

!-------------------------------------------------------------
! --- NEW: Build elem_nghbr (neighbour elements for each element)

! allocate temporaries with safe upper bounds
allocate(ind_el(elem2d))
ind_el = 0
allocate(aux_el(elem2d))
aux_el = -999

allocate(elem_nghbr(elem2d)) ! allocate container for element neighbours

! loop over all elements
do el = 1, elem2d
    counter = 0 ! number of neighbours found for this element

    ! for each node of element el, collect incident elements
    do k = 1, 3 ! Loop over this element's nodes (elem2d_nodes(k,e) gives node index)
        a = elem2d_nodes(k, el)             ! node ID of the k-th node of element e

        ! For this node, loop over all elements that contain it
        do m = 1, nod_in_elem2D(a)%nmb
            other_el = nod_in_elem2D(a)%addresses(m)

            if (other_el == el) cycle       ! skip the element itself

            if (ind_el(other_el) == 0) then ! check if this element was already added
                ind_el(other_el) = 1        ! mark as added
                counter = counter + 1       ! increase neighbour count
                aux_el(counter) = other_el  !  store neighbour
            end if
        end do
    end do

    ! store count and allocate storage for this element's neighbours
    elem_nghbr(el)%nmb = counter
    if (counter > 0) then
        allocate(elem_nghbr(el)%addresses(counter))
    else
        nullify(elem_nghbr(el)%addresses)
    end if

    ! sort the temporary list aux_el(1:counter) into elem_nghbr(el)%addresses
    ! using the same maxloc approach as for nodes: take the maximum repeatedly
     
    do m = counter, 1, -1
        ml = maxloc(aux_el(1:counter))
        a = ml(1)
        elem_nghbr(el)%addresses(m) = aux_el(a)
        ind_el(aux_el(a)) = 0    ! reset marker so it's ready for next element
        aux_el(a) = -999         ! clear the aux slot
    end do

end do

! free temporaries
deallocate(ind_el)
deallocate(aux_el)

!-------------------------------------------------------------

  !-------------------------------------------------------------
  ! save neighbor file

  write(*,*) 'writing output ...'
  write(*,*) 'n2d: ' ,nod2d

  open(13,file='neighbors_rotated.out')
  do j=1,nod2d
     counter=nghbr_nod2D(j)%nmb
     write(13,'(i8)', advance="no") counter
     do m=1,counter
        write(13,'(i8)', advance="no") nghbr_nod2D(j)%addresses(m)
     end do
     do m=counter+1,9
        write(13,'(i8)', advance="no") -999
     end do
     write(13,'(i8)') -999
  end do
  close(13)
  write(*,*) 'The neighbor file is saved.'


  write(*,*) 'writing output elements ...'

  ! --- WRITE TO FILE in same style as neighbors.out ---
open(unit=14, file='elem_neighbours_rotated.out')
do el = 1, elem2d
   counter = elem_nghbr(el)%nmb

   ! write the count (same field width i8 as your other file)
   write(14,'(i8)', advance="no") counter

   ! write the neighbour element ids
   do m = 1, counter
      write(14,'(i8)', advance="no") elem_nghbr(el)%addresses(m)
   end do

   ! pad to 9 neighbour fields with -999 (same approach as node neighbours)
   do m = counter+1, 15
      write(14,'(i8)', advance="no") -999
   end do

   ! trailing -999 and newline
   write(14,'(i8)') -999
end do
close(14)

  write(*,*) "Written element neighbours to elem_neighbours.out"
  
end program main

  !
  !----------------------------------------------------------------
  !

subroutine r2g(lon, lat, rlon, rlat)

  implicit none

  real(kind=8)        :: rotate_matrix(3,3)
  real(kind=8)      :: al, be, ga, rad
  real(kind=8)      :: xr, yr, zr, xg, yg, zg
  real(kind=8), intent(out)      :: lon, lat
  real(kind=8), intent(in)       :: rlon, rlat

  real(kind=8)         	:: alphaEuler=50. 		![degree] Euler angles, convention:
  real(kind=8)         	:: betaEuler=15.  		![degree] first around z, then around new x,
  real(kind=8)		:: gammaEuler=-90.		![degree] then around new z.
  
  rad=3.141592653589793/180.0

  al=alphaEuler*rad
  be=betaEuler*rad
  ga=gammaEuler*rad

  ! rotation matrix
  rotate_matrix(1,1)=cos(ga)*cos(al)-sin(ga)*cos(be)*sin(al)
  rotate_matrix(1,2)=cos(ga)*sin(al)+sin(ga)*cos(be)*cos(al)
  rotate_matrix(1,3)=sin(ga)*sin(be)
  rotate_matrix(2,1)=-sin(ga)*cos(al)-cos(ga)*cos(be)*sin(al)
  rotate_matrix(2,2)=-sin(ga)*sin(al)+cos(ga)*cos(be)*cos(al)
  rotate_matrix(2,3)=cos(ga)*sin(be)
  rotate_matrix(3,1)=sin(be)*sin(al) 
  rotate_matrix(3,2)=-sin(be)*cos(al)  
  rotate_matrix(3,3)=cos(be)


  ! Rotated Cartesian coordinates:
  xr=cos(rlat)*cos(rlon)
  yr=cos(rlat)*sin(rlon)
  zr=sin(rlat)

  ! Geographical Cartesian coordinates:
  xg=rotate_matrix(1,1)*xr + rotate_matrix(2,1)*yr + rotate_matrix(3,1)*zr
  yg=rotate_matrix(1,2)*xr + rotate_matrix(2,2)*yr + rotate_matrix(3,2)*zr  
  zg=rotate_matrix(1,3)*xr + rotate_matrix(2,3)*yr + rotate_matrix(3,3)*zr  

  ! Geographical coordinates:
  lat=asin(zg)
  if(yg==0. .and. xg==0.) then
     lon=0.0     ! exactly at the poles
  else
     lon=atan2(yg,xg)
  end if
end subroutine r2g

