!###################################################################
!#
!#
!###################################################################

module o_param
  implicit none
  save
  !    
  ! *** Fixed parameters ***
  real(kind=8), parameter  	:: pi=3.141592653589793, rad=pi/180.0 
  real(kind=8), parameter       :: rpol=6356912.0_8,  req=6378388.0_8 ! radius pole and equator

 ! real(kind=8), parameter  	:: omega=2.0*pi/(24.0*60.0*60.0)
 ! real(kind=8), parameter  	:: g=9.81                       ![m/s^2]
 ! real(kind=8), parameter  	:: r_earth=6.3675e6             ![m]
 !real(kind=8), parameter  	:: rho0=1030.                   ![kg/m^3]
 !real(kind=8), parameter 	:: rho0r=1.0/rho0 
 ! real(kind=8), parameter  	:: vcpw=4.2e6                   ![J/m^3/K]volum. heat cap. of water
 ! real(kind=8), parameter	:: small=1.0e-8                 !small value
end module o_param


!----------------------------------------------------------------------

module g_config
	use o_param
	implicit none
	save

	logical			:: rotated_grid=.true.	!option only valid for coupled model case now
	real(kind=8)		:: alphaEuler=50.*rad	![rad] Euler angles, convention:
	real(kind=8)		:: betaEuler=15.*rad	![rad] first around z, then around new x,
	real(kind=8)		:: gammaEuler=-90.*rad	![rad] then around new z.
end module g_config


!----------------------------------------------------------------------


module o_DATA_TYPES
  implicit none
  save
  !
  type addresstype
     integer                                :: nmb
     integer(KIND=4), dimension(:), pointer :: addresses
  end type addresstype

  type addresstype2
     integer                               :: nmb
     integer(KIND=4), dimension(:),pointer :: addresses
  end type addresstype2
  !
end module o_DATA_TYPES

!--------------------------------------------------------------------

module o_mesh
  
  use o_DATA_TYPES
  implicit none
  save
  

  
  integer                                      :: nod2D        
  real(kind=8), allocatable, dimension(:,:)    :: coord_nod2D  
 ! integer, allocatable, dimension(:)           :: index_nod2D  
  integer                                      :: nod3D        
  real(kind=8), allocatable, dimension(:,:)    :: coord_nod3D  
 ! integer(KIND=4), allocatable, dimension(:)   :: index_nod3D  
  real(kind=8), allocatable, dimension(:)      :: glat, glon, lat3d, lon3d, lat3d_rad, lon3d_rad
  real(kind=8), allocatable, dimension(:)      :: depth3d


  
  integer                                      :: z_layers
  real, allocatable, dimension(:)              :: depths
  integer, allocatable, dimension(:)           :: num_layers_below_nod2D
  integer(KIND=4), allocatable, dimension(:,:) :: nod3D_below_nod2D  
  integer(KIND=4), allocatable, dimension(:)   :: nod2D_corresp_to_nod3D 

  type(addresstype), allocatable, dimension(:) :: nod_in_elem2D
  type(addresstype), allocatable, dimension(:) :: nghbr_nod2D
  type(addresstype2), allocatable, dimension(:,:) :: tetra_in_prism

end module o_mesh

!----------------------------------------------------------------------------

module o_elements
  implicit none
  save
  integer                                      :: elem2D
  integer(KIND=4), allocatable, dimension(:,:) :: elem2D_nodes 
  !integer(KIND=4), allocatable, dimension(:,:) :: elem2D_nghbrs 
  integer                                      :: elem3D
  integer(KIND=4), allocatable, dimension(:,:) :: elem3D_nodes 
 ! integer(KIND=4), allocatable, dimension(:,:) :: elem3D_nghbrs  
  integer(KIND=4), allocatable, dimension(:)   :: elem2D_corresp_to_elem3D, elem3d_layer
  !

  
end module o_elements


!
