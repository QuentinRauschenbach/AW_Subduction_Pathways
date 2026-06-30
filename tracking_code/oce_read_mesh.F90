module g_read_mesh
  ! for grid rotation, north pole to Greenland
  !
  ! Coded by Qiang Wang
  ! Reviewed by ??
  !----------------------------------------------------------------

  use g_config
  use o_PARAM
  use o_MESH
  use o_ELEMENTS
  use o_DATA_TYPES
  use g_rotate_grid

  implicit none
  save


  contains


subroutine read_2Dmesh(MeshPath)

  implicit none
!  save

  integer           :: i, n, ind
  real(kind=8)      :: x, y
  character(100)    :: MeshPath

  open (20,file=trim(MeshPath)//'nod2d.out',  status='old')
  open (21,file=trim(MeshPath)//'elem2d.out', status='old')
  write(*,*) '2D mesh:  opened'

  read(20,*) nod2D 
  print*,'nod2d: ',nod2d
  allocate(coord_nod2D(2, nod2D))
  do i=1, nod2D
     read(20,*) n, x, y, ind
     coord_nod2D(1, i)=x
     coord_nod2D(2, i)=y
  end do
  close(20)

  read(21,*)  elem2D    
  print*,'elem2d: ',elem2d  
  allocate(elem2D_nodes(3,elem2D))  
  do n=1, elem2D
     read(21,*) elem2D_nodes(:,n)     
  end do
  close(21)
  
  write(*,*) 'read_2Dmesh: DONE'

end subroutine read_2Dmesh
!
!---------------------------------------------------------------------------
!
subroutine read_3Dmesh(MeshPath)

  implicit none
  
  integer           :: i, j, n, node, ind, m
  real(kind=8)      :: x, y, z
  character(100)    :: MeshPath

  open(10,file=trim(MeshPath)//'nod3d.out', status='old')
  open(11,file=trim(MeshPath)//'elem3d.out',status='old')
  open(12,file=trim(MeshPath)//'aux3d.out', status='old')
  open(13, file=trim(MeshPath)//'m3d.ini', status='old')

  write(*,*) '3D mesh: opened'

  ! Read node data
  read(10,*) nod3D
  allocate(coord_nod3D(3, nod3D),depth3d(nod3D))
  do i=1, nod3D 
     read(10,*) n, x, y, z, ind
     coord_nod3D(1,n)=x
     coord_nod3D(2,n)=y
     coord_nod3D(3,n)=z
  end do   
  close(10)

  depth3d = -1*coord_nod3D(3,:)
  
  ! Read the element data
  read(11, *)  elem3D   
  allocate(elem3D_nodes(4,elem3D)) 
  do i=1,elem3D 
     read(11,*) elem3D_nodes(:,i)
  end do
  close(11)

  ! Read auxilliary data

  read(12,*) z_layers
  allocate(nod3D_below_nod2D(z_layers,nod2D))
  read(12,*) nod3D_below_nod2D
  
  allocate(nod2D_corresp_to_nod3D(nod3D)) 

  do i=1, nod3D
     read(12, *) n
     nod2D_corresp_to_nod3D(i)=n
  end do

  allocate(elem2D_corresp_to_elem3D(elem3D)) 

  do i=1,elem3D
     read(12,*) n
     elem2D_corresp_to_elem3D(i)=n
  end do      
  close(12)


!Lesen der Tiefendaten: welche Ebene liegt in welcher Tiefe      !neu
 
  allocate(depths(z_layers))
  do i = 1, z_layers
   read(13, *) depths(i)
  enddo
  close(13)
  print*, 'finished reading m3d.ini'
  
  write(*,*) 'read_3Dmesh: DONE'  
end subroutine read_3Dmesh
!
!---------------------------------------------------------------------------
!

subroutine mesh_scaling  
  !
  ! Transforms degrees in rad in coord_nod2D(2,nod2D)
  ! Constructs the arrays cos_elem2D(elem2D)
  ! Constructs num_layers_below_nod2D, 
  ! and does transform to rad in coord_nod3D


  implicit none

  integer         :: i,j,ind2,ind3
  integer         :: n,  node, nodeup
  real(kind=8)    :: lon,lat,rlon,rlat

  ! =======================
  !  Lon and lat to radians
  ! =======================

  coord_nod2D(1,:)=coord_nod2D(1,:)*rad
  coord_nod2D(2,:)=coord_nod2D(2,:)*rad

  coord_nod3D(1,:)=coord_nod3D(1,:)*rad
  coord_nod3D(2,:)=coord_nod3D(2,:)*rad

  ! =======================
  ! Mean cos on 2D elements
  ! This sets spherical geometry!
  ! =======================
!  allocate(cos_elem2D(elem2D))
!  do i=1, elem2D  
 !    cos_elem2D(i)=sum(cos(coord_nod2D(2,elem2D_nodes(:,i))))/3.0
 ! end do
 ! if(cartesian) cos_elem2D=1.0  
  ! =======================
  ! number of layers 
  ! =======================
 ! allocate(num_layers_below_nod2D(nod2D))
 ! num_layers_below_nod2D=-1
 ! do n=1, nod2D
 !    do j=1,max_num_layers
 !       node=nod3D_below_nod2D(j,n)
 !       if (node > 0) then
 !          num_layers_below_nod2D(n)=num_layers_below_nod2D(n) + 1
 !       else
 !          exit
 !       end if
 !    end do
 ! end do
  ! ========================
  ! Lon and lat to radians for 3D nodes
  ! ========================
 ! do n=1,nod2D
 !    !   coord_nod3D correction:
 !    do j=1,max_num_layers
 !       node=nod3D_below_nod2D(j,n)
 !       if (node < 1) exit
 !       coord_nod3D(1,node)=coord_nod2D(1,n) ! equal x and y coords
 !       coord_nod3D(2,node)=coord_nod2D(2,n) ! equal x and y coords
 !    end do
 ! end do
  ! ========================
  ! setup geolat, which contains geographic latitude
  ! ========================
  allocate(glat(nod2d),glon(nod2d))
  allocate(lat3d(nod3d),lon3d(nod3d),lat3d_rad(nod3d),lon3d_rad(nod3d))

 !print*, 'rotated_grid: ',rotated_grid

!  if(rotated_grid) then
     do i=1,nod2d
        rlon=coord_nod2d(1,i)
        rlat=coord_nod2d(2,i)
        call r2g(lon,lat,rlon,rlat)
        glat(i)=lat
        glon(i)=lon
     end do
     do i=1,nod3d
        rlon=coord_nod3d(1,i)
        rlat=coord_nod3d(2,i)
        call r2g(lon,lat,rlon,rlat)
        lat3d(i)=lat
        lon3d(i)=lon
     end do
!  else
!     glat=coord_nod2d(2,:)
!     glon=coord_nod2d(1,:)
!     lat3d=coord_nod3d(2,:)
!     lon3d=coord_nod3d(1,:)
!  end if
  
 lon3d_rad=lon3d
 lat3d_rad=lat3d
 lon3d=lon3d/rad
 lat3d=lat3d/rad


end subroutine mesh_scaling
!
!---------------------------------------------------------------------------
!
subroutine find_layer_elem3d
  !find the layer number of 3d elements
  !
  ! Coded by Qiang Wang
  ! Reviewed by ??
  !-------------------------------------------------------------
  

  implicit none

  integer                   :: i, j, k, elem3, tetra_nodes(4)
  integer, allocatable      :: auxind(:)

  allocate(elem3d_layer(elem3D))
  allocate(auxind(nod3d))	

  do i=1,nod2d
     do k=1,z_layers
        j=nod3d_below_nod2d(k,i)  
        if(j<0) exit
        auxind(j)=k
     end do
  end do
  do elem3=1, elem3D
     tetra_nodes=elem3d_nodes(:,elem3)
     elem3d_layer(elem3)=minval(auxind(tetra_nodes))
  end do

  deallocate(auxind)
end subroutine find_layer_elem3d

!
!==========================================================================
!
subroutine find_tetra_in_prism
  !find the three tetrahedra that define prism
  !
  ! Coded by Claudia
  ! Reviewed by ??
  !-------------------------------------------------------------
  

  implicit none

  integer                              :: i, j, lay, el
  integer, allocatable, dimension(:,:) :: ind

  allocate(tetra_in_prism(z_layers,elem2D))
  allocate(ind(z_layers,elem2D))

  ind=0
  do j=1,elem3d
     lay=elem3d_layer(j) ! upper layer
     el=elem2D_corresp_to_elem3D(j)
     ind(lay,el)=ind(lay,el)+1
  enddo

  do i=1,z_layers
     do j=1,elem2D
        tetra_in_prism(i,j)%nmb=ind(i,j)
        allocate(tetra_in_prism(i,j)%addresses(ind(i,j)))
     enddo
  enddo
  ind=0
  do j=1,elem3d
     lay=elem3d_layer(j)
     el=elem2D_corresp_to_elem3D(j)
     ind(lay,el)=ind(lay,el)+1
     tetra_in_prism(lay,el)%addresses(ind(lay,el))=j
  enddo
 
  deallocate(ind)

end subroutine find_tetra_in_prism
!
! ----------------------------------------------------------------------------
!
end module g_read_mesh
