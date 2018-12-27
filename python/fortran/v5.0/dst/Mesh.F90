#include "KIM.h"
!-------------------------------------------------------------------------------
   module mesh
!-------------------------------------------------------------------------------
!
!  abstract: 
!
!  history log :
!    201x-0x-xx  your   name    code comment
!
!  structure:
!
!-------------------------------------------------------------------------------
   use kiapsbase, only: iulog=>kim_iu_log, long_kind=>kim_long_kind, real_kind=>kim_real8_kind, longdouble_kind=>kim_longdouble_kind
   use physicalconstants, only: dd_pi
   use control, only: max_file_len
   use netcdf ! _external
   implicit none
!
   integer, parameter :: mxstln = 32
! ===============================
! Public methods for Mesh
! ===============================
   logical, public :: meshusemeshfile = .false.
   public :: meshopen ! must be called first
   public :: meshcubeedgecount ! called anytime afer meshopen
   public :: meshcubeelemcount ! called anytime afer meshopen
   public :: meshcubetopology ! called afer meshopen
   public :: meshsetcoordinates ! called after meshcubetopology
   public :: meshprint ! show the contents of the mesh after it has been loaded into the module
   public :: meshclose
! ===============================
! Private members
! ===============================
   integer, private, parameter :: nfaces = 6 ! number of faces on the cube
   integer, private, parameter :: ninnerelemedge = 8 ! number of edges for an interior element
   character(len=max_file_len), private :: p_mesh_file_name
   integer, private :: p_ncid
   integer, private :: p_number_elements
   integer, private :: p_number_elements_per_face
   integer, private :: p_number_blocks
   integer, private :: p_number_nodes
   integer, private :: p_number_dimensions
   integer, private :: p_number_neighbor_edges
   real(real_kind), private, allocatable :: p_node_coordinates(:,:)
   integer, private, allocatable :: p_connectivity(:,:)
! Not used will eliminate later
   integer, private :: p_elem_block_ids
! ===============================
! Private methods
! ===============================
   private :: find_side_neighbor_information
   private :: find_side_neighbors
   private :: find_corner_neighbors
   private :: find_nodal_neighbors
   private :: get_node_coordinates
   private :: get_2d_sub_coordinate_indexes
   private :: mesh_connectivity
   private :: cube_face_element_centroids
   private :: smallest_diameter_element
   private :: cube_to_cube_coordinates
   private :: sphere_to_cube_coordinates
   private :: initialize_space_filling_curve
   private :: handle_error
   private :: open_mesh_file
   private :: close_mesh_file
   private :: get_number_of_elements
   private :: get_number_of_dimensions
   private :: get_number_of_elements_per_face
   private :: get_number_of_nodes
   private :: get_number_of_element_blocks
   private :: get_node_multiplicity
   private :: get_face_connectivity
! Not used will eliminate later
   private :: get_block_ids
!
   contains
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine handle_error(status, file, line)
!-------------------------------------------------------------------------------
   use kiapsparallel, only: abortpar
   implicit none
!
   integer, intent(in   ) :: status
   character(len=*), intent(in   ) :: file
   integer, intent(in   ) :: line
!
   print*,file,':',line,':',trim(nf90_strerror(status))
   call abortpar(message = "terminating program due to netcdf error while obtaining mesh information,please see message in standard output.")
!
   end subroutine handle_error
!> Open the netcdf file containing the mesh.
!! Assign the holder to the file to p_ncid so everyone else knows
!! how to use it without passing the argument around.
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine open_mesh_file()
!-------------------------------------------------------------------------------
   implicit none
!
   integer :: status
!
   status = nf90_open(p_mesh_file_name,nf90_nowrite,p_ncid)
   if (status/= nf90_noerr) call handle_error(status,__file__,__line__)
   meshusemeshfile = .true.
!
   end subroutine open_mesh_file
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine close_mesh_file()
!-------------------------------------------------------------------------------
   implicit none
!
   integer :: status
!
   status = nf90_close(p_ncid)
   if (status/= nf90_noerr) call handle_error(status,__file__,__line__)
!
   end subroutine close_mesh_file
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   function get_number_of_dimensions() result(number_dimensions)
!-------------------------------------------------------------------------------
   implicit none
!
   integer :: number_dimensions
! local variables
   integer :: status, number_of_dim_id
! Get the id of 'num_elem', if such dimension is not there panic and quit :P
!
   status = nf90_inq_dimid(p_ncid,"num_dim",number_of_dim_id)
   if (status/= nf90_noerr) call handle_error(status,__file__,__line__)
! How many values for 'num_elem' are there?
   status = nf90_inquire_dimension(p_ncid,number_of_dim_id,len=number_dimensions)
   if (status/= nf90_noerr) call handle_error(status,__file__,__line__)
!
   end function get_number_of_dimensions
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   function get_number_of_elements() result(number_elements)
!-------------------------------------------------------------------------------
   implicit none
!
   integer :: number_elements
! local variables
   integer :: status, number_of_elements_id
! Get the id of 'num_elem', if such dimension is not there panic and quit :P
!
   status = nf90_inq_dimid(p_ncid,"num_elem",number_of_elements_id)
   if (status/= nf90_noerr) call handle_error(status,__file__,__line__)
! How many values for 'num_elem' are there?
   status = nf90_inquire_dimension(p_ncid,number_of_elements_id,len=number_elements)
   if (status/= nf90_noerr) call handle_error(status,__file__,__line__)
!
   end function get_number_of_elements
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   function get_number_of_nodes() result(number_nodes)
!-------------------------------------------------------------------------------
   implicit none
!
   integer :: number_nodes
! local variables
   integer :: status, number_of_nodes_id
! Get the id of 'num_nodes', if such dimension is not there panic and quit :P
!
   status = nf90_inq_dimid(p_ncid,"num_nodes",number_of_nodes_id)
   if (status/= nf90_noerr) call handle_error(status,__file__,__line__)
! How many values for 'num_nodes' are there?
   status = nf90_inquire_dimension(p_ncid,number_of_nodes_id,len=number_nodes)
   if (status/= nf90_noerr) call handle_error(status,__file__,__line__)
!
   end function get_number_of_nodes
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   function get_number_of_element_blocks() result(number_element_blocks)
!-------------------------------------------------------------------------------
   use kiapsparallel, only: abortpar
   implicit none
!
   integer :: number_element_blocks
! local variables
   integer :: status, number_of_element_blocks_id
! Get the id of 'num_el_blk', if such dimension is not there panic and quit :P
!
   status = nf90_inq_dimid(p_ncid,"num_el_blk",number_of_element_blocks_id)
   if (status/= nf90_noerr) call handle_error(status,__file__,__line__)
! How many values for 'num_el_blk' are there?
   status = nf90_inquire_dimension(p_ncid,number_of_element_blocks_id,len=number_element_blocks)
   if (status/= nf90_noerr) call handle_error(status,__file__,__line__)
   if (number_element_blocks/= 1) then
     if (number_element_blocks/= 6) then
       call abortpar(message = 'Reading cube-sphere from input file is not supported')
     else
       call abortpar(message = 'Number of elements blocks not exactly 1(sphere) or 6(cube) ')
     endif
   endif
!
   end function get_number_of_element_blocks
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   function get_number_of_elements_per_face() result(number_elements_per_face)
!-------------------------------------------------------------------------------
   use kiapsparallel, only: abortpar
   implicit none
!
   integer :: number_elements_per_face
   integer :: face_num ! for each of the face, we get the information
   character(len=mxstln) :: element_type ! each face is composed of elements of certain type
   integer :: number_elements_in_face ! how many elements in this face
   integer :: num_nodes_per_elem ! how many nodes in each element
   integer :: number_of_attributes ! how many attributes in the face
   integer :: status, dimension_id
!
   if (p_number_blocks==0) then
     call abortpar(message = 'get_number_of_elements_per_face called before meshopen')
   elseif (p_number_blocks==1) then ! we are in the presence of a sphere
! First we get sure the number of nodes per element is four
     status = nf90_inq_dimid(p_ncid,"num_nod_per_el1",dimension_id)
     if (status/= nf90_noerr) call handle_error(status,__file__,__line__)
     status = nf90_inquire_dimension(p_ncid,dimension_id,len=num_nodes_per_elem)
     if (status/= nf90_noerr) call handle_error(status,__file__,__line__)
     if (num_nodes_per_elem/= 4) call abortpar(message = 'Number of nodes per element is not four')
! now we check how many elements there are in the face
     status = nf90_inq_dimid(p_ncid,"num_el_in_blk1",dimension_id)
     if (status/= nf90_noerr) call handle_error(status,__file__,__line__)
     status = nf90_inquire_dimension(p_ncid,dimension_id,len=number_elements_in_face)
     if (status/= nf90_noerr) call handle_error(status,__file__,__line__)
     number_elements_per_face = number_elements_in_face
   elseif (p_number_blocks==6) then ! we are in the presence of a cube-sphere
     call abortpar(message = 'Reading a mesh for a cube-sphere is not supported')
   else
     call abortpar(message = 'Number of elements blocks not exactly 1(sphere) or 6(cube) ')
   endif
!
   end function get_number_of_elements_per_face
! This function is used to set the value of p_elem_block_ids  but such variable is never used
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   function get_block_ids(idexo) result(block_ids)
!-------------------------------------------------------------------------------
   use kiapsparallel, only: abortpar
   implicit none
!
   integer(long_kind), intent(in   ) :: idexo
   integer(long_kind) :: block_ids(p_number_blocks)
!
   end function get_block_ids
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine get_face_connectivity()
!-------------------------------------------------------------------------------
   use kiapsparallel, only: abortpar
   implicit none
!
   integer :: var_id, status
!
   status = nf90_inq_varid(p_ncid,"connect1",var_id)
   if (status/= nf90_noerr) call handle_error(status,__file__,__line__)
   status = nf90_get_var(p_ncid,var_id,p_connectivity)
   if (status/= nf90_noerr) call handle_error(status,__file__,__line__)
!
   end subroutine get_face_connectivity
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine get_node_multiplicity(node_multiplicity)
!-------------------------------------------------------------------------------
   use kiapsparallel, only: abortpar
   use dimensions, only: max_elements_attached_to_node
   implicit none
!
   integer, intent(  out) :: node_multiplicity(:)
   integer :: node_num(4)
   integer :: k, number_nodes
!
   node_multiplicity(:) = 0
   number_nodes = size(node_multiplicity)
! check this external buffer was allocated correctly
   if (number_nodes/= p_number_nodes) call abortpar(message = 'Number of nodes does not matches size of node multiplicity array')
! for each node,we have for four other nodes
   if (minval(p_connectivity)<1.or.number_nodes<maxval(p_connectivity)) then
     call abortpar(message = 'get_node_multiplicity:node number less than 1 or greater than max.')
   endif
   do k = 1,p_number_elements_per_face
     node_num = p_connectivity(:,k)
     node_multiplicity(node_num) = node_multiplicity(node_num)+1
   enddo
   if (minval(node_multiplicity)<3.or.max_elements_attached_to_node<maxval(node_multiplicity)) then
     print*,'minval(node_multiplicity) ',minval(node_multiplicity)
     print*,'maxval(node_multiplicity) ',maxval(node_multiplicity),' and max_elements_attached_to_node ',max_elements_attached_to_node
     call abortpar(message = 'get_node_multiplicity:number of elements attached to node less than 3 or greater than maximum.')
   endif
!
   end subroutine get_node_multiplicity
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine get_node_coordinates()
!-------------------------------------------------------------------------------
   use coordinatesystems, only: cartesian3d_t
   use kiapsparallel, only: abortpar
   implicit none
!
   integer :: var_id, status
!
   status = nf90_inq_varid(p_ncid,"coord",var_id)
   if (status/= nf90_noerr) call handle_error(status,__file__,__line__)
   status = nf90_get_var(p_ncid,var_id,p_node_coordinates)
   if (status/= nf90_noerr) call handle_error(status,__file__,__line__)
!
   end subroutine get_node_coordinates
! ================================================================================
!
! -----------------Internal private routines that do not use netCDF IO -----------
!
! ================================================================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine get_2d_sub_coordinate_indexes(x, y, sgnx, sgny, face_no)
!-------------------------------------------------------------------------------
   implicit none
!
   integer, intent(in   ) :: face_no
   integer, intent(  out) :: x, y
   integer, intent(  out) :: sgnx, sgny
!
   if (face_no==1.or.face_no==3) then
     x = 2
     y = 3
   elseif (face_no==2.or.face_no==4) then
     x = 1
     y = 3
   else
     x = 2
     y = 1
   endif
!
   if (face_no==1.or.face_no==4.or.face_no==5) then
     sgnx = 1
     sgny = 1
   elseif (face_no==2.or.face_no==3) then
     sgnx =-1
     sgny = 1
   else
     sgnx = 1
     sgny =-1
   endif
!
   end subroutine get_2d_sub_coordinate_indexes
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   function find_side_neighbor_information(me, nodes, nel, element_nodes, direction) result(neighbor)
!-------------------------------------------------------------------------------
   use kiapsparallel, only: abortpar
!
   integer, intent(in   ) :: me
   integer, intent(in   ) :: nodes(2)
   integer, intent(in   ) :: nel
   integer, intent(  out) :: direction
   integer, intent(in   ) :: element_nodes(nel, 4)
   integer :: i, j
   integer :: neighbor
   integer :: side_nodes(2)
!
   neighbor = 0
   do i = 1,nel
     if (i/= me) then
       do j = 1, 4
         side_nodes(1) = element_nodes(i,j)
         side_nodes(2) = element_nodes(i,mod(j,4)+1)
         if ((side_nodes(1)==nodes(2).and.side_nodes(2)==nodes(1)).or.(side_nodes(1)==nodes(1).and.side_nodes(2)==nodes(2))) then
           neighbor = i
           direction = j
           return
         endif
       enddo
       endif
   enddo
   if (neighbor==0) call abortpar(message = 'find_side_neighbor_information:neighbor not found! every side should have a neighbor.')
!
   end function find_side_neighbor_information
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   function find_nodal_neighbors(node, element_nodes, elements) result(multiplicity)
!-------------------------------------------------------------------------------
   use kiapsparallel, only: abortpar
   use dimensions, only: max_elements_attached_to_node
   implicit none
! parameters
!
   integer, intent(in   ) :: node
   integer, intent(in   ) :: element_nodes(p_number_elements, 4)
   integer, intent(  out) :: elements(2*max_elements_attached_to_node)
! return value
   integer :: multiplicity
! local variables
   integer :: i, j
!
   multiplicity = 0
   do i = 1,p_number_elements
     do j = 1, 4
       if (element_nodes(i, j)==node) then
         multiplicity = multiplicity+1
         elements(multiplicity) = i
       endif
       enddo
   enddo
   if (multiplicity<3.or.max_elements_attached_to_node<multiplicity) then
     call abortpar(message = 'find_nodal_neighbors:number of elements attached to node less than 3 or greater than maximum.')
   endif
!
   end function find_nodal_neighbors
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine mesh_connectivity(connect)
!-------------------------------------------------------------------------------
   use kiapsparallel, only: abortpar
   implicit none
!
   integer, intent(  out) :: connect(p_number_elements, 4)
   integer :: k, j
!
   if (0==p_number_blocks) call abortpar(message = 'mesh_connectivity called before meshopen')
   j = 0
   do k = 1,p_number_elements_per_face
     j = j+1
     connect(j,:) = p_connectivity(:,k)
   enddo
   if (j/= p_number_elements) call abortpar(message = 'mesh_connectivity:number of elements in side sets not equal to total elements')
   if (minval(connect)<1.or.maxval(connect)>p_number_nodes) then
     call abortpar(message = 'mesh_connectivity:node number out of bounds')
   endif
!
   end subroutine mesh_connectivity
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine find_side_neighbors(gridvertex, normal_to_homme_ordering, element_nodes)
!-------------------------------------------------------------------------------
   use coordinatesystems, only: cartesian3d_t
   use gridgraph, only: gridvertex_t
   use kiapsparallel, only: abortpar
   implicit none
!
   integer, intent(in   ) :: normal_to_homme_ordering(8)
   integer, intent(in   ) :: element_nodes(p_number_elements, 4)
   type(gridvertex_t), intent(inout) :: gridvertex(:)
   integer :: side_nodes(2)
   integer :: neighbor, direction
   integer :: j, k, l
!
   if (0==p_number_blocks) call abortpar(message = 'find_side_neighbors called before meshopen')
   do k = 1,p_number_elements
     do l = 1, 4
       allocate(gridvertex(k)%nbrs(l)%n(1))
       allocate(gridvertex(k)%nbrs(l)%f(1))
       gridvertex(k)%nbrs(l)%n(1) = 0
     enddo
   enddo
   do k = 1,p_number_elements
     do l = 1, 4
       j = normal_to_homme_ordering(l)
       if (gridvertex(k)%nbrs(j)%n(1)==0) then
         side_nodes(1) = element_nodes(k,l)
         side_nodes(2) = element_nodes(k,mod(l,4)+1)
         neighbor = find_side_neighbor_information(k,side_nodes,p_number_elements,element_nodes,direction)
         j = normal_to_homme_ordering(l)
         gridvertex(k)%nbrs(j)%n(1) = neighbor
         j = normal_to_homme_ordering(direction)
         gridvertex(neighbor)%nbrs(j)%n(1) = k
       endif
     enddo
   enddo
   do k = 1,p_number_elements
     do l = 1, 4
       if (0==gridvertex(k)%nbrs(l)%n(1)) then
         call abortpar(message = 'Found one side of one element witout a neighbor. bummer!')
       endif
       enddo
   enddo
!
   end subroutine find_side_neighbors
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   function smallest_diameter_element(element_nodes) result(min_diameter)
!-------------------------------------------------------------------------------
   use kiapsparallel, only: abortpar
   implicit none
!
   integer, intent(in   ) :: element_nodes(:,:)
   integer :: i, j
   integer :: node_numbers(4)
   real(real_kind) :: coordinates(4, 3)
   real :: x(3), y(3), r(3), d, min_diameter
!
   if (size(element_nodes, dim=1)/= p_number_elements) then
     call abortpar(message = 'smallest_diameter_element:element count check failed in exodus_mesh. connectivity array length not equal to number of elements.')
   endif
!
   if (p_number_elements_per_face/= p_number_elements) then
     call abortpar(message = 'smallest_diameter_element:element count check failed in exodus_mesh. element array length not equal to sum of face.')
   endif
!
   min_diameter = 9999999.
   do i = 1,p_number_elements
     node_numbers = element_nodes(i,:)
     coordinates = p_node_coordinates(node_numbers,:)
! smallest side length
     do j = 1,4
       x = coordinates(j,:)
       y = coordinates(1+mod(j,4),:)
       r = x-y
       d = dot_product(r,r)
       if (d<min_diameter) then
         min_diameter = d
       endif
     enddo
! smallest diameter length
     do j = 1,2
       x = coordinates(j,:)
       y = coordinates(2+mod(j,4),:)
       r = x-y
       d = dot_product(r,r)
       if (d<min_diameter) then
         min_diameter = d
       endif
     enddo
   enddo
   min_diameter = sqrt(min_diameter)
!
   end function smallest_diameter_element
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine cube_to_cube_coordinates(cube_coor, node_coor, face_number)
!-------------------------------------------------------------------------------
   use kiapsparallel, only: abortpar
   use physicalconstants, only: dd_pi
   implicit none
!
   real(real_kind), intent(in   ) :: node_coor(4, 3)
   integer, intent(in   ) :: face_number
   real(real_kind), intent(  out) :: cube_coor(4, 2)
   real(real_kind) :: test_coor(4, 2)
   integer :: i, j, x_index, y_index, sgnx, sgny
!
   call get_2d_sub_coordinate_indexes(x_index,y_index,sgnx,sgny,face_number)
   cube_coor(:,1) = sgnx*node_coor(:,x_index)
   cube_coor(:,2) = sgny*node_coor(:,y_index)
!
   end subroutine cube_to_cube_coordinates
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine sphere_to_cube_coordinates(cube_coor, node_coor, face_number)
!-------------------------------------------------------------------------------
   use coordinatesystems, only: cartesian3d_t, cartesian2d_t, spherical_polar_t, change_coordinates, sphere2cubedsphere
   implicit none
!
   real(real_kind), intent(in   ) :: node_coor(4, 3)
   integer, intent(in   ) :: face_number
   real(real_kind), intent(  out) :: cube_coor(4, 2)
   integer :: i, l
   type(cartesian2d_t) :: cart(4)
!
   do i = 1, 4
     cart(i) = sphere2cubedsphere(change_coordinates(node_coor(i,:)),face_number)
   enddo
!
   cube_coor(:,1) = cart(:)%x
   cube_coor(:,2) = cart(:)%y
!
   end subroutine sphere_to_cube_coordinates
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine cube_face_element_centroids(centroids, face_numbers, element_nodes)
!-------------------------------------------------------------------------------
   use kiapsparallel, only: abortpar
   implicit none
!
   integer, intent(in   ) :: element_nodes(:,:)
   integer, intent(in   ) :: face_numbers(p_number_elements)
   real, intent(  out) :: centroids(p_number_elements, 2)
   real(real_kind) :: coordinates(4, 3)
   real(real_kind) :: cube_coor(4, 2)
   integer :: i, j, node_numbers(4)
!
   if (0==p_number_blocks) call abortpar(message = 'cube_face_element_centroids called before meshopen')
   if (size(element_nodes,dim=1)/= p_number_elements) then
     call abortpar(message = 'cube_face_element_centroids:element count check failed in exodus_mesh. connectivity array length not equal to number of elements.')
   endif
   if (p_number_elements_per_face/= p_number_elements) then
     call abortpar(message = 'cube_face_element_centroids:element count check failed in exodus_mesh. element array length not equal to sum of face.')
   endif
   do i = 1,p_number_elements
     node_numbers = element_nodes(i,:)
     coordinates = p_node_coordinates(node_numbers,:)
     if (6==p_number_blocks) then
       call cube_to_cube_coordinates(cube_coor,coordinates,face_numbers(i))
     else
       call sphere_to_cube_coordinates(cube_coor,coordinates,face_numbers(i))
     endif
     centroids(i,:) = sum(cube_coor,dim=1)/4.0
   enddo
!
   end subroutine cube_face_element_centroids
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine initialize_space_filling_curve(gridvertex, element_nodes)
!-------------------------------------------------------------------------------
   use gridgraph, only: gridvertex_t
   use kiapsparallel, only: abortpar
   use spacecurve, only: genspacecurve
   implicit none
!
   type(gridvertex_t), intent(inout) :: gridvertex(:)
   integer, intent(in   ) :: element_nodes(:,:)
   integer, allocatable :: mesh2(:,:), mesh2_map(:,:), sfcij(:,:)
   real :: centroids(p_number_elements, 2)
   integer :: face_numbers(p_number_elements)
   real :: x, y, h
   integer :: i, j, i2, j2, ne, ne2
   integer :: sfc_index, face, nelem
!
   if (size(gridvertex)/= p_number_elements) then
     call abortpar(message = 'initialize_space_filling_curve:element count check failed in exodus_mesh. vertex array length not equal to number of elements.')
   endif
!
   if (size(element_nodes, dim=1)/= p_number_elements) then
     call abortpar(message = 'initialize_space_filling_curve:element count check failed in exodus_mesh. connectivity array length not equal to number of elements.')
   endif
!
   face_numbers(:) = gridvertex(:)%face_number
   h = smallest_diameter_element(element_nodes)
   call cube_face_element_centroids(centroids,face_numbers,element_nodes)
   if (h<.00001) call abortpar(message = 'initialize_space_filling_curve:unreasonably small element found. less than .00001')
   ne = ceiling(0.5*dd_pi/(h/2)) ;
! find the smallest ne2 which is a power of 2 and ne2>ne
   ne2 = 2**ceiling(log(real(ne))/log(2d0))
   if (ne2<ne) call abortpar(message = 'initialize_space_filling_curve:fatel sfc error')
   allocate(mesh2(ne2,ne2))
   allocate(mesh2_map(ne2,ne2))
   allocate(sfcij(0:ne2*ne2,2))
! create a reverse index array for Mesh2
! j = Mesh2(i,j)
! (i,j) = (sfcij(j,1),sfci(j,2))
   call genspacecurve(mesh2) ! sfc partition for ne2
   do j2 = 1,ne2
     do i2 = 1, ne2
       j = mesh2(i2,j2)
       sfcij(j,1) = i2
       sfcij(j,2) = j2
     enddo
   enddo
   gridvertex(:)%spacecurve =-1
   sfc_index = 0
   do face = 1,nfaces
! associate every element on the ne x ne mesh (Mesh)
! with its closest element on the ne2 x ne2 mesh (Mesh2)
! Store this as a map from Mesh2 -> Mesh in Mesh2_map.
! elements in Mesh2 which are not mapped get assigned a value of 0
     mesh2_map = 0
     do i = 1,p_number_elements
       if (face_numbers(i)==face) then
         x = centroids(i,1)
         y = centroids(i,2)
! map this element to an (i2,j2) element
! [ -DD_PI/4,DD_PI/4 ]  -> [ 0,ne2 ]
         i2 = nint((0.5+2.0*x/dd_pi)*ne2+.5)
         j2 = nint((0.5+2.0*y/dd_pi)*ne2+.5)
         if (face==4.or.face==6) i2 = ne2-i2+1
         if (face==1.or.face==2.or.face==6) j2 = ne2-j2+1
         if (i2<1) i2 = 1
         if (i2>ne2) i2 = ne2
         if (j2<1) j2 = 1
         if (j2>ne2) j2 = ne2
         mesh2_map(i2,j2) = i
       endif
     enddo
! generate a SFC for Mesh with the same ordering as the
! elements in Mesh2 which map to Mesh.
     do j = 0,ne2*ne2-1
       i2 = sfcij(j,1)
       j2 = sfcij(j,2)
       i = mesh2_map(i2,j2)
       if (i/= 0) then
! (i2,j2) element maps to element
         gridvertex(i)%spacecurve = sfc_index
         sfc_index = sfc_index+1
       endif
     enddo
   enddo
   deallocate(mesh2)
   deallocate(mesh2_map)
   deallocate(sfcij)
   if (minval(gridvertex(:)%spacecurve)==-1) then
     do i = 1, p_number_elements
       if (-1==gridvertex(i)%spacecurve) then
         write(*,*) " error in projecting element ",i," to space filling curve."
         write(*,*) " face:",face_numbers(i)
         write(*,*) " centroid:",centroids(i,:)
       endif
       enddo
         call abortpar(message = 'initialize_space_filling_curve:vertex not on spacecurve')
   endif
!
   end subroutine initialize_space_filling_curve
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine find_corner_neighbors(gridvertex, normal_to_homme_ordering, element_nodes)
!-------------------------------------------------------------------------------
   use kiapsparallel, only: abortpar
   use gridgraph, only: gridvertex_t
   use dimensions, only: max_elements_attached_to_node
   implicit none
!
   type(gridvertex_t), intent(inout) :: gridvertex(:)
   integer, intent(in   ) :: normal_to_homme_ordering(8)
   integer, intent(in   ) :: element_nodes(p_number_elements, 4)
   integer :: node_elements(2*max_elements_attached_to_node)
   integer :: elem_neighbor(2*max_elements_attached_to_node)
   integer :: multiplicity
   integer :: i, j, k, l
!
   node_elements(:) = 0
   elem_neighbor(:) = 0
   do i = 1,p_number_elements
     do j = 1, 4
       multiplicity = find_nodal_neighbors(element_nodes(i,j),element_nodes,node_elements)
       k = 0
       do l = 1,multiplicity
         if (i/= node_elements(l).and.gridvertex(i)%nbrs(1)%n(1)/= node_elements(l).and.gridvertex(i)%nbrs(2)%n(1)/= node_elements(l).and.gridvertex(i)%nbrs(3)%n(1)/= node_elements(l).and.gridvertex(i)%nbrs(4)%n(1)/= node_elements(l)) then
           k = k+1
           elem_neighbor(k) = node_elements(l)
         endif
       enddo
       if (0<k) then
         l = normal_to_homme_ordering(j+4)
         allocate(gridvertex(i)%nbrs(l)%n(k))
         allocate(gridvertex(i)%nbrs(l)%f(k))
         gridvertex(i)%nbrs(l)%n(:) = 0
         gridvertex(i)%nbrs(l)%f(:) = 0
         gridvertex(i)%nbrs(l)%n(1:k) = elem_neighbor(1:k)
         gridvertex(i)%nbrs(l)%f(1:k) = gridvertex(elem_neighbor(1:k))%face_number
       endif
     enddo
   enddo
!
   end subroutine find_corner_neighbors
! ================================================================================
!
! -------------------------------Public Methods-----------------------------------
!
! ================================================================================
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine meshopen(mesh_file_name, par)
!-------------------------------------------------------------------------------
   use kiapsparallel, only: parallel_t, abortpar
   implicit none
!
   character(len=*), intent(in   ) :: mesh_file_name
   type(parallel_t), intent(in   ) :: par
   integer :: ncid
   integer, allocatable :: node_multiplicity(:)
   integer :: k
!
   p_mesh_file_name = mesh_file_name
   call open_mesh_file()
   p_number_elements = get_number_of_elements()
   p_number_nodes = get_number_of_nodes()
   p_number_blocks = get_number_of_element_blocks()
   p_number_dimensions = get_number_of_dimensions()
   if (p_number_dimensions/= 3) then
     call abortpar(message = 'The number of dimensions must be 3,otherwise the mesh algorithms will not work')
   endif
! Only spheres are allowed in input files.
   if (par%ismasterproc) then
     if (p_number_blocks==1) then
       write(iulog,*) "since the mesh file has only one block,it is assumed to be a sphere."
     endif
   endif
   if (p_number_blocks/= 1) then
     call abortpar(message = 'Number of elements blocks not exactly 1(sphere) ')
   endif
   p_number_elements_per_face = get_number_of_elements_per_face()
! Because all elements are in one face,this value must match  p_number_elements
   if (p_number_elements/= p_number_elements_per_face) then
     call abortpar(message = 'The value of the total number of elements does not match all the elements found in face 1')
   endif
   allocate(p_connectivity(4,p_number_elements_per_face))
   p_connectivity(:,:) = 0
! extract the connectivity from the netcdf file
   call get_face_connectivity()
   allocate(node_multiplicity(p_number_nodes))
   call get_node_multiplicity(node_multiplicity)
! tricky:  For each node with multiplicity n,there are n(n-1) neighbor links
! created.  But this counts each edge twice,so:  n(n-1) -n
! Should be the same as SUM(SIZE(GridVertex(i)%nbrs(j)%n),i=1:p_number_elements,j=1:8)
! p_number_neighbor_edges = dot_product(mult,mult) - 2*sum(mult)
   p_number_neighbor_edges = 0
   do k = 1,p_number_nodes
     p_number_neighbor_edges = p_number_neighbor_edges+node_multiplicity(k)*(node_multiplicity(k)-2)
   enddo
   deallocate(node_multiplicity)
! allocate the space for the coordinates,this is used in many functions
   allocate(p_node_coordinates(p_number_nodes,p_number_dimensions))
   call get_node_coordinates()
   if (p_number_elements_per_face/= p_number_elements) then
     call abortpar(message = 'MeshOpen:total number of elements not equal to the number of elements on face 1!')
   endif
!
   end subroutine meshopen
! This routine acts as a destructor cleaning the memory allocated in MeshOpen
! which acts as a constructor allocated dynamical memory for the nodes coordinates.
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine meshclose
! release memory
!
   deallocate(p_node_coordinates)
   deallocate(p_connectivity)
! let the file go
   call close_mesh_file()
!
   end subroutine meshclose
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine meshprint(par)
!-------------------------------------------------------------------------------
   use kiapsparallel, only: parallel_t, abortpar
   implicit none
!
   type(parallel_t), intent(in   ) :: par
!
   if (par%ismasterproc) then
     print*,'This are the values for file ',trim(p_mesh_file_name)
     print*,'The value for the number of dimensions(num_dim) is ',p_number_dimensions
     print*,'The number of elements in the mesh file is ',p_number_elements
     print*,'The number of nodes in the mesh file is ',p_number_nodes
     print*,'The number of blocks in the mesh file is ',p_number_blocks
     print*,'The number of elements in the face 1(sphere) is ',p_number_elements_per_face
     if (p_number_elements==p_number_elements) then
       print*,'The value of the total number of elements does match all the elements found in face 1(the only face) '
     else
       print*,'The value of the total number of elements does not match all the elements found in face 1'
       print*,'This message should not be appearing,there is something wrong in the code'
     endif
     print*,'The number of neighbor edges ',p_number_neighbor_edges
!print *,'The node connectivity are (compare with ncdump -v connect1) ',p_connectivity
!print *,' ========================================================='
!print *,'The node coordinates are (compare with ncdump -v coord) ',p_node_coordinates
   endif
!
   end subroutine meshprint
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine meshcubetopology(gridedge, gridvertex)
!-------------------------------------------------------------------------------
   use kiapsparallel, only: abortpar
   use dimensions, only: np
   use coordinatesystems, only: cartesian3d_t, cube_face_number_from_cart, cube_face_number_from_sphere
   use gridgraph, only: gridvertex_t
   use gridgraph, only: gridedge_t
   use cube, only: cubesetupedgeindex
   use gridgraph, only: initgridedge
   use control, only: north, south, east, west, neast, seast, swest, nwest
   implicit none
!
   type(gridedge_t), intent(  out) :: gridedge(:)
   type(gridvertex_t), intent(  out) :: gridvertex(:)
   real(real_kind) :: coordinates(4, 3)
   real(real_kind) :: centroid(3)
   type(cartesian3d_t) :: face_center
   integer :: i, j, k, l, m
   integer :: element_nodes(p_number_elements, 4)
   integer :: edgewgtp, cornerwgt
   integer :: normal_to_homme_ordering(8)
   integer :: node_numbers(4)
!
   normal_to_homme_ordering(1) = south
   normal_to_homme_ordering(2) = east
   normal_to_homme_ordering(3) = north
   normal_to_homme_ordering(4) = west
   normal_to_homme_ordering(5) = swest
   normal_to_homme_ordering(6) = seast
   normal_to_homme_ordering(7) = neast
   normal_to_homme_ordering(8) = nwest
   if (size(gridvertex)/= p_number_elements) then
     call abortpar(message = 'MeshCubeTopology:element count check failed in exodus_mesh. vertex array length not equal to number of elements.')
   endif
   if (p_number_elements_per_face/= p_number_elements) then
     call abortpar(message = 'MeshCubeTopology:element count check failed in exodus_mesh. element array length not equal to sum of face.')
   endif
   edgewgtp = np
   cornerwgt = 1
   call mesh_connectivity(element_nodes)
   do i = 1,p_number_elements
     gridvertex(i)%number = i
     gridvertex(i)%face_number = 0
     gridvertex(i)%processor_number = 0
     gridvertex(i)%spacecurve = 0
     gridvertex(i)%wgtp(:) = 0
     gridvertex(i)%wgtp_ghost(:) = 1
     do j = 1,8
       nullify(gridvertex(i)%nbrs(j)%n)
       nullify(gridvertex(i)%nbrs(j)%f)
     enddo
   enddo
! side neighbors
   call find_side_neighbors(gridvertex,normal_to_homme_ordering,element_nodes)
! All side and corner weights
   do i = 1,p_number_elements
     gridvertex(i)%wgtp(:) = edgewgtp
   enddo
! vertex faces
   do i = 1,p_number_elements
     node_numbers = element_nodes(i,:)
     coordinates = p_node_coordinates(node_numbers,:)
     centroid = sum(coordinates,dim=1)/4.0
     face_center%x = centroid(1)
     face_center%y = centroid(2)
     face_center%z = centroid(3)
     gridvertex(i)%face_number = cube_face_number_from_cart(face_center)
   enddo
! neighbor faces
   do i = 1,p_number_elements
     do j = 1, 4
       k = normal_to_homme_ordering(j)
       l = gridvertex(i)%nbrs(k)%n(1)
       gridvertex(i)%nbrs(k)%f(1) = gridvertex(l)%face_number
     enddo
   enddo
! corner neighbor and faces
   call find_corner_neighbors(gridvertex,normal_to_homme_ordering,element_nodes)
   do i = 1,p_number_elements
     do j = 5, 8
       if (associated(gridvertex(i)%nbrs(j)%n)) gridvertex(i)%wgtp(j) = cornerwgt
     enddo
   enddo
   call initgridedge(gridedge,gridvertex)
   do i = 1,size(gridedge)
     call cubesetupedgeindex(gridedge(i))
   enddo
   call initialize_space_filling_curve(gridvertex,element_nodes)
!
   end subroutine meshcubetopology
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine meshsetcoordinates(elem)
!-------------------------------------------------------------------------------
   use element, only: element_t
   use kiapsparallel, only: abortpar
   use coordinatesystems, only: cartesian3d_t, cartesian2d_t, spherical_polar_t, change_coordinates, sphere2cubedsphere
   implicit none
!
   type(element_t), intent(inout) :: elem(:)
   integer :: connectivity(p_number_elements, 4)
   integer :: node_multiplicity(p_number_nodes)
   integer :: face_no, i, k, l
   integer :: number
   integer :: node_num(4)
   real(real_kind) :: coordinates(4, 3)
   real(real_kind) :: cube_coor(4, 2)
   real(real_kind) :: x_double
   real :: x_real
   type(cartesian2d_t) :: cart2
!
   connectivity = 0
   node_multiplicity = 0
   call mesh_connectivity(connectivity)
   do k = 1,p_number_elements
     node_num = connectivity(k,:)
     node_multiplicity(node_num(:)) = node_multiplicity(node_num(:))+1
   enddo
   do k = 1,size(elem)
     number = elem(k)%vertex%number
     face_no = elem(k)%vertex%face_number
     node_num = connectivity(number,:)
     coordinates = p_node_coordinates(node_num,:)
     if (6==p_number_blocks) then
       call cube_to_cube_coordinates(cube_coor,coordinates,face_no)
     else
       call sphere_to_cube_coordinates(cube_coor,coordinates,face_no)
     endif
     elem(k)%node_numbers = node_num
     elem(k)%node_multiplicity(:) = node_multiplicity(node_num(:))
     elem(k)%corners(:)%x = cube_coor(:,1)
     elem(k)%corners(:)%y = cube_coor(:,2)
   enddo
!
   end subroutine meshsetcoordinates
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   function meshcubeedgecount() result(nedge)
!-------------------------------------------------------------------------------
   use kiapsparallel, only: abortpar
   implicit none
!
   integer :: nedge
!
   if (0==p_number_blocks) call abortpar(message = 'MeshCubeEdgeCount called before meshopenmesh')
   if (meshusemeshfile) then
! should be the same as SUM(SIZE(GridVertex(i)%nbrs(j)%n),i=1:p_number_elements,j=1:nInnerElemEdge)
! the total number of neighbors.
     nedge = p_number_neighbor_edges
   else
     call abortpar(message = 'Error in meshcubeedgecount:should not call for non-exodus mesh file.')
   endif
!
   end function meshcubeedgecount
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   function meshcubeelemcount() result(nelem)
!-------------------------------------------------------------------------------
   use kiapsparallel, only: abortpar
   implicit none
!
   integer :: nelem
!
   if (0==p_number_blocks) call abortpar(message = 'MeshCubeElemCount called before meshopenmesh')
   if (meshusemeshfile) then
     nelem = p_number_elements
   else
     call abortpar(message = 'Error in meshcubeelemcount:should not call for non-exodus mesh file.')
   endif
!
   end function meshcubeelemcount
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine test_private_methods
!-------------------------------------------------------------------------------
   implicit none
!
   integer :: element_nodes(p_number_elements, 4)
!
   call mesh_connectivity(element_nodes)
!
   end subroutine test_private_methods
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   end module mesh
!-------------------------------------------------------------------------------
