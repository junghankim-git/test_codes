!-------------------------------------------------------------------------------
   program test
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
   use kinds,        only: i4, l4, r4, r8
   use cubed_sphere, only: cubed_sphere_t, cs_initialize, cs_finalize, cs_run, &
                           cs_print_coordinates, cs_print_directions,          &
                           cs_write_cube_infos, cs_write_lonslats,             &
                           cs_write_scrip
!
   type(cubed_sphere_t) :: cs1, cs2
   integer(i4)          :: np1, ne1, nprocs1
   integer(i4)          :: np2, ne2, nprocs2
   logical(l4)          :: rotated1, rotated2
   logical(l4)          :: dogenremap = .false.
   character(len=512)   :: filename1, filename2
   character(len= 16)   :: lonname1, latname1, lonname2, latname2
   character(len=512)   :: scripfile1, scripfile2
   character(len=512)   :: remapfile12, remapfile21
!
   call readnl()
!
! make the SCRIP inputs
!
   call cs_initialize(cs1,np1,ne1,nprocs1,rotated1,.true.,filename1,lonname1,latname1)
   call cs_run(cs1)
   !if (.not.dogenremap) then
     call cs_write_cube_infos(cs1,'./out1.nc')
     call cs_write_lonslats(cs1,'./lonslats1.nc')
   !endif
   call cs_write_scrip(cs1,trim(scripfile1))
   call cs_finalize(cs1)
!
   if (dogenremap) then
     call cs_initialize(cs2,np2,ne2,nprocs2,rotated2,.true.,filename2,lonname2,latname2)
     call cs_run(cs2)
     !if (.not.dogenremap) then
       call cs_write_cube_infos(cs2,'./out2.nc')
       call cs_write_lonslats(cs2,'./lonslats2.nc')
     !endif
     call cs_write_scrip(cs2,trim(scripfile2))
     call cs_finalize(cs2)
!
! make namelist for remap matrix (with SCRIP)
!
     call make_scrip_input()
   endif
!
   contains
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine readnl()
!-------------------------------------------------------------------------------
   implicit none
!
   namelist /inputs/ np1,ne1,nprocs1,filename1,lonname1,latname1,rotated1,     &
                     np2,ne2,nprocs2,filename2,lonname2,latname2,rotated2,     &
                     dogenremap,scripfile1,scripfile2,remapfile12,remapfile21
   open(31,file='inputs.nl',status='old')
   read(31,inputs)
   close(31)
!
   end subroutine readnl
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine make_scrip_input()
!-------------------------------------------------------------------------------
   implicit none
!
   character(len=16) :: str_rotated1, str_rotated2
   character(len=64) :: descrip1, descrip2
!
   if (rotated1) then
     str_rotated1 = 'rotated'
   else
     str_rotated1 = 'unrotated'
   endif
!
   if (rotated2) then
     str_rotated2 = 'rotated'
   else
     str_rotated2 = 'unrotated'
   endif
!
   write(descrip1,'(a,i3.3,a)') 'Cubed-sphere on ne',ne1,' and np4('//trim(str_rotated1)//') '
   write(descrip2,'(a,i3.3,a)') 'Cubed-sphere on ne',ne2,' and np4('//trim(str_rotated2)//') '
   open(unit=41,file='scrip_in',status='unknown')
   write(41,'(a)') "&remap_inputs"
   write(41,'(a)') " num_maps        = 2"
   write(41,'(a)') " grid1_file      = '"//trim(scripfile1)//"'"
   write(41,'(a)') " grid2_file      = '"//trim(scripfile2)//"'"
   write(41,'(a)') " interp_file1    = '"//trim(remapfile12)//"'"
   write(41,'(a)') " interp_file2    = '"//trim(remapfile21)//"'"
   write(41,'(a)') " map1_name       = '"//trim(descrip1)//"'"
   write(41,'(a)') " map2_name       = '"//trim(descrip2)//"'"
   write(41,'(a)') " "
   write(41,'(a)') " map_method      = 'conservative'"
   write(41,'(a)') " normalize_opt   = 'fracarea'"
   write(41,'(a)') " output_opt      = 'scrip'"
   write(41,'(a)') " restrict_type   = 'latlon'"
   write(41,'(a)') " num_srch_bins   = 90"
   write(41,'(a)') " luse_grid1_area = .false."
   write(41,'(a)') " luse_grid2_area = .false."
   write(41,'(a)') "/"
   close(unit=41)
!
   end subroutine make_scrip_input
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   end program test
!-------------------------------------------------------------------------------
