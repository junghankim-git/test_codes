!
   program remapmain
!-------------------------------------------------------------------------------
!
!  abstract : 
!
!  history log :
!    2016-03-18  junghan kim    initial setup
!    2017-06-02  junghan kim    apply the remapping matrix of kim
!
!  structure : 
!
!-------------------------------------------------------------------------------
   use kinds
   use remapper, only: remap_t
   use remapper, only: remap_initialize, remap_finalize, remap_driver
!-------------------------------------------------------------------------------
   implicit none
!
   integer(i4),       parameter            :: max_ndims = 5, max_nvars = 12
!
   type(remap_t)                           :: remap
   character(len=512)                      :: remapfile, infile, oufile
   integer(i4)                             :: ndims
   character(len=16), dimension(max_ndims) :: dimnames
   integer(i4)                             :: nvars
   character(len=16), dimension(max_nvars) :: varnames
   integer(i4),       dimension(max_nvars) :: ndims_vars
   character(len=16), dimension(max_ndims,max_nvars) :: dimnames_vars
!
   call read_namelist()
   call remap_initialize(remap,trim(remapfile),infile,oufile,                  &
                        ndims,dimnames,nvars,varnames,ndims_vars,dimnames_vars)
   call remap_driver(remap)
   call remap_finalize(remap)
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine read_namelist()
!-------------------------------------------------------------------------------
   implicit none
!
   character(len=16), dimension(max_ndims) :: dimnames_var01, dimnames_var02,  &
              dimnames_var03, dimnames_var04, dimnames_var05, dimnames_var06,  &
              dimnames_var07, dimnames_var08, dimnames_var09, dimnames_var10,  &
              dimnames_var11, dimnames_var12
   namelist/inputs/                                remapfile, infile, oufile,  &
                                ndims, dimnames, nvars, varnames, ndims_vars,  &
              dimnames_var01, dimnames_var02, dimnames_var03, dimnames_var04,  &
              dimnames_var05, dimnames_var06, dimnames_var07, dimnames_var08,  &
              dimnames_var09, dimnames_var10, dimnames_var11, dimnames_var12
! local variables
   integer(i4) :: id, iv
!
   do id = 1,max_ndims
     dimnames_var01(id) = 'none'; dimnames_var02(id) = 'none'
     dimnames_var03(id) = 'none'; dimnames_var04(id) = 'none'
     dimnames_var05(id) = 'none'; dimnames_var06(id) = 'none'
     dimnames_var07(id) = 'none'; dimnames_var08(id) = 'none'
     dimnames_var09(id) = 'none'; dimnames_var10(id) = 'none'
     dimnames_var11(id) = 'none'; dimnames_var12(id) = 'none'
   enddo
!
   read(*,nml=inputs)
!
   do id = 1,max_ndims
     dimnames_vars(id,01) = dimnames_var01(id)
     dimnames_vars(id,02) = dimnames_var02(id)
     dimnames_vars(id,03) = dimnames_var03(id)
     dimnames_vars(id,04) = dimnames_var04(id)
     dimnames_vars(id,05) = dimnames_var05(id)
     dimnames_vars(id,06) = dimnames_var06(id)
     dimnames_vars(id,07) = dimnames_var07(id)
     dimnames_vars(id,08) = dimnames_var08(id)
     dimnames_vars(id,09) = dimnames_var09(id)
     dimnames_vars(id,10) = dimnames_var10(id)
     dimnames_vars(id,11) = dimnames_var11(id)
     dimnames_vars(id,12) = dimnames_var12(id)
   enddo
!
   return
   end subroutine read_namelist
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end program remapmain
!-------------------------------------------------------------------------------
