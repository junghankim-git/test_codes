!-------------------------------------------------------------------------------
   program remapmain
!-------------------------------------------------------------------------------
!
!  abstract: 
!
!  history log :
!    2016-03-18  junghan kim    first written
!
!  structure:
!
!-------------------------------------------------------------------------------
   use kinds, only: i4, r8, l4
   use remapper, only: remap_t
   use remapper, only: remap_initialize, remap_finalize, remap_driver
   implicit none
!
   type(remap_t) :: remap
   character(len=512) :: remapfile, oufile, srcfile
   character(len=16) :: dimname, varname
   integer(i4),       dimension(1,1) :: ndims_vars
   character(len=16), dimension(1,1) :: dimnames_vars
!
   remapfile = '/home/jhkim/study/library/main/scrip/2.remap_matrix/cs_cs/cs_ne030_unrotated_to_cs_ne120_unrotated.nc'
   oufile    = './result.nc'
   srcfile   = '/scratch/jhkim/testbed/kim/output/2.3/ne30/gnu/10h/sw_g/101/up-20110725220000-000011.nc'
   dimname   = 'ncol'
   varname   = 'u10m'
   call readnl()
!   call remap_initialize(remap,trim(remapfile),trim(oufile),trim(srcfile),trim(dimname),trim(varname))
!   call remap_driver(remap)
!   call writeresult(remap)
   ndims_vars    = 1
   dimnames_vars = trim(dimname)
   call remap_initialize(remap,trim(remapfile))
   call remap_driver(remap,trim(srcfile),trim(oufile),1,(/trim(dimname)/),1,(/trim(varname)/), &
                                                                      ndims_vars,dimnames_vars)
   call remap_finalize(remap)
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
   namelist/inputs/remapfile,oufile,srcfile,dimname,varname
! local
   integer(i4) :: err
!
! open
   open(21,file='./inputs_remap.nl',status='old',iostat=err)
   if (err.ne.0) then
     print*,'some error in open namelist...',err
     stop
   endif
!
! read
   read(21,nml=inputs,iostat=err)
   if (err.ne.0) then
     print*,'some error in read namelist...',err
     stop
   endif
!
! close
   close(21)
   print*,' '
!
   end subroutine readnl
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   end program remapmain
