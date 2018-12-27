!-------------------------------------------------------------------------------
   module wav_main
!-------------------------------------------------------------------------------
!
!  abstract : wave main
! 
!  history log : 
!    2016-03-18   junghan kim     initial setup
!    2016-04-15   junghan kim     the version for application to the KIM 
!                                 (KIM-CPL v1.0)
!
!-------------------------------------------------------------------------------
   use kindscpl,    only : i4, r8, l4
   use parallelcpl, only : parallel_t, copycommunicator, decompose1d, createdof
   use coupler,     only : ncomps
   use auxiliary,   only : createx, createvar, printglobalvar
   use netcdf_cpl,  only : netcdf_t, createncfile, closencfile, adddimension,  &
                           defvariable, addvariable
!
   implicit none
!
   private
!
   type(parallel_t) :: par
   integer(i4) :: gvarsize, lvarsize
   logical(l4) :: usecontinuous
   character(len=512) :: nlfilename
   integer(i4), dimension(:), allocatable :: dofs
   real(r8),    dimension(:), allocatable :: wav_u10m, wav_v10m
   type(netcdf_t) :: my_nc
!
   public :: set, ini, run, fin
   public :: lvarsize, gvarsize, dofs
   public :: wav_u10m, wav_v10m
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine set(par_in,nlfilename_in)
!-------------------------------------------------------------------------------
!
   implicit none
!
   type(parallel_t),intent(in   ) :: par_in
   character(len=*),intent(in   ) :: nlfilename_in
!
   call copycommunicator(par_in,par)
   nlfilename = nlfilename_in
!
   return
   end subroutine set
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine read_nl()
!-------------------------------------------------------------------------------
!
   implicit none
!
   namelist /wav/ gvarsize,usecontinuous
!
   open(21,file=trim(nlfilename),status='old')
   read(21,nml=wav)
   close(21)
!
   return
   end subroutine read_nl
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine ini()
!-------------------------------------------------------------------------------
!
   implicit none
!
   integer(i4) :: ista,iend,i
!
   call read_nl()
   lvarsize = decompose1d(1,gvarsize,par%nprocs,par%rank,ista,iend)
!
   if (usecontinuous) then
     dofs = createdof(ista,iend)
   else
     dofs = createdof(gvarsize,par%rank,par%nprocs)
   endif
!
   allocate(wav_u10m(lvarsize))
   allocate(wav_v10m(lvarsize))
   wav_u10m(:) = 0.0_r8
   wav_v10m(:) = 0.0_r8
!
   call createncfile(my_nc,par,"file_wav.nc")
   call adddimension(my_nc,"ncols",gvarsize)
   call defvariable(my_nc,"wav_u10m")
!
   return
   end subroutine ini
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine run()
!-------------------------------------------------------------------------------
!
   implicit none
!
   integer(i4) :: i
!
   call addvariable(my_nc,"wav_u10m",wav_u10m,dofs)
!   do i = 1,lvarsize
!     wav_u(i)  = wav_u(i)+wav_du(i)
!     wav_v(i)  = wav_v(i)+wav_dv(i)
!     wav_du(i) = -1.0_r8*wav_du(i)
!     wav_dv(i) = -1.0_r8*wav_dv(i)
!   end do
!
   return
   end subroutine run
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine fin()
!-------------------------------------------------------------------------------
!
   implicit none
!
   deallocate(dofs)
   deallocate(wav_u10m)
   deallocate(wav_v10m)
!
   call closencfile(my_nc)
!
   return
   end subroutine fin
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end module wav_main
