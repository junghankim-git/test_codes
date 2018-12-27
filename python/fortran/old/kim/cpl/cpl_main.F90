!-------------------------------------------------------------------------------
   module cpl_main
!-------------------------------------------------------------------------------
!
!  abstract : coupler main
! 
!  history log : 
!    2016-03-18   junghan kim     initial setup
!    2016-04-15   junghan kim     the version for application to the KIM 
!                                 (KIM-CPL v1.0)
!
!-------------------------------------------------------------------------------
   use kindscpl,    only : i4, r8, i4
   use parallelcpl, only : parallel_t, copycommunicator, decompose1d
   use coupler,     only : id_atm, id_wav
   use auxiliary,   only : createx, createvar
!
   implicit none
!
   private
!
   type(parallel_t) :: par
!
   public :: set,ini,run,fin
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine set(par_in)
!-------------------------------------------------------------------------------
!
   implicit none
!
   type(parallel_t),intent(in) :: par_in
!
   call copycommunicator(par_in,par)
!
   return
   end subroutine set
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine ini()
!-------------------------------------------------------------------------------
!
   implicit none
!
   integer(i4) :: ista,iend,i,lvarsize,gvarsize
!
#if 0
!
! atm
!
   gvarsize = nsize(id_atm)
   lvarsize = decompose1d(1,gvarsize,par%nprocs,par%rank,ista,iend)
   comp_vars(id_atm)%gvarsize = nsize(id_atm)
   comp_vars(id_atm)%lvarsize = decompose1d(1,gvarsize,par%nprocs,par%rank,ista,iend)
   allocate(comp_vars(id_atm)%dofs(lvarsize))
   allocate(comp_vars(id_atm)%cpl_x(lvarsize))
   allocate(comp_vars(id_atm)%cpl_var(lvarsize))
   comp_vars(id_atm)%cpl_x = createx(gvarsize,ista,iend)
   comp_vars(id_atm)%cpl_var = createvar(gvarsize,1,ista,iend)
   do i = ista,iend
     comp_vars(id_atm)%dofs(i-ista+1) = i
   enddo
!
! wav
!
   if (id_wav.gt.0) then
     gvarsize = nsize(id_wav)
     lvarsize = decompose1d(1,gvarsize,par%nprocs,par%rank,ista,iend)
     comp_vars(id_wav)%gvarsize = nsize(id_wav)
     comp_vars(id_wav)%lvarsize = decompose1d(1,gvarsize,par%nprocs,par%rank,ista,iend)
     allocate(comp_vars(id_wav)%dofs(lvarsize))
     allocate(comp_vars(id_wav)%cpl_x(lvarsize))
     allocate(comp_vars(id_wav)%cpl_var(lvarsize))
     comp_vars(id_wav)%cpl_x = createx(gvarsize,ista,iend)
     comp_vars(id_wav)%cpl_var = createvar(gvarsize,1,ista,iend)
     do i = ista,iend
       comp_vars(id_wav)%dofs(i-ista+1) = i
     enddo
   endif
#endif
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
#if 0
   deallocate(comp_vars(id_atm)%dofs)
   deallocate(comp_vars(id_atm)%cpl_x)
   deallocate(comp_vars(id_atm)%cpl_var)
   if (id_wav.gt.0) then
     deallocate(comp_vars(id_wav)%dofs)
     deallocate(comp_vars(id_wav)%cpl_x)
     deallocate(comp_vars(id_wav)%cpl_var)
   endif
#endif
!
   return
   end subroutine fin
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end module cpl_main
