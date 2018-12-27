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
   use kinds,    only: i4, l4, r8
   use parallel, only: par, par_initialize, par_finalize, onesided_comm
   use infra,    only: map, var, gvar, infra_initialize, infra_finalize, check_var
   implicit none
!
   integer(i4) :: gnx = 10, gny = 3, nx
   logical(l4) :: issame
!
   call par_initialize()
   call read_nl()
!
   nx = infra_initialize(par%nprocs,par%rank,gnx,gny,.true.)
   if (par%ismaster) then
     print *, 'gnx        = ', gnx
     print *, 'gny        = ', gny
     print *, 'gsize (MB) = ', gnx*gny*8.0_r8/1000./1000.
     print *, 'lsize (MB) = ',  nx*gny*8.0_r8/1000./1000.
   endif
!
   call onesided_comm(gnx,gny,nx,map,var,gvar)
!
   if (par%ismaster) then
     issame = check_var(gnx,gny,gvar)
     if (issame) then
       print *, '* same !'
     else
       print *, '* diff !'
       print *, gvar(nx-3:nx+3,1)
     endif
   endif
   if (par%ismaster) then
     print *, 'gsize (MB) = ', gnx*gny*8.0_r8/1000./1000.
     print *, 'lsize (MB) = ',  nx*gny*8.0_r8/1000./1000.
   endif
!
   call infra_finalize()
   call par_finalize()
!
   contains
!-------------------------------------------------------------------------------
!
!
!
!-------------------------------------------------------------------------------
   subroutine read_nl()
!-------------------------------------------------------------------------------
   implicit none
!
   namelist/inputs/gnx, gny
!
   open(21,file='inputs.nl',status='old')
   read(21,nml=inputs)
   close(21)
!
   return
   end subroutine read_nl
!
   end program test
!-------------------------------------------------------------------------------
