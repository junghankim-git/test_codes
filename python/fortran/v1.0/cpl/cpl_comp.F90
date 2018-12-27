!-------------------------------------------------------------------------------
   module cpl_comp
!-------------------------------------------------------------------------------
!
!  abstract : coupler model component
! 
!  history log : 
!    2016-03-18   junghan kim     initial setup
!    2016-04-14   junghan kim     the version for application to the KIM 
!                                 (KIM-CPL v1.0)
!
!-------------------------------------------------------------------------------
   use kindscpl,    only : i4, r8, i4
   use controls,    only : lunits
   use parallelcpl, only : decompose1d, syncpar, getstartcount
   use loggercpl,   only : elev, llev, printlog, tostring,                     &
                           createlogfile, closelogfile
   use coupler,     only : ncomps, id_cpl, id_atm, id_wav, id_ocn,             &
                           natms, ids_atm, component_t, copycomponent,         &
                           coupler_t, coupledvariables_t, remapmatrix_t,       &
                           setmap, global_size, setglobalsize, printglobalsize,&
                           createcoupler, deletecoupler,                       &
                           inicouplerworld, fincouplerworld,                   &
                           inienvironment, finenvironment,                     &
                           exchangevariables, putvariable, getvariable,        &
                           inicouplervariables, fincouplervariables,           &
                           iniremappingmatrix, finremappingmatrix
   use cpl_main,    only : set_cpl_main=>set, ini_cpl_main=>ini,               &
                           run_cpl_main=>run, fin_cpl_main=>fin
!
   implicit none
!
   private
!
   type(component_t) :: comp
   type(coupler_t), dimension(:), allocatable :: cpls
   integer(i4) :: myid
   integer(i4) :: log_unit
   type(remapmatrix_t), dimension(:), allocatable :: matrix
!
   public :: set,ini,inicpl,exchangecpl,run,fin
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine set(comp_in,id_in)
!-------------------------------------------------------------------------------
!
   implicit none
!
   type(component_t), intent(in   ) :: comp_in
   integer(i4),       intent(in   ) :: id_in
!
! local variables
!
   character(len=512) :: logfilename
!
   call copycomponent(comp_in,comp)
   myid     = id_in
   log_unit = lunits(myid)
   write(logfilename,'(a)'), 'cpl.log'
!
   if (comp%iscoupling) then
     call inicouplerworld(comp)
   endif
!
   call createlogfile(comp%par,trim(logfilename),log_unit)
   call printlog(comp%par,'set cpl_comp',elev,llev,un=log_unit)
!
   call set_cpl_main(comp%par)
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
   call printlog(comp%par,'ini cpl_comp',elev,llev,un=log_unit)
!
   call ini_cpl_main()
!
   call setglobalsize(comp)
!
   return
   end subroutine ini
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine inicpl()
!-------------------------------------------------------------------------------
!
   implicit none
!
! local variables
!
   integer(i4) :: i,ista,iend
   integer(i4) :: lvarsize
   integer(i4), dimension(:), allocatable :: starts,counts
!
   call printlog(comp%par,'ini coupling',elev,llev,un=log_unit)
!
   call createcoupler(comp,cpls)
!
! atm
!
   if (id_atm.gt.0) then
     if (natms.eq.1) then
       call printlog(comp%par,' global size of atm = '//tostring(global_size(id_atm)),elev,llev,un=log_unit)
       lvarsize = decompose1d(1,global_size(id_atm),comp%par%nprocs,comp%par%rank,ista,iend)
       call getstartcount(ista,iend,starts,counts)
       call setmap(comp,cpls,id_atm,starts,counts)
       deallocate(starts,counts)
     else
       do i = 1,natms
       call printlog(comp%par,' global size of atm = '//tostring(global_size(ids_atm(i))),elev,llev,un=log_unit)
       lvarsize = decompose1d(1,global_size(ids_atm(i)),comp%par%nprocs,comp%par%rank,ista,iend)
       call getstartcount(ista,iend,starts,counts)
       call setmap(comp,cpls,ids_atm(i),starts,counts)
       deallocate(starts,counts)
       enddo
     endif
   endif
!
! wav
!
   if (id_wav.gt.0) then
     call printlog(comp%par,' global size of wav = '//tostring(global_size(id_wav)),elev,llev,un=log_unit)
     lvarsize = decompose1d(1,global_size(id_wav),comp%par%nprocs,comp%par%rank,ista,iend)
     call getstartcount(ista,iend,starts,counts)
     call setmap(comp,cpls,id_wav,starts,counts)
     deallocate(starts,counts)
   endif
!
! ocn
!
   if (id_ocn.gt.0) then
     call printlog(comp%par,' global size of ocn = '//tostring(global_size(id_ocn)),elev,llev,un=log_unit)
     lvarsize = decompose1d(1,global_size(id_ocn),comp%par%nprocs,comp%par%rank,ista,iend)
     call getstartcount(ista,iend,starts,counts)
     call setmap(comp,cpls,id_ocn,starts,counts)
     deallocate(starts,counts)
   endif
!
   if (allocated(starts)) deallocate(starts)
   if (allocated(counts)) deallocate(counts)
!
   call inienvironment(comp,cpls)
!
   call iniremappingmatrix(comp,cpls,matrix)
!
   return
   end subroutine inicpl
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine exchangecpl(itime)
!-------------------------------------------------------------------------------
!
   implicit none
!
   integer(i4),intent(in) :: itime
!
   call printlog(comp%par,'do coupling('//tostring(itime)//') ',elev,llev,un=log_unit)
!
   call exchangevariables(comp,cpls,itime,matrix)
!
   return
   end subroutine exchangecpl
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine run(itime)
!-------------------------------------------------------------------------------
!
   implicit none
!
   integer(i4),intent(in) :: itime
!
   call printlog(comp%par,'run cpl_comp('//tostring(itime)//') ',elev,llev,un=log_unit)
!
   call run_cpl_main()
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
   call printlog(comp%par,'fin cpl_comp',elev,llev,un=log_unit)
!
   call fin_cpl_main()
!
   call finremappingmatrix(comp,cpls,matrix)
   call finenvironment(comp,cpls)
   call deletecoupler(comp,cpls)
   call closelogfile(comp%par,log_unit)
!
   return
   end subroutine fin
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end module cpl_comp
