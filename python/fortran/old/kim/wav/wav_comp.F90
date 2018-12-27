!-------------------------------------------------------------------------------
   module wav_comp
!-------------------------------------------------------------------------------
!
!  abstract : wave model component
! 
!  history log : 
!    2016-03-18   junghan kim     initial setup
!    2016-04-15   junghan kim     the version for application to the KIM 
!                                 (KIM-CPL v1.0)
!
!-------------------------------------------------------------------------------
   use kindscpl,    only : i4, r8, i4
   use controls,    only : lunits
   use parallelcpl, only : decompose1d, syncpar, getstartcount
   use loggercpl,   only : elev, llev, printlog, tostring,                     &
                           createlogfile, closelogfile
   use coupler,     only : ncomps, id_cpl, id_atm, id_wav, id_ocn,             &
                           wav_nlfile, usestridedmap,                          &
                           natms, ids_atm, component_t, copycomponent,         &
                           coupler_t, coupledvariables_t, remapmatrix_t,       &
                           setmap, global_size, setglobalsize, printglobalsize,&
                           createcoupler, deletecoupler,                       &
                           inicouplerworld, fincouplerworld,                   &
                           inienvironment, finenvironment,                     &
                           exchangevariables, putvariable, getvariable,        &
                           inicouplervariables, fincouplervariables,           &
                           iniremappingmatrix, finremappingmatrix
   use wav_main,    only : lvarsize, gvarsize, dofs,                           &
                           set_wav_main=>set, ini_wav_main=>ini,               &
                           run_wav_main=>run, fin_wav_main=>fin,               &
                           wav_u10m, wav_v10m
!
   implicit none
!
   private
!
   integer(i4)       :: myid
   integer(i4)       :: log_unit
   type(component_t) :: comp
   type(coupler_t), dimension(:), allocatable :: cpls
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
   write(logfilename,'(a)') 'wav.log'
!
   if (comp%iscoupling) then
     call inicouplerworld(comp)
   endif
!
   call createlogfile(comp%par,trim(logfilename),log_unit)
   call printlog(comp%par,'set wavecomp',elev,llev,un=log_unit)
!
   call setwavemain(comp%par,wav_nlfile)
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
   call printlog(comp%par,'ini wavecomp',elev,llev,un=log_unit)
!
   call ini_wav_main()
!
   if (comp%iscoupling) then
     call setglobalsize(comp,myid,gvarsize)
   endif
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
   integer(i4), dimension(:), allocatable :: starts,counts
!
   if (comp%iscoupling) then
     call printlog(comp%par,'ini coupling',elev,llev,un=log_unit)
     call createcoupler(comp,cpls)
     call getstartcount(dofs,starts,counts,usestridedmap)
     call setmap(comp,cpls,1,starts,counts)
     deallocate(starts,counts)
     call inienvironment(comp,cpls)
   endif
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
   integer(i4), intent(in   ) :: itime
!
   if (comp%iscoupling) then
     call printlog(comp%par,'do coupling('//tostring(itime)//') ',elev,llev,un=log_unit)
! put (user)
!     call putvariable(comp,cpls,1,id_atm,'du_wav',wav_du,itime)
!     call putvariable(comp,cpls,1,id_atm,'dv_wav',wav_dv,itime)
! exchange variable
     call exchangevariables(comp,cpls,itime)
   endif
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
   call printlog(comp%par,'run wavecomp('//tostring(itime)//') ',elev,llev,un=log_unit)
!
   if (comp%iscoupling) then
     call waitexchangevariables(comp,cpls,itime)
! get
     call getvariable(comp,cpls,1,id_atm,'u10m',wav_u10m,itime)
     call getvariable(comp,cpls,1,id_atm,'v10m',wav_v10m,itime)
   endif
!
   call run_wav_main()
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
   call printlog(comp%par,'fin wavecomp',elev,llev,un=log_unit)
   call fin_wav_main()
!
   if (comp%iscoupling) then
     call finenvironment(comp,cpls)
     call deletecoupler(comp,cpls)
   endif
!
   call closelogfile(comp%par,log_unit)
!
   return
   end subroutine fin
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end module wav_comp
