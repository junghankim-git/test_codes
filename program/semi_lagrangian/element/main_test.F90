!-------------------------------------------------------------------------------
   program testmain
!-------------------------------------------------------------------------------
   use kinds,           only: i4, l4, r8
   use logger,          only: printlog, tostring
   use io_ascii,        only: open_file, close_file, write_ascii, write_ascii_lines
   use controls,        only: np, ne, ntimelevs, dt, velocity, iexp, amp, mu, sigma
   use element,         only: element_t, nups
   use element,         only: elem_initialize, elem_finalize, elem_set_position, elem_set_state
   use element,         only: elem2up_state
   use semi_lagrangian, only: semilagrangian_t, sl_initialize, sl_finalize 
   use semi_lagrangian, only: sl_gen_ctl_volume, sl_run, sl_get_mass
!-------------------------------------------------------------------------------
   implicit none
!
! output unit number
   integer(i4) :: unitnum
! logger
   integer(i4) :: elev, llev
! domain and steps
   real(r8) :: xmin, xmax
   integer(i4) :: nstep, istep
! element
   type(element_t), dimension(:), allocatable :: elem, elem_out
! semi-lagrangian
   type(semilagrangian_t) :: sl
   logical(l4) :: mono
   integer(i4) :: bndry
! output
   real(r8), dimension(:), allocatable :: output
! logger setting
   elev = 0
   llev = 1
   unitnum = 91
! domain
   xmin =-1.0_r8
   xmax = 1.0_r8
! np,ne
   ne = 1
   np = 4
   dt = 1.0_r8
   mono = .false.
   bndry = 2
!
   call printlog('START : read namelist',elev,llev)
   call read_namelist('inputs.nl')
   call printlog('END : read namelist',elev,llev)
   call printlog(' ',elev,5)
!
   call printlog('START : initialize element and semilagrangian modules',elev,llev)
   call elem_initialize(elem)
   call sl_initialize(sl,dt,mono,bndry)
   call printlog('END : initialize element and semilagrangian modules',elev,llev)
   call printlog(' ',elev,5)
!
   call printlog('START : set positions and states,generate control volume',elev,llev)
   call elem_set_position(elem,xmin,xmax)
   call elem_set_state(elem,iexp)
   call sl_gen_ctl_volume(sl,elem)
   call printlog('END : set positions and states,generate control volume',elev,llev)
   call printlog(' ',elev,5)
!
   call printlog('START : open output file',elev,llev)
   allocate(output(nups))
   call open_file('Result.dat',unitnum)
   call write_ascii(unitnum,np)
   call write_ascii(unitnum,ne)
   call write_ascii(unitnum,nups)
   call write_ascii(unitnum,nstep)
   call write_ascii(unitnum,dt)
   call write_ascii(unitnum,velocity)
   call write_ascii(unitnum,iexp)
   call write_ascii(unitnum,amp)
   call write_ascii(unitnum,mu)
   call write_ascii(unitnum,sigma)
   call write_ascii(unitnum,nups,sl%grid)
   call write_ascii(unitnum,nups+1,sl%cell)
   call printlog('END : open output file',elev,llev)
!
   call printlog('START : save initial output',elev,llev)
   call elem2up_state(elem,output,2)
! Check
!call elem_initialize(elem_out)
!call up2elem_state(output,elem_out,2)
!call elem_difference(elem,elem_out)
!call elem_finalize(elem_out)
! Check
   call write_ascii(unitnum,nups,output)
   !call printlog('mass = '//tostring(sl_get_mass(nups,sl%cell,output)),elev,llev)
   print*,'mass = ',sl_get_mass(nups,sl%cell,output)
   call write_ascii(unitnum,sl_get_mass(nups,sl%cell,output))
!
   call printlog('END : save initial output',elev,llev)
   call printlog(' ',elev,5)
!
   do istep = 1,nstep
     call printlog('istep = '//tostring(istep),elev,llev)
     call printlog('START : run semi-lagrangian,istep = '//tostring(istep),elev,llev)
     call sl_run(sl,elem)
     call printlog('END : run semi-lagrangian,istep = '//tostring(istep),elev,llev)
     call printlog(' ',elev,5)
     !
     call printlog('START : save output',elev,llev)
     call elem2up_state(elem,output,2)
     call write_ascii(unitnum,nups,output)
     !call printlog('mass = '//tostring(sl%mass),elev,llev)
     print*,'mass = ',sl_get_mass(nups,sl%cell,output)
     call write_ascii(unitnum,sl%mass)
     call printlog('END : save output',elev,llev)
     call printlog(' ',elev,5)
   enddo
!
   call printlog('START : close output file',elev,llev)
   call close_file(unitnum)
   deallocate(output)
   call printlog('END : close output file',elev,llev)
   call printlog(' ',elev,5)
!
   call sl_finalize(sl)
   call elem_finalize(elem)
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine read_namelist(nlfilename)
!-------------------------------------------------------------------------------
   implicit none
!
   character(len=*), optional,                       intent(in   ) :: nlfilename
   namelist/inputs/xmin, xmax, nstep, np, ne, velocity, dt, mono, bndry, iexp, amp, mu, sigma
!
   if (present(nlfilename)) then
     open(21,file=trim(nlfilename),status='old')
     read(21,nml=inputs)
     close(21)
   else
     read(*,nml=inputs)
   endif
!
   call printlog('xmin = '//tostring(xmin),elev,llev+1)
   call printlog('xmax = '//tostring(xmax),elev,llev+1)
   call printlog('nstep = '//tostring(nstep),elev,llev+1)
   call printlog('ne = '//tostring(ne),elev,llev+1)
   call printlog('np = '//tostring(np),elev,llev+1)
   call printlog('velocity = '//tostring(velocity),elev,llev+1)
   call printlog('dt = '//tostring(dt),elev,llev+1)
   if (mono) then
     call printlog('mono = '//tostring(1),elev,llev+1)
   else
     call printlog('mono = '//tostring(0),elev,llev+1)
   endif
   call printlog('bndry = '//tostring(bndry),elev,llev+1)
   call printlog('iexp = '//tostring(iexp),elev,llev+1)
   call printlog('amp = '//tostring(amp),elev,llev+1)
   call printlog('mu = '//tostring(mu),elev,llev+1)
   call printlog('sigma = '//tostring(sigma),elev,llev+1)
!
   return
   end subroutine read_namelist
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end program testmain
!-------------------------------------------------------------------------------
