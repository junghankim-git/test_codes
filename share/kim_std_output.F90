!-------------------------------------------------------------------------------
   module kim_std_output
!-------------------------------------------------------------------------------
!
!  abstract :  bi-linear interpolation module
!
!  history log :
!    201?-??-??  junghan kim    initial setup
!    2017-02-15  junghan kim    code cleanup
!
!  structure : 
!
!-------------------------------------------------------------------------------
   use kinds,        only:i4, l4, r4, r8
   use statistics,   only:get_lnorm, get_lerror
   use kim_std_vars, only:len_filename, len_varname, max_nfiles
   use kim_std_vars, only:kim_var_t, kim_file_t, nc_check
   use kim_std_vars, only:inikimvars, setreskimfile, finkimvars
   use kim_std_vars, only:inikimfiles, finkimfiles, inikimfile, finkimfile
   use kim_std_vars, only:writekimfile, readkimfile
   use kim_std_vars, only:putvarlevel1, putvarlevels_nfiles, analysisvars
   use kim_std_vars, only:printfilestatus, printfileanalysis, printfileminmax
   use time, only:time_t, timeinterval_t, string2second
   use time, only:time_create, time_interval_create, time_inc, time_get, time_print
   use netcdf
!-------------------------------------------------------------------------------
   implicit none
!
   private
!
   ! kim output
   type kim_out_t
     integer(i4) :: nfiles
     type(kim_file_t), dimension(:), allocatable :: files ! nfiles
     logical(l4) :: usedb
     type(kim_file_t) :: filedb
     real(r8), dimension(:,:,:), allocatable :: lerr !(l1, l2, linf), nlevs, nvars
   end type kim_out_t
!
   ! kim db
   logical(l4) :: useinitial = .true.
   type kim_db_t
     integer(i4) :: nfiles
     character(len=len_filename), dimension(max_nfiles) :: filenames ! nfiles
     type(kim_file_t) :: file
   end type kim_db_t
!
   logical(l4) :: issame
!
   public :: kim_out_t, kim_db_t ! derived type
   public :: kimout_initialize, kimout_finalize, kimout_run, kimout_print ! kim out
   public :: kimdb_initialize, kimdb_add_file, kimdb_print, kimdb_finalize, kimdb_generate, kimdb_write_file ! gen db
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine kimout_initialize(kim_std, filename1, filename2, nvars, varnames, usedb)
!-------------------------------------------------------------------------------
   implicit none
!
   type(kim_out_t),                          intent(inout) :: kim_std
   character(len=*),                         intent(in   ) :: filename1
   character(len=*), optional,               intent(in   ) :: filename2
   integer(i4), optional,                    intent(in   ) :: nvars
   character(len=*), dimension(:), optional, intent(in   ) :: varnames
   logical(l4), optional,                    intent(in   ) :: usedb
! local variables
   logical(l4) :: l_usedb
   integer(i4) :: l_nfiles, l_nvars, ivar, ivar_std
   character(len=len_filename), dimension(:), allocatable :: filenames
!
   if (present(filename2)) then
     allocate(filenames(2))
     filenames(1) = filename1
     filenames(2) = filename2
     kim_std%nfiles = 2
   else
     allocate(filenames(1))
     filenames(1) = filename1
     kim_std%nfiles = 1
   endif
   l_nfiles = kim_std%nfiles
!
   if (present(usedb)) then
     l_usedb = usedb
   else
     l_usedb = .false.
   endif
   kim_std%usedb = l_usedb
!
   if (present(nvars).and.present(varnames)) then
     call inikimfiles(kim_std%files,filenames,nvars,varnames)
   else
     call inikimfiles(kim_std%files,filenames)
   endif
   if (l_usedb) then
     call inikimfile(kim_std%filedb,'/home/jhkim/study/library/main/stdkim/testdb.nc',usedb = .true.)
     call readkimfile(kim_std%filedb)
   endif
!
   deallocate(filenames)
!
   return
   end subroutine kimout_initialize
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine kimout_finalize(kim_std)
!-------------------------------------------------------------------------------
   implicit none
!
   type(kim_out_t), intent(inout) :: kim_std
!
   call finkimfiles(kim_std%files)
   deallocate(kim_std%lerr)
!
   if (kim_std%usedb) then
     call finkimfile(kim_std%filedb)
   endif
!
   return
   end subroutine kimout_finalize
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine kimout_run(kim_std)
!-------------------------------------------------------------------------------
   implicit none
!
   type(kim_out_t), intent(inout) :: kim_std
! local variables
   integer(i4) :: err_nc, nfiles, ifile, nvars, ivar, ilev, nhoriz, nlevs
   integer(i4), dimension(2) :: fileid, dimid, varid
   integer(i4), dimension(:,:), allocatable :: vars_i4_1d
   integer(i4), dimension(:,:,:), allocatable :: vars_i4_2d
   real(r8), dimension(:,:), allocatable :: vars_r8_1d
   real(r8), dimension(:,:,:), allocatable :: vars_r8_2d
   real(r8), dimension(:,:,:), allocatable :: vars_r8_2d_i
   character(len=512) :: string
!
   nfiles = kim_std%nfiles
   nvars = kim_std%files(1)%nvars
!
   ! open file
   print*,'Open files... '
   if (kim_std%nfiles<=2) then
     do ifile = 1,nfiles
       err_nc = nf90_open(trim(kim_std%files(ifile)%filename),nf90_nowrite,fileid(ifile))
       call nc_check(err_nc)
     enddo
   else
     print*,'check 1 in processing'
     stop
   endif
!
   ! dimension
   print*,'Read dimension...'
   err_nc = nf90_inq_dimid(fileid(1),'ncol',dimid(1))
   call nc_check(err_nc)
   err_nc = nf90_inquire_dimension(fileid(1),dimid(1),len=nhoriz)
   call nc_check(err_nc)
   err_nc = nf90_inq_dimid(fileid(1),'lev',dimid(1))
   if (err_nc.ne.0) then
     err_nc = nf90_inq_dimid(fileid(1),'nlev',dimid(1))
   endif
   call nc_check(err_nc)
   err_nc = nf90_inquire_dimension(fileid(1),dimid(1),len=nlevs)
   call nc_check(err_nc)
   call setreskimfile(kim_std%files(1),nhoriz,nlevs)
!
   allocate(kim_std%lerr(3,0:nlevs+1,kim_std%files(1)%nvars))
   kim_std%lerr =-9999.0
!
   if (nfiles.gt.1) then
     err_nc = nf90_inq_dimid(fileid(2),'ncol',dimid(2))
     call nc_check(err_nc)
     err_nc = nf90_inquire_dimension(fileid(2),dimid(2),len=nhoriz)
     call nc_check(err_nc)
     err_nc = nf90_inq_dimid(fileid(2),'lev',dimid(2))
     if (err_nc.ne.0) then
       err_nc = nf90_inq_dimid(fileid(2),'nlev',dimid(2))
     endif
     call nc_check(err_nc)
     err_nc = nf90_inquire_dimension(fileid(2),dimid(2),len=nlevs)
     !
     if (nhoriz.ne.kim_std%files(1)%nhoriz.or.nlevs.ne.kim_std%files(1)%nlevs) then
       print*,'check 2 in processing'
       stop
     endif
     call setreskimfile(kim_std%files(2),nhoriz,nlevs)
   endif
!
   allocate(vars_i4_1d(nhoriz,nfiles))
   allocate(vars_i4_2d(nhoriz,nlevs,nfiles))
   allocate(vars_r8_1d(nhoriz,nfiles))
   allocate(vars_r8_2d(nhoriz,nlevs,nfiles))
   allocate(vars_r8_2d_i(nhoriz,nlevs+1,nfiles))
!
   ! variables
   print*,'Read variables...'
   if (nfiles.gt.1) then
   write(string,'(a) ') '======================================================================'
   print*,trim(string)
   write(string,'(a) ') '| file 1 = '//trim(kim_std%files(1)%filename)
   print*,trim(string)
   write(string,'(a) ') '| file 2 = '//trim(kim_std%files(2)%filename)
   print*,trim(string)
   write(string,'(a) ') '----------------------------------------------------------------------'
   print*,trim(string)
   endif
!
   issame = .true.
!
   do ivar = 1,nvars
     ! read variable
     do ifile = 1,nfiles
       err_nc = nf90_inq_varid(fileid(ifile),trim(kim_std%files(ifile)%vars(ivar)%varname),varid(ifile))
       if (err_nc.eq.nf90_noerr) then
         kim_std%files(ifile)%vars(ivar)%isvar = .true.
         if (kim_std%files(ifile)%vars(ivar)%varndim==1) then
           err_nc = nf90_get_var(fileid(ifile),varid(ifile),vars_r8_1d(:,ifile))
           if (kim_std%usedb) call analysisvars(kim_std%files(ifile)%vars(ivar),kim_std%filedb,vars_r8_1d(:,ifile))
         elseif (kim_std%files(ifile)%vars(ivar)%varndim==2) then
           err_nc = nf90_get_var(fileid(ifile),varid(ifile),vars_r8_2d(:,:,ifile))
           if (kim_std%usedb) call analysisvars(kim_std%files(ifile)%vars(ivar),kim_std%filedb,vars_r8_2d(:,:,ifile))
         elseif (kim_std%files(ifile)%vars(ivar)%varndim==3) then
           err_nc = nf90_get_var(fileid(ifile),varid(ifile),vars_r8_2d_i(:,:,ifile))
           if (kim_std%usedb) call analysisvars(kim_std%files(ifile)%vars(ivar),kim_std%filedb,vars_r8_2d_i(:,:,ifile))
         endif
         call nc_check(err_nc)
       else
         kim_std%files(ifile)%vars(ivar)%isvar = .false.
         vars_r8_1d(:,ifile) =-99999.9
         vars_r8_2d(:,:,ifile) =-99999.9
         vars_r8_2d_i(:,:,ifile) =-99999.9
       endif
     enddo
  !
     ! analysis variable
     if (nfiles==2) then
       if (kim_std%files(1)%vars(ivar)%varndim==1) then
         kim_std%lerr(1,1,ivar) = get_lerror(nhoriz,vars_r8_1d(:,1),vars_r8_1d(:,2),'L1')
         kim_std%lerr(2,1,ivar) = get_lerror(nhoriz,vars_r8_1d(:,1),vars_r8_1d(:,2),'L2')
         kim_std%lerr(3,1,ivar) = get_lerror(nhoriz,vars_r8_1d(:,1),vars_r8_1d(:,2),'Linf')
         kim_std%lerr(1,0,ivar) = kim_std%lerr(1,1,ivar)
         kim_std%lerr(2,0,ivar) = kim_std%lerr(2,1,ivar)
         kim_std%lerr(3,0,ivar) = kim_std%lerr(3,1,ivar)
       elseif (kim_std%files(1)%vars(ivar)%varndim==2) then
         do ilev = 1,nlevs
           kim_std%lerr(1,ilev,ivar) = get_lerror(nhoriz,vars_r8_2d(:,ilev,1),vars_r8_2d(:,ilev,2),'L1')
           kim_std%lerr(2,ilev,ivar) = get_lerror(nhoriz,vars_r8_2d(:,ilev,1),vars_r8_2d(:,ilev,2),'L2')
           kim_std%lerr(3,ilev,ivar) = get_lerror(nhoriz,vars_r8_2d(:,ilev,1),vars_r8_2d(:,ilev,2),'Linf')
         enddo
         kim_std%lerr(1,0,ivar) = get_lerror(nhoriz,nlevs,vars_r8_2d(:,:,1),vars_r8_2d(:,:,2),'L1')
         kim_std%lerr(2,0,ivar) = get_lerror(nhoriz,nlevs,vars_r8_2d(:,:,1),vars_r8_2d(:,:,2),'L2')
         kim_std%lerr(3,0,ivar) = get_lerror(nhoriz,nlevs,vars_r8_2d(:,:,1),vars_r8_2d(:,:,2),'Linf')
       elseif (kim_std%files(1)%vars(ivar)%varndim==3) then
         do ilev = 1,nlevs+1
           kim_std%lerr(1,ilev,ivar) = get_lerror(nhoriz,vars_r8_2d_i(:,ilev,1),vars_r8_2d_i(:,ilev,2),'L1')
           kim_std%lerr(2,ilev,ivar) = get_lerror(nhoriz,vars_r8_2d_i(:,ilev,1),vars_r8_2d_i(:,ilev,2),'L2')
           kim_std%lerr(3,ilev,ivar) = get_lerror(nhoriz,vars_r8_2d_i(:,ilev,1),vars_r8_2d_i(:,ilev,2),'Linf')
         enddo
         kim_std%lerr(1,0,ivar) = get_lerror(nhoriz,nlevs+1,vars_r8_2d_i(:,:,1),vars_r8_2d_i(:,:,2),'L1')
         kim_std%lerr(2,0,ivar) = get_lerror(nhoriz,nlevs+1,vars_r8_2d_i(:,:,1),vars_r8_2d_i(:,:,2),'L2')
         kim_std%lerr(3,0,ivar) = get_lerror(nhoriz,nlevs+1,vars_r8_2d_i(:,:,1),vars_r8_2d_i(:,:,2),'Linf')
       endif
       !
       if (kim_std%files(1)%vars(ivar)%isvar.and.kim_std%files(2)%vars(ivar)%isvar) then
         if (kim_std%lerr(2, 0, ivar).eq.0.0_r8) then
           write(string,'(a1,x,a14,a,f7.2,a) ') '|',trim(kim_std%files(1)%vars(ivar)%varname),' : same...,',kim_std%lerr(2,0,ivar)*100.0_r8,' %'
         else
           write(string,'(a1,x,a14,a,f7.2,a) ') '|',trim(kim_std%files(1)%vars(ivar)%varname),' : diff...,',kim_std%lerr(2,0,ivar)*100.0_r8,' %'
           issame = .false.
         endif
       else
         write(string,'(a1,x,a14,a,a7,a) ') '|',trim(kim_std%files(1)%vars(ivar)%varname),' : no var.,','-',' %'
       endif
       print*,trim(string)
       write(string,'(a) ') '----------------------------------------------------------------------'
       print*,trim(string)
     endif ! if nfiles==2
  !
     do ifile = 1,nfiles
       if (kim_std%files(ifile)%vars(ivar)%varndim==1) then
         call putvarlevel1(kim_std%files(ifile)%vars(ivar),vars_r8_1d(:,ifile))
       elseif (kim_std%files(ifile)%vars(ivar)%varndim==2) then
         call putvarlevel1(kim_std%files(ifile)%vars(ivar),vars_r8_2d(:,:,ifile))
       elseif (kim_std%files(ifile)%vars(ivar)%varndim==3) then
         call putvarlevel1(kim_std%files(ifile)%vars(ivar),vars_r8_2d_i(:,:,ifile))
       else
         print*,'some error'
         stop
       endif
     enddo
!
   enddo ! ivar = 1,std_nvars
!
   write(string,'(a) ') ' '
   print*,trim(string)
!
   ! close file
   if (kim_std%nfiles<=1) then
     do ifile = 1,nfiles
       err_nc = nf90_close(fileid(ifile))
       call nc_check(err_nc)
     enddo
   endif
!
   deallocate(vars_i4_1d)
   deallocate(vars_r8_1d)
   deallocate(vars_i4_2d)
   deallocate(vars_r8_2d)
   deallocate(vars_r8_2d_i)
!
   return
   end subroutine kimout_run
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine kimout_print(kim_std, mylev)
!-------------------------------------------------------------------------------
   implicit none
!
   type(kim_out_t),       intent(inout) :: kim_std
   integer(i4), optional, intent(in   ) :: mylev
! local variables
   integer(i4) :: nfiles, nvars, nlevs, ivar, ilev
   character(len=64) :: fmts
   character(len=512) :: string
!
   nfiles = kim_std%nfiles
   nvars = kim_std%files(1)%nvars
   nlevs = kim_std%files(1)%nlevs
!
   if (present(mylev)) then
     if (mylev.lt.1) then
       ilev = 0
     else
       ilev = mylev
     endif
   else
     ilev = 0
   endif
!
   ! l error
   if (nfiles.gt.1) then
   write(string,'(a,i4)') ' statistics:l1, l2, linf norm error, ilev = ',ilev
   print*,trim(string)
   write(string,'(a)') '======================================================================'
   print*,trim(string)
   write(string,'(a)') '| file 1 = '//trim(kim_std%files(1)%filename)
   print*,trim(string)
   write(string,'(a)') '| file 2 = '//trim(kim_std%files(2)%filename)
   print*,trim(string)
   write(string,'(a)') '----------------------------------------------------------------------'
   print*,trim(string)
   write(string,'(a2,a14,a1,a15,a,a15,a,a15,a6)') '| ', 'varname', '|', 'L1', '|', 'L2', '|', 'Linf', '|same|'
   print*,trim(string)
   write(string,'(a) ') '======================================================================'
   print*,trim(string)
   do ivar = 1,nvars
   if (kim_std%files(1)%vars(ivar)%isvar.and.kim_std%files(2)%vars(ivar)%isvar) then
   if (kim_std%lerr(2,0,ivar).eq.0.0_r8) then
   write(string,'(a2,a14,a1,e15.6,a,e15.6,a,e15.6,a6)') '| ',kim_std%files(1)%vars(ivar)%varname,'|',kim_std%lerr(1,ilev,ivar),'|',kim_std%lerr(2,ilev,ivar),'|',kim_std%lerr(3,ilev,ivar),'|  O |'
   else
   write(string,'(a2,a14,a1,e15.6,a,e15.6,a,e15.6,a6)') '| ',kim_std%files(1)%vars(ivar)%varname,'|',kim_std%lerr(1,ilev,ivar),'|',kim_std%lerr(2,ilev,ivar),'|',kim_std%lerr(3,ilev,ivar),'|  X |'
   endif
   else
   write(string,'(a2,a14,a1,a15,a,a15,a,a15,a6)') '| ',kim_std%files(1)%vars(ivar)%varname,'|','-    ','|','-    ','|','-    ','|  - |'
   endif
   print*,trim(string)
   write(string,'(a) ') '----------------------------------------------------------------------'
   print*,trim(string)
   enddo
   write(string, '(a) ') ' '
   print*, trim(string)
   write(string, '(a) ') ' '
   print*, trim(string)
   endif
!
   call printfilestatus(kim_std%files,ilev)
   if (kim_std%usedb) then
     call printfileanalysis(kim_std%files(1),kim_std%filedb)
     call printfileminmax(kim_std%files(1),kim_std%filedb)
   endif
!
   if (kim_std%nfiles.gt.1) then
     if (issame) then
       print*,'# two files are same...'
     else
       print*,'# two files are diff...'
     endif
   endif
!
   return
   end subroutine kimout_print
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine kimdb_initialize(kim_db, filename)
!-------------------------------------------------------------------------------
   implicit none
!
   type(kim_db_t),   intent(inout) :: kim_db
   character(len=*), intent(in   ) :: filename
! local variables
   integer(i4) :: nvars, ivar, std_ivar
!
   kim_db%nfiles = 0
!
   call inikimfile(kim_db%file,trim(filename),usedb=.true.)
!
   return
   end subroutine kimdb_initialize
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine kimdb_finalize(kim_db)
!-------------------------------------------------------------------------------
   implicit none
!
   type(kim_db_t), intent(inout) :: kim_db
!
   call finkimfile(kim_db%file)
!
   return
   end subroutine kimdb_finalize
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine kimdb_add_file(kim_db, filepath, startdate, ftimes_s, outfreq_s)
!-------------------------------------------------------------------------------
   implicit none
!
   type(kim_db_t),   intent(inout) :: kim_db
   character(len=*), intent(in   ) :: filepath
   character(len=*), intent(in   ) :: startdate
   character(len=*), intent(in   ) :: ftimes_s
   character(len=*), intent(in   ) :: outfreq_s
! local variables
   integer(i4) :: nfiles, ifile, nfiles_db
   character(len=len_filename), dimension(:), allocatable :: filenames
!
   filenames = generatefilenames(filepath,startdate,ftimes_s,outfreq_s)
   nfiles = size(filenames)
!
   nfiles_db = kim_db%nfiles
   kim_db%nfiles = kim_db%nfiles+nfiles
!
   do ifile = 1,nfiles
     kim_db%filenames(nfiles_db+ifile) = filenames(ifile)
   enddo
!
   if (allocated(filenames)) deallocate(filenames)
!
   return
   end subroutine kimdb_add_file
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   function generatefilenames(filepath, startdate, ftimes_s, outfreq_s) result(filenames)
!-------------------------------------------------------------------------------
   implicit none
!
   character(len=*), intent(in   ) :: filepath
   character(len=*), intent(in   ) :: startdate
   character(len=*), intent(in   ) :: ftimes_s
   character(len=*), intent(in   ) :: outfreq_s
   character(len=len_filename), dimension(:), allocatable :: filenames
! local variables
   type(time_t) :: time
   integer(i4) :: year0, month0, day0, hour0, minute0, second0, ftimes, outfreq
   integer(i4) :: year, month, day, hour, minute, second, nfiles, ifile
   integer(i4) :: ii
!
   if (len(startdate).ne.14) then
     print*,'check startdate in genvariabledb...'
     stop
   endif
!
   ! initial time
   read(startdate(1:4),'(i4)') year0
   read(startdate(5:6),'(i2)') month0
   read(startdate(7:8),'(i2)') day0
   read(startdate(9:10),'(i2)') hour0
   read(startdate(11:12),'(i2)') minute0
   read(startdate(13:14),'(i2)') second0
!
   time = time_create(year0,month0,day0,hour0)
!
   ftimes = string2second(ftimes_s)
   outfreq = string2second(outfreq_s)
   nfiles = ftimes/outfreq
   time%interval = time_interval_create(hour = 0,minute = 0,second = outfreq)
!
   if (allocated(filenames)) deallocate(filenames)
   if (useinitial) then
     allocate(filenames(nfiles+1))
   else
     allocate(filenames(nfiles))
   endif
   !
   ii = 0
   do ifile = 1,nfiles
!
     if (useinitial.and.ifile==1) then
       ii = ii+1
       call time_get(time,year,month,day,hour,minute,second)
       write(filenames(ii),'(a,a,i4.4,i2.2,i2.2,i2.2,i2.2,i2.2,a,i6.6,a) ') trim(filepath),'/UP-',year,month,day,hour,minute,second,'-',ifile,'.nc'
     endif
     !
     ii = ii+1
     call time_inc(time)
     call time_get(time,year,month,day,hour,minute,second)
     write(filenames(ii),'(a,a,i4.4,i2.2,i2.2,i2.2,i2.2,i2.2,a,i6.6,a) ') trim(filepath),'/UP-',year,month,day,hour,minute,second,'-',ifile+1,'.nc'
!
   enddo
!
   end function generatefilenames
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine kimdb_print(kim_db)
!-------------------------------------------------------------------------------
   implicit none
!
   type(kim_db_t), intent(in   ) :: kim_db
! local variables
   integer(i4) :: ifile, nfiles, ivar, nvars
!
   nfiles = kim_db%nfiles
   nvars = kim_db%file%nvars
   print*,' '
!
   do ivar = 1,nvars
     print '(i2,x,a32,x,i2) ',ivar,trim(kim_db%file%vars(ivar)%varname),kim_db%file%vars(ivar)%varndim
   enddo
   print*,' '
!
   return
   end subroutine kimdb_print
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine kimdb_generate(kim_db)
!-------------------------------------------------------------------------------
   implicit none
!
   type(kim_db_t), intent(inout) :: kim_db
! local variables
   integer(i4) :: ifile, nfiles, ivar, nvars, nhoriz, nlevs, idx_dim
   integer(i4), dimension(:), allocatable :: fileids
   integer(i4) :: ierr, dimid, varid
   real(r8), dimension(:), allocatable :: plane
   real(r8), dimension(:,:), allocatable :: var_1d
   real(r8), dimension(:,:,:), allocatable :: var_2d
   real(r8), dimension(:,:,:), allocatable :: var_2d_i
!
   nfiles = kim_db%nfiles
   nvars = kim_db%file%nvars
!
   ! open file
   allocate(fileids(nfiles))
   do ifile = 1,nfiles
     print*,'Loading... ',trim(kim_db%filenames(ifile))
     ierr = nf90_open(trim(kim_db%filenames(ifile)),nf90_nowrite,fileids(ifile))
     call nc_check(ierr)
   enddo
!
   print*,'Get dimensions and set variables'
   do ifile = 1,nfiles
!
     ! dimensions
     ierr = nf90_inq_dimid(fileids(ifile),'ncol',dimid)
     call nc_check(ierr)
     ierr = nf90_inquire_dimension(fileids(ifile),dimid,len=nhoriz)
     call nc_check(ierr)
     ierr = nf90_inq_dimid(fileids(ifile),'lev',dimid)
     call nc_check(ierr)
     ierr = nf90_inquire_dimension(fileids(ifile),dimid,len=nlevs)
     call nc_check(ierr)
     if (ifile.eq.1) then
       call setreskimfile(kim_db%file,nhoriz,nlevs)
     else
       if (nhoriz.ne.kim_db%file%nhoriz.or.nlevs.ne.kim_db%file%nlevs) then
         print*,'error(dimension) in kimdb_generate'
         stop
       endif
     endif
   enddo
!
   ! coordinates
   ierr = nf90_inq_varid(fileids(1),'lon',varid)
   call nc_check(ierr)
   ierr = nf90_get_var(fileids(1),varid,kim_db%file%coord(:)%lon)
   call nc_check(ierr)
   ierr = nf90_inq_varid(fileids(1),'lat',varid)
   call nc_check(ierr)
   ierr = nf90_get_var(fileids(1),varid,kim_db%file%coord(:)%lat)
   call nc_check(ierr)
   ! variables
   allocate(plane(nhoriz*nfiles))
   allocate(var_1d(nhoriz,nfiles))
   allocate(var_2d(nhoriz,nlevs,nfiles))
   allocate(var_2d_i(nhoriz,nlevs+1,nfiles))
!
   print*,'Get variables and stat...'
   do ivar = 1,nvars
     idx_dim=kim_db%file%vars(ivar)%varndim
     print*,'-varname:'//trim(kim_db%file%vars(ivar)%varname)
     do ifile = 1,nfiles
       ierr = nf90_inq_varid(fileids(ifile),trim(kim_db%file%vars(ivar)%varname),varid)
       call nc_check(ierr)
       if (idx_dim.eq.1) then
         ierr = nf90_get_var(fileids(ifile),varid,var_1d(:,ifile))
       elseif (idx_dim.eq.2) then
         ierr = nf90_get_var(fileids(ifile),varid,var_2d(:,:,ifile))
       elseif (idx_dim.eq.3) then
         ierr = nf90_get_var(fileids(ifile),varid,var_2d_i(:,:,ifile))
       endif
       call nc_check(ierr)
     enddo
     ! levels
     if (idx_dim.eq.1) then
       call putvarlevels_nfiles(kim_db%file%vars(ivar),var_1d(:,:))
     elseif (idx_dim.eq.2) then
       call putvarlevels_nfiles(kim_db%file%vars(ivar),var_2d(:,:,:))
     elseif (idx_dim.eq.3) then
       call putvarlevels_nfiles(kim_db%file%vars(ivar),var_2d_i(:,:,:))
     endif
   enddo
!
   deallocate(plane)
   deallocate(var_1d)
   deallocate(var_2d)
   deallocate(var_2d_i)
   ! close file
   do ifile = 1,nfiles
     print*,'Close... ',trim(kim_db%filenames(ifile))
     ierr = nf90_close(fileids(ifile))
   enddo
!
   return
   end subroutine kimdb_generate
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine kimdb_write_file(kim_db)
!-------------------------------------------------------------------------------
   implicit none
!
   type(kim_db_t), intent(in   ) :: kim_db
!
   call writekimfile(kim_db%file)
!
   return
   end subroutine kimdb_write_file
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end module kim_std_output
!-------------------------------------------------------------------------------
