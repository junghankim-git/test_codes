!-------------------------------------------------------------------------------
   program test
   use kinds, only : i4, r4, r8
   use kim_std_output, only : kim_db_t
   use kim_std_output, only : kimdb_initialize, kimdb_add_file, kimdb_print, kimdb_finalize, kimdb_generate, kimdb_write_file
!-------------------------------------------------------------------------------
   implicit none
!
   type(kim_db_t) :: kimdb
!
   ! files
   call kimdb_initialize(kimdb,'./testdb.nc')
!
   !call kimdb_add_file(kimdb,'/scratch/jhkim/testbed/kim/output/2.2.14/ne30/gnu/10h/sw_g/101','20110725120000','2h','1h')
   !call kimdb_add_file(kimdb,'/scratch/jhkim/testbed/kim/output/2.2.15/ne30/gnu/10h/sw_g/101','20110725120000','10h','1h')
   call kimdb_add_file(kimdb,'/scratch/jhkim/testbed/kim/output/2.2.16/ne30/gnu/10h/sw_g/101','20110725120000','10h','1h')
!
   call kimdb_print(kimdb)
   call kimdb_generate(kimdb)
   call kimdb_write_file(kimdb)
!
   call kimdb_finalize(kimdb)
!
   end
!-------------------------------------------------------------------------------
