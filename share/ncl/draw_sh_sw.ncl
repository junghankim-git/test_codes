; Compare snapshots on target day
;
; Contents : Baroclinic Instability Test
; Date     : 20150514
; Author   : Junghan Kim
;
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "/home/jhkim/Study/Library/Shared/NCL/Lib_Math.ncl"
load "/home/jhkim/Study/Library/Shared/NCL/Lib_Plot.ncl"

begin

  onUP      = False

  
  ;varname    = "ps"
  varname    = "wind"

  if (varname .eq. "ps") then
    var_name1  = "ps_dyn"
    var_name2  = "ps"
    var_lname = "Surface Pressure [hPa]"
  else
    var_name1  = "u_dyn"
    var_name2  = "u"
    var_lname = "wind [m/s]"
  end if




  ;;;; File configuration

  if (onUP) then
    filename = "UP-20110725120000-000001.nc"
  else
    filename = "LL-20110725120000-000001.nc"
  end if

  filename1 = "LL-20110725120000-000001.nc"
  filename2 = "LL-000001.nc"

  filedir1 = "/scratch/jhkim/TestBed/vsSW"
  filedir2 = "/scratch/sjchoi/KIAPSGM_GRIMS.v0.9393_wophy_jhkim/2011072512"

  infile1 = filedir1+"/"+filename1
  infile2 = filedir2+"/"+filename2

  print("loading file 1: "+infile1)
  print("loading file 2: "+infile2)

  history_nc_1 = addfile(infile1,"r")
  history_nc_2 = addfile(infile2,"r")



  ;;;; Read Variables
  lats = history_nc_1->lat
  lons = history_nc_1->lon
  nlats = dimsizes(lats)
  nlons = dimsizes(lons)
  print("nlats    = "+tostring(nlats))
  print("nlons    = "+tostring(nlons))
  
; print(lons)
; print(lats(140))
; print(lons(90))

  times = history_nc_1->time
  levs  = history_nc_1->lev
  nlevs_t = dimsizes(times)
  nlevs_v = dimsizes(levs)
  print("nlevs_t  = "+tostring(nlevs_t))
  print("nlevs_v  = "+tostring(nlevs_v))

  ilev_t = 120
  ilev_v = 0
  dlev_t = 4
  do ilev_t = 0, nlevs_t, dlev_t

    if (ilev_t .gt. nlevs_t-1) then
       print("check ilev_t")
       exit()
    end if
  
    if (varname .eq. "ps") then

       var1  = history_nc_1->$var_name1$(ilev_t,:,:)
       var2  = history_nc_2->$var_name2$(ilev_t,:,:)
       var1 = var1/100.0
       var2 = var2/100.0
       var1@units = "hPa"
       var2@units = "hPa"
     
       minvar  =  600.0
       maxvar  = 1040.0
       intvar  =   40.0
       diffminvar = -10.0
       diffmaxvar =  10.0
       diffintvar =   1.0

    else

       var1_u  = history_nc_1->u_dyn(ilev_t,nlevs_v-1,:,:)
       var1_v  = history_nc_1->v_dyn(ilev_t,nlevs_v-1,:,:)
       var2_u  = history_nc_2->u(ilev_t,0,:,:)
       var2_v  = history_nc_2->v(ilev_t,0,:,:)
       var1 = sqrt(var1_u^2+var1_v^2)
       var2 = sqrt(var2_u^2+var2_v^2)
       copy_VarCoords(var1_u, var1)
       copy_VarCoords(var2_u, var2)
;       var1@units = "m/s"
;       var2@units = "m/s"
     
       minvar  = -55.0
       maxvar  =  55.0
       intvar  =  10.0
       diffminvar = -16.5
       diffmaxvar =  16.5
       diffintvar =   2.5

    end if


    loc_x1 = 105.0
    loc_x2 = 145.0
    loc_y1 =  25.0
    loc_y2 =  50.0
    dt = 60.0*36.0
  
    outfilename = DiffPlotLL(varname, "KIM-SH", "KIM-SW", ilev_v, ilev_t, dt, var1, var2, minvar, maxvar, intvar, diffminvar, diffmaxvar, diffintvar, lats, lons, loc_y1, loc_y2, loc_x1, loc_x2, True)


  end do


end
