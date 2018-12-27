
;load "/home/jhkim/Study/Library/Shared/NCL/Lib_LowPlot.ncl"

std_kim_vars =  (/"u", "v", "T", "q", "p", "pint", "ps", "psl", "omega", "ncrain", "crain", "tsfc", "weasd", "canopy", "soilm1", "soilm2", "soilm3", "soilm4", "soilt1", "soilt2", "soilt3", "soilt4", "seaice", "landfrac", "topo", "tg3", "znt", "rh", "hgt", "u10m", "v10m", "t2m", "q2m", "rh2m"/)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; FUNCTION DiffPlotLL
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function DiffPlotLL(varname, org_name, new_name, ilev_v, ilev_t, dt, var1[*][*], var2[*][*], plot_min_in, plot_max_in, plot_int_in, plot_diffmin_in, plot_diffmax_in, plot_diffint_in, lats, lons, local_y1, local_y2, local_x1, local_x2, doPlot)

begin

  diff_var = var2-var1
  copy_VarCoords(var1, diff_var)

  nlats = dimsizes(lats)
  nlons = dimsizes(lons)
  print("# nlats = "+tostring(nlats))
  print("# nlons = "+tostring(nlons))

  nlats_loc = 0
  nlons_loc = 0
  slat_loc = -99
  slon_loc = -99
  elat_loc = -99
  elon_loc = -99

;  print(lats)

  do i = 0, nlats-1
     isSta = False
     if (lats(i) .ge. local_y1 .and. lats(i) .le. local_y2) then
       if (.not. isSta) then
         slat_loc = i
         isSta = True
       end if
       nlats_loc = nlats_loc + 1
     end if
  end do
  do i = 0, nlons-1
     isSta = False
     if (lons(i) .ge. local_x1 .and. lons(i) .le. local_x2) then
       if (.not. isSta) then
         slon_loc = i
         isSta = True
       end if
       nlons_loc = nlons_loc + 1
     end if
  end do
  print("# lat local = "+" "+tostring(slat_loc)+" "+tostring(nlats_loc))
  print("# lon local = "+" "+tostring(slon_loc)+" "+tostring(nlons_loc))
  
  var1_loc = var1(slat_loc:slat_loc+nlats_loc, slon_loc:slon_loc+nlons_loc)
  var2_loc = var2(slat_loc:slat_loc+nlats_loc, slon_loc:slon_loc+nlons_loc)



  Err_l1     = GetLError2D(var1(:,:),var2(:,:),"L1")
  Err_l2     = GetLError2D(var1(:,:),var2(:,:),"L2")
  Err_linf   = GetLError2D(var1(:,:),var2(:,:),"Linf")

  Err_loc_l1     = GetLError2D(var1_loc(:,:),var2_loc(:,:),"L1")
  Err_loc_l2     = GetLError2D(var1_loc(:,:),var2_loc(:,:),"L2")
  Err_loc_linf   = GetLError2D(var1_loc(:,:),var2_loc(:,:),"Linf")


  print(" ")
  print("- Rel. L1   = "+Err_l1)
  print("- Rel. L2   = "+Err_l2)
  print("- Rel. Linf = "+Err_linf)
  print(" ")
  print(" ")
  print("- Rel. L1 loc   = "+Err_loc_l1)
  print("- Rel. L2 loc   = "+Err_loc_l2)
  print("- Rel. Linf loc = "+Err_loc_linf)
  print(" ")
  print(" ")
  print("- var1 for      = "+var1(140,90))
  print("- var2 for      = "+var2(140,90))
  print(" ")

  outdataname = "data/"+varname+".txt"
  if (ilev_t .eq. 0) then
;    write_table(outdataname, "w", [/ilev_v, ilev_t, Err_l2, Err_loc_l2, var1(140, 90), var2(140, 90)/], "%03i %03i %6.3e %6.3e %6.3e %6.3e")
  else
;    write_table(outdataname, "a", [/ilev_v, ilev_t, Err_l2, Err_loc_l2, var1(140, 90), var2(140, 90)/], "%03i %03i %6.3e %6.3e %6.3e %6.3e")
  end if




  ;;;;  Plot

  if(doPlot) then
    format = "%6.3e"
    tx_err_l1     = sprintf("L1   = "+format,Err_l1)
    tx_err_l2     = sprintf("L2   = "+format,Err_l2)
    tx_err_linf   = sprintf("Linf = "+format,Err_linf)


    wks_type = "ps"

    if(wks_type .eq."x11".or.wks_type.eq."png") then 
       wks_type@wkWidth = 1000 
       wks_type@wkHeight = 1000
    end if
    if(wks_type.eq."ps".or.wks_type.eq."eps".or.wks_type.eq."epsi".or.wks_type.eq."pdf") then  
       ;wks_type@wkOrientation = "landscape" 
       wks_type@wkPaperSize = "A3" 
    end if


    outfilename = "Figs/"+varname+"_"+org_name+"_"+new_name+"_"+sprinti("%02d",ilev_v)+"_"+sprinti("%04d",ilev_t)
    print("outfilename = "+outfilename)
    wks  = gsn_open_wks(wks_type,outfilename)
;    gsn_merge_colormaps(wks,"BlGrYeOrReVi200","hotcold_18lev");
;    gsn_merge_colormaps(wks,"rainbow","gsdtol");
    gsn_define_colormap(wks,"rainbow")


;    gsn_draw_colormap(wks)
    res                        = True                   ; plot mods desired
    res@gsnDraw                = False
    res@gsnFrame               = False

    ; for unique point grid
    if (isUP) then
      res@sfXArray              = lon;
      res@sfYArray              = lat;
    end if
    res@gsnSpreadColors = True;
    res@mpCenterLonF = 202.0;

    res@cnLevelSelectionMode = "ManualLevels";

;    res@gsnSpreadColorStart = 2;
;    res@gsnSpreadColorEnd = 201;

    ;res@tiYAxisString = "Latitude"
    ;res@tiXAxisString = "Longitude"

    res@cnFillOn = True;     
    res@cnInfoLabelOn = False
    res@lbLabelBarOn  = True
    res@lbLabelStride = 2



    if (.not. isUP) then
       plot_var1 = var1(:,:)
       plot_var2 = var2(:,:)
       plot_diff_var = diff_var(:,:)
    else
       plot_var1 = var1(:)
       plot_var2 = var2(:)
       plot_diff_var = diff_var(:)
    end if


    ; Time Level
    plot_min = plot_min_in
    plot_max = plot_max_in
    plot_int = plot_int_in
    plot_diffmin = plot_diffmin_in
    plot_diffmax = plot_diffmax_in
    plot_diffint = plot_diffint_in


    xpos = new(4, double)
    ypos = new(3, double)

    x_margin = 0.03
    y_margin = -0.04

    width_g  = 0.22
    width_l  = 0.165 ;0.16 was small, 0.17 was large

    xpos(0)  = 0.05
    xpos(1)  = xpos(0) + width_g + x_margin
    xpos(2)  = xpos(1) + width_g + x_margin
    xpos(3)  = xpos(2) + width_g + x_margin

    ypos(0)  = 0.95
    ypos(1)  = ypos(0) - width_g - y_margin
    ypos(2)  = ypos(1) - width_g - y_margin
    adjust_l = 0.0275 ; 0.027 more down


    ;local_x1 = 105.0
    ;local_x2 = 150.0
    ;local_y1 =  20.0
    ;local_y2 =  50.0
    ;local_x1 = 105.0
    ;local_x2 = 145.0
    ;local_y1 =  25.0
    ;local_y2 =  50.0


    ; ######################## Box and Line ########################
    xndc_beg = new(1,double)
    yndc_beg = new(1,double)
    xndc_end = new(1,double)
    yndc_end = new(1,double)

    x_box = (/local_x1, local_x1, local_x2, local_x2, local_x1/)
    y_box = (/local_y1, local_y2, local_y2, local_y1, local_y1/)

    bres             = True
    ;bres@gsFillColor = "black"

    lnres                  = True
    lnres@gsLineThicknessF = 1.5  ; default is 1.5
    lnres@gsLineColor      = "black"


  
    plots = new (6, graphic)

    ; [[[[[[[[[[[[[[[[[[[[[[[[      Ps       ]]]]]]]]]]]]]]]]]]]]]]]]
    ; ######################## Global ########################
    ; draw contour line
    res@cnLinesOn      = False
    res@cnLineLabelsOn = False

    res@mpMinLonF =   0.0
    res@mpMaxLonF = 360.0
    res@mpMinLatF = -90.0
    res@mpMaxLatF =  90.0

    ; for Ps
    res@cnMaxLevelValF  = plot_max
    res@cnMinLevelValF  = plot_min
    res@cnLevelSpacingF = plot_int
;    res@cnLevelCount    = "5"

    res@vpXF = xpos(0)
    res@vpYF = ypos(0)
    res@vpHeightF = width_g
    res@vpWidthF  = width_g
    res@tiMainString = org_name+"("+sprintf("%05.2f",dt*ilev_t/60./60.)+" hour)"
    plots(0) = gsn_csm_contour_map_ce(wks,plot_var1,res)
    tmp = gsn_add_polyline(wks,plots(0),x_box,y_box,bres)
    draw(plots(0))

    res@vpXF = xpos(0)
    res@vpYF = ypos(1)
    res@vpHeightF = width_g
    res@vpWidthF  = width_g
    res@tiMainString = new_name+"("+sprintf("%05.2f",dt*ilev_t/60./60.)+" hour)"
    plots(1) = gsn_csm_contour_map_ce(wks,plot_var2,res)
    tmp = gsn_add_polyline(wks,plots(1),x_box,y_box,bres)
    draw(plots(1))

    res@cnMaxLevelValF  = plot_diffmax
    res@cnMinLevelValF  = plot_diffmin
    res@cnLevelSpacingF = plot_diffint

    res@vpXF = xpos(0)
    res@vpYF = ypos(2)
    res@vpHeightF = width_g
    res@vpWidthF  = width_g
    res@tiMainString = "Difference"
    plots(2) = gsn_csm_contour_map_ce(wks,plot_diff_var,res)
    tmp = gsn_add_polyline(wks,plots(2),x_box,y_box,bres)

    ; error text
    pos_lon  = 30
    pos_lat1 = 80
    pos_lat2 = 70
    pos_lat3 = 60
    txres = True
    txres@txFont                = "helvetica-bold"
    txres@txFontQuality         = "High"
    txres@txJust                = "CenterLeft"
    txres@txFontHeightF         = 0.004
    txres@txBackgroundFillColor = 0 ;default -1
    txres@txFontColor           = "black"
    ;txres@txFontColor           = "red"
    txid = gsn_add_text(wks, plots(2), tx_err_l1,     pos_lon, pos_lat1, txres)
    txid = gsn_add_text(wks, plots(2), tx_err_l2,     pos_lon, pos_lat2, txres)
    txid = gsn_add_text(wks, plots(2), tx_err_linf,   pos_lon, pos_lat3, txres)

    draw(plots(2))


    ; ########################  Local ########################
    ; draw contour line
    if (varname .eq. "ps") then
     res@cnLinesOn      = True
     res@cnLineLabelsOn = True
    end if

    ; Local Range
    res@mpMinLonF = local_x1
    res@mpMaxLonF = local_x2
    res@mpMinLatF = local_y1
    res@mpMaxLatF = local_y2

     ; Time Level
     ; Surface Pressure Range
     ;plot_int =  local_y1
     ;plot_min =  700.0;-plot_int
     ;plot_max = 1100.0;+plot_int
     ;plot_diffmin = -1
     ;plot_diffmax =  1
     ;plot_diffint = 0.1

    ; for Ps
    res@cnMaxLevelValF  = plot_max
    res@cnMinLevelValF  = plot_min
    res@cnLevelSpacingF = plot_int
;    res@cnLevelCount    = "5"

    res@vpXF = xpos(1)
    res@vpYF = ypos(0)-adjust_l
    res@vpHeightF = width_l
    res@vpWidthF  = width_l
    res@tiMainString = org_name+"("+sprintf("%05.2f",dt*ilev_t/60./60.)+" hour)"
    plots(3) = gsn_csm_contour_map_ce(wks,plot_var1,res)
    draw(plots(3))
    res@vpXF = xpos(1)
    res@vpYF = ypos(1)-adjust_l
    res@vpHeightF = width_l
    res@vpWidthF  = width_l
    res@tiMainString = new_name+"("+sprintf("%05.2f",dt*ilev_t/60./60.)+" hour)"
    plots(4) = gsn_csm_contour_map_ce(wks,plot_var2,res)
    draw(plots(4))

    res@cnMaxLevelValF  = plot_diffmax
    res@cnMinLevelValF  = plot_diffmin
    res@cnLevelSpacingF = plot_diffint

    ;res@cnInfoLabelOn = False
    res@cnLinesOn      = False
    res@cnLineLabelsOn = False
    res@vpXF = xpos(1)
    res@vpYF = ypos(2)-adjust_l
    res@vpHeightF = width_l
    res@vpWidthF  = width_l
    res@tiMainString = "Difference"
    plots(5) = gsn_csm_contour_map_ce(wks,plot_diff_var,res)
    draw(plots(5))


    ; Draw line

    datatondc(plots(0), local_x2, local_y1, xndc_beg, yndc_beg)
    datatondc(plots(3), local_x1, local_y1, xndc_end, yndc_end)
    gsn_polyline_ndc(wks, (/xndc_beg,xndc_end/),(/yndc_beg,yndc_end/), lnres)
    datatondc(plots(0), local_x2, local_y2, xndc_beg, yndc_beg)
    datatondc(plots(3), local_x1, local_y2, xndc_end, yndc_end)
    gsn_polyline_ndc(wks, (/xndc_beg,xndc_end/),(/yndc_beg,yndc_end/), lnres)

    datatondc(plots(1), local_x2, local_y1, xndc_beg, yndc_beg)
    datatondc(plots(4), local_x1, local_y1, xndc_end, yndc_end)
    gsn_polyline_ndc(wks, (/xndc_beg,xndc_end/),(/yndc_beg,yndc_end/), lnres)
    datatondc(plots(1), local_x2, local_y2, xndc_beg, yndc_beg)
    datatondc(plots(4), local_x1, local_y2, xndc_end, yndc_end)
    gsn_polyline_ndc(wks, (/xndc_beg,xndc_end/),(/yndc_beg,yndc_end/), lnres)

    datatondc(plots(2), local_x2, local_y1, xndc_beg, yndc_beg)
    datatondc(plots(5), local_x1, local_y1, xndc_end, yndc_end)
    gsn_polyline_ndc(wks, (/xndc_beg,xndc_end/),(/yndc_beg,yndc_end/), lnres)
    datatondc(plots(2), local_x2, local_y2, xndc_beg, yndc_beg)
    datatondc(plots(5), local_x1, local_y2, xndc_end, yndc_end)
    gsn_polyline_ndc(wks, (/xndc_beg,xndc_end/),(/yndc_beg,yndc_end/), lnres)


    frame(wks)


    onSys = True
    if (onSys) then
       system("mkdir -p ./Figs")
       psname   = outfilename+".ps"
       pngname  = outfilename+".png"
       pngsname = outfilename+"_small.png"

       if (wks_type .eq. "ps") then
          system("ps2png "+psname+" -o "+pngname)
          system("ps2png -s 50 "+psname+" -o "+pngsname)
       else
          system("convert -s 50 "+pngname+" "+pngsname)
       end if
    end if

    ;pres = True
    ;pres@gsnPanelRowSpec = True
    ;gsn_panel(wks,plots,(/3,2/),pres)



  end if ; if (True)


  return(outfilename)

end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; END FUNCTION DiffPlotLL
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;








;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; FUNCTION PlotStdOutKIM
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function PlotStdOutKIM(filename_in, filepath_in, varnames_in, reverse_in, itime_in, ilev_in, isLocal_in)

begin

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Filenames and variables
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  if (filename_in .eq. "") then
    filename   = "UP-20110725220000-000011.nc"
  else
    filename   = filename_in
  end if

  if (filepath_in .eq. "") then
    filepath  = "/scratch/jhkim/TestBed/KIM/Output/0.24.02/ne30/gnu/10hours/SH_KPPACK/403"
  else
    filepath  = filepath_in
  end if


  if (varnames_in(0) .eq. "") then
    varnames   = std_kim_vars
  else
    varnames   = varnames_in
  end if

  if (isLocal_in) then
    varnames   = (/"pint", "rain"/)
  end if

  reverse      = reverse_in

  print("### Step0: Input arguments... : ")
  print(" ")



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Loading input file and Coordinate
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  infilename   = filepath+"/"+filename
  infile       = addfile(infilename,"r")

  print("### Step1: Loading input file... : "+infilename)
  print(" ")

  lats   = infile->lat
  lons   = infile->lon
  times  = infile->time
  levs   = infile->lev
  ilevs  = infile->ilev
  hyam   = infile->hyam
  hybm   = infile->hybm
  hyai   = infile->hyai
  hybi   = infile->hybi

  nlats  = dimsizes(lats)
  nlons  = dimsizes(lons)
  ntimes = dimsizes(times)
  nlevs  = dimsizes(levs)


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Levels
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  nvars        = dimsizes(varnames)
  itime        = itime_in
  if (reverse) then
    ilev       = nlevs-1-ilev_in
  else
    ilev       = ilev_in
  end if

  isUP         = True


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Generate Positions
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print("### Step2 - START : Generate Position... : ")
  ; Positions
  ; - input
  nplots  = nvars
  ; - outputs
  nx      = 0
  ny      = 0
  xpos    = new(nplots, double) 
  ypos    = new(nplots, double) 
  icanvas = new(nplots, integer) 
  xwidth  = 0.0d0
  ywidth  = 0.0d0
  ncanvas = 0
  nmod    = 0
  err = GetCanvasPosition(nplots, nx, ny, xpos, ypos, xwidth, ywidth, ncanvas, nmod, icanvas)

  print("### Step2 -  END  : Generate Position... : ")
  print(" ")



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Make Canvas
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  wks_type = "ps"
  if(wks_type .eq."x11".or.wks_type.eq."png") then
     wks_type@wkWidth  = 1000
     wks_type@wkHeight = 1000
  end if
  if(wks_type.eq."ps".or.wks_type.eq."eps".or.wks_type.eq."epsi".or.wks_type.eq."pdf") then
     ;wks_type@wkOrientation = "landscape"
     wks_type@wkPaperSize    = "A4"
  end if

  outfilename = new(ncanvas, string)
  do icanv = 0, ncanvas-1
    outfilename(icanv) = "output"+tostring(icanv)
  end do

;  wks         = new(ncanvas, graphic)
;  do icanv = 0, ncanvas-1
;    outfilename(icanv) = "output"+tostring(icanv)
;    wks(icanv)         = gsn_open_wks(wks_type,outfilename(icanv))
;    gsn_define_colormap(wks(icanv), "rainbow")
;    print(wks(icanv))
;  end do
;    wks(0)         = gsn_open_wks(wks_type,outfilename(0))
;    print(wks)

  
 

  print("### Step3 - START : Variable Loop... : ")
  icanv = -1
  do ivar = 0, nvars-1

     iplot = ivar
     if (icanvas(iplot) .gt. icanv) then
        if (icanv .ne. -1) then
           frame(wks)
        end if
        icanv  = icanv + 1
        wks    = gsn_open_wks(wks_type,outfilename(icanvas(iplot)))
     end if
     print(" - icanvas = "+tostring(icanv)+", for ivar = "+tostring(ivar))

     print(" - Reading variables... : "+varnames(ivar))
  
     vardim = dimsizes(infile->$varnames(ivar)$)
     ndims  = dimsizes(vardim)
   
     if (ndims .eq. 1) then
       var   = infile->$varnames(ivar)$(:)
     else if (ndims .eq. 2) then
       var   = infile->$varnames(ivar)$(itime, :)
     else if (ndims .eq. 3) then
       var   = infile->$varnames(ivar)$(itime, ilev, :)
     end if
     end if
     end if
   



   ; plot range

     range = False
 
   
     ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     ;; Plotting
     ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     res = PlotSingle(wks, xpos(iplot), ypos(iplot), xwidth, ywidth, varnames(iplot), var, range, lats, lons, isUP, isLocal_in, False)
     print(" - Plotting... : ")
     print(" ")

  
     delete(vardim) ; for reassignment

  end do
  frame(wks)

  print("### Step3 -   END : Variable Loop... : ")
  print(" ")

  delete(infile)

  return(0)

end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; END FUNCTION PlotStdOutKIM
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;






;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; FUNCTION PlotDiffStdOutKIM
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function PlotDiffStdOutKIM(filename1_in, filename2_in, filepath1_in, filepath2_in, varnames_in, reverse1_in, reverse2_in, itime_in, ilev_in, isPlot_in, isLocal_in)

begin

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Filenames and variables
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  if (filename1_in .eq. "") then
    filename1   = "UP-20110725220000-000011.nc"
  else
    filename1   = filename1_in
  end if
  if (filename2_in .eq. "") then
    filename2   = "UP-20110725220000-000011.nc"
  else
    filename2   = filename2_in
  end if

  if (filepath1_in .eq. "") then
    filepath1  = "/scratch/jhkim/TestBed/KIM/Output/0.24.02/ne30/gnu/10hours/SH_KPPACK/403"
  else
    filepath1  = filepath1_in
  end if

  if (filepath2_in .eq. "") then
    filepath2  = "/scratch/jhkim/TestBed/Data/history"
  else
    filepath2  = filepath2_in
  end if

  if (varnames_in(0) .eq. "") then
    varnames   = std_kim_vars
  else
    varnames   = varnames_in
  end if
;print(varnames)


  reverse1     = reverse1_in
  reverse2     = reverse2_in
  isPlot       = isPlot_in

  print("### Step0: Input arguments... : ")
  print(" ")



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Loading input file and Coordinate
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  infilename1   = filepath1+"/"+filename1
  infilename2   = filepath2+"/"+filename2
  infile1      = addfile(infilename1,"r")
  infile2      = addfile(infilename2,"r")

  print("### Step1: Loading input file... : "+infilename1)
  print("### Step1: Loading input file... : "+infilename2)
  print(" ")

  lats   = infile1->lat
  lons   = infile1->lon
  times  = infile1->time
  levs   = infile1->lev
  ilevs  = infile1->ilev
  hyam   = infile1->hyam
  hybm   = infile1->hybm
  hyai   = infile1->hyai
  hybi   = infile1->hybi

  nlats  = dimsizes(lats)
  nlons  = dimsizes(lons)
  ntimes = dimsizes(times)
  nlevs  = dimsizes(levs)


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Levels
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  nvars        = dimsizes(varnames)
  itime        = itime_in
  if (reverse1) then
    ilev1      = nlevs-1-ilev_in
  else
    ilev1      = ilev_in
  end if

  if (reverse2) then
    ilev2      = nlevs-1-ilev_in
  else
    ilev2      = ilev_in
  end if

  isUP         = True


  if (isPlot) then
   
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;; Generate Positions
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   
     print("### Step2 - START : Generate Position... : ")
     ; Positions
     ; - input
     nplots  = 3*nvars
     ; - outputs
     nx      = 0
     ny      = 0
     xpos    = new(nplots, double) 
     ypos    = new(nplots, double) 
     icanvas = new(nplots, integer) 
     xwidth  = 0.0d0
     ywidth  = 0.0d0
     ncanvas = 0
     nmod    = 0
     err = GetCanvasPosition(nplots, nx, ny, xpos, ypos, xwidth, ywidth, ncanvas, nmod, icanvas)
   
     print("### Step2 -  END  : Generate Position... : ")
     print(" ")
   
   
   
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;; Make Canvas
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     wks_type = "ps"
     if(wks_type .eq."x11".or.wks_type.eq."png") then
        wks_type@wkWidth  = 1000
        wks_type@wkHeight = 1000
     end if
     if(wks_type.eq."ps".or.wks_type.eq."eps".or.wks_type.eq."epsi".or.wks_type.eq."pdf") then
        ;wks_type@wkOrientation = "landscape"
        wks_type@wkPaperSize    = "A4"
     end if
   
     outfilename = new(ncanvas, string)
     do icanv = 0, ncanvas-1
       outfilename(icanv) = "output"+tostring(icanv)
     end do
   
   ;  wks         = new(ncanvas, graphic)
   ;  do icanv = 0, ncanvas-1
   ;    outfilename(icanv) = "output"+tostring(icanv)
   ;    wks(icanv)         = gsn_open_wks(wks_type,outfilename(icanv))
   ;    gsn_define_colormap(wks(icanv), "rainbow")
   ;    print(wks(icanv))
   ;  end do
   ;    wks(0)         = gsn_open_wks(wks_type,outfilename(0))
   ;    print(wks)

  end if ; isPlot
  
 

  print("### Step3 - START : Variable Loop... : ")


  icanv = -1
  do ivar = 0, nvars-1

     if (isPlot) then
        iplot = nx*ivar
        if (icanvas(iplot) .gt. icanv) then
           if (icanv .ne. -1) then
              frame(wks)
           end if
           icanv  = icanv + 1
           wks    = gsn_open_wks(wks_type,outfilename(icanvas(iplot)))
        end if
        print(" - icanvas = "+tostring(icanv)+", for ivar = "+tostring(ivar))
     end if

     print(" - Reading variables... : "+varnames(ivar))
  
     vardim1 = dimsizes(infile1->$varnames(ivar)$)
     ndims1  = dimsizes(vardim1)
     vardim2 = dimsizes(infile2->$varnames(ivar)$)
     ndims2  = dimsizes(vardim2)
   
     if (ndims1 .eq. 1) then
       var1   = infile1->$varnames(ivar)$(:)
     else if (ndims1 .eq. 2) then
       var1   = infile1->$varnames(ivar)$(itime, :)
     else if (ndims1 .eq. 3) then
       var1   = infile1->$varnames(ivar)$(itime, ilev1, :)
     end if
     end if
     end if

     if (ndims2 .eq. 1) then
       var2   = infile2->$varnames(ivar)$(:)
     else if (ndims2 .eq. 2) then
       var2   = infile2->$varnames(ivar)$(itime, :)
     else if (ndims2 .eq. 3) then
       var2   = infile2->$varnames(ivar)$(itime, ilev2, :)
     end if
     end if
     end if
   
   
     ;-----------------------
     ; Statistic
     ;-----------------------
     Err_l1     = GetLError(var1(:),var2(:),"L1")
     Err_l2     = GetLError(var1(:),var2(:),"L2")
     Err_linf   = GetLError(var1(:),var2(:),"Linf")
     print(" - Rel. L1   = "+Err_l1)
     print(" - Rel. L2   = "+Err_l2)
     print(" - Rel. Linf = "+Err_linf)

     ;-----------------------
     ; Difference
     ;-----------------------
     diff_var = var2 - var1
     copy_VarCoords(var1, diff_var)
     
     if (isPlot) then
       diff_var@long_name = "L2 = "+tostring(Err_l2)
       if (isdefined("diff_var@units")) then
          diff_var@units     = var1@units
       else
          diff_var@units     = "unknown"
       end if
     end if


     if (isPlot) then
   
      ; plot range
        minvar = min(var1)
        maxvar = max(var1)
        minp = (toint(minvar/10.0d0))*10.0d0
        maxp = (toint(maxvar/10.0d0))*10.0d0
        if (minp .eq. 0.0d0) then
          if (minvar .lt. 0.0d0) then
            minp = -1.0d0
          else
            minp =  0.0d0
          end if
        end if
        if (maxp .eq. 0.0d0) then
          maxp = 1.0d0
        end if
        intp = (maxp-minp)/10.0d0
    
      
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        ;; Plotting
        ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        range = True
        range@min = minp
        range@max = maxp
        range@int = intp

        fraction = 1.0d0/20.0d0
        diff_range = True
        diff_range@min = -(maxp-minp)*fraction
        diff_range@max =  (maxp-minp)*fraction
        diff_range@int =  (maxp-minp)*fraction/5.0d0

        res = PlotSingle(wks, xpos(iplot+0), ypos(iplot+0), xwidth, ywidth, varnames(ivar)+"-1", var1, range, lats, lons, isUP, isLocal_in, False)
        res = PlotSingle(wks, xpos(iplot+1), ypos(iplot+1), xwidth, ywidth, varnames(ivar)+"-2", var2, range, lats, lons, isUP, isLocal_in, False)
        res = PlotSingle(wks, xpos(iplot+2), ypos(iplot+2), xwidth , ywidth, "difference", diff_var, diff_range, lats, lons, isUP, isLocal_in, False)
      

        print(" - Plotting... : ")
        print(" ")

     end if
  
     delete(vardim1) ; for reassignment
     delete(vardim2) ; for reassignment

  end do

  if (isPlot) then
     frame(wks)
  end if

  print("### Step3 -   END : Variable Loop... : ")
  print(" ")

  delete(infile1)
  delete(infile2)

  return(0)

end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; FUNCTION PlotDiffStdOutKIM
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;






;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; FUNCTION PlotDiffRainStdOutKIM
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function PlotDiffRainStdOutKIM(filename1_in, filename2_in, filepath1_in, filepath2_in, reverse1_in, reverse2_in, itime_in, ilev_in, isPlot_in)

begin

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Filenames and variables
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  filename1 = filename1_in
  filename2 = filename2_in
  nfiles1 = dimsizes(filename1)
  nfiles2 = dimsizes(filename2)


  if (filepath1_in .eq. "") then
    filepath1  = "/scratch/jhkim/TestBed/KIM/Output/0.24.02/ne30/gnu/10hours/SH_KPPACK/403"
  else
    filepath1  = filepath1_in
  end if

  if (filepath2_in .eq. "") then
    filepath2  = "/scratch/jhkim/TestBed/Data/history"
  else
    filepath2  = filepath2_in
  end if

  reverse1     = reverse1_in
  reverse2     = reverse2_in
  isPlot       = isPlot_in

  print("### Step0: Input arguments... : ")
  print(" ")



  ; [[ Local ]]
  local_x1 = 105.0
  local_x2 = 150.0
  local_y1 =  20.0
  local_y2 =  50.0


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Levels
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  nvars        = 2
  itime        = itime_in
  if (reverse1) then
    ilev1      = nlevs-1-ilev_in
  else
    ilev1      = ilev_in
  end if

  if (reverse2) then
    ilev2      = nlevs-1-ilev_in
  else
    ilev2      = ilev_in
  end if

  isUP         = True


   
  print("### Step3 - START : Variable Loop... : ")

  varname = (/"ps", "prec."/)

  do ivar = 0, 1

     if (ivar .eq. 0) then
       nvart    = 1
       varnamet = (/"pint"/)
     else
       nvart    = 2
       varnamet = (/"crain","ncrain"/)
     end if

     print(" - Reading variables... : "+varnamet)

     infile1 = addfile(filepath1+"/"+filename1(0), "r")
     infile2 = addfile(filepath2+"/"+filename2(0), "r")
     lats   = infile1->lat
     lons   = infile1->lon

     vardim1 = dimsizes(infile1->$varnamet(0)$)
     ndims1  = dimsizes(vardim1)
     vardim2 = dimsizes(infile2->$varnamet(0)$)
     ndims2  = dimsizes(vardim2)

     ; ps
     if (ndims1 .eq. 1) then

       var1 = infile1->$varnamet(0)$(:)
       do k = 1, nvart-1
          var1(:) = var1(:) + infile1->$varnamet(k)$(:)
       end do
       do i = 1, nfiles1-1
         infile1 = addfile(filepath1+"/"+filename1(i), "r")
       do k = 0, nvart-1
          var1(:) = var1(:) + infile1->$varnamet(k)$(:)
       end do
       end do

     else if (ndims1 .eq. 2) then

       var1 = infile1->$varnamet(0)$(itime,:)
       do k = 1, nvart-1
          var1(:) = var1(:) + infile1->$varnamet(k)$(itime,:)
       end do
       do i = 1, nfiles1-1
         infile1 = addfile(filepath1+"/"+filename1(i), "r")
       do k = 0, nvart-1
          var1(:) = var1(:) + infile1->$varnamet(k)$(itime,:)
       end do
       end do

     else if (ndims1 .eq. 3) then

       var1 = infile1->$varnamet(0)$(itime,ilev1,:)
       do k = 1, nvart-1
          var1(:) = var1(:) + infile1->$varnamet(k)$(itime,ilev1,:)
       end do
       do i = 1, nfiles1-1
         infile1 = addfile(filepath1+"/"+filename1(i), "r")
       do k = 0, nvart-1
          var1(:) = var1(:) + infile1->$varnamet(k)$(itime,ilev1,:)
       end do
       end do

     end if
     end if
     end if
     if (ivar .eq. 0) then
        var1 = var1/todouble(nfiles1)
     end if

     if (ndims2 .eq. 1) then

       var2 = infile2->$varnamet(0)$(:)
       do k = 1, nvart-1
          var2(:) = var2(:) + infile2->$varnamet(k)$(:)
       end do
       do i = 1, nfiles2-1
         infile2 = addfile(filepath2+"/"+filename2(i), "r")
       do k = 0, nvart-1
          var2(:) = var2(:) + infile2->$varnamet(k)$(:)
       end do
       end do

     else if (ndims2 .eq. 2) then

       var2 = infile2->$varnamet(0)$(itime, :)
       do k = 1, nvart-1
          var2(:) = var2(:) + infile2->$varnamet(k)$(itime,:)
       end do
       do i = 1, nfiles2-1
         infile2 = addfile(filepath2+"/"+filename2(i), "r")
       do k = 0, nvart-1
          var2(:) = var2(:) + infile2->$varnamet(k)$(itime,:)
       end do
       end do

     else if (ndims2 .eq. 3) then
       var2 = infile2->$varnamet(0)$(itime,ilev2, :)
       do k = 1, nvart-1
          var2(:) = var2(:) + infile2->$varnamet(k)$(itime,ilev2,:)
       end do
       do i = 1, nfiles2-1
         infile2 = addfile(filepath2+"/"+filename2(i), "r")
       do k = 0, nvart-1
          var2(:) = var2(:) + infile2->$varnamet(k)$(itime,ilev2,:)
       end do
       end do
     end if
     end if
     end if
     if (ivar .eq. 0) then
        var2 = var2/todouble(nfiles2)
     end if

     if (ivar .eq. 0) then
        ps_1 = var1
        ps_2 = var2
     else
        prec_1 = var1
        prec_2 = var2
     end if

     delete(vardim1) ; for reassignment
     delete(vardim2) ; for reassignment
     delete(varnamet)



  end do



  ps_1 = ps_1/100.0d0
  ps_2 = ps_2/100.0d0
  ps_1@units = "hPa"
  ps_2@units = "hPa"
  ps_1@long_name = "Surface pressure"
  ps_2@long_name = "Surface pressure"
  
  diff_ps = ps_2-ps_1
  copy_VarCoords(ps_1, diff_ps)

  prec_2@long_name = "Precipitation"
  diff_prec = prec_2-prec_1
  copy_VarCoords(prec_2, diff_prec)


  min_ps_1    = min(ps_1)
  max_ps_1    = max(ps_1)
  min_ps_2    = min(ps_2)
  max_ps_2    = max(ps_2)
  min_diff_ps = min(diff_ps)
  max_diff_ps = max(diff_ps)

  min_prec_1    = min(prec_1)
  max_prec_1    = max(prec_1)
  min_prec_2    = min(prec_2)
  max_prec_2    = max(prec_2)
  min_diff_prec = min(diff_prec)
  max_diff_prec = max(diff_prec)


  l_ps_1      = CutLocal(ps_1, lons, lats, local_x1, local_x2, local_y1, local_y2)
  l_ps_2      = CutLocal(ps_2, lons, lats, local_x1, local_x2, local_y1, local_y2)
  diff_l_ps   = CutLocal(diff_ps, lons, lats, local_x1, local_x2, local_y1, local_y2)
  l_prec_1    = CutLocal(prec_1, lons, lats, local_x1, local_x2, local_y1, local_y2)
  l_prec_2    = CutLocal(prec_2, lons, lats, local_x1, local_x2, local_y1, local_y2)
  diff_l_prec = CutLocal(diff_prec, lons, lats, local_x1, local_x2, local_y1, local_y2)

  min_l_ps_1    = min(l_ps_1)
  max_l_ps_1    = max(l_ps_1)
  min_l_ps_2    = min(l_ps_2)
  max_l_ps_2    = max(l_ps_2)
  min_diff_l_ps = min(diff_l_ps)
  max_diff_l_ps = max(diff_l_ps)

  min_l_prec_1    = min(l_prec_1)
  max_l_prec_1    = max(l_prec_1)
  min_l_prec_2    = min(l_prec_2)
  max_l_prec_2    = max(l_prec_2)
  min_diff_l_prec = min(diff_l_prec)
  max_diff_l_prec = max(diff_l_prec)

  

  ;print(ps_1)
  if (isUP) then
     Pserr_l1     = GetLError(ps_1(:),ps_2(:),"L1")
     Pserr_l2     = GetLError(ps_1(:),ps_2(:),"L2")
     Pserr_linf   = GetLError(ps_1(:),ps_2(:),"Linf")
     precerr_l1   = GetLError(prec_1(:),prec_2(:),"L1")
     precerr_l2   = GetLError(prec_1(:),prec_2(:),"L2")
     precerr_linf = GetLError(prec_1(:),prec_2(:),"Linf")

     l_Pserr_l2   = GetLError(l_ps_1(:),l_ps_2(:),"L2")
     l_precerr_l2 = GetLError(l_prec_1(:),l_prec_2(:),"L2")
  else
     Pserr_l1     = GetLError2D(ps_1(:,:),ps_2(:,:),"L1")
     Pserr_l2     = GetLError2D(ps_1(:,:),ps_2(:,:),"L2")
     Pserr_linf   = GetLError2D(ps_1(:,:),ps_2(:,:),"Linf")
     precerr_l1   = GetLError2D(prec_1(:,:),prec_2(:,:),"L1")
     precerr_l2   = GetLError2D(prec_1(:,:),prec_2(:,:),"L2")
     precerr_linf = GetLError2D(prec_1(:,:),prec_2(:,:),"Linf")

     l_Pserr_l2   = GetLError(l_ps_1(:,:),l_ps_2(:,:),"L2")
     l_precerr_l2 = GetLError(l_prec_1(:,:),l_prec_2(:,:),"L2")
  end if

  print(" ")
  print("- Ps   L1   = "+Pserr_l1)
  print("- Ps   L2   = "+Pserr_l2)
  print("- Ps   Linf = "+Pserr_linf)
  print("- Prec L1   = "+precerr_l1)
  print("- Prec L2   = "+precerr_l2)
  print("- Prec Linf = "+precerr_linf)
  print(" ")
  print("- Ps   L2 (local)  = "+l_Pserr_l2)
  print("- Prec L2 (local)  = "+l_precerr_l2)
  print(" ")

  format = "%6.3e"
  tx_pserr_l1     = sprintf("L1   = "+format,Pserr_l1)
  tx_pserr_l2     = sprintf("L2   = "+format,Pserr_l2)
  tx_pserr_linf   = sprintf("Linf = "+format,Pserr_linf)
  tx_precerr_l1   = sprintf("L1   = "+format,precerr_l1)
  tx_precerr_l2   = sprintf("L2   = "+format,precerr_l2)
  tx_precerr_linf = sprintf("Linf = "+format,precerr_linf)

  tx_l_pserr_l2   = sprintf("L2   = "+format,l_Pserr_l2)
  tx_l_precerr_l2 = sprintf("L2   = "+format,l_precerr_l2)



  ;;;;  Plot

  if(True) then

    wks_type = "ps"

    if(wks_type .eq."x11".or.wks_type.eq."png") then 
       wks_type@wkWidth = 1000 
       wks_type@wkHeight = 1000
    end if
    if(wks_type.eq."ps".or.wks_type.eq."eps".or.wks_type.eq."epsi".or.wks_type.eq."pdf") then  
       ;wks_type@wkOrientation = "landscape" 
       wks_type@wkPaperSize = "A3" 
    end if


    wks  = gsn_open_wks(wks_type,"output")
;    gsn_merge_colormaps(wks,"BlGrYeOrReVi200","hotcold_18lev");
;    gsn_merge_colormaps(wks,"rainbow","gsdtol");
    gsn_define_colormap(wks,"rainbow")


;    gsn_draw_colormap(wks)
    res                        = True                   ; plot mods desired
    res@gsnDraw                = False
    res@gsnFrame               = False

    ; for unique point grid
    if (isUP) then
      res@sfXArray              = lons;
      res@sfYArray              = lats;
    end if
    res@gsnSpreadColors = True;
    res@mpCenterLonF = 202.0d0;

    res@cnLevelSelectionMode = "ManualLevels";

;    res@gsnSpreadColorStart = 2;
;    res@gsnSpreadColorEnd = 201;

    ;res@tiYAxisString = "Latitude"
    ;res@tiXAxisString = "Longitude"

    res@cnFillOn = True;     
    res@cnInfoLabelOn = False
    res@lbLabelBarOn  = True
    res@lbLabelStride = 2



    if (.not. isUP) then
       plot_ps1 = ps_1(:,:)
       plot_ps2 = ps_2(:,:)
       plot_diff_ps = diff_ps(:,:)

       plot_prec1 = prec_1(:,:)
       plot_prec2 = prec_2(:,:)
       plot_diff_prec = diff_prec(:,:)
    else
       plot_ps1 = ps_1(:)
       plot_ps2 = ps_2(:)
       plot_diff_ps = diff_ps(:)

       plot_prec1 = prec_1(:)
       plot_prec2 = prec_2(:)
       plot_diff_prec = diff_prec(:)
    end if


    ; Time Level
    psmin =  600.0d0;-psint
    psmax = 1040.0d0;+psint
    psint =   40.0d0
    diffpsmin = -2.0d0
    diffpsmax =  2.0d0
    diffpsint = 0.2d0

    precint =  0.02d0
    precmin =  0.0d0
    precmax =  0.2d0
    diffprecmin = -0.1d0
    diffprecmax =  0.1d0
    diffprecint =  0.02d0

    precint =  1.5d0
    precmin =  0.0d0
    precmax = 15.0d0
    diffprecmin = -15.0d0
    diffprecmax =  15.0d0
    diffprecint =   3.0d0


    xpos = new(4, double)
    ypos = new(3, double)

    x_margin = 0.03d0
    y_margin = -0.04d0

    width_g  = 0.22d0
    width_l  = 0.165d0 ;0.16 was small, 0.17 was large

    xpos(0)  = 0.05d0
    xpos(1)  = xpos(0) + width_g + x_margin
    xpos(2)  = xpos(1) + width_g + x_margin
    xpos(3)  = xpos(2) + width_g + x_margin

    ypos(0)  = 0.95d0
    ypos(1)  = ypos(0) - width_g - y_margin
    ypos(2)  = ypos(1) - width_g - y_margin
    adjust_l = 0.0275d0 ; 0.027 more down


    ;local_x1 = 105.0
    ;local_x2 = 150.0
    ;local_y1 =  20.0
    ;local_y2 =  50.0


    ; ######################## Box and Line ########################
    xndc_beg = new(1,float)
    yndc_beg = new(1,float)
    xndc_end = new(1,float)
    yndc_end = new(1,float)

    x_box = (/local_x1, local_x1, local_x2, local_x2, local_x1/)
    y_box = (/local_y1, local_y2, local_y2, local_y1, local_y1/)

    bres             = True
    ;bres@gsFillColor = "black"

    lnres                  = True
    lnres@gsLineThicknessF = 1.5d0  ; default is 1.5
    lnres@gsLineColor      = "black"


  
    plots = new (12, graphic)

    ; [[[[[[[[[[[[[[[[[[[[[[[[      Ps       ]]]]]]]]]]]]]]]]]]]]]]]]
    ; ######################## Global ########################
    ; draw contour line
    res@cnLinesOn      = False
    res@cnLineLabelsOn = False

    res@mpMinLonF =   0.0d0
    res@mpMaxLonF = 360.0d0
    res@mpMinLatF = -90.0d0
    res@mpMaxLatF =  90.0d0

    ; for Ps
    res@cnMaxLevelValF  = psmax
    res@cnMinLevelValF  = psmin
    res@cnLevelSpacingF = psint
;    res@cnLevelCount    = "5"

    res@vpXF = xpos(0)
    res@vpYF = ypos(0)
    res@vpHeightF = width_g
    res@vpWidthF  = width_g
    res@tiMainString    = "Ps_1"
    res@gsnCenterString = "min/max = "+min_ps_1+"/"+max_ps_1
    res@gsnLeftString   = plot_ps1@long_name
    plots(0) = gsn_csm_contour_map_ce(wks,plot_ps1,res)
    tmp = gsn_add_polyline(wks,plots(0),x_box,y_box,bres)
    draw(plots(0))

    res@vpXF = xpos(0)
    res@vpYF = ypos(1)
    res@vpHeightF = width_g
    res@vpWidthF  = width_g
    res@tiMainString    = "Ps_2"
    res@gsnCenterString = "min/max = "+min_ps_2+"/"+max_ps_2
    res@gsnLeftString   = plot_ps2@long_name
    plots(1) = gsn_csm_contour_map_ce(wks,plot_ps2,res)
    tmp = gsn_add_polyline(wks,plots(1),x_box,y_box,bres)
    draw(plots(1))

    res@cnMaxLevelValF  = diffpsmax
    res@cnMinLevelValF  = diffpsmin
    res@cnLevelSpacingF = diffpsint

    res@vpXF = xpos(0)
    res@vpYF = ypos(2)
    res@vpHeightF = width_g
    res@vpWidthF  = width_g
    res@tiMainString    = "Diff."
    ;res@gsnCenterString = "min/max = "+min_diff_ps+"/"+max_diff_ps
    ;res@gsnLeftString   = ""
    res@gsnCenterString = min_diff_ps+"/"+max_diff_ps
    res@gsnLeftString   = tx_pserr_l2
    plots(2) = gsn_csm_contour_map_ce(wks,plot_diff_ps,res)
    tmp = gsn_add_polyline(wks,plots(2),x_box,y_box,bres)

    ; error text
    pos_lon  = 30
    pos_lat1 = 80
    pos_lat2 = 70
    pos_lat3 = 60
    txres = True
    txres@txFont                = "helvetica-bold"
    txres@txFontQuality         = "High"
    txres@txJust                = "CenterLeft"
    txres@txFontHeightF         = 0.004d0
    txres@txBackgroundFillColor = 0 ;default -1
    txres@txFontColor           = "black"
    ;txres@txFontColor           = "red"
    txid = gsn_add_text(wks, plots(2), tx_pserr_l1,     pos_lon, pos_lat1, txres)
    txid = gsn_add_text(wks, plots(2), tx_pserr_l2,     pos_lon, pos_lat2, txres)
    txid = gsn_add_text(wks, plots(2), tx_pserr_linf,   pos_lon, pos_lat3, txres)

    draw(plots(2))


    ; ########################  Local ########################
    ; draw contour line
    res@cnLinesOn      = True
    res@cnLineLabelsOn = True

    ; Local Range
    res@mpMinLonF = local_x1
    res@mpMaxLonF = local_x2
    res@mpMinLatF = local_y1
    res@mpMaxLatF = local_y2

     ; Time Level
     ; Surface Pressure Range
     psint =  local_y1
     psmin =  700.0d0;-psint
     psmax = 1100.0d0;+psint
     ;diffpsmin = -1
     ;diffpsmax =  1
     ;diffpsint = 0.1d0

    ; for Ps
    res@cnMaxLevelValF  = psmax
    res@cnMinLevelValF  = psmin
    res@cnLevelSpacingF = psint
;    res@cnLevelCount    = "5"

    res@vpXF = xpos(1)
    res@vpYF = ypos(0)-adjust_l
    res@vpHeightF = width_l
    res@vpWidthF  = width_l
    res@tiMainString    = "Ps_1"
    res@gsnCenterString = "min/max = "+min_l_ps_1+"/"+max_l_ps_1
    ;res@gsnLeftString   = plot_ps1@long_name
    res@gsnLeftString   = ""
    plots(3) = gsn_csm_contour_map_ce(wks,plot_ps1,res)
    draw(plots(3))
    res@vpXF = xpos(1)
    res@vpYF = ypos(1)-adjust_l
    res@vpHeightF = width_l
    res@vpWidthF  = width_l
    res@tiMainString    = "Ps_2"
    res@gsnCenterString = "min/max = "+min_l_ps_2+"/"+max_l_ps_2
    ;res@gsnLeftString   = plot_ps2@long_name
    res@gsnLeftString   = ""
    plots(4) = gsn_csm_contour_map_ce(wks,plot_ps2,res)
    draw(plots(4))

    res@cnMaxLevelValF  = diffpsmax
    res@cnMinLevelValF  = diffpsmin
    res@cnLevelSpacingF = diffpsint

    ;res@cnInfoLabelOn = False
    res@cnLinesOn      = False
    res@cnLineLabelsOn = False
    res@vpXF = xpos(1)
    res@vpYF = ypos(2)-adjust_l
    res@vpHeightF = width_l
    res@vpWidthF  = width_l
    res@tiMainString    = "Diff."
    ;res@gsnCenterString = tx_l_pserr_l2
    ;res@gsnLeftString   = ""
    res@gsnCenterString = min_diff_l_ps+"/"+max_diff_l_ps
    res@gsnLeftString   = l_Pserr_l2
    plots(5) = gsn_csm_contour_map_ce(wks,plot_diff_ps,res)
    draw(plots(5))


    ; Draw line


    datatondc(plots(0), local_x2, local_y1, xndc_beg, yndc_beg)
    datatondc(plots(3), local_x1, local_y1, xndc_end, yndc_end)
    gsn_polyline_ndc(wks, (/xndc_beg,xndc_end/),(/yndc_beg,yndc_end/), lnres)
    datatondc(plots(0), local_x2, local_y2, xndc_beg, yndc_beg)
    datatondc(plots(3), local_x1, local_y2, xndc_end, yndc_end)
    gsn_polyline_ndc(wks, (/xndc_beg,xndc_end/),(/yndc_beg,yndc_end/), lnres)

    datatondc(plots(1), local_x2, local_y1, xndc_beg, yndc_beg)
    datatondc(plots(4), local_x1, local_y1, xndc_end, yndc_end)
    gsn_polyline_ndc(wks, (/xndc_beg,xndc_end/),(/yndc_beg,yndc_end/), lnres)
    datatondc(plots(1), local_x2, local_y2, xndc_beg, yndc_beg)
    datatondc(plots(4), local_x1, local_y2, xndc_end, yndc_end)
    gsn_polyline_ndc(wks, (/xndc_beg,xndc_end/),(/yndc_beg,yndc_end/), lnres)

    datatondc(plots(2), local_x2, local_y1, xndc_beg, yndc_beg)
    datatondc(plots(5), local_x1, local_y1, xndc_end, yndc_end)
    gsn_polyline_ndc(wks, (/xndc_beg,xndc_end/),(/yndc_beg,yndc_end/), lnres)
    datatondc(plots(2), local_x2, local_y2, xndc_beg, yndc_beg)
    datatondc(plots(5), local_x1, local_y2, xndc_end, yndc_end)
    gsn_polyline_ndc(wks, (/xndc_beg,xndc_end/),(/yndc_beg,yndc_end/), lnres)


    ; [[[[[[[[[[[[[[[[[[[[[[[[ precipitation ]]]]]]]]]]]]]]]]]]]]]]]]


    ; ######################## Global ########################
    gsn_define_colormap(wks,"MPL_Blues")
    ; for Precipitation
    res@cnMaxLevelValF  = precmax
    res@cnMinLevelValF  = precmin
    res@cnLevelSpacingF = precint
;    res@cnLevelCount    = "5"


    res@mpMinLonF =   0.0d0
    res@mpMaxLonF = 360.0d0
    res@mpMinLatF = -90.0d0
    res@mpMaxLatF =  90.0d0

    res@vpXF = xpos(2)
    res@vpYF = ypos(0)
    res@vpHeightF = width_g
    res@vpWidthF  = width_g
    res@tiMainString    = "Prec.-1"
    res@gsnCenterString = "min/max = "+min_prec_1+"/"+max_prec_1
    res@gsnLeftString   = plot_prec1@long_name
    plots(6) = gsn_csm_contour_map_ce(wks,plot_prec1,res)
    tmp = gsn_add_polyline(wks,plots(6),x_box,y_box,bres)
    draw(plots(6))

    res@vpXF = xpos(2)
    res@vpYF = ypos(1)
    res@vpHeightF = width_g
    res@vpWidthF  = width_g
    res@tiMainString    = "Prec.-2"
    res@gsnCenterString = "min/max = "+min_prec_2+"/"+max_prec_2
    res@gsnLeftString   = plot_prec2@long_name
    plots(7) = gsn_csm_contour_map_ce(wks,plot_prec2,res)
    tmp = gsn_add_polyline(wks,plots(7),x_box,y_box,bres)
    draw(plots(7))

    gsn_define_colormap(wks,"rainbow")

    res@cnMaxLevelValF  = diffprecmax
    res@cnMinLevelValF  = diffprecmin
    res@cnLevelSpacingF = diffprecint

    res@vpXF = xpos(2)
    res@vpYF = ypos(2)
    res@vpHeightF = width_g
    res@vpWidthF  = width_g
    res@tiMainString    = "Diff."
    ;res@gsnCenterString = "min/max = "+min_diff_prec+"/"+max_diff_prec
    ;res@gsnLeftString   = ""
    res@gsnCenterString = min_diff_prec+"/"+max_diff_prec
    res@gsnLeftString   = tx_precerr_l2
    plots(8) = gsn_csm_contour_map_ce(wks,plot_diff_prec,res)
    tmp = gsn_add_polyline(wks,plots(8),x_box,y_box,bres)


    ; error text
    pos_lon  = 30
    pos_lat1 = 80
    pos_lat2 = 70
    pos_lat3 = 60
    txres = True
    txres@txFont                = "helvetica-bold"
    txres@txFontQuality         = "High"
    txres@txJust                = "CenterLeft"
    txres@txFontHeightF         = 0.004d0
    txres@txBackgroundFillColor = 0 ;default -1
    txres@txFontColor           = "black"
    ;txres@txFontColor           = "red"
    txid = gsn_add_text(wks, plots(8), tx_precerr_l1,     pos_lon, pos_lat1, txres)
    txid = gsn_add_text(wks, plots(8), tx_precerr_l2,     pos_lon, pos_lat2, txres)
    txid = gsn_add_text(wks, plots(8), tx_precerr_linf,   pos_lon, pos_lat3, txres)
    draw(plots(8))



    ; ########################  Local ########################
;    gsn_define_colormap(wks,"MPL_Blues")
    ; for Precipitation
    res@cnMaxLevelValF  = 100.0d0;precmax
    res@cnMinLevelValF  = 0.0d0;precmin
    res@cnLevelSpacingF = 10;precint
    gsn_define_colormap(wks,"amwg")
    res@cnLevelSelectionMode = "ExplicitLevels" ; default "AutomaticLevels"
    res@cnLevels        = (/ 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0 /)
 
;    res@cnLevelCount    = "5"


    ; Local Range
    res@mpMinLonF = local_x1
    res@mpMaxLonF = local_x2
    res@mpMinLatF = local_y1
    res@mpMaxLatF = local_y2

    res@vpXF = xpos(3)
    res@vpYF = ypos(0)-adjust_l
    res@vpHeightF = width_l
    res@vpWidthF  = width_l
    res@tiMainString    = "Prec.-1"
    res@gsnCenterString = "min/max = "+min_l_prec_1+"/"+max_l_prec_1
    ;res@gsnLeftString   = plot_prec1@long_name
    res@gsnLeftString   = ""
    plots(9) = gsn_csm_contour_map_ce(wks,plot_prec1,res)
    draw(plots(9))

    res@vpXF = xpos(3)
    res@vpYF = ypos(1)-adjust_l
    res@vpHeightF = width_l
    res@vpWidthF  = width_l
    res@tiMainString    = "Prec.-2"
    res@gsnCenterString = "min/max = "+min_l_prec_2+"/"+max_l_prec_2
    ;res@gsnLeftString   = plot_prec2@long_name
    res@gsnLeftString   = ""
    plots(10) = gsn_csm_contour_map_ce(wks,plot_prec2,res)
    draw(plots(10))

    gsn_define_colormap(wks,"rainbow")
    res@cnLevelSelectionMode = "ManualLevels"

    res@cnMaxLevelValF  = diffprecmax
    res@cnMinLevelValF  = diffprecmin
    res@cnLevelSpacingF = diffprecint

    res@vpXF = xpos(3)
    res@vpYF = ypos(2)-adjust_l
    res@vpHeightF = width_l
    res@vpWidthF  = width_l
    res@tiMainString    = "Diff."
    res@gsnCenterString = min_diff_l_prec+"/"+max_diff_l_prec
    res@gsnLeftString   = l_precerr_l2
    plots(11) = gsn_csm_contour_map_ce(wks,plot_diff_prec,res)


    ; Draw line
    draw(plots(11))


    datatondc(plots(6), local_x2, local_y1, xndc_beg, yndc_beg)
    datatondc(plots(9), local_x1, local_y1, xndc_end, yndc_end)
    gsn_polyline_ndc(wks, (/xndc_beg,xndc_end/),(/yndc_beg,yndc_end/), lnres)
    datatondc(plots(6), local_x2, local_y2, xndc_beg, yndc_beg)
    datatondc(plots(9), local_x1, local_y2, xndc_end, yndc_end)
    gsn_polyline_ndc(wks, (/xndc_beg,xndc_end/),(/yndc_beg,yndc_end/), lnres)

    datatondc(plots(7), local_x2, local_y1, xndc_beg, yndc_beg)
    datatondc(plots(10), local_x1, local_y1, xndc_end, yndc_end)
    gsn_polyline_ndc(wks, (/xndc_beg,xndc_end/),(/yndc_beg,yndc_end/), lnres)
    datatondc(plots(7), local_x2, local_y2, xndc_beg, yndc_beg)
    datatondc(plots(10), local_x1, local_y2, xndc_end, yndc_end)
    gsn_polyline_ndc(wks, (/xndc_beg,xndc_end/),(/yndc_beg,yndc_end/), lnres)

    datatondc(plots(8), local_x2, local_y1, xndc_beg, yndc_beg)
    datatondc(plots(11), local_x1, local_y1, xndc_end, yndc_end)
    gsn_polyline_ndc(wks, (/xndc_beg,xndc_end/),(/yndc_beg,yndc_end/), lnres)
    datatondc(plots(8), local_x2, local_y2, xndc_beg, yndc_beg)
    datatondc(plots(11), local_x1, local_y2, xndc_end, yndc_end)
    gsn_polyline_ndc(wks, (/xndc_beg,xndc_end/),(/yndc_beg,yndc_end/), lnres)


    frame(wks)


    pres = True
    ;pres@gsnPanelRowSpec = True
;    gsn_panel(wks,plots,(/3,3/),pres)


  end if ; if (True)



  return(0)

end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; END FUNCTION PlotDiffRainStdOutKIM
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;








;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; FUNCTION PlotRainStdOutKIM
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function PlotRainStdOutKIM(names_in, filenames_in, filepaths_in)

begin

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Filenames and variables
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  names     = names_in
  filenames = filenames_in

  ncases    = dimsizes(names)
  nfiles    = dimsizes(filenames)


  itime = 0;itime_in
  ilev  = 0;ilev_in

  filepaths = filepaths_in
  nfp       = dimsizes(filepaths)
  if (ncases .ne. nfp) then
     print("check number of cases...")
     exit
  end if
  

  print("### Step0: Input arguments... : ")
  print(" ")

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Levels
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  isUP         = True


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Generate Positions
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print("### Step2 - START : Generate Position... : ")
  nplots  = ncases
  ; - outputs
  nx      = 0
  ny      = 0
  xpos    = new(nplots, double) 
  ypos    = new(nplots, double) 
  icanvas = new(nplots, integer) 
  xwidth  = 0.0d0
  ywidth  = 0.0d0
  ncanvas = 0
  nmod    = 0

;  err = GetCanvasPosition(nplots, nx, ny, xpos, ypos, xwidth, ywidth, ncanvas, nmod, icanvas)
;print(xpos)
;print(ypos)

  xpos(0) = 0.05d0
  xpos(1) = 0.43d0
  xpos(2) = 0.05d0
  xpos(3) = 0.43d0

  ypos(0) = 0.98d0
  ypos(1) = 0.98d0
  ypos(2) = 0.66d0
  ypos(3) = 0.66d0

  xwidth  = 0.30d0
  ywidth  = 0.30d0

  print("### Step2 -  END  : Generate Position... : ")
  print(" ")


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Make Canvas
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  wks_type = "ps"
  if(wks_type .eq."x11".or.wks_type.eq."png") then
     wks_type@wkWidth  = 1000
     wks_type@wkHeight = 1000
  end if
  if(wks_type.eq."ps".or.wks_type.eq."eps".or.wks_type.eq."epsi".or.wks_type.eq."pdf") then
     ;wks_type@wkOrientation = "landscape"
     wks_type@wkPaperSize    = "A4"
  end if

  outfilename = "output"
  wks        = gsn_open_wks(wks_type,outfilename)

  
  print("### Step3 - START : Variable Loop... : ")


  varname = (/"prec."/)

;  do ivar = 0, 0

    
  nvart    = 2
  varnamet = (/"crain","ncrain"/)
  
  print(" - Reading variables... : "+varnamet(0))

  infilenames = new(ncases, string)
  do i = 0, ncases-1
    infilenames(i) = filepaths(i)+"/"+filenames(0)
  end do


  infiles = addfiles(infilenames, "r")

  lats   = infiles[0]->lat
  lons   = infiles[0]->lon

  ndims   = new(ncases, integer)
  nhoriz  = new(ncases, integer)
  do i = 0, ncases-1
    vardims   = dimsizes(infiles[i]->$varnamet(0)$)
    ndims(i)  = dimsizes(vardims)
    nhoriz(i) = vardims(ndims(i)-1)
  end do


 
  var      = new((/ncases,nhoriz(0)/),double)
  var(:,:) = 0.0

  do icase = 0, ncases-1

    do ifile = 0, nfiles-1

      infile  = addfile(filepaths(icase)+"/"+filenames(ifile), "r")
  
      if (ndims(icase) .eq. 1) then
        do k = 0, nvart-1
          var(icase,:) = var(icase,:) + infile->$varnamet(k)$(:)
        end do
      else if (ndims(icase) .eq. 2) then
        do k = 0, nvart-1
          var(icase,:) = var(icase,:) + infile->$varnamet(k)$(itime,:)
        end do
      else if (ndims(icase) .eq. 3) then
        do k = 0, nvart-1
          var(icase,:) = var(icase,:) + infile->$varnamet(k)$(itime,ilev,:)
        end do
      end if
      end if
      end if

      delete(infile)

    end do

  end do
  

 

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; Plotting
  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  range = False

  res = PlotSingle(wks, xpos(0), ypos(0), xwidth, ywidth, names(0), var(0,:), range, lats, lons, isUP, True, False)
  res = PlotSingle(wks, xpos(1), ypos(1), xwidth, ywidth, names(1), var(1,:), range, lats, lons, isUP, True, False)
  res = PlotSingle(wks, xpos(2), ypos(2), xwidth, ywidth, names(2), var(2,:), range, lats, lons, isUP, True, False)
  res = PlotSingle(wks, xpos(3), ypos(3), xwidth, ywidth, names(3), var(3,:), range, lats, lons, isUP, True, False)
 

  print(" - Plotting... : ")
  print(" ")

  
  delete(varnamet)

  frame(wks)

  print("### Step3 -   END : Variable Loop... : ")
  print(" ")


  return(0)

end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; END FUNCTION PlotRainStdOutKIM
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;



































; backup






;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; FUNCTION PlotDiffRainStdOutKIM_old
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function PlotDiffRainStdOutKIM_old(filename1_in, filename2_in, filepath1_in, filepath2_in, reverse1_in, reverse2_in, itime_in, ilev_in, isPlot_in)

begin

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Filenames and variables
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  ;filename1 = new(nfiles1_d, string)
  ;filename2 = new(nfiles2_d, string)
  filename1 = filename1_in
  filename2 = filename2_in
  nfiles1 = dimsizes(filename1)
  nfiles2 = dimsizes(filename2)


  if (filepath1_in .eq. "") then
    filepath1  = "/scratch/jhkim/TestBed/KIM/Output/0.24.02/ne30/gnu/10hours/SH_KPPACK/403"
  else
    filepath1  = filepath1_in
  end if

  if (filepath2_in .eq. "") then
    filepath2  = "/scratch/jhkim/TestBed/Data/history"
  else
    filepath2  = filepath2_in
  end if

  reverse1     = reverse1_in
  reverse2     = reverse2_in
  isPlot       = isPlot_in

  print("### Step0: Input arguments... : ")
  print(" ")



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;; Loading input file and Coordinate
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;
;  infile1  = new(nfiles1, file)
;  infile2  = new(nfiles2, file)
;
;  do i = 0, nfiles1-1
;    infilename1   = filepath1+"/"+filename1(i)
;    infile1(i)      = addfile(infilename1,"r")
;    infilename2   = filepath2+"/"+filename2(i)
;    infile2(i)      = addfile(infilename2,"r")
;  end do
;
;  print("### Step1: Loading input file... : "+infilename1)
;  print("### Step1: Loading input file... : "+infilename2)
;  print(" ")

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Levels
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  nvars        = 2
  itime        = itime_in
  if (reverse1) then
    ilev1      = nlevs-1-ilev_in
  else
    ilev1      = ilev_in
  end if

  if (reverse2) then
    ilev2      = nlevs-1-ilev_in
  else
    ilev2      = ilev_in
  end if

  isUP         = True


   
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Generate Positions
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print("### Step2 - START : Generate Position... : ")
  ; Positions
  ; - input
  nplots  = 3*nvars
  ; - outputs
  nx      = 0
  ny      = 0
  xpos    = new(nplots, double) 
  ypos    = new(nplots, double) 
  icanvas = new(nplots, integer) 
  xwidth  = 0.0
  ywidth  = 0.0
  ncanvas = 0
  nmod    = 0
  err = GetCanvasPosition(nplots, nx, ny, xpos, ypos, xwidth, ywidth, ncanvas, nmod, icanvas)

  ypos(3) = ypos(3) - 0.1
  ypos(4) = ypos(4) - 0.1
  ypos(5) = ypos(5) - 0.1

  print("### Step2 -  END  : Generate Position... : ")
  print(" ")


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Make Canvas
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  wks_type = "ps"
  if(wks_type .eq."x11".or.wks_type.eq."png") then
     wks_type@wkWidth  = 1000
     wks_type@wkHeight = 1000
  end if
  if(wks_type.eq."ps".or.wks_type.eq."eps".or.wks_type.eq."epsi".or.wks_type.eq."pdf") then
     ;wks_type@wkOrientation = "landscape"
     wks_type@wkPaperSize    = "A4"
  end if

  outfilename = "output"
  wks        = gsn_open_wks(wks_type,outfilename)
  gsn_define_colormap(wks, "rainbow")

  
 

  print("### Step3 - START : Variable Loop... : ")


  varname = (/"ps", "prec."/)

  do ivar = 0, 1

    
     if (ivar .eq. 0) then
       nvart    = 1
       varnamet = (/"pint"/)
     else
       nvart    = 2
       varnamet = (/"crain","ncrain"/)
     end if
  
     print(" - Reading variables... : "+varnamet(0))

     infile1 = addfile(filepath1+"/"+filename1(0), "r")
     infile2 = addfile(filepath2+"/"+filename2(0), "r")
     lats   = infile1->lat
     lons   = infile1->lon

     vardim1 = dimsizes(infile1->$varnamet(0)$)
     ndims1  = dimsizes(vardim1)
     vardim2 = dimsizes(infile2->$varnamet(0)$)
     ndims2  = dimsizes(vardim2)
  
     ; ps
     if (ndims1 .eq. 1) then
       var1 = infile1->$varnamet(0)$(:)
       do k = 1, nvart-1
          var1(:) = var1(:) + infile1->$varnamet(k)$(:)
       end do
       do i = 1, nfiles1-1
         infile1 = addfile(filepath1+"/"+filename1(i), "r")
         var1(:) = var1(:) + infile1->$varnamet(0)$(:)
       do k = 1, nvart-1
          var1(:) = var1(:) + infile1->$varnamet(k)$(:)
       end do
       end do
     else if (ndims1 .eq. 2) then
       var1 = infile1->$varnamet(0)$(itime,:)
       do k = 1, nvart-1
          var1(:) = var1(:) + infile1->$varnamet(k)$(itime,:)
       end do
       do i = 1, nfiles1-1
         infile1 = addfile(filepath1+"/"+filename1(i), "r")
         var1(:) = var1(:) + infile1->$varnamet(0)$(itime,:)
       end do
       do k = 1, nvart-1
          var1(:) = var1(:) + infile1->$varnamet(k)$(itime,:)
       end do
     else if (ndims1 .eq. 3) then
       var1 = infile1->$varnamet(0)$(itime,ilev1,:)
       do k = 1, nvart-1
          var1(:) = var1(:) + infile1->$varnamet(k)$(itime,ilev1,:)
       end do
       do i = 1, nfiles1-1
         infile1 = addfile(filepath1+"/"+filename1(i), "r")
         var1(:) = var1(:) + infile1->$varnamet(0)$(itime,ilev1,:)
       end do
       do k = 1, nvart-1
          var1(:) = var1(:) + infile1->$varnamet(k)$(itime,ilev1,:)
       end do
     end if
     end if
     end if
     if (ivar .eq. 0) then
        var1 = var1/todouble(nfiles1)
;     else
;        var1 = var1*todouble(nfiles1)
     end if
   
     if (ndims2 .eq. 1) then
       var2 = infile2->$varnamet(0)$(:)
       do k = 1, nvart-1
          var2(:) = var2(:) + infile2->$varnamet(k)$(:)
       end do
       do i = 1, nfiles1-1
         infile2 = addfile(filepath2+"/"+filename2(i), "r")
         var2(:) = var2(:) + infile2->$varnamet(0)$(:)
       do k = 1, nvart-1
          var2(:) = var2(:) + infile2->$varnamet(k)$(:)
       end do
       end do
     else if (ndims2 .eq. 2) then
       var2 = infile2->$varnamet(0)$(itime, :)
       do k = 1, nvart-1
          var2(:) = var2(:) + infile2->$varnamet(k)$(itime,:)
       end do
       do i = 1, nfiles1-1
         infile2 = addfile(filepath2+"/"+filename2(i), "r")
         var2(:) = var2(:) + infile2->$varnamet(0)$(itime,:)
       do k = 1, nvart-1
          var2(:) = var2(:) + infile2->$varnamet(k)$(itime,:)
       end do
       end do
     else if (ndims2 .eq. 3) then
       var2 = infile2->$varnamet(0)$(itime, ilev2, :)
       do k = 1, nvart-1
          var2(:) = var2(:) + infile2->$varnamet(k)$(itime,ilev2,:)
       end do
       do i = 1, nfiles1-1
         infile2 = addfile(filepath2+"/"+filename2(i), "r")
         var2(:) = var2(:) + infile2->$varnamet(0)$(itime,ilev2,:)
       end do
       do k = 1, nvart-1
          var2(:) = var2(:) + infile2->$varnamet(k)$(itime,ilev2,:)
       end do
     end if
     end if
     end if
     if (ivar .eq. 0) then
        var2 = var2/todouble(nfiles2)
;     else
;        var2 = var2*todouble(nfiles1)
     end if
   
   
     ;-----------------------
     ; Statistic
     ;-----------------------
     Err_l1     = GetLError(var1(:),var2(:),"L1")
     Err_l2     = GetLError(var1(:),var2(:),"L2")
     Err_linf   = GetLError(var1(:),var2(:),"Linf")
     print(" - Rel. L1   = "+Err_l1)
     print(" - Rel. L2   = "+Err_l2)
     print(" - Rel. Linf = "+Err_linf)
   
     ;-----------------------
     ; Difference
     ;-----------------------
     diff_var = var2 - var1
     copy_VarCoords(var1, diff_var)
     
     diff_var@long_name = "L2 = "+tostring(Err_l2)
     if (isdefined("diff_var@units")) then
        diff_var@units     = var1@units
     else
        diff_var@units     = "unknown"
     end if
   
   
   ; plot range
     minvar = min(var1)
     maxvar = max(var1)
     minp = (toint(minvar/10.0))*10.0
     maxp = (toint(maxvar/10.0))*10.0
     if (minp .eq. 0.0) then
       minp = -1.0
     end if
     if (maxp .eq. 0.0) then
       maxp = 1.0
     end if
     intp = (maxp-minp)/10.0
    
   
     ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     ;; Plotting
     ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     fraction = 1.0/20.0
     ;res = PlotSingle(wks, xpos(iplot+0), ypos(iplot+0), xwidth, ywidth, varnames(ivar)+"-1", var1, minp, maxp, intp, lats, lons, isUP, isLocal_in)
     ;res = PlotSingle(wks, xpos(iplot+1), ypos(iplot+1), xwidth, ywidth, varnames(ivar)+"-2", var2, minp, maxp, intp, lats, lons, isUP, isLocal_in)
     ;res = PlotSingle(wks, xpos(iplot+2), ypos(iplot+2), xwidth , ywidth, "difference", diff_var, -(maxp-minp)*fraction, (maxp-minp)*fraction, (maxp-minp)*fraction/5.0, lats, lons, isUP, isLocal_in)
     print(" - Plotting... : ")
     print(" ")
   
     
     delete(vardim1) ; for reassignment
     delete(vardim2) ; for reassignment
     delete(varnamet)

  end do

  frame(wks)

  print("### Step3 -   END : Variable Loop... : ")
  print(" ")

  delete(infile1)
  delete(infile2)

  return(0)

end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; END FUNCTION PlotDiffRainStdOutKIM_old
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


































;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; FUNCTION PlotRainStdOutKIM_old
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function PlotRainStdOutKIM_old(names_in, filenames_in, filepaths_in)

begin

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Filenames and variables
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  names     = names_in
  filenames = filenames_in

  ncases    = dimsizes(names)
  nfiles    = dimsizes(filenames)


  itime = 0;itime_in
  ilev  = 0;ilev_in

  filepaths = filepaths_in
  nfp       = dimsizes(filepaths)
  if (ncases .ne. nfp) then
     print("check number of cases...")
     exit
  end if
  

  print("### Step0: Input arguments... : ")
  print(" ")

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Levels
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


  isUP         = True


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Generate Positions
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

  print("### Step2 - START : Generate Position... : ")
  nplots  = ncases
  ; - outputs
  nx      = 0
  ny      = 0
  xpos    = new(nplots, double) 
  ypos    = new(nplots, double) 
  icanvas = new(nplots, integer) 
  xwidth  = 0.0d0
  ywidth  = 0.0d0
  ncanvas = 0
  nmod    = 0

;  err = GetCanvasPosition(nplots, nx, ny, xpos, ypos, xwidth, ywidth, ncanvas, nmod, icanvas)
;print(xpos)
;print(ypos)

  xpos(0) = 0.05d0
  xpos(1) = 0.43d0
  xpos(2) = 0.05d0
  xpos(3) = 0.43d0

  ypos(0) = 0.98d0
  ypos(1) = 0.98d0
  ypos(2) = 0.66d0
  ypos(3) = 0.66d0

  xwidth  = 0.30d0
  ywidth  = 0.30d0

  print("### Step2 -  END  : Generate Position... : ")
  print(" ")


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; Make Canvas
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  wks_type = "ps"
  if(wks_type .eq."x11".or.wks_type.eq."png") then
     wks_type@wkWidth  = 1000
     wks_type@wkHeight = 1000
  end if
  if(wks_type.eq."ps".or.wks_type.eq."eps".or.wks_type.eq."epsi".or.wks_type.eq."pdf") then
     ;wks_type@wkOrientation = "landscape"
     wks_type@wkPaperSize    = "A4"
  end if

  outfilename = "output"
  wks        = gsn_open_wks(wks_type,outfilename)
;    gsn_draw_colormap(wks)
  res                        = True                   ; plot mods desired
  res@gsnDraw                = False
  res@gsnFrame               = False

  res@gsnSpreadColors = True;
  res@mpCenterLonF = 202.0d0;

;    res@gsnSpreadColorStart = 2;
;    res@gsnSpreadColorEnd = 201;

  ;res@tiYAxisString = "Latitude"
  ;res@tiXAxisString = "Longitude"

  res@cnFillOn = True;     
  res@cnLinesOn      = False
  res@cnLineLabelsOn = False
  res@cnInfoLabelOn = False
  res@lbLabelBarOn  = True
  res@lbLabelStride = 2

  res@cnMaxLevelValF  = 100.0d0;precmax
  res@cnMinLevelValF  = 0.0d0;precmin
  res@cnLevelSpacingF = 10;precint
  gsn_define_colormap(wks,"amwg")
  res@cnLevelSelectionMode = "ExplicitLevels" ; default "AutomaticLevels"
  res@cnLevels        = (/ 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0 /)

  
  print("### Step3 - START : Variable Loop... : ")


  varname = (/"prec."/)

;  do ivar = 0, 0

    
     nvart    = 2
     varnamet = (/"crain","ncrain"/)
  
     print(" - Reading variables... : "+varnamet(0))

     infilenames = new(ncases, string)
     do i = 0, ncases-1
       infilenames(i) = filepaths(i)+"/"+filenames(0)
     end do


     infiles = addfiles(infilenames, "r")

     lats   = infiles[0]->lat
     lons   = infiles[0]->lon

    ; for unique point grid
     if (isUP) then
       res@sfXArray              = lons;
       res@sfYArray              = lats;
     end if


     ndims   = new(ncases, integer)
     nhoriz  = new(ncases, integer)
     do i = 0, ncases-1
       vardims   = dimsizes(infiles[i]->$varnamet(0)$)
       ndims(i)  = dimsizes(vardims)
       nhoriz(i) = vardims(ndims(i)-1)
     end do


    
     var      = new((/ncases,nhoriz(0)/),double)
     var(:,:) = 0.0

     do icase = 0, ncases-1

       do ifile = 0, nfiles-1

         infile  = addfile(filepaths(icase)+"/"+filenames(ifile), "r")
  
         if (ndims(icase) .eq. 1) then
           do k = 0, nvart-1
             var(icase,:) = var(icase,:) + infile->$varnamet(k)$(:)
           end do
         else if (ndims(icase) .eq. 2) then
           do k = 0, nvart-1
             var(icase,:) = var(icase,:) + infile->$varnamet(k)$(itime,:)
           end do
         else if (ndims(icase) .eq. 3) then
           do k = 0, nvart-1
             var(icase,:) = var(icase,:) + infile->$varnamet(k)$(itime,ilev,:)
           end do
         end if
         end if
         end if

         delete(infile)

       end do
     end do
  


   
    
   
     ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
     ;; Plotting
     ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;





     plots = new (ncases, graphic)

     ; Local Range
     local_x1 = 105.0
     local_x2 = 150.0
     local_y1 =  20.0
     local_y2 =  50.0


     l_var1 = CutLocal(var(0,:), lons, lats, local_x1, local_x2, local_y1, local_y2)
     l_var2 = CutLocal(var(1,:), lons, lats, local_x1, local_x2, local_y1, local_y2)
     l_var3 = CutLocal(var(2,:), lons, lats, local_x1, local_x2, local_y1, local_y2)
     l_var4 = CutLocal(var(3,:), lons, lats, local_x1, local_x2, local_y1, local_y2)
     min_l_var1 = min(l_var1)
     max_l_var1 = max(l_var1)
     min_l_var2 = min(l_var2)
     max_l_var2 = max(l_var2)
     min_l_var3 = min(l_var3)
     max_l_var3 = max(l_var3)
     min_l_var4 = min(l_var4)
     max_l_var4 = max(l_var4)
  


     res@mpMinLonF = local_x1
     res@mpMaxLonF = local_x2
     res@mpMinLatF = local_y1
     res@mpMaxLatF = local_y2
 

     res@vpXF = xpos(0)
     res@vpYF = ypos(0)
     res@vpWidthF  = xwidth
     res@vpHeightF = ywidth
     res@tiMainString    = names(0)
     res@gsnCenterString = "min/max = "+min_l_var1+"/"+max_l_var1
     res@gsnLeftString   = ""
     plots(0) = gsn_csm_contour_map_ce(wks,var(0,:),res)
     draw(plots(0))
 

     res@vpXF = xpos(1)
     res@vpYF = ypos(1)
     res@vpWidthF  = xwidth
     res@vpHeightF = ywidth
     res@tiMainString    = names(1)
     res@gsnCenterString = "min/max = "+min_l_var2+"/"+max_l_var2
     res@gsnLeftString   = ""
     plots(1) = gsn_csm_contour_map_ce(wks,var(1,:),res)
     draw(plots(1))
 

     res@vpXF = xpos(2)
     res@vpYF = ypos(2)
     res@vpWidthF  = xwidth
     res@vpHeightF = ywidth
     res@tiMainString    = names(2)
     res@gsnCenterString = "min/max = "+min_l_var3+"/"+max_l_var3
     res@gsnLeftString   = ""
     plots(2) = gsn_csm_contour_map_ce(wks,var(2,:),res)
     draw(plots(2))
 

     res@vpXF = xpos(3)
     res@vpYF = ypos(3)
     res@vpWidthF  = xwidth
     res@vpHeightF = ywidth
     res@tiMainString    = names(3)
     res@gsnCenterString = "min/max = "+min_l_var4+"/"+max_l_var4
     res@gsnLeftString   = ""
     plots(3) = gsn_csm_contour_map_ce(wks,var(3,:),res)
     draw(plots(3))
 

     print(" - Plotting... : ")
     print(" ")
   
     
     delete(varnamet)

;  end do

  frame(wks)

  print("### Step3 -   END : Variable Loop... : ")
  print(" ")


  return(0)

end
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; END FUNCTION PlotRainStdOutKIM_old
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
