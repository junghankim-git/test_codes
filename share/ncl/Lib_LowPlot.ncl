;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; FUNCTION GetCanvasPosition
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function GetCanvasPosition(nplots_in, nx, ny, xpos, ypos, xwidth, ywidth, ncanvas, nmod, icanvas)

begin

  ; canvas size
  ; org: 3 x 5
  nx          = 3
  ny          = 5
  xpos_c      = new(nx, double)  ; position x in canvas
  ypos_c      = new(ny, double)  ; position y in canvas

  margin_x    = 0.02d0
  margin_y    = 0.02d0
  width_tx    = 1.0d0-(2.0d0*margin_x)
  width_ty    = 1.0d0-(2.0d0*margin_y)

  width_x     = (width_tx-margin_x)/todouble(nx)
  width_y     = width_ty/todouble(ny-1.0d0)           ; hand-control (old: 0.9)

  xpos_c(0)   = margin_x
  ypos_c(0)   = 1.0d0 - margin_y

  back_margin_x = width_tx*0.000d0                 ; hand-control
  back_margin_y = width_ty*0.070d0                 ; hand-control (old: 0.054)
  do i = 1, nx-1
    xpos_c(i) = xpos_c(i-1)+width_x-back_margin_x
  end do
  do i = 1, ny-1
    ypos_c(i) = ypos_c(i-1)-width_y+back_margin_y
  end do


  do iplot = 0, nplots_in-1
    nn = mod(iplot, nx*ny)
    ii = mod(nn, nx)
    jj = nn/nx
    xpos(iplot)    = xpos_c(ii)
    ypos(iplot)    = ypos_c(jj)
    icanvas(iplot) = iplot/(nx*ny)
  end do
  xwidth = width_x
  ywidth = width_y


  ; Canvas
  ncanvas = nplots_in/(nx*ny)
  nmod    = mod(nplots_in, nx*ny)
  if (nmod .gt. 0) then
    ncanvas = ncanvas + 1
  end if

  return(0)

end








;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; FUNCTION PlotSingle
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function PlotSingle(wks_in:graphic, xpos_in, ypos_in, xwidth_in, ywidth_in, varname_in, var_in, range_in, lats_in, lons_in, isUP, isLocal, isPrec)

begin

  ;-----------------------
  ; Set Attributes
  ;-----------------------
  gsn_define_colormap(wks_in, "rainbow")

  minvar = min(var_in)
  maxvar = max(var_in)

  att_in = True
; UP or LL
  if (isUP) then
    att_in@sfXArray     = lons_in
    att_in@sfYArray     = lats_in
  end if

; draw contour line
  att_in@cnLinesOn      = False
  att_in@cnLineLabelsOn = False
  ;att_in@cnLevelSelectionMode = ""
  att_in@cnFillOn        = True
  att_in@cnInfoLabelOn   = True

; plot range
  ncolors  = 9.0


  if (range_in) then

     l_minp = range_in@min
     l_maxp = range_in@max
     l_intp = range_in@int

  else

     if (abs(minvar) .gt. 10.0d) then
        l_minp = (toint(minvar/10.0d0))*10.0d0
     else if (abs(minvar) .gt. 1.0d) then
        l_minp = (toint(minvar))*1.0d0
     else if (abs(minvar) .gt. 0.1d) then
        l_minp = (toint(minvar*10.0))/10.0d0
     else
        l_minp = -0.001
     end if
     end if
     end if

     if (abs(maxvar) .gt. 10.0d) then
        l_maxp = (toint(maxvar/10.0d0))*10.0d0
     else if (abs(maxvar) .gt. 1.0d) then
        l_maxp = (toint(maxvar))*1.0d0
     else if (abs(maxvar) .gt. 0.1d) then
        l_maxp = (toint(maxvar*10.0))/10.0d0
     else
        l_maxp = 0.001
     end if
     end if
     end if

     l_intp = (l_maxp-l_minp)/10.0d0

  end if

  att_in@cnLevelSelectionMode = "ManualLevels" ; "AutomaticLevels", "ManualLevels", "ExplicitLevels"
  att_in@cnMinLevelValF  = l_minp
  att_in@cnMaxLevelValF  = l_maxp
  att_in@cnLevelSpacingF = l_intp
  ;att_in@cnLevelCount    = "6"
  att_in@cnInfoLabelOn   = False
  att_in@lbOrientation   = "Vertical"          ; "Horizontal", "Vertical"

  range_x1 =   0.0
  range_x2 = 360.0
  range_y1 = -90.0
  range_y2 =  90.0

  if (isLocal) then
     range_x1 = 105.0
     range_x2 = 150.0
     range_y1 =  20.0
     range_y2 =  50.0
     l_var = CutLocal(var_in(:), lons_in, lats_in, range_x1, range_x2, range_y1, range_y2)
     minvar = min(l_var)
     maxvar = max(l_var)
     
     if (isPrec) then
        att_in@cnLevelSelectionMode = "ExplicitLevels"
        att_in@cnMinLevelValF  = 0.0
        att_in@cnMaxLevelValF  = 100.0
        att_in@cnLevelSpacingF = 10
        att_in@cnLevels        = (/ 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0 /)
        gsn_define_colormap(wks_in,"amwg")
        att_in@lbOrientation   = "Horizontal"          ; "Horizontal", "Vertical"
     end if

  end if



; canvas range & position
  ;att_in@mpCenterLonF    = 202.0        ; default "0"
  att_in@mpMinLonF       = range_x1
  att_in@mpMaxLonF       = range_x2
  att_in@mpMinLatF       = range_y1
  att_in@mpMaxLatF       = range_y2

; Label & Draw options
  ; gsnDraw, gsnFrame = False => need frame(wks) for multi plot
  att_in@gsnDraw         = False       ; False: no plot
  att_in@gsnFrame        = False       ; same as gsnDraw?
  att_in@gsnSpreadColors = True        ; maybe fixed as True
  att_in@gsnCenterString = "min/max = "+minvar+"/"+maxvar          ; default: none
  ;att_in@gsnLeftString   = ""          ; default: "long_name"
  ;att_in@gsnRightString  = ""          ; default: "units"
  att_in@tiMainString    = varname_in
  att_in@tiYAxisString   = ""          ; default: X's "long_name"
  att_in@tiXAxisString   = ""          ; default: Y's "long_name"


; plot position
  att_in@vpXF            = xpos_in
  att_in@vpYF            = ypos_in
  att_in@vpHeightF       = xwidth_in
  att_in@vpWidthF        = ywidth_in


  ;-----------------------
  ; Plot
  ;-----------------------
  plot = gsn_csm_contour_map_ce(wks_in, var_in, att_in)
  draw(plot)

  return(0)

end






;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; FUNCTION PlotDifference
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function PlotDifference(wks_in:graphic, xpos_in, ypos_in, xwidth_in, ywidth_in, varname_in, var1_in, var2_in, diff_in, min_plot, max_plot, int_plot, lats_in, lons_in, isUP, isLocal)

begin

  range = True
  range@min = min_plot
  range@max = max_plot
  range@int = int_plot

  fraction = 1.0d0/20.0d0
  diff_range = True
  diff_range@min = -(max_plot-min_plot)*fraction
  diff_range@max =  (max_plot-min_plot)*fraction
  diff_range@int =  (max_plot-min_plot)*fraction/5.0d0
  fraction = 1.0/20.0

  ;-----------------------
  ; Plotting
  ;-----------------------
  res = PlotSingle(wks_in, xpos_in(0), ypos_in(0), xwidth_in, ywidth_in, varname_in+"-1", var1_in, range, lats_in, lons_in, isUP, isLocal, False)
  res = PlotSingle(wks_in, xpos_in(1), ypos_in(1), xwidth_in, ywidth_in, varname_in+"-2", var2_in, range, lats_in, lons_in, isUP, isLocal, False)
  res = PlotSingle(wks_in, xpos_in(2), ypos_in(2), xwidth_in , ywidth_in, "difference", diff_in, diff_range, lats_in, lons_in, isUP, isLocal, False)

  return(0)

end



