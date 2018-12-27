
function AveDev(var[*])
local dims,i,ave,dev
begin

  dims = dimsizes(var)
        
  output = new((/2/),double)


  ave = 0.0d0
  dev = 0.0d0
  do i = 0, dims-1
     ave = ave + var(i)
  end do
  ave = ave/dims


  do i = 0, dims-1
     dev = dev + (var(i)-ave)*(var(i)-ave)
  end do
  dev = sqrt(dev/(dims))
  ;dev = sqrt(dev/(dims-1))

  output(0) = ave
  output(1) = dev

  return(output)
end



function GetDev(var[*])
local dims,i,ave,dev
begin

  dims = dimsizes(var)
        
  output = -1.0


  ave = 0.0d0
  dev = 0.0d0
  do i = 0, dims-1
     ave = ave + var(i)
  end do
  ave = ave/dims


  do i = 0, dims-1
     dev = dev + (var(i)-ave)*(var(i)-ave)
  end do
  dev = sqrt(dev/(dims))

  output = dev

  return(output)
end



; 1 Dimension


function Get1LNorm(var1[*],pnorm:string)
local dims1,diff,i
begin

  if (pnorm.ne."L1").and.(pnorm.ne."L2").and.(pnorm.ne."Linf") then
     print("check pnorm...")
     return(-1)
  end if

  dims1  = dimsizes(var1)
        
  output = 0.0d0
  do i = 0, dims1-1
     if (pnorm .eq. "L1") then
        output = output + abs(var1(i))
     else if (pnorm .eq. "L2") then
        output = output + var1(i)*var1(i)
     end if
     end if
  end do

  if (pnorm .eq. "L2") then
     output = sqrt(output)
  else if (pnorm .eq. "Linf") then
     diff   = abs(var1)
     output = max(diff)
  end if
  end if

  return(output)
end


function GetLNorm(var1[*],var2[*],pnorm:string)
local dims1,dims2,diff,i
begin

  if (pnorm.ne."L1").and.(pnorm.ne."L2").and.(pnorm.ne."Linf") then
     print("check pnorm...")
     return(-1)
  end if

  dims1  = dimsizes(var1)
  dims2  = dimsizes(var2)
  if (dims1 .ne. dims2) then
     print("var1, var2 was different size of dimension.")
     return(-1)
  end if
        
  if (pnorm .eq. "Linf") then
     diff = var2-var1
  end if

  output = 0.0d0
  do i = 0, dims1-1
     if (pnorm .eq. "L1") then
        output = output + abs(var1(i)-var2(i))
     else if (pnorm .eq. "L2") then
        output = output + (var1(i)-var2(i))*(var1(i)-var2(i))
     else if (pnorm .eq. "Linf") then
        diff = abs(diff)
     end if
     end if
     end if
  end do

  if (pnorm .eq. "L2") then
     output = sqrt(output)
  else if (pnorm .eq. "Linf") then
     output = max(diff)
  end if
  end if

  return(output)
end



function GetLError(var1[*],var2[*],pnorm:string)
local dims1,dims2,tmp,i,j
begin
  if (pnorm.ne."L1").and.(pnorm.ne."L2").and.(pnorm.ne."Linf") then
     print("check pnorm...")
     return(-1)
  end if

  dims1  = dimsizes(var1)
  dims2  = dimsizes(var2)
  if (dims1 .ne. dims2) then
     print("var1, var2 was different size of dimension.")
     return(-1)
  end if
        
  output = 0.0d0
  tmp    = 0.0d0
  output = GetLNorm(var1, var2, pnorm)
  tmp    = Get1LNorm(var1, pnorm)
  if (tmp .eq. 0.0) then
     output = 0.0
  else
     output = output/tmp
  end if

  return(output)
end





function GetL1Norm(var1[*],var2[*])
local dims1,dims2,i
begin

  dims1  = dimsizes(var1)
  dims2  = dimsizes(var2)
  if (dims1 .ne. dims2) then
     print("var1, var2 was different size of dimension.")
     return(-1)
  end if
        
  output = 0.0d0
  do i = 0, dims1-1
     output = output + abs(var1(i)-var2(i))
  end do

  return(output)
end



function GetL2Norm(var1[*],var2[*])
local dims1,dims2,i
begin

  dims1  = dimsizes(var1)
  dims2  = dimsizes(var2)
  if (dims1 .ne. dims2) then
     print("var1, var2 was different size of dimension.")
     return(-1)
  end if
        
  output = 0.0d0
  do i = 0, dims1-1
     output = output + (var1(i)-var2(i))*(var1(i)-var2(i))
  end do
  output = sqrt(output)

  return(output)
end



function GetLinfNorm(var1[*],var2[*])
local dims1,dims2,diff,i
begin

  dims1  = dimsizes(var1)
  dims2  = dimsizes(var2)
  if (dims1 .ne. dims2) then
     print("var1, var2 was different size of dimension.")
     return(-1)
  end if
        
  output = 0.0d0
  diff = var2-var1
  do i = 0, dims1-1
     diff = abs(diff)
  end do
  output = max(diff)

  return(output)
end



function GetL1Error(var1[*],var2[*])
local dims1,dims2,tmp,i
begin

  dims1  = dimsizes(var1)
  dims2  = dimsizes(var2)
  if (dims1 .ne. dims2) then
     print("var1, var2 was different size of dimension.")
     return(-1)
  end if
        
  output = 0.0d0
  tmp    = 0.0d0
  output = GetL1Norm(var1, var2)
  do i = 0, dims1-1
     tmp = tmp + abs(var1(i))
  end do
  output = output/tmp

  return(output)
end



function GetL2Error(var1[*],var2[*])
local dims1,dims2,tmp,i
begin

  dims1  = dimsizes(var1)
  dims2  = dimsizes(var2)
  if (dims1 .ne. dims2) then
     print("var1, var2 was different size of dimension.")
     return(-1)
  end if
        
  output = 0.0d0
  tmp    = 0.0d0
  output = GetL2Norm(var1, var2)
  do i = 0, dims1-1
     tmp = tmp + var1(i)*var1(i)
  end do
  tmp = sqrt(tmp)
  output = output/tmp

  return(output)
end



function GetLinfError(var1[*],var2[*])
local dims1,dims2,tmp,diff,i
begin

  dims1  = dimsizes(var1)
  dims2  = dimsizes(var2)
  if (dims1 .ne. dims2) then
     print("var1, var2 was different size of dimension.")
     return(-1)
  end if
        
  output = 0.0d0
  tmp    = 0.0d0
  output = GetLinfNorm(var1, var2)
  diff   = abs(var1)
  tmp    = max(diff)
  output = output/tmp

  return(output)
end






; 2 Dimension


function Get1LNorm2D(var1[*][*],pnorm:string)
local dims1,diff,i
begin

  if (pnorm.ne."L1").and.(pnorm.ne."L2").and.(pnorm.ne."Linf") then
     print("check pnorm...")
     return(-1)
  end if

  dims1  = dimsizes(var1)
        
  output = 0.0d0
  do i = 0, dims1(0)-1
   do j = 0, dims1(1)-1
     if (pnorm .eq. "L1") then
        output = output + abs(var1(i,j))
     else if (pnorm .eq. "L2") then
        output = output + var1(i,j)*var1(i,j)
     end if
     end if
   end do
  end do

  if (pnorm .eq. "L2") then
     output = sqrt(output)
  else if (pnorm .eq. "Linf") then
     diff   = abs(var1)
     output = max(diff)
  end if
  end if

  return(output)
end


function GetLNorm2D(var1[*][*],var2[*][*],pnorm:string)
local dims1,dims2,diff,i,j
begin

  if (pnorm.ne."L1").and.(pnorm.ne."L2").and.(pnorm.ne."Linf") then
     print("check pnorm...")
     return(-1)
  end if

  dims1  = dimsizes(var1)
  dims2  = dimsizes(var2)
  if ((dims1(0).ne.dims2(0)) .or. (dims1(1).ne.dims2(1))) then
     print("var1, var2 was different size of dimension.")
     return(-1)
  end if
        
  if (pnorm .eq. "Linf") then
     diff = var2-var1
  end if
        
  output = 0.0d0
  do i = 0, dims1(0)-1
   do j = 0, dims1(1)-1
     if (pnorm .eq. "L1") then
        output = output + abs(var1(i,j)-var2(i,j))
     else if (pnorm .eq. "L2") then
        output = output + (var1(i,j)-var2(i,j))*(var1(i,j)-var2(i,j))
     else if (pnorm .eq. "Linf") then
        diff = abs(diff)
     end if
     end if
     end if
   end do
  end do

  if (pnorm .eq. "L2") then
     output = sqrt(output)
  else if (pnorm .eq. "Linf") then
     output = max(diff)
  end if
  end if


  return(output)
end



function GetLError2D(var1[*][*],var2[*][*],pnorm:string)
local dims1,dims2,tmp,i,j
begin
  if (pnorm.ne."L1").and.(pnorm.ne."L2").and.(pnorm.ne."Linf") then
     print("check pnorm...")
     return(-1)
  end if

  dims1  = dimsizes(var1)
  dims2  = dimsizes(var2)
  if ((dims1(0).ne.dims2(0)) .or. (dims1(1).ne.dims2(1))) then
     print("var1, var2 was different size of dimension.")
     return(-1)
  end if
        
  output = 0.0d0
  tmp    = 0.0d0
  output = GetLNorm2D(var1, var2, pnorm)
  tmp    = Get1LNorm2D(var1, pnorm)
  output = output/tmp

  return(output)
end



function CutLocal(var_in, lons_in, lats_in, min_lon, max_lon, min_lat, max_lat)
begin

  nvars     = dimsizes(var_in)
  nvars_out = 0
  
  do i = 0, nvars-1
    if (lons_in(i) .ge. min_lon .and. lons_in(i) .le. max_lon) then
    if (lats_in(i) .ge. min_lat .and. lats_in(i) .le. max_lat) then
;       var_out = array_append_record(var_out, var_in, 0)
       nvars_out = nvars_out + 1
    end if
    end if
  end do

  ivar    = 0
  var_out = new(nvars_out, typeof(var_in))
  do i = 0, nvars-1
    if (lons_in(i) .ge. min_lon .and. lons_in(i) .le. max_lon) then
    if (lats_in(i) .ge. min_lat .and. lats_in(i) .le. max_lat) then
       var_out(ivar) = var_in(i)
       ivar = ivar + 1
    end if
    end if
  end do

  return(var_out)
end
