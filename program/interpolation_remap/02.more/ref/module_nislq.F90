#include <define.h> 
#ifndef KIM
!-------------------------------------------------------------------------------
   module nislq
!-------------------------------------------------------------------------------
!
!  abstract :
!
!  history log :
!    2016-09-25   myung-seo koo      clean-up
!    2017-12-05   myung-seo koo      optional remap (iopt)
!
!  variable :
!    iopt: remap option
!
!      [ ppm(-1) in the horizontal ]
!         = -1 for piecewise parabolic method (PPM)    (henry juang)
!            0 for parabolic spline    method (PSM)    (CAM)
!            1 for piecewise parabolic method (PPM) v1 (CAM)
!            2 for piecewise parabolic method (PPM) v2 (CAM)
!            3 for piecewise parabolic method (PPM) v3 (jh.kim)
!            4 for piecewise linear    method (PLM)    (jh.kim)
!            5 for piecewise constant  method (PCM)    (jh.kim)
!            6 for parabolic spline    method (PSM)    (jh.kim)
!            7 for quartic   spline    method (QSM)    (jh.kim)
!
!      [ in the both direction]
!         = 16 for parabolic spline    method (PSM)    (jh.kim)
!         = 17 for quartic   spline    method (QSM)    (jh.kim)... not yet
!
!-------------------------------------------------------------------------------
   use comio          , only : iope
#ifdef NISLQ_NEW
   use piecewise_remap, only : piecewise_interface
   use spline_remap   , only : spline_interface,spline_cyclic_interface
#endif
!
   implicit none
!
! Common variables
!
#ifdef NISLQ_NEW
   integer, parameter                    ::  iopt = 17
#else
   integer, parameter                    ::  iopt = -1
#endif
   integer                               ::  nx         ,& ! global longitude
                                             my         ,& ! global latitude
                                             my_max     ,& ! local  latitude
                                             lev        ,& ! global level
                                             ncld       ,& ! # of hydrometeor.
                                             nlevs      ,&
                                             nlevsp     ,&
!                                            jlistnum   ,& ! myrank lat. length
!                                            lat_s      ,& ! myrank lat. start
!                                            lat_e      ,& ! myrank lat. end
                                             lonfull    ,& ! nx
                                             latfull    ,& ! my*2
                                             lonhalf    ,& ! nx/2
                                             lathalf    ,& ! my
                                             lonpart    ,& ! lonhalf/nsize+1
                                             latpart    ,& ! lathalf/nsize+1
                                             mylonlen   ,& ! myrank lon. length
                                             mylatlen      ! myrank lat. length
   integer, allocatable, dimension(:)    ::  lonstr     ,& ! allpe lon. start
                                             lonlen     ,& ! allpe lon. length
                                             latstr     ,& ! allpe lat. start
                                             latlen     ,& ! allpe lat. length
                                             truej      ,&
                                             shflj      ,&
                                             jlist1
   real   , allocatable, dimension(:)    ::  glat       ,& ! gaussian lat.(my)
                                             gglat      ,& ! gaussian lat.(my*2)
                                             gglati     ,&
                                             gglon      ,& ! gaussian lon. (nx)
                                             ggloni     ,&
                                             racos      ,& ! 1./(Re*cos)
                                             racos2        ! 1./(Re.cos^2)
   real   , allocatable, dimension(:,:,:) :: slq_q1     ,& ! moisture at n-1
                                             slq_q2     ,& ! moisture at n
                                             slq_q3        ! moisture at n+1
   real   , allocatable, dimension(:,:)   :: slq_p2
   real   , allocatable, dimension(:,:,:) :: slq_u2,slq_v2,slq_w2
                                             
   contains
!-------------------------------------------------------------------------------
   subroutine nislq_init(nsize,myrank,colrad,rbs2)
!-------------------------------------------------------------------------------
   use paramodel, only : lonf_,latg_,latg2_,levs_,ntotal_
   use paramodel, only : LEVSS
   use comio    , only : iope
#ifdef DFS
   use dfsvar   , only : ib,jbw,levh,latdef
#else
   use paramodel, only : LONF2S,LATG2S,levh_
#ifdef MP
   use commpi   , only : latdef
#else
   use comfgrid , only : latdef
#endif
#endif /* DFS end */
   use constant , only : pi_,rrerth_
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!
! passing variables
!
   integer,intent(in)                    :: nsize, myrank
   real   ,intent(in), dimension(latg2_) :: colrad,rbs2
!
! local variables
!
   integer                               :: myh,my2,i,j,ii,jj,js,k
   real                                  :: hfpi,twopi,dlat,dlon
   integer                               :: n,nm,nr,lonp
   real, dimension(latg2_)               :: colatrad
   integer                               :: j1,j2
!
! assign parameter
!
   nx      = lonf_
   my      = latg_
   lev     = levs_
!sldp   ncld    = 1+ntotal_
   ncld    = ntotal_
   nlevs   = ncld*lev
   nlevsp  = ncld*LEVSS
!
   lonfull = nx
   lonhalf = nx / 2             ! nx has to be even
   lonpart = lonhalf/nsize+1
   latfull = my * 2
   lathalf = my 
!  latpart = lathalf/nsize+1
   latpart = 2*((lathalf/2)/nsize+1)
#ifdef MP
   my_max  = latpart
#else
   my_max  = latg_
#endif
!
! constant
!
   hfpi  = pi_ * 0.5
   twopi = pi_ * 2.0
   myh   = my/2
   my2   = my*2
!
! humidity mixing ratio and other for nislq
!
#ifdef DFS
   allocate(  slq_q1 (ib,jbw,levh)                                            ,&
              slq_q2 (ib,jbw,levh)                                            ,&
              slq_q3 (ib,jbw,levh)                                            ,&
              slq_p2 (ib,jbw)                                                 ,&
              slq_u2 (ib,jbw,levs_)                                           ,&
              slq_v2 (ib,jbw,levs_)                                           ,&
              slq_w2 (ib,jbw,levs_+1)                                          )
#else /* SPH */
   allocate(  slq_q1 (LONF2S,levh_  ,LATG2S)                                  ,&
              slq_q2 (LONF2S,levh_  ,LATG2S)                                  ,&
              slq_q3 (LONF2S,levh_  ,LATG2S)                                  ,&
              slq_p2 (LONF2S        ,LATG2S)                                  ,&
              slq_u2 (LONF2S,levs_  ,LATG2S)                                  ,&
              slq_v2 (LONF2S,levs_  ,LATG2S)                                  ,&
              slq_w2 (LONF2S,levs_+1,LATG2S)                                   )
#endif /* DFS end */
!
! initialize slq_q
!
   slq_q1=0.
   slq_q2=0.
   slq_q3=0.
   slq_p2=0.
   slq_u2=0.
   slq_v2=0.
   slq_w2=0.
!
! informations for latitude
!
   allocate(jlist1(latg2_))
   allocate(racos(latg2_))
   allocate(racos2(latg2_))
!
! define latitude
!
   do j = 1,latg2_
#ifdef MP
     jlist1(j)=j
#else
     jlist1(j)=j
#endif
   enddo
!
! define colatitude in radian
!
   do j = 1,latg2_
     colatrad(j)=colrad(j)
   enddo
!
! ------------------- gaussian latitude 
!
   allocate ( glat(my) )
   allocate ( gglat(my2), gglati(my2+1) )
!
   do j = 1,myh
     glat(j)=hfpi-colatrad(j)
     glat(my+1-j)=-glat(j)
   enddo
!
   do j = 1,myh       ! co-latitude
     gglat(j) = hfpi - glat(j)
     gglat(my+1-j) = hfpi + glat(j)
   enddo
   do j = my+1,my2
     gglat(j) = twopi - gglat(my2+1-j)
   enddo
!
   gglati(myh+1) = hfpi
   dlat = (gglati(myh+1)-gglat(myh))*2.0
   do j = myh,2,-1
     gglati(j) = gglati(j+1) - dlat
     dlat      = (gglati(j) - gglat(j-1))*2.0
   enddo
   gglati(1) = 0.0
   do j = myh+2,my
     gglati(j) = pi_ - gglati(my+2-j)
   enddo
   gglati(my+1) = pi_
   do j = my+2,my2
     gglati(j) = twopi - gglati(my2+2-j)
   enddo
   gglati(my2+1) = twopi
#ifdef SLDBG
!
   if( iope ) then
     print *,' j interface =1    gglati=',gglati(1)
     do j = 1,my2
       print *,'               j mean =',j,'  gglat= ',gglat(j)
       print *,' j interface =',j+1,' gglati= ',gglati(j+1)
     enddo
   endif
#endif
!
!  determin the longitude with full grid
!
   allocate( gglon(nx), ggloni(nx+1) )
   dlon = twopi / nx
   do i = 1,nx
     gglon(i)=(i-1)*dlon
   enddo
   do i = 2,nx
     ggloni(i)=0.5*(gglon(i-1)+gglon(i))
   enddo
   ggloni(     1)=ggloni(   2)-dlon
   ggloni(nx+1)=ggloni(nx)+dlon
!
#ifdef SLDBG
   if( iope ) then
     print *,' ------ total edge nx=',nx,' -------'
     print *,' i edge number=1    ggloni=',ggloni(1)
     do i = 1,nx
       print *,'               i cell number=',i,'  gglon= ',gglon(i)
       print *,' i edge number=',i+1,' ggloni= ',ggloni(i+1)
     enddo
   endif
#endif
!
! --------------------- for parallel --------------------
!
! in nisl, we do great circle, so nx and my2 are full
! transpose will between (lonfull,lev,latpart)  by (nx  ,lev,my/nsize+1)
!                    and (latfull,lev,lonpart)  by (my*2,lev,nx/2/nsize+1)
!
   allocate( lonstr(nsize), lonlen(nsize) )
   allocate( latstr(nsize), latlen(nsize) )
!
! equally distribute len, no location is considered
!
   lonlen(1:nsize) = 0
   i=1
   do ii = 1,lonpart
     do n = 1,nsize
       if(i.le.lonhalf) then
         lonlen(n) = lonlen(n)+1
         i=i+1
       endif
     enddo
   enddo
!
! sequential location for longitude
!
   nm=1
   do n = 1,nsize
     lonstr(n) = nm
     nm = nm + lonlen(n)
   enddo
!
! check make_list  to have consistent latitude number for each pe here
!
! refer equdis
   latlen(1:nsize) = 0
   i=1
   n=1
   do jj = 1,myh
     latlen(n)=latlen(n)+1
     n=n+i
     if(n.eq.nsize+1) then
       i=-1
       n=n+i
     endif
     if(n.eq.0) then
       i=1
       n=n+i
     endif
   enddo
   latlen=latlen*2
#ifdef SLDBG
   if(iope) then
     do n = 1,nsize
       print *,'n,latlen=',n,latlen(n)
     enddo
   endif
#endif
!
! sequential location for latitude
!
   nm=1
   do n = 1,nsize
     latstr(n) = nm
     nm = nm + latlen(n)
   enddo
!
!   jlistnum=latlen(myrank+1)
!   lat_s=latstr(myrank+1)
!   lat_e=latstr(myrank+1)+latlen(myrank+1)-1
#ifdef SLDBG
!
! check by print
!
   if( iope ) then
     do n = 1,nsize
       print *,' pe lonstr lonlen ',n,lonstr(n),lonlen(n)
       print *,' pe latstr latlen ',n,latstr(n),latlen(n)
     enddo
   endif
#endif
!
! make true latitude index follow make_list
!
   allocate( truej(my), shflj(my) )
!!   j=1
!!   do jj=1,my/nsize+1
!!     js = jj
!!     do n=1,nsize
!!       if( j.le.my ) then
!!         truej (js) = j
!!         shflj (j ) = js
!!         j = j + 1
!!       endif
!!       js = js + latlen(n)
!!     enddo
!!   enddo
   do jj = 1,latg2_
     j1=jj*2-1
     j2=jj*2
     truej(j1)=latdef(jj)
     truej(j2)=latg_-latdef(jj)+1
   enddo
#ifdef SLDBG
!
! check by print
!
   if(iope) then
!     do j=1,my
!       print *,' true j=',j,' to shaffled j=',shflj(j),   &
!                 ' back to true j=',truej(shflj(j))
!     enddo
     do j = 1,my
       print *,' shfl j=',j,' to true j=',truej(j)                            ,&
                 ' back to shfl j=',shflj(truej(j))
     enddo
     j=0
     do n = 1,nsize
       print *,' --- start pe =',n-1
       do jj = 1,latlen(n)
         j=j+1
         print *,' shaffled j=',j,' to true j=',truej(j)
       enddo
     enddo
   endif
#endif
!
   mylonlen = lonlen(myrank+1)
   mylatlen = latlen(myrank+1)
! debug
!  if( n.eq.n ) then
!    call mpe_finalize
!    stop
!  endif
!
   return
   end subroutine nislq_init
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine cyclic_cell_massadvx(levs,ncld,delt,uc,qq,mass)
!-------------------------------------------------------------------------------
!
! compute local positive advection with mass conservation
! qq is advected by uc which is in radiance/sec from past to next position
!
! author: hann-ming henry juang 2008
!
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   real   , parameter                    ::  fa1 = 9./16.                     ,&
                                             fa2 = 1./16.
   integer                               ::  levs,ncld,mass
   integer                               ::  i,k,im,n
   real                                  ::  delt, rm2, dist, sc
   real   , dimension(lonfull,levs)      ::  uc
   real   , dimension(lonfull,levs,ncld) ::  qq
   real   , dimension(lonfull)           ::  past,next,da,dxfact
   real   , dimension(lonfull+1)         ::  xpast,xnext,uint
!
! preparations ---------------------------
!
! x is equal grid spacing, so location can be specified by grid point number
!
   im = lonfull
!
   do k = 1,levs
!
! 4th order interpolation from mid point to cell interfaces
!
     do i = 3,im-1
       uint(i)=fa1*(uc(i,k)+uc(i-1,k))-fa2*(uc(i+1,k)+uc(i-2,k))
     enddo
     uint(2)=fa1*(uc(2,k)+uc(1 ,k))-fa2*(uc(3,k)+uc(im  ,k))
     uint(1)=fa1*(uc(1,k)+uc(im,k))-fa2*(uc(2,k)+uc(im-1,k))
     uint(im+1)=uint(1)
     uint(im  )=fa1*(uc(im,k)+uc(im-1,k)) -fa2*(uc(1,k)+uc(im-2,k))
!
! compute past and next positions of cell interfaces
!
     do i = 1,im+1
       dist     = uint(i) * delt
       xpast(i) = ggloni(i) - dist
       xnext(i) = ggloni(i) + dist
     enddo
!      
     if( mass.eq.1 ) then
       do i = 1,im
         dxfact(i) = (xpast(i+1)-xpast(i)) / (xnext(i+1)-xnext(i))
       enddo
     endif
!
!  mass positive advection
!
     sc=ggloni(im+1)-ggloni(1)
     do n = 1,ncld
       past(1:im) = qq(1:im,k,n)
#ifdef NISLQ_NEW
       if ( iopt.eq.16 .or. iopt.eq.17 ) then
         call spline_cyclic_interface(im,ggloni,past,im,xpast ,da  ,2,.true.,2)
         call spline_cyclic_interface(im,xnext ,da  ,im,ggloni,next,2,.true.,2)
       else
#endif
         call cyclic_cell_ppm_intp(ggloni,past,xpast,da,im,im,im,sc)
         if(mass.eq.1) da(1:im) = da(1:im) * dxfact(1:im)
         call cyclic_cell_ppm_intp(xnext,da,ggloni,next,im,im,im,sc)
#ifdef NISLQ_NEW
       endif
#endif
       qq(1:im,k,n) = next(1:im)
     enddo
!
   enddo
!
   return
   end subroutine cyclic_cell_massadvx
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine cyclic_cell_massadvy(levs,ncld,delt,vc,qq,mass)
!-------------------------------------------------------------------------------
!
! compute local positive advection with mass conserving
! qq will be advect by vc from past to next location with 2*delt
!
! author: hann-ming henry juang 2007
!
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   real   , parameter                    ::  fa1 = 9./16.                     ,&
                                             fa2 = 1./16.
   integer                               ::  levs,ncld,mass                   ,&
                                             n,k,j,jm,jm2
   real                                  ::  delt,sc
   real   , dimension(latfull,levs)      ::  vc
   real   , dimension(latfull,levs,ncld) ::  qq
   real   , dimension(latfull)           ::  var,past,da,next,dyfact
   real   , dimension(latfull+1)         ::  ypast,ynext,dist
!
! preparations ---------------------------
!
   jm   = lathalf
   jm2  = latfull
!
   do k = 1,levs
!
     do j = 1,jm
       var(j)      =-vc(j   ,k) * delt
       var(j+jm)   = vc(j+jm,k) * delt
     enddo
!
     do j = 3,jm2-1
       dist(j)=fa1*(var(j)+var(j-1))-fa2*(var(j+1)+var(j-2))
     enddo
     dist(2)=fa1*(var(2)+var(1  ))-fa2*(var(3)+var(jm2  ))
     dist(1)=fa1*(var(1)+var(jm2))-fa2*(var(2)+var(jm2-1))
     dist(jm2+1)=dist(1)
     dist(jm2  )=fa1*(var(jm2)+var(jm2-1))-fa2*(var(1)+var(jm2-2))
!
     do j = 1,jm2+1
       ypast(j) = gglati(j) - dist(j)
       ynext(j) = gglati(j) + dist(j)
     enddo
!
     if( mass.eq.1 ) then
       do j = 1,jm2
         dyfact(j) = (ypast(j+1)-ypast(j)) / (ynext(j+1)-ynext(j))
       enddo
     endif
!
! advection all in y
!
     sc=gglati(jm2+1)-gglati(1)
     do n = 1,ncld
       past(1:jm2) = qq(1:jm2,k,n)
#ifdef NISLQ_NEW
       if ( iopt.eq.16 .or. iopt.eq.17 ) then
         call spline_cyclic_interface(jm2,gglati,past,jm2,ypast ,da  ,2,.true.,2)
         call spline_cyclic_interface(jm2,ynext ,da  ,jm2,gglati,next,2,.true.,2)
       else
#endif
         call cyclic_cell_ppm_intp(gglati,past,ypast,da,jm2,jm2,jm2,sc)
         if( mass.eq.1 ) da(1:jm2) = da(1:jm2) * dyfact(1:jm2)
         call cyclic_cell_ppm_intp(ynext,da,gglati,next,jm2,jm2,jm2,sc)
#ifdef NISLQ_NEW
       endif
#endif
       qq(1:jm2,k,n) = next(1:jm2)
     enddo
! 
   enddo
!
   return
   end subroutine cyclic_cell_massadvy
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine cyclic_cell_intpx(levs,imp,imf,qq)
!-------------------------------------------------------------------------------
!
! do  mass conserving interpolation from different grid at given latitude
!
! author: hann-ming henry juang 2008
!
!-------------------------------------------------------------------------------
   use constant, only : pi_
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer                          ::  levs, imp, imf, i,k,im
   real                             ::  two_pi,dxp,dxf,hfdxp,hfdxf,sc
   real   , dimension(lonfull,levs) ::  qq
   real   , dimension(lonfull)      ::  past,next
   real   , dimension(lonfull+1)    ::  xpast,xnext
!
   im = lonfull
! ..................................  
   if( imp.ne.imf ) then
! ..................................
     two_pi = 2.*pi_
     dxp = two_pi / imp
     dxf = two_pi / imf
     hfdxp = 0.5 * dxp
     hfdxf = 0.5 * dxf
!
     do i = 1,imp+1
       xpast(i) = (i-1) * dxp - hfdxp
     enddo
!
     do i = 1,imf+1
       xnext(i) = (i-1) * dxf - hfdxf
     enddo
!
     sc=two_pi
     do k = 1,levs
       do i = 1,imp
         past(i)=qq(i,k)
       enddo
       call cyclic_cell_ppm_intp(xpast,past,xnext,next,im,imp,imf,sc)
       do i = 1,imf
         qq(i,k)=next(i)
       enddo
     enddo
! .................       
   endif
! .................
   return
   end subroutine cyclic_cell_intpx
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine cyclic_cell_ppm_intp(pp,qq,pn,qn,lons,lonp,lonn,sc)
!-------------------------------------------------------------------------------
!
! mass conservation in cyclic bc interpolation: interpolate a group
! of grid point  coordiante call pp at interface with quantity qq at
! cell averaged to a group of new grid point coordinate call pn at
! interface with quantity qn at cell average with ppm spline.
! in horizontal with mass conservation is under the condition that
! pp(1)= pp(lons+1)=pn(lons+1)
!
! pp    location at interfac level as input
! qq    quantity at averaged-cell as input
! pn    location at interface of new grid structure as input
! qn    quantity at averaged-cell as output
! lons  numer of cells for dimension
! lonp  numer of cells for input
! lonn  numer of cells for output
! lev   number of vertical layers
! mono  monotonicity o:no, 1:yes
!
! author : henry.juang@noaa.gov
!
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer, parameter          ::  mono=1
   integer                     ::  lons,lonp,lonn             ,&
                                   ik,le,kstr,kend            ,&
                                   i,k, kl, kh, kk, kkl, kkh
   real                        ::  length,sc                  ,&
                                   dqi,dqimax,dqimin          ,&
                                   tl,tl2,tl3,qql,dql         ,&
                                   th,th2,th3,qqh,dqh         ,&
                                   dpp,dpq,c1,c2              ,&
                                   dpph,dppl,dpqh,dpql
   real   , dimension(lons)    ::  qq,qn
   real   , dimension(lons+1)  ::  pp,pn
   real   , dimension(3*lonp)  ::  locs,mass,hh,dqmono,qmi,qpi
!
!  length = pp(lonp+1) - pp(1)
   length=sc
!
! arrange input array cover output location with cyclic boundary condition
!
   locs(lonp+1:2*lonp) = pp(1:lonp)
   do i = 1,lonp
     locs(i) = locs(i+lonp) - length
     locs(i+2*lonp) = locs(i+lonp) + length
   enddo
!
   find_kstr : do i = 1,3*lonp
     if( pn(1).ge.locs(i) .and. pn(1).lt.locs(i+1) ) then
       kstr = i
       exit find_kstr
     else
       cycle find_kstr
     endif
   enddo find_kstr
   kstr=max(1,kstr)
!
   mass(lonp+1:2*lonp) = qq(1:lonp)
   do i = 1,lonp
     mass(i) = mass(i+lonp)
     mass(i+2*lonp) = mass(i+lonp)
   enddo
!
! prepare grid spacing
!
   do i = lonp+1,2*lonp
     hh(i) = locs(i+1)-locs(i)
   enddo
   do i = 1,lonp
     hh(i) = hh(i+lonp)
     hh(i+2*lonp) = hh(i+lonp)
   enddo
!
! prepare location with monotonic concerns
!
   do i = lonp+1,2*lonp
     dqi = 0.25*(mass(i+1)-mass(i-1))
     dqimax = max(mass(i-1),mass(i),mass(i+1)) - mass(i)
     dqimin = mass(i) - min(mass(i-1),mass(i),mass(i+1))
     dqmono(i) = sign( min( abs(dqi), dqimin, dqimax ), dqi)
   enddo
   do i = 1,lonp
     dqmono(i) = dqmono(i+lonp)
     dqmono(i+2*lonp) = dqmono(i+lonp)
   enddo
!
! compute value at interface with monotone
!
   do i = lonp+1,2*lonp
     qmi(i)=(mass(i-1)*hh(i)+mass(i)*hh(i-1))/(hh(i)+hh(i-1))                  &
             +(dqmono(i-1)-dqmono(i))/3.0
!    qmi(i)=(mass(i-1)*hh(i)+mass(i)*hh(i-1))/(hh(i)+hh(i-1))
   enddo
   qmi(2*lonp+1)=qmi(lonp+1)
   do i = lonp+1,2*lonp
     qpi(i)=qmi(i+1)
   enddo
!
! do less diffusive
!
   do i = lonp+1,2*lonp
     qmi(i)=mass(i)                                                            &
            -sign(min(abs(2.*dqmono(i)),abs(qmi(i)-mass(i))),2.*dqmono(i))
   enddo
   do i = lonp+1,2*lonp
     qpi(i)=mass(i)                                                            &
            +sign(min(abs(2.*dqmono(i)),abs(qpi(i)-mass(i))),2.*dqmono(i))
   enddo
!
! do monotonicity
!
   if( mono.eq.1 ) then
     do i = lonp+1,2*lonp
       c1=qpi(i)-mass(i)
       c2=mass(i)-qmi(i)
       if( c1*c2.le.0.0 ) then
         qmi(i)=mass(i)
         qpi(i)=mass(i)
       endif
     enddo
     do i = lonp+1,2*lonp
       c1=(qpi(i)-qmi(i))*(mass(i)-0.5*(qpi(i)+qmi(i)))
       c2=(qpi(i)-qmi(i))*(qpi(i)-qmi(i))/6.
       if( c1.gt.c2 ) then
         qmi(i)=3.*mass(i)-2.*qpi(i)
       else if( c1.lt.-c2 ) then
         qpi(i)=3.*mass(i)-2.*qmi(i)
       endif
     enddo
   endif
!
! extend array with cyclic condition
!
   do i = 1,lonp
     qmi(i)        = qmi(i+lonp)
     qmi(i+2*lonp) = qmi(i+lonp)
     qpi(i)        = qpi(i+lonp)
     qpi(i+2*lonp) = qpi(i+lonp)
   enddo
!
! start interpolation by integral of ppm spline
!
   kkl = kstr
   do i = 1,lonn
     kl = i
     kh = i + 1
! find kkh
     do kk = kkl+1,3*lonp
       if( pn(kh).lt.locs(kk) ) then
         kkh = kk-1
         go to 100
       endif
     enddo
100  continue
! mass interpolate
     if( kkh.eq.kkl ) then
!      print *,' condition 000000000000 ',kl,kh,kkl,kkh
!      print *,' pl ph ll lh ',pn(kl),pn(kh),locs(kkl),locs(kkl+1)
       tl=(pn(kl)-locs(kkl))/hh(kkl)
       tl2=tl*tl
       tl3=tl2*tl
       th=(pn(kh)-locs(kkl))/hh(kkl)
       th2=th*th
       th3=th2*th
       qqh=(th3-th2)*qpi(kkl)                                                  &
           +(th3-2.*th2+th)*qmi(kkl)+(-2.*th3+3.*th2)*mass(kkl)
       qql=(tl3-tl2)*qpi(kkl)                                                  &
           +(tl3-2.*tl2+tl)*qmi(kkl)+(-2.*tl3+3.*tl2)*mass(kkl)
       qn(i) = (qqh-qql)/(th-tl)
     else if( kkh.gt.kkl ) then
       tl=(pn(kl)-locs(kkl))/hh(kkl)
       tl2=tl*tl
       tl3=tl2*tl
       qql=(tl3-tl2)*qpi(kkl)                                                  &
           +(tl3-2.*tl2+tl)*qmi(kkl)+(-2.*tl3+3.*tl2)*mass(kkl)
       dql = mass(kkl)-qql
       th=(pn(kh)-locs(kkh))/hh(kkh)
       th2=th*th
       th3=th2*th
       dqh=(th3-th2)*qpi(kkh)                                                  &
           +(th3-2.*th2+th)*qmi(kkh)+(-2.*th3+3.*th2)*mass(kkh)
       dpp  = (1.-tl)*hh(kkl) + th*hh(kkh)
       dpq  = dql*hh(kkl) + dqh*hh(kkh)
       if( kkh-kkl.gt.1 ) then
!        print *,' condition 2222222222 ',kl,kh,kkl,kkh
!        print *,' pl ph ll lh ',pn(kl),pn(kh),locs(kkl),locs(kkh)
         do kk = kkl+1,kkh-1
           dpp = dpp + hh(kk)
           dpq = dpq + mass(kk)*hh(kk)
         enddo
       endif
       qn(i) = dpq / dpp
     else
       print *,' Error in cyclic_cell_ppm_intp location not found '
#ifdef MP
       call mpabort
#else
       call abort
#endif
     endif
! next one
     kkl = kkh
   enddo
!
   return
   end subroutine cyclic_cell_ppm_intp
!-------------------------------------------------------------------------------
!
!
#ifdef NISLQ_NEW
!-------------------------------------------------------------------------------
   subroutine vertical_cell_advect(lons,londim,levs,nvars,delt,ssi,wwi,qql,mass)
!-------------------------------------------------------------------------------
!
!  Abstract: Vertical Lagrangian advection for moisture and tracers
!
!  Histor:
!    2012-01-01  henry juang     initial version
!    
!-------------------------------------------------------------------------------
!
   implicit none
!
! passing variables
!
   integer                              , intent(in   )  ::  londim,levs,nvars,&
                                                             lons,mass
   real                                 , intent(in   )  ::  delt
   real   , dimension(londim,levs+1)    , intent(in   )  ::  ssi,wwi
   real   , dimension(londim,levs,nvars), intent(inout)  ::  qql
!
! local variables
!
   integer                         ::  i,k,n,isp
   real                            ::  sstmp, dpdt, check
   real   , dimension(levs)        ::  dsfact,dpd,dpi,dpa
   real   , dimension(levs+1)      ::  ssii,ssid,ssia
   real   , dimension(levs,nvars)  ::  rqmm,rqnn,rqda
!
! i-loop
!
   do i = 1,lons

     do k = 1,levs+1
       ssii(k) = ssi(i,k)
       ssid(k) = ssi(i,k)-wwi(i,k)*delt
       ssia(k) = ssi(i,k)+wwi(i,k)*delt
     enddo
!
! vertical remapping
!
     if ( iopt.eq.-1 ) then  ! original (bottom to top)
!
       if ( mass.eq.1 ) then
         do k = 1,levs
           dsfact(k) = (ssid(k)-ssid(k+1))/(ssia(k)-ssia(k+1))
         enddo
       endif
       do n = 1,nvars
         rqmm(1:levs,n) = qql(i,1:levs,n)
       enddo
       call vertical_cell_ppm_intp(ssii,rqmm,ssid,rqda,levs,nvars)
       if ( mass.eq.1 ) then
         do n = 1,nvars
           rqda(1:levs,n) = rqda(1:levs,n) * dsfact(1:levs)
         enddo
       endif
       call vertical_cell_ppm_intp(ssia,rqda,ssii,rqnn,levs,nvars)
       do n = 1,nvars
         qql(i,1:levs,n) = rqnn(1:levs,n)
       enddo

     elseif ( iopt.ge.0 .and. iopt.le.2 ) then  ! cam remap (top to bottom)

       do k = 1,levs
         dpi(levs+1-k) = ssii(k)-ssii(k+1)
         dpd(levs+1-k) = ssid(k)-ssid(k+1)
         dpa(levs+1-k) = ssia(k)-ssia(k+1)
       enddo
       do n = 1,nvars
         do k = 1,levs
           rqmm(k,n) = qql(i,levs+1-k,n)
         enddo
       enddo
       call remap1(rqmm,1,levs,nvars,dpi,dpd,iopt)
       call remap1(rqmm,1,levs,nvars,dpa,dpi,iopt)
       do n = 1,nvars
         do k = 1,levs
           qql(i,levs+1-k,n) = rqmm(k,n)
         enddo
       enddo
!
     elseif ( iopt.ge.3 ) then  ! kim remap
!
       if ( iopt.eq.3 ) isp = 2  ! ppm
       if ( iopt.eq.4 ) isp = 1  ! plm
       if ( iopt.eq.5 ) isp = 0  ! pcm
       if ( iopt.eq.6 .or. iopt.eq.16 ) isp = 2  ! psm
       if ( iopt.eq.7 .or. iopt.eq.17 ) isp = 4  ! qsm

       do n = 1,nvars
         dpi(1:levs) = qql(i,1:levs,n)
         if ( iopt.ge.6 ) then
           call spline_interface(levs,ssii,dpi,levs,ssid,dpd,isp,.true.,1)
           call spline_interface(levs,ssia,dpd,levs,ssii,dpa,isp,.true.,1)
         else
           call piecewise_interface(levs,ssii,dpi,levs,ssid,dpd,isp,.true.,1)
           call piecewise_interface(levs,ssia,dpd,levs,ssii,dpa,isp,.true.,1)
         endif
         qql(i,1:levs,n) = dpa(1:levs)
       enddo

     endif
!
   enddo
!
   end subroutine vertical_cell_advect
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine vertical_cell_ppm_intp(pp,qq,pn,qn,levs,nvars)
!-------------------------------------------------------------------------------
!
! mass conservation in vertical interpolation: interpolate a group
! of grid point  coordiante call pp at interface with quantity qq at
! cell averaged to a group of new grid point coordinate call pn at
! interface with quantity qn at cell average with ppm spline.
! in vertical with mass conservation is under the condition that
! pp(1)=pn(1), pp(levs+1)=pn(levs+1)
!
! pp    pressure at interfac level as input
! qq    quantity at layer as input
! pn    pressure at interface of new grid structure as input
! qn    quantity at layer as output
! levs  numer of verical layers
!
! author : henry.juang@noaa.gov
!
!-------------------------------------------------------------------------------
   implicit none
!
   integer, parameter :: mono=1
!
! passing variables
!
   integer                         , intent(in   )  ::  levs,nvars
   real   , dimension(levs+1)      , intent(in   )  ::  pp,pn
   real   , dimension(levs  ,nvars), intent(in   )  ::  qq
   real   , dimension(levs  ,nvars), intent(  out)  ::  qn
!
! local variables
!
   integer                           ::  i,k, kl, kh, kk, kkl, kkh,n
   real                              ::  massm,massc,massp,massbot,masstop,    &
                                         dqi,dqimax,dqimin,                    &
                                         tl,tl2,tl3,tlp,tlm,tlc,               &
                                         th,th2,th3,thp,thm,thc,               &
                                         dpp,dqq,c1,c2
   real   , dimension(levs)          ::  hh
   real   , dimension       (nvars)  ::  dql,dqh
   real   , dimension(levs  ,nvars)  ::  qmi,qpi,dqmono
!
      if( pp(1).ne.pn(1) .or. pp(levs+1).ne.pn(levs+1) ) then
        print *,' Error in vertical_cell_ppm_intp for domain values '
        print *,' i pp1 pn1 ppt pnt ',i,pp(1),pn(1),pp(levs+1),pn(levs+1)
#ifdef MP
        call mpabort
#else
        call abort
#endif
      endif
!
! prepare thickness for grid
!
      do k=1,levs
        hh(k) = pp(k+1)-pp(k)
      enddo
!
! prepare location with monotonic concerns
!
      do n=1,nvars
        massbot=(3.*hh(1)+hh(2))*qq(1,n)-2.*hh(1)*qq(2,n)
        massm = massbot/(hh(1)+hh(2))
        massc = qq(1  ,n)
        massp = qq(1+1,n)
        dqi = 0.25*(massp-massm)
        dqimax = max(massm,massc,massp) - massc
        dqimin = massc - min(massm,massc,massp)
        dqmono(1,n) = sign( min( abs(dqi), dqimin, dqimax ), dqi)
        do k=2,levs-1
          massp = qq(k+1,n)
          massc = qq(k  ,n)
          massm = qq(k-1,n)
          dqi = 0.25*(massp-massm)
          dqimax = max(massm,massc,massp) - massc
          dqimin = massc - min(massm,massc,massp)
          dqmono(k,n) = sign( min( abs(dqi), dqimin, dqimax ), dqi)
        enddo
        masstop=(3.*hh(levs)+hh(levs-1))*qq(levs,n) -2.*hh(levs)*qq(levs-1,n)
        massp = masstop/(hh(levs)+hh(levs-1))
        massc = qq(levs  ,n)
        massm = qq(levs-1,n)
        dqi = 0.25*(massp-massm)
        dqimax = max(massm,massc,massp) - massc
        dqimin = massc - min(massm,massc,massp)
        dqmono(levs,n) = sign( min( abs(dqi), dqimin, dqimax ), dqi)
!
! compute value at interface with momotone
!
        do k=2,levs
          qmi(k,n)=(qq(k-1,n)*hh(k)+qq(k,n)*hh(k-1))/(hh(k)+hh(k-1))           &
            +(dqmono(k-1,n)-dqmono(k,n))/3.0
        enddo
        do k=1,levs-1
          qpi(k,n)=qmi(k+1,n)
        enddo
        qmi(1,n)=qq(1,n)
        qpi(1,n)=qq(1,n)
        qmi(levs,n)=qq(levs,n)
        qpi(levs,n)=qq(levs,n)
      enddo
!
! do monotonicity
!
      if( mono.eq.1 ) then
        do n=1,nvars
        do k=1,levs
          c1=qpi(k,n)-qq(k,n)
          c2=qq(k,n)-qmi(k,n)
          if( c1*c2.le.0.0 ) then
            qmi(k,n)=qq(k,n)
            qpi(k,n)=qq(k,n)
          endif
        enddo
        do k=1,levs
          c1=(qpi(k,n)-qmi(k,n))*(qq(k,n)-0.5*(qpi(k,n)+qmi(k,n)))
          c2=(qpi(k,n)-qmi(k,n))*(qpi(k,n)-qmi(k,n))/6.
          if( c1.gt.c2 ) then
            qmi(k,n)=3.*qq(k,n)-2.*qpi(k,n)
          else if( c1.lt.-c2 ) then
            qpi(k,n)=3.*qq(k,n)-2.*qmi(k,n)
          endif
        enddo
        enddo
      endif
!
! start interpolation by integral of ppm spline
!
      kkl = 1
      tl=0
      do n=1,nvars
        dql(n)=0.0
      enddo

      do k=1,levs

        kl = k
        kh = k + 1
! find kkh
        do kk=kkl+1,levs+1
          if( pn(kh).ge.pp(kk) ) then
            kkh = kk-1
            go to 100
          endif
        enddo
        print *,' Error in vertical_cell_ppm_intp for no lev found '
        print *,' kh kl kkl',kh,kl,kkl
        print *,' pn ',(pn(kk),kk=1,levs+1)
        print *,' pp ',(pp(kk),kk=1,levs+1)
#ifdef MP
        call mpabort
#else
        call abort
#endif
 100    continue
        th=(pn(kh)-pp(kkh))/hh(kkh)
        th2=th*th
        th3=th2*th
        thp = th3-th2
        thm = th3-2.*th2+th
        thc = -2.*th3+3.*th2
        do n=1,nvars
          dqh(n)=thp*qpi(kkh,n)+thm*qmi(kkh,n)+thc*qq(kkh,n)
        enddo
! mass interpolate
        if( kkh.eq.kkl ) then
          do n=1,nvars
            qn(k,n) = (dqh(n)-dql(n))/(th-tl)
          enddo
        else if( kkh.gt.kkl ) then
          dpp  = (1.-tl)*hh(kkl) + th*hh(kkh)
          do kk=kkl+1,kkh-1
            dpp = dpp + hh(kk)
          enddo
          do n=1,nvars
            dql(n) = qq(kkl,n)-dql(n)
            dqq  = dql(n)*hh(kkl) + dqh(n)*hh(kkh)
            do kk=kkl+1,kkh-1
              dqq = dqq + qq(kk,n)*hh(kk)
            enddo
            qn(k,n) = dqq / dpp
          enddo
        else
          print *,' Error in vertical_cell_ppm_intp for lev messed up '
          print *,' i kh kl ',i,kh,kl
          print *,' pn ',(pn(kk),kk=1,levs+1)
          print *,' pp ',(pp(kk),kk=1,levs+1)
#ifdef MP
          call mpabort
#else
          call abort
#endif
        endif
! next one
        kkl = kkh
        tl  = th
        do n=1,nvars
          dql(n) = dqh(n)
        enddo

      enddo ! end of k loop
!
   return
   end subroutine vertical_cell_ppm_intp
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine remap1(Qdp,nx,nlev,qsize,dp1,dp2,vert_remap_q_alg)
!-------------------------------------------------------------------------------
  ! remap 1 field
  ! input:  Qdp   field to be remapped (NOTE: MASS, not MIXING RATIO)
  !         dp1   layer thickness (source)
  !         dp2   layer thickness (target)
  !
  ! output: remaped Qdp, conserving mass, monotone on Q=Qdp/dp
  !
  implicit none
  integer, intent(in) :: nx,nlev,qsize,vert_remap_q_alg
  real, intent(inout) :: Qdp(nx,nlev,qsize)
  real, intent(in) :: dp1(nx,nlev),dp2(nx,nlev)
  ! ========================
  ! Local Variables
  ! ========================

  real, dimension(nlev+1)    :: rhs,lower_diag,diag,upper_diag,q_diag,zgam,z1c,z2c,zv
  real, dimension(nlev)      :: h,Qcol,dy,za0,za1,za2,zarg,zhdp,dp_star,dp_np1
  real  :: f_xm,level1,level2,level3,level4,level5, &
           peaks_min,peaks_max,tmp_cal,xm,xm_d,zv1,zv2, &
           zero = 0,one = 1,tiny = 1e-12,qmax = 1d50
  integer :: zkr(nlev+1),filter_code(nlev),peaks,im1,im2,im3,ip1,ip2, &
             lt1,lt2,lt3,t0,t1,t2,t3,t4,tm,tp,ie,i,ilev,jk,k,q
  logical :: abort=.false.
!mskoo  integer  ::  vert_remap_q_alg = 0
!mskoo  call t_startf('remap1')

  if (vert_remap_q_alg == 1 .or. vert_remap_q_alg == 2) then
     call remap_Q_ppm(qdp,nx,nlev,qsize,dp1,dp2,vert_remap_q_alg)
!mskoo     call t_stopf('remap1')
     return
  endif

  do q=1,qsize
    do i=1,nx
      z1c(1)=0 ! source grid
      z2c(1)=0 ! target grid
      do k=1,nlev
         z1c(k+1)=z1c(k)+dp1(i,k)
         z2c(k+1)=z2c(k)+dp2(i,k)
      enddo

      zv(1)=0
      do k=1,nlev
        Qcol(k)=Qdp(i,k,q)*(z1c(k+1)-z1c(k)) !  input is mass
!mskoo        Qcol(k)=Qdp(i,k,q)!  *(z1c(k+1)-z1c(k)) input is mass
        zv(k+1) = zv(k)+Qcol(k)
      enddo

      if (ABS(z2c(nlev+1)-z1c(nlev+1)).GE.0.000001) then
        write(6,*) 'SURFACE PRESSURE IMPLIED BY ADVECTION SCHEME'
        write(6,*) 'NOT CORRESPONDING TO SURFACE PRESSURE IN    '
        write(6,*) 'DATA FOR MODEL LEVELS'
        write(6,*) 'PLEVMODEL=',z2c(nlev+1)
        write(6,*) 'PLEV     =',z1c(nlev+1)
        write(6,*) 'DIFF     =',z2c(nlev+1)-z1c(nlev+1)
        abort=.true.
      endif

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! quadratic splies with UK met office monotonicity constraints  !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      zkr  = 99
      ilev = 2
      zkr(1) = 1
      zkr(nlev+1) = nlev
      kloop: do k = 2,nlev
        do jk = ilev,nlev+1
          if (z1c(jk).ge.z2c(k)) then
            ilev      = jk
            zkr(k)   = jk-1
            cycle kloop
          endif
        enddo
      enddo kloop

      zgam  = (z2c(1:nlev+1)-z1c(zkr)) / (z1c(zkr+1)-z1c(zkr))
      zgam(1)      = 0.0
      zgam(nlev+1) = 1.0
      zhdp = z1c(2:nlev+1)-z1c(1:nlev)


      h = 1/zhdp
      zarg = Qcol * h
      rhs = 0
      lower_diag = 0
      diag = 0
      upper_diag = 0

      rhs(1)=3*zarg(1)
      rhs(2:nlev) = 3*(zarg(2:nlev)*h(2:nlev) + zarg(1:nlev-1)*h(1:nlev-1))
      rhs(nlev+1)=3*zarg(nlev)

      lower_diag(1)=1
      lower_diag(2:nlev) = h(1:nlev-1)
      lower_diag(nlev+1)=1

      diag(1)=2
      diag(2:nlev) = 2*(h(2:nlev) + h(1:nlev-1))
      diag(nlev+1)=2

      upper_diag(1)=1
      upper_diag(2:nlev) = h(2:nlev)
      upper_diag(nlev+1)=0

      q_diag(1)=-upper_diag(1)/diag(1)
      rhs(1)= rhs(1)/diag(1)

      do k=2,nlev+1
        tmp_cal    =  1/(diag(k)+lower_diag(k)*q_diag(k-1))
        q_diag(k) = -upper_diag(k)*tmp_cal
        rhs(k) =  (rhs(k)-lower_diag(k)*rhs(k-1))*tmp_cal
      enddo
      do k=nlev,1,-1
        rhs(k)=rhs(k)+q_diag(k)*rhs(k+1)
      enddo

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!  monotonicity modifications  !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      filter_code = 0
      dy(1:nlev-1) = zarg(2:nlev)-zarg(1:nlev-1)
      dy(nlev) = dy(nlev-1)

      dy = merge(zero, dy, abs(dy) < tiny )

      do k=1,nlev
        im1=MAX(1,k-1)
        im2=MAX(1,k-2)
        im3=MAX(1,k-3)
        ip1=MIN(nlev,k+1)
        t1 = merge(1,0,(zarg(k)-rhs(k))*(rhs(k)-zarg(im1)) >= 0)
        t2 = merge(1,0,dy(im2)*(rhs(k)-zarg(im1)) > 0 .AND. dy(im2)*dy(im3) > 0 &
             .AND. dy(k)*dy(ip1) > 0 .AND. dy(im2)*dy(k) < 0 )
        t3 = merge(1,0,ABS(rhs(k)-zarg(im1)) > ABS(rhs(k)-zarg(k)))

        filter_code(k) = merge(0,1,t1+t2 > 0)
        rhs(k) = (1-filter_code(k))*rhs(k)+filter_code(k)*(t3*zarg(k)+(1-t3)*zarg(im1))
        filter_code(im1) = MAX(filter_code(im1),filter_code(k))
      enddo

      rhs = merge(qmax,rhs,rhs > qmax)
      rhs = merge(zero,rhs,rhs < zero)

      za0 = rhs(1:nlev)
      za1 = -4*rhs(1:nlev) - 2*rhs(2:nlev+1) + 6*zarg
      za2 =  3*rhs(1:nlev) + 3*rhs(2:nlev+1) - 6*zarg

      dy(1:nlev) = rhs(2:nlev+1)-rhs(1:nlev)
      dy = merge(zero, dy, abs(dy) < tiny )

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! Compute the 3 quadratic spline coeffients {za0, za1, za2}                 !!
      !! knowing the quadratic spline parameters {rho_left,rho_right,zarg}         !!
      !! Zerroukat et.al., Q.J.R. Meteorol. Soc., Vol. 128, pp. 2801-2820 (2002).   !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      h = rhs(2:nlev+1)

      do k=1,nlev
        xm_d = merge(one,2*za2(k),abs(za2(k)) < tiny)
        xm = merge(zero,-za1(k)/xm_d, abs(za2(k)) < tiny)
        f_xm = za0(k) + za1(k)*xm + za2(k)*xm**2

        t1 = merge(1,0,ABS(za2(k)) > tiny)
        t2 = merge(1,0,xm <= zero .OR. xm >= 1)
        t3 = merge(1,0,za2(k) > zero)
        t4 = merge(1,0,za2(k) < zero)
        tm = merge(1,0,t1*((1-t2)+t3) .EQ. 2)
        tp = merge(1,0,t1*((1-t2)+(1-t3)+t4) .EQ. 3)

        peaks=0
        peaks = merge(-1,peaks,tm .EQ. 1)
        peaks = merge(+1,peaks,tp .EQ. 1)
        peaks_min = merge(f_xm,MIN(za0(k),za0(k)+za1(k)+za2(k)),tm .EQ. 1)
        peaks_max = merge(f_xm,MAX(za0(k),za0(k)+za1(k)+za2(k)),tp .EQ. 1)

        im1=MAX(1,k-1)
        im2=MAX(1,k-2)
        ip1=MIN(nlev,k+1)
        ip2=MIN(nlev,k+2)

        t1 = merge(abs(peaks),0,(dy(im2)*dy(im1) <= tiny) .OR. &
             (dy(ip1)*dy(ip2) <= tiny) .OR. (dy(im1)*dy(ip1) >= tiny) .OR. &
             (dy(im1)*float(peaks) <= tiny))

        filter_code(k) = merge(1,t1+(1-t1)*filter_code(k),(rhs(k) >= qmax) .OR. &
             (rhs(k) <= zero) .OR. (peaks_max > qmax) .OR. (peaks_min < tiny))

        if (filter_code(k) > 0) then
          level1 = rhs(k)
          level2 = (2*rhs(k)+h(k))/3
          level3 = 0.5*(rhs(k)+h(k))
          level4 = (1/3d0)*rhs(k)+2*(1/3d0)*h(k)
          level5 = h(k)

          t1 = merge(1,0,h(k) >= rhs(k))
          t2 = merge(1,0,zarg(k) <= level1 .OR.  zarg(k) >= level5)
          t3 = merge(1,0,zarg(k) >  level1 .AND. zarg(k) <  level2)
          t4 = merge(1,0,zarg(k) >  level4 .AND. zarg(k) <  level5)

          lt1 = t1*t2
          lt2 = t1*(1-t2+t3)
          lt3 = t1*(1-t2+1-t3+t4)

          za0(k) = merge(zarg(k),za0(k),lt1 .EQ. 1)
          za1(k) = merge(zero,za1(k),lt1 .EQ. 1)
          za2(k) = merge(zero,za2(k),lt1 .EQ. 1)

          za0(k) = merge(rhs(k),za0(k),lt2 .EQ. 2)
          za1(k) = merge(zero,za1(k),lt2 .EQ. 2)
          za2(k) = merge(3*(zarg(k)-rhs(k)),za2(k),lt2 .EQ. 2)

          za0(k) = merge(-2*h(k)+3*zarg(k),za0(k),lt3 .EQ. 3)
          za1(k) = merge(+6*h(k)-6*zarg(k),za1(k),lt3 .EQ. 3)
          za2(k) = merge(-3*h(k)+3*zarg(k),za2(k),lt3 .EQ. 3)

          t2 = merge(1,0,zarg(k) >= level1 .OR.  zarg(k) <= level5)
          t3 = merge(1,0,zarg(k) <  level1 .AND. zarg(k) >  level2)
          t4 = merge(1,0,zarg(k) <  level4 .AND. zarg(k) >  level5)

          lt1 = (1-t1)*t2
          lt2 = (1-t1)*(1-t2+t3)
          lt3 = (1-t1)*(1-t2+1-t3+t4)

          za0(k) = merge(zarg(k),za0(k),lt1 .EQ. 1)
          za1(k) = merge(zero,za1(k),lt1 .EQ. 1)
          za2(k) = merge(zero,za2(k),lt1 .EQ. 1)

          za0(k) = merge(rhs(k),za0(k),lt2 .EQ. 2)
          za1(k) = merge(zero,za1(k),lt2 .EQ. 2)
          za2(k) = merge(3*(zarg(k)-rhs(k)),za2(k),lt2 .EQ. 2)

          za0(k) = merge(-2*h(k)+3*zarg(k),za0(k),lt3 .EQ. 3)
          za1(k) = merge(+6*h(k)-6*zarg(k),za1(k),lt3 .EQ. 3)
          za2(k) = merge(-3*h(k)+3*zarg(k),za2(k),lt3 .EQ. 3)
        endif
      enddo

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !! start iteration from top to bottom of atmosphere !!
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      zv1 = 0
      do k=1,nlev
        if (zgam(k+1)>1d0) then
          WRITE(*,*) 'r not in [0:1]', zgam(k+1)
          abort=.true.
        endif
        zv2 = zv(zkr(k+1))+(za0(zkr(k+1))*zgam(k+1)+(za1(zkr(k+1))/2)*(zgam(k+1)**2)+ &
             (za2(zkr(k+1))/3)*(zgam(k+1)**3))*zhdp(zkr(k+1))
        Qdp(i,k,q) = (zv2 - zv1) / (z2c(k+1)-z2c(k) ) ! dont convert back to mixing ratio
!mskoo        Qdp(i,k,q) = (zv2 - zv1) ! / (z2c(k+1)-z2c(k) ) dont convert back to mixing ratio
        zv1 = zv2
      enddo
    enddo
  enddo ! q loop
  if (abort) call mpabort('Bad levels in remap1.  usually CFL violatioin')
!mskoo  call t_stopf('remap1')
!-------------------------------------------------------------------------------
   end subroutine remap1
!-------------------------------------------------------------------------------
!
!
!This uses the exact same model and reference grids and data as remap_Q, but it interpolates
!using PPM instead of splines.
!-------------------------------------------------------------------------------
   subroutine remap_Q_ppm(Qdp,nx,nlev,qsize,dp1,dp2,vert_remap_q_alg)
!-------------------------------------------------------------------------------
  ! remap 1 field
  ! input:  Qdp   field to be remapped (NOTE: MASS, not MIXING RATIO)
  !         dp1   layer thickness (source)
  !         dp2   layer thickness (target)
  !
  ! output: remaped Qdp, conserving mass
  !
  implicit none
  integer,intent(in) :: nx,nlev,qsize,vert_remap_q_alg
  real, intent(inout) :: Qdp(nx,nlev,qsize)
  real, intent(in) :: dp1(nx,nlev),dp2(nx,nlev)
  ! Local Variables
  integer, parameter :: gs = 2                              !Number of cells to place in the ghost region
  real, dimension(       nlev+2 ) :: pio    !Pressure at interfaces for old grid
  real, dimension(       nlev+1 ) :: pin    !Pressure at interfaces for new grid
  real, dimension(       nlev+1 ) :: masso  !Accumulate mass up to each interface
  real, dimension(  1-gs:nlev+gs) :: ao     !Tracer value on old grid
  real, dimension(  1-gs:nlev+gs) :: dpo    !change in pressure over a cell for old grid
  real, dimension(  1-gs:nlev+gs) :: dpn    !change in pressure over a cell for old grid
  real, dimension(3,     nlev   ) :: coefs  !PPM coefficients within each cell
  real, dimension(       nlev   ) :: z1, z2
  real :: ppmdx(10,0:nlev+1)  !grid spacings
  real :: mymass, massn1, massn2
  integer :: i, k, q, kk, kid(nlev)

!mskoo  call t_startf('remap_Q_ppm')
    do i = 1 , nx

      pin(1)=0
      pio(1)=0
      do k=1,nlev
         dpn(k)=dp2(i,k)
         dpo(k)=dp1(i,k)
         pin(k+1)=pin(k)+dpn(k)
         pio(k+1)=pio(k)+dpo(k)
      enddo



      pio(nlev+2) = pio(nlev+1) + 1.  !This is here to allow an entire block of k threads to run in the remapping phase.
                                      !It makes sure there's an old interface value below the domain that is larger.
      pin(nlev+1) = pio(nlev+1)       !The total mass in a column does not change.
                                      !Therefore, the pressure of that mass cannot either.
      !Fill in the ghost regions with mirrored values. if vert_remap_q_alg is defined, this is of no consequence.
      do k = 1 , gs
        dpo(1   -k) = dpo(       k)
        dpo(nlev+k) = dpo(nlev+1-k)
      enddo

      !Compute remapping intervals once for all tracers. Find the old grid cell index in which the
      !k-th new cell interface resides. Then integrate from the bottom of that old cell to the new
      !interface location. In practice, the grid never deforms past one cell, so the search can be
      !simplified by this. Also, the interval of integration is usually of magnitude close to zero
      !or close to dpo because of minimial deformation.
      !Numerous tests confirmed that the bottom and top of the grids match to machine precision, so
      !I set them equal to each other.
      do k = 1 , nlev
        kk = k  !Keep from an order n^2 search operation by assuming the old cell index is close.
        !Find the index of the old grid cell in which this new cell's bottom interface resides.
        do while ( pio(kk) <= pin(k+1) )
          kk = kk + 1
        enddo
        kk = kk - 1                   !kk is now the cell index we're integrating over.
        if (kk == nlev+1) kk = nlev   !This is to keep the indices in bounds.
                                      !Top bounds match anyway, so doesn't matter what coefficients are used
        kid(k) = kk                   !Save for reuse
        z1(k) = -0.5D0                !This remapping assumes we're starting from the left interface of an old grid cell
                                      !In fact, we're usually integrating very little or almost all of the cell in question
        z2(k) = ( pin(k+1) - ( pio(kk) + pio(kk+1) ) * 0.5 ) / dpo(kk)  !PPM interpolants are normalized to an independent
                                                                        !coordinate domain [-0.5,0.5].
      enddo

      !This turned out a big optimization, remembering that only parts of the PPM algorithm depends on the data, namely the
      !limiting. So anything that depends only on the grid is pre-computed outside the tracer loop.
      ppmdx(:,:) = compute_ppm_grids( vert_remap_q_alg,nlev,dpo)

      !From here, we loop over tracers for only those portions which depend on tracer data, which includes PPM limiting and
      !mass accumulation
      do q = 1 , qsize
        !Accumulate the old mass up to old grid cell interface locations to simplify integration
        !during remapping. Also, divide out the grid spacing so we're working with actual tracer
        !values and can conserve mass. The option for ifndef ZEROHORZ I believe is there to ensure
        !tracer consistency for an initially uniform field. I copied it from the old remap routine.
        masso(1) = 0.
        do k = 1 , nlev
          ao(k) = Qdp(i,k,q)
          masso(k+1) = masso(k) + ao(k)*dpo(k) !Accumulate the old mass. This will simplify the remapping
!mskoo          masso(k+1) = masso(k) + ao(k) !Accumulate the old mass. This will simplify the remapping
!mskoo          ao(k) = ao(k) / dpo(k)        !Divide out the old grid spacing because we want the tracer mixing ratio, not mass.
        enddo
        !Fill in ghost values. Ignored if vert_remap_q_alg == 2
        do k = 1 , gs
          ao(1   -k) = ao(       k)
          ao(nlev+k) = ao(nlev+1-k)
        enddo
        !Compute monotonic and conservative PPM reconstruction over every cell
        coefs(:,:) = compute_ppm( ao , ppmdx ,nlev, vert_remap_q_alg)
        !Compute tracer values on the new grid by integrating from the old cell bottom to the new
        !cell interface to form a new grid mass accumulation. Taking the difference between
        !accumulation at successive interfaces gives the mass inside each cell. Since Qdp is
        !supposed to hold the full mass this needs no normalization.
        massn1 = 0.
        do k = 1 , nlev
          kk = kid(k)
          massn2 = masso(kk) + integrate_parabola( coefs(:,kk) , z1(k) , z2(k) ) * dpo(kk)
!mskoo          Qdp(i,k,q) = massn2 - massn1
          Qdp(i,k,q) = ( massn2 - massn1 ) / dpn(k)
          massn1 = massn2
        enddo
      enddo
    enddo
!mskoo  call t_stopf('remap_Q_ppm')
!-------------------------------------------------------------------------------
   end subroutine remap_Q_ppm
!-------------------------------------------------------------------------------
!
!
!THis compute grid-based coefficients from Collela & Woodward 1984.
!-------------------------------------------------------------------------------
   function compute_ppm_grids( vert_remap_q_alg, nlev,dx )   result(rslt)
!-------------------------------------------------------------------------------
  implicit none
  integer, intent(in)  ::  vert_remap_q_alg,nlev
  real, intent(in) :: dx(-1:nlev+2)  !grid spacings
  real             :: rslt(10,0:nlev+1)  !grid spacings
  integer :: j
  integer :: indB, indE

  !Calculate grid-based coefficients for stage 1 of compute_ppm
  if (vert_remap_q_alg == 2) then
    indB = 2
    indE = nlev-1
  else
    indB = 0
    indE = nlev+1
  endif
  do j = indB , indE
    rslt( 1,j) = dx(j) / ( dx(j-1) + dx(j) + dx(j+1) )
    rslt( 2,j) = ( 2.*dx(j-1) + dx(j) ) / ( dx(j+1) + dx(j) )
    rslt( 3,j) = ( dx(j) + 2.*dx(j+1) ) / ( dx(j-1) + dx(j) )
  enddo

  !Caculate grid-based coefficients for stage 2 of compute_ppm
  if (vert_remap_q_alg == 2) then
    indB = 2
    indE = nlev-2
  else
    indB = 0
    indE = nlev
  endif
  do j = indB , indE
    rslt( 4,j) = dx(j) / ( dx(j) + dx(j+1) )
    rslt( 5,j) = 1. / sum( dx(j-1:j+2) )
    rslt( 6,j) = ( 2. * dx(j+1) * dx(j) ) / ( dx(j) + dx(j+1 ) )
    rslt( 7,j) = ( dx(j-1) + dx(j  ) ) / ( 2. * dx(j  ) + dx(j+1) )
    rslt( 8,j) = ( dx(j+2) + dx(j+1) ) / ( 2. * dx(j+1) + dx(j  ) )
    rslt( 9,j) = dx(j  ) * ( dx(j-1) + dx(j  ) ) / ( 2.*dx(j  ) +    dx(j+1) )
    rslt(10,j) = dx(j+1) * ( dx(j+1) + dx(j+2) ) / (    dx(j  ) + 2.*dx(j+1) )
  enddo
end function compute_ppm_grids

!This computes a limited parabolic interpolant using a net 5-cell stencil, but the stages of computation are broken up into 3 stages
function compute_ppm( a , dx ,nlev, vert_remap_q_alg)    result(coefs)
  implicit none
  integer, intent(in) ::  nlev,vert_remap_q_alg
  real, intent(in) :: a    (    -1:nlev+2)  !Cell-mean values
  real, intent(in) :: dx   (10,  0:nlev+1)  !grid spacings
  real ::             coefs(0:2,   nlev  )  !PPM coefficients (for parabola)
  real :: ai (0:nlev  )                     !fourth-order accurate, then limited interface values
  real :: dma(0:nlev+1)                     !An expression from Collela's '84 publication
  real :: da                                !Ditto
  ! Hold expressions based on the grid (which are cumbersome).
  real :: dx1, dx2, dx3, dx4, dx5, dx6, dx7, dx8, dx9, dx10
  real :: al, ar                            !Left and right interface values for cell-local limiting
  integer :: j
  integer :: indB, indE

  ! Stage 1: Compute dma for each cell, allowing a 1-cell ghost stencil below and above the domain
  if (vert_remap_q_alg == 2) then
    indB = 2
    indE = nlev-1
  else
    indB = 0
    indE = nlev+1
  endif
  do j = indB , indE
    da = dx(1,j) * ( dx(2,j) * ( a(j+1) - a(j) ) + dx(3,j) * ( a(j) - a(j-1) ) )
    dma(j) = minval( (/ abs(da) , 2. * abs( a(j) - a(j-1) ) , 2. * abs( a(j+1) - a(j) ) /) ) * sign(1.D0,da)
    if ( ( a(j+1) - a(j) ) * ( a(j) - a(j-1) ) <= 0. ) dma(j) = 0.
  enddo

  ! Stage 2: Compute ai for each cell interface in the physical domain (dimension nlev+1)
  if (vert_remap_q_alg == 2) then
    indB = 2
    indE = nlev-2
  else
    indB = 0
    indE = nlev
  endif
  do j = indB , indE
    ai(j) = a(j) + dx(4,j) * ( a(j+1) - a(j) ) + dx(5,j) * ( dx(6,j) * ( dx(7,j) - dx(8,j) ) &
         * ( a(j+1) - a(j) ) - dx(9,j) * dma(j+1) + dx(10,j) * dma(j) )
  enddo

  ! Stage 3: Compute limited PPM interpolant over each cell in the physical domain
  ! (dimension nlev) using ai on either side and ao within the cell.
  if (vert_remap_q_alg == 2) then
    indB = 3
    indE = nlev-2
  else
    indB = 1
    indE = nlev
  endif
  do j = indB , indE
    al = ai(j-1)
    ar = ai(j  )
    if ( (ar - a(j)) * (a(j) - al) <= 0. ) then
      al = a(j)
      ar = a(j)
    endif
    if ( (ar - al) * (a(j) - (al + ar)/2.) >  (ar - al)**2/6. ) al = 3.*a(j) - 2. * ar
    if ( (ar - al) * (a(j) - (al + ar)/2.) < -(ar - al)**2/6. ) ar = 3.*a(j) - 2. * al
    !Computed these coefficients from the edge values and cell mean in Maple. Assumes normalized coordinates: xi=(x-x0)/dx
    coefs(0,j) = 1.5 * a(j) - ( al + ar ) / 4.
    coefs(1,j) = ar - al
    coefs(2,j) = -6. * a(j) + 3. * ( al + ar )
  enddo

  !If we're not using a mirrored boundary condition, then make the two cells bordering the top and bottom
  !material boundaries piecewise constant. Zeroing out the first and second moments, and setting the zeroth
  !moment to the cell mean is sufficient to maintain conservation.
  if (vert_remap_q_alg == 2) then
    coefs(0,1:2) = a(1:2)
    coefs(1:2,1:2) = 0.
    coefs(0,nlev-1:nlev) = a(nlev-1:nlev)
    coefs(1:2,nlev-1:nlev) = 0.D0
  endif
!
   end function compute_ppm
!-------------------------------------------------------------------------------
!
!
!Simple function computes the definite integral of a parabola in normalized coordinates, xi=(x-x0)/dx,
!given two bounds. Make sure this gets inlined during compilation.
!-------------------------------------------------------------------------------
   function integrate_parabola( a , x1 , x2 )    result(mass)
!-------------------------------------------------------------------------------
  implicit none
!
  real, intent(in) :: a(0:2)  !Coefficients of the parabola
  real, intent(in) :: x1      !lower domain bound for integration
  real, intent(in) :: x2      !upper domain bound for integration
  real             :: mass
  mass = a(0) * (x2 - x1) + a(1) * (x2 ** 2 - x1 ** 2) / 0.2D1 + a(2) * (x2 ** 3 - x1 ** 3) / 0.3D1
!
   end function integrate_parabola
!-------------------------------------------------------------------------------
#else /* ~NISLQ_NEW */
!-------------------------------------------------------------------------------
   subroutine vertical_cell_advect(lons,londim,levs,ncld,deltim,               &
                                   ppi,wwi,qql,mass)
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer                               ::  lons,londim,levs,ncld,i,k,n,mass
   real                                  ::  deltim
   real   , dimension(londim,levs+1)     ::  ppi,wwi
   real   , dimension(londim,levs,ncld)  ::  qql
   real   , dimension(lons,levs)         ::  dsfact,rqmm,rqnn,rqda
   real   , dimension(lons,levs+1)       ::  ppii,ppid,ppia
!
   do k = 1,levs+1
     do i = 1,lons
       ppii(i,k)=ppi(i,k)
       ppid(i,k)=ppi(i,k)-wwi(i,k)*deltim
       ppia(i,k)=ppi(i,k)+wwi(i,k)*deltim
     enddo
   enddo
!
   if( mass.eq.1) then
     do k = 1,levs
       do i = 1,lons
         dsfact(i,k)=(ppid(i,k)-ppid(i,k+1))/(ppia(i,k)-ppia(i,k+1))
       enddo
     enddo
   endif
!
   do n = 1,ncld                              !hmhj nisl
     do k = 1,levs
       do i = 1,lons
         rqmm(i,k) = qql(i,k,n)
       enddo
     enddo
     call vertical_cell_ppm_intp(ppii,rqmm,ppid,rqda,lons,levs)
     if( mass.eq.1 ) then
       do k = 1,levs
         do i = 1,lons
           rqda(i,k) = rqda(i,k) * dsfact(i,k)
         enddo
       enddo
     endif
     call vertical_cell_ppm_intp(ppia,rqda,ppii,rqnn,lons,levs)
     do k = 1,levs
       do i = 1,lons
         qql(i,k,n)=rqnn(i,k)
       enddo
     enddo
   enddo
!
   return
   end subroutine vertical_cell_advect
!-------------------------------------------------------------------------------
! 
!-------------------------------------------------------------------------------
   subroutine vertical_cell_ppm_intp(pp,qq,pn,qn,lons,levs)
!-------------------------------------------------------------------------------
!
! mass conservation in vertical interpolation: interpolate a group
! of grid point  coordiante call pp at interface with quantity qq at
! cell averaged to a group of new grid point coordinate call pn at
! interface with quantity qn at cell average with ppm spline.
! in vertical with mass conservation is under the condition that
! pp(1)=pn(1), pp(lev+1)=pn(lev+1)
!
! pp    pressure at interfac level as input
! qq    quantity at layer as input
! pn    pressure at interface of new grid structure as input
! qn    quantity at layer as output
! lev  numer of verical layers
!
! author : henry.juang@noaa.gov
!
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer                           ::  lons,levs                            ,&
                                         i,k, kl, kh, kk, kkl, kkh
   real   , dimension(lons,levs)   ::  qq,qn
   real   , dimension(lons,levs+1) ::  pp,pn
!
   real                              ::  dqi,dqimax,dqimin    ,&
                                         tl,tl2,tl3,qql,dql   ,&
                                         th,th2,th3,qqh,dqh   ,&
                                         dpp,dpq,c1,c2
   real   , dimension(levs)          ::  hh,dqmono,qmi,qpi
   real   , dimension(0:levs+1)      ::  mass
   integer, parameter                ::  mono=1
!
   do i = 1,lons
!
     if( pp(i,1).ne.pn(i,1) .or. pp(i,levs+1).ne.pn(i,levs+1) ) then
       print *,' Error in vertical_cell_ppm_intp for domain values '
       print *,' i,pp1 pn1 ppt pnt ',i,pp(i,1),pn(i,1),pp(i,levs+1),pn(i,levs+1)
#ifdef MP
       call mpabort
#else
       call abort
#endif
     endif
!
! prepare thickness for uniform grid
!
     do k = 1,levs
       hh(k) = pp(i,k+1)-pp(i,k)         ! (top to bottm) ??
     enddo
!
! prepare location with monotonic concerns
!
     mass(1:levs)=qq(i,1:levs)
!
     mass(0)=(3.*hh(1)+hh(2))*mass(1)-2.*hh(1)*mass(2)
     mass(0)=mass(0)/(hh(1)+hh(2))
     mass(levs+1)=(3.*hh(levs)+hh(levs-1))*mass(levs)-2.*hh(levs)*mass(levs-1)
     mass(levs+1)=mass(levs+1)/(hh(levs)+hh(levs-1))
     do k = 1,levs
       dqi = 0.25*(mass(k+1)-mass(k-1))
       dqimax = max(mass(k-1),mass(k),mass(k+1)) - mass(k)
       dqimin = mass(k) - min(mass(k-1),mass(k),mass(k+1))
       dqmono(k) = sign( min( abs(dqimin), dqimin, dqimax ), dqi)
     enddo
!
! compute value at interface with momotone
!
     do k = 2,levs
       qmi(k)=(mass(k-1)*hh(k)+mass(k)*hh(k-1))/(hh(k)+hh(k-1))                &
             +(dqmono(k-1)-dqmono(k))/3.0
     enddo
     do k = 1,levs-1
       qpi(k)=qmi(k+1)
     enddo
     qmi(1)=mass(1)
     qpi(1)=mass(1)
     qmi(levs)=mass(levs)
     qpi(levs)=mass(levs)
!
! do monotonicity
!
     if( mono.eq.1 ) then
       do k = 1,levs
         c1=qpi(k)-mass(k)
         c2=mass(k)-qmi(k)
         if( c1*c2.le.0.0 ) then
           qmi(k)=mass(k)
           qpi(k)=mass(k)
         endif
       enddo
       do k = 1,levs
         c1=(qpi(k)-qmi(k))*(mass(k)-0.5*(qpi(k)+qmi(k)))
         c2=(qpi(k)-qmi(k))*(qpi(k)-qmi(k))/6.
         if( c1.gt.c2 ) then
           qmi(k)=3.*mass(k)-2.*qpi(k)
         else if( c1.lt.-c2 ) then
           qpi(k)=3.*mass(k)-2.*qmi(k)
         endif
       enddo
     endif
!
! start interpolation by integral of ppm spline
!
     kkl = 1
     do k = 1,levs
       kl = k
       kh = k + 1
! find kkh
       do kk = kkl+1,levs+1
!         if( pn(i,kh).ge.pp(i,kk) ) then
         if( pn(i,kh).le.pp(i,kk) ) then        ! top to bottom (hwangso)
           kkh = kk-1
           go to 100
         endif
       enddo
! mass interpolate
 100   if( kkh.eq.kkl ) then
         tl=(pn(i,kl)-pp(i,kkl))/hh(kkl)
         tl2=tl*tl
         tl3=tl2*tl
         th=(pn(i,kh)-pp(i,kkl))/hh(kkl)
         th2=th*th
         th3=th2*th
         qqh=(th3-th2)*qpi(kkl)+(th3-2.*th2+th)*qmi(kkl)   &
             +(-2.*th3+3.*th2)*mass(kkl)
         qql=(tl3-tl2)*qpi(kkl)+(tl3-2.*tl2+tl)*qmi(kkl)   &
             +(-2.*tl3+3.*tl2)*mass(kkl)
         qn(i,k) = (qqh-qql)/(th-tl)
       else if( kkh.gt.kkl ) then
         tl=(pn(i,kl)-pp(i,kkl))/hh(kkl)
         tl2=tl*tl
         tl3=tl2*tl
         qql=(tl3-tl2)*qpi(kkl)+(tl3-2.*tl2+tl)*qmi(kkl)   &
             +(-2.*tl3+3.*tl2)*mass(kkl)
         dql = qq(i,kkl)-qql
         th=(pn(i,kh)-pp(i,kkh))/hh(kkh)
         th2=th*th
         th3=th2*th
         dqh=(th3-th2)*qpi(kkh)+(th3-2.*th2+th)*qmi(kkh)   &
             +(-2.*th3+3.*th2)*mass(kkh)
         dpp= (1.-tl)*hh(kkl) + th*hh(kkh)
         dpq= dql*hh(kkl) + dqh*hh(kkh)
         if( kkh-kkl.gt.1 ) then
           do kk=kkl+1,kkh-1
             dpp = dpp + hh(kk)
             dpq = dpq + qq(i,kk)*hh(kk)
           enddo
         endif
         qn(i,k) = dpq / dpp
       else
         print *,' Error in vertical_cell_ppm_intp for no lev found '
         print *,' i kh kl ',i,kh,kl
         print *,' pn ',(pn(i,kk),kk=1,levs+1)
         print *,' pp ',(pp(i,kk),kk=1,levs+1)
#ifdef MP
         call mpabort
#else
         call abort
#endif
       endif
! next one
       kkl = kkh
     enddo
!
   enddo
!
   return
   end subroutine vertical_cell_ppm_intp
!-------------------------------------------------------------------------------
#endif /* NISLQ_NEW end */
!
!-------------------------------------------------------------------------------
   end module nislq
!-------------------------------------------------------------------------------
#endif /* ~KIM end */
