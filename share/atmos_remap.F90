!-------------------------------------------------------------------------------
   module atmos_remap
!-------------------------------------------------------------------------------
!
!  abstract :  
!
!  history log :
!    201?-??-??  junghan kim    initial setup
!    2017-02-13  junghan kim    code cleanup
!
!  structure : 
!
!-------------------------------------------------------------------------------
   use kinds, only : i4, r4, r8, l4
!-------------------------------------------------------------------------------
   implicit none
!
   private
!
   real(r8), parameter :: rd    = 287.04d0
   real(r8), parameter :: grav  = 9.80616d0
   real(r8), parameter :: ginv  = 1.d0/grav
   real(r8), parameter :: alpha = 0.0065d0*rd*ginv
!
   public :: remap_vertical
!
   contains
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   subroutine remap_vertical(npts, nlev_i, nlev_o, pres_i, pres_o, varflg, ps, ts, topo, dat_i, dat_o)
!-------------------------------------------------------------------------------
   implicit none
!
   integer,                           intent(in   ) :: npts, nlev_i, nlev_o
   real(r8), dimension(npts, nlev_i), intent(in   ) :: pres_i
   real(r8), dimension(nlev_o),       intent(in   ) :: pres_o
   integer,                           intent(in   ) :: varflg
   real(r8), dimension(npts),         intent(in   ) :: ps
   real(r8), dimension(npts),         intent(in   ) :: ts
   real(r8), dimension(npts),         intent(in   ) :: topo
   real(r8), dimension(npts, nlev_i), intent(in   ) :: dat_i
   real(r8), dimension(npts, nlev_o), intent(  out) :: dat_o
! local variables
   real(r8), dimension(nlevi) :: plevi
   real(r8), dimension(npts) :: phis
   real(r8) :: a2ln, a1
   real(r8) :: tstar, hgt, alnp, t0, tplat, tprime0, alpha, alph, psfcmb
   integer :: j, k, kp, kpi
   logical :: valid_ps, vert_interp
!
!-------------------------------------------------------------------------------
! Begins
!-------------------------------------------------------------------------------
!
! Statement function for double LOG interpolation on
! pressure surfaces presumes pressure is in mb
! --------------------------------------------------
!
   a2ln(a1) = log(log(a1+2.72d0))
   phis = topo/ginv
!
! Main loop
! ---------
!
   ma : do j = 1,npts
!
! In case surface pressure is undefined,vertical
! vertical interpolation is impossible
! -----------------------------------------------
!
     valid_ps = .true.
     if (ps(j)==spvl) then
       do k = 1,nlev_o
         dato(j,k) = spvl
       enddo
       valid_ps = .false.
     endif
!
! Only in case surface pressure is valid,interpolation or
! extrapolation is computed
! --------------------------------------------------------
!
     if (valid_ps) then
!
! Get pressure values for hybrid surfaces for current spatial point
! -----------------------------------------------------------------
!
       do k = 1,nlevi
         plevi(k) = pres_i(j,k)*0.01d0
       enddo
!
! Perfoms vertical interpolation/extrapolation
! --------------------------------------------
!
       vt : do k = 1,nlev_o
!
         vert_interp = .true.
!
! Output pressure is lower than minimum input pressure
! ----------------------------------------------------
!
         if (pres_o(k)<=plevi(1)) then ! higher than top pres
!
           kp = 1
           vert_interp = .true.
!
! Output pressure is larger than maximum input pressure
! -----------------------------------------------------
!
         elseif (pres_o(k)>plevi(nlev_i)) then ! lower than bottom pres
    
           if (kxtrp==0) then ! no extrapolation
        
             dato(j,k) = spvl
             vert_interp = .false.
        
           else ! extrapolation
           !
             if (varflg>0) then ! for temperature
          
               psmb = ps(j)*0.01d0
               tstar = dati(j,nlev_i)*(1.d0+alpha*(psmb/plevi(nlev_i)-1.d0))
               hgt = topo(j)
            
               if (hgt<2000.d0) then
                 alnp = alpha*log(pres_o(k)/psmb)
               else
                 t0 = tstar+0.0065d0*hgt
                 tplat = min(t0,298.d0)
                 if (hgt<=2500.d0) then
                   tprime0 = 0.002d0*((2500.d0-hgt)*t0+(hgt-2000.d0)*tplat)
                 else
                   tprime0 = tplat
                 endif
                 if (tprime0<tstar) then
                   alnp = 0.d0
                 else
                   alnp = rd*(tprime0-tstar)/phis(j)*log(pres_o(k)/psmb)
                 endif
               endif
            
               dato(j,k) = tstar*(1.d0+alnp+.5d0*alnp**2+1.d0/6.d0*alnp**3)
               vert_interp = .false.
          
             elseif (varflg<0) then ! for geopotential
          
               psmb = ps(j)*0.01d0
               hgt = topo(j)
               tstar = tbot(j)*(1.d0+alpha*(psmb/plevi(nlev_i)-1.d0))
               t0 = tstar+0.0065d0*hgt
            
               if (tstar<=290.5d0.and.t0>290.5d0) then
                 alph = rd/phis(j)*(290.5d0-tstar)
               elseif (tstar>290.5d0.and.t0>290.5d0) then
                 alph = 0
                 tstar = 0.5d0*(290.5d0+tstar)
               else
                 alph = alpha
               endif
            
               if (tstar<255.d0) then
                 tstar = 0.5d0*(tstar+255.d0)
               endif
               alnp = alph*log(pres_o(k)/psmb)
               dato(j,k) = hgt-rd*tstar*ginv*log(pres_o(k)/psmb)*  (1.d0+.5d0*alnp+1.d0/6.d0*alnp**2)
               vert_interp = .false.
          
             else ! for any other variables
          
               dato(j,k) = dati(j,nlev_i)
               vert_interp = .false.
          
             endif ! varflg
        
           endif ! kxtrp
    
         elseif (pres_o(k)>= plevi(nlev_i-1)) then
    
           kp = nlev_i-1
           vert_interp = .true.
    
         else
    
           kp = 0
           do
             kp = kp+1
             if (pres_o(k)<=plevi(kp+1)) then
               vert_interp = .true.
               exit
             endif
             if (kp>nlev_i) then
             write(6,'(a,2i5) ') 'KP>klevi. kp,klevi = ',kp,nlev_i
             endif
           enddo
    
         endif
!  
         if (vert_interp) then
           if (intyp==1) then
             dato(j,k) = dati(j,kp)+(dati(j,kp+1)-dati(j,kp))*(pres_o(k)-plevi(kp))/(plevi(kp+1)-plevi(kp))
           elseif (intyp==2) then
             dato(j,k) = dati(j,kp)+(dati(j,kp+1)-dati(j,kp))*log(pres_o(k)/plevi(kp))/log(plevi(kp+1)/plevi(kp))
           elseif (intyp==3) then
             dato(j,k) = dati(j,kp)+(dati(j,kp+1)-dati(j,kp))*(a2ln(pres_o(k))-a2ln(plevi(kp)))/(a2ln(plevi(kp+1))-a2ln(plevi(kp)))
           endif
         endif
!
       enddo vt
!
     endif ! valid_ps
!
   enddo ma
!
   return
   end subroutine remap_vertical
!-------------------------------------------------------------------------------
!
!
!-------------------------------------------------------------------------------
   end module atmos_remap
