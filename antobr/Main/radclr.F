
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine radclr(coszrs  ,trayoslp,pflx    ,abh2o   ,abo3    ,
     $                  abco2   ,abo2    ,uth2o   ,uto3    ,utco2   ,
     $                  uto2    ,   tauxar_mix_css, tauasc_mix_css    , 
     $                  gtota_mix_css, ftota_mix_css  ,
     $                  nloop   ,is      ,ie      ,rdir    ,rdif    ,
     $                  tdir    ,tdif    ,explay  ,exptdn  ,rdndif  ,
     $                  tottrn  )
C-----------------------------------------------------------------------
C
C Delta-Eddington solution for special clear sky computation
C
C Computes total reflectivities and transmissivities for two atmospheric
C layers: an overlying purely ozone absorbing layer, and the rest of the
C column below.
C
C For more details , see Briegleb, Bruce P., 1992: Delta-Eddington
C Approximation for Solar Radiation in the NCAR Community Climate Model,
C Journal of Geophysical Research, Vol 97, D7, pp7603-7612).
C
C---------------------------Code history--------------------------------
C
C Original version:  B. Briegleb
C Standardized:      J. Rosinski, June 1992
C Reviewed:          J. Kiehl, B. Briegleb, August 1992
C
C-----------------------------------------------------------------------
c
c $Id: radclr.F,v 1.2 1995/03/17 18:54:08 ccm2 Exp $
c $Author: ccm2 $
c
      implicit none
c
C------------------------------Parameters-------------------------------
c
#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
c
      include 'Commons/parrad.cb'
C
C-----------------------------------------------------------------------
C
C Minimum total transmission below which no layer computation are done:
C
      real*8   trmin,          ! Minimum total transmission allowed
     $      wray,           ! Rayleigh single scatter albedo
     $      gray,           ! Rayleigh asymetry parameter
     $      fray            ! Rayleigh forward scattered fraction
      parameter (trmin = 1.e-3,
     $           wray = 0.999999,
     $           gray = 0.0,
     $           fray = 0.1)
C
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real*8  coszrs(plond),         ! Cosine zenith angle
     $     trayoslp,              ! Tray/sslp
     $     pflx(plond,0:plevp),   ! Interface pressure
     $     abh2o,                 ! Absorption coefficiant for h2o
     $     abo3 ,                 ! Absorption coefficiant for o3
     $     abco2,                 ! Absorption coefficiant for co2
     $     abo2 ,                 ! Absorption coefficiant for o2
     $     uth2o(plond),          ! Total column absorber amount of h2o
     $     uto3(plond),           ! Total column absorber amount of  o3
     $     utco2(plond),          ! Total column absorber amount of co2
     $     uto2(plond)            ! Total column absorber amount of  o2

      real*8  tauxar_mix_css(plond),         ! Total column aerosol extinction
     $     tauasc_mix_css (plond),          ! Aerosol single scattering albedo
     $     gtota_mix_css(plond),           ! Aerosol asymmetry parameter
     $     ftota_mix_css(plond)            ! Aerosol forward scattering fraction
      integer nloop,              ! Number of loops (1 or 2)
     $        is(2),              ! Starting index for 1 or 2 loops
     $        ie(2)               ! Ending index for 1 or 2 loops
C
C Input/Output arguments
C
C Following variables are defined for each layer; note, we use layer 0
C to refer to the entire atmospheric column:
C
      real*8  rdir(plond,0:plev),    ! Layer reflectivity to direct rad
     $     rdif(plond,0:plev),    ! Layer refflectivity to diffuse rad
     $     tdir(plond,0:plev),    ! Layer transmission to direct rad
     $     tdif(plond,0:plev),    ! Layer transmission to diffuse rad
     $     explay(plond,0:plev)   ! Solar beam exp transmn for layer
C
C Note that the following variables are defined on interfaces, with
C the index k referring to the top interface of the kth layer:
C exptdn,rdndif,tottrn; for example, tottrn(k=5) refers to the total
C transmission to the top interface of the 5th layer.
C
      real*8  exptdn(plond,0:plevp), ! Solar beam exp down transmn from top
     $     rdndif(plond,0:plevp), ! Added dif ref for layers above
     $     tottrn(plond,0:plevp)  ! Total transmission for layers above
C
c      external  resetr,     ! Resets array elements to zero
c     $          whenfgt     ! Collect indices for greater than condition
C
C---------------------------Local variables-----------------------------
C
      integer i,            ! Longitude index
     $        k,            ! Level index
     $        nn,           ! Index of longitude loops (max=nloop)
     $        ii,           ! Longitude index
     $        nval,         ! Number of long values satisfying criteria
     $        index(plond)  ! Array of longitude indices
C
      real*8  taugab(plond),   ! Total column gas absorption optical depth
     $     tauray(plond),   ! Column rayleigh optical depth
     $     tautot       ,   ! Total column optical depth
     $       wtot       ,   ! Total column single scatter albedo
     $       gtot       ,   ! Total column asymmetry parameter
     $       ftot           ! Total column forward scatter fraction
      real*8  ts,            ! Column scaled extinction optical depth
     $       ws,            ! Column scaled single scattering albedo
     $       gs             ! Column scaled asymmetry parameter
      real*8  rdenom,          ! Mulitiple scattering term
     $     rdirexp,         ! Layer direct ref times exp transmission
     $     tdnmexp          ! Total transmission minus exp transmission
C
C---------------------------Statement functions-------------------------
C
C Statement functions for delta-Eddington solution; for detailed
C explanation of individual terms, see the routine 'radded'.
C
      real*8  alpha,gamma,el,taus,omgs,asys,u,n,lm,ne
      real*8  w,uu,g,e,f,t,et
C
C Intermediate terms for delta-Eddington solution
C
      real*8  alp,gam,ue,arg,extins,amg,apg
C
      alpha(w,uu,g,e) = .75*w*uu*((1. + g*(1-w))/(1. - e*e*uu*uu))
      gamma(w,uu,g,e) = .50*w*((3.*g*(1.-w)*uu*uu + 1.)/(1.-e*e*uu*uu))
      el(w,g)         = dsqrt(3.*(1-w)*(1. - w*g))
      taus(w,f,t)     = (1. - w*f)*t
      omgs(w,f)       = (1. - f)*w/(1. - w*f)
      asys(g,f)       = (g - f)/(1. - f)
      u(w,g,e)        = 1.5*(1. - w*g)/e
      n(uu,et)        = ((uu+1.)*(uu+1.)/et ) - ((uu-1.)*(uu-1.)*et)
C
C-----------------------------------------------------------------------
C
C Initialize all total transmimission values to 0, so that nighttime
C values from previous computations are not used:
C
      call resetr(tottrn,plond*2,0.d0)
c      print*,'dans radclr', maxval(tottrn)
      do k=0,plevp
      do i=1,plond
         tottrn(i,k)=0.
      enddo
      enddo
C
C Compute total direct beam transmission, total transmission, and
C reflectivity for diffuse radiation (from below) for all layers
C above each interface by starting from the top and adding layers
C down:
C
C The top layer is assumed to be a purely absorbing ozone layer, and
C that the mean diffusivity for diffuse transmission is 1.66:
C
      do nn=1,nloop
         do i=is(nn),ie(nn)
C
            taugab(i) = abo3*uto3(i)
C
C Limit argument of exponential to 25, in case coszrs is very small:
C
            arg         = dmin1(taugab(i)/coszrs(i),25.d0)
            explay(i,0) = dexp(-arg)
            tdir(i,0)   = explay(i,0)
C
C Same limit for diffuse transmission:
C
            arg         = dmin1(1.66*taugab(i),25.d0)
            tdif(i,0)   = dexp(-arg)
C
            rdir(i,0)   = 0.0
            rdif(i,0)   = 0.0
C
C Initialize top interface of extra layer:
C
            exptdn(i,0) =   1.0
            rdndif(i,0) =   0.0
            tottrn(i,0) =   1.0
C
            rdndif(i,1) = rdif(i,0)
            tottrn(i,1) = tdir(i,0)
C
         end do
      end do
C
C Now, complete the rest of the column; if the total transmission
C through the top ozone layer is less than trmin, then no
C delta-Eddington computation for the underlying column is done:
C
      do 200 k=1,1
C
C Initialize current layer properties to zero;only if total transmission
C to the top interface of the current layer exceeds the minimum, will
C these values be computed below:
C
         do nn=1,nloop
            do i=is(nn),ie(nn)
C
               rdir(i,k)   =  0.0
               rdif(i,k)   =  0.0
               tdir(i,k)   =  0.0
               tdif(i,k)   =  0.0
               explay(i,k) =  0.0
C
C Calculates the solar beam transmission, total transmission, and
C reflectivity for diffuse radiation from below at the top of the
C current layer:
C
               exptdn(i,k) = exptdn(i,k-1)*explay(i,k-1)
               rdenom      = 1./(1. - rdif(i,k-1)*rdndif(i,k-1))
               rdirexp     = rdir(i,k-1)*exptdn(i,k-1)
               tdnmexp     = tottrn(i,k-1) - exptdn(i,k-1)
               tottrn(i,k) = exptdn(i,k-1)*tdir(i,k-1) + tdif(i,k-1)*
     $                      (tdnmexp + rdndif(i,k-1)*rdirexp)*rdenom
               rdndif(i,k) = rdif(i,k-1)  +
     $                (rdndif(i,k-1)*tdif(i,k-1))*(tdif(i,k-1)*rdenom)
C
            end do
         end do
C
C Compute next layer delta-Eddington solution only if total transmission
C of radiation to the interface just above the layer exceeds trmin.
C
         call whenfgt(plon,tottrn(1,k),1,trmin,index,nval)
         if(nval.gt.0) then
CCDIR$ IVDEP
            do 100 ii=1,nval
               i=index(ii)
C
C Remember, no ozone absorption in this layer:
C
               tauray(i) = trayoslp*pflx(i,plevp)
               taugab(i) = abh2o*uth2o(i) +
     $                     abco2*utco2(i) + abo2*uto2(i)
C
               tautot    = tauray(i) + taugab(i) +  tauxar_mix_css(i)
C
               wtot      = (wray*tauray(i) + tauasc_mix_css(i) )
     $                             /tautot
C
               gtot      = (gray*wray*tauray(i)  +
     $                     gtota_mix_css(i) )
     $                             / (wtot*tautot)
C
               ftot      = (fray*wray*tauray(i)  +
     $                     ftota_mix_css (i) )
     $                             / (wtot*tautot)
C
               ts        = taus(wtot,ftot,tautot)
               ws        = omgs(wtot,ftot)
               gs        = asys(gtot,ftot)
               lm        = el(ws,gs)
               alp       = alpha(ws,coszrs(i),gs,lm)
               gam       = gamma(ws,coszrs(i),gs,lm)
               ue        = u(ws,gs,lm)
C
C Limit argument of exponential to 25, in case lm very large:
C
               arg       = dmin1(lm*ts,25.d0)
               extins    = dexp(-arg)
               ne        = n(ue,extins)
C
               rdif(i,k) = (ue+1.)*(ue-1.)*(1./extins - extins)/ne
               tdif(i,k) =   4.*ue/ne
C
C Limit argument of exponential to 25, in case coszrs is very small:
C
               arg       = dmin1(ts/coszrs(i),25.d0)
               explay(i,k) = dexp(-arg)
C
               apg       = alp + gam
               amg       = alp - gam
               rdir(i,k) = amg*(tdif(i,k)*explay(i,k) - 1.) +
     $                     apg*rdif(i,k)
               tdir(i,k) = apg*tdif(i,k) +
     $                     (amg*rdif(i,k) - (apg-1.))*explay(i,k)
C
C Under rare conditions, reflectivies and transmissivities can be
C negative; zero out any negative values
C
               rdir(i,k) = dmax1(rdir(i,k),0.d0)
               tdir(i,k) = dmax1(tdir(i,k),0.d0)
               rdif(i,k) = dmax1(rdif(i,k),0.d0)
               tdif(i,k) = dmax1(tdif(i,k),0.d0)
C
  100       continue
         end if
C
  200 continue
C
C Compute total direct beam transmission, total transmission, and
C reflectivity for diffuse radiation (from below) for both layers
C above the surface:
C
      k = 2
      do nn=1,nloop
         do i=is(nn),ie(nn)
            exptdn(i,k) = exptdn(i,k-1)*explay(i,k-1)
            rdenom      = 1./(1. - rdif(i,k-1)*rdndif(i,k-1))
            rdirexp     = rdir(i,k-1)*exptdn(i,k-1)
            tdnmexp     = tottrn(i,k-1) - exptdn(i,k-1)
            tottrn(i,k) = exptdn(i,k-1)*tdir(i,k-1) + tdif(i,k-1)*
     $                   (tdnmexp + rdndif(i,k-1)*rdirexp)*rdenom
            rdndif(i,k) = rdif(i,k-1)  +
     $                  (rdndif(i,k-1)*tdif(i,k-1))*(tdif(i,k-1)*rdenom)
         end do
      end do
C
      return
      end
