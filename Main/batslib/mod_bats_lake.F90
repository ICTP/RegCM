!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    ICTP RegCM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with ICTP RegCM.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

!
!     LAKE MODEL
!
module mod_bats_lake
!
  use mod_realkinds
  use mod_dynparam
  use mod_bats_common
  use mod_bats_internal
  use mod_bats_mppio
!
  private
!
  public :: allocate_mod_bats_lake , lakesav_i, lakesav_o
  public :: initlake , lakescatter , lakegather , lakedrv
  public :: dhlake1
!
  real(dp) , pointer , dimension(:,:,:) :: dhlake1
  integer , pointer , dimension(:,:,:) :: idep2d
  real(dp) , pointer , dimension(:,:,:) :: eta2d
  real(dp) , pointer , dimension(:,:,:) :: hi2d
  real(dp) , pointer , dimension(:,:,:) :: aveice2d
  real(dp) , pointer , dimension(:,:,:) :: hsnow2d
  real(dp) , pointer , dimension(:,:,:,:) :: tlak3d
!
  real(dp) , dimension(ndpmax) :: de , dnsty , tt
!
  ! surface thickness
  real(dp) , parameter :: surf = d_one
  ! vertical grid spacing in m
  real(dp) , parameter :: dz = surf
  ! minimum ice depth in mm: less that this is removed
  real(dp) , parameter :: iceminh = d_10
  ! reference hgt in mm for latent heat removal from ice
  real(dp) , parameter :: href = d_two * iceminh
  ! steepness factor of latent heat removal
  real(dp) , parameter :: steepf = 1.0D0  ! Tuning needed !
!
  contains
!
!-----------------------------------------------------------------------
!
  subroutine allocate_mod_bats_lake(ilake)
    implicit none
    integer , intent(in) :: ilake
    if ( ilake == 1 ) then
      call getmem3d(dhlake1,1,nnsg,1,iy,1,jxp,'mod_lake:dhlake1')
      call getmem3d(idep2d,1,nnsg,1,iy,1,jxp,'mod_lake:idep2d')
      call getmem3d(eta2d,1,nnsg,1,iy,1,jxp,'mod_lake:eta2d')
      call getmem3d(hi2d,1,nnsg,1,iy,1,jxp,'mod_lake:hi2d')
      call getmem3d(aveice2d,1,nnsg,1,iy,1,jxp,'mod_lake:aveice2d')
      call getmem3d(hsnow2d,1,nnsg,1,iy,1,jxp,'mod_lake:hsnow2d')
      call getmem4d(tlak3d,1,ndpmax,1,nnsg,1,iy,1,jxp,'mod_lake:tlak3d')
    end if
  end subroutine allocate_mod_bats_lake

  subroutine initlake(jstart,jend,istart,iend)
#ifndef IBM
  use mpi
#endif
  implicit none
  integer , intent(in) :: jstart , jend , istart , iend
#ifdef IBM
  include 'mpif.h'
#endif
! 
  integer :: i, j, n
  integer :: ierr
!
  hi2d     = 0.01D0
  aveice2d = d_zero
  hsnow2d  = d_zero
  eta2d    = d_half
  tlak3d   = 6.0D0
  idep2d   = 0

  do i = istart , iend
    do j = jstart , jend
      do n = 1 , nnsg

!     ******  initialize hostetler lake model
        if ( (lndcat1(n,j,i) > 13.9D0 .and.   &
              lndcat1(n,j,i) < 14.1D0) .and.  &
             dhlake1(n,i,j) > d_one) then
          idep2d(n,i,j) = idint(dmax1(d_two,dmin1(dhlake1(n,i,j), &
                                dble(ndpmax)))/dz)
          if ( ocld2d(n,j,i) == 2 ) then
            tlak3d(1,n,i,j) = 1.78D0
            tlak3d(2,n,i,j) = 1.78D0
            aveice2d(n,i,j) = d_1000
            hi2d(n,i,j) = d_one
            hsnow2d(n,i,j) = d_zero
          end if
          ! Azar Zarrin: Fixed unrealistic high ice tickness and
          ! high water temperatures during warm months.
          if (idep2d(n,i,j) < 50) then
            eta2d(n,i,j) = 0.5D0
          else if (idep2d(n,i,j) > 100) then
            eta2d(n,i,j) = 0.1D0
          else
            eta2d(n,i,j) = 0.3D0
          end if
        else
          idep2d(n,i,j) = 0
        end if
        if (idep2d(n,i,j) == 0) then
          hi2d(n,i,j)     = dmissval
          aveice2d(n,i,j) = dmissval
          hsnow2d(n,i,j)  = dmissval
          eta2d(n,i,j)    = dmissval
          tlak3d(:,n,i,j) = dmissval
        else if (idep2d(n,i,j) < ndpmax) then
          tlak3d(idep2d(n,i,j)+1:,n,i,j) = dmissval
        end if
      end do
    end do
  end do

  call mpi_gather(idep2d,   nnsg*iym1*jxp,mpi_integer, &
                  idep2d_io,nnsg*iym1*jxp,mpi_integer, &
                  0, mycomm,ierr)
  end subroutine initlake
!
  subroutine lakedrv(jstart,jend,istart,iend)
  implicit none
  integer , intent(in) :: jstart , jend , istart , iend
!
  real(dp) :: flwx , fswx , hsen , prec , &
             ql , tgl , tl , vl , zl , xl , evp , toth
  integer :: i , j , n
!
  do i = istart , iend
    do j = jstart , jend
      do n = 1 , nnsg
        if ( idep2d(n,i,j) > 1 ) then
          tl = sts(n,j,i)
          vl = dsqrt(usw(j,i)**d_two+vsw(j,i)**d_two)
          zl = zh(n,j,i)
          ql = qs(n,j,i)
          fswx = fsw(j,i)
          flwx = -d_one*flw(j,i)
          prec = prcp(n,j,i)*dtbat
          hsen = -d_one*sent(n,j,i)
          evp = evpr(n,j,i)
          if (nnsg == 1) then
            xl = xlat(i,j)
          else
            xl = xlat1(n,j,i)
          end if

          call lake( dtlake,tl,vl,zl,ql,fswx,flwx,hsen,xl,    &
                     tgl,prec,idep2d(n,i,j),eta2d(n,i,j),  &
                     hi2d(n,i,j),aveice2d(n,i,j),          &
                     hsnow2d(n,i,j),evp,tlak3d(:,n,i,j) )

          ! Feed back ground temperature
          tgrd(n,j,i) = tgl
          tgbrd(n,j,i) = tgl

          if ( aveice2d(n,i,j) <= iceminh ) then
            ocld2d(n,j,i) = 0 
            ldimsk(n,j,i) = 0
            lveg(n,j,i) = 14
            sfice(n,j,i) = d_zero
            sncv(n,j,i) = d_zero
            snag(n,j,i) = d_zero
          else
            ocld2d(n,j,i) = 2 
            ldimsk(n,j,i) = 2
            lveg(n,j,i) = 12
            sfice(n,j,i) = aveice2d(n,i,j)  !  units of ice = mm
            sncv(n,j,i)  = hsnow2d(n,i,j)   !  units of snw = mm
            evpr(n,j,i) = evp                 !  units of evp = mm/sec
            ! Reduce sensible heat flux for ice presence
            toth = sfice(n,j,i) + sncv(n,j,i)
            if ( toth > href ) then
              sent(n,j,i) = sent(n,j,i) * (href/toth)**steepf
            end if
          end if
        end if
      end do
    end do
  end do
 
  end subroutine lakedrv
!
!-----------------------------------------------------------------------
!
  subroutine lake(dtlake,tl,vl,zl,ql,fsw,flw,hsen,xl,tgl,  &
                  prec,ndpt,eta,hi,aveice,hsnow,evl,tprof)
 
  implicit none
!
  real(dp) :: dtlake , evl , aveice , hsen , hsnow , flw , &
             prec , ql , fsw , tl , tgl , vl , zl , eta , hi , xl
  real(dp) , dimension(ndpmax) :: tprof
  integer :: ndpt
  intent (in) hsen , ql , tl , vl , zl
  intent (in) ndpt , eta
  intent (out) tgl
  intent (inout) evl , aveice , hsnow
  intent (inout) tprof
!
  real(dp) :: ai , ea , ev , hs , ld , lu , qe , qh , tac , tk , u2
!
!***  dtlake:  time step in seconds
!***  zo:      surface roughness length
!
  real(dp) , parameter :: zo = 0.001D0
  real(dp) , parameter :: z2 = d_two
  real(dp) , parameter :: tcutoff = -0.001D0
  real(dp) , parameter :: twatui = 1.78D0
  logical , parameter :: lfreeze = .false.
  integer , parameter :: kmin = 1
!
!     interpolate winds at z1 m to 2m via log wind profile
  u2 = vl*dlog(z2/zo)/dlog(zl/zo)
  if ( u2 < d_half ) u2 = d_half
 
!     ****** Check if conditions not exist for lake ice
  if ( (aveice < 1.0D-8) .and. (tprof(1) > tcutoff) ) then
 
    ! Graziano: removed hlat. It is calculated from evaporation
    qe = -d_one*evl*wlhv
    qh = hsen

!       ******    Calculate eddy diffusivities
    call eddy(ndpt,dtlake,u2,xl,tprof)
 
!       ******    Lake temperature calc using sensible and latent heats
    call temp(ndpt,dtlake,fsw,flw,qe,qh,eta,tprof)
 
!       ******    Convective mixer
    call mixer(kmin,ndpt,tprof)

    hi     = 0.01D0
    aveice = d_zero
    hsnow  = d_zero

!     ****** Lake ice
  else
 
!       convert mixing ratio to air vapor pressure
    ea  = ql*88.0D0/(ep2+0.378D0*ql)
    tac = tl - tzero
    tk  = tzero + tprof(1)
    lu  = -emsw*sigm*tk**d_four
    ld  = flw - lu
    ev  = evl*secph         ! convert to mm/hr
    ai  = aveice / d_1000   ! convert to m
    hs  = hsnow / d_100     ! convert to m

    call ice(dtlake,fsw,ld,tac,u2,ea,hs,hi,ai,ev,prec,tprof)
    if ( .not. lfreeze ) tprof(1) = twatui

    evl    = ev/secph       ! convert evl  from mm/hr to mm/sec
    aveice = ai*d_1000      ! convert ice  from m to mm
    hsnow  = hs*d_100       ! convert snow from m depth to mm h20
    if (aveice < dlowval) aveice = d_zero
    if (hsnow < dlowval) hsnow = d_zero
 
  end if
 
  tgl = tprof(1) + tzero
 
  end subroutine lake
!
!-----------------------------------------------------------------------
!
  subroutine eddy(ndpt,dtlake,u2,xl,tprof)
 
! Computes density and eddy diffusivity
 
  implicit none
!
  integer , intent (in) :: ndpt
  real(dp) , intent (in) :: dtlake , u2 , xl
  real(dp) , dimension(ndpmax) , intent (in) :: tprof
!
  real(dp) :: demax , demin , dpdz , ks , n2 , po
  real(dp) :: zmax , rad , ri , ws , z
  integer :: k
!
!     demin molecular diffusion of heat in water
  demin = hdmw
!
!     Added to keep numerical stability of code
  demax = .50D0*dz**d_two/dtlake
  demax = .99D0*demax
!
  do k = 1 , ndpt
    dnsty(k) = d_1000*(d_one-1.9549D-05 * &
                  (dabs((tprof(k)+tzero)-277.0D0))**1.68D0)
  end do
! 
! Compute eddy diffusion profile
!
! Reference:
!
! B. Henderson-Sellers
!  New formulation of eddy diffusion thermocline models.
!  Appl. Math. Modelling, 1985, Vol. 9 December, pp. 441-446
!
 
!     Decay constant of shear velocity - Ekman profile parameter
  ks = 6.6D0*dsqrt(dsin(xl*degrad))*u2**(-1.84D0)

!     Ekman layer depth where eddy diffusion happens
  zmax = dble(ceiling(surf+40.0D0/(vonkar*ks)))

!     Surface shear velocity
  ws = 0.0012D0*u2

!     Inverse of turbulent Prandtl number
  po = d_one
 
  do k = 1 , ndpt - 1

!       Actual depth from surface
    z = surf + dble(k-1)*dz
    if (z >= zmax) then
      de(k) = demin
      cycle
    end if

    if ( k == 1 ) then
      dpdz = (dnsty(k+1)-dnsty(k))/surf
    else
      dpdz = (dnsty(k+1)-dnsty(k))/dz
    end if

!       Brunt Vaisala frequency squared : we do not mind stability,
!       we just look for energy here.
!        n2 = dabs((dpdz/dnsty(k))*egrav)
    n2 = (dpdz/dnsty(k))*egrav
    if (dabs(n2) < dlowval) then
      de(k) = demin
      cycle
    end if

!       Richardson number estimate
    rad = d_one+40.0D0*n2*((vonkar*z)/(ws*dexp(-ks*z)))**d_two
    if (rad < d_zero) rad = d_zero
    ri = (-d_one+dsqrt(rad))/20.0D0

!       Total diffusion coefficient for heat: molecular + eddy (Eqn 42)
    de(k) = demin + vonkar*ws*z*po*dexp(-ks*z) / &
                    (d_one+37.0D0*ri**d_two)
    if ( de(k) < demin ) de(k) = demin
    if ( de(k) > demax ) de(k) = demax

  end do
  de(ndpt) = demin
 
  end subroutine eddy
!
!-----------------------------------------------------------------------
!
  subroutine temp(ndpt,dtlake,fsw,flw,qe,qh,eta,tprof)
!
!*****************BEGIN SUBROUTINE TEMP********************
!             COMPUTES TEMPERATURE PROFILE                *
!**********************************************************
!
  implicit none
!
  integer , intent(in) :: ndpt
  real(dp) , intent(in) :: dtlake , eta , flw , qe , qh , fsw
  real(dp) , dimension(ndpmax) , intent(inout) :: tprof
!
  real(dp) :: bot , dt1 , dt2 , top
  integer :: k
 
!******    solve differential equations of heat transfer

  tt(1:ndpt) = tprof(1:ndpt)
 
  dt1 = (fsw*(d_one-dexp(-eta*surf))+(flw+qe+qh)) / &
          (surf*dnsty(1)*cpw)
  dt2 = -de(1)*(tprof(1)-tprof(2))/surf
  tt(1) = tt(1) + (dt1+dt2)*dtlake
 
  do k = 2 , ndpt - 1
    top = (surf+(k-2)*dz)
    bot = (surf+(k-1)*dz)
    dt1 = fsw*(dexp(-eta*top)-dexp(-eta*bot))/(dz*dnsty(k)*cpw)
    dt2 = (de(k-1)*(tprof(k-1)-tprof(k))    -    &
           de(k)  *(tprof(k)  -tprof(k+1))) / dz
    tt(k) = tt(k) + (dt1+dt2)*dtlake
  end do
 
  top = (surf+(ndpt-2)*dz)
  dt1 = fsw*dexp(-eta*top)/(dz*dnsty(ndpt)*cpw)
  dt2 = de(ndpt-1)*(tprof(ndpt-1)-tprof(ndpt))/dz
  tt(ndpt) = tt(ndpt) + (dt1+dt2)*dtlake
 
  do k = 1 , ndpt
    tprof(k) = tt(k)
    dnsty(k) = d_1000*(d_one-1.9549D-05 * &
               (dabs((tprof(k)+tzero)-277.0D0))**1.68D0)
  end do

  end subroutine temp
!
!-----------------------------------------------------------------------
!
  subroutine mixer(kmin,ndpt,tprof)
!
! Simulates convective mixing
!
  implicit none
!
  integer , intent(in) :: ndpt , kmin
  real(dp) , intent(inout) , dimension(ndpmax) :: tprof
!
  real(dp) :: avet , avev , tav , vol
  integer :: k , k2
! 
  tt(kmin:ndpt) = tprof(kmin:ndpt)
 
  do k = kmin , ndpt - 1
    avet = d_zero
    avev = d_zero
 
    if ( dnsty(k) > dnsty(k+1) ) then
 
      do k2 = kmin , k + 1
        if ( k2 == 1 ) then
          vol = surf
        else
          vol = dz
        end if
        avet = avet + tt(k2)*vol
        avev = avev + vol
      end do
 
      tav = avet/avev

      do k2 = kmin , k + 1
        tt(k2) = tav
        dnsty(k2) = d_1000*(d_one-1.9549D-05 * &
                    (dabs((tav+tzero)-277.0D0))**1.68D0)
      end do
    end if
 
  end do ! K loop
 
  tprof(kmin:ndpt) = tt(kmin:ndpt)
 
  end subroutine mixer
!
!-----------------------------------------------------------------------
!
  subroutine ice(dtx,fsw,ld,tac,u2,ea,hs,hi,aveice,evl,prec,tprof)

  implicit none
  real(dp) :: ea , evl , hi , aveice , hs , fsw , &
             ld , prec , tac , u2 , dtx
  real(dp) , dimension(ndpmax) :: tprof
  intent (in) dtx , ea , ld , prec , tac , u2
  intent (out) evl
  intent (inout) hi , aveice , hs , fsw , tprof
!
  real(dp) :: di , ds , f0 , f1 , khat , psi , q0 , qpen , t0 , t1 , &
             t2 , tf , theta , rho , xlexpc
  real(dp) :: xea , xeb , xec
  integer :: nits
!
  real(dp) , parameter :: isurf = 0.6D0
  ! attenuation coeff for ice in visible band (m-1)
  real(dp) , parameter :: lami1 = 1.5D0
  ! attenuation coeff for ice in infrared band (m-1)
  real(dp) , parameter :: lami2 = 20.0D0
  ! attenuation coeff for snow in visible band (m-1)
  real(dp) , parameter :: lams1 = 6.0D0
  ! attenuation coeff for snow in infrared band (m-1)
  real(dp) , parameter :: lams2 = 20.0D0
  ! thermal conductivity of ice (W/m/C)
  real(dp) , parameter :: ki = 2.3D0
  ! thermal conductivity of snow (W/m/C)
  real(dp) , parameter :: ks = 0.31D0
  ! standard atmospheric pressure (hPa) ????
  real(dp) , parameter :: atm = 950.0D0
  ! heat flux from water to ice (w/m2) ???
  real(dp) , parameter :: qw = 1.389D0
  ! latent heat of fusion (J/kg)
  real(dp) , parameter :: li = 334.0D03
  ! drag coefficient for the turbulent momentum flux.
  real(dp) , parameter :: cd = 0.001D0
  ! Maximum exponent
  real(dp) , parameter :: minexp = -25.0D0
!
!
!****************************SUBROUINE ICE*****************************
!     SIMULATES LAKE ICE                           
!**********************************************************************
 
  if ( (tac <= d_zero) .and. (aveice > d_zero) ) &
    hs = hs + prec*d_r100  ! convert prec(mm) to depth(m)
  if ( hs < dlowval ) hs = d_zero
 
  ! temperature of ice/snow surface
  t0 = tprof(1)
  ! freezing temp of water
  tf = d_zero
  ! approximate density of air (1 kg/m3)
  rho = rhoh2o*d_r1000
 
  khat = (ki*hs+ks*hi)/(ki*ks)
  theta = cpd*rho*cd*u2
  psi = wlhv*rho*cd*u2*ep2/atm
  evl = d_100*psi*(eomb(t0)-ea)/(wlhv*rho)
  ! amount of radiation that penetrates through the ice (W/m2)
  xea = -lams1*hs
  xeb = -lami1*hi
  xec = -lami2*hi
  if ( xea > minexp ) then
    xea = dexp(xea)
  else
    xea = d_zero
  end if
  if ( xeb > minexp ) then
    xeb = dexp(xeb)
  else
    xeb = d_zero
  end if
  if ( xec > minexp ) then
    xec = dexp(xec)
  else
    xec = d_zero
  end if

  qpen = fsw*0.7D0*((d_one-xea)/(ks*lams1) +            &
                    (xea*(d_one-xeb)/(ki*lami1))) +     &
         fsw*0.3D0*((d_one-dexp(-lams2))/(ks*lams2)+    &
                    (-lams2*hs)*(d_one-xec)/(ki*lami2))
  ! radiation absorbed at the ice surface
  fsw = fsw - qpen
 
  ! test qpen sensitivity
  !qpen = qpen * 0.5

  nits = 0
  t1 = -50.0D0
  f0 = f(t0)
  f1 = f(t1)
  do
    nits = nits + 1
    t2 = t1 - (t1-t0)*f1/(f1-f0)
    if ( dabs((t2-t1)/t1) >= 0.001D0 ) then
      t0 = t1
      t1 = t2
      f0 = f1
      f1 = f(t1)
      cycle
    end if
 
    t0 = t2
    if ( t0 >= tf ) then
 
      if ( hs > d_zero ) then
        ds = dtx*                                            &
             ((-ld+0.97D0*sigm*t4(tf)+psi*(eomb(tf)-ea)+     &
              theta*(tf-tac)-fsw)-d_one/khat*(tf-t0+qpen)) / &
             (rhosnow*li)
        if ( ds > d_zero ) ds = d_zero
        hs = hs + ds
        if ( hs < d_zero ) then
          hs = d_zero
          tprof(1) = (aveice*t0+(isurf-aveice)*tprof(2))/isurf
        end if
      end if
      if ( (dabs(hs) < dlowval) .and. (aveice > d_zero) ) then
        di = dtx*                                        &
            ((-ld+0.97D0*sigm*t4(tf)+psi*(eomb(tf)-ea) + &
             theta*(tf-tac)-fsw)-d_one/khat*(tf-t0+qpen))/ &
             (rhoice*li)
        if ( di > d_zero ) di = d_zero
        hi = hi + di
      end if
 
    else if ( t0 < tf ) then
 
      q0 = -ld + 0.97D0*sigm*t4(t0) + psi*(eomb(t0)-ea) + &
           theta*(t0-tac) - fsw
      xlexpc = -(lams1*hs+lami1*hi)
      ! Graziano : limit exponential
      if (xlexpc > minexp ) then
        qpen = fsw*0.7D0*(d_one-dexp(-(lams1*hs+lami1*hi))) + &
               fsw*0.3D0*(d_one-dexp(-(lams2*hs+lami2*hi)))
      else
        qpen = fsw
      end if
      di = dtx*(q0-qw-qpen)/(rhoice*li)
 
      hi = hi + di
    end if
 
    if ( hi <= 0.01D0 ) then
      hi = 0.01D0
      aveice = d_zero
      hs = d_zero
      tprof(1) = (hi*t0+(isurf-hi)*tprof(2))/isurf
    else
      aveice = hi
      tprof(1) = t0
    end if
    exit
  end do
 
  contains

  function t4(x)
    implicit none
    real(dp) :: t4
    real(dp) , intent(in) :: x
    t4 = (x+tzero)**d_four
  end function t4
  ! Computes air vapor pressure as a function of temp (in K)
  function tr1(x)
    implicit none
    real(dp) :: tr1
    real(dp) , intent(in) :: x
    tr1 = d_one - (tboil/(x+tzero))
  end function tr1
  function eomb(x)
    implicit none
    real(dp) :: eomb
    real(dp) , intent(in) :: x
    eomb = stdpmb*dexp(13.3185D0*tr1(x)-1.976D0*tr1(x)**d_two   &
           -0.6445D0*tr1(x)**d_three- 0.1299D0*tr1(x)**d_four)
   end function eomb
  function f(x)
    implicit none
    real(dp) :: f
    real(dp) , intent(in) :: x
    f = (-ld+0.97D0*sigm*t4(x)+psi*(eomb(x)-ea)+theta*(x-tac)-fsw)  &
        - d_one/khat*(qpen+tf-x)
  end function f
 
  end subroutine ice
!
!-----------------------------------------------------------------------
!
  subroutine lakegather

#ifndef IBM
  use mpi
#endif
  implicit none
#ifdef IBM
  include 'mpif.h'
#endif
!
  integer :: ierr
!
  call mpi_gather(eta2d,   nnsg*iym1*jxp,mpi_real8, &
                  eta2d_io,nnsg*iym1*jxp,mpi_real8, &
                  0, mycomm,ierr)
  call mpi_gather(hi2d,   nnsg*iym1*jxp,mpi_real8, &
                  hi2d_io,nnsg*iym1*jxp,mpi_real8, &
                  0, mycomm,ierr)
  call mpi_gather(aveice2d,   nnsg*iym1*jxp,mpi_real8, &
                  aveice2d_io,nnsg*iym1*jxp,mpi_real8, &
                  0, mycomm,ierr)
  call mpi_gather(hsnow2d,   nnsg*iym1*jxp,mpi_real8, &
                  hsnow2d_io,nnsg*iym1*jxp,mpi_real8, &
                  0, mycomm,ierr)
  call mpi_gather(tlak3d,   ndpmax*nnsg*iym1*jxp,mpi_real8, &
                  tlak3d_io,ndpmax*nnsg*iym1*jxp,mpi_real8, &
                  0, mycomm,ierr)

  end subroutine lakegather
!
!-----------------------------------------------------------------------
!
  subroutine lakescatter

#ifndef IBM
  use mpi
#endif
  implicit none
#ifdef IBM
  include 'mpif.h'
#endif
!
  integer :: ierr
!
  call mpi_scatter(idep2d_io,nnsg*iym1*jxp,mpi_integer, &
                   idep2d,   nnsg*iym1*jxp,mpi_integer, &
                   0, mycomm,ierr)
  call mpi_scatter(eta2d_io,nnsg*iym1*jxp,mpi_real8, &
                   eta2d,   nnsg*iym1*jxp,mpi_real8, &
                   0, mycomm,ierr)
  call mpi_scatter(hi2d_io,nnsg*iym1*jxp,mpi_real8, &
                   hi2d,   nnsg*iym1*jxp,mpi_real8, &
                   0, mycomm,ierr)
  call mpi_scatter(aveice2d_io,nnsg*iym1*jxp,mpi_real8, &
                   aveice2d,   nnsg*iym1*jxp,mpi_real8, &
                   0, mycomm,ierr)
  call mpi_scatter(hsnow2d_io,nnsg*iym1*jxp,mpi_real8, &
                   hsnow2d,   nnsg*iym1*jxp,mpi_real8, &
                   0, mycomm,ierr)
  call mpi_scatter(tlak3d_io,ndpmax*nnsg*iym1*jxp,mpi_real8, &
                   tlak3d,   ndpmax*nnsg*iym1*jxp,mpi_real8, &
                   0, mycomm,ierr)

  end subroutine lakescatter
!
!-----------------------------------------------------------------------
!
  subroutine lakesav_o(iutl)

  implicit none
  integer :: iutl
  intent (in) iutl
!
  integer :: i , j , k , n
!
#ifdef BAND
  write (iutl) (((idep2d_io(n,i,j),n=1,nnsg),i=2,iym1),j=1,jx)
#else
  write (iutl) (((idep2d_io(n,i,j),n=1,nnsg),i=2,iym1),j=2,jxm1)
#endif
#ifdef BAND
  do j = 1 , jx
#else
  do j = 2 , jxm1
#endif
    do i = 2 , iym1
      do n = 1 , nnsg
        if ( idep2d_io(n,i,j) > 1 ) then
          write(iutl) eta2d_io(n,i,j), hi2d_io(n,i,j), &
               aveice2d_io(n,i,j), hsnow2d_io(n,i,j),  &
               (tlak3d_io(k,n,i,j),k=1,idep2d_io(n,i,j))  
        end if
      end do
    end do
  end do

  end subroutine lakesav_o
!
!-----------------------------------------------------------------------
!
  subroutine lakesav_i(iutl)

  implicit none
  integer :: iutl
  intent (in) iutl
!
  integer :: i , j , k , n
!
  idep2d_io   = 0
  hi2d_io     = 0.01D0
  aveice2d_io = d_zero
  hsnow2d_io  = d_zero
  eta2d_io    = d_half
  tlak3d_io   = 6.0D0
!
#ifdef BAND
  read (iutl) (((idep2d_io(n,i,j),n=1,nnsg),i=2,iym1),j=1,jx)
#else
  read (iutl) (((idep2d_io(n,i,j),n=1,nnsg),i=2,iym1),j=2,jxm1)
#endif
#ifdef BAND
  do j = 1 , jx
#else
  do j = 2 , jxm1
#endif
    do i = 2 , iym1
      do n = 1 , nnsg
        if ( idep2d_io(n,i,j) > 1 ) then
          read(iutl) eta2d_io(n,i,j), hi2d_io(n,i,j), &
               aveice2d_io(n,i,j), hsnow2d_io(n,i,j), &
               (tlak3d_io(k,n,i,j),k=1,idep2d_io(n,i,j))  
        end if
      end do
    end do
  end do

#ifdef BAND
  do j = 1 , jx
#else
  do j = 2 , jxm1 
#endif
    do i = 2 , iym1
      do n = 1 , nnsg
        if (idep2d_io(n,i,j) == 0) then
          hi2d_io(n,i,j)     = dmissval
          aveice2d_io(n,i,j) = dmissval
          hsnow2d_io(n,i,j)  = dmissval
          eta2d_io(n,i,j)    = dmissval
          tlak3d_io(:,n,i,j) = dmissval
        else if (idep2d_io(n,i,j) < ndpmax) then
          tlak3d_io(idep2d_io(n,i,j)+1:,n,i,j) = dmissval
        end if
      end do
    end do
  end do

  end subroutine lakesav_i
!
end module mod_bats_lake
