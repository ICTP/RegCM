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

module mod_ocn_lake

  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_dynparam
  use mod_service
  use mod_runparams , only : rcmtimer , iocnrough
  use mod_ocn_internal

  implicit none

  private

  public :: allocate_mod_ocn_lake , initlake , lakedrv
  public :: lake_fillvar

  integer(ik4) :: nlakep = 0

  real(rkx) , dimension(ndpmax) :: de , dnsty , tt

  ! surface thickness
  real(rkx) , parameter :: surf = d_one
  ! vertical grid spacing in m
  real(rkx) , parameter :: dz = surf
  ! reference hgt in mm for latent heat removal from ice
  real(rkx) , parameter :: href = d_two * iceminh
  ! steepness factor of latent heat removal
  real(rkx) , parameter :: steepf = 1.0_rkx  ! Tuning needed !

  real(rkx) , pointer , dimension(:,:) :: tlak
  real(rkx) , pointer , dimension(:) :: hi , eta
  logical , pointer , public , dimension(:) :: lakmsk
  integer(ik4) , pointer , dimension(:) :: idep , ilp

  integer , public , parameter :: var_eta    = 1
  integer , public , parameter :: var_hi     = 2
  integer , public , parameter :: var_tlak   = 5

  ! Lake Malawi bottom temperature can be as high as 22.75 Celsius
  ! with surface as hot as 25.5 degrees.
  real(rkx) , parameter :: slake_trop = 26.0_rkx
  real(rkx) , parameter :: blake_trop = 22.0_rkx
  real(rkx) , parameter :: int_trop = slake_trop - blake_trop

  interface lake_fillvar
    module procedure lake_fillvar_real8_1d
    module procedure lake_fillvar_real8_2d
  end interface lake_fillvar

  contains

  subroutine allocate_mod_ocn_lake
    implicit none
    integer :: i , lp
    nlakep = 0
    call getmem1d(lakmsk,1,nocnp,'ocn::initlake::lakmsk')
    lakmsk = (ilake == 1)
    nlakep = count(lakmsk)
#ifdef DEBUG
    write(ndebug,*) 'NUMBER OF LAKE POINTS : ',nlakep
#endif
    if ( nlakep == 0 ) return
    call getmem1d(idep,1,nlakep,'ocn::initlake::idep')
    call getmem1d(ilp,1,nlakep,'ocn::initlake::ilp')
    call getmem1d(hi,1,nlakep,'ocn::initlake::hi')
    call getmem1d(eta,1,nlakep,'ocn::initlake::eta')
    call getmem2d(tlak,1,nlakep,1,ndpmax,'ocn::initlake::tlak')
    lp = 1
    do i = iocnbeg , iocnend
      if ( lakmsk(i) ) then
        ilp(lp) = i
        lp = lp + 1
      end if
    end do
  end subroutine allocate_mod_ocn_lake

  subroutine initlake
    implicit none
    integer(ik4) :: i , k , lp , n
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'initlake'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    ! initialize hostetler lake model

    if ( nlakep == 0 ) then
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
      return
    end if

    if ( rcmtimer%integrating( ) ) then
      do lp = 1 , nlakep
        i = ilp(lp)
        idep(lp) = int(max(d_two,min(dhlake(i),real(ndpmax,rkx)))/dz)
      end do
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
      return
    end if

    hi(:)     = dmissval
    eta(:)    = dmissval
    tlak(:,:) = dmissval

    do lp = 1 , nlakep
      i = ilp(lp)
      idep(lp) = int(max(d_two,min(dhlake(i),real(ndpmax,rkx)))/dz)
      hi(lp) = iceminh
      ! Azar Zarrin: Fixed unrealistic high ice tickness and
      ! high water temperatures during warm months.
      ! Graziano: Take a data driven approach.
      ! The attenuation coefficient is function of the turbidity
      ! of the water. We assume here that the more deep, the less
      ! turbid the water is. Real data from:
      !
      ! http://www.waterontheweb.org/under/lakeecology/04_light.html
      !
      if ( idep(lp) < 5 ) then
        ! A High eutrophic lake with suspended sediments can
        ! reach even value of vertical extinction coefficient of 4 !!!
        ! This means euphotic zone very limited.
        eta(lp) = -1.20_rkx
      else if ( idep(lp) > 5 .and. idep(lp) < 10) then
        ! Something like Mesotrophic
        eta(lp) = -0.80_rkx
      else if ( idep(lp) >= 10 .and. idep(lp) < 40) then
        eta(lp) = -0.60_rkx
      else if ( idep(lp) >= 40 .and. idep(lp) < 100) then
        eta(lp) = -0.40_rkx
      else
        ! This is a mean value for Great Lakes.
        ! Oligotropic lake value.
        eta(lp) = -0.20_rkx
      end if
      ! Put winter surface water a bit colder and summer or tropical
      ! surface water a little warmer to ease spinup nudging in the
      ! correct direction the profile.
      if ( abs(lat(i)) > 25.0_rkx ) then
        tlak(lp,1) = max(min(tgrd(i)-tzero+d_one,20.0_rkx),4.0_rkx)
        tlak(lp,2) = tlak(lp,1)
        if ( idep(lp) >= 3 ) then
          do k = 3 , idep(lp)
            tlak(lp,k) = min(max(tlak(lp,k-1)-0.1,4.0_rkx),20.0_rkx)
          end do
        end if
      else
        ! This needs tuning for tropical lakes.
        tlak(lp,1) = slake_trop
        tlak(lp,2) = slake_trop
        if ( idep(lp) >= 3 ) then
          if ( idep(lp) <= 20 ) then
            tlak(lp,3:idep(lp) ) = slake_trop
          else
            tlak(lp,3:20) = slake_trop
            do n = 21 , min(39,idep(lp))
              tlak(lp,n) = slake_trop - real(n-20,rkx)/20.0_rkx*int_trop
            end do
            if ( idep(lp) >= 40 ) then
              tlak(lp,40:idep(lp) ) = blake_trop
            end if
          end if
        end if
        !tlak(lp,1) = min(max(tgrd(i)-tzero+d_one,18.0_rkx),25.0_rkx)
        !tlak(lp,2) = min(max(tlak(lp,1)-d_half,18.0_rkx),25.0_rkx)
        !if ( idep(lp) >= 3 ) then
        !  do k = 3 , idep(lp)
        !    tlak(lp,k) = min(max(tlak(lp,k-1)-d_half,18.0_rkx),25.0_rkx)
        !  end do
        !end if
      end if
      sfice(i) = d_zero
      sncv(i) = d_zero
    end do
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine initlake

  subroutine lakedrv
    implicit none
    real(rkx) :: flwx , fswx , hsen , prec , qs , tgl , tl , vl , zl
    real(rkx) :: xl , toth
    real(rkx) :: age , age1 , age2 , arg , arg2 , cdr , cdrmin , cdrn
    real(rkx) :: cdrx , clead , dela , dela0 , delq , dels , delt
    real(rkx) :: fact , factuv , qgrd , qgrnd , qice , rhosw , rhosw3
    real(rkx) :: ribd , ribl , ribn
    real(rkx) :: sold , vspda , u1 , tc , visa , rho
    real(rkx) :: wt1 , wt2
    real(rkx) , dimension(ndpmax) :: tp

    integer(ik4) :: lp , i
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'lakedrv'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    if ( nlakep == 0 ) then
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
      return
    end if

    wt1 = (threedays-dtlake)/threedays
    wt2 = dtlake/threedays

    do lp = 1 , nlakep
      i = ilp(lp)
      tl = tatm(i)
      tc = tl - tzero
      sold = sncv(i)
      vl = sqrt(usw(i)**2+vsw(i)**2)
      zl = ht(i)
      qs = qv(i)
      fswx = rswf(i)
      flwx = -d_one*rlwf(i)
      prec = prcp(i)*dtlake
      hsen = -d_one*sent(i)
      xl = lat(i)
      tp = tlak(lp,:)
      rho = rhox(i)
      visa = 1.326e-5_rkx*(d_one + 6.542e-3_rkx * tc +     &
                                   8.301e-6_rkx * tc*tc -  &
                                   4.840e-9_rkx * tc*tc*tc)

      call lake( dtlake,tl,vl,zl,qs,fswx,flwx,hsen,xl, &
                 tgl,prec,idep(lp),eta(lp),hi(lp),sfice(i), &
                 sncv(i),evpr(i),tp,sfps(i),rho )

      tlak(lp,:) = tp
      tgrd(i)  = tgl
      tgbrd(i) = tgl
      qgrd = pfqsat(tgrd(i),sfps(i))
      delt = tatm(i) - tgrd(i)
      ! Move to specific humidities
      qs = qs/(d_one+qs)
      qgrd = qgrd/(d_one+qgrd)
      delq = (qs - qgrd)

      if ( sfice(i) <= iceminh ) then
        mask(i) = 3
        tgb(i)  = tgl
        sfice(i) = d_zero
        sncv(i) = d_zero
        snag(i) = d_zero
        sm(i) = d_zero
        ribd = usw(i)**2 + vsw(i)**2 + wtur**2
        vspda = sqrt(ribd)
        cdrn = (vonkar/log(ht(i)/zoce))**2
        ribn = ht(i)*egrav*(delt/tatm(i))
        br(i) = ribn/ribd
        if ( br(i) < d_zero ) then
          cdrx = cdrn*(d_one+24.5_rkx*sqrt(-cdrn*br(i)))
        else
          cdrx = cdrn/(d_one+11.5_rkx*br(i))
        end if
        cdrmin = max(0.25_rkx*cdrn,6.0e-4_rkx)
        if ( cdrx < cdrmin ) cdrx = cdrmin
        rah1(i) = d_one/cdrx
        ram1(i) = d_one/cdrx
        drag(i) = cdrx*vspda*rhox(i)
        evpr(i) = -drag(i)*delq
        sent(i) = -drag(i)*cpd*delt
      else
        mask(i) = 4
        tgbrd(i) = tzero - d_two
        toth = sfice(i) + sncv(i)
        if ( sncv(i) < dlowval ) then
          sncv(i) = d_zero
          snag(i) = d_zero
        else
          arg = 5.0e3_rkx*(d_one/tzero-d_one/tgrd(i))
          age1 = exp(arg)
          arg2 = min(d_zero,d_10*arg)
          age2 = exp(arg2)
          age = age1 + age2 + age3
          dela0 = 1.0e-6_rkx*dtlake
          dela = dela0*age
          dels = d_r10*max(d_zero,sncv(i)-sold)
          snag(i) = (snag(i)+dela)*(d_one-dels)
          if ( snag(i) < dlowval ) snag(i) = d_zero
          if ( sncv(i) > 800.0_rkx ) snag(i) = d_zero
        end if
        age = (d_one-d_one/(d_one+snag(i)))
        cdrn = (vonkar/log(ht(i)/zlnd))**2
        if ( delt < d_zero ) then
          u1 = wtur + d_two*sqrt(-delt)
        else
          u1 = wtur
        end if
        ribd = usw(i)**2 + vsw(i)**2 + u1**2
        vspda = sqrt(ribd)
        ribn = ht(i)*egrav*(delt/tatm(i))
        br(i) = ribn/ribd
        if ( br(i) < d_zero ) then
          cdr = cdrn*(d_one+24.5_rkx*sqrt(-cdrn*br(i)))
        else
          cdr = cdrn/(d_one+11.5_rkx*br(i))
        end if
        cdrmin = max(0.25_rkx*cdrn,6.0e-4_rkx)
        if ( cdr < cdrmin ) cdr = cdrmin
        rhosw = 0.10_rkx*(d_one+d_three*age)
        rhosw3 = rhosw**3
        cdrn = (vonkar/log(ht(i)/zsno))**2
        ribl = (d_one-271.5_rkx/tatm(i))*ht(i)*egrav/ribd
        if ( ribl < d_zero ) then
          clead = cdrn*(d_one+24.5_rkx*sqrt(-cdrn*ribl))
        else
          clead = cdrn/(d_one+11.5_rkx*br(i))
        end if
        cdrx = (d_one-aarea)*cdr + aarea*clead
        rah1(i) = d_one/cdrx
        ram1(i) = d_one/cdrx
        drag(i) = cdrx*vspda*rhox(i)
        qice = 3.3e-3_rkx * stdp/sfps(i)
        qgrnd = ((d_one-aarea)*cdr*qgrd + aarea*clead*qice)/cdrx
        tgb(i) = ((d_one-aarea)*cdr*tgrd(i) + aarea*clead*(tzero-1.8_rkx))/cdrx
        delt = tatm(i) - tgb(i)
        delq = qs - qgrnd
        evpr(i) = -drag(i)*delq
        sent(i) = -drag(i)*cpd*delt
        ! Reduce sensible heat flux for ice presence
        if ( toth > href ) then
          sent(i) = sent(i) * (href/toth)**steepf
        end if
      end if
      if ( abs(sent(i)) < dlowval ) sent(i) = d_zero
      if ( abs(evpr(i)) < dlowval ) evpr(i) = d_zero
      fact = log(ht(i)*d_half)/log(ht(i)/zoce)
      factuv = log(ht(i)*d_r10)/log(ht(i)/zoce)
      u10m(i) = usw(i)*(d_one-factuv)
      v10m(i) = vsw(i)*(d_one-factuv)
      rhoa(i) = rhox(i)
      um10(i) = um10(i) * wt1 + sqrt(u10m(i)**2+v10m(i)**2) * wt2
      ustr(i) = sqrt(sqrt((u10m(i)*drag(i))**2 + &
                          (v10m(i)*drag(i))**2)/rhoa(i))
      ustr(i) = max(ustr(i),1.0e-5_rkx)
      call ocnrough(zoo(i),ustr(i),um10(i),vl,visa)
      taux(i) = drag(i) * (u10m(i)/usw(i))
      tauy(i) = drag(i) * (v10m(i)/vsw(i))
      t2m(i) = tatm(i) - delt*fact
      q2m(i) = qs - delq*fact
    end do

#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif

    contains

#include <pfesat.inc>
#include <pfqsat.inc>
    !
    ! our formulation for zo
    !
    subroutine ocnrough(zo,ustar,um10,wc,visa)
      implicit none
      real(rkx) , intent (in) :: ustar , um10 , wc , visa
      real(rkx) , intent (out) :: zo
      real(rkx) :: cp , charnockog , alph
      ! if surface roughness not provided by wave model
      ! Wave age. The wind here is the mean last N days wind
      cp = 1.2_rkx*um10
      ! Smith et al. (1992), Carlsson et al. (2009)
      ! Charnock parameter as power function of the wave age
      ! We consider here dominant wind sea waves
      ! Swell dominated sea would require a wave model...
      charnockog = regrav*0.063_rkx*(cp/ustar)**(-0.4_rkx)
      if ( iocnrough == 1 ) then
        zo = 0.0065_rkx*regrav*ustar*ustar
      else if ( iocnrough == 2 ) then
        zo = 0.013_rkx*regrav*ustar*ustar + 0.11_rkx*visa/ustar
      else if ( iocnrough == 3 ) then
        zo = 0.017_rkx*regrav*ustar*ustar
      else if ( iocnrough == 4 ) then
        ! C.H. Huang, 2012
        ! Modification of the Charnock Wind Stress Formula
        ! to Include the Effects of Free Convection and Swell
        ! Advanced Methods for Practical Applications in Fluid Mechanics
        zo = charnockog*(ustar*ustar*ustar+0.11_rkx*wc*wc*wc)**twot
      else if ( iocnrough == 5 ) then
        if ( um10 < 10.0_rkx ) then
          alph = 0.011_rkx
        else if ( um10 > 10.0_rkx .and. um10 < 18.0_rkx ) then
          alph = 0.011_rkx + 0.000875*(um10-10.0_rkx)
        else if ( um10 > 18.0_rkx .and. um10 < 25.0_rkx ) then
          alph = 0.018_rkx
        else
          alph = max(2.e-3_rkx,0.018_rkx / &
               (d_one+0.050_rkx*(ustar-0.02_rkx)**2 - &
                0.018_rkx*(ustar-0.02_rkx)**1.6_rkx))
        end if
        zo = alph*regrav*ustar*ustar + 0.11_rkx*visa/ustar
      else
        zo = charnockog*ustar*ustar
      end if
      zo = max(zo,1.0e-8_rkx)
     end subroutine ocnrough

  end subroutine lakedrv

  subroutine lake(dtlake,tl,vl,zl,ql,fsw,flw,hsen,xl,tgl,  &
                  prec,ndpt,eta,hi,aveice,hsnow,evl,tprof,ps,dens)
    implicit none
    real(rkx) , intent(in) :: dtlake , hsen , flw , &
               prec , ql , fsw , tl , vl , zl , eta , xl , ps , dens
    real(rkx) , intent(out) :: tgl
    real(rkx) , intent(inout) :: hi , evl , aveice , hsnow
    real(rkx) , dimension(ndpmax) , intent(inout) :: tprof
    integer(ik4) , intent(in) :: ndpt
    real(rkx) :: ai , ea , ev , hs , ld , lu , qe , qh , tac , tk , u2
    ! zo: surface roughness length
    real(rkx) , parameter :: zo = 0.001_rkx
    real(rkx) , parameter :: z2 = d_two
    real(rkx) , parameter :: tcutoff = -0.001_rkx
    real(rkx) , parameter :: twatui = 1.78_rkx
    logical , parameter :: lfreeze = .true.
    integer(ik4) , parameter :: kmin = 1
    integer(ik4) , parameter :: kmax = 200

    ! interpolate winds at z1 m to 2m via log wind profile
    u2 = vl*log(z2/zo)/log(zl/zo)
    if ( u2 < d_half ) u2 = d_half

    ! Check if conditions not exist for lake ice
    if ( (aveice < iceminh) .and. (tprof(1) > tcutoff) ) then

      ! Graziano: removed hlat. It is calculated from evaporation
      qe = -d_one*evl*wlhv
      qh = hsen

      ! Calculate eddy diffusivities
      call lakeeddy(ndpt,dtlake,u2,xl,tprof)

      ! Lake temperature calc using sensible and latent heats
      call laketemp(ndpt,dtlake,fsw,flw,qe,qh,eta,tprof)

      ! Convective mixer
      call lakemixer(kmin,kmax,ndpt,tprof)

      hi     = iceminh
      aveice = d_zero
      hsnow = d_zero

    else

      ! Calculate eddy diffusivities
      call lakeeddy(ndpt,dtlake,u2,xl,tprof)

      ! Lake ice
      ! convert mixing ratio to air vapor pressure
      ea  = ql*88.0_rkx/(ep2+0.378_rkx*ql)
      tac = tl - tzero
      tk  = tzero + tprof(1)
      lu  = -emsw*sigm*tk**4
      ld  = flw - lu
      ev  = evl*secph          ! convert to mm/hr
      ai  = aveice
      ! RegCM snow is mm h2o eq. Show depth is this * 10.0.
      ! Then we have from mm to m, i.e. multiply by 10E-3.
      ! Effect is multiply by 10E-2
      hs  = hsnow * d_r100    ! convert to m.

      call lakeice(dtlake,fsw,ld,tac,u2,ea,hs,hi,ai,ev,prec,ps,tprof,dens)

      ! Convective mixer
      call lakemixer(kmin,kmax,ndpt,tprof)

      if ( .not. lfreeze ) tprof(1) = twatui

      evl    = ev/secph      ! convert evl  from mm/hr to mm/sec
      aveice = ai
      hsnow  = hs*d_100      ! convert back. See Above.
      if ( aveice < dlowval ) then
        aveice = d_zero
        hsnow = d_zero
      end if
    end if
    tgl = tprof(1) + tzero

  end subroutine lake

  subroutine lakeeddy(ndpt,dtlake,u2,xl,tprof)

    ! Computes density and eddy diffusivity

    implicit none

    integer(ik4) , intent (in) :: ndpt
    real(rkx) , intent (in) :: dtlake , u2 , xl
    real(rkx) , dimension(ndpmax) , intent (in) :: tprof

    real(rkx) :: demax , demin , dpdz , ks , n2 , po
    real(rkx) :: zmax , rad , ri , ws , z
    integer(ik4) :: k

    ! demin molecular diffusion of heat in water
    demin = hdmw

    ! Added to keep numerical stability of code
    demax = 0.50_rkx*dz**2/dtlake
    demax = 0.99_rkx*demax

    do k = 1 , ndpt
      dnsty(k) = d_1000*(d_one-1.9549e-5_rkx*(abs(tprof(k)-4.0_rkx))**1.68_rkx)
    end do

    ! Compute eddy diffusion profile
    !
    ! Reference:
    !
    ! B. Henderson-Sellers
    !  New formulation of eddy diffusion thermocline models.
    !  Appl. Math. Modelling, 1985, Vol. 9 December, pp. 441-446
    !
    ! Decay constant of shear velocity - Ekman profile parameter
    if ( xl > 25.0_rkx ) then
      ks = 6.6_rkx*sqrt(sin(xl*degrad))*u2**(-1.84_rkx)
    else
      ks = 0.001_rkx
    end if

    ! Ekman layer depth where eddy diffusion happens
    zmax = real(ceiling(surf+40.0_rkx/(vonkar*ks)),rkx)

    ! Surface shear velocity
    ws = 0.0012_rkx*u2

    ! Inverse of turbulent Prandtl number
    po = d_one

    do k = 1 , ndpt - 1

      ! Actual depth from surface
      z = surf + real(k-1,rkx)*dz
      if (z >= zmax) then
        de(k) = demin
        cycle
      end if

      if ( k == 1 ) then
        dpdz = (dnsty(k+1)-dnsty(k))/surf
      else
        dpdz = (dnsty(k+1)-dnsty(k))/dz
      end if

      ! Brunt Vaisala frequency squared : we do not mind stability,
      ! we just look for energy here.
      ! n2 = abs((dpdz/dnsty(k))*egrav)
      n2 = (dpdz/dnsty(k))*egrav
      if (abs(n2) < dlowval) then
        de(k) = demin
        cycle
      end if

      ! Richardson number estimate
      ! Total diffusion coefficient for heat: molecular + eddy (Eqn 42)
      if ( ks*z > 12.0_rkx ) then
        de(k) = demin
      else
        rad = max(d_zero,d_one+40.0_rkx*n2*((vonkar*z)/(ws*exp(-ks*z)))**2)
        ri = (-d_one+sqrt(rad))/20.0_rkx
        de(k) = demin + vonkar*ws*z*po*exp(-ks*z) / (d_one+37.0_rkx*ri**2)
      end if
      if ( de(k) < demin ) de(k) = demin
      if ( de(k) > demax ) de(k) = demax
    end do
    de(ndpt) = demin
  end subroutine lakeeddy

  subroutine laketemp(ndpt,dtlake,fsw,flw,qe,qh,eta,tprof)
    implicit none
    integer(ik4) , intent(in) :: ndpt
    real(rkx) , intent(in) :: dtlake , eta , flw , qe , qh , fsw
    real(rkx) , dimension(ndpmax) , intent(inout) :: tprof
    real(rkx) :: bot , dt1 , dt2 , top
    integer(ik4) :: k

    ! Computes temperature profile
    ! solve differential equations of heat transfer

    tt(1:ndpt) = tprof(1:ndpt)

    dt1 = (fsw*(d_one-exp(eta*surf))+(flw+qe+qh)) / &
            (surf*dnsty(1)*cpw)
    dt2 = -de(1)*(tprof(1)-tprof(2))/surf
    tt(1) = tt(1) + (dt1+dt2)*dtlake

    do k = 2 , ndpt - 1
      top = (surf+(k-2)*dz)
      bot = (surf+(k-1)*dz)
      dt1 = fsw*(exp(eta*top)-exp(eta*bot))/(dz*dnsty(k)*cpw)
      dt2 = (de(k-1)*(tprof(k-1)-tprof(k))    -    &
             de(k)  *(tprof(k)  -tprof(k+1))) / dz
      tt(k) = tt(k) + (dt1+dt2)*dtlake
    end do

    top = (surf+(ndpt-2)*dz)
    dt1 = fsw*exp(eta*top)/(dz*dnsty(ndpt)*cpw)
    dt2 = de(ndpt-1)*(tprof(ndpt-1)-tprof(ndpt))/dz
    tt(ndpt) = tt(ndpt) + (dt1+dt2)*dtlake

    do k = 1 , ndpt
      tprof(k) = tt(k)
      dnsty(k) = d_1000*(d_one-1.9549e-5_rkx*(abs(tprof(k)-4.0_rkx))**1.68_rkx)
    end do
  end subroutine laketemp

  subroutine lakemixer(kmin,kmax,ndpt,tprof)
    implicit none
    integer(ik4) , intent(in) :: ndpt , kmin , kmax
    real(rkx) , intent(inout) , dimension(ndpmax) :: tprof
    real(rkx) :: avet , avev , tav , vol
    integer(ik4) :: k , k2

    ! Simulates convective mixing
    tt(kmin:ndpt) = tprof(kmin:ndpt)

    do k = max(kmin,2) , min(kmax,ndpt-1)
      avet = d_zero
      avev = d_zero

      if ( dnsty(k) > dnsty(k+1) ) then

        do k2 = k - 1 , k + 1
          if ( k2 == 1 ) then
            vol = surf
          else
            vol = dz
          end if
          avet = avet + tt(k2)*vol
          avev = avev + vol
        end do

        tav = avet/avev

        do k2 = k - 1 , k + 1
          tt(k2) = tav
          dnsty(k2) = d_1000*(d_one-1.9549e-5_rkx*(abs(tav-4.0_rkx))**1.68_rkx)
        end do
      end if

    end do ! K loop

    tprof(kmin:ndpt) = tt(kmin:ndpt)
  end subroutine lakemixer

  subroutine lakeice(dtx,fsw,ld,tac,u2,ea,hs,hi,aveice,evl,prec,ps,tprof,dens)
    implicit none
    real(rkx) , intent(in) :: ea , ld , prec , tac , u2 , dtx , ps , dens , fsw
    real(rkx) , intent(out) :: evl
    real(rkx) , intent(inout) :: hi , aveice , hs
    real(rkx) , dimension(ndpmax) , intent(inout) :: tprof
    real(rkx) :: di , ds , f0 , f1 , khat , psi , q0 , qpen , t0 , t1 , &
                 t2 , tf , theta , rho
    real(rkx) , parameter :: isurf = 0.6_rkx
    ! attenuation coeff for ice in visible band (m-1)
    real(rkx) , parameter :: lami1 = 1.5_rkx
    ! attenuation coeff for ice in infrared band (m-1)
    real(rkx) , parameter :: lami2 = 20.0_rkx
    ! attenuation coeff for snow in visible band (m-1)
    real(rkx) , parameter :: lams1 = 6.0_rkx
    ! attenuation coeff for snow in infrared band (m-1)
    real(rkx) , parameter :: lams2 = 20.0_rkx
    ! thermal conductivity of ice (W/m/C)
    real(rkx) , parameter :: ki = 2.3_rkx
    ! thermal conductivity of snow (W/m/C)
    ! 0.078 W m-1 K-1 for new snow
    ! 0.290 W m-1 K-1 for an ubiquitous wind slab
    ! the average bulk value for the full snowpack is 0.14 W m-1 K-1
    ! Sturm, Perovich and Holmgren
    ! Journal of Geophysical Research: Oceans (1978-2012)
    !real(rkx) , parameter :: ks = 0.14_rkx
    real(rkx) , parameter :: ks = 0.31_rkx
    ! heat flux from water to ice (w m-2)
    real(rkx) , parameter :: qw = 1.389_rkx
    ! latent heat of fusion (J/kg)
    real(rkx) , parameter :: li = 334.0e3_rkx
    ! drag coefficient for the turbulent momentum flux.
    real(rkx) , parameter :: cd = 0.001_rkx
    integer(ik4) :: icount
    integer(ik4) , parameter :: maxiter = 10

    if ( (tac <= d_zero) .and. (aveice > d_zero) ) then
      ! National Weather Service indicates the average snowfall is in a
      ! ratio of 10 inches of snow to 1 inch of equivalent rainfall
      ! So this is convert to mm and to snowfall ( prec*10/1000 )
      hs = hs + prec*d_r100  ! convert prec(mm) to depth(m)
      ! Maximum snow depth suggested by Patterson and Hamblin (1988)
      ds = hi*(rhoh2o-rhoice)/rhosnow
      if ( hs > ds ) then
        ds = hs - ds
        hi = hi + ds
        hs = hs - ds
      end if
    end if
    if ( hs < dlowval ) hs = d_zero

    ! temperature of ice/snow surface
    t0 = tprof(1)
    ! temperature of the water under the snow
    tf = 0.0_rkx
    ! density of air (1 kg/m3)
    rho = dens

    khat = (ki*hs+ks*hi)/(ki*ks)
    theta = cpw0*rho*cd*u2
    psi = wlhs*rho*cd*u2*ep2/(ps*d_r100)
    evl = psi*(eomb(t0)-ea)/(wlhs*rho)
    ! amount of radiation that penetrates through the ice (W/m2)
    qpen = fsw * 0.7_rkx *                                      &
       ((d_one-exp(-lams1*hs))                  / (ks*lams1)  + &
        (exp(-lams1*hs))*(d_one-exp(-lami1*hi)) / (ki*lami1)) + &
           fsw * 0.3_rkx *                                      &
       ((d_one-exp(-lams2))                     / (ks*lams2)  + &
        ((-lams2*hs))*(d_one-exp(-lami2*hi))    / (ki*lami2))
    t1 = -50.0_rkx
    f0 = f(t0)
    f1 = f(t1)
    icount = 0
    do
      if ( abs(f1-f0) < 1.0e-8_rkx ) then
        t2 = t1
        exit
      end if
      t2 = t1 - (t1-t0)*f1/(f1-f0)
      if ( abs(t2-t1) < 0.001_rkx .or. t2 > 0.0_rkx .or. &
           icount == maxiter ) then
        exit
      end if
      t0 = t1
      t1 = t2
      f0 = f1
      f1 = f(t1)
      icount = icount + 1
    end do

    t0 = t2

    if ( t0 >= tf ) then
      if ( hs > d_zero ) then
        ds = dtx *                                                  &
             ( (-ld + emsw*sigm*t4(tf) + psi * (eomb(tf)-ea) +      &
                theta*(tf-tac)-fsw) - d_one/khat * (tf-t0+qpen) ) / &
              (rhosnowp*li)
        if ( ds > d_zero ) ds = d_zero
        hs = hs + ds * 10.0_rkx
        if ( hs < d_zero ) then
          hs = d_zero
        end if
      end if
      if ( (abs(hs) < dlowval) .and. (aveice > d_zero) ) then
        di = dtx *                                                &
            ( (-ld + emsw*sigm*t4(tf) + psi * (eomb(tf)-ea) +     &
              theta*(tf-tac)-fsw) - d_one/khat * (tf-t0+qpen) ) / &
             (rhoice*li)
        if ( di > d_zero ) di = d_zero
        hi = hi + di
      end if
    else
      q0 = -ld + emsw*sigm*t4(t0) + psi*(eomb(t0)-ea) + &
           theta*(t0-tac) - fsw
      qpen = fsw*0.7_rkx*(d_one-exp(-(lams1*hs+lami1*hi))) + &
             fsw*0.3_rkx*(d_one-exp(-(lams2*hs+lami2*hi)))
      di = dtx*(q0-qw-qpen)/(rhoice*li)
      hi = hi + di
    end if

    if ( hi <= iceminh ) then
      hi = iceminh
      aveice = d_zero
      hs = d_zero
      tprof(1) = (hi*t0+(isurf-hi)*tprof(2))/isurf
    else
      aveice = hi
      tprof(1) = min(t0,0.0_rkx)
    end if

    contains

    pure real(rkx) function t4(x)
      implicit none
      real(rkx) , intent(in) :: x
      t4 = (x+tzero)**4
    end function t4

    ! Computes air vapor pressure as a function of temp (in K)
    pure real(rkx) function eomb(x)
      implicit none
      real(rkx) , intent(in) :: x
      real(rkx) :: tr1
      tr1 = d_one - (tboil/(x+tzero))
      eomb = stdpmb*exp(13.3185_rkx*tr1 - 1.976_rkx*tr1**2 - &
                         0.6445_rkx*tr1**3 - 0.1299_rkx*tr1**4)
     end function eomb

    pure real(rkx) function f(x)
      implicit none
      real(rkx) , intent(in) :: x
      f = (-ld + emsw*sigm*t4(x) + psi*(eomb(x)-ea) + &
           theta*(x-tac)-fsw) - d_one/khat*(qpen+tf-x)
    end function f

  end subroutine lakeice

  subroutine lake_fillvar_real8_1d(ivar,rvar,idir)
    implicit none
    integer(ik4) , intent(in) :: ivar , idir
    real(rkx) , intent(inout) , pointer , dimension(:) :: rvar
    real(rkx) , pointer , dimension(:) :: p

    if ( nlakep == 0 ) then
      if ( idir == 0 ) rvar = dmissval
      return
    end if
    select case (ivar)
      case (var_eta)
        p => eta
      case (var_hi)
        p => hi
      case default
        return
    end select
    if ( idir == 1 ) then
      p = pack(rvar,lakmsk)
    else
      rvar = dmissval
      rvar = unpack(p,lakmsk,rvar)
    end if
  end subroutine lake_fillvar_real8_1d

  subroutine lake_fillvar_real8_2d(ivar,rvar,idir)
    implicit none
    integer(ik4) , intent(in) :: ivar , idir
    real(rkx) , intent(inout) , pointer , dimension(:,:) :: rvar
    real(rkx) , pointer , dimension(:,:) :: p
    integer(ik4) :: k

    if ( nlakep == 0 ) then
      if ( idir == 0 ) rvar = dmissval
      return
    end if
    select case (ivar)
      case (var_tlak)
        p => tlak
      case default
        return
    end select
    if ( idir == 1 ) then
      do k = 1 , ndpmax
        p(:,k) = pack(rvar(:,k),lakmsk)
      end do
    else
      rvar = dmissval
      do k = 1 , ndpmax
        rvar(:,k) = unpack(p(:,k),lakmsk,rvar(:,k))
      end do
    end if
  end subroutine lake_fillvar_real8_2d

end module mod_ocn_lake
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
