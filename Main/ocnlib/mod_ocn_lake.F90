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
module mod_ocn_lake

  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_service
  use mod_runparams , only : xmonth , ktau
  use mod_ocn_internal

  implicit none

  private

  public :: allocate_mod_ocn_lake , initlake , lakedrv
  public :: lake_fillvar

  integer(ik4) :: nlakep = 0

  real(rk8) , dimension(ndpmax) :: de , dnsty , tt

  ! surface thickness
  real(rk8) , parameter :: surf = d_one
  ! vertical grid spacing in m
  real(rk8) , parameter :: dz = surf
  ! reference hgt in mm for latent heat removal from ice
  real(rk8) , parameter :: href = d_two * iceminh
  ! steepness factor of latent heat removal
  real(rk8) , parameter :: steepf = 1.0D0  ! Tuning needed !
!
  real(rk8) , pointer , dimension(:,:) :: tlak
  real(rk8) , pointer , dimension(:) :: hi , eta
  logical , pointer , public , dimension(:) :: lakmsk
  integer(ik4) , pointer , dimension(:) :: idep , ilp

  integer , public , parameter :: var_eta    = 1
  integer , public , parameter :: var_hi     = 2
  integer , public , parameter :: var_tlak   = 5

  interface lake_fillvar
    module procedure lake_fillvar_real8_1d
    module procedure lake_fillvar_real8_2d
  end interface lake_fillvar

  contains
!
!-----------------------------------------------------------------------
!
  subroutine allocate_mod_ocn_lake
    implicit none
    integer :: i , lp
    nlakep = 0
    call getmem1d(lakmsk,1,nocnp,'ocn::initlake::lakmsk')
    lakmsk = (mask >= 3)
    nlakep = count(lakmsk)
#ifdef DEBUG
    write(ndebug+myid,*) 'NUMBER OF LAKE POINTS : ',nlakep
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
    integer(ik4) :: i , k , lp
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

    if ( ktau /= 0 ) then
      do lp = 1 , nlakep
        i = ilp(lp)
        idep(lp) = idint(dmax1(d_two,dmin1(dhlake(i),dble(ndpmax)))/dz)
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
      idep(lp) = idint(dmax1(d_two,dmin1(dhlake(i),dble(ndpmax)))/dz)
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
        eta(lp) = -1.20D0
      else if ( idep(lp) > 5 .and. idep(lp) < 10) then
        ! Something like Mesotrophic
        eta(lp) = -0.80D0
      else if ( idep(lp) >= 10 .and. idep(lp) < 40) then
        eta(lp) = -0.60D0
      else if ( idep(lp) >= 40 .and. idep(lp) < 100) then
        eta(lp) = -0.40D0
      else
        ! This is a mean value for Great Lakes.
        ! Oligotropic lake value.
        eta(lp) = -0.20D0
      end if
      ! Put winter surface water a bit colder and summer or tropical
      ! surface water a little warmer to ease spinup nudging in the
      ! correct direction the profile.
      if ( abs(lat(i)) > 25.0 ) then
        tlak(lp,1) = min(max(tgrd(i)-tzero,4.0D0),20.0D0)
        tlak(lp,2) = min(max(tlak(lp,1)-d_half,4.0D0),20.0D0)
        if ( idep(lp) >= 3 ) then
          do k = 3 , idep(lp)
            tlak(lp,k) = min(max(tlak(lp,k-1)-d_half,4.0D0),20.0D0)
          end do
        end if
      else
        ! This needs tuning for tropical lakes.
        ! Lake Malawi bottom temperature can be as high as 22.75 Celsius
        ! with surface as hot as 25.5 degrees.
        tlak(lp,1) = min(max(tgrd(i)-tzero,18.0D0),25.0D0)
        tlak(lp,2) = min(max(tlak(lp,1)-d_half,18.0D0),25.0D0)
        if ( idep(lp) >= 3 ) then
          do k = 3 , idep(lp)
            tlak(lp,k) = min(max(tlak(lp,k-1)-d_half,18.0D0),25.0D0)
          end do
        end if
      end if
      ! Ice over lake from the start
      if ( mask(i) == 4 ) then
        tlak(lp,1) = -2.0D0
        tlak(lp,2) = -1.5D0
        hi(lp) = 0.20D0 ! Set this to 20 cm of ice
        sfice(i) = hi(lp)
        sncv(i) = 1.0D0 ! Put 1 cm of snow on top.
      else
        sfice(i) = d_zero
        sncv(i) = d_zero
      end if
    end do
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine initlake
!
  subroutine lakedrv
    implicit none
    real(rk8) :: flwx , fswx , hsen , prec , qs , tgl , tl , vl , zl
    real(rk8) :: xl , toth , rhs
    real(rk8) :: age , age1 , age2 , arg , arg2 , cdr , cdrmin , cdrn
    real(rk8) :: cdrx , clead , dela , dela0 , delq , dels , delt
    real(rk8) :: fact , factuv , qgrd , qgrnd , qice , rhosw , rhosw3
    real(rk8) :: rib , ribd , ribl , ribn , scrat , tage , tgrnd
    real(rk8) :: sold , vspda , u1
    real(rk8) , dimension(ndpmax) :: tp

    integer(ik4) :: lp , i
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'lakedrv'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
!
    if ( nlakep == 0 ) then
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
      return
    end if

    do lp = 1 , nlakep
      i = ilp(lp)
      tl = sts(i)
      sold = sncv(i)
      vl = dsqrt(usw(i)**2+vsw(i)**2)
      zl = ht(i)
      qs = qv(i)
      fswx = rswf(i)
      flwx = -d_one*rlwf(i)
      prec = prcp(i)*dtlake
      hsen = -d_one*sent(i)
      xl = lat(i)
      tp = tlak(lp,:)

      call lake( dtlake,tl,vl,zl,qs,fswx,flwx,hsen,xl, &
                 tgl,prec,idep(lp),eta(lp),hi(lp),sfice(i), &
                 sncv(i),evpr(i),tp,sfps(i) )

      tlak(lp,:) = tp
      tgb(i)   = tgl
      tgrd(i)  = tgl
      tgbrd(i) = tgl
      rhs = sfps(i)/(rgas*sts(i))
      qgrd = pfqsat(sfps(i),tgrd(i))
      delt = sts(i) - tgrd(i)
      delq = (qs - qgrd)

      if ( sfice(i) <= iceminh ) then
        mask(i) = 3
        sfice(i) = d_zero
        sncv(i) = d_zero
        snag(i) = d_zero
        sm(i) = d_zero
        ribd = usw(i)**2 + vsw(i)**2 + wtur**2
        vspda = dsqrt(ribd)
        cdrn = (vonkar/dlog(ht(i)/zoce))**2
        ribn = ht(i)*egrav*(delt/sts(i))
        rib = ribn/ribd
        if ( rib < d_zero ) then
          cdrx = cdrn*(d_one+24.5D0*dsqrt(-cdrn*rib))
        else
          cdrx = cdrn/(d_one+11.5D0*rib)
        end if
        cdrmin = dmax1(0.25D0*cdrn,6.0D-4)
        if ( cdrx < cdrmin ) cdrx = cdrmin
        drag(i) = cdrx*vspda*rhs
        evpr(i) = -drag(i)*delq
        sent(i) = -drag(i)*cpd*delt
      else
        mask(i) = 4 
        tgbrd(i) = tzero - d_two
        ! Reduce sensible heat flux for ice presence
        toth = sfice(i) + sncv(i)
        if ( toth > href ) then
          sent(i) = sent(i) * (href/toth)**steepf
        end if

        if ( sncv(i) < dlowval ) then
          sncv(i) = d_zero
          snag(i) = d_zero
        else
          arg = 5.0D3*(d_one/tzero-d_one/tgrd(i))
          age1 = dexp(arg)
          arg2 = dmin1(d_zero,d_10*arg)
          age2 = dexp(arg2)
          tage = age1 + age2 + age3
          dela0 = 1.0D-6*dtocn
          dela = dela0*tage
          dels = d_r10*dmax1(d_zero,sncv(i)-sold)
          snag(i) = (snag(i)+dela)*(d_one-dels)
          if ( snag(i) < dlowval ) snag(i) = d_zero
        end if
        if ( sncv(i) > 800.0D0 ) snag(i) = d_zero
        age = (d_one-d_one/(d_one+snag(i)))
        scrat = sncv(i)*0.01D0/(d_one+d_three*age)
        scvk(i) = scrat/(0.1D0+scrat)
        cdrn = (vonkar/dlog(ht(i)/zlnd))**2
        if ( delt < d_zero ) then
          u1 = wtur + d_two*dsqrt(-delt)
        else
          u1 = wtur
        end if
        ribd = usw(i)**2 + vsw(i)**2 + u1**2
        vspda = dsqrt(ribd)
        ribn = ht(i)*egrav*(delt/sts(i))
        rib = ribn/ribd
        if ( rib < d_zero ) then
          cdr = cdrn*(d_one+24.5D0*dsqrt(-cdrn*rib))
        else
          cdr = cdrn/(d_one+11.5D0*rib)
        end if
        cdrmin = dmax1(0.25D0*cdrn,6.0D-4)
        if ( cdr < cdrmin ) cdr = cdrmin
        rhosw = 0.10D0*(d_one+d_three*age)
        rhosw3 = rhosw**3
        cdrn = (vonkar/dlog(ht(i)/zsno))**2
        ribl = (d_one-271.5D0/sts(i))*ht(i)*egrav/ribd
        if ( ribl < d_zero ) then
          clead = cdrn*(d_one+24.5D0*dsqrt(-cdrn*ribl))
        else
          clead = cdrn/(d_one+11.5D0*rib)
        end if
        cdrx = (d_one-aarea)*cdr + aarea*clead
        rhs = sfps(i)/(rgas*sts(i))
        drag(i) = cdrx*vspda*rhs
        qice = 3.3D-3 * stdp/sfps(i)
        qgrnd = ((d_one-aarea)*cdr*qgrd + aarea*clead*qice)/cdrx
        tgrnd = ((d_one-aarea)*cdr*tgrd(i) + aarea*clead*(tzero-1.8D0))/cdrx
        fact = -drag(i)
        delt = sts(i) - tgrnd
        delq = qs - qgrnd
      end if
      if ( dabs(sent(i)) < dlowval ) sent(i) = d_zero
      if ( dabs(evpr(i)) < dlowval ) evpr(i) = d_zero
      fact = dlog(ht(i)*d_half)/dlog(ht(i)/zoce)
      factuv = dlog(ht(i)*d_r10)/dlog(ht(i)/zoce)
      u10m(i) = usw(i)*(d_one-factuv)
      v10m(i) = vsw(i)*(d_one-factuv)
      taux(i) = dmissval
      tauy(i) = dmissval
      t2m(i) = sts(i) - delt*fact
      q2m(i) = qs - delq*fact
    end do

#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine lakedrv
!
!-----------------------------------------------------------------------
!
  subroutine lake(dtlake,tl,vl,zl,ql,fsw,flw,hsen,xl,tgl,  &
                  prec,ndpt,eta,hi,aveice,hsnow,evl,tprof,ps)
    implicit none
    real(rk8) :: dtlake , evl , aveice , hsen , hsnow , flw , &
               prec , ql , fsw , tl , tgl , vl , zl , eta , hi , xl , ps
    real(rk8) , dimension(ndpmax) :: tprof
    integer(ik4) :: ndpt
    intent (in) hsen , ql , tl , vl , zl , ps
    intent (in) ndpt , eta
    intent (out) tgl
    intent (inout) tprof , evl , aveice , hsnow
    real(rk8) :: ai , ea , ev , hs , ld , lu , qe , qh , tac , tk , u2
    ! zo: surface roughness length
    real(rk8) , parameter :: zo = 0.001D0
    real(rk8) , parameter :: z2 = d_two
    real(rk8) , parameter :: tcutoff = -0.001D0
    real(rk8) , parameter :: twatui = 1.78D0
    logical , parameter :: lfreeze = .true.
    integer(ik4) , parameter :: kmin = 1
    integer(ik4) , parameter :: kmax = 200
!
    ! interpolate winds at z1 m to 2m via log wind profile
    u2 = vl*dlog(z2/zo)/dlog(zl/zo)
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

    else
      ! Lake ice
      ! convert mixing ratio to air vapor pressure
      ea  = ql*88.0D0/(ep2+0.378D0*ql)
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

      call lakeice(dtlake,fsw,ld,tac,u2,ea,hs,hi,ai,ev,prec,ps,tprof)
      if ( .not. lfreeze ) tprof(1) = twatui

      evl    = ev/secph       ! convert evl  from mm/hr to mm/sec
      aveice = ai
      hsnow  = hs*d_100      ! convert back. See Above.
      if (aveice < dlowval) then
        aveice = d_zero
        hsnow = d_zero
      end if
    end if
    tgl = tprof(1) + tzero
  end subroutine lake
!
!-----------------------------------------------------------------------
!
  subroutine lakeeddy(ndpt,dtlake,u2,xl,tprof)

    ! Computes density and eddy diffusivity

    implicit none
!
    integer(ik4) , intent (in) :: ndpt
    real(rk8) , intent (in) :: dtlake , u2 , xl
    real(rk8) , dimension(ndpmax) , intent (in) :: tprof
!
    real(rk8) :: demax , demin , dpdz , ks , n2 , po
    real(rk8) :: zmax , rad , ri , ws , z
    integer(ik4) :: k
!
    ! demin molecular diffusion of heat in water
    demin = hdmw
!
    ! Added to keep numerical stability of code
    demax = .50D0*dz**2/dtlake
    demax = .99D0*demax
!
    do k = 1 , ndpt
      dnsty(k) = d_1000*(d_one-1.9549D-05*(dabs(tprof(k)-4.0D0))**1.68D0)
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
    ! Decay constant of shear velocity - Ekman profile parameter
    if ( xl > 25.0D0 ) then
      ks = 6.6D0*dsqrt(dsin(xl*degrad))*u2**(-1.84D0)
    else
      ks = 0.001
    end if

    ! Ekman layer depth where eddy diffusion happens
    zmax = dble(ceiling(surf+40.0D0/(vonkar*ks)))

    ! Surface shear velocity
    ws = 0.0012D0*u2

    ! Inverse of turbulent Prandtl number
    po = d_one

    do k = 1 , ndpt - 1

      ! Actual depth from surface
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

      ! Brunt Vaisala frequency squared : we do not mind stability,
      ! we just look for energy here.
      ! n2 = dabs((dpdz/dnsty(k))*egrav)
      n2 = (dpdz/dnsty(k))*egrav
      if (dabs(n2) < dlowval) then
        de(k) = demin
        cycle
      end if

      ! Richardson number estimate
      rad = d_one+40.0D0*n2*((vonkar*z)/(ws*dexp(-ks*z)))**2
      if (rad < d_zero) rad = d_zero
      ri = (-d_one+dsqrt(rad))/20.0D0

      ! Total diffusion coefficient for heat: molecular + eddy (Eqn 42)
      de(k) = demin + vonkar*ws*z*po*dexp(-ks*z) / (d_one+37.0D0*ri**2)
      if ( de(k) < demin ) de(k) = demin
      if ( de(k) > demax ) de(k) = demax

    end do
    de(ndpt) = demin
  end subroutine lakeeddy
!
!-----------------------------------------------------------------------
!
  subroutine laketemp(ndpt,dtlake,fsw,flw,qe,qh,eta,tprof)
    implicit none
    integer(ik4) , intent(in) :: ndpt
    real(rk8) , intent(in) :: dtlake , eta , flw , qe , qh , fsw
    real(rk8) , dimension(ndpmax) , intent(inout) :: tprof
    real(rk8) :: bot , dt1 , dt2 , top
    integer(ik4) :: k

    ! Computes temperature profile
    ! solve differential equations of heat transfer

    tt(1:ndpt) = tprof(1:ndpt)

    dt1 = (fsw*(d_one-dexp(eta*surf))+(flw+qe+qh)) / &
            (surf*dnsty(1)*cpw)
    dt2 = -de(1)*(tprof(1)-tprof(2))/surf
    tt(1) = tt(1) + (dt1+dt2)*dtlake

    do k = 2 , ndpt - 1
      top = (surf+(k-2)*dz)
      bot = (surf+(k-1)*dz)
      dt1 = fsw*(dexp(eta*top)-dexp(eta*bot))/(dz*dnsty(k)*cpw)
      dt2 = (de(k-1)*(tprof(k-1)-tprof(k))    -    &
             de(k)  *(tprof(k)  -tprof(k+1))) / dz
      tt(k) = tt(k) + (dt1+dt2)*dtlake
    end do

    top = (surf+(ndpt-2)*dz)
    dt1 = fsw*dexp(eta*top)/(dz*dnsty(ndpt)*cpw)
    dt2 = de(ndpt-1)*(tprof(ndpt-1)-tprof(ndpt))/dz
    tt(ndpt) = tt(ndpt) + (dt1+dt2)*dtlake

    do k = 1 , ndpt
      tprof(k) = tt(k)
      dnsty(k) = d_1000*(d_one-1.9549D-05*(dabs(tprof(k)-4.0D0))**1.68D0)
    end do
  end subroutine laketemp
!
!-----------------------------------------------------------------------
!
  subroutine lakemixer(kmin,kmax,ndpt,tprof)
    implicit none
    integer(ik4) , intent(in) :: ndpt , kmin , kmax
    real(rk8) , intent(inout) , dimension(ndpmax) :: tprof
    real(rk8) :: avet , avev , tav , vol
    integer(ik4) :: k , k2
! 
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
          dnsty(k2) = d_1000*(d_one-1.9549D-05*(dabs(tav-4.0D0))**1.68D0)
        end do
      end if

    end do ! K loop

    tprof(kmin:ndpt) = tt(kmin:ndpt)
  end subroutine lakemixer
!
!-----------------------------------------------------------------------
!
  subroutine lakeice(dtx,fsw,ld,tac,u2,ea,hs,hi,aveice,evl,prec,ps,tprof)
    implicit none
    real(rk8) , intent(in) :: ea , ld , prec , tac , u2 , dtx , ps
    real(rk8) , intent(out) :: evl
    real(rk8) , intent(inout) :: hi , aveice , hs , fsw
    real(rk8) , dimension(ndpmax) , intent(inout) :: tprof
    real(rk8) :: di , ds , f0 , f1 , khat , psi , q0 , qpen , t0 , t1 , &
                 t2 , tf , theta , rho
    real(rk8) , parameter :: isurf = 0.6D0
    ! attenuation coeff for ice in visible band (m-1)
    real(rk8) , parameter :: lami1 = 1.5D0
    ! attenuation coeff for ice in infrared band (m-1)
    real(rk8) , parameter :: lami2 = 20.0D0
    ! attenuation coeff for snow in visible band (m-1)
    real(rk8) , parameter :: lams1 = 6.0D0
    ! attenuation coeff for snow in infrared band (m-1)
    real(rk8) , parameter :: lams2 = 20.0D0
    ! thermal conductivity of ice (W/m/C)
    real(rk8) , parameter :: ki = 2.3D0
    ! thermal conductivity of snow (W/m/C)
    ! 0.078 W m−1 K−1 for new snow 
    ! 0.290 W m−1 K−1 for an ubiquitous wind slab
    ! the average bulk value for the full snowpack is 0.14 W m−1 K−1
    ! Sturm, Perovich and Holmgren
    ! Journal of Geophysical Research: Oceans (1978–2012)
    real(rk8) , parameter :: ks = 0.14D0
    ! heat flux from water to ice (w m-2)
    real(rk8) , parameter :: qw = 1.389D0
    ! latent heat of fusion (J/kg)
    real(rk8) , parameter :: li = 334.0D03
    ! drag coefficient for the turbulent momentum flux.
    real(rk8) , parameter :: cd = 0.001D0

    if ( (tac <= d_zero) .and. (aveice > d_zero) ) then
      ! National Weather Service indicates the average snowfall is in a
      ! ratio of 10 inches of snow to 1 inch of equivalent rainfall
      ! So this is convert to mm and to snowfall ( prec*10/1000 )
      hs = hs + prec*d_r100  ! convert prec(mm) to depth(m)
    end if
    if ( hs < dlowval ) hs = d_zero

    ! temperature of ice/snow surface
    t0 = tprof(1)
    ! temperature of the water under the snow
    tf = -d_two
    ! approximate density of air (1 kg/m3)
    rho = rhoh2o*d_r1000

    khat = (ki*hs+ks*hi)/(ki*ks)
    theta = cpd*rho*cd*u2
    psi = wlhv*rho*cd*u2*ep2/ps
    evl = psi*(eomb(t0)-ea)/(wlhv*rho)
    ! amount of radiation that penetrates through the ice (W/m2)
    qpen = fsw * 0.7D0 * &
       ((d_one-dexp(-lams1*hs))                   / (ks*lams1) + &
        (dexp(-lams1*hs))*(d_one-dexp(-lami1*hi)) / (ki*lami1)) + &
           fsw * 0.3D0 * &
       ((d_one-dexp(-lams2*hs))                   / (ks*lams2) + &
        (dexp(-lams2*hs))*(d_one-dexp(-lami2*hi)) / (ki*lami2))
    ! radiation absorbed at the ice surface
    fsw = fsw - qpen
    t1 = -20.0D0
    f0 = f(t0)
    f1 = f(t1)
    do
      t2 = t1 - (t1-t0)*f1/(f1-f0)
      if ( (t2-t1)/t1 >= 0.001D0 ) then
        t0 = t1
        t1 = t2
        f0 = f1
        f1 = f(t1)
        cycle
      end if
      exit
    end do

    t0 = t2
    if ( t0 >= tf ) then
      if ( hs > d_zero ) then
        ds = dtx *                                                  &
             ( (-ld + 0.97D0*sigm*t4(tf) + psi * (eomb(tf)-ea) +    &
                theta*(tf-tac)-fsw) - d_one/khat * (tf-t0+qpen) ) / &
              (rhosnow*li)
        if ( ds > d_zero ) ds = d_zero
        hs = hs + ds * 10.0D0
        if ( hs < d_zero ) then
          hs = d_zero
        end if
      end if
      if ( (dabs(hs) < dlowval) .and. (aveice > d_zero) ) then
        di = dtx *                                                  &
            ( (-ld + 0.97D0*sigm*t4(tf) + psi * (eomb(tf)-ea)     + &
              theta*(tf-tac)-fsw) - d_one/khat * (tf-t0+qpen) ) /   &
             (rhoice*li)
        if ( di > d_zero ) di = d_zero
        hi = hi + di
      end if
    else
      q0 = -ld + 0.97D0*sigm*t4(t0) + psi*(eomb(t0)-ea) + &
           theta*(t0-tac) - fsw
      qpen = fsw*0.7D0*(d_one-dexp(-(lams1*hs+lami1*hi))) + &
             fsw*0.3D0*(d_one-dexp(-(lams2*hs+lami2*hi)))
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
      tprof(1) = t0
    end if

    contains

    pure real(rk8) function t4(x)
      implicit none
      real(rk8) , intent(in) :: x
      t4 = (x+tzero)**4
    end function t4

    ! Computes air vapor pressure as a function of temp (in K)
    pure real(rk8) function eomb(x)
      implicit none
      real(rk8) , intent(in) :: x
      real(rk8) :: tr1
      tr1 = d_one - (tboil/(x+tzero))
      eomb = stdpmb*dexp(13.3185D0*tr1 - 1.976D0*tr1**2 - &
                         0.6445D0*tr1**3 - 0.1299D0*tr1**4)
     end function eomb
    pure real(rk8) function f(x)
      implicit none
      real(rk8) , intent(in) :: x
      f = (-ld + 0.97D0*sigm*t4(x) + psi*(eomb(x)-ea) + &
           theta*(x-tac)-fsw) - d_one/khat*(qpen+tf-x)
    end function f

  end subroutine lakeice
!
!-----------------------------------------------------------------------
!
  subroutine lake_fillvar_real8_1d(ivar,rvar,idir)
    implicit none
    integer(ik4) , intent(in) :: ivar , idir
    real(rk8) , intent(inout) , pointer , dimension(:) :: rvar
    real(rk8) , pointer , dimension(:) :: p

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
!
  subroutine lake_fillvar_real8_2d(ivar,rvar,idir)
    implicit none
    integer(ik4) , intent(in) :: ivar , idir
    real(rk8) , intent(inout) , pointer , dimension(:,:) :: rvar
    real(rk8) , pointer , dimension(:,:) :: p
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
        p(:,k) = pack(rvar(:,k),lakmsk)-tzero
      end do
    else
      rvar = dmissval
      do k = 1 , ndpmax
        rvar(:,k) = unpack(p(:,k),lakmsk,rvar(:,k))+tzero
      end do
    end if
  end subroutine lake_fillvar_real8_2d
!
end module mod_ocn_lake
