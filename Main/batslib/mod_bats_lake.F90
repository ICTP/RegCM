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
  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_service
  use mod_bats_common
  use mod_bats_internal
  use mod_bats_mppio
!
  private
!
  public :: initlake , lakedrv
  public :: lakesav_i , lakesav_o
!
  real(rk8) , dimension(ndpmax) :: de , dnsty , tt , ttx
!
  ! surface thickness
  real(rk8) , parameter :: surf = d_one
  ! vertical grid spacing in m
  real(rk8) , parameter :: dz = surf
  ! minimum ice depth in mm: less that this is removed
  real(rk8) , parameter :: iceminh = 1.0D0
  ! reference hgt in mm for latent heat removal from ice
  real(rk8) , parameter :: href = d_two * iceminh
  ! steepness factor of latent heat removal
  real(rk8) , parameter :: steepf = 1.0D0  ! Tuning needed !
!
  contains
!
!-----------------------------------------------------------------------
!
  subroutine initlake
    implicit none
    integer(ik4) :: i, j, n
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'initlake'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    ! initialize hostetler lake model

    idep(:,:,:)   = 0
    hi(:,:,:)     = dmissval
    aveice(:,:,:) = dmissval
    hsnow(:,:,:)  = dmissval
    eta(:,:,:)    = dmissval
    tlak(:,:,:,:) = dmissval

    do i = ici1 , ici2
      do j = jci1 , jci2
        do n = 1 , nnsg
          if ( iveg1(n,j,i) == 14 ) then
            idep(n,j,i) = idint(dmax1(d_two,dmin1(dhlake1(n,j,i), &
                                  dble(ndpmax)))/dz)
            if ( ldmsk1(n,j,i) == 2 ) then
              tlak(n,j,i,1) = -2.0D0
              tlak(n,j,i,2) = -2.0D0
              aveice(n,j,i) = d_10
              hi(n,j,i) = d_one
              hsnow(n,j,i) = d_zero
            end if
            ! Azar Zarrin: Fixed unrealistic high ice tickness and
            ! high water temperatures during warm months.
            if (idep(n,j,i) < 50) then
              eta(n,j,i) = 0.5D0
            else if (idep(n,j,i) > 100) then
              eta(n,j,i) = 0.1D0
            else
              eta(n,j,i) = 0.3D0
            end if
            tlak(n,j,i,1:idep(n,j,i)) = 6.0D0
            hi(n,j,i) = 0.01D0
            aveice(n,j,i) = d_zero
            hsnow(n,j,i)  = d_zero
          end if
        end do
      end do
    end do
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine initlake
!
  subroutine lakedrv
    implicit none
    real(rk8) :: flwx , fswx , hsen , prec , ql , tgl , tl , vl , zl , &
                xl , evp , toth
    integer(ik4) :: i , j , n
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'lakedrv'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
!
    do i = ici1 , ici2
      do j = jci1 , jci2
        do n = 1 , nnsg
          if ( idep(n,j,i) > 1 ) then
            tl = sts(n,j,i)
            vl = dsqrt(usw(j,i)**2+vsw(j,i)**2)
            zl = zh(n,j,i)
            ql = qs(n,j,i)
            fswx = fsw(j,i)
            flwx = -d_one*flw(j,i)
            prec = prcp(n,j,i)*dtbat
            hsen = -d_one*sent(n,j,i)
            evp = evpr(n,j,i)
            if (nnsg == 1) then
              xl = xlat(j,i)
            else
              xl = xlat1(n,j,i)
            end if

            ttx(:) = tlak(n,j,i,:)
            call lake( dtlake,tl,vl,zl,ql,fswx,flwx,hsen,xl, &
                       tgl,prec,idep(n,j,i),eta(n,j,i),  &
                       hi(n,j,i),aveice(n,j,i),          &
                       hsnow(n,j,i),evp,ttx )
            tlak(n,j,i,:) = ttx(:)

            ! Feed back ground temperature
            tgrd(n,j,i) = tgl
            tgbrd(n,j,i) = tgl

            if ( aveice(n,j,i) <= iceminh ) then
              ldmsk1(n,j,i) = 0 
              lveg(n,j,i) = 14
              sfice(n,j,i) = d_zero
              sncv(n,j,i) = d_zero
              snag(n,j,i) = d_zero
            else
              ldmsk1(n,j,i) = 2 
              lveg(n,j,i) = 12
              sfice(n,j,i) = aveice(n,j,i)  !  units of ice = mm
              sncv(n,j,i)  = hsnow(n,j,i)   !  units of snw = mm
              evpr(n,j,i) = evp                 !  units of evp = mm/sec
              ! Reduce sensible heat flux for ice presence
              toth = sfice(n,j,i) + sncv(n,j,i)
              if ( toth > href ) then
                sent(n,j,i) = sent(n,j,i) * (href/toth)**steepf
              end if
              if ( dabs(sent(n,j,i)) < dlowval ) sent(n,j,i) = d_zero
              if ( dabs(evpr(n,j,i)) < dlowval ) evpr(n,j,i) = d_zero
            end if
          end if
        end do
      end do
    end do
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine lakedrv
!
!-----------------------------------------------------------------------
!
  subroutine lake(dtlake,tl,vl,zl,ql,fsw,flw,hsen,xl,tgl,  &
                  prec,ndpt,eta,hi,aveice,hsnow,evl,tprof)
    implicit none
    real(rk8) :: dtlake , evl , aveice , hsen , hsnow , flw , &
               prec , ql , fsw , tl , tgl , vl , zl , eta , hi , xl
    real(rk8) , dimension(ndpmax) :: tprof
    integer(ik4) :: ndpt
    intent (in) hsen , ql , tl , vl , zl
    intent (in) ndpt , eta
    intent (out) tgl
    intent (inout) evl , aveice , hsnow
    intent (inout) tprof
    real(rk8) :: ai , ea , ev , hs , ld , lu , qe , qh , tac , tk , u2
    ! zo: surface roughness length
    real(rk8) , parameter :: zo = 0.001D0
    real(rk8) , parameter :: z2 = d_two
    real(rk8) , parameter :: tcutoff = -0.001D0
    real(rk8) , parameter :: twatui = 1.78D0
    logical , parameter :: lfreeze = .false.
    integer(ik4) , parameter :: kmin = 1
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'lake'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
!
    ! interpolate winds at z1 m to 2m via log wind profile
    u2 = vl*dlog(z2/zo)/dlog(zl/zo)
    if ( u2 < d_half ) u2 = d_half

    ! Check if conditions not exist for lake ice
    if ( (aveice < 1.0D-8) .and. (tprof(1) > tcutoff) ) then

      ! Graziano: removed hlat. It is calculated from evaporation
      qe = -d_one*evl*wlhv
      qh = hsen

      ! Calculate eddy diffusivities
      call lakeeddy(ndpt,dtlake,u2,xl,tprof)

      ! Lake temperature calc using sensible and latent heats
      call laketemp(ndpt,dtlake,fsw,flw,qe,qh,eta,tprof)

      ! Convective mixer
      call lakemixer(kmin,ndpt,tprof)

      hi     = 0.01D0
      aveice = d_zero
      hsnow  = d_zero

    else
      ! Lake ice
      ! convert mixing ratio to air vapor pressure
      ea  = ql*88.0D0/(ep2+0.378D0*ql)
      tac = tl - tzero
      tk  = tzero + tprof(1)
      lu  = -emsw*sigm*tk**4
      ld  = flw - lu
      ev  = evl*secph         ! convert to mm/hr
      ai  = aveice / d_1000   ! convert to m
      hs  = hsnow / d_100     ! convert to m

      call lakeice(dtlake,fsw,ld,tac,u2,ea,hs,hi,ai,ev,prec,tprof)
      if ( .not. lfreeze ) tprof(1) = twatui

      evl    = ev/secph       ! convert evl  from mm/hr to mm/sec
      aveice = ai*d_1000      ! convert ice  from m to mm
      hsnow  = hs*d_100       ! convert snow from m depth to mm h20
      if (aveice < dlowval) aveice = d_zero
      if (hsnow < dlowval) hsnow = d_zero

    end if
    tgl = tprof(1) + tzero
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
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
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'lakeeddy'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
!
    ! demin molecular diffusion of heat in water
    demin = hdmw
!
    ! Added to keep numerical stability of code
    demax = .50D0*dz**2/dtlake
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
    ! Decay constant of shear velocity - Ekman profile parameter
    ks = 6.6D0*dsqrt(dsin(xl*degrad))*u2**(-1.84D0)

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
      de(k) = demin + vonkar*ws*z*po*dexp(-ks*z) / &
                      (d_one+37.0D0*ri**2)
      if ( de(k) < demin ) de(k) = demin
      if ( de(k) > demax ) de(k) = demax

    end do
    de(ndpt) = demin
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
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
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'laketemp'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    ! Computes temperature profile
    ! solve differential equations of heat transfer

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
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine laketemp
!
!-----------------------------------------------------------------------
!
  subroutine lakemixer(kmin,ndpt,tprof)
    implicit none
    integer(ik4) , intent(in) :: ndpt , kmin
    real(rk8) , intent(inout) , dimension(ndpmax) :: tprof
    real(rk8) :: avet , avev , tav , vol
    integer(ik4) :: k , k2
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'lakemixer'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
! 
    ! Simulates convective mixing
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
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine lakemixer
!
!-----------------------------------------------------------------------
!
  subroutine lakeice(dtx,fsw,ld,tac,u2,ea,hs,hi,aveice,evl,prec,tprof)
    implicit none
    real(rk8) :: ea , evl , hi , aveice , hs , fsw , &
               ld , prec , tac , u2 , dtx
    real(rk8) , dimension(ndpmax) :: tprof
    intent (in) dtx , ea , ld , prec , tac , u2
    intent (out) evl
    intent (inout) hi , aveice , hs , fsw , tprof
    real(rk8) :: di , ds , f0 , f1 , khat , psi , q0 , qpen , t0 , t1 , &
               t2 , tf , theta , rho , xlexpc
    real(rk8) :: xea , xeb , xec
    integer(ik4) :: nits
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
    real(rk8) , parameter :: ks = 0.31D0
    ! standard atmospheric pressure (hPa) ????
    real(rk8) , parameter :: atm = 950.0D0
    ! heat flux from water to ice (w/m2) ???
    real(rk8) , parameter :: qw = 1.389D0
    ! latent heat of fusion (J/kg)
    real(rk8) , parameter :: li = 334.0D03
    ! drag coefficient for the turbulent momentum flux.
    real(rk8) , parameter :: cd = 0.001D0
    ! Maximum exponent
    real(rk8) , parameter :: minexp = -25.0D0
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'lakeice'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

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
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif

    contains

    function t4(x)
      implicit none
      real(rk8) :: t4
      real(rk8) , intent(in) :: x
      t4 = (x+tzero)**4
    end function t4
    ! Computes air vapor pressure as a function of temp (in K)
    function tr1(x)
      implicit none
      real(rk8) :: tr1
      real(rk8) , intent(in) :: x
      tr1 = d_one - (tboil/(x+tzero))
    end function tr1
    function eomb(x)
      implicit none
      real(rk8) :: eomb
      real(rk8) , intent(in) :: x
      eomb = stdpmb*dexp(13.3185D0*tr1(x)-1.976D0*tr1(x)**2   &
             -0.6445D0*tr1(x)**3- 0.1299D0*tr1(x)**4)
     end function eomb
    function f(x)
      implicit none
      real(rk8) :: f
      real(rk8) , intent(in) :: x
      f = (-ld+0.97D0*sigm*t4(x)+psi*(eomb(x)-ea)+theta*(x-tac)-fsw)  &
          - d_one/khat*(qpen+tf-x)
    end function f

  end subroutine lakeice
!
!-----------------------------------------------------------------------
!
  subroutine lakesav_o(iutl)
    implicit none
    integer(ik4) , intent(in) :: iutl
    integer(ik4) :: i , j , k , n
    write (iutl) idep_io
    do i = iout1 , iout2
      do j = jout1 , jout2
        do n = 1 , nnsg
          if ( idep_io(n,j,i) > 0 ) then
            write(iutl) eta_io(n,j,i), hi_io(n,j,i), &
                 aveice_io(n,j,i), hsnow_io(n,j,i),  &
                 (tlak_io(n,j,i,k),k=1,idep_io(n,j,i))  
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
    integer(ik4) , intent(in) :: iutl
    integer(ik4) :: i , j , k , n
    idep_io   = 0
    hi_io     = dmissval
    aveice_io = dmissval
    hsnow_io  = dmissval
    eta_io    = dmissval
    tlak_io   = dmissval
    read (iutl) idep_io
    do i = iout1 , iout2
      do j = jout1 , jout2
        do n = 1 , nnsg
          if ( idep_io(n,j,i) > 0 ) then
            read(iutl) eta_io(n,j,i), hi_io(n,j,i), &
                 aveice_io(n,j,i), hsnow_io(n,j,i), &
                 (tlak_io(n,j,i,k),k=1,idep_io(n,j,i))  
          end if
        end do
      end do
    end do
  end subroutine lakesav_i
!
end module mod_bats_lake
