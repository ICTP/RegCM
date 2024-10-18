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

! Version 1.0 released by Yi-Chuan Lu on May 20, 2022.
!
! When using this code, please cite:
!
! @article{20heatindex,
!   Title   = {Extending the Heat Index},
!   Author  = {Yi-Chuan Lu and David M. Romps},
!   Journal = {Journal of Applied Meteorology and Climatology},
!   Year    = {2022},
!   Volume  = {61},
!   Number  = {10},
!   Pages   = {1367--1383},
!   Year    = {2022},
! }
!
! This headindex function returns the Heat Index in Kelvin. The inputs are:
! - T, the temperature in Kelvin
! - RH, the relative humidity, which is a value from 0 to 1

module mod_heatindex

  use mod_realkinds
  use mod_stdio
  use mod_message

  implicit none

  private

  ! w/m^2/k^4 , boltzmann constant
  real(rkx) , parameter :: sigma    = 5.67e-8_rkx
  ! k         , vapor temperature at triple point
  real(rkx) , parameter :: ttrip    = 273.16_rkx
  ! pa        , vapor pressure at triple point
  real(rkx) , parameter :: ptrip    = 611.65_rkx
  real(rkx) , parameter :: e0v      = 2.3740e6_rkx  ! J/kg
  real(rkx) , parameter :: e0s      = 0.3337e6_rkx  ! J/kg
  real(rkx) , parameter :: rgasa    = 287.04_rkx    ! J/kg/K
  real(rkx) , parameter :: rgasv    = 461.0_rkx     ! J/kg/K
  real(rkx) , parameter :: cva      = 719.0_rkx     ! J/kg/K
  real(rkx) , parameter :: cvv      = 1418.0_rkx    ! J/kg/K
  real(rkx) , parameter :: cvl      = 4119.0_rkx    ! J/kg/K
  real(rkx) , parameter :: cvs      = 1861.0_rkx    ! J/kg/K
  real(rkx) , parameter :: cpa      = cva + rgasa
  real(rkx) , parameter :: cpv      = cvv + rgasv
  !           , emissivity of surface
  real(rkx) , parameter :: reps  = 0.97_rkx
  ! kg        , mass of average us adults, FRYAR2018
  real(rkx) , parameter :: m        = 83.6_rkx
  ! m         , height of average us adults, FRYAR2018
  real(rkx) , parameter :: h        = 1.69_rkx
  ! m^2       , dubois formula, parson2014
  real(rkx) , parameter :: a        = 0.202_rkx*(m**0.425_rkx)*(h**0.725_rkx)
  ! J/kg/K    , specific heat capacity of core, Gagge1972
  real(rkx) , parameter :: cpc      = 3492.0_rkx
  !           , heat capacity of core
  real(rkx) , parameter :: c        = m*cpc/a
  ! Pa/K      , zf/rf
  real(rkx) , parameter :: r        = 124.0_rkx
  ! W/m^2     , metabolic rate per skin area
  real(rkx) , parameter :: q        = 180.0_rkx
  !           , vapor saturation pressure level of saline solution
  real(rkx) , parameter :: phi_salt = 0.9_rkx
  ! K         , core temeprature
  real(rkx) , parameter :: tc        = 310.0_rkx
  !           , core vapor pressure, phi_salt*pvstar(tc)
  real(rkx) , parameter :: pc        = 5612.299124343398_rkx
  ! Pa        , atmospheric pressure
  real(rkx) , parameter :: p         = 1.013e5_rkx
  ! kg/J      , "inhaled mass" / "metabolic rate"
  real(rkx) , parameter :: eta       = 1.43e-6_rkx
  ! Pa        , reference air vapor pressure in regions III, IV, V, VI,
  !             chosen by Steadman
  real(rkx) , parameter :: pa0       = 1.6e3_rkx
  ! Pa m^2/W  , mass transfer resistance through air, exposed part of skin
  real(rkx) , parameter :: za        = 60.6_rkx/17.4_rkx
  ! Pa m^2/W  , mass transfer resistance through air, clothed part of skin
  real(rkx) , parameter :: za_bar    = 60.6_rkx/11.6_rkx
  ! Pa m^2/W  , mass transfer resistance through air, when being naked
  real(rkx) , parameter :: za_un     = 60.6_rkx/12.3_rkx
  !           , tolorence of the root solver for heatindex
  real(rkx) , parameter :: errt      = 1.0e-8_rkx
  !           , tolorence of the root solver
  real(rkx) , parameter :: err       = 1.0e-8_rkx

  ! maximum number of iteration of the root solver
  integer, parameter :: maxiter = 100

  integer , parameter :: eqvar_indx_1 = 1 ! phi
  integer , parameter :: eqvar_indx_2 = 2 ! rf
  integer , parameter :: eqvar_indx_3 = 3 ! rs
  integer , parameter :: eqvar_indx_4 = 4 ! rs*
  integer , parameter :: eqvar_indx_5 = 5 ! dtcdt

  public  :: heatindex

  type eqvar
    integer :: indx
    real(rkx) , dimension(4) :: var
  end type eqvar

  contains

  pure real(rkx) function pvstar(t)
    implicit none
    real(rkx) , intent(in) :: t
    if ( t <= 0.0_rkx ) then
      pvstar = 0.0_rkx
    else if ( t < ttrip ) then
      pvstar = ptrip * (t/ttrip)**((cpv-cvs)/rgasv) * &
        exp( (e0v + e0s -(cvv-cvs)*ttrip)/rgasv * (1.0_rkx/ttrip - 1.0_rkx/t) )
    else
      pvstar = ptrip * (t/ttrip)**((cpv-cvl)/rgasv) * &
        exp( (e0v       -(cvv-cvl)*ttrip)/rgasv * (1.0_rkx/ttrip - 1.0_rkx/t) )
    end if
  end function pvstar

  pure real(rkx) function le(t)
    implicit none
    real(rkx) , intent(in):: t
    le = (e0v + (cvv-cvl)*(t-ttrip) + rgasv*t)
  end function le

  pure real(rkx) function qv(ta,pa) ! respiratory heat loss, w/m^2
    implicit none
    real(rkx) , intent(in) :: ta , pa
    qv = eta * q *(cpa*(tc-ta) + le(tc)*rgasa/(p*rgasv) * ( pc-pa ) )
  end function qv

  ! mass transfer resistance through skin, pa m^2/w
  pure real(rkx) function zs(rs)
    implicit none
    real(rkx) , intent(in) :: rs
    if ( rs == 0.0387_rkx ) then
      zs = 52.1_rkx
    else
      zs = 6.0e8_rkx * rs**5
    end if
  end function zs

  ! heat transfer resistance through air, exposed part of skin, k m^2/w
  pure real(rkx) function ra(ts,ta)
    implicit none
    real(rkx) , intent(in) :: ts , ta
    real(rkx) :: hc , phi_rad , hr
    hc      = 17.4_rkx
    phi_rad = 0.85_rkx
    hr      = reps * phi_rad * sigma* (ts**2 + ta**2)*(ts + ta)
    ra      = 1.0_rkx/(hc+hr)
  end function ra

  ! heat transfer resistance through air, clothed part of skin, k m^2/w
  pure real(rkx) function ra_bar(tf,ta)
    implicit none
    real(rkx), intent(in) :: tf , ta
    real(rkx) :: hc , phi_rad, hr
    hc      = 11.6_rkx
    phi_rad = 0.79_rkx
    hr      = reps * phi_rad * sigma* (tf**2 + ta**2)*(tf + ta)
    ra_bar  = 1.0_rkx/(hc+hr)
  end function ra_bar

  ! heat transfer resistance through air, when being naked, k m^2/w
  pure real(rkx) function ra_un(ts,ta)
    implicit none
    real(rkx) , intent(in) :: ts , ta
    real(rkx) :: hc , phi_rad , hr
    hc      = 12.3_rkx
    phi_rad = 0.80_rkx
    hr      = reps * phi_rad * sigma* (ts**2 + ta**2)*(ts + ta)
    ra_un   = 1.0_rkx/(hc+hr)
  end function ra_un

#ifdef DEBUG
  type(eqvar) function initial_find_eqvar(ta,rh) result(res)
#else
  pure type(eqvar) function initial_find_eqvar(ta,rh) result(res)
#endif
    implicit none
    real(rkx) , intent(in) :: ta , rh
    real(rkx) :: pa , phi , rf , rs , dtcdt , ts , ts_bar
    real(rkx) :: tf , ps , flux1 , flux2 , flux3 , m , m_bar

    pa    = rh*pvstar(ta)
    phi   = 0.84_rkx
    rs    = 0.0387_rkx
    dtcdt = 0.0_rkx

    m     = (pc-pa)/(zs(rs)+za)
    m_bar = (pc-pa)/(zs(rs)+za_bar)
    ts = solve1(ta,pa,rs,max(0.0_rkx,min(tc,ta)-rs*abs(m)), &
                max(tc,ta)+rs*abs(m),    err,maxiter)
    tf = solve2(ta,pa,rs,max(0.0_rkx,min(tc,ta)-rs*abs(m_bar)), &
                max(tc,ta)+rs*abs(m_bar),err,maxiter)
    ! c*dtc/dt when rf=zf=\inf
    flux1 = q-qv(ta,pa)-(1.0_rkx-phi)*(tc-ts)/rs
    ! c*dtc/dt when rf=zf=0
    flux2 = q-qv(ta,pa)-(1.0_rkx-phi)*(tc-ts)/rs - phi*(tc-tf)/rs
    if ( flux1 <= 0.0_rkx ) then ! region I
      phi = 1.0_rkx-(q-qv(ta,pa))*rs/(tc-ts)
      rf  = huge(0.0_rkx)
      res%indx = eqvar_indx_1
    else if ( flux2 <= 0.0_rkx ) then ! region II&III
      ts_bar = tc - (q-qv(ta,pa))*rs/phi + (1.0_rkx/phi -1.0_rkx)*(tc-ts)
      tf = solve3(ta,pa,rs,ts_bar,ta,ts_bar,err,maxiter)
      rf = ra_bar(tf,ta)*(ts_bar-tf)/(tf-ta)
      res%indx = eqvar_indx_2
    else ! region IV,V,VI
      rf = 0.0_rkx
      flux3 = q-qv(ta,pa)-(tc-ta)/ra_un(tc,ta)-(phi_salt*pvstar(tc)-pa)/za_un
      if ( flux3 <= 0.0_rkx ) then ! region IV,V
        ts = solve4(ta,pa,0.0_rkx,tc,err,maxiter)
        rs = (tc-ts)/(q-qv(ta,pa))
        ps = pc - (pc-pa)* zs(rs)/( zs(rs)+za_un)
        res%indx = eqvar_indx_3
        if ( ps > phi_salt * pvstar(ts) ) then  ! region V
          ts = solve5(ta,pa,0.0_rkx,tc,err,maxiter)
          rs = (tc-ts)/(q-qv(ta,pa))
          res%indx = eqvar_indx_4
        end if
      else ! region VI
        rs = 0.0_rkx
        dtcdt = (1.0_rkx/c)* flux3
        res%indx = eqvar_indx_5
      end if
    end if
    res%var = (/ phi,rf,rs,dtcdt /)
  end function initial_find_eqvar

#ifdef DEBUG
  function find_eqvar(ta,rh)
#else
  pure function find_eqvar(ta,rh)
#endif
    implicit none
    real(rkx) , dimension(4) :: find_eqvar
    real(rkx) , intent(in) :: ta , rh
    real(rkx) :: pa , phi , rf , rs , dtcdt , ts , ts_bar
    real(rkx) :: tf , ps , flux1 , flux2 , flux3 , m , m_bar

    pa    = rh*pvstar(ta)
    phi   = 0.84_rkx
    rs    = 0.0387_rkx
    dtcdt = 0.0_rkx

    m     = (pc-pa)/(zs(rs)+za)
    m_bar = (pc-pa)/(zs(rs)+za_bar)
    ts = solve1(ta,pa,rs,max(0.0_rkx,min(tc,ta)-rs*abs(m)), &
                max(tc,ta)+rs*abs(m),    err,maxiter)
    tf = solve2(ta,pa,rs,max(0.0_rkx,min(tc,ta)-rs*abs(m_bar)), &
                max(tc,ta)+rs*abs(m_bar),err,maxiter)
    ! c*dtc/dt when rf=zf=\inf
    flux1 = q-qv(ta,pa)-(1.0_rkx-phi)*(tc-ts)/rs
    ! c*dtc/dt when rf=zf=0
    flux2 = q-qv(ta,pa)-(1.0_rkx-phi)*(tc-ts)/rs - phi*(tc-tf)/rs
    if ( flux1 <= 0.0_rkx ) then ! region I
      phi = 1.0_rkx-(q-qv(ta,pa))*rs/(tc-ts)
      rf  = huge(0.0_rkx)
    else if ( flux2 <= 0.0_rkx ) then ! region II&III
      ts_bar = tc - (q-qv(ta,pa))*rs/phi + (1.0_rkx/phi -1.0_rkx)*(tc-ts)
      tf = solve3(ta,pa,rs,ts_bar,ta,ts_bar,err,maxiter)
      rf = ra_bar(tf,ta)*(ts_bar-tf)/(tf-ta)
    else ! region IV,V,VI
      rf = 0.0_rkx
      flux3 = q-qv(ta,pa)-(tc-ta)/ra_un(tc,ta)-(phi_salt*pvstar(tc)-pa)/za_un
      if ( flux3 <= 0.0_rkx ) then ! region IV,V
        ts = solve4(ta,pa,0.0_rkx,tc,err,maxiter)
        rs = (tc-ts)/(q-qv(ta,pa))
        ps = pc - (pc-pa)* zs(rs)/( zs(rs)+za_un)
        if ( ps > phi_salt * pvstar(ts) ) then  ! region V
          ts = solve5(ta,pa,0.0_rkx,tc,err,maxiter)
          rs = (tc-ts)/(q-qv(ta,pa))
        end if
      else ! region VI
        rs = 0.0_rkx
        dtcdt = (1.0_rkx/c)* flux3
      end if
    end if
    find_eqvar = (/ phi,rf,rs,dtcdt /)
  end function find_eqvar

#ifdef DEBUG
  real(rkx) function find_t(eqvar_indx, eqvar)
#else
  pure real(rkx) function find_t(eqvar_indx, eqvar)
#endif
    implicit none
    integer , intent(in) :: eqvar_indx
    real(rkx) , intent(in) :: eqvar
    real(rkx) :: t

    select case ( eqvar_indx )
      case ( eqvar_indx_1 )
        ! Region 1
        t = solvei(eqvar)
      case ( eqvar_indx_2 )
        ! Region 2 if pa0 > pvstar(t)
        ! Region 3 if pa0 < pvstar(t)
        t = solveii(eqvar)
      case ( eqvar_indx_3 )
        ! Region 4
        t = solveiii(eqvar)
      case ( eqvar_indx_4 )
        ! Region 5
        t = solveiii(eqvar)
      case default
        ! Region 6
        t = solveiv(eqvar)
    end select
    find_t = t
  end function find_t

#ifdef DEBUG
  real(rkx) function heatindex(ta,rh)
#else
  pure real(rkx) function heatindex(ta,rh)
#endif
    implicit none
    real(rkx), intent(in) :: ta , rh
    type(eqvar) :: initial

    if ( ta == 0.0_rkx ) then
      heatindex = 0.0_rkx
    else
      initial = initial_find_eqvar(ta, rh)
      select case ( initial%indx )
        case ( eqvar_indx_1 )
          heatindex = find_t(initial%indx, initial%var(1))
        case ( eqvar_indx_2 )
          heatindex = find_t(initial%indx, initial%var(2))
        case ( eqvar_indx_3, eqvar_indx_4 )
          heatindex = find_t(initial%indx, initial%var(3))
        case default
          heatindex = find_t(initial%indx, initial%var(4))
      end select
    end if
  end function heatindex

  !!!!!!!!!!!!!!!!! root solvers !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef DEBUG
  real(rkx) function solve1(ta,pa,rs,x1,x2,err,maxiter)
#else
  pure real(rkx) function solve1(ta,pa,rs,x1,x2,err,maxiter)
#endif
    implicit none
    real(rkx) , intent(in) :: ta , pa , rs , x1 , x2 , err
    integer , intent(in) :: maxiter
    integer :: iter
    real(rkx) :: a , b , c , fa , fb , fc

    a = x1
    b = x2
    fa = (a-ta)/ra(a,ta) + (pc-pa)/(zs(rs)+za) - (tc-a)/rs
    fb = (b-ta)/ra(b,ta) + (pc-pa)/(zs(rs)+za) - (tc-b)/rs
#ifdef DEBUG
    if ( fa*fb > 0.0_rkx ) then
      write(stderr,*) 'solve1 : ta, pa, rs : ', ta, pa, rs
      call fatal(__FILE__,__LINE__,'error interval, solve1')
    end if
#endif
    do iter = 1 , maxiter
      c = (a+b) * 0.5_rkx
      fc = (c-ta)/ra(c,ta) + (pc-pa)/(zs(rs)+za) - (tc-c)/rs
      if ( abs(b-a) < err ) exit
      if ( fb*fc > 0.0_rkx ) then
        b = c
        fb = fc
      else
        a = c
      end if
    end do
#ifdef DEBUG
    if ( iter == maxiter+1 ) write(stderr,*) "maxiter, solve1,"
#endif
    solve1 = c
  end function solve1

#ifdef DEBUG
  real(rkx) function solve2(ta,pa,rs,x1,x2,err,maxiter)
#else
  pure real(rkx) function solve2(ta,pa,rs,x1,x2,err,maxiter)
#endif
    implicit none
    real(rkx) , intent(in) :: ta , pa , rs , x1 , x2 , err
    integer , intent(in) :: maxiter
    integer :: iter
    real(rkx) :: a , b , c , fa , fb , fc

    a = x1
    b = x2
    fa = (a-ta)/ra_bar(a,ta) + (pc-pa)/(zs(rs)+za_bar) - (tc-a)/rs
    fb = (b-ta)/ra_bar(b,ta) + (pc-pa)/(zs(rs)+za_bar) - (tc-b)/rs
#ifdef DEBUG
    if ( fa*fb > 0.0_rkx ) then
      write(stderr,*) 'solve2'
      call fatal(__FILE__,__LINE__,'error interval, solve2')
    end if
#endif
    do iter = 1 , maxiter
      c = (a+b) * 0.5_rkx
      fc = (c-ta)/ra_bar(c,ta) + (pc-pa)/(zs(rs)+za_bar) - (tc-c)/rs
      if ( abs(b-a) < err ) exit
      if ( fb*fc > 0.0_rkx ) then
        b = c
        fb = fc
      else
        a = c
      end if
    end do
#ifdef DEBUG
    if ( iter == maxiter+1 ) write(stderr,*) "maxiter, solve2"
#endif
    solve2 = c
  end function solve2

#ifdef DEBUG
  real(rkx) function solve3(ta,pa,rs,ts_bar,x1,x2,err,maxiter)
#else
  pure real(rkx) function solve3(ta,pa,rs,ts_bar,x1,x2,err,maxiter)
#endif
    implicit none
    real(rkx) , intent(in) :: ta , pa , rs , ts_bar , x1 , x2 , err
    integer , intent(in) :: maxiter
    integer :: iter
    real(rkx) :: a , b , c , fa , fb , fc

    a = x1
    b = x2
    fa = (a-ta)/ra_bar(a,ta) + &
      (pc-pa)*(a-ta)/((a-ta)*(zs(rs)+za_bar)+r*ra_bar(a,ta)*(ts_bar-a)) - &
      (tc-ts_bar)/rs
    fb = (b-ta)/ra_bar(b,ta) + &
      (pc-pa)*(b-ta)/((b-ta)*(zs(rs)+za_bar)+r*ra_bar(b,ta)*(ts_bar-b)) - &
      (tc-ts_bar)/rs
#ifdef DEBUG
    if ( fa*fb > 0.0_rkx ) then
      write(stderr,*) 'solve3 : ta, pa, rs, ts_bar, a, b, fa, fb : ', &
        ta, pa, rs, ts_bar, a, b, fa, fb
      call fatal(__FILE__,__LINE__,'error interval, solve3')
    end if
#endif
    do iter = 1 , maxiter
      c = (a+b) * 0.5_rkx
      fc = (c-ta)/ra_bar(c,ta) + (pc-pa)*(c-ta)/((c-ta)*(zs(rs)+za_bar) + &
        r*ra_bar(c,ta)*(ts_bar-c)) - (tc-ts_bar)/rs
      if ( abs(b-a) < err ) exit
      if ( fb*fc > 0.0_rkx ) then
        b = c
        fb = fc
      else
        a = c
      end if
    end do
#ifdef DEBUG
    if ( iter == maxiter+1 ) write(stderr,*) "maxiter, solve3"
#endif
    solve3 = c
  end function solve3

#ifdef DEBUG
  real(rkx) function solve4(ta,pa,x1,x2,err,maxiter)
#else
  pure real(rkx) function solve4(ta,pa,x1,x2,err,maxiter)
#endif
    implicit none
    real(rkx) , intent(in) :: ta , pa , x1 , x2 , err
    integer , intent(in) :: maxiter
    integer :: iter
    real(rkx) :: a , b , c , fa , fb , fc

    a = x1
    b = x2
    fa = (a-ta)/ra_un(a,ta) + &
      (pc-pa)/(zs((tc-a)/(q-qv(ta,pa)))+za_un)-(q-qv(ta,pa))
    fb = (b-ta)/ra_un(b,ta) + &
      (pc-pa)/(zs((tc-b)/(q-qv(ta,pa)))+za_un)-(q-qv(ta,pa))
#ifdef DEBUG
    if ( fa*fb > 0.0_rkx ) then
      write(stderr,*) 'solve4'
      call fatal(__FILE__,__LINE__,'error interval, solve4')
    end if
#endif
    do iter = 1 , maxiter
      c = (a+b) * 0.5_rkx
      fc = (c-ta)/ra_un(c,ta) + &
        (pc-pa)/(zs((tc-c)/(q-qv(ta,pa)))+za_un)-(q-qv(ta,pa))
      if ( abs(b-a) < err ) exit
      if ( fb*fc > 0.0_rkx ) then
        b = c
        fb = fc
      else
        a = c
      end if
    end do
#ifdef DEBUG
    if ( iter == maxiter+1 ) write(stderr,*) "maxiter, solve4"
#endif
    solve4 = c
  end function solve4

#ifdef DEBUG
  real(rkx) function solve5(ta,pa,x1,x2,err,maxiter)
#else
  pure real(rkx) function solve5(ta,pa,x1,x2,err,maxiter)
#endif
    implicit none
    real(rkx) , intent(in) :: ta , pa , x1 , x2 , err
    integer , intent(in) :: maxiter
    integer :: iter
    real(rkx) :: a , b , c , fa , fb , fc

    a = x1
    b = x2
    fa = (a-ta)/ra_un(a,ta) + (phi_salt*pvstar(a)-pa)/za_un -(q-qv(ta,pa))
    fb = (b-ta)/ra_un(b,ta) + (phi_salt*pvstar(b)-pa)/za_un -(q-qv(ta,pa))
#ifdef DEBUG
    if ( fa*fb > 0.0_rkx ) then
      write(stderr,*) 'solve5'
      call fatal(__FILE__,__LINE__,'error interval, solve5')
    end if
#endif
    do iter = 1 , maxiter
      c = (a+b) * 0.5_rkx
      fc = (c-ta)/ra_un(c,ta) + (phi_salt*pvstar(c)-pa)/za_un -(q-qv(ta,pa))
      if ( abs(b-a) < err ) exit
      if ( fb*fc > 0.0_rkx ) then
        b = c
        fb = fc
      else
        a = c
      end if
    end do
#ifdef DEBUG
    if ( iter == maxiter+1 ) write(stderr,*) "maxiter, solve5"
#endif
    solve5 = c
  end function solve5

#ifdef DEBUG
  real(rkx) function solvei(eqvar)
#else
  pure real(rkx) function solvei(eqvar)
#endif
    implicit none
    real(rkx) , intent(in) :: eqvar
    integer :: iter
    real(rkx) :: a , b , c , fa , fb , fc
    real(rkx) , dimension(4) :: tmp

    a = 0.0_rkx
    b = 240.0_rkx
    tmp = find_eqvar(a,1.0_rkx)
    fa = tmp(1) - eqvar
    tmp = find_eqvar(b,1.0_rkx)
    fb = tmp(1) - eqvar
#ifdef DEBUG
    if ( fa*fb > 0.0_rkx ) then
      write(stderr,*) 'solvei'
      call fatal(__FILE__,__LINE__,'error interval, solvei')
    end if
#endif
    do iter = 1 , maxiter
      c = (a+b) * 0.5_rkx
      tmp = find_eqvar(c,1.0_rkx)
      fc = tmp(1) - eqvar
      if ( abs(b-a) < errt .or. abs(fc) < tiny(fc) ) exit
      if ( fb*fc > 0.0_rkx ) then
        b = c
        fb = fc
      else
        a = c
      end if
    end do
#ifdef DEBUG
    if ( iter == maxiter+1 ) write(stderr,*) "maxiter, solvei"
#endif
    solvei = c
  end function solvei

#ifdef DEBUG
  real(rkx) function solveii(eqvar)
#else
  pure real(rkx) function solveii(eqvar)
#endif
    implicit none
    real(rkx) , intent(in) :: eqvar
    integer :: iter
    real(rkx) :: a , b , c , fa , fb , fc
    real(rkx) , dimension(4) :: tmp

    a = 230.0_rkx
    b = 300.0_rkx
    tmp = find_eqvar(a,min(1.0_rkx,pa0/pvstar(a)))
    fa = tmp(2) - eqvar
    tmp = find_eqvar(b,min(1.0_rkx,pa0/pvstar(b)))
    fb = tmp(2) - eqvar
#ifdef DEBUG
    if ( fa*fb > 0.0_rkx ) then
      write(stderr,*) 'solveii : a, b, fa, fb, eqvar : ', &
        a, b, fa, fb, eqvar
      call fatal(__FILE__,__LINE__,'error interval, solveii')
    end if
#endif
    do iter = 1 , maxiter
      c = (a+b) * 0.5_rkx
      tmp = find_eqvar(c,min(1.0_rkx,pa0/pvstar(c)))
      fc = tmp(2) - eqvar
      if ( abs(b-a) < errt .or. abs(fc) < tiny(fc) ) exit
      if ( fb*fc > 0.0_rkx ) then
        b = c
        fb = fc
      else
        a = c
      end if
    end do
#ifdef DEBUG
    if ( iter == maxiter+1 ) write(stderr,*) "maxiter, solveii"
#endif
    solveii = c
  end function solveii

#ifdef DEBUG
  real(rkx) function solveiii(eqvar)
#else
  pure real(rkx) function solveiii(eqvar)
#endif
    implicit none
    real(rkx) , intent(in) :: eqvar
    integer :: iter
    real(rkx) :: a , b , c , fa , fb , fc
    real(rkx) , dimension(4) :: tmp

    a = 295.0_rkx
    b = 350.0_rkx
    tmp = find_eqvar(a,pa0/pvstar(a))
    fa = tmp(3) - eqvar
    tmp = find_eqvar(b,pa0/pvstar(b))
    fb = tmp(3) - eqvar
#ifdef DEBUG
    if ( fa*fb > 0.0_rkx ) then
      write(stderr,*) 'solveiii'
      call fatal(__FILE__,__LINE__,'error interval, solveiii')
    end if
#endif
    do iter = 1 , maxiter
      c = (a+b) * 0.5_rkx
      tmp = find_eqvar(c,pa0/pvstar(c))
      fc = tmp(3) - eqvar
      if ( abs(b-a) < errt .or. abs(fc) < tiny(fc) ) exit
      if ( fb*fc > 0.0_rkx ) then
        b = c
        fb = fc
      else
        a = c
      end if
    end do
#ifdef DEBUG
    if ( iter == maxiter+1 ) write(stderr,*) "maxiter, solveiii"
#endif
    solveiii = c
  end function solveiii

#ifdef DEBUG
  real(rkx) function solveiv(eqvar)
#else
  pure real(rkx) function solveiv(eqvar)
#endif
    implicit none
    real(rkx) , intent(in) :: eqvar
    integer :: iter
    real(rkx) :: a , b , c , fa , fb , fc
    real(rkx) , dimension(4) :: tmp

    a = 340.0_rkx
    b = 800.0_rkx
    tmp = find_eqvar(a,pa0/pvstar(a))
    fa = tmp(4) - eqvar
    tmp = find_eqvar(b,pa0/pvstar(b))
    fb = tmp(4) - eqvar
#ifdef DEBUG
    if ( fa*fb > 0.0_rkx ) then
      write(stderr,*) 'solveiv : eqvar : ', eqvar
      call fatal(__FILE__,__LINE__,'error interval, solveiv')
    end if
#endif
    do iter = 1 , maxiter
      c = (a+b) * 0.5_rkx
      tmp = find_eqvar(c,pa0/pvstar(c))
      fc = tmp(4) - eqvar
      if ( abs(b-a) < errt .or. abs(fc) < tiny(fc) ) exit
      if ( fb*fc > 0.0_rkx ) then
        b = c
        fb = fc
      else
        a = c
      end if
    end do
#ifdef DEBUG
    if ( iter == maxiter+1 ) write(stderr,*) "maxiter, solveiv"
#endif
    solveiv = c
  end function solveiv

end module mod_heatindex

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
