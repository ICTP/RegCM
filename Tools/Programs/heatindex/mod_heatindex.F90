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

  use iso_fortran_env

  implicit none

  private

  integer , parameter :: wp = real64

  ! w/m^2/k^4 , boltzmann constant
  real(wp) , parameter :: sigma    = 5.67e-8_wp
  ! k         , vapor temperature at triple point
  real(wp) , parameter :: ttrip    = 273.16_wp
  ! pa        , vapor pressure at triple point
  real(wp) , parameter :: ptrip    = 611.65_wp
  real(wp) , parameter :: e0v      = 2.3740e6_wp  ! J/kg
  real(wp) , parameter :: e0s      = 0.3337e6_wp  ! J/kg
  real(wp) , parameter :: rgasa    = 287.04_wp    ! J/kg/K
  real(wp) , parameter :: rgasv    = 461.0_wp     ! J/kg/K
  real(wp) , parameter :: cva      = 719.0_wp     ! J/kg/K
  real(wp) , parameter :: cvv      = 1418.0_wp    ! J/kg/K
  real(wp) , parameter :: cvl      = 4119.0_wp    ! J/kg/K
  real(wp) , parameter :: cvs      = 1861.0_wp    ! J/kg/K
  real(wp) , parameter :: cpa      = cva + rgasa
  real(wp) , parameter :: cpv      = cvv + rgasv
  !           , emissivity of surface
  real(wp) , parameter :: reps  = 0.97_wp
  ! kg        , mass of average us adults, FRYAR2018
  real(wp) , parameter :: m        = 83.6_wp
  ! m         , height of average us adults, FRYAR2018
  real(wp) , parameter :: h        = 1.69_wp
  ! m^2       , dubois formula, parson2014
  real(wp) , parameter :: a        = 0.202_wp*(m**0.425_wp)*(h**0.725_wp)
  ! J/kg/K    , specific heat capacity of core, Gagge1972
  real(wp) , parameter :: cpc      = 3492.0_wp
  !           , heat capacity of core
  real(wp) , parameter :: c        = m*cpc/a
  ! Pa/K      , zf/rf
  real(wp) , parameter :: r        = 124.0_wp
  ! W/m^2     , metabolic rate per skin area
  real(wp) , parameter :: q        = 180.0_wp
  !           , vapor saturation pressure level of saline solution
  real(wp) , parameter :: phi_salt = 0.9_wp
  ! K         , core temeprature
  real(wp) , parameter :: tc        = 310.0_wp
  !           , core vapor pressure, phi_salt*pvstar(tc)
  real(wp) , parameter :: pc        = 5612.299124343398_wp
  ! Pa        , atmospheric pressure
  real(wp) , parameter :: p         = 1.013e5_wp
  ! kg/J      , "inhaled mass" / "metabolic rate"
  real(wp) , parameter :: eta       = 1.43e-6_wp
  ! Pa        , reference air vapor pressure in regions III, IV, V, VI,
  !             chosen by Steadman
  real(wp) , parameter :: pa0       = 1.6e3_wp
  ! Pa m^2/W  , mass transfer resistance through air, exposed part of skin
  real(wp) , parameter :: za        = 60.6_wp/17.4_wp
  ! Pa m^2/W  , mass transfer resistance through air, clothed part of skin
  real(wp) , parameter :: za_bar    = 60.6_wp/11.6_wp
  ! Pa m^2/W  , mass transfer resistance through air, when being naked
  real(wp) , parameter :: za_un     = 60.6_wp/12.3_wp
  !           , tolorence of the root solver for heatindex
  real(wp) , parameter :: errt      = 1.0e-8_wp
  !           , tolorence of the root solver
  real(wp) , parameter :: err       = 1.0e-8_wp

  ! maximum number of iteration of the root solver
  integer, parameter :: maxiter = 100

  integer , parameter :: region_1 = 1
  integer , parameter :: region_2 = 2
  integer , parameter :: region_3 = 3
  integer , parameter :: region_4 = 4
  integer , parameter :: region_5 = 5
  integer , parameter :: region_6 = 6

  integer , parameter :: eqvar_indx_1 = 1 ! phi
  integer , parameter :: eqvar_indx_2 = 2 ! rf
  integer , parameter :: eqvar_indx_3 = 3 ! rs
  integer , parameter :: eqvar_indx_4 = 4 ! rs*
  integer , parameter :: eqvar_indx_5 = 5 ! dtcdt

  public  :: heatindex

  contains

  pure real(wp) function pvstar(t)
    implicit none
    real(wp) , intent(in) :: t
    if ( t <= 0.0_wp ) then
      pvstar = 0.0_wp
    else if ( t < ttrip ) then
      pvstar = ptrip * (t/ttrip)**((cpv-cvs)/rgasv) * &
        exp( (e0v + e0s -(cvv-cvs)*ttrip)/rgasv * (1.0_wp/ttrip - 1.0_wp/t) )
    else
      pvstar = ptrip * (t/ttrip)**((cpv-cvl)/rgasv) * &
        exp( (e0v       -(cvv-cvl)*ttrip)/rgasv * (1.0_wp/ttrip - 1.0_wp/t) )
    end if
  end function pvstar

  pure real(wp) function le(t)
    implicit none
    real(wp) , intent(in):: t
    le = (e0v + (cvv-cvl)*(t-ttrip) + rgasv*t)
  end function le

  pure real(wp) function qv(ta,pa) ! respiratory heat loss, w/m^2
    implicit none
    real(wp) , intent(in) :: ta , pa
    qv = eta * q *(cpa*(tc-ta) + le(tc)*rgasa/(p*rgasv) * ( pc-pa ) )
  end function qv

  ! mass transfer resistance through skin, pa m^2/w
  pure real(wp) function zs(rs)
    implicit none
    real(wp) , intent(in) :: rs
    if ( rs == 0.0387_wp ) then
      zs = 52.1_wp
    else
      zs = 6.0e8_wp * rs**5
    end if
  end function zs

  ! heat transfer resistance through air, exposed part of skin, k m^2/w
  pure real(wp) function ra(ts,ta)
    implicit none
    real(wp) , intent(in) :: ts , ta
    real(wp) :: hc , phi_rad , hr
    hc      = 17.4_wp
    phi_rad = 0.85_wp
    hr      = reps * phi_rad * sigma* (ts**2 + ta**2)*(ts + ta)
    ra      = 1.0_wp/(hc+hr)
  end function ra

  ! heat transfer resistance through air, clothed part of skin, k m^2/w
  pure real(wp) function ra_bar(tf,ta)
    implicit none
    real(wp), intent(in) :: tf , ta
    real(wp) :: hc , phi_rad, hr
    hc      = 11.6_wp
    phi_rad = 0.79_wp
    hr      = reps * phi_rad * sigma* (tf**2 + ta**2)*(tf + ta)
    ra_bar  = 1.0_wp/(hc+hr)
  end function ra_bar

  ! heat transfer resistance through air, when being naked, k m^2/w
  pure real(wp) function ra_un(ts,ta)
    implicit none
    real(wp) , intent(in) :: ts , ta
    real(wp) :: hc , phi_rad , hr
    hc      = 12.3_wp
    phi_rad = 0.80_wp
    hr      = reps * phi_rad * sigma* (ts**2 + ta**2)*(ts + ta)
    ra_un   = 1.0_wp/(hc+hr)
  end function ra_un

  function initial_find_eqvar(ta,rh,eqvar_indx)
    implicit none
    real(wp) , dimension(4) :: initial_find_eqvar
    real(wp) , intent(in) :: ta , rh
    integer , intent(out) :: eqvar_indx
    real(wp) :: pa , phi , rf , rs , dtcdt , ts , ts_bar
    real(wp) :: tf , ps , flux1 , flux2 , flux3 , m , m_bar

    pa    = rh*pvstar(ta)
    phi   = 0.84_wp
    rs    = 0.0387_wp
    dtcdt = 0.0_wp

    m     = (pc-pa)/(zs(rs)+za)
    m_bar = (pc-pa)/(zs(rs)+za_bar)
    ts = solve1(ta,pa,rs,max(0.0_wp,min(tc,ta)-rs*abs(m)), &
                max(tc,ta)+rs*abs(m),    err,maxiter)
    tf = solve2(ta,pa,rs,max(0.0_wp,min(tc,ta)-rs*abs(m_bar)), &
                max(tc,ta)+rs*abs(m_bar),err,maxiter)
    ! c*dtc/dt when rf=zf=\inf
    flux1 = q-qv(ta,pa)-(1.0_wp-phi)*(tc-ts)/rs
    ! c*dtc/dt when rf=zf=0
    flux2 = q-qv(ta,pa)-(1.0_wp-phi)*(tc-ts)/rs - phi*(tc-tf)/rs
    if ( flux1 <= 0.0_wp ) then ! region I
      eqvar_indx = eqvar_indx_1
      phi = 1.0_wp-(q-qv(ta,pa))*rs/(tc-ts)
      rf  = huge(0.0_wp)
    else if ( flux2 <= 0.0_wp ) then ! region II&III
      eqvar_indx = eqvar_indx_2
      ts_bar = tc - (q-qv(ta,pa))*rs/phi + (1.0_wp/phi -1.0_wp)*(tc-ts)
      tf = solve3(ta,pa,rs,ts_bar,ta,ts_bar,err,maxiter)
      rf = ra_bar(tf,ta)*(ts_bar-tf)/(tf-ta)
    else ! region IV,V,VI
      rf = 0.0_wp
      flux3 = q-qv(ta,pa)-(tc-ta)/ra_un(tc,ta)-(phi_salt*pvstar(tc)-pa)/za_un
      if ( flux3 <= 0.0_wp ) then ! region IV,V
        ts = solve4(ta,pa,0.0_wp,tc,err,maxiter)
        rs = (tc-ts)/(q-qv(ta,pa))
        eqvar_indx = eqvar_indx_3
        ps = pc - (pc-pa)* zs(rs)/( zs(rs)+za_un)
        if ( ps > phi_salt * pvstar(ts) ) then  ! region V
          ts = solve5(ta,pa,0.0_wp,tc,err,maxiter)
          rs = (tc-ts)/(q-qv(ta,pa))
          eqvar_indx = eqvar_indx_3
        end if
      else ! region VI
        rs = 0.0_wp
        eqvar_indx = eqvar_indx_4
        dtcdt = (1.0_wp/c)* flux3
      end if
    end if
    initial_find_eqvar = (/ phi,rf,rs,dtcdt /)
  end function initial_find_eqvar

  function find_eqvar(ta,rh)
    implicit none
    real(wp) , dimension(4) :: find_eqvar
    real(wp) , intent(in) :: ta , rh
    real(wp) :: pa , phi , rf , rs , dtcdt , ts , ts_bar
    real(wp) :: tf , ps , flux1 , flux2 , flux3 , m , m_bar

    pa    = rh*pvstar(ta)
    phi   = 0.84_wp
    rs    = 0.0387_wp
    dtcdt = 0.0_wp

    m     = (pc-pa)/(zs(rs)+za)
    m_bar = (pc-pa)/(zs(rs)+za_bar)
    ts = solve1(ta,pa,rs,max(0.0_wp,min(tc,ta)-rs*abs(m)), &
                max(tc,ta)+rs*abs(m),    err,maxiter)
    tf = solve2(ta,pa,rs,max(0.0_wp,min(tc,ta)-rs*abs(m_bar)), &
                max(tc,ta)+rs*abs(m_bar),err,maxiter)
    ! c*dtc/dt when rf=zf=\inf
    flux1 = q-qv(ta,pa)-(1.0_wp-phi)*(tc-ts)/rs
    ! c*dtc/dt when rf=zf=0
    flux2 = q-qv(ta,pa)-(1.0_wp-phi)*(tc-ts)/rs - phi*(tc-tf)/rs
    if ( flux1 <= 0.0_wp ) then ! region I
      phi = 1.0_wp-(q-qv(ta,pa))*rs/(tc-ts)
      rf  = huge(0.0_wp)
    else if ( flux2 <= 0.0_wp ) then ! region II&III
      ts_bar = tc - (q-qv(ta,pa))*rs/phi + (1.0_wp/phi -1.0_wp)*(tc-ts)
      tf = solve3(ta,pa,rs,ts_bar,ta,ts_bar,err,maxiter)
      rf = ra_bar(tf,ta)*(ts_bar-tf)/(tf-ta)
    else ! region IV,V,VI
      rf = 0.0_wp
      flux3 = q-qv(ta,pa)-(tc-ta)/ra_un(tc,ta)-(phi_salt*pvstar(tc)-pa)/za_un
      if ( flux3 <= 0.0_wp ) then ! region IV,V
        ts = solve4(ta,pa,0.0_wp,tc,err,maxiter)
        rs = (tc-ts)/(q-qv(ta,pa))
        ps = pc - (pc-pa)* zs(rs)/( zs(rs)+za_un)
        if ( ps > phi_salt * pvstar(ts) ) then  ! region V
          ts = solve5(ta,pa,0.0_wp,tc,err,maxiter)
          rs = (tc-ts)/(q-qv(ta,pa))
        end if
      else ! region VI
        rs = 0.0_wp
        dtcdt = (1.0_wp/c)* flux3
      end if
    end if
    find_eqvar = (/ phi,rf,rs,dtcdt /)
  end function find_eqvar

  real(wp) function find_t(eqvar_indx, eqvar)
    implicit none
    integer , intent(in) :: eqvar_indx
    real(wp) , intent(in) :: eqvar
    real(wp) :: t

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

  real(wp) function heatindex(ta,rh)
    implicit none
    real(wp), intent(in) :: ta , rh
    real(wp) :: eqvars(4)
    integer :: eqvar_indx

    if ( ta == 0.0_wp ) then
      heatindex = 0.0_wp
    else
      eqvars = initial_find_eqvar(ta, rh, eqvar_indx)
      select case ( eqvar_indx )
        case ( eqvar_indx_1 )
          heatindex = find_t(eqvar_indx, eqvars(1))
        case ( eqvar_indx_2 )
          heatindex = find_t(eqvar_indx, eqvars(2))
        case ( eqvar_indx_3, eqvar_indx_4 )
          heatindex = find_t(eqvar_indx, eqvars(3))
        case default
          heatindex = find_t(eqvar_indx, eqvars(4))
      end select
    end if
  end function heatindex

  !!!!!!!!!!!!!!!!! root solvers !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(wp) function solve1(ta,pa,rs,x1,x2,err,maxiter)
    implicit none
    real(wp) , intent(in) :: ta , pa , rs , x1 , x2 , err
    integer , intent(in) :: maxiter
    integer :: iter
    real(wp) :: a , b , c , fa , fb , fc

    a = x1
    b = x2
    fa = (a-ta)/ra(a,ta) + (pc-pa)/(zs(rs)+za) - (tc-a)/rs
    fb = (b-ta)/ra(b,ta) + (pc-pa)/(zs(rs)+za) - (tc-b)/rs
    if ( fa*fb > 0.0_wp ) then
      write(error_unit,*) 'solve1 : ta, pa, rs : ', ta, pa, rs
      stop 'SOLVER_ERROR'
    end if
    do iter = 1 , maxiter
      c = (a+b) * 0.5_wp
      fc = (c-ta)/ra(c,ta) + (pc-pa)/(zs(rs)+za) - (tc-c)/rs
      if ( abs(b-a) < err ) exit
      if ( fb*fc > 0.0_wp ) then
        b = c
        fb = fc
      else
        a = c
        fa = fc
      end if
    end do
    if ( iter == maxiter+1 ) write(error_unit,*) "maxiter, solve1"
    solve1 = c
  end function solve1

  real(wp) function solve2(ta,pa,rs,x1,x2,err,maxiter)
    implicit none
    real(wp) , intent(in) :: ta , pa , rs , x1 , x2 , err
    integer , intent(in) :: maxiter
    integer :: iter
    real(wp) :: a , b , c , fa , fb , fc

    a = x1
    b = x2
    fa = (a-ta)/ra_bar(a,ta) + (pc-pa)/(zs(rs)+za_bar) - (tc-a)/rs
    fb = (b-ta)/ra_bar(b,ta) + (pc-pa)/(zs(rs)+za_bar) - (tc-b)/rs
    if ( fa*fb > 0.0_wp ) then
      write(error_unit,*) 'solve2'
      stop 'SOLVER_ERROR'
    end if
    do iter = 1 , maxiter
      c = (a+b) * 0.5_wp
      fc = (c-ta)/ra_bar(c,ta) + (pc-pa)/(zs(rs)+za_bar) - (tc-c)/rs
      if ( abs(b-a) < err ) exit
      if ( fb*fc > 0.0_wp ) then
        b = c
        fb = fc
      else
        a = c
        fa = fc
      end if
    end do
    if ( iter == maxiter+1 ) write(error_unit,*) "maxiter, solve2"
    solve2 = c
  end function solve2

  real(wp) function solve3(ta,pa,rs,ts_bar,x1,x2,err,maxiter)
    implicit none
    real(wp) , intent(in) :: ta , pa , rs , ts_bar , x1 , x2 , err
    integer , intent(in) :: maxiter
    integer :: iter
    real(wp) :: a , b , c , fa , fb , fc

    a = x1
    b = x2
    fa = (a-ta)/ra_bar(a,ta) + &
      (pc-pa)*(a-ta)/((a-ta)*(zs(rs)+za_bar)+r*ra_bar(a,ta)*(ts_bar-a)) - &
      (tc-ts_bar)/rs
    fb = (b-ta)/ra_bar(b,ta) + &
      (pc-pa)*(b-ta)/((b-ta)*(zs(rs)+za_bar)+r*ra_bar(b,ta)*(ts_bar-b)) - &
      (tc-ts_bar)/rs
    if ( fa*fb > 0.0_wp ) then
      write(error_unit,*) 'solve3 : ta, pa, rs, ts_bar, a, b, fa, fb : ', &
        ta, pa, rs, ts_bar, a, b, fa, fb
      stop 'SOLVER_ERROR'
    end if
    do iter = 1 , maxiter
      c = (a+b) * 0.5_wp
      fc = (c-ta)/ra_bar(c,ta) + (pc-pa)*(c-ta)/((c-ta)*(zs(rs)+za_bar) + &
        r*ra_bar(c,ta)*(ts_bar-c)) - (tc-ts_bar)/rs
      if ( abs(b-a) < err ) exit
      if ( fb*fc > 0.0_wp ) then
        b = c
        fb = fc
      else
        a = c
        fa = fc
      end if
    end do
    if ( iter == maxiter+1 ) write(error_unit,*) "maxiter, solve3"
    solve3 = c
  end function solve3

  real(wp) function solve4(ta,pa,x1,x2,err,maxiter)
    implicit none
    real(wp) , intent(in) :: ta , pa , x1 , x2 , err
    integer , intent(in) :: maxiter
    integer :: iter
    real(wp) :: a , b , c , fa , fb , fc

    a = x1
    b = x2
    fa = (a-ta)/ra_un(a,ta) + &
      (pc-pa)/(zs((tc-a)/(q-qv(ta,pa)))+za_un)-(q-qv(ta,pa))
    fb = (b-ta)/ra_un(b,ta) + &
      (pc-pa)/(zs((tc-b)/(q-qv(ta,pa)))+za_un)-(q-qv(ta,pa))
    if ( fa*fb > 0.0_wp ) then
      write(error_unit,*) 'solve4'
      stop 'SOLVER_ERROR'
    end if
    do iter = 1 , maxiter
      c = (a+b) * 0.5_wp
      fc = (c-ta)/ra_un(c,ta) + &
        (pc-pa)/(zs((tc-c)/(q-qv(ta,pa)))+za_un)-(q-qv(ta,pa))
      if ( abs(b-a) < err ) exit
      if ( fb*fc > 0.0_wp ) then
        b = c
        fb = fc
      else
        a = c
        fa = fc
      end if
    end do
    if ( iter == maxiter+1 ) write(error_unit,*) "maxiter, solve4"
    solve4 = c
  end function solve4

  real(wp) function solve5(ta,pa,x1,x2,err,maxiter)
    implicit none
    real(wp) , intent(in) :: ta , pa , x1 , x2 , err
    integer , intent(in) :: maxiter
    integer :: iter
    real(wp) :: a , b , c , fa , fb , fc

    a = x1
    b = x2
    fa = (a-ta)/ra_un(a,ta) + (phi_salt*pvstar(a)-pa)/za_un -(q-qv(ta,pa))
    fb = (b-ta)/ra_un(b,ta) + (phi_salt*pvstar(b)-pa)/za_un -(q-qv(ta,pa))
    if ( fa*fb > 0.0_wp ) then
      write(error_unit,*) 'solve5'
      stop 'SOLVER_ERROR'
    end if
    do iter = 1 , maxiter
      c = (a+b) * 0.5_wp
      fc = (c-ta)/ra_un(c,ta) + (phi_salt*pvstar(c)-pa)/za_un -(q-qv(ta,pa))
      if ( abs(b-a) < err ) exit
      if ( fb*fc > 0.0_wp ) then
        b = c
        fb = fc
      else
        a = c
        fa = fc
      end if
    end do
    if ( iter == maxiter+1 ) write(error_unit,*) "maxiter, solve5"
    solve5 = c
  end function solve5

  real(wp) function solvei(eqvar)
    implicit none
    real(wp) , intent(in) :: eqvar
    integer :: iter
    real(wp) :: a , b , c , fa , fb , fc
    real(wp) , dimension(4) :: tmp

    a = 0.0_wp
    b = 240.0_wp
    tmp = find_eqvar(a,1.0_wp)
    fa = tmp(1) - eqvar
    tmp = find_eqvar(b,1.0_wp)
    fb = tmp(1) - eqvar
    if ( fa*fb > 0.0_wp ) then
      write(error_unit,*) 'solvei'
      stop 'SOLVER_ERROR'
    end if
    do iter = 1 , maxiter
      c = (a+b) * 0.5_wp
      tmp = find_eqvar(c,1.0_wp)
      fc = tmp(1) - eqvar
      if ( abs(b-a) < errt ) exit
      if ( fb*fc > 0.0_wp ) then
        b = c
        fb = fc
      else
        a = c
        fa = fc
      end if
    end do
    if ( iter == maxiter+1 ) write(error_unit,*) "maxiter, solvei"
    solvei = c
  end function solvei

  real(wp) function solveii(eqvar)
    implicit none
    real(wp) , intent(in) :: eqvar
    integer :: iter
    real(wp) :: a , b , c , fa , fb , fc
    real(wp) , dimension(4) :: tmp
    integer :: region

    a = 230.0_wp
    b = 300.0_wp
    tmp = find_eqvar(a,min(1.0_wp,pa0/pvstar(a)))
    fa = tmp(2) - eqvar
    tmp = find_eqvar(b,min(1.0_wp,pa0/pvstar(b)))
    fb = tmp(2) - eqvar
    if ( fa*fb > 0.0_wp ) then
      write(error_unit,*) 'solveii : a, b, fa, fb, eqvar : ', &
        a, b, fa, fb, eqvar
      stop 'SOLVER_ERROR'
    end if
    do iter = 1 , maxiter
      c = (a+b) * 0.5_wp
      tmp = find_eqvar(c,min(1.0_wp,pa0/pvstar(c)))
      fc = tmp(2) - eqvar
      if ( abs(b-a) < errt ) exit
      if ( fb*fc > 0.0_wp ) then
        b = c
        fb = fc
      else
        a = c
        fa = fc
      end if
    end do
    if ( iter == maxiter+1 ) write(error_unit,*) "maxiter, solveii"
    solveii = c
  end function solveii

  real(wp) function solveiii(eqvar)
    implicit none
    real(wp) , intent(in) :: eqvar
    integer :: iter
    real(wp) :: a , b , c , fa , fb , fc
    real(wp) , dimension(4) :: tmp
    integer :: region

    a = 295.0_wp
    b = 350.0_wp
    tmp = find_eqvar(a,pa0/pvstar(a))
    fa = tmp(3) - eqvar
    tmp = find_eqvar(b,pa0/pvstar(b))
    fb = tmp(3) - eqvar
    if ( fa*fb > 0.0_wp ) then
      write(error_unit,*) 'solveiii'
      stop 'SOLVER_ERROR'
    end if
    do iter = 1 , maxiter
      c = (a+b) * 0.5_wp
      tmp = find_eqvar(c,pa0/pvstar(c))
      fc = tmp(3) - eqvar
      if ( abs(b-a) < errt ) exit
      if ( fb*fc > 0.0_wp ) then
        b = c
        fb = fc
      else
        a = c
        fa = fc
      end if
    end do
    if ( iter == maxiter+1 ) write(error_unit,*) "maxiter, solveiii"
    solveiii = c
  end function solveiii

  real(wp) function solveiv(eqvar)
    implicit none
    real(wp) , intent(in) :: eqvar
    integer :: iter
    real(wp) :: a , b , c , fa , fb , fc
    real(wp) , dimension(4) :: tmp
    integer :: region

    a = 340.0_wp
    b = 800.0_wp
    tmp = find_eqvar(a,pa0/pvstar(a))
    fa = tmp(4) - eqvar
    tmp = find_eqvar(b,pa0/pvstar(b))
    fb = tmp(4) - eqvar
    if ( fa*fb > 0.0_wp ) then
      write(error_unit,*) 'solveiv : eqvar : ', eqvar
      stop 'SOLVER_ERROR'
    end if
    do iter = 1 , maxiter
      c = (a+b) * 0.5_wp
      tmp = find_eqvar(c,pa0/pvstar(c))
      fc = tmp(4) - eqvar
      if ( abs(b-a) < errt ) exit
      if ( fb*fc > 0.0_wp ) then
        b = c
        fb = fc
      else
        a = c
        fa = fc
      end if
    end do
    if ( iter == maxiter+1 ) write(error_unit,*) "maxiter, solveiv"
    solveiv = c
  end function solveiv

end module mod_heatindex

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
