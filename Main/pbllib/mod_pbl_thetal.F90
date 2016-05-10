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

module mod_pbl_thetal

  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_mpmessage
  use mod_stdio

  implicit none

  private

  real(rkx) :: myp , mythetal , myqt , myexner

  real(rkx) , parameter :: delta = 1.0e-8_rkx
  real(rkx) , parameter :: mindt = 1.0e-3_rkx
  real(rkx) , parameter :: mindq = 1.0e-3_rkx

  public :: solve_for_t

  contains

#include <pfesat.inc>
#include <pfqsat.inc>

  real(rkx) function zerofunc(t,bderiv)
    implicit none
    real(rkx) , intent(in) :: t
    logical , intent(in) :: bderiv
    real(rkx) :: qc , qv , es

    es = pfesat(t)
    qv = max(ep2/(myp/es-d_one),d_zero)
    if ( myqt > qv ) then
      qc = myqt - qv
    else
      qv = myqt
      qc = d_zero
    end if
    zerofunc = myexner*mythetal + wlhvocp*qc - t
    ! If necessary, find the minimum of the zerofuncion, which should be
    ! the place where it equals zero if no solution is available otherwise
    if ( bderiv ) then
      zerofunc = -(cpd/rwat)*(myp/es)*(wlhvocp*qv/t)**2 - d_one
    end if
  end function zerofunc

  ! Determine temperature from the liquid water potential temperature,
  ! total water, and pressure.

  real(rkx) function solve_for_t(thetal,qt,p,tprev,qtprev,   &
                                 qcprev,thlprev,imax,imethod, &
                                 outqv,outqc)
    implicit none
    real(rkx) , intent(in) :: thetal , qt , p , tprev , qtprev , &
                              qcprev , thlprev
    integer(ik4) , intent(in) :: imethod , imax
    real(rkx) , intent(out) :: outqc , outqv
    real(rkx) :: a , b , c , d , s , dum
    real(rkx) :: fa , fb , fc , fs
    real(rkx) :: smb , bmc , cmd
    real(rkx) :: temps , templ , rvls
    real(rkx) :: dthetal , dqt , qvprev , deltat , es , qcthresh
    real(rkx) :: tempqv , tempqc , tempt , tempes , ddq
    logical :: mflag
    integer(ik4) :: iteration , itqt
    integer(ik4) , parameter :: itmax = 6

    ! set some interal module variables
    mythetal = thetal
    myqt = qt
    myp = p
    myexner = (myp/1.0e5_rkx)**rovcp
    templ = mythetal*myexner

    ! Use the Brent 1973 method to determine t from liquid water pot.
    ! temp., pressure, and total water
    ! (see http://en.wikipedia.org/wiki/brent's_method)

    if ( imethod == 1 ) then

      ! set the bounds on t based on the change in thetal and qw:
      ! the bounds are: all of the change in qw goes in to changing
      ! cloud water at one end, and none of the change in qw goes in to it
      ! in the other
      a = tprev + myexner*(mythetal - thlprev)
      b = a + wlhvocp*(myqt-qtprev)

      ! if a and b are essentially the same, then the total water did not
      ! change, and so t is given only by the change in liquid water
      ! potential temperature; return with this as the solution
      if ( abs(a-b) <= delta ) then
        solve_for_t = b
        call getqvqc(solve_for_t,es,outqv,outqc)
        return
      end if

      ! make sure that the dynamic boundaries above don't conflict with the
      ! thermodynamic boundaries
      if ( a > b ) call swap(a,b)
      dum = myexner*mythetal
      a = max(a, dum)
      b = min(b, dum + wlhvocp*myqt)
      ! run these through the solution-finder
      fa = zerofunc(a,.false.)
      fb = zerofunc(b,.false.)
      fs = zerofunc(tprev,.false.)

      ! check if we already have the solution:
      if ( fa == 0 ) then
        b = a
        fb = fa
      end if
      ! check if the previous timestep works as a solution
      if ( abs(fs) < dlowval ) then
        b = tprev
        fb = fs
      end if
      if ( abs(fb) < dlowval ) then
        solve_for_t = b
        call getqvqc(solve_for_t,es,outqv,outqc)
        return
      end if

      ! if the initial guesses are out of range, return with a bogus
      ! temperature to flag a failure
      if ( fa*fb >= d_zero ) then
        solve_for_t = dum
        solve_for_t = tprev
        call getqvqc(solve_for_t,es,outqv,outqc)
      end if

      ! put a and b in the proper order
      if ( abs(fa) < abs(fb) ) then
        call swap(a,b)
        call swap(fa,fb)
      end if

      ! set c and f(c)
      c = a
      fc = fa
      ! set the bisection flag
      mflag = .true.
      d = d_zero

      !*************** the main iteration loop *******************

      do iteration = 1 , imax
        if ( fa /= fc .and. fb /= fc ) then
          ! inverse quadratic interpolation
          s =   a*fb*fc/( (fa-fb)*(fa-fc) ) &
              + b*fa*fc/( (fb-fa)*(fb-fc) ) &
              + c*fa*fb/( (fc-fa)*(fc-fb) )
        else
          ! secant rule
          s = b - fb*(b-a)/(fb-fa)
        end if

        dum = (d_three*a+b)*d_rfour
        smb = abs(s-b)
        bmc = abs(b-c)
        cmd = abs(c-d)
        if ( (.not. (s > dum .and. s < b)) .or.              &
             (mflag .and. (smb >= (bmc/d_two))) .or.         &
             ((.not. mflag) .and. (smb >= (cmd/d_two))) .or. &
             (mflag .and. (bmc < delta)) .or.                &
             ((.not. mflag) .and. (cmd < delta)) ) then
          s = (a+b)/d_two
          mflag = .true.
        else
          mflag = .false.
        end if

        fs = zerofunc(s,.false.)
        d = c
        c = b

        if ( fa*fs < d_zero ) then
          b = s
        else
          a = s
        end if

        ! check for convergence; if we converged, return from this function
        if ( abs(fs) < dlowval ) then
          b = s
          fb = fs
        end if

        if ( (abs(fa-fb) < mindq) .or. (abs(fb) < dlowval) ) then
          solve_for_t = b
          call getqvqc(solve_for_t,es,outqv,outqc)
          return
        end if
        if ( abs(fa) < abs(fb) ) call swap(a,b)
        fa = zerofunc(a,.false.)
        fb = zerofunc(b,.false.)
        fc = zerofunc(c,.false.)
      end do

      ! if we ran out of iterations, return what we have and hope for the best.
      if ( fb <= fs ) then
        solve_for_t = b
        call getqvqc(solve_for_t,es,outqv,outqc)
      else
        solve_for_t = s
        call getqvqc(solve_for_t,es,outqv,outqc)
      end if

    else if ( imethod == 2 ) then
      ! Bretherton's iterative method
      ! set the dummy temperature temps (which will eventually
      ! be the actual temperature)
      temps = templ

      !*******************************************************
      !******* condense any water, if possible ***************
      !*******************************************************

      itqt = 0
      bigloop : &
      do
        ! calculate the saturation mixing ratio
        rvls = ep2/(myp/pfesat(temps)-d_one)
        ! go through 3 iterations of calculating the saturation
        ! humidity and updating the temperature

        ! increase the number of iterations for high total water
        ! content; the algorithm takes longer to converge for high liquid
        ! water content.  for too few iterations, really high liquid water
        ! contents (and correspondingly high temperatures) result.

        do iteration = 1 , itmax ! condenseliquid
          ! update the dummy temperature
          deltat = ((templ-temps)*cpowlhv + myqt-rvls)/   &
                    (cpowlhv+ep2*wlhv*rvls/rgas/temps/temps)
          temps = temps + deltat
          ! re-calculate the saturation mixing ratio
          rvls = ep2/(myp/pfesat(temps)-d_one)
          if ( abs(deltat) < mindt ) exit
        end do ! condenseliquid
        ! cloud water is the total minus the saturation
        ddq = max(myqt-rvls,d_zero)
        outqc = ddq
        ! water vapor is anything left over (should be =rvls)
        outqv = myqt - ddq
        ! add the enthalpy of condensation to templ and
        ! convert to potential temperature
        solve_for_t = templ + ddq*wlhvocp
        ! if there is cloud, check that qv is consistent with t
        if ( ddq > d_zero ) then
          dum = pfesat(solve_for_t)
          rvls = ep2/(myp/dum-d_one)
          dum = rvls - outqv
          ! if the solution did not converge, add three more iterations
          if ( abs(dum) > mindq ) then
            ! if the solution really is not converging, issue a warning
            if ( itqt > imax ) then
              write(stderr,*) '(mod_thetal) warning: non-convergence of ', &
                         'temperature solution'
              call fatal(__FILE__,__LINE__, &
                         'model stops for UW PBL error')
            end if
            itqt = itqt + 1
            cycle bigloop
          end if
        end if
        exit
      end do bigloop

    else

      ! O'brien's finite difference method

      myqt = qt
      dthetal = thetal - thlprev
      dqt = qt - qtprev
      qvprev = myqt - qcprev
      tempt = tprev
      temps = templ
      tempqv = qvprev
      tempqc = qcprev
      myqt = qtprev
      call getqvqc(temps,tempes,tempqv,tempqc)

      ! set the threshold for (approximately) zero cloud water
      ! have it be consistent with the threshold for determining
      ! when tprev+deltat is consistent enough with thetal
      ! qcthresh = mindt*myexner/wlhvocp
      qcthresh = minqx

      do iteration = 1 , imax
        ! if the air is undersaturated, then thetal is just the normal
        ! potential temperature, so the delta t is given by dthetal
        ! if(tempqc < qcthresh) then
        if ( abs(qt-qvprev) <= qcthresh ) then
          deltat = myexner*dthetal
          ! otherwise, add a contribution from the change in cloud water.
        else
          dum = (myp/ep2/tempes)*(cpd/rwat)*      &
                (wlhvocp*tempqv/(cpd*tempt))**2 + d_one
          deltat = (myexner*dthetal + wlhvocp*dqt)/dum
        end if

        if ( abs(deltat) < mindt ) then
          tempt = tprev
          tempqv = qvprev
          tempqc = qcprev
          exit
        end if

        ! update the dummy temperature
        tempt = tprev + deltat
        ! and get qc and qv from it
        call getqvqc(tempt,tempes,tempqv,tempqc)
        ! check if the solution has converged
        dum = thetal - (tempt - wlhvocp*tempqc)/myexner
        if (abs(dum) < mindt ) exit
      end do

      ! return the number of iterations that it took to converge (or return
      ! -999 if we failed to converge)
      solve_for_t = tempt
      outqv = tempqv
      outqc = tempqc
    end if

  end function solve_for_t

  ! returns the saturation vapor pressure over water in units (cb)
  ! given the input pressure in cb and temperature in k
  ! modified from buck (1981), j. app. met. v 20
  function esatw(p,t)
    implicit none
    real(rkx) , intent(in) :: p , t
    real(rkx) :: esatw
    real(rkx) :: dum , arg , tdum
    ! Limit t to reasonable values.  I believe that this is necessary because
    ! in the iterative calculation of t and qv from the liquid water
    ! temperature, the temperature can take on some crazy values before it
    ! converges. --TAO 01/05/2011
    tdum = max(100.0_rkx,min(399.99_rkx,t))
    dum = 1.0007 + 3.46e-5_rkx*p
    !arg = 17.502*(t-273.15)/(t-32.18)
    arg = 17.502_rkx*(tdum-tzero)/(tdum-32.18_rkx)
    esatw = dum*0.61121_rkx*exp(arg)
  end function esatw

  ! returns the saturation vapor pressure over ice in units (cb)
  ! given the input pressure in cb and temperature in k
  ! modified from buck (1981), j. app. met. v 20
  function esati(p,t)
    implicit none
    real(rkx) , intent(in) :: p , t
    real(rkx) :: esati
    real(rkx) :: dum , arg , tdum
    ! Limit t to reasonable values.  I believe that this is necessary because
    ! in the iterative calculation of t and qv from the liquid water
    ! temperature, the temperature can take on some crazy values before it
    ! converges. --TAO 01/05/2011
    tdum = max(100.0_rkx,min(399.99_rkx,t))
    dum = 1.0003_rkx + 4.18e-5_rkx*p
    !arg = 22.452*(t-273.15)/(t-0.6)
    arg = 22.452_rkx*(tdum-tzero)/(tdum-0.6_rkx)
    esati = dum*0.61115_rkx*exp(arg)
  end function esati

  subroutine swap(a,b)
    implicit none
    real(rkx) , intent(inout) :: a , b
    real(rkx) :: c
    c = a
    a = b
    b = c
  end subroutine swap

  subroutine getqvqc(t,es,qv,qc)
    implicit none
    real(rkx) , intent(in) :: t
    real(rkx) , intent(out) :: es , qv , qc
    es = pfesat(t)
    qv = max(ep2/(myp/es-d_one),d_zero)
    if ( myqt > qv ) then
      qc = myqt - qv
    else
      qv = myqt
    end if
  end subroutine getqvqc

end module mod_pbl_thetal

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
