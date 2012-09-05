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

  private

  real(rk8) :: myp , mythetal , myqt , myexner , mylovcp

  logical :: bderiv

  real(rk8) , parameter :: mymaxdiff = 1.0d-3
  real(rk8) , parameter :: delta = 1.0d-8
  real(rk8) , parameter :: trange = 20.0D0

  public :: solve_for_t

  contains

  real(rk8) function zerofunc(t,isice)
    implicit none
    real(rk8) , intent(in) :: t
    real(rk8) :: qc , qv , locp , es
    integer(ik4) , intent(in) :: isice

    if ( isice == 0 ) then
      es = esatw(myp,t)
      locp = wlhvocp
    else
      es = esati(myp,t)
      locp = wlhsocp
    end if
    qv = dmax1(ep2/(myp/es-d_one),d_zero)

    if ( myqt > qv ) then
      qc = myqt - qv
    else
      qv = myqt
      qc = d_zero
    end if

    zerofunc = myexner*mythetal + locp*qc - t

    ! If necessary, find the minimum of the zerofuncion, which should be
    ! the place where it equals zero if no solution is available otherwise
    if ( bderiv ) then
      zerofunc = -(cpd/rwat)*(myp/es)*(locp*qv/t)**d_two - d_one
    end if
  end function zerofunc

  ! Determine temperature from the liquid water potential temperature,
  ! total water, and pressure.

  real(rk8) function solve_for_t(thetal,qt,p,tprev,qtprev, &
                                qcprev,thlprev,imax,      &
                                imethod,isice,outqc,outqv)
    implicit none
    real(rk8) , intent(in) :: thetal , qt , p , tprev , qtprev , &
                             qcprev , thlprev
    integer(ik4) , intent(in) :: imethod , isice
    integer(ik4) , intent(inout) :: imax
    real(rk8) , intent(out) :: outqc , outqv
    real(rk8) :: a , b , c , d , s , dum
    real(rk8) :: fa , fb , fc , fs
    real(rk8) :: smb , bmc , cmd
    real(rk8) :: temps , templ , rvls
    real(rk8) :: dthetal , dqt , qvprev , deltat , es , qcthresh
    real(rk8) :: tempqv , tempqc , tempt , tempes
    logical :: mflag
    integer(ik4) :: i , iteration , itqt , itqtsupp
    integer(ik4) , parameter :: itmax = 6
    integer(ik4) , parameter :: itmin = 3

    ! set some interal module variables
    mythetal = thetal
    myqt = qt
    myp = p
    myexner = (myp/d_100)**rovcp


!*******************************************************************************
!*******************************************************************************
!************************ the brent method *************************************
!*******************************************************************************
!*******************************************************************************

    ! Use the Brent 1973 method to determine t from liquid water pot.
    ! temp., pressure, and total water
    ! (see http://en.wikipedia.org/wiki/brent's_method)

    if ( imethod == 1 ) then
      bderiv = .false.
      if ( isice == 0 ) then
        mylovcp = wlhvocp
      else
        mylovcp = wlhsocp
      end if

      ! set the bounds on t based on the change in thetal and qw:
      ! the bounds are: all of the change in qw goes in to changing
      ! cloud water at one end, and none of the change in qw goes in to it
      ! in the other
      a = tprev + myexner*(mythetal - thlprev)
      b = a + mylovcp*(myqt-qtprev)

      ! if a and b are essentially the same, then the total water did not
      ! change, and so t is given only by the change in liquid water
      ! potential temperature; return with this as the solution
      if ( abs(a-b) <= delta ) then
        solve_for_t = b
        call getqvqc(solve_for_t,es,outqv,outqc,isice)
        imax = 0
        return
      end if

      ! make sure that the dynamic boundaries above don't conflict with the
      ! thermodynamic boundaries
      if ( a > b ) call swap(a,b)
      dum = myexner*mythetal
      a = dmax1(a,dum)
      b = dmin1(b,dum+mylovcp*myqt)
      ! run these through the solution-finder
      fa = zerofunc(a,isice)
      fb = zerofunc(b,isice)
      fs = zerofunc(tprev,isice)

      ! check if we already have the solution:
      if ( fa == 0 ) then
        b = a
        fb = fa
      end if
      ! check if the previous timestep works as a solution
      if ( dabs(fs) < dlowval ) then
        b = tprev
        fb = fs
      end if
      if ( dabs(fb) < dlowval ) then
        solve_for_t = b
        call getqvqc(solve_for_t,es,outqv,outqc,isice)
        imax = 0
        return
      end if

      ! if the initial guesses are out of range, return with a bogus
      ! temperature to flag a failure
      if ( fa*fb >= d_zero ) then
        solve_for_t = dum
        solve_for_t = tprev
        call getqvqc(solve_for_t,es,outqv,outqc,isice)
        imax = -999
!       return
!       bderiv = .true.
!       fa = zerofunc(a,isice)
!       fb = zerofunc(b,isice)
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

      !***********************************************************
      !***********************************************************
      !*************** the main iteration loop *******************
      !***********************************************************
      !***********************************************************

      do i = 1 , imax
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

        fs = zerofunc(s,isice)
        d = c
        c = b

        if ( fa*fs < d_zero ) then
          b = s
        else
          a = s
        end if

        ! check for convergence; if we converged, return from this function
        if ( dabs(fs) < dlowval ) then
          b = s
          fb = fs
        end if

        if ( (abs(fa-fb) < mymaxdiff) .or. (abs(fb) < dlowval) ) then
          solve_for_t = b
          call getqvqc(solve_for_t,es,outqv,outqc,isice)
          imax = i
          return
        end if


        if ( abs(fa) < abs(fb) ) call swap(a,b)

        fa = zerofunc(a,isice)
        fb = zerofunc(b,isice)
        fc = zerofunc(c,isice)

      end do

      ! if we ran out of iterations, return what we have and hope for the best.
      if ( fb <= fs ) then
        solve_for_t = b
        call getqvqc(solve_for_t,es,outqv,outqc,isice)
      else
        solve_for_t = s
        call getqvqc(solve_for_t,es,outqv,outqc,isice)
      end if

!*******************************************************************************
!*******************************************************************************
!************************ Bretherton's iterative method*************************
!*******************************************************************************
!*******************************************************************************

    else if ( imethod == 2 ) then

      templ = mythetal*myexner
      ! set the dummy temperature temps (which will eventually
      ! be the actual temperature)
      temps = templ

      !*******************************************************
      !******* condense any water, if possible ***************
      !*******************************************************

bigloop : &
      do
        if ( isice == 0 ) then

          mylovcp = wlhvocp

          ! calculate the saturation mixing ratio

          rvls = ep2/(myp/esatw(myp,temps)-d_one)

          ! go through 3 iterations of calculating the saturation
          ! humidity and updating the temperature

          ! increase the number of iterations for high total water
          ! content; the algorithm takes longer to converge for high liquid
          ! water content.  for too few iterations, really high liquid water
          ! contents (and correspondingly high temperatures) result.

          itqtsupp = 0
          itqt = max(min(int(d_two*(dble(itmax)+log(myqt))),itmax),itmin) + &
                 itqtsupp
          do iteration = 1, itqt ! condenseliquid
            ! update the dummy temperature
            temps = temps  + ((templ-temps)*cpowlhv + myqt-rvls)/   &
                              (cpowlhv+ep2*wlhv*rvls/rgas/temps/temps)
            ! re-calculate the saturation mixing ratio
            rvls = ep2/(myp/esatw(myp,temps)-d_one)
          end do ! condenseliquid
          ! cloud water is the total minus the saturation
          outqc = dmax1(myqt-rvls,d_zero)
          ! water vapor is anything left over (should be =rvls)
          outqv= myqt - outqc
          ! add the enthalpy of condensation to templ and
          ! convert to potential temperature
          solve_for_t = (templ+outqc*wlhvocp)
        else
          mylovcp = wlhsocp
          ! calculate the saturation mixing ratio (wrt ice)
          rvls = ep2/(myp/esati(myp,temps)-d_one)
          do iteration = 1, itqt ! condenseice
            ! update the dummy temperature
            temps = temps +                                 &
                    ((templ-temps)*cpowlhs + myqt -rvls)/   &
                     (cpowlhs +ep2*wlhs*rvls/rgas/temps/temps)
            ! re-calculate the saturation mixing ratio
            rvls = ep2/(myp/esati(myp,temps)-d_one)
          end do ! condenseice
          outqc = dmax1(myqt-rvls,d_zero)
          outqv = myqt - outqc
          solve_for_t = (templ+outqc*wlhsocp)
        end if

        ! if there is cloud, check that qv is consistent with t
        if (outqc > d_zero ) then
          if ( isice == 0 ) then
            dum = esatw(myp,solve_for_t)
          else
            dum = esati(myp,solve_for_t)
          end if
          rvls = ep2/(myp/dum-d_one)
          dum = rvls-outqv
          ! if the solution did not converge, add three more iterations
          if ( abs(dum) > mymaxdiff ) then
            itqtsupp = itqtsupp + 3
            ! if the solution really is not converging, issue a warning
            if ( itqtsupp > 99 ) then
              write(*,*) '(mod_thetal) warning: non-convergence of ', &
                         'temperature solution'
            end if
            cycle bigloop
          end if
        end if
        exit
      end do bigloop


!     ! check if the solution has converged properly
!     ! check that thetal calculates properly
!     dum = mythetal - (solve_for_t - mylovcp*outqc)/myexner
!     if ( abs(dum) > mymaxdiff ) then
!       imax = -999
!     end if
!     ! check that total water is conserved
!     dum = myqt - (outqv + outqc)
!     if ( abs(dum) > mymaxdiff*myexner/mylovcp ) then
!       imax = -999
!     end if



!*******************************************************************************
!*******************************************************************************
!************************ O'brien's finite difference method *******************
!*******************************************************************************
!*******************************************************************************

    else
      dthetal = thetal-thlprev
      dqt = qt - qtprev
      qvprev = qt-qcprev

      tempt = tprev

      myqt = qtprev
      call getqvqc(temps,tempes,tempqv,tempqc,isice)
      myqt = qt

      ! set the threshold for (approximately) zero cloud water
      ! have it be consistent with the threshold for determining
      ! when tprev+deltat is consistent enough with thetal
      ! qcthresh = mymaxdiff*myexner/mylovcp
      qcthresh = d_zero

      do iteration = 1 , imax
        ! if the air is undersaturated, then thetal is just the normal
        ! potential temperature, so the delta t is given by dthetal
        ! if(tempqc < qcthresh) then 
        if ( abs(qt-qvprev) <= qcthresh ) then 
          deltat = myexner*dthetal
          ! otherwise, add a contribution from the change in cloud water.
        else
          dum = (myp/ep2/tempes)*(cpd/rwat)*      &
                (mylovcp*tempqv/(cpd*tempt))**d_two + d_one
          deltat = (myexner*dthetal + mylovcp*dqt)/dum
        end if

        if ( abs(deltat) < mymaxdiff ) then
          tempt = tprev
          tempqc = qcprev
          tempqv = qvprev
          exit
        end if

        ! update the dummy temperature
        tempt = tprev + deltat
        ! and get qc and qv from it
        call getqvqc(tempt,tempes,tempqv,tempqc,isice)

        ! check if the solution has converged
        dum = thetal - (tempt - mylovcp*tempqc)/myexner
        if (abs(dum) < mymaxdiff ) exit
      end do

      ! return the number of iterations that it took to converge (or return
      ! -999 if we failed to converge)
      if ( iteration >= imax ) then
        imax = -999
        dum = sqrt(float(imax))
      else
        imax = iteration
      end if

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
    real(rk8) , intent(in) :: p , t
    real(rk8) :: esatw
    real(rk8) :: dum , arg , tdum
    ! Limit t to reasonable values.  I believe that this is necessary because
    ! in the iterative calculation of t and qv from the liquid water
    ! temperature, the temperature can take on some crazy values before it
    ! converges. --TAO 01/05/2011
    tdum = dmax1(100.0D0,dmin1(399.99D0,t))
    dum = 1.0007 + 3.46D-5*p
    !arg = 17.502*(t-273.15)/(t-32.18)
    arg = 17.502D0*(tdum-tzero)/(tdum-32.18D0)
    esatw = dum*0.61121D0*dexp(arg)
  end function esatw

  ! returns the saturation vapor pressure over ice in units (cb)
  ! given the input pressure in cb and temperature in k
  ! modified from buck (1981), j. app. met. v 20
  function esati(p,t)
    implicit none
    real(rk8) , intent(in) :: p , t
    real(rk8) :: esati
    real(rk8) :: dum , arg , tdum
    ! Limit t to reasonable values.  I believe that this is necessary because
    ! in the iterative calculation of t and qv from the liquid water
    ! temperature, the temperature can take on some crazy values before it
    ! converges. --TAO 01/05/2011
    tdum = dmax1(100.0D0,dmin1(399.99D0,t))
    dum = 1.0003D0 + 4.18D-5*p
    !arg = 22.452*(t-273.15)/(t-0.6)
    arg = 22.452D0*(tdum-tzero)/(tdum-0.6D0)
    esati = dum*0.61115D0*dexp(arg)
  end function esati

  subroutine swap(a,b)
    implicit none
    real(rk8) , intent(inout) :: a , b
    real(rk8) :: c
    c = a
    a = b
    b = c
  end subroutine swap

  subroutine getqvqc(t,es,qv,qc,isice)
    implicit none
    real(rk8) , intent(in) :: t
    real(rk8) , intent(out) :: es , qv , qc
    integer(ik4) , intent(in) :: isice
    if ( isice == 0 ) then
      es = esatw(myp,t)
      mylovcp = wlhvocp
    else
      es = esati(myp,t)
      mylovcp = wlhsocp
    end if
    qv = dmax1(ep2/(myp/es-d_one),d_zero)
    if ( myqt > qv ) then
      qc = myqt - qv
    else
      qv = myqt
      qc = d_zero
    end if
  end subroutine getqvqc

end module mod_pbl_thetal
