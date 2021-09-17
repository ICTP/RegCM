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

module mod_massck
  !
  ! Compute the total dry air and water substance within the domain
  ! and compares with the initial values.
  ! The unit used in all the calculation is "kg".
  !
  use mod_intkinds
  use mod_realkinds
  use mod_stdio
  use mod_dynparam
  use mod_constants
  use mod_runparams
  use mod_mppparam
  use mod_atm_interface

  implicit none

  private

  public :: massck

#ifndef QUAD_PRECISION
  integer , parameter :: wrkp = rk8
#else
  integer , parameter :: wrkp = rk16
#endif

  real(wrkp) , parameter :: q_zero = 0.0_wrkp
  real(wrkp) , public :: dryini = q_zero
  real(wrkp) , public :: watini = q_zero
  real(rk8) , public :: dryerror = d_zero
  real(rk8) , public :: waterror = d_zero

  contains

  subroutine massck
    implicit none
    real(wrkp) :: tttmp
    real(wrkp) :: tdrym , tdadv , tqmass , tqadv
    real(wrkp) :: tcrai , tncrai , tqeva
    real(wrkp) :: drymass , dryadv , qmass , qadv , craim , ncraim , evapm
    real(wrkp) :: north , south , east , west
    real(wrkp) , save :: mcrai = 0.0_wrkp
    real(wrkp) , save :: mncrai = 0.0_wrkp
    real(wrkp) , save :: mevap = 0.0_wrkp
    real(wrkp) , save :: mdryadv = 0.0_wrkp
    real(wrkp) , save :: mqadv = 0.0_wrkp
    real(wrkp) :: w1 , w2
    integer(ik4) :: i , j , k , n

    w2 = real(dt,wrkp)/86400.0_wrkp
    w1 = 1.0_wrkp - w2

    tdrym = q_zero
    tdadv = q_zero
    drymass = q_zero
    dryadv = q_zero
    !
    ! Internal dry air mass
    !
    if ( idynamic == 3 ) then
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            tdrym = tdrym + dxsq * mo_atm%dz(j,i,k) * mo_atm%rho(j,i,k)
          end do
        end do
      end do
      tdadv = q_zero
      if ( ma%has_bdyleft ) then
        do k = 1 , kz
          do i = ice1 , ice2
            tdadv = tdadv + mo_atm%u(jde1,i,k) * dt * dx * &
              mo_atm%dz(jce1,i,k) * mo_atm%rho(jce1,i,k)
          end do
        end do
      end if
      if ( ma%has_bdyright ) then
        do k = 1 , kz
          do i = ice1 , ice2
            tdadv = tdadv - mo_atm%u(jde2,i,k) * dt * dx * &
              mo_atm%dz(jce2,i,k) * mo_atm%rho(jce2,i,k)
          end do
        end do
      end if
      if ( ma%has_bdybottom ) then
        do k = 1 , kz
          do j = jci1 , jci2
            tdadv = tdadv + mo_atm%v(j,ide1,k) * dt * dx * &
              mo_atm%dz(j,ice1,k) * mo_atm%rho(j,ice1,k)
          end do
        end do
      end if
      if ( ma%has_bdytop ) then
        do k = 1 , kz
          do j = jci1 , jci2
            tdadv = tdadv - mo_atm%v(j,ide2,k) * dt * dx * &
              mo_atm%dz(j,ice2,k) * mo_atm%rho(j,ice2,k)
          end do
        end do
      end if
      !
      ! Moisture
      !
      tqmass = q_zero
      tqadv = q_zero
      tcrai = q_zero
      tncrai = q_zero
      tqeva = q_zero
      qmass = q_zero
      qadv = q_zero
      craim = q_zero
      ncraim = q_zero
      evapm = q_zero
      do i = ici1 , ici2
        do j = jci1 , jci2
          tcrai = tcrai + crrate(j,i)*dxsq*dt
          tncrai = tncrai + ncrrate(j,i)*dxsq*dt
          tqeva = tqeva + sfs%qfx(j,i)*dxsq*dt
        end do
      end do
      do n = 1 , nqx
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
              tqmass = tqmass + mo_atm%qx(j,i,k,n) * &
                dxsq * mo_atm%dz(j,i,k) * mo_atm%rho(j,i,k)
            end do
          end do
        end do
        !
        ! Boundary input
        !
        if ( ma%has_bdyleft ) then
          do k = 1 , kz
            do i = ice1 , ice2
              tqadv = tqadv + mo_atm%qx(jce1,i,k,n) * &
                mo_atm%u(jde1,i,k) * dt * dx * &
                mo_atm%dz(jce1,i,k) * mo_atm%rho(jce1,i,k)
            end do
          end do
        end if
        if ( ma%has_bdyright ) then
          do k = 1 , kz
            do i = ice1 , ice2
              tqadv = tqadv - mo_atm%qx(jce2,i,k,n) * &
                mo_atm%u(jde2,i,k) * dt * dx * &
                mo_atm%dz(jce2,i,k) * mo_atm%rho(jce2,i,k)
            end do
          end do
        end if
        if ( ma%has_bdybottom ) then
          do k = 1 , kz
            do j = jci1 , jci2
              tqadv = tqadv + mo_atm%qx(j,ice1,k,n) * &
                mo_atm%v(j,ide1,k) * dt * dx * &
                mo_atm%dz(j,ice1,k) * mo_atm%rho(j,ice1,k)
            end do
          end do
        end if
        if ( ma%has_bdytop ) then
          do k = 1 , kz
            do j = jci1 , jci2
              tqadv = tqadv - mo_atm%qx(j,ice2,k,n) * &
                mo_atm%v(j,ide2,k) * dt * dx * &
                mo_atm%dz(j,ice2,k) * mo_atm%rho(j,ice2,k)
            end do
          end do
        end if
      end do
    else
      tttmp = q_zero
      do i = ici1 , ici2
        do j = jci1 , jci2
          tttmp = tttmp + sfs%psa(j,i)
        end do
      end do
      do k = 1 , kz
        tdrym = tdrym + tttmp*dsigma(k)
      end do
      !
      ! Boundary input
      !
      tdadv = q_zero
      if ( ma%has_bdyleft ) then
        do k = 1 , kz
          do i = ice1 , ice2
            west = (atm1%u(jde1,i+1,k)+atm1%u(jde1,i,k))
            tdadv = tdadv + dt*3.e4_wrkp*dsigma(k)*dx*regrav*west
          end do
        end do
      end if
      if ( ma%has_bdyright ) then
        do k = 1 , kz
          do i = ice1 , ice2
            east = (atm1%u(jde2,i+1,k)+atm1%u(jde2,i,k))
            tdadv = tdadv - dt*3.e4_wrkp*dsigma(k)*dx*regrav*east
          end do
        end do
      end if
      if ( ma%has_bdybottom ) then
        do k = 1 , kz
          do j = jci1 , jci2
            south = (atm1%v(j+1,ide1,k)+atm1%v(j,ide1,k))
            tdadv = tdadv + dt*3.e4_wrkp*dsigma(k)*dx*regrav*south
          end do
        end do
      end if
      if ( ma%has_bdytop ) then
        do k = 1 , kz
          do j = jci1 , jci2
            north = (atm1%v(j+1,ide2,k)+atm1%v(j,ide2,k))
            tdadv = tdadv - dt*3.e4_wrkp*dsigma(k)*dx*regrav*north
          end do
        end do
      end if
      tdrym = tdrym*dxsq*1000.0_wrkp*regrav
      !
      ! Moisture
      !
      tqmass = q_zero
      tqadv = q_zero
      tcrai = q_zero
      tncrai = q_zero
      tqeva = q_zero
      qmass = q_zero
      qadv = q_zero
      craim = q_zero
      ncraim = q_zero
      evapm = q_zero
      do i = ici1 , ici2
        do j = jci1 , jci2
          tcrai = tcrai + crrate(j,i)*dxsq*dt
          tncrai = tncrai + ncrrate(j,i)*dxsq*dt
          tqeva = tqeva + sfs%qfx(j,i)*dxsq*dt
        end do
      end do
      do n = 1 , nqx
        do k = 1 , kz
          tttmp = q_zero
          do i = ici1 , ici2
            do j = jci1 , jci2
              tttmp = tttmp + atm1%qx(j,i,k,n)
            end do
          end do
          tqmass = tqmass + tttmp*dsigma(k)
        end do
        !
        ! Boundary input
        !
        if ( ma%has_bdyleft ) then
          do k = 1 , kz
            do i = ice1 , ice2
              west = (atm1%u(jde1,i+1,k)+atm1%u(jde1,i,k)) * &
                   (atm1%qx(jce1,i,k,n)/sfs%psa(jce1,i))
              tqadv = tqadv + dt*3.e4_wrkp*dsigma(k)*dx*regrav*west
            end do
          end do
        end if
        if ( ma%has_bdyright ) then
          do k = 1 , kz
            do i = ice1 , ice2
              east = (atm1%u(jde2,i+1,k)+atm1%u(jde2,i,k)) * &
                   (atm1%qx(jce2,i,k,n)/sfs%psa(jce2,i))
              tqadv = tqadv - dt*3.e4_wrkp*dsigma(k)*dx*regrav*east
            end do
          end do
        end if
        if ( ma%has_bdybottom ) then
          do k = 1 , kz
            do j = jci1 , jci2
              south = (atm1%v(j+1,ide1,k)+atm1%v(j,ide1,k)) * &
                   (atm1%qx(j,ice1,k,n)/sfs%psa(j,ice1))
              tqadv = tqadv + dt*3.e4_wrkp*dsigma(k)*dx*regrav*south
            end do
          end do
        end if
        if ( ma%has_bdytop ) then
          do k = 1 , kz
            do j = jci1 , jci2
              north = (atm1%v(j+1,ide2,k)+atm1%v(j,ide2,k)) * &
                   (atm1%qx(j,ice2,k,n)/sfs%psa(j,ice2))
              tqadv = tqadv - dt*3.e4_wrkp*dsigma(k)*dx*regrav*north
            end do
          end do
        end if
      end do
      tqmass = tqmass*dxsq*1000.0_wrkp*regrav
    end if

    call sumall(tdrym,drymass)
    call sumall(tqmass,qmass)

    if ( rcmtimer%start( ) ) then
      dryerror = d_zero
      waterror = d_zero
      if ( myid == italk ) then
        dryini = drymass
        watini = qmass
      end if
      return
    end if

    call sumall(tqadv,qadv)
    call sumall(tdadv,dryadv)
    call sumall(tcrai,craim)
    call sumall(tncrai,ncraim)
    call sumall(tqeva,evapm)

    mcrai = w1*mcrai + w2*craim
    mncrai = w1*mncrai + w2*ncraim
    mevap = w1*mevap + w2*evapm
    mdryadv = w1*mdryadv + w2*dryadv
    mqadv = w1*mqadv + w2*qadv

    if ( myid == italk ) then
      drymass = drymass - dryadv
      qmass = qmass + craim + ncraim - qadv - evapm
      if ( dryini < dlowval ) dryini = drymass
      if ( watini < dlowval ) watini = qmass
      dryerror = dryerror + &
        (real((drymass-dryini)/dryini,rk8) * d_100) * dt/86400.0_rk8
      waterror = waterror + &
        (real((qmass-watini)/watini,rk8) * d_100) * dt/86400.0_rk8
      if ( alarm_day%act( ) .or. syncro_dbg%act( ) ) then
        write(stdout,'(a)') &
            ' ********************* MASS CHECK ********************'
        write(stdout,*) 'At ', trim(rcmtimer%str( ))
        write(stdout,'(a,e12.5,a,f9.5,a)') ' Total dry air   =', drymass, &
                   ' kg, error = ', dryerror, ' %'
        write(stdout,'(a,e12.5,a,f9.5,a)') ' Total water     =', qmass, &
                   ' kg, error = ', waterror, ' %'
        write(stdout,'(a)') ' Mean values over past 24 hours :'
        write(stdout,'(a,e12.5,a)') ' Dry air boundary    = ', mdryadv , ' kg.'
        write(stdout,'(a,e12.5,a)') ' Water boundary      = ', mqadv, ' kg.'
        write(stdout,'(a,e12.5,a)') ' Convective rain     = ', mcrai, ' kg.'
        write(stdout,'(a,e12.5,a)') ' Nonconvective rain  = ', mncrai, ' kg.'
        write(stdout,'(a,e12.5,a)') ' Ground Evaporation  = ', mevap, ' kg.'
        write(stdout,'(a)') &
            ' *****************************************************'
        dryerror = 0.0_rk8
        waterror = 0.0_rk8
      end if
      dryini = drymass
      watini = qmass
    end if

  end subroutine massck

end module mod_massck

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
