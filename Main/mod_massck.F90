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
  use mod_dynparam
  use mod_constants
  use mod_runparams
  use mod_mppparam
  use mod_atm_interface

  implicit none

  private

  public :: massck

  real(rk16) , parameter :: q_zero = 0.0_rk16
  real(rk16) , public :: dryini , watini

  contains

  subroutine massck
    implicit none
    real(rkx) :: error1 , error2
    real(rk16) :: tttmp
    real(rk16) :: tdrym , tdadv , tqmass , tqadv
    real(rk16) :: tcrai , tncrai , tqeva
    real(rk16) :: drymass , dryadv , qmass , qadv , craim , ncraim , evapm
    real(rk16) :: north , south , east , west
    integer(ik4) :: i , j , k , n
    character (len=32) :: appdat

    tdrym = q_zero
    tdadv = q_zero
    !
    ! Internal dry air mass
    !
    tttmp = q_zero
    do i = ice1 , ice2
      do j = jce1 , jce2
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
          tdadv = tdadv + dt*3.e4_rkx*dsigma(k)*dx*regrav*west
        end do
      end do
    end if
    if ( ma%has_bdyright ) then
      do k = 1 , kz
        do i = ice1 , ice2
          east = (atm1%u(jde2,i+1,k)+atm1%u(jde2,i,k))
          tdadv = tdadv - dt*3.e4_rkx*dsigma(k)*dx*regrav*east
        end do
      end do
    end if
    if ( ma%has_bdybottom ) then
      do k = 1 , kz
        do j = jce1 , jce2
          south = (atm1%v(j+1,ide1,k)+atm1%v(j,ide1,k))
          tdadv = tdadv + dt*3.e4_rkx*dsigma(k)*dx*regrav*south
        end do
      end do
    end if
    if ( ma%has_bdytop ) then
      do k = 1 , kz
        do j = jce1 , jce2
          north = (atm1%v(j+1,ide2,k)+atm1%v(j,ide2,k))
          tdadv = tdadv - dt*3.e4_rkx*dsigma(k)*dx*regrav*north
        end do
      end do
    end if
    tdrym = tdrym*dxsq*1000.0_rkx*regrav
    !
    ! Moisture
    !
    tqmass = q_zero
    tqadv = q_zero
    tcrai = q_zero
    tncrai = q_zero
    tqeva = q_zero
    do i = ici1 , ici2
      do j = jci1 , jci2
        tcrai = tcrai + pptnc(j,i)*dxsq*dt
        tncrai = tncrai + pptc(j,i)*dxsq*dt
        tqeva = tqeva + sfs%qfx(j,i)*dxsq*dt
      end do
    end do
    do n = 1 , nqx
      do k = 1 , kz
        tttmp = q_zero
        do i = ice1 , ice2
          do j = jce1 , jce2
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
            tqadv = tqadv + dt*3.e4_rkx*dsigma(k)*dx*regrav*west
          end do
        end do
      end if
      if ( ma%has_bdyright ) then
        do k = 1 , kz
          do i = ice1 , ice2
            east = (atm1%u(jde2,i+1,k)+atm1%u(jde2,i,k)) * &
                 (atm1%qx(jce2,i,k,n)/sfs%psa(jce2,i))
            tqadv = tqadv - dt*3.e4_rkx*dsigma(k)*dx*regrav*east
          end do
        end do
      end if
      if ( ma%has_bdybottom ) then
        do k = 1 , kz
          do j = jce1 , jce2
            south = (atm1%v(j+1,ide1,k)+atm1%v(j,ide1,k)) * &
                 (atm1%qx(j,ice1,k,n)/sfs%psa(j,ice1))
            tqadv = tqadv + dt*3.e4_rkx*dsigma(k)*dx*regrav*south
          end do
        end do
      end if
      if ( ma%has_bdytop ) then
        do k = 1 , kz
          do j = jce1 , jce2
            north = (atm1%v(j+1,ide2,k)+atm1%v(j,ide2,k)) * &
                 (atm1%qx(j,ice2,k,n)/sfs%psa(j,ice2))
            tqadv = tqadv - dt*3.e4_rkx*dsigma(k)*dx*regrav*north
          end do
        end do
      end if
    end do
    tqmass = tqmass*dxsq*1000.0_rkx*regrav

    call sumall(tdrym,drymass)
    call sumall(tqmass,qmass)

    if ( ktau == 0 ) then
      return
    end if

    call sumall(tdadv,dryadv)
    call sumall(tqadv,qadv)
    call sumall(tcrai,craim)
    call sumall(tncrai,ncraim)
    call sumall(tqeva,evapm)

    if ( myid == italk ) then
      drymass = drymass - dryadv
      qmass = qmass + tcrai + tncrai - qadv - evapm
      error1 = error1 + &
        (real((drymass-dryini)/dryini,rkx) * d_100) * dt/86400.0_rkx
      error2 = error2 + &
        (real((qmass-watini)/watini,rkx) * d_100) * dt/86400.0_rkx
      if ( mod(ktau,kday) == 0 ) then
        appdat = tochar(idatex)
        write(stdout,'(a)') &
            ' ********************* MASS CHECK ********************'
        write(stdout,'(a,a23,a,i16)') ' At ', appdat, ', ktau   = ', ktau
        write(stdout,'(a,e12.5,a,f9.4,a)') ' Total dry air   =', drymass, &
                   ' kg, error = ', error1, ' %'
        write(stdout,'(a,e12.5,a,f9.4,a)') ' Total water     =', qmass, &
                   ' kg, error = ', error2, ' %'
        write(stdout,'(a,e12.5,a)') ' Dry air boundary    = ', dryadv , ' kg.'
        write(stdout,'(a,e12.5,a)') ' Water boundary      = ', qadv, ' kg.'
        write(stdout,'(a,e12.5,a)') ' Convective rain     = ', tcrai, ' kg.'
        write(stdout,'(a,e12.5,a)') ' Nonconvective rain  = ', tncrai, ' kg.'
        write(stdout,'(a,e12.5,a)') ' Ground Evaporation  = ', evapm, ' kg.'
        write(stdout,'(a)') &
            ' *****************************************************'
        error1 = d_zero
        error2 = d_zero
      end if
      dryini = drymass
      watini = qmass
    end if

  end subroutine massck

end module mod_massck

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
