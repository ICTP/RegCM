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

  real(rkx) , public :: dryini , watini

  contains

  subroutine massck
    implicit none
    real(rkx) :: error1 , error2 , tttmp
    real(rkx) :: tdrym , tdadv , tqmass , tqadv
    real(rkx) :: tcrai , tncrai , tqeva
    real(rkx) :: drymass , dryadv , qmass , qadv , craim , ncraim , evapm
    real(rkx) :: north , south , east , west
    integer(ik4) :: i , j , k , n
    character (len=32) :: appdat

    tdrym = d_zero
    tdadv = d_zero
    !
    ! Internal dry air mass
    !
    tttmp = d_zero
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
    tdadv = d_zero
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
    tqmass = d_zero
    tqadv = d_zero
    tcrai = d_zero
    tncrai = d_zero
    tqeva = d_zero
    do i = ici1 , ici2
      do j = jci1 , jci2
        tcrai = tcrai + sfs%rainc(j,i)*dxsq
        tncrai = tncrai + sfs%rainnc(j,i)*dxsq
        tqeva = tqeva + sfs%qfx(j,i)*dxsq*dt
      end do
    end do
    do n = 1 , nqx
      do k = 1 , kz
        tttmp = d_zero
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
      dryini = drymass
      watini = qmass
      return
    end if

    call sumall(tdadv,dryadv)
    call sumall(tqadv,qadv)
    call sumall(tcrai,craim)
    call sumall(tncrai,ncraim)
    call sumall(tqeva,evapm)

    drymass = drymass - dryadv
    qmass = qmass + tcrai + tncrai - qadv - evapm

    if ( myid == italk .and. mod(ktau,krep) == 0 ) then
      error1 = (drymass-dryini)/dryini * d_100
      error2 = (qmass-watini)/watini * d_100
      appdat = tochar(idatex)
      write(stdout,'(a,a23,a,i16)') ' *** ', appdat, ', ktau   = ', ktau
      write(stdout,'(a,e12.5,a,f5.2,a)') '   total air   =', drymass, &
                 ' kg, error = ', error1, ' %'
      write(stdout,'(a,e12.5,a)') ' horizontal advection   = ', dryadv , ' kg.'
      write(stdout,'(a,e12.5,a,f5.2,a)') '   total water =', qmass, &
                 ' kg, error = ', error2, ' %'
      write(stdout,'(a,e12.5,a)') ' horizontal advection    = ', qadv , ' kg.'
      write(stdout,'(a,e12.5,a)') ' convective rainfall     = ', tcrai , ' kg.'
      write(stdout,'(a,e12.5,a)') ' nonconvective rainfall  = ', tncrai , ' kg.'
      write(stdout,'(a,e12.5,a)') ' evaporation from ground = ', evapm , ' kg.'
    end if

    dryini = drymass
    watini = qmass

  end subroutine massck

end module mod_massck

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
