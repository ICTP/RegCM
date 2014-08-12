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

module mod_sladvection

  ! To compute the horizontal advection using semi-lagrangian

  use mod_realkinds
  use mod_intkinds
  use mod_constants
  use mod_dynparam
  use mod_runparams
  use mod_atm_interface
  use mod_sldepparam
  use mod_mpmessage
  use mod_mppparam
  use mod_service

  implicit none

  private

  public :: trajcalc_x , trajcalc_d
  public :: slhadv_x, slhadv_d , hdvg_x , hdvg_d

  interface slhadv_x
    module procedure slhadv_x3d
    module procedure slhadv_x4d
  end interface

  interface hdvg_x
    module procedure hdvg_x3d
    module procedure hdvg_x4d
  end interface hdvg_x

  contains

  subroutine adv_velocity(ldot)
    implicit none
    logical , intent(in) :: ldot
    integer(ik4) :: i , j , k
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'adv_velocity'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    if ( ldot ) then
      do k = 1, kz
        do i = idi1 , idi2
          do j = jdi1 , jdi2
            uadvx_d(j,i,k)  = atmx%u(j,i,k)/mddom%msfd(j,i)
            vadvy_d(j,i,k)  = atmx%v(j,i,k)/mddom%msfd(j,i)
            uadxp1_d(j,i,k) = atmx%u(j+1,i,k)/mddom%msfd(j+1,i)
            uadxm1_d(j,i,k) = atmx%u(j-1,i,k)/mddom%msfd(j-1,i)
            vadyp1_d(j,i,k) = atmx%v(j,i+1,k)/mddom%msfd(j,i+1)
            vadym1_d(j,i,k) = atmx%v(j,i-1,k)/mddom%msfd(j,i-1)
          end do
        end do
      end do
    else
      do k = 1, kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            uadvx_x(j,i,k) = 0.25D0 * (atmx%u(j,i,k) + &
                 atmx%u(j,i+1,k) + atmx%u(j+1,i+1,k) + &
                 atmx%u(j+1,i,k)) / mddom%msfx(j,i)
            vadvy_x(j,i,k) = 0.25D0 * (atmx%v(j,i,k) + &
                 atmx%v(j,i+1,k) + atmx%v(j+1,i+1,k) + &
                 atmx%v(j+1,i,k)) / mddom%msfx(j,i)
            uadxp1_x(j,i,k) = 0.25D0 * (atmx%u(j+1,i,k) + &
                 atmx%u(j+1,i+1,k) + atmx%u(j+2,i+1,k)  + &
                 atmx%u(j+2,i,k)) / mddom%msfx(j+1,i)
            uadxm1_x(j,i,k) = 0.25D0 * (atmx%u(j,i,k) + &
                 atmx%u(j,i+1,k) + atmx%u(j-1,i+1,k)  + &
                 atmx%u(j-1,i,k)) / mddom%msfx(j-1,i)
            vadyp1_x(j,i,k) = 0.25D0 * (atmx%v(j,i+1,k) + &
                 atmx%v(j+1,i+1,k) + atmx%v(j+1,i+2,k)  + &
                 atmx%v(j,i+2,k)) / mddom%msfx(j,i+1)
            vadym1_x(j,i,k) = 0.25D0 * (atmx%v(j,i,k) + &
                 atmx%v(j,i-1,k) + atmx%v(j+1,i,k)    + &
                 atmx%v(j+1,i,k)) / mddom%msfx(j,i-1)
          end do
        end do
      end do
    end if
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine adv_velocity

  subroutine trajcalc_x
    implicit none
    real(rk8) :: ux , uxx , xdis , xn , alfax , vy , vyy , ydis , &
                 yn , betay , ddx , ddy
    integer(ik4) :: i , j , k , xnp , xsn , ynp , ysn
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'trajcalc_x'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    ! get the advective velocity

    call adv_velocity(.false.)

    ddx = dx
    ddy = dx

    do k = 1, kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          ux = 0.5D0 * (uadxp1_x(j,i,k) - uadxm1_x(j,i,k))/ddx
          uxx = (uadxp1_x(j,i,k) - &
                  2.0D0*uadvx_x(j,i,k) + uadxm1_x(j,i,k))/(ddx*ddx)
          xdis = - uadvx_x(j,i,k)*dt + 0.5D0 * (dtsq*uadvx_x(j,i,k)*ux) - &
                   (dtcb*uadvx_x(j,i,k))*(ux*ux+uadvx_x(j,i,k)*uxx)/6.0D0
          xn = xdis/ddx
          xnp = int(xn)
          if ( abs(xnp) >= 3 ) then
            write(stderr,*) 'SL Advection problem in WE direction !'
            write(stderr,*) 'Cannot compute DP at J = ',global_dot_jstart+j
            write(stderr,*) '                     I = ',global_dot_istart+i
            call fatal(__FILE__,__LINE__,'SLADVECTION')
          end if
          alfax = dabs((xnp*ddx - xdis)/ddx)
          xsn = int(sign(d_one,xn))
          xndp_x(j,i,k) = j + xnp
          xnnm1dp_x(j,i,k) = xndp_x(j,i,k) + xsn
          xnnm2dp_x(j,i,k) = xnnm1dp_x(j,i,k) + xsn
          xnnp1dp_x(j,i,k) = xndp_x(j,i,k) - xsn
          if ( ma%has_bdyleft ) then
            if ( xndp_x(j,i,k) < jce1 ) xndp_x(j,i,k) = jce1
            if ( xnnm1dp_x(j,i,k) < jce1 ) xnnm1dp_x(j,i,k) = jce1
            if ( xnnm2dp_x(j,i,k) < jce1 ) xnnm2dp_x(j,i,k) = jce1
            if ( xnnp1dp_x(j,i,k) < jce1 ) xnnp1dp_x(j,i,k) = jce1
          end if
          if ( ma%has_bdyright ) then
            if ( xndp_x(j,i,k) > jce2 ) xndp_x(j,i,k) = jce2
            if ( xnnm1dp_x(j,i,k) > jce2 ) xnnm1dp_x(j,i,k) = jce2
            if ( xnnm2dp_x(j,i,k) > jce2 ) xnnm2dp_x(j,i,k) = jce2
            if ( xnnp1dp_x(j,i,k) > jce2 ) xnnp1dp_x(j,i,k) = jce2
          end if

          ! for the y (meridional) direction

          vy = 0.5D0*(vadyp1_x(j,i,k) - vadym1_x(j,i,k))/ddy
          vyy = (vadyp1_x(j,i,k) - &
                  2.0D0*vadvy_x(j,i,k) + vadym1_x(j,i,k))/(ddy*ddy)
          ydis = - vadvy_x(j,i,k)*dt + 0.5D0*(dtsq*vadvy_x(j,i,k)*vy) - &
                  (dtcb*vadvy_x(j,i,k))*(vy*vy +vadvy_x(j,i,k)*vyy)/6.0D0
          ! GTD ydis = - vadvy_x(j,i,k)*dt
          yn = ydis/ddy
          ynp = int(yn)
          if ( abs(ynp) >= 3 ) then
            write(stderr,*) 'SL Advection problem in SN direction !'
            write(stderr,*) 'Cannot compute DP at J = ',global_dot_jstart+j
            write(stderr,*) '                     I = ',global_dot_istart+i
            call fatal(__FILE__,__LINE__,'SLADVECTION')
          end if
          betay = dabs((ynp*ddy - ydis)/ddy)
          ysn = int(sign(d_one,yn))
          yndp_x(j,i,k) = i + ynp
          ynnm1dp_x(j,i,k) = yndp_x(j,i,k) + ysn
          ynnm2dp_x(j,i,k) = ynnm1dp_x(j,i,k) + ysn
          ynnp1dp_x(j,i,k) = yndp_x(j,i,k) - ysn
          if ( ma%has_bdybottom ) then
            if ( yndp_x(j,i,k) < ice1 ) yndp_x(j,i,k) = ice1
            if ( ynnm1dp_x(j,i,k) < ice1 ) ynnm1dp_x(j,i,k) = ice1
            if ( ynnm2dp_x(j,i,k) < ice1 ) ynnm2dp_x(j,i,k) = ice1
            if ( ynnp1dp_x(j,i,k) < ice1 ) ynnp1dp_x(j,i,k) = ice1
          end if
          if ( ma%has_bdytop ) then
            if ( yndp_x(j,i,k) > ice2 ) yndp_x(j,i,k) = ice2
            if ( ynnm1dp_x(j,i,k) > ice2 ) ynnm1dp_x(j,i,k) = ice2
            if ( ynnm2dp_x(j,i,k) > ice2 ) ynnm2dp_x(j,i,k) = ice2
            if ( ynnp1dp_x(j,i,k) > ice2 ) ynnp1dp_x(j,i,k) = ice2
          end if

          ! to estimate the value at the departure point at t=tau -1
          ! using cubic interpolation

          ! weighting coefficients

          alffbl_x(j,i,k) = alfax
          alfm2dp_x(j,i,k) = -(alfax*(1.0D0 - alfax*alfax))/6.0D0
          alfm1dp_x(j,i,k) = (alfax*(1.0D0 + alfax)*(2.0D0 - alfax))/2.0D0
          alfdp_x(j,i,k) = ((1.0D0 - alfax*alfax)*(2.0D0 - alfax))/2.0D0
          alfp1dp_x(j,i,k) = -(alfax*(1.0D0 - alfax)*(2.0D0 - alfax))/6.0D0
          betm2dp_x(j,i,k) = -(betay*(1.0D0 - betay*betay))/6.0D0
          betm1dp_x(j,i,k) = (betay*(1.0D0 + betay)*(2.0D0 - betay))/2.0D0
          betdp_x(j,i,k) = ((1.0D0 - betay*betay)*(2.0D0 - betay))/2.0D0
          betp1dp_x(j,i,k) = -(betay*(1.0D0 - betay)*(2.0D0 - betay))/6.0D0
        end do
      end do
    end do
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine trajcalc_x

  subroutine trajcalc_d
    implicit none
    real(rk8) :: ux , uxx , xdis , xn , alfax , vy , vyy , ydis , &
                 yn , betay , ddx , ddy
    integer(ik4) :: i , j , k , xnp , xsn , ynp , ysn
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'trajcalc_d'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    ! get the advective velocity

    call adv_velocity(.true.)

    ddx = dx
    ddy = dx

    do k = 1, kz
      do i = idi1 , idi2
        do j = jdi1 , jdi2
          ux = 0.5D0 * (uadxp1_d(j,i,k) - uadxm1_d(j,i,k))/ddx
          uxx = (uadxp1_d(j,i,k) - &
                  2.0D0*uadvx_d(j,i,k) + uadxm1_d(j,i,k))/(ddx*ddx)
          xdis = - uadvx_d(j,i,k)*dt + 0.5D0* (dtsq*uadvx_d(j,i,k)*ux) - &
                   (dtcb*uadvx_d(j,i,k))*(ux*ux+uadvx_d(j,i,k)*uxx)/6.0D0
          xn = xdis/ddx
          xnp = int(xn)
          if ( abs(xnp) >= 3 ) then
            write(stderr,*) 'SL Advection problem in WE direction.'
            call fatal(__FILE__,__LINE__,'SLADVECTION')
          end if
          alfax = abs((xnp*ddx -xdis)/ddx)
          xsn = int(sign(d_one,xn))
          xndp_d(j,i,k) = j + xnp
          xnnm1dp_d(j,i,k) = xndp_d(j,i,k) + xsn
          xnnm2dp_d(j,i,k) = xnnm1dp_d(j,i,k) + xsn
          xnnp1dp_d(j,i,k) = xndp_d(j,i,k) - xsn
          if ( ma%has_bdyleft ) then
            if ( xndp_d(j,i,k) < jde1 ) xndp_x(j,i,k) = jde1
            if ( xnnm1dp_d(j,i,k) < jde1 ) xnnm1dp_x(j,i,k) = jde1
            if ( xnnm2dp_d(j,i,k) < jde1 ) xnnm2dp_x(j,i,k) = jde1
            if ( xnnp1dp_d(j,i,k) < jde1 ) xnnp1dp_x(j,i,k) = jde1
          end if
          if ( ma%has_bdyright ) then
            if ( xndp_d(j,i,k) > jde2 ) xndp_x(j,i,k) = jde2
            if ( xnnm1dp_d(j,i,k) > jde2 ) xnnm1dp_x(j,i,k) = jde2
            if ( xnnm2dp_d(j,i,k) > jde2 ) xnnm2dp_x(j,i,k) = jde2
            if ( xnnp1dp_d(j,i,k) > jde2 ) xnnp1dp_x(j,i,k) = jde2
          end if

          ! for the y (meridional) direction

          vy = 0.5D0*(vadyp1_d(j,i,k) - vadym1_d(j,i,k))/ddy
          vyy = (vadyp1_d(j,i,k) - 2.0D0*vadvy_d(j,i,k) + &
                 vadym1_d(j,i,k))/(ddy*ddy)
          ydis = - vadvy_d(j,i,k)*dt + 0.5D0*(dtsq*vadvy_d(j,i,k)*vy) - &
                 (dtcb*vadvy_d(j,i,k))*(vy*vy +vadvy_d(j,i,k)*vyy)/6.0D0
          yn = ydis/ddy
          ynp = int(yn)
          if ( abs(ynp) >= 3 ) then
            write(stderr,*) 'SL Advection problem in SN direction.'
            call fatal(__FILE__,__LINE__,'SLADVECTION')
          end if
          betay = abs((ynp*ddy - ydis)/ddy)
          ysn = int(sign(1.0D0,yn))
          yndp_d(j,i,k) = i + ynp
          ynnm1dp_d(j,i,k) = yndp_d(j,i,k) + ysn
          ynnm2dp_d(j,i,k) = ynnm1dp_d(j,i,k) + ysn
          ynnp1dp_d(j,i,k) = yndp_d(j,i,k) - ysn
          if ( ma%has_bdybottom ) then
            if ( yndp_x(j,i,k) < ide1 ) yndp_x(j,i,k) = ide1
            if ( ynnm1dp_x(j,i,k) < ide1 ) ynnm1dp_x(j,i,k) = ide1
            if ( ynnm2dp_x(j,i,k) < ide1 ) ynnm2dp_x(j,i,k) = ide1
            if ( ynnp1dp_x(j,i,k) < ide1 ) ynnp1dp_x(j,i,k) = ide1
          end if
          if ( ma%has_bdytop ) then
            if ( yndp_x(j,i,k) > ide2 ) yndp_x(j,i,k) = ide2
            if ( ynnm1dp_x(j,i,k) > ide2 ) ynnm1dp_x(j,i,k) = ide2
            if ( ynnm2dp_x(j,i,k) > ide2 ) ynnm2dp_x(j,i,k) = ide2
            if ( ynnp1dp_x(j,i,k) > ide2 ) ynnp1dp_x(j,i,k) = ide2
          end if

          ! to estimate the value at the departure point at t=tau -1
          ! using cubic interpolation

          ! weighting coefficients
          alffbl_d(j,i,k) = alfax
          alfm2dp_d(j,i,k) = -(alfax*(1.0D0 - alfax*alfax))/6.0D0
          alfm1dp_d(j,i,k) = (alfax*(1.0D0 + alfax)*(2.0D0 - alfax))/2.0D0
          alfdp_d(j,i,k) = ((1.0D0 - alfax*alfax)*(2.0D0 - alfax))/2.0D0
          alfp1dp_d(j,i,k) = -(alfax*(1.0D0 - alfax)*(2.0D0 - alfax))/6.0D0
          betm2dp_d(j,i,k) = -(betay*(1.0D0 - betay*betay))/6.0D0
          betm1dp_d(j,i,k) = (betay*(1.0D0 + betay)*(2.0D0 - betay))/2.0D0
          betdp_d(j,i,k) = ((1.0D0 - betay*betay)*(2.0D0 - betay))/2.0D0
          betp1dp_d(j,i,k) = -(betay*(1.0D0 - betay)*(2.0D0 - betay))/6.0D0
        end do
      end do
    end do
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine trajcalc_d

  subroutine slhadv_x3d(ften,var)
    implicit none
    real(rk8) , pointer , intent(inout) , dimension(:,:,:) :: ften
    real(rk8) , pointer , intent(in) , dimension(:,:,:) :: var
    real(rk8) :: tbadp , tbmax , tbmin , tsla , bl1 , bl2 , cb1 , cb2
    integer(ik4) :: i , j , k
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'slhadv_x3d'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    do k = 1, kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          ! GTD bilinear for the two outermost zonal grid points
          bl1 = alffbl_x(j,i,k)*var(xnnm1dp_x(j,i,k),ynnp1dp_x(j,i,k),k) + &
             (d_one - alffbl_x(j,i,k))*var(xndp_x(j,i,k),ynnp1dp_x(j,i,k),k)
          bl2 = alffbl_x(j,i,k)*var(xnnm1dp_x(j,i,k),ynnm2dp_x(j,i,k),k) + &
             (d_one - alffbl_x(j,i,k))*var(xndp_x(j,i,k),ynnm2dp_x(j,i,k),k)
          ! GTD Cubic for the two innermost zonal grid points
          cb1 = alfm2dp_x(j,i,k)*var(xnnm2dp_x(j,i,k),yndp_x(j,i,k),k) + &
                alfm1dp_x(j,i,k)*var(xnnm1dp_x(j,i,k),yndp_x(j,i,k),k) + &
                alfdp_x(j,i,k)*var(xndp_x(j,i,k),yndp_x(j,i,k),k) + &
                alfp1dp_x(j,i,k)*var(xnnp1dp_x(j,i,k),yndp_x(j,i,k),k)
          cb2 = alfm2dp_x(j,i,k)*var(xnnm2dp_x(j,i,k),ynnm1dp_x(j,i,k),k) + &
                alfm1dp_x(j,i,k)*var(xnnm1dp_x(j,i,k),ynnm1dp_x(j,i,k),k) + &
                alfdp_x(j,i,k)*var(xndp_x(j,i,k),ynnm1dp_x(j,i,k),k) + &
                alfp1dp_x(j,i,k)*var(xnnp1dp_x(j,i,k),ynnm1dp_x(j,i,k),k)
          ! GTD finally a Cubic interpolation on
          !              cb1,cb2,bl1 and bl2 pts in y dirn
          tbadp = betm2dp_x(j,i,k)*bl2 + betm1dp_x(j,i,k)*cb2 + &
                  betdp_x(j,i,k)*cb1 + betp1dp_x(j,i,k)*bl1
          ! to get the maximum and minimum value
          tbmax = max(var(xndp_x(j,i,k),yndp_x(j,i,k),k),    &
                      var(xndp_x(j,i,k),ynnm1dp_x(j,i,k),k), &
                      var(xnnm1dp_x(j,i,k),yndp_x(j,i,k),k), &
                      var(xnnm1dp_x(j,i,k),ynnm1dp_x(j,i,k),k))
          tbmin = min(var(xndp_x(j,i,k),yndp_x(j,i,k),k) ,   &
                      var(xndp_x(j,i,k),ynnm1dp_x(j,i,k),k), &
                      var(xnnm1dp_x(j,i,k),yndp_x(j,i,k),k), &
                      var(xnnm1dp_x(j,i,k),ynnm1dp_x(j,i,k),k))
          ! for the quasi monotonic sl
          if ( tbadp > tbmax ) then
            tsla = tbmax
          else if ( tbadp < tbmin ) then
            tsla = tbmin
          else
            tsla = tbadp
          end if
          ! to compute the tendency
          ften(j,i,k) =  (tsla - var(j,i,k))/dt
        end do
      end do
    end do
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine slhadv_x3d

  subroutine slhadv_x4d(ften,var,m,p)
    implicit none
    real(rk8) , pointer , intent(inout) , dimension(:,:,:,:) :: ften
    real(rk8) , pointer , intent(in) , dimension(:,:,:,:) :: var
    integer(ik4) , optional , intent(in) :: m , p
    real(rk8) :: tbadp , tbmax , tbmin , tsla , bl1 , bl2 , cb1 , cb2
    integer(ik4) :: i , j , k , n , n1 , n2
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'slhadv_x4d'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    if ( present(m) ) then
      if ( present(p) ) then
        n1 = m
        n2 = p
      else
        n1 = m
        n2 = m
      end if
    else
      n1 = lbound(var,4)
      n2 = ubound(var,4)
    end if
    do n = n1 , n2
      do k = 1, kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            ! GTD bilinear for the two outermost zonal grid points
            bl1 = alffbl_x(j,i,k)*var(xnnm1dp_x(j,i,k),ynnp1dp_x(j,i,k),k,n) + &
               (d_one - alffbl_x(j,i,k))*var(xndp_x(j,i,k),ynnp1dp_x(j,i,k),k,n)
            bl2 = alffbl_x(j,i,k)*var(xnnm1dp_x(j,i,k),ynnm2dp_x(j,i,k),k,n) + &
               (d_one - alffbl_x(j,i,k))*var(xndp_x(j,i,k),ynnm2dp_x(j,i,k),k,n)
            ! GTD Cubic for the two innermost zonal grid points
            cb1 = alfm2dp_x(j,i,k)*var(xnnm2dp_x(j,i,k),yndp_x(j,i,k),k,n) + &
                  alfm1dp_x(j,i,k)*var(xnnm1dp_x(j,i,k),yndp_x(j,i,k),k,n) + &
                  alfdp_x(j,i,k)*var(xndp_x(j,i,k),yndp_x(j,i,k),k,n) + &
                  alfp1dp_x(j,i,k)*var(xnnp1dp_x(j,i,k),yndp_x(j,i,k),k,n)
            cb2 = alfm2dp_x(j,i,k)*var(xnnm2dp_x(j,i,k),ynnm1dp_x(j,i,k),k,n) +&
                  alfm1dp_x(j,i,k)*var(xnnm1dp_x(j,i,k),ynnm1dp_x(j,i,k),k,n) +&
                  alfdp_x(j,i,k)*var(xndp_x(j,i,k),ynnm1dp_x(j,i,k),k,n) + &
                  alfp1dp_x(j,i,k)*var(xnnp1dp_x(j,i,k),ynnm1dp_x(j,i,k),k,n)
            ! GTD finally a Cubic interpolation on
            !              cb1,cb2,bl1 and bl2 pts in y dirn
            tbadp = betm2dp_x(j,i,k)*bl2 + betm1dp_x(j,i,k)*cb2 + &
                    betdp_x(j,i,k)*cb1 + betp1dp_x(j,i,k)*bl1
            ! to get the maximum and minimum value
            tbmax = max(var(xndp_x(j,i,k),yndp_x(j,i,k),k,n),    &
                        var(xndp_x(j,i,k),ynnm1dp_x(j,i,k),k,n), &
                        var(xnnm1dp_x(j,i,k),yndp_x(j,i,k),k,n), &
                        var(xnnm1dp_x(j,i,k),ynnm1dp_x(j,i,k),k,n))
            tbmin = min(var(xndp_x(j,i,k),yndp_x(j,i,k),k,n) ,   &
                        var(xndp_x(j,i,k),ynnm1dp_x(j,i,k),k,n), &
                        var(xnnm1dp_x(j,i,k),yndp_x(j,i,k),k,n), &
                        var(xnnm1dp_x(j,i,k),ynnm1dp_x(j,i,k),k,n))
            ! for the quasi monotonic sl
            if ( tbadp > tbmax ) then
              tsla = tbmax
            else if ( tbadp < tbmin ) then
              tsla = tbmin
            else
              tsla = tbadp
            end if
            ! to compute the tendency
            ften(j,i,k,n) =  (tsla - var(j,i,k,n))/dt
          end do
        end do
      end do
    end do
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine slhadv_x4d

  subroutine slhadv_d(ften,var)
    implicit none
    real(rk8) , pointer , intent(inout) , dimension(:,:,:) :: ften
    real(rk8) , pointer , intent(in) , dimension(:,:,:) :: var

    real(rk8) :: tbadp, tbmax, tbmin, tsla, bl1, bl2, cb1, cb2
    integer(ik4) :: i , j , k
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'slhadv_d'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    do k = 1, kz
      do i = idi1 , idi2
        do j = jdi1 , jdi2
          ! GTD bilinear for the two outermost zonal grid points
          bl1 = alffbl_d(j,i,k)*var(xnnm1dp_d(j,i,k),ynnp1dp_d(j,i,k),k) + &
             (d_one - alffbl_d(j,i,k))*var(xndp_d(j,i,k),ynnp1dp_d(j,i,k),k)
          bl2 = alffbl_d(j,i,k)*var(xnnm1dp_d(j,i,k),ynnm2dp_d(j,i,k),k) + &
             (d_one - alffbl_d(j,i,k))*var(xndp_d(j,i,k),ynnm2dp_d(j,i,k),k)
          ! GTD Cubic for the two innermost zonal grid points
          cb1 = alfm2dp_d(j,i,k)*var(xnnm2dp_d(j,i,k),yndp_d(j,i,k),k) + &
                alfm1dp_d(j,i,k)*var(xnnm1dp_d(j,i,k),yndp_d(j,i,k),k) + &
                alfdp_d(j,i,k)*var(xndp_d(j,i,k),yndp_d(j,i,k),k) +      &
                alfp1dp_d(j,i,k)*var(xnnp1dp_d(j,i,k),yndp_d(j,i,k),k)
          cb2 = alfm2dp_d(j,i,k)*var(xnnm2dp_d(j,i,k),ynnm1dp_d(j,i,k),k) + &
                alfm1dp_d(j,i,k)*var(xnnm1dp_d(j,i,k),ynnm1dp_d(j,i,k),k) + &
                alfdp_d(j,i,k)*var(xndp_d(j,i,k),ynnm1dp_d(j,i,k),k) +      &
                alfp1dp_d(j,i,k)*var(xnnp1dp_d(j,i,k),ynnm1dp_d(j,i,k),k)
          ! GTD finally a Cubic interpolation on
          !   cb1,cb2,bl1 and bl2 pts in y dirn
          tbadp = betm2dp_d(j,i,k)*bl2 + betm1dp_d(j,i,k)*cb2 + &
                  betdp_d(j,i,k)*cb1 + betp1dp_d(j,i,k)*bl1
          ! to get the maximum and minimum value
          tbmax = max(var(xndp_d(j,i,k),yndp_d(j,i,k),k) ,   &
                      var(xndp_d(j,i,k),ynnm1dp_d(j,i,k),k), &
                      var(xnnm1dp_d(j,i,k),yndp_d(j,i,k),k), &
                      var(xnnm1dp_d(j,i,k),ynnm1dp_d(j,i,k),k))
          tbmin = min(var(xndp_d(j,i,k),yndp_d(j,i,k),k),    &
                      var(xndp_d(j,i,k),ynnm1dp_d(j,i,k),k), &
                      var(xnnm1dp_d(j,i,k),yndp_d(j,i,k),k), &
                      var(xnnm1dp_d(j,i,k),ynnm1dp_d(j,i,k),k))
          ! for the quasi monotonic sl
          if ( tbadp > tbmax ) then
            tsla = tbmax
          else if ( tbadp < tbmin ) then
            tsla = tbmin
          else
            tsla = tbadp
          end if
          ! to compute the tendency
          ften(j,i,k) =  (tsla - var(j,i,k))/dt
        end do
      end do
    end do
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine slhadv_d

  subroutine hdvg_x3d(ften,var)
    implicit none
    real(rk8) , pointer , intent(inout) , dimension(:,:,:) :: ften
    real(rk8) , pointer , intent(in) , dimension(:,:,:) :: var

    real(rk8) :: ucapf_x, ucapi_x , vcapf_x , vcapi_x
    real(rk8) :: ducapdx , dvcapdy , hdvg , tatotdvtrm
    integer(ik4) :: i , j , k
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'hdvg_x3d'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          ucapf_x = (atmx%u(j+1,i+1,k)*mddom%msfd(j+1,i+1) + &
                     atmx%u(j+1,i,k)*mddom%msfd(j+1,i))*d_half
          ucapi_x = (atmx%u(j,i+1,k)*mddom%msfd(j,i+1) + &
                     atmx%u(j,i,k)*mddom%msfd(j,i))*d_half
          ducapdx = (ucapf_x - ucapi_x) / dx
          vcapf_x = (atmx%v(j+1,i+1,k)*mddom%msfd(j+1,i+1) + &
                     atmx%v(j,i+1,k)*mddom%msfd(j,i+1))*d_half
          vcapi_x = (atmx%v(j+1,i,k)*mddom%msfd(j+1,i) + &
                     atmx%v(j,i,k)*mddom%msfd(j,i))*d_half
          dvcapdy = (vcapf_x - vcapi_x) / dx
          hdvg = (ducapdx + dvcapdy)/(mddom%msfx(j,i)*mddom%msfx(j,i))
          tatotdvtrm = var(j,i,k)*hdvg
          ften(j,i,k) = ften(j,i,k) - tatotdvtrm
          !GTD ucapf_x = 0.25D0*(atmx%u(j+1,i,k) + &
          !GTD   atmx%u(j+1,i+1,k) + atmx%u(j+2,i+1,k) + &
          !GTD   atmx%u(j+2,i,k))*mddom%msfx(j+1,i)
          !GTD ucapi_x = 0.25D0*(atmx%u(j,i,k) + &
          !GTD   atmx%u(j,i+1,k) + atmx%u(j-1,i+1,k) + &
          !GTD   atmx%u(j-1,i,k))*mddom%msfx(j-1,i)
          !GTD vcapf_x = 0.25D0*(atmx%v(j,i+1,idyp1_x,k) + &
          !GTD   atmx%v(j+1,i+1,k+ atmx%v(j+1,i+2,k) + &
          !GTD   atmx%v(j,i+2,k))*mddom%msfx(j,i+1)
          !GTD vcapi_x = 0.25D0*(atmx%v(j,i,k) + &
          !GTD   atmx%v(j,i-1,k) + atmx%v(j+1,i-1,k) + &
          !GTD   atmx%v(j+1,i,k))*mddom%msfx(j,i-1)
          !GTD ducapdx = (ucapf_x - ucapi_x) / (2.0D0*dx)
          !GTD dvcapdy = (vcapf_x - vcapi_x) / (2.0D0*dx)
          !GTD hdvg = (ducapdx + dvcapdy)/(mddom%msfx(j,i)*mddom%msfx(j,i))
          !GTD tatotdvtrm = var(j,i,k)*hdvg
          !GTD ften(j,i,k) = ften(j,i,k) - tatotdvtrm
        end do
      end do
    end do
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine hdvg_x3d

  subroutine hdvg_x4d(ften,var,m,p)
    implicit none
    real(rk8) , pointer , intent(inout) , dimension(:,:,:,:) :: ften
    real(rk8) , pointer , intent(in) , dimension(:,:,:,:) :: var
    integer(ik4) , optional , intent(in) :: m , p
    real(rk8) :: ucapf_x, ucapi_x , vcapf_x , vcapi_x
    real(rk8) :: ducapdx , dvcapdy , hdvg , tatotdvtrm
    integer(ik4) :: i , j , k , n , n1 , n2
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'hdvg_x4d'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    if ( present(m) ) then
      if ( present(p) ) then
        n1 = m
        n2 = p
      else
        n1 = m
        n2 = m
      end if
    else
      n1 = lbound(var,4)
      n2 = ubound(var,4)
    end if
    do n = n1 , n2
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            ucapf_x = (atmx%u(j+1,i+1,k)*mddom%msfd(j+1,i+1) + &
                       atmx%u(j+1,i,k)*mddom%msfd(j+1,i))*d_half
            ucapi_x = (atmx%u(j,i+1,k)*mddom%msfd(j,i+1) + &
                       atmx%u(j,i,k)*mddom%msfd(j,i))*d_half
            ducapdx = (ucapf_x - ucapi_x) / dx
            vcapf_x = (atmx%v(j+1,i+1,k)*mddom%msfd(j+1,i+1) + &
                       atmx%v(j,i+1,k)*mddom%msfd(j,i+1))*d_half
            vcapi_x = (atmx%v(j+1,i,k)*mddom%msfd(j+1,i) + &
                       atmx%v(j,i,k)*mddom%msfd(j,i))*d_half
            dvcapdy = (vcapf_x - vcapi_x) / dx
            hdvg = (ducapdx + dvcapdy)/(mddom%msfx(j,i)*mddom%msfx(j,i))
            tatotdvtrm = var(j,i,k,n)*hdvg
            ften(j,i,k,n) = ften(j,i,k,n) - tatotdvtrm
            !GTD ucapf_x = 0.25D0*(atmx%u(j+1,i,k) + &
            !GTD   atmx%u(j+1,i+1,k) + atmx%u(j+2,i+1,k) + &
            !GTD   atmx%u(j+2,i,k))*mddom%msfx(j+1,i)
            !GTD ucapi_x = 0.25D0*(atmx%u(j,i,k) + &
            !GTD   atmx%u(j,i+1,k) + atmx%u(j-1,i+1,k) + &
            !GTD   atmx%u(j-1,i,k))*mddom%msfx(j-1,i)
            !GTD vcapf_x = 0.25D0*(atmx%v(j,i+1,idyp1_x,k) + &
            !GTD   atmx%v(j+1,i+1,k+ atmx%v(j+1,i+2,k) + &
            !GTD   atmx%v(j,i+2,k))*mddom%msfx(j,i+1)
            !GTD vcapi_x = 0.25D0*(atmx%v(j,i,k) + &
            !GTD   atmx%v(j,i-1,k) + atmx%v(j+1,i-1,k) + &
            !GTD   atmx%v(j+1,i,k))*mddom%msfx(j,i-1)
            !GTD ducapdx = (ucapf_x - ucapi_x) / (2.0D0*dx)
            !GTD dvcapdy = (vcapf_x - vcapi_x) / (2.0D0*dx)
            !GTD hdvg = (ducapdx + dvcapdy)/(mddom%msfx(j,i)*mddom%msfx(j,i))
            !GTD tatotdvtrm = var(j,i,k,n)*hdvg
            !GTD ften(j,i,k,n) = ften(j,i,k,n) - tatotdvtrm
          end do
        end do
      end do
    end do
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine hdvg_x4d

  subroutine hdvg_d(ften,var,dx)
    implicit none
    real(rk8) , intent(in) :: dx
    real(rk8) , pointer , intent(inout) , dimension(:,:,:) :: ften
    real(rk8) , pointer , intent(in) , dimension(:,:,:) :: var
    real(rk8) :: ucapf, ucapi , ducapdx , vcapf , vcapi ,   &
                 dvcapdy , hdvg , tatotdvtrm
    integer(ik4) :: i , j , k
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'hdvg_x'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    do k = 1 , kz
      do i = idi1 , idi2
        do j = jdi1 , jdi2
          ucapf = atmx%u(j+1,i,k)*mddom%msfd(j+1,i)
          ucapi = atmx%u(j-1,i,k)*mddom%msfd(j-1,i)
          ducapdx = (ucapf - ucapi) / (2.0D0*dx)
          vcapf = atmx%v(j,i+1,k)*mddom%msfd(j,i+1)
          vcapi = atmx%v(j,i-1,k)*mddom%msfd(j,i-1)
          dvcapdy = (vcapf - vcapi) / (2.0D0*dx)
          hdvg = (ducapdx + dvcapdy) / (mddom%msfd(j,i)*mddom%msfd(j,i))
          tatotdvtrm = var(j,i,k)*hdvg/mddom%msfd(j,i)
          ften(j,i,k) = ften(j,i,k) - tatotdvtrm
        end do
      end do
    end do
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine hdvg_d

end module mod_sladvection
