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

module mod_diffusion
!
! Diffusion calculations
!
  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_memutil
  use mod_dynparam
  use mod_runparams
  use mod_mppparam
  use mod_service

  implicit none

  private

  real(rkx) , pointer , dimension(:,:,:) :: xkc , xkd , xkcf
  real(rkx) , public , pointer , dimension(:,:) :: hgfact

  real(rkx) :: dydc , xkhmax

  interface diffu_x
    module procedure diffu_x3df
    module procedure diffu_x3d
    module procedure diffu_x4d
    module procedure diffu_x4d3d
  end interface

  public :: allocate_mod_diffusion
  public :: initialize_diffusion
  public :: calc_coeff
  public :: diffu_d
  public :: diffu_x

  real(rkx) , pointer , dimension(:,:,:) :: ud , vd , wx
  real(rkx) , pointer , dimension(:,:) :: pc , pd

  !
  ! Use 9-point laplacian as in LeVeque,
  !   Finite Difference Methods for Differential Equations , Eq. 3.17
  !
  real(rkx) , parameter :: o4_c1 =   4.0_rkx/6.0_rkx
  real(rkx) , parameter :: o4_c2 =   1.0_rkx/6.0_rkx
  real(rkx) , parameter :: o4_c3 = -20.0_rkx/6.0_rkx

  contains

  subroutine allocate_mod_diffusion
    implicit none
    call getmem2d(hgfact,jce1ga,jce2ga,ice1ga,ice2ga,'diffusion:hgfact')
    call getmem3d(xkc,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,'diffusion:xkc')
    call getmem3d(xkd,jdi1,jdi2,idi1,idi2,1,kz,'diffusion:xkd')
    call getmem3d(xkcf,jci1,jci2,ici1,ici2,1,kzp1,'diffusion:xkcf')
  end subroutine allocate_mod_diffusion

  subroutine initialize_diffusion
    use mod_atm_interface , only : mddom , sfs , atms
    implicit none
    integer(ik4) :: i , j
    real(rkx) :: hg1 , hg2 , hg3 , hg4
    real(rkx) :: hgmax , xkhz , minxkh , maxxkh
    !
    ! Diffusion coefficients: for non-hydrostatic, follow the MM5
    ! The hydrostatic diffusion is following the RegCM3 formulation
    !
    xkhmax = dxsq/(64.0_rkx*dtsec)  ! Computation stability
    dydc = adyndif*vonkar*vonkar*dx*d_rfour ! Deformation term coefficient
    if ( idynamic == 1 ) then
      xkhz = ckh * 1.5e-3_rkx*dxsq/dtsec
    else
      ! (Xu et al., MWR, 2001, 502-516)
      xkhz = ckh * dx
      ! xkhz = ckh * 1.5e-3_rkx*dxsq/dtsec
      ! xkhz = ckh * 3.0e-3_rkx*dxsq/dtsec
      xkhmax = d_two*xkhmax
    end if
    if ( myid == 0 ) then
      write(stdout,'(a,e13.6,a)') &
        ' Constant hor. diff. coef. = ',xkhz,' m^2 s-1'
      write(stdout,'(a,e13.6,a)') &
        ' Maximum  hor. diff. coef. = ',xkhmax,' m^2 s-1'
    end if
    !
    ! Calculate topographical correction to diffusion coefficient
    !
    do i = ice1ga , ice2ga
      do j = jce1ga , jce2ga
        hgfact(j,i) = xkhz
      end do
    end do
    if ( diffu_hgtf == 1 ) then
      ! Should we have a vertical profile for this?
      do concurrent ( j = jci1ga:jci2ga , i = ici1ga:ici2ga )
        hg1 = abs((mddom%ht(j,i)-mddom%ht(j,i-1))/dx)
        hg2 = abs((mddom%ht(j,i)-mddom%ht(j,i+1))/dx)
        hg3 = abs((mddom%ht(j,i)-mddom%ht(j-1,i))/dx)
        hg4 = abs((mddom%ht(j,i)-mddom%ht(j+1,i))/dx)
        hgmax = max(hg1,hg2,hg3,hg4)*regrav*1.0e3_rkx
        hgfact(j,i) = xkhz/(d_one+hgmax**2)
      end do
      call maxall(maxval(hgfact),maxxkh)
      call minall(minval(hgfact),minxkh)
      if ( myid == 0 ) then
        write(stdout,'(a,e13.6,a,e13.6,a)') &
          ' Computed hor. diff. coef. range = ',minxkh,' -',maxxkh,' m^2 s-1'
      end if
    end if
    call assignpnt(sfs%psb,pc)
    call assignpnt(sfs%psdotb,pd)
    call assignpnt(atms%ubd3d,ud)
    call assignpnt(atms%vbd3d,vd)
    call assignpnt(atms%wb3d,wx)
  end subroutine initialize_diffusion

  subroutine calc_coeff
    implicit none
    real(rkx) :: dudx , dvdx , dudy , dvdy , dwdz , duv
    integer(ik4) :: i , j , k

    xkc(:,:,:)  = d_zero
    xkd(:,:,:)  = d_zero
    xkcf(:,:,:) = d_zero
    !
    ! compute the horizontal diffusion coefficient and stored in xkc:
    ! the values are calculated at cross points, but they also used
    ! for dot-point variables.
    !
    if ( idynamic == 1 ) then
      do concurrent ( j = jce1:jce2 , i = ice1:ice2 , k = 1:kz )
        ! Following Smagorinsky et al, 1965 for eddy viscosity
        dudx = ud(j+1,i,k) + ud(j+1,i+1,k) - &
               ud(j,i,k)   - ud(j,i+1,k)
        dvdx = vd(j+1,i,k) + vd(j+1,i+1,k) - &
               vd(j,i,k)   - vd(j,i+1,k)
        dudy = ud(j,i+1,k) + ud(j+1,i+1,k) - &
               ud(j,i,k)   - ud(j+1,i,k)
        dvdy = vd(j,i+1,k) + vd(j+1,i+1,k) - &
               vd(j,i,k)   - vd(j+1,i,k)
        duv = sqrt((dudx-dvdy)*(dudx-dvdy)+(dvdx+dudy)*(dvdx+dudy))
        xkc(j,i,k) = min((hgfact(j,i) + dydc*duv),xkhmax)
      end do
    else
      do concurrent ( j = jce1:jce2 , i = ice1:ice2 , k = 1:kz )
        ! Following Smagorinsky et al, 1965 for eddy viscosity
        dudx = ud(j+1,i,k) + ud(j+1,i+1,k) - &
               ud(j,i,k)   - ud(j,i+1,k)
        dvdx = vd(j+1,i,k) + vd(j+1,i+1,k) - &
               vd(j,i,k)   - vd(j,i+1,k)
        dudy = ud(j,i+1,k) + ud(j+1,i+1,k) - &
               ud(j,i,k)   - ud(j+1,i,k)
        dvdy = vd(j,i+1,k) + vd(j+1,i+1,k) - &
               vd(j,i,k)   - vd(j+1,i,k)
        dwdz = wx(j,i,k) - wx(j,i,k+1)
        duv = sqrt(max((dudx-dvdy)*(dudx-dvdy) + &
                       (dvdx+dudy)*(dvdx+dudy) - dwdz*dwdz,d_zero))
        xkc(j,i,k) = min((hgfact(j,i) + dydc*duv),xkhmax)
      end do
    end if
    xkcf(:,:,1) = xkc(jci1:jci2,ici1:ici2,1)
    do concurrent ( j = jci1:jci2 , i = ici1:ici2 , k = 1:kz )
      xkcf(j,i,k+1) = xkc(j,i,k)
    end do
    call exchange(xkc,1,jce1,jce2,ice1,ice2,1,kz)
    do concurrent ( j = jdi1:jdi2 , i = idi1:idi2 , k = 1:kz )
      xkd(j,i,k) = min(d_rfour*(xkc(j,i,k)+xkc(j-1,i-1,k) + &
                                xkc(j-1,i,k)+xkc(j,i-1,k)),xkhmax)
    end do
    do concurrent ( j = jci1:jci2 , i = ici1:ici2 , k = 1:kz )
      xkc(j,i,k) = xkc(j,i,k) * rdxsq * pc(j,i)
    end do
    do concurrent ( j = jci1:jci2 , i = ici1:ici2 , k = 1:kzp1 )
      xkcf(j,i,k) = xkcf(j,i,k) * rdxsq * pc(j,i)
    end do
    do concurrent ( j = jdi1:jdi2 , i = idi1:idi2 , k = 1:kz )
      xkd(j,i,k) = xkd(j,i,k) * rdxsq * pd(j,i)
    end do
  end subroutine calc_coeff
  !
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !                                                                     c
  !     These subroutines computes the diffusion term for decoupled     c
  !     variable on constant sigma surface.                             c
  !                                                                     c
  !     ften    : tendency for variable                                 c
  !                                                                     c
  !     xk      : horizontal diffusion coefficient                      c
  !                                                                     c
  !     f       : variable f on cross/dot points                        c
  !                                                                     c
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  subroutine diffu_d(uten,vten,u,v)
    implicit none
    real(rkx) , pointer , dimension(:,:,:) , intent(in) :: u , v
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: uten , vten
    integer(ik4) :: i , j , k
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'diffu_d'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    !
    ! fourth-order scheme
    !
    do concurrent ( j = jdi1:jdi2 , i = idi1:idi2 , k = 1:kz )
      uten(j,i,k) = uten(j,i,k) + xkd(j,i,k) * &
             (o4_c1*(u(j+1,i,k)+u(j-1,i,k)+u(j,i+1,k)+u(j,i-1,k)) + &
              o4_c2*(u(j+1,i+1,k)+u(j-1,i-1,k)+u(j-1,i+1,k)+u(j+1,i-1,k)) + &
              o4_c3*(u(j,i,k)))
      vten(j,i,k) = vten(j,i,k) + xkd(j,i,k) * &
              (o4_c1*(v(j+1,i,k)+v(j-1,i,k)+v(j,i+1,k)+v(j,i-1,k)) + &
               o4_c2*(v(j+1,i+1,k)+v(j-1,i-1,k)+v(j-1,i+1,k)+v(j+1,i-1,k)) + &
               o4_c3*(v(j,i,k)))
    end do
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine diffu_d

  subroutine diffu_x3df(ften,f,fac)
    implicit none
    real(rkx) , pointer , dimension(:,:,:) , intent(in) :: f
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: ften
    real(rkx) , intent(in) :: fac
    integer(ik4) :: i , j , k
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'diffu_x3df'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    !
    ! fourth-order scheme
    !
    do concurrent ( j = jci1:jci2 , i = ici1:ici2 , k = 1:kzp1 )
      ften(j,i,k) = ften(j,i,k) + fac * xkcf(j,i,k) * &
           (o4_c1*(f(j+1,i,k)+f(j-1,i,k)+f(j,i+1,k)+f(j,i-1,k)) +   &
            o4_c2*(f(j+1,i+1,k)+f(j-1,i-1,k)+f(j-1,i+1,k)+f(j+1,i-1,k)) + &
            o4_c3*(f(j,i,k)))
    end do
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine diffu_x3df

  subroutine diffu_x3d(ften,f)
    implicit none
    real(rkx) , pointer , dimension(:,:,:) , intent(in) :: f
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: ften
    integer(ik4) :: i , j , k
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'diffu_x3d'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    !
    ! fourth-order scheme
    !
    do concurrent ( j = jci1:jci2 , i = ici1:ici2 , k = 1:kz )
      ften(j,i,k) = ften(j,i,k) + xkc(j,i,k) * &
           (o4_c1*(f(j+1,i,k)+f(j-1,i,k)+f(j,i+1,k)+f(j,i-1,k)) +   &
            o4_c2*(f(j+1,i+1,k)+f(j-1,i-1,k)+f(j-1,i+1,k)+f(j+1,i-1,k)) + &
            o4_c3*(f(j,i,k)))
    end do
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine diffu_x3d

  subroutine diffu_x4d(ften,f,n1,n2,fac)
    implicit none
    integer(ik4) , intent(in) :: n1 , n2
    real(rkx) , pointer , dimension(:,:,:,:) , intent(in) :: f
    real(rkx) , pointer , dimension(:,:,:,:) , intent(inout) :: ften
    real(rkx) , intent(in) :: fac
    integer(ik4) ::  n
    !
    do n = n1 , n2
      call diffu_x4d3d(ften,f,n,fac)
    end do
  end subroutine diffu_x4d

  subroutine diffu_x4d3d(ften,f,n,fac)
    implicit none
    integer(ik4) , intent(in) :: n
    real(rkx) , pointer , dimension(:,:,:,:) , intent(in) :: f
    real(rkx) , pointer , dimension(:,:,:,:) , intent(inout) :: ften
    real(rkx) , intent(in) :: fac
    integer(ik4) :: i , j , k
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'diffu_x4d3d'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    !
    ! fourth-order scheme
    !
    do concurrent ( j = jci1:jci2 , i = ici1:ici2 , k = 1:kz )
      ften(j,i,k,n) = ften(j,i,k,n) + fac * xkc(j,i,k) * &
           (o4_c1*(f(j+1,i,k,n)+f(j-1,i,k,n) + &
                   f(j,i+1,k,n)+f(j,i-1,k,n)) +   &
            o4_c2*(f(j+1,i+1,k,n)+f(j-1,i-1,k,n) + &
                   f(j-1,i+1,k,n)+f(j+1,i-1,k,n)) + &
            o4_c3*f(j,i,k,n))
    end do
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine diffu_x4d3d

end module mod_diffusion
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
