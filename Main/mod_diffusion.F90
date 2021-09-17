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
  use mod_stdio
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
  real(rkx) , pointer , dimension(:,:) :: pc , pd , mpd
  !
  ! Use 9-point laplacian as in LeVeque,
  !   Finite Difference Methods for Differential Equations , Eq. 3.17
  !
  real(rkx) , parameter :: o4_c1 =   4.0_rkx/6.0_rkx
  real(rkx) , parameter :: o4_c2 =   1.0_rkx/6.0_rkx
  real(rkx) , parameter :: o4_c3 = -20.0_rkx/6.0_rkx
  !
  ! Weights
  !
  real(rkx) , parameter :: z4_c1 =   1.0_rkx
  real(rkx) , parameter :: z4_c2 =  -4.0_rkx
  real(rkx) , parameter :: z4_c3 =  12.0_rkx
  !
  ! Weights 6th order
  !
  real(rkx) , parameter :: h4_c1 =  10.0_rkx
  real(rkx) , parameter :: h4_c2 =  -5.0_rkx
  real(rkx) , parameter :: h4_c3 =   1.0_rkx
  real(rkx) , parameter :: diff_6th_factor = 0.12_rkx
  real(rkx) :: diff_6th_coef

  contains

  subroutine allocate_mod_diffusion
    implicit none
    if ( idiffu < 3 ) then
      call getmem2d(hgfact,jce1ga,jce2ga,ice1ga,ice2ga,'diffusion:hgfact')
    end if
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
    if ( idiffu < 3 ) then
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
        do i = ici1ga , ici2ga
          do j = jci1ga , jci2ga
            hg1 = abs((mddom%ht(j,i)-mddom%ht(j,i-1))/dx)
            hg2 = abs((mddom%ht(j,i)-mddom%ht(j,i+1))/dx)
            hg3 = abs((mddom%ht(j,i)-mddom%ht(j-1,i))/dx)
            hg4 = abs((mddom%ht(j,i)-mddom%ht(j+1,i))/dx)
            hgmax = max(hg1,hg2,hg3,hg4)*regrav*1.0e3_rkx
            hgfact(j,i) = xkhz/(d_one+hgmax**2)
          end do
        end do
        call maxall(maxval(hgfact),maxxkh)
        call minall(minval(hgfact),minxkh)
        if ( myid == 0 ) then
          write(stdout,'(a,e13.6,a,e13.6,a)') &
            ' Computed background diff. coeff. = ',minxkh,' -',maxxkh,' m^2 s-1'
        end if
      else
        if ( myid == 0 ) then
          write(stdout,'(a,e13.6,a)') &
            ' Computed background diff. coeff. = ', xkhz, ' m^2 s-1'
        end if
      end if
    else
      diff_6th_coef = diff_6th_factor * 0.015625_rkx / ( 2.0_rkx * dtsec )
      if ( myid == 0 ) then
        write(stdout,'(a,e13.6,a)') &
            ' Constant hor. diff. coef. = ',diff_6th_coef,' m^2 s-1'
      end if
    end if

    call assignpnt(mddom%msfd,mpd)
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

    if ( idiffu == 3 ) then
      do concurrent ( j = jci1:jci2 , i = ici1:ici2 , k = 1:kz )
        xkc(j,i,k) = diff_6th_coef * pc(j,i)
      end do
      do concurrent ( j = jci1:jci2 , i = ici1:ici2 , k = 1:kzp1 )
        xkcf(j,i,k) = diff_6th_coef * pc(j,i)
      end do
      do concurrent ( j = jdi1:jdi2 , i = idi1:idi2 , k = 1:kz )
        xkd(j,i,k) = diff_6th_coef * pd(j,i)
      end do
    else
      xkc(:,:,:)  = d_zero
      xkd(:,:,:)  = d_zero
      xkcf(:,:,:) = d_zero
      !
      ! compute the horizontal diffusion coefficient and stored in xkc:
      ! the values are calculated at cross points, but they also used
      ! for dot-point variables.
      !
      if ( idynamic == 1 ) then
        do k = 1 , kz
          do i = ice1 , ice2
            do j = jce1 , jce2
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
          end do
        end do
      else
        do k = 1 , kz
          do i = ice1 , ice2
            do j = jce1 , jce2
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
          end do
        end do
      end if
      xkcf(:,:,1) = xkc(jci1:jci2,ici1:ici2,1)
      do concurrent ( j = jci1:jci2 , i = ici1:ici2 , k = 1:kz )
        xkcf(j,i,k+1) = xkc(j,i,k)
      end do
      call exchange(xkc,1,jce1,jce2,ice1,ice2,1,kz)
      do concurrent ( j = jdi1:jdi2 , i = idi1:idi2 , k = 1:kz )
        xkd(j,i,k) = d_rfour*(xkc(j,i,k)+xkc(j-1,i-1,k) + &
                              xkc(j-1,i,k)+xkc(j,i-1,k))
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
    end if
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
    integer(ik4) :: im1 , im2 , im3
    integer(ik4) :: jm1 , jm2 , jm3
    integer(ik4) :: ip1 , ip2 , ip3
    integer(ik4) :: jp1 , jp2 , jp3
    real(rkx) :: dflux_x_p0 , dflux_x_p1 , dflux_y_p0 , dflux_y_p1
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'diffu_d'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    if ( idiffu == 1 ) then
      !
      ! fourth-order scheme for interior:
      !
      do concurrent ( j = jdii1:jdii2 , i = idii1:idii2 , k = 1:kz )
        uten(j,i,k) = uten(j,i,k) - xkd(j,i,k) * &
               (z4_c1*(u(j+2,i,k)/mpd(j+2,i) +   &
                       u(j-2,i,k)/mpd(j-2,i) +   &
                       u(j,i+2,k)/mpd(j,i+2) +   &
                       u(j,i-2,k)/mpd(j,i-2)) +  &
                z4_c2*(u(j+1,i,k)/mpd(j+1,i) +   &
                       u(j-1,i,k)/mpd(j-1,i) +   &
                       u(j,i+1,k)/mpd(j,i+1) +   &
                       u(j,i-1,k)/mpd(j,i-1)) +  &
                z4_c3*(u(j,i,k)/mpd(j,i)))
        vten(j,i,k) = vten(j,i,k) - xkd(j,i,k) * &
               (z4_c1*(v(j+2,i,k)/mpd(j+2,i) +   &
                       v(j-2,i,k)/mpd(j-2,i) +   &
                       v(j,i+2,k)/mpd(j,i+2) +   &
                       v(j,i-2,k)/mpd(j,i-2)) +  &
                z4_c2*(v(j+1,i,k)/mpd(j+1,i) +   &
                       v(j-1,i,k)/mpd(j-1,i) +   &
                       v(j,i+1,k)/mpd(j,i+1) +   &
                       v(j,i-1,k)/mpd(j,i-1)) +  &
                z4_c3*(v(j,i,k)/mpd(j,i)))
      end do
      !
      ! second-order scheme for boundaries:
      !
      if ( ma%has_bdyleft ) then
        j = jdi1
        do k = 1 , kz
          do i = idi1 , idi2
            uten(j,i,k) = uten(j,i,k) + xkd(j,i,k) * &
              (z4_c1*(u(j+1,i,k)/mpd(j+1,i) +  &
                      u(j-1,i,k)/mpd(j-1,i) +  &
                      u(j,i+1,k)/mpd(j,i+1) +  &
                      u(j,i-1,k)/mpd(j,i-1)) + &
               z4_c2*(u(j,i,k)/mpd(j,i)))
            vten(j,i,k) = vten(j,i,k) + xkd(j,i,k) * &
              (z4_c1*(v(j+1,i,k)/mpd(j+1,i) +  &
                      v(j-1,i,k)/mpd(j-1,i) +  &
                      v(j,i+1,k)/mpd(j,i+1) +  &
                      v(j,i-1,k)/mpd(j,i-1)) + &
               z4_c2*(v(j,i,k)/mpd(j,i)))
          end do
        end do
      end if
      if ( ma%has_bdyright ) then
        j = jdi2
        do k = 1 , kz
          do i = idi1 , idi2
            uten(j,i,k) = uten(j,i,k) + xkd(j,i,k) * &
              (z4_c1*(u(j+1,i,k)/mpd(j+1,i) +  &
                      u(j-1,i,k)/mpd(j-1,i) +  &
                      u(j,i+1,k)/mpd(j,i+1) +  &
                      u(j,i-1,k)/mpd(j,i-1)) + &
               z4_c2*(u(j,i,k)/mpd(j,i)))
            vten(j,i,k) = vten(j,i,k) + xkd(j,i,k) * &
              (z4_c1*(v(j+1,i,k)/mpd(j+1,i) +  &
                      v(j-1,i,k)/mpd(j-1,i) +  &
                      v(j,i+1,k)/mpd(j,i+1) +  &
                      v(j,i-1,k)/mpd(j,i-1)) + &
               z4_c2*(v(j,i,k)/mpd(j,i)))
          end do
        end do
      end if
      if ( ma%has_bdybottom ) then
        i = idi1
        do k = 1 , kz
          do j = jdi1 , jdi2
            uten(j,i,k) = uten(j,i,k) + xkd(j,i,k) * &
              (z4_c1*(u(j+1,i,k)/mpd(j+1,i) +  &
                      u(j-1,i,k)/mpd(j-1,i) +  &
                      u(j,i+1,k)/mpd(j,i+1) +  &
                      u(j,i-1,k)/mpd(j,i-1)) + &
               z4_c2*(u(j,i,k)/mpd(j,i)))
            vten(j,i,k) = vten(j,i,k) + xkd(j,i,k) * &
              (z4_c1*(v(j+1,i,k)/mpd(j+1,i) +  &
                      v(j-1,i,k)/mpd(j-1,i) +  &
                      v(j,i+1,k)/mpd(j,i+1) +  &
                      v(j,i-1,k)/mpd(j,i-1)) + &
               z4_c2*(v(j,i,k)/mpd(j,i)))
          end do
        end do
      end if
      if ( ma%has_bdytop ) then
        i = idi2
        do k = 1 , kz
          do j = jdi1 , jdi2
            uten(j,i,k) = uten(j,i,k) + xkd(j,i,k) * &
              (z4_c1*(u(j+1,i,k)/mpd(j+1,i) +  &
                      u(j-1,i,k)/mpd(j-1,i) +  &
                      u(j,i+1,k)/mpd(j,i+1) +  &
                      u(j,i-1,k)/mpd(j,i-1)) + &
               z4_c2*(u(j,i,k)/mpd(j,i)))
            vten(j,i,k) = vten(j,i,k) + xkd(j,i,k) * &
              (z4_c1*(v(j+1,i,k)/mpd(j+1,i) +  &
                      v(j-1,i,k)/mpd(j-1,i) +  &
                      v(j,i+1,k)/mpd(j,i+1) +  &
                      v(j,i-1,k)/mpd(j,i-1)) + &
               z4_c2*(v(j,i,k)/mpd(j,i)))
          end do
        end do
      end if
    else if ( idiffu == 2 ) then
      !
      ! fourth-order scheme
      !
      do concurrent ( j = jdi1:jdi2 , i = idi1:idi2 , k = 1:kz )
        uten(j,i,k) = uten(j,i,k) + xkd(j,i,k) *    &
               (o4_c1*(u(j+1,i,k)/mpd(j+1,i) +      &
                       u(j-1,i,k)/mpd(j-1,i) +      &
                       u(j,i+1,k)/mpd(j,i+1) +      &
                       u(j,i-1,k)/mpd(j,i-1)) +     &
                o4_c2*(u(j+1,i+1,k)/mpd(j+1,i+1) +  &
                       u(j-1,i-1,k)/mpd(j-1,i-1) +  &
                       u(j-1,i+1,k)/mpd(j-1,i+1) +  &
                       u(j+1,i-1,k)/mpd(j+1,i-1)) + &
                o4_c3*(u(j,i,k)/mpd(j,i)))
        vten(j,i,k) = vten(j,i,k) + xkd(j,i,k) * &
               (o4_c1*(v(j+1,i,k)/mpd(j+1,i) +      &
                       v(j-1,i,k)/mpd(j-1,i) +      &
                       v(j,i+1,k)/mpd(j,i+1) +      &
                       v(j,i-1,k)/mpd(j,i-1)) +     &
                o4_c2*(v(j+1,i+1,k)/mpd(j+1,i+1) +  &
                       v(j-1,i-1,k)/mpd(j-1,i-1) +  &
                       v(j-1,i+1,k)/mpd(j-1,i+1) +  &
                       v(j+1,i-1,k)/mpd(j+1,i-1)) + &
                o4_c3*(v(j,i,k)/mpd(j,i)))
      end do
    else if ( idiffu == 3 ) then
      do k = 1 , kz
        do i = idi1 , idi2
          im1 = max(i-1,1)
          im2 = max(i-2,1)
          im3 = max(i-3,1)
          ip1 = min(i+1,iy)
          ip2 = min(i+2,iy)
          ip3 = min(i+3,iy)
          do j = jdi2 , jdi2
            jm1 = max(j-1,1)
            jm2 = max(j-2,1)
            jm3 = max(j-3,1)
            jp1 = min(j+1,jx)
            jp2 = min(j+2,jx)
            jp3 = min(j+3,jx)
            dflux_x_p0 = ( h4_c1 * (u(j,i,k)/mpd(j,i)     -   &
                                    u(jm1,i,k)/mpd(jm1,i) ) + &
                           h4_c2 * (u(jp1,i,k)/mpd(jp1,i) -   &
                                    u(jm2,i,k)/mpd(jm2,i) ) + &
                           h4_c3 * (u(jp2,i,k)/mpd(jp2,i) -   &
                                    u(jm3,i,k)/mpd(jm3,i) ) )
            if  ( dflux_x_p0 * (u(j,i,k)/mpd(j,i) - &
                                u(jm1,i,k)/mpd(jm1,i)) <= d_zero ) then
              dflux_x_p0 = d_zero
            end if
            dflux_x_p1 = ( h4_c1 * (u(jp1,i,k)/mpd(jp1,i) -   &
                                    u(j,i,k)/mpd(j,i) )     + &
                           h4_c2 * (u(jp2,i,k)/mpd(jp2,i) -   &
                                    u(jm1,i,k)/mpd(jm1,i) ) + &
                           h4_c3 * (u(jp3,i,k)/mpd(jp3,i) -   &
                                    u(jm2,i,k)/mpd(jm2,i) ) )
            if  ( dflux_x_p1 * (u(jp1,i,k)/mpd(jp1,i) - &
                                u(j,i,k)/mpd(j,i)) <= d_zero ) then
              dflux_x_p1 = d_zero
            end if
            dflux_y_p0 = ( h4_c1 * (u(j,i,k)/mpd(j,i)     -   &
                                    u(j,im1,k)/mpd(j,im1) ) + &
                           h4_c2 * (u(j,ip1,k)/mpd(j,ip1) -   &
                                    u(j,im2,k)/mpd(j,im2) ) + &
                           h4_c3 * (u(j,ip2,k)/mpd(j,ip2) -   &
                                    u(j,im3,k)/mpd(j,im3) ) )
            if  ( dflux_y_p0 * (u(j,i,k)/mpd(j,i) - &
                                u(j,im1,k)/mpd(j,im1)) <= d_zero ) then
              dflux_y_p0 = d_zero
            end if
            dflux_y_p1 = ( h4_c1 * (u(j,ip1,k)/mpd(j,ip1) -   &
                                    u(j,i,k)/mpd(j,i) )     + &
                           h4_c2 * (u(j,ip2,k)/mpd(j,ip2) -   &
                                    u(j,im1,k)/mpd(j,im1) ) + &
                           h4_c3 * (u(j,ip3,k)/mpd(j,ip3) -   &
                                    u(j,im2,k)/mpd(j,im2) ) )
            if  ( dflux_y_p1 * (u(j,ip1,k)/mpd(j,ip1) - &
                                u(j,i,k)/mpd(j,i)) <= d_zero ) then
              dflux_y_p1 = d_zero
            end if
            uten(j,i,k) = uten(j,i,k) + xkd(j,i,k) * &
                     ( ( dflux_x_p1 - dflux_x_p0 ) + &
                       ( dflux_y_p1 - dflux_y_p0 ) )
            dflux_x_p0 = ( h4_c1 * (v(j,i,k)/mpd(j,i)     -   &
                                    v(jm1,i,k)/mpd(jm1,i) ) + &
                           h4_c2 * (v(jp1,i,k)/mpd(jp1,i) -   &
                                    v(jm2,i,k)/mpd(jm2,i) ) + &
                           h4_c3 * (v(jp2,i,k)/mpd(jp2,i) -   &
                                    v(jm3,i,k)/mpd(jm3,i) ) )
            if  ( dflux_x_p0 * (v(j,i,k)/mpd(j,i) - &
                                v(jm1,i,k)/mpd(jm1,i)) <= d_zero ) then
              dflux_x_p0 = d_zero
            end if
            dflux_x_p1 = ( h4_c1 * (v(jp1,i,k)/mpd(jp1,i) -   &
                                    v(j,i,k)/mpd(j,i) )     + &
                           h4_c2 * (v(jp2,i,k)/mpd(jp2,i) -   &
                                    v(jm1,i,k)/mpd(jm1,i) ) + &
                           h4_c3 * (v(jp3,i,k)/mpd(jp3,i) -   &
                                    v(jm2,i,k)/mpd(jm2,i) ) )
            if  ( dflux_x_p1 * (v(jp1,i,k)/mpd(jp1,i) - &
                                v(j,i,k)/mpd(j,i)) <= d_zero ) then
              dflux_x_p1 = d_zero
            end if
            dflux_y_p0 = ( h4_c1 * (v(j,i,k)/mpd(j,i)     -   &
                                    v(j,im1,k)/mpd(j,im1) ) + &
                           h4_c2 * (v(j,ip1,k)/mpd(j,ip1) -   &
                                    v(j,im2,k)/mpd(j,im2) ) + &
                           h4_c3 * (v(j,ip2,k)/mpd(j,ip2) -   &
                                    v(j,im3,k)/mpd(j,im3) ) )
            if  ( dflux_y_p0 * (v(j,i,k)/mpd(j,i) - &
                                v(j,im1,k)/mpd(j,im1)) <= d_zero ) then
              dflux_y_p0 = d_zero
            end if
            dflux_y_p1 = ( h4_c1 * (v(j,ip1,k)/mpd(j,ip1) -   &
                                    v(j,i,k)/mpd(j,i) )     + &
                           h4_c2 * (v(j,ip2,k)/mpd(j,ip2) -   &
                                    v(j,im1,k)/mpd(j,im1) ) + &
                           h4_c3 * (v(j,ip3,k)/mpd(j,ip3) -   &
                                    v(j,im2,k)/mpd(j,im2) ) )
            if  ( dflux_y_p1 * (v(j,ip1,k)/mpd(j,ip1) - &
                                v(j,i,k)/mpd(j,i)) <= d_zero ) then
              dflux_y_p1 = d_zero
            end if
            vten(j,i,k) = vten(j,i,k) + xkd(j,i,k) * &
                     ( ( dflux_x_p1 - dflux_x_p0 ) + &
                       ( dflux_y_p1 - dflux_y_p0 ) )
          end do
        end do
      end do
    end if
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
    integer(ik4) :: im1 , im2 , im3
    integer(ik4) :: jm1 , jm2 , jm3
    integer(ik4) :: ip1 , ip2 , ip3
    integer(ik4) :: jp1 , jp2 , jp3
    real(rkx) :: dflux_x_p0 , dflux_x_p1 , dflux_y_p0 , dflux_y_p1
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'diffu_x3df'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    if ( idiffu == 1 ) then
      !
      ! fourth-order scheme for interior:
      !
      do concurrent ( j = jcii1:jcii2 , i = icii1:icii2 , k = 1:kzp1 )
        ften(j,i,k) = ften(j,i,k) - fac * xkcf(j,i,k) * &
             (z4_c1*(f(j+2,i,k)+f(j-2,i,k)+f(j,i+2,k)+f(j,i-2,k)) +   &
              z4_c2*(f(j+1,i,k)+f(j-1,i,k)+f(j,i+1,k)+f(j,i-1,k)) + &
              z4_c3*(f(j,i,k)))
      end do
      !
      ! second-order scheme for boundaries:
      !
      if ( ma%has_bdyleft ) then
        j = jci1
        do k = 1 , kzp1
          do i = ici1 , ici2
            ften(j,i,k) = ften(j,i,k) + fac * xkcf(j,i,k) * &
              (z4_c1*(f(j+1,i,k)+f(j-1,i,k)+f(j,i+1,k)+f(j,i-1,k)) + &
               z4_c2*(f(j,i,k)))
          end do
        end do
      end if
      if ( ma%has_bdyright ) then
        j = jci2
        do k = 1 , kzp1
          do i = ici1 , ici2
            ften(j,i,k) = ften(j,i,k) + fac * xkcf(j,i,k) * &
              (z4_c1*(f(j+1,i,k)+f(j-1,i,k)+f(j,i+1,k)+f(j,i-1,k)) + &
               z4_c2*(f(j,i,k)))
          end do
        end do
      end if
      if ( ma%has_bdybottom ) then
        i = ici1
        do k = 1 , kzp1
          do j = jci1 , jci2
            ften(j,i,k) = ften(j,i,k) + fac * xkcf(j,i,k) * &
              (z4_c1*(f(j+1,i,k)+f(j-1,i,k)+f(j,i+1,k)+f(j,i-1,k)) + &
               z4_c2*(f(j,i,k)))
          end do
        end do
      end if
      if ( ma%has_bdytop ) then
        i = ici2
        do k = 1 , kzp1
          do j = jci1 , jci2
            ften(j,i,k) = ften(j,i,k) + fac * xkcf(j,i,k) * &
              (z4_c1*(f(j+1,i,k)+f(j-1,i,k)+f(j,i+1,k)+f(j,i-1,k)) + &
               z4_c2*(f(j,i,k)))
          end do
        end do
      end if
    else if ( idiffu == 2 ) then
      !
      ! fourth-order scheme
      !
      do concurrent ( j = jci1:jci2 , i = ici1:ici2 , k = 1:kzp1 )
        ften(j,i,k) = ften(j,i,k) + fac * xkcf(j,i,k) * &
             (o4_c1*(f(j+1,i,k)+f(j-1,i,k)+f(j,i+1,k)+f(j,i-1,k)) +   &
              o4_c2*(f(j+1,i+1,k)+f(j-1,i-1,k)+f(j-1,i+1,k)+f(j+1,i-1,k)) + &
              o4_c3*(f(j,i,k)))
      end do
    else if ( idiffu == 3 ) then
      do k = 1 , kzp1
        do i = ici1 , ici2
          im1 = max(i-1,1)
          im2 = max(i-2,1)
          im3 = max(i-3,1)
          ip1 = min(i+1,iym1)
          ip2 = min(i+2,iym1)
          ip3 = min(i+3,iym1)
          do j = jci2 , jci2
            jm1 = max(j-1,1)
            jm2 = max(j-2,1)
            jm3 = max(j-3,1)
            jp1 = min(j+1,jxm1)
            jp2 = min(j+2,jxm1)
            jp3 = min(j+3,jxm1)
            dflux_x_p0 = ( h4_c1 * (f(j,i,k)   - f(jm1,i,k) ) + &
                           h4_c2 * (f(jp1,i,k) - f(jm2,i,k) ) + &
                           h4_c3 * (f(jp2,i,k) - f(jm3,i,k) ) )
            if  ( dflux_x_p0 * (f(j,i,k)/mpd(j,i) - &
                                f(jm1,i,k)/mpd(jm1,i)) <= d_zero ) then
              dflux_x_p0 = d_zero
            end if
            dflux_x_p1 = ( h4_c1 * (f(jp1,i,k) - f(j,i,k) )   + &
                           h4_c2 * (f(jp2,i,k) - f(jm1,i,k) ) + &
                           h4_c3 * (f(jp3,i,k) - f(jm2,i,k) ) )
            if  ( dflux_x_p1 * (f(jp1,i,k)/mpd(jp1,i) - &
                                f(j,i,k)/mpd(j,i)) <= d_zero ) then
              dflux_x_p1 = d_zero
            end if
            dflux_y_p0 = ( h4_c1 * (f(j,i,k)   - f(j,im1,k) ) + &
                           h4_c2 * (f(j,ip1,k) - f(j,im2,k) ) + &
                           h4_c3 * (f(j,ip2,k) - f(j,im3,k) ) )
            if  ( dflux_y_p0 * (f(j,i,k)/mpd(j,i) - &
                                f(j,im1,k)/mpd(j,im1)) <= d_zero ) then
              dflux_y_p0 = d_zero
            end if
            dflux_y_p1 = ( h4_c1 * (f(j,ip1,k) - f(j,i,k) )   + &
                           h4_c2 * (f(j,ip2,k) - f(j,im1,k) ) + &
                           h4_c3 * (f(j,ip3,k) - f(j,im2,k) ) )
            if  ( dflux_y_p1 * (f(j,ip1,k)/mpd(j,ip1) - &
                                f(j,i,k)/mpd(j,i)) <= d_zero ) then
              dflux_y_p1 = d_zero
            end if
            ften(j,i,k) = ften(j,i,k) + fac * xkc(j,i,k) * &
                     ( ( dflux_x_p1 - dflux_x_p0 ) + &
                       ( dflux_y_p1 - dflux_y_p0 ) )
          end do
        end do
      end do
    end if
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine diffu_x3df

  subroutine diffu_x3d(ften,f)
    implicit none
    real(rkx) , pointer , dimension(:,:,:) , intent(in) :: f
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: ften
    integer(ik4) :: i , j , k
    integer(ik4) :: im1 , im2 , im3
    integer(ik4) :: jm1 , jm2 , jm3
    integer(ik4) :: ip1 , ip2 , ip3
    integer(ik4) :: jp1 , jp2 , jp3
    real(rkx) :: dflux_x_p0 , dflux_x_p1 , dflux_y_p0 , dflux_y_p1
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'diffu_x3d'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    if ( idiffu == 1 ) then
      !
      ! fourth-order scheme for interior:
      !
      do concurrent ( j = jcii1:jcii2 , i = icii1:icii2 , k = 1:kz )
        ften(j,i,k) = ften(j,i,k) - xkc(j,i,k) * &
             (z4_c1*(f(j+2,i,k)+f(j-2,i,k)+f(j,i+2,k)+f(j,i-2,k)) +   &
              z4_c2*(f(j+1,i,k)+f(j-1,i,k)+f(j,i+1,k)+f(j,i-1,k)) + &
              z4_c3*(f(j,i,k)))
      end do
      !
      ! second-order scheme for boundaries:
      !
      if ( ma%has_bdyleft ) then
        j = jci1
        do k = 1 , kz
          do i = ici1 , ici2
            ften(j,i,k) = ften(j,i,k) + xkc(j,i,k) * &
              (z4_c1*(f(j+1,i,k)+f(j-1,i,k)+f(j,i+1,k)+f(j,i-1,k)) + &
               z4_c2*(f(j,i,k)))
          end do
        end do
      end if
      if ( ma%has_bdyright ) then
        j = jci2
        do k = 1 , kz
          do i = ici1 , ici2
            ften(j,i,k) = ften(j,i,k) + xkc(j,i,k) * &
              (z4_c1*(f(j+1,i,k)+f(j-1,i,k)+f(j,i+1,k)+f(j,i-1,k)) + &
               z4_c2*(f(j,i,k)))
          end do
        end do
      end if
      if ( ma%has_bdybottom ) then
        i = ici1
        do k = 1 , kz
          do j = jci1 , jci2
            ften(j,i,k) = ften(j,i,k) + xkc(j,i,k) * &
              (z4_c1*(f(j+1,i,k)+f(j-1,i,k)+f(j,i+1,k)+f(j,i-1,k)) + &
               z4_c2*(f(j,i,k)))
          end do
        end do
      end if
      if ( ma%has_bdytop ) then
        i = ici2
        do k = 1 , kz
          do j = jci1 , jci2
            ften(j,i,k) = ften(j,i,k) + xkc(j,i,k) * &
              (z4_c1*(f(j+1,i,k)+f(j-1,i,k)+f(j,i+1,k)+f(j,i-1,k)) + &
               z4_c2*(f(j,i,k)))
          end do
        end do
      end if
    else if ( idiffu == 2 ) then
      !
      ! fourth-order scheme
      !
      do concurrent ( j = jci1:jci2 , i = ici1:ici2 , k = 1:kz )
        ften(j,i,k) = ften(j,i,k) + xkc(j,i,k) * &
             (o4_c1*(f(j+1,i,k)+f(j-1,i,k)+f(j,i+1,k)+f(j,i-1,k)) +   &
              o4_c2*(f(j+1,i+1,k)+f(j-1,i-1,k)+f(j-1,i+1,k)+f(j+1,i-1,k)) + &
              o4_c3*(f(j,i,k)))
      end do
    else if ( idiffu == 3 ) then
      do k = 1 , kz
        do i = ici1 , ici2
          im1 = max(i-1,1)
          im2 = max(i-2,1)
          im3 = max(i-3,1)
          ip1 = min(i+1,iym1)
          ip2 = min(i+2,iym1)
          ip3 = min(i+3,iym1)
          do j = jci2 , jci2
            jm1 = max(j-1,1)
            jm2 = max(j-2,1)
            jm3 = max(j-3,1)
            jp1 = min(j+1,jxm1)
            jp2 = min(j+2,jxm1)
            jp3 = min(j+3,jxm1)
            dflux_x_p0 = ( h4_c1 * (f(j,i,k)   - f(jm1,i,k) ) + &
                           h4_c2 * (f(jp1,i,k) - f(jm2,i,k) ) + &
                           h4_c3 * (f(jp2,i,k) - f(jm3,i,k) ) )
            if  ( dflux_x_p0 * (f(j,i,k)/mpd(j,i) - &
                                f(jm1,i,k)/mpd(jm1,i)) <= d_zero ) then
              dflux_x_p0 = d_zero
            end if
            dflux_x_p1 = ( h4_c1 * (f(jp1,i,k) - f(j,i,k) )   + &
                           h4_c2 * (f(jp2,i,k) - f(jm1,i,k) ) + &
                           h4_c3 * (f(jp3,i,k) - f(jm2,i,k) ) )
            if  ( dflux_x_p1 * (f(jp1,i,k)/mpd(jp1,i) - &
                                f(j,i,k)/mpd(j,i)) <= d_zero ) then
              dflux_x_p1 = d_zero
            end if
            dflux_y_p0 = ( h4_c1 * (f(j,i,k)   - f(j,im1,k) ) + &
                           h4_c2 * (f(j,ip1,k) - f(j,im2,k) ) + &
                           h4_c3 * (f(j,ip2,k) - f(j,im3,k) ) )
            if  ( dflux_y_p0 * (f(j,i,k)/mpd(j,i) - &
                                f(j,im1,k)/mpd(j,im1)) <= d_zero ) then
              dflux_y_p0 = d_zero
            end if
            dflux_y_p1 = ( h4_c1 * (f(j,ip1,k) - f(j,i,k) )   + &
                           h4_c2 * (f(j,ip2,k) - f(j,im1,k) ) + &
                           h4_c3 * (f(j,ip3,k) - f(j,im2,k) ) )
            if  ( dflux_y_p1 * (f(j,ip1,k)/mpd(j,ip1) - &
                                f(j,i,k)/mpd(j,i)) <= d_zero ) then
              dflux_y_p1 = d_zero
            end if
            ften(j,i,k) = ften(j,i,k) + xkc(j,i,k) * &
                     ( ( dflux_x_p1 - dflux_x_p0 ) + &
                       ( dflux_y_p1 - dflux_y_p0 ) )
          end do
        end do
      end do
    end if
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
    integer(ik4) :: im1 , im2 , im3
    integer(ik4) :: jm1 , jm2 , jm3
    integer(ik4) :: ip1 , ip2 , ip3
    integer(ik4) :: jp1 , jp2 , jp3
    real(rkx) :: dflux_x_p0 , dflux_x_p1 , dflux_y_p0 , dflux_y_p1
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'diffu_x4d3d'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    if ( idiffu == 1 ) then
      !
      ! fourth-order scheme for interior:
      !
      do concurrent ( j = jcii1:jcii2 , i = icii1:icii2 , k = 1:kz )
        ften(j,i,k,n) = ften(j,i,k,n) - fac * xkc(j,i,k) * &
             (z4_c1*(f(j+2,i,k,n)+f(j-2,i,k,n) +  &
                     f(j,i+2,k,n)+f(j,i-2,k,n)) + &
              z4_c2*(f(j+1,i,k,n)+f(j-1,i,k,n) +  &
                     f(j,i+1,k,n)+f(j,i-1,k,n)) + &
              z4_c3*f(j,i,k,n))
      end do
      !
      ! second-order scheme for boundaries:
      !
      if ( ma%has_bdyleft ) then
        j = jci1
        do k = 1 , kz
          do i = ici1 , ici2
            ften(j,i,k,n) = ften(j,i,k,n) + fac * xkc(j,i,k) * &
             (z4_c1*(f(j+1,i,k,n)+f(j-1,i,k,n) +  &
                     f(j,i+1,k,n)+f(j,i-1,k,n)) + &
              z4_c2*f(j,i,k,n))
          end do
        end do
      end if
      if ( ma%has_bdyright ) then
        j = jci2
        do k = 1 , kz
          do i = ici1 , ici2
            ften(j,i,k,n) = ften(j,i,k,n) + fac * xkc(j,i,k) * &
             (z4_c1*(f(j+1,i,k,n)+f(j-1,i,k,n) +  &
                     f(j,i+1,k,n)+f(j,i-1,k,n)) + &
              z4_c2*f(j,i,k,n))
          end do
        end do
      end if
      if ( ma%has_bdybottom ) then
        i = ici1
        do k = 1 , kz
          do j = jci1 , jci2
            ften(j,i,k,n) = ften(j,i,k,n) + fac * xkc(j,i,k) * &
             (z4_c1*(f(j+1,i,k,n)+f(j-1,i,k,n) +  &
                     f(j,i+1,k,n)+f(j,i-1,k,n)) + &
              z4_c2*f(j,i,k,n))
          end do
        end do
      end if
      if ( ma%has_bdytop ) then
        i = ici2
        do k = 1 , kz
          do j = jci1 , jci2
            ften(j,i,k,n) = ften(j,i,k,n) + fac * xkc(j,i,k) * &
             (z4_c1*(f(j+1,i,k,n)+f(j-1,i,k,n) +  &
                     f(j,i+1,k,n)+f(j,i-1,k,n)) + &
              z4_c2*f(j,i,k,n))
          end do
        end do
      end if
    else if ( idiffu == 2 ) then
      !
      ! fourth-order scheme
      !
      do concurrent ( j = jci1:jci2 , i = ici1:ici2 , k = 1:kz )
        ften(j,i,k,n) = ften(j,i,k,n) + fac * xkc(j,i,k) * &
             (o4_c1*(f(j+1,i,k,n)+f(j-1,i,k,n) +      &
                     f(j,i+1,k,n)+f(j,i-1,k,n)) +     &
              o4_c2*(f(j+1,i+1,k,n)+f(j-1,i-1,k,n) +  &
                     f(j-1,i+1,k,n)+f(j+1,i-1,k,n)) + &
              o4_c3*f(j,i,k,n))
      end do
    else if ( idiffu == 3 ) then
      do k = 1 , kz
        do i = ici1 , ici2
          im1 = max(i-1,1)
          im2 = max(i-2,1)
          im3 = max(i-3,1)
          ip1 = min(i+1,iym1)
          ip2 = min(i+2,iym1)
          ip3 = min(i+3,iym1)
          do j = jci2 , jci2
            jm1 = max(j-1,1)
            jm2 = max(j-2,1)
            jm3 = max(j-3,1)
            jp1 = min(j+1,jxm1)
            jp2 = min(j+2,jxm1)
            jp3 = min(j+3,jxm1)
            dflux_x_p0 = ( h4_c1 * (f(j,i,k,n)   - f(jm1,i,k,n) ) + &
                           h4_c2 * (f(jp1,i,k,n) - f(jm2,i,k,n) ) + &
                           h4_c3 * (f(jp2,i,k,n) - f(jm3,i,k,n) ) )
            if  ( dflux_x_p0 * (f(j,i,k,n)/mpd(j,i) - &
                                f(jm1,i,k,n)/mpd(jm1,i)) <= d_zero ) then
              dflux_x_p0 = d_zero
            end if
            dflux_x_p1 = ( h4_c1 * (f(jp1,i,k,n) - f(j,i,k,n) )   + &
                           h4_c2 * (f(jp2,i,k,n) - f(jm1,i,k,n) ) + &
                           h4_c3 * (f(jp3,i,k,n) - f(jm2,i,k,n) ) )
            if  ( dflux_x_p1 * (f(jp1,i,k,n)/mpd(jp1,i) - &
                                f(j,i,k,n)/mpd(j,i)) <= d_zero ) then
              dflux_x_p1 = d_zero
            end if
            dflux_y_p0 = ( h4_c1 * (f(j,i,k,n)   - f(j,im1,k,n) ) + &
                           h4_c2 * (f(j,ip1,k,n) - f(j,im2,k,n) ) + &
                           h4_c3 * (f(j,ip2,k,n) - f(j,im3,k,n) ) )
            if  ( dflux_y_p0 * (f(j,i,k,n)/mpd(j,i) - &
                                f(j,im1,k,n)/mpd(j,im1)) <= d_zero ) then
              dflux_y_p0 = d_zero
            end if
            dflux_y_p1 = ( h4_c1 * (f(j,ip1,k,n) - f(j,i,k,n) )   + &
                           h4_c2 * (f(j,ip2,k,n) - f(j,im1,k,n) ) + &
                           h4_c3 * (f(j,ip3,k,n) - f(j,im2,k,n) ) )
            if  ( dflux_y_p1 * (f(j,ip1,k,n)/mpd(j,ip1) - &
                                f(j,i,k,n)/mpd(j,i)) <= d_zero ) then
              dflux_y_p1 = d_zero
            end if
            ften(j,i,k,n) = ften(j,i,k,n) + xkc(j,i,k) * &
                     ( ( dflux_x_p1 - dflux_x_p0 ) + &
                       ( dflux_y_p1 - dflux_y_p0 ) )
          end do
        end do
      end do
    end if
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine diffu_x4d3d

end module mod_diffusion
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
