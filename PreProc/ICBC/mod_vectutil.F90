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

module mod_vectutil

  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_message

  private

  public :: crs2dot , dot2crs , top2btm , btm2top , meandiv

  contains

  subroutine crs2dot(pd,px,ni,nj,iband)
    implicit none
    integer(ik4) , intent(in) :: ni , nj
    integer(ik4) , intent(in) :: iband
    real(rk8) , intent(in) , dimension(ni,nj) :: px
    real(rk8) , intent(out) , dimension(ni,nj) :: pd

    integer(ik4) :: i , j , ni1 , nj1 , im1
    !
    ! THIS ROUTINE DETERMINES P(.) FROM P(X) BY A 4-POINT INTERPOLATION.
    ! ON THE X-GRID, A P(X) POINT OUTSIDE THE GRID DOMAIN IS ASSUMED TO
    ! SATISFY P(0,J) = P(1,J); P(NI,J) = P(NI-1,J); AND SIMILARLY FOR THE
    ! I'S.
    !
    if ( iband == 1 ) then
      nj1 = nj - 1
      do j = 2 , nj1
        do i = 1 , ni
          im1 = i-1
          if (im1 == 0) im1 = ni
          pd(i,j) = d_rfour*(px(i,j)+px(im1,j)+px(i,j-1)+px(im1,j-1))
        end do
      end do
      do i = 1 , ni
        im1 = i-1
        if (im1 == 0) im1 = ni
        pd(i,1) = d_half*(px(i,1)+px(im1,1))
        pd(i,nj) = d_half*(px(i,nj1)+px(im1,nj1))
      end do
    else
      ni1 = ni - 1
      nj1 = nj - 1
      do j = 2 , nj1
        do i = 2 , ni1
          pd(i,j) = d_rfour*(px(i,j)+px(i-1,j)+px(i,j-1)+px(i-1,j-1))
        end do
      end do
      do i = 2 , ni1
        pd(i,1) = d_half*(px(i,1)+px(i-1,1))
        pd(i,nj) = d_half*(px(i,nj1)+px(i-1,nj1))
      end do
      do j = 2 , nj1
        pd(1,j) = d_half*(px(1,j)+px(1,j-1))
        pd(ni,j) = d_half*(px(ni1,j)+px(ni1,j-1))
      end do
      pd(1,1) = px(1,1)
      pd(1,nj) = px(1,nj1)
      pd(ni,1) = px(ni1,1)
      pd(ni,nj) = px(ni1,nj1)
    end if
  end subroutine crs2dot
  !
  !-----------------------------------------------------------------------
  !
  subroutine dot2crs(px,pd,ni,nj,iband)
    implicit none
    integer(ik4) , intent(in) :: ni , nj
    real(rk8) , intent(in) , dimension(ni,nj) :: pd
    real(rk8) , intent(out) , dimension(ni,nj) :: px
    integer(ik4) , intent(in) :: iband
    integer(ik4) :: i , j
    do j = 1 , nj - 1
      do i = 1 , ni - 1
        px(i,j) = ( pd(i  ,j  ) + &
                    pd(i+1,j  ) + &
                    pd(i  ,j+1) + &
                    pd(i+1,j+1) ) * 0.25D0
      end do
    end do
    if ( iband == 1 ) then
      do j = 1 , nj - 1
        px(ni,j) = ( pd(ni,j  ) + &
                     pd(1 ,j  ) + &
                     pd(ni,j+1) + &
                     pd(1, j+1) ) * 0.25D0
      end do
    end if
  end subroutine dot2crs
  !
  !-----------------------------------------------------------------------
  !
  subroutine top2btm(x,nlon1,nlat1,nlev1)
    implicit none
    integer(ik4) , intent(in) :: nlat1 , nlev1 , nlon1
    real(rk8) , intent(inout) , dimension(nlon1,nlat1,nlev1) :: x

    integer(ik4) :: i , j , k , kr
    real(rk8) , dimension(nlev1) :: work

    do j = 1 , nlat1
      do i = 1 , nlon1
        do k = 1 , nlev1
          work(k) = x(i,j,k)
        end do
        do k = 1 , nlev1
          kr = nlev1 - k + 1
          x(i,j,k) = work(kr)
        end do
      end do
    end do
  end subroutine top2btm
  !
  !-----------------------------------------------------------------------
  !
  subroutine btm2top(x,nlon1,nlat1,nlev1)
    implicit none
    integer(ik4) , intent(in) :: nlat1 , nlev1 , nlon1
    real(rk8) , intent(inout) , dimension(nlon1,nlat1,nlev1) :: x
    call top2btm(x,nlon1,nlat1,nlev1)
  end subroutine btm2top
  !
  !-----------------------------------------------------------------------
  !
  subroutine relax(chi,ff,rd,imx,jmx,ie,je,ds)
    implicit none
    integer(ik4) , intent(in) :: imx , jmx , ie , je
    real(rk8) , intent(inout) , dimension(imx,jmx) :: chi
    real(rk8) , intent(inout) , dimension(imx,jmx) :: rd
    real(rk8) , intent(inout) , dimension(imx,jmx) :: ff
    real(rk8) , intent(in) :: ds
    real(rk8) , parameter :: smallres = 1.0D-9
    integer(ik4) , parameter :: mm = 20000
    real(rk8) , parameter :: alpha = 1.8D0
    real(rk8) , parameter :: alphaov4 = alpha * d_rfour
    integer(ik4) :: i , j , iter , mi
    real(rk8) , dimension(jmx) :: chimx
    real(rk8) , dimension(jmx) :: rdmax
    real(rk8) :: epx , fac
    logical :: converged
    converged = .false.
    fac = d_two * ds * ds
    rd = d_zero
    do j = 1 , je + 1
      do i = 1 , ie + 1
        ff(i,j) = fac * ff(i,j)
      end do
    end do
    iter_loop : do iter = 1 , mm
      mi = iter
      chimx = d_zero
      do j = 2 , je
        do i = 2 , ie
          chimx(j) = max(abs(chi(i,j)),chimx(j))
        end do
      end do
      epx = maxval(chimx) * smallres * d_four / alpha
      do j = 2 , je
        do i = 2 , ie
          rd(i,j) = chi(i+1,j+1) + chi(i-1,j+1) + &
                    chi(i+1,j-1) + chi(i-1,j-1) - &
                    4.0 * chi(i,j) - ff(i,j)
          chi(i,j) = chi(i,j) + rd(i,j) * alphaov4
        end do
      end do
      rdmax = d_zero
      do j = 2 , je
        do i = 2 , ie
          rdmax(j) = max(abs(rd(i,j)),rdmax(j))
        end do
      end do
      if ( maxval(rdmax) < epx) then
        converged = .true.
        exit iter_loop
      end if
    end do iter_loop
    if ( .not. converged ) then
      call fatal(__FILE__,__LINE__,'Relaxation did not converge !')
    end if
  end subroutine relax
  !
  !-----------------------------------------------------------------------
  !
  subroutine fill(f,ix,jx,imx,jmx,ifirst,ilast,jfirst,jlast)
    implicit none
    integer(ik4) , intent(in) :: ix , jx , imx , jmx , &
                                 ifirst , ilast , jfirst , jlast
    real(rk8) , intent(inout) , dimension(ix,jx) :: f
    integer(ik4) :: i , j
    do j = jfirst , jlast
      do i = 1 , ifirst - 1
        f(i,j) = f(ifirst,j)
      end do
      do i = ilast + 1 , imx
        f(i,j) = f(ilast,j)
      end do
    end do
    do j = 1 , jfirst - 1
      f(:,j) = f(:,jfirst)
    end do
    do j = jlast + 1 , jmx
      f(:,j) = f(:,jlast)
    end do
  end subroutine fill
  !
  !-----------------------------------------------------------------------
  !
  subroutine meandiv(u,v,psd,dm,sigh,imx,jmx,kxs,ds,imxm,jmxm)
    implicit none
    integer(ik4) , intent(in) :: imx , jmx , kxs , imxm , jmxm
    real(rk8) , intent(inout) , dimension(imx,jmx,kxs) :: u , v
    real(rk8) , intent(in) , dimension(imx,jmx) :: psd
    real(rk8) , intent(in) , dimension(imx,jmx) :: dm
    real(rk8) , intent(in) , dimension(kxs) :: sigh
    real(rk8) , intent(in) :: ds
    integer(ik4) :: i , j , k
    real(rk8) , dimension(imx,jmx) :: chi , div , f
    real(rk8) , dimension(imx,jmx) :: dudx , dvdy
    real(rk8) , dimension(imx,jmx) :: udiverg , vdiverg
    real(rk8) , dimension(imx,jmx) :: uslb , vslb
    real(rk8) , dimension(kxs) :: dsg , weight
    real(rk8) , dimension(kxs+1) :: sigf
    real(rk8) :: oneov2ds

    oneov2ds = d_one / (d_two * ds)
    !
    ! Integrate p* v/m, compute div, to dot point, to (x,y) format
    !
    sigf(1) = d_zero
    do k = 1, kxs
      sigf(k+1) = d_two * sigh(k) - sigf(k)
      dsg(k) = sigf(k+1) - sigf(k)
    end do
    do j = 1 , jmx
      do i = 1 , imx
        uslb(i,j) = d_zero
        vslb(i,j) = d_zero
      end do
    end do
    do j = 1 , jmx
      do i = 1 , imx
        do k = 1 , kxs
          uslb(i,j) = uslb(i,j) + u(i,j,k) * dsg(k)
          vslb(i,j) = vslb(i,j) + v(i,j,k) * dsg(k)
        end do
      end do
    end do
    do j = 1 , jmx
      do i = 1 , imx
        uslb(i,j) = uslb(i,j) * psd(i,j) / dm(i,j)
        vslb(i,j) = vslb(i,j) * psd(i,j) / dm(i,j)
      end do
    end do
    do j = 1 , jmxm
      do i = 1 , imxm
        dudx(i,j) = uslb(i+1,j+1) - uslb(i+1,j) + uslb(i,j+1) - uslb(i,j)
        dvdy(i,j) = vslb(i+1,j+1) - vslb(i,j+1) + vslb(i+1,j) - vslb(i,j)
      end do
    end do
    do j = 1 , jmxm
      do i = 1 , imxm
        div(i,j) = oneov2ds * (dudx(i,j) + dvdy(i,j))
      end do
    end do
    !
    ! Iteratively solve laplacian from good first guess.
    !
    do j = 1 , jmx
      do i = 1 , imx
        chi(i,j) = d_zero
      end do
    end do
    call relax(chi,div,f,imx,jmx,imx-2,jmx-2,ds)
    !
    ! Get divergent component of wind, 2d field on dot points.
    !
    do j = 2 , jmxm
      do i = 2 , imxm
         udiverg(i,j) = (chi(i,j) - chi(i,j-1) + &
                         chi(i-1,j) - chi(i-1,j-1)) * oneov2ds
      end do
    end do
    do j = 2 , jmxm
      do i = 2 , imxm
        vdiverg(i,j) = (chi(i,j) - chi(i-1,j) + &
                        chi(i,j-1) - chi(i-1,j-1)) * oneov2ds
      end do
    end do
    call fill(udiverg,imx,jmx,imx,jmx,2,imxm,2,jmxm)
    call fill(vdiverg,imx,jmx,imx,jmx,2,imxm,2,jmxm)
    !
    ! Remove mean divergent component
    !
    do j = 1 , jmx
      do i = 1 , imx
        udiverg(i,j) = udiverg(i,j) * dm(i,j) / psd(i,j)
        vdiverg(i,j) = vdiverg(i,j) * dm(i,j) / psd(i,j)
      end do
    end do
    do k = 1 , kxs
      weight(k) = d_two * (d_one - sigh(k))
    end do
    do j = 1 , jmx
      do i = 1 , imx
        do k = 1 , kxs
          u(i,j,k) = u(i,j,k) - weight(k) * udiverg(i,j)
          v(i,j,k) = v(i,j,k) - weight(k) * vdiverg(i,j)
        end do
      end do
    end do
  end subroutine meandiv

end module mod_vectutil

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
