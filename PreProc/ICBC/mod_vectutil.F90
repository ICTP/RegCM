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

  public :: crs2dot , dot2crs , top2btm , btm2top , meandiv , meandivf
  public :: ucrs2dot , vcrs2dot

  interface ucrs2dot
    module procedure ucrs2dot_2d
    module procedure ucrs2dot_3d
  end interface ucrs2dot

  interface vcrs2dot
    module procedure vcrs2dot_2d
    module procedure vcrs2dot_3d
  end interface vcrs2dot

  interface crs2dot
    module procedure crs2dot_2d
    module procedure crs2dot_3d
  end interface crs2dot

  contains

  subroutine ucrs2dot_3d(pd,px,ni,nj,nk,iband)
    implicit none
    integer(ik4) , intent(in) :: ni , nj , nk
    integer(ik4) , intent(in) :: iband
    real(rkx) , intent(in) , dimension(ni,nj,nk) :: px
    real(rkx) , intent(out) , dimension(ni,nj,nk) :: pd

    integer(ik4) :: i , j , k , ni1 , im1 , im2 , ip1
    !
    ! THIS ROUTINE DETERMINES P(.) FROM P(X) BY A 4-POINT INTERPOLATION.
    ! ON THE X-GRID, A P(X) POINT OUTSIDE THE GRID DOMAIN IS ASSUMED TO
    ! SATISFY P(0,J) = P(1,J); P(NI,J) = P(NI-1,J); AND SIMILARLY FOR THE
    ! I'S.
    !
    if ( iband == 1 ) then
      do k = 1 , nk
        do j = 1 , nj
          do i = 1 , ni
            im1 = i-1
            im2 = i-2
            ip1 = i+1
            if (ip1 > ni) ip1 = 1
            if (im1 == 0) im1 = ni
            if (im2 == 0) im2 = ni
            if (im2 == -1) im2 = ni-1
            pd(i,j,k) = 0.5625_rkx*(px(i,j,k)+px(im1,j,k)) - &
                        0.0625_rkx*(px(ip1,j,k)+px(im2,j,k))
          end do
        end do
      end do
    else
      do k = 1 , nk
        ni1 = ni - 1
        do j = 1 , nj
          do i = 3 , ni1
            pd(i,j,k) = 0.5625_rkx*(px(i,j,k)+px(i-1,j,k)) - &
                        0.0625_rkx*(px(i+1,j,k)+px(i-2,j,k))
          end do
        end do
        pd(1,:,k) = px(1,:,k)
        pd(2,:,k) = 0.5_rkx * (px(1,:,k)+px(2,:,k))
        pd(ni,:,k) = 0.5_rkx * (px(ni-1,:,k)+px(ni,:,k))
      end do
    end if
  end subroutine ucrs2dot_3d

  subroutine vcrs2dot_3d(pd,px,ni,nj,nk,icrs)
    implicit none
    integer(ik4) , intent(in) :: ni , nj , nk
    integer(ik4) , intent(in) :: icrs
    real(rkx) , intent(in) , dimension(ni,nj,nk) :: px
    real(rkx) , intent(out) , dimension(ni,nj,nk) :: pd

    integer(ik4) :: i , j , k , nj1 , jm1 , jm2 , jp1
    !
    ! THIS ROUTINE DETERMINES P(.) FROM P(X) BY A 4-POINT INTERPOLATION.
    ! ON THE X-GRID, A P(X) POINT OUTSIDE THE GRID DOMAIN IS ASSUMED TO
    ! SATISFY P(0,J) = P(1,J); P(NI,J) = P(NI-1,J); AND SIMILARLY FOR THE
    ! I'S.
    !
    if ( icrs == 1 ) then
      do k = 1 , nk
        do j = 1 , nj
          jm1 = j-1
          jm2 = j-2
          jp1 = j+1
          if (jp1 > ni) jp1 = 1
          if (jm1 == 0) jm1 = nj
          if (jm2 == 0) jm2 = nj
          if (jm2 == -1) jm2 = nj-1
          do i = 1 , ni
            pd(i,j,k) = 0.5625_rkx*(px(i,j,k)+px(i,jm1,k)) - &
                        0.0625_rkx*(px(i,jp1,k)+px(i,jm2,k))
          end do
        end do
      end do
    else
      do k = 1 , nk
        nj1 = nj - 1
        do j = 3 , nj1
          do i = 1 , ni
            pd(i,j,k) = 0.5625_rkx*(px(i,j,k)+px(i,j-1,k)) - &
                        0.0625_rkx*(px(i,j+1,k)+px(i,j-2,k))
          end do
        end do
        pd(:,1,k) = px(:,1,k)
        pd(:,2,k) = 0.5_rkx * (px(:,1,k)+px(:,2,k))
        pd(:,nj,k) = 0.5_rkx * (px(:,nj-1,k)+px(:,nj,k))
      end do
    end if
  end subroutine vcrs2dot_3d

  subroutine ucrs2dot_2d(pd,px,ni,nj,iband)
    implicit none
    integer(ik4) , intent(in) :: ni , nj
    integer(ik4) , intent(in) :: iband
    real(rkx) , intent(in) , dimension(ni,nj) :: px
    real(rkx) , intent(out) , dimension(ni,nj) :: pd

    integer(ik4) :: i , j , ni1 , im1 , im2 , ip1
    !
    ! THIS ROUTINE DETERMINES P(.) FROM P(X) BY A 4-POINT INTERPOLATION.
    ! ON THE X-GRID, A P(X) POINT OUTSIDE THE GRID DOMAIN IS ASSUMED TO
    ! SATISFY P(0,J) = P(1,J); P(NI,J) = P(NI-1,J); AND SIMILARLY FOR THE
    ! I'S.
    !
    if ( iband == 1 ) then
      do j = 1 , nj
        do i = 1 , ni
          im1 = i-1
          im2 = i-2
          ip1 = i+1
          if (ip1 > ni) ip1 = 1
          if (im1 == 0) im1 = ni
          if (im2 == 0) im2 = ni
          if (im2 == -1) im2 = ni-1
          pd(i,j) = 0.5625_rkx*(px(i,j)+px(im1,j)) - &
                    0.0625_rkx*(px(ip1,j)+px(im2,j))
        end do
      end do
    else
      ni1 = ni - 1
      do j = 1 , nj
        do i = 3 , ni1
          pd(i,j) = 0.5625_rkx*(px(i,j)+px(i-1,j)) - &
                    0.0625_rkx*(px(i+1,j)+px(i-2,j))
        end do
      end do
      pd(1,:) = px(1,:)
      pd(2,:) = 0.5_rkx * (px(1,:)+px(2,:))
      pd(ni,:) = 0.5_rkx * (px(ni-1,:)+px(ni,:))
    end if
  end subroutine ucrs2dot_2d

  subroutine vcrs2dot_2d(pd,px,ni,nj,icrs)
    implicit none
    integer(ik4) , intent(in) :: ni , nj
    integer(ik4) , intent(in) :: icrs
    real(rkx) , intent(in) , dimension(ni,nj) :: px
    real(rkx) , intent(out) , dimension(ni,nj) :: pd

    integer(ik4) :: i , j , nj1 , jm1 , jm2 , jp1
    !
    ! THIS ROUTINE DETERMINES P(.) FROM P(X) BY A 4-POINT INTERPOLATION.
    ! ON THE X-GRID, A P(X) POINT OUTSIDE THE GRID DOMAIN IS ASSUMED TO
    ! SATISFY P(0,J) = P(1,J); P(NI,J) = P(NI-1,J); AND SIMILARLY FOR THE
    ! I'S.
    !
    if ( icrs == 1 ) then
      do j = 1 , nj
        jm1 = j-1
        jm2 = j-2
        jp1 = j+1
        if (jp1 > ni) jp1 = 1
        if (jm1 == 0) jm1 = nj
        if (jm2 == 0) jm2 = nj
        if (jm2 == -1) jm2 = nj-1
        do i = 1 , ni
          pd(i,j) = 0.5625_rkx*(px(i,j)+px(i,jm1)) - &
                    0.0625_rkx*(px(i,jp1)+px(i,jm2))
        end do
      end do
    else
      nj1 = nj - 1
      do j = 3 , nj1
        do i = 1 , ni
          pd(i,j) = 0.5625_rkx*(px(i,j)+px(i,j-1)) - &
                    0.0625_rkx*(px(i,j+1)+px(i,j-2))
        end do
      end do
      pd(:,1) = px(:,1)
      pd(:,2) = 0.5_rkx * (px(:,1)+px(:,2))
      pd(:,nj) = 0.5_rkx * (px(:,nj-1)+px(:,nj))
    end if
  end subroutine vcrs2dot_2d

  subroutine crs2dot_2d(pd,px,ni,nj,iband,icrm)
    implicit none
    integer(ik4) , intent(in) :: ni , nj
    integer(ik4) , intent(in) :: iband , icrm
    real(rkx) , intent(in) , dimension(ni,nj) :: px
    real(rkx) , intent(out) , dimension(ni,nj) :: pd

    integer(ik4) :: i , j , ni1 , nj1 , im1 , jm1
    !
    ! THIS ROUTINE DETERMINES P(.) FROM P(X) BY A 4-POINT INTERPOLATION.
    ! ON THE X-GRID, A P(X) POINT OUTSIDE THE GRID DOMAIN IS ASSUMED TO
    ! SATISFY P(0,J) = P(1,J); P(NI,J) = P(NI-1,J); AND SIMILARLY FOR THE
    ! I'S.
    !
    if ( iband == 1 ) then
      if ( icrm == 1 ) then
        do j = 1 , nj
          do i = 1 , ni
            im1 = i-1
            jm1 = j-1
            if (im1 == 0) im1 = ni
            if (jm1 == 0) jm1 = nj
            pd(i,j) = d_rfour*(px(i,j)+px(im1,j)+px(i,jm1)+px(im1,jm1))
          end do
        end do
      else
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
      end if
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
  end subroutine crs2dot_2d

  subroutine crs2dot_3d(pd,px,ni,nj,nk,iband,icrm)
    implicit none
    integer(ik4) , intent(in) :: ni , nj , nk
    integer(ik4) , intent(in) :: iband , icrm
    real(rkx) , intent(in) , dimension(ni,nj,nk) :: px
    real(rkx) , intent(out) , dimension(ni,nj,nk) :: pd

    integer(ik4) :: i , j , k , ni1 , nj1 , im1 , jm1
    !
    ! THIS ROUTINE DETERMINES P(.) FROM P(X) BY A 4-POINT INTERPOLATION.
    ! ON THE X-GRID, A P(X) POINT OUTSIDE THE GRID DOMAIN IS ASSUMED TO
    ! SATISFY P(0,J) = P(1,J); P(NI,J) = P(NI-1,J); AND SIMILARLY FOR THE
    ! I'S.
    !
    if ( iband == 1 ) then
      if ( icrm == 1 ) then
        do k = 1 , nk
          do j = 1 , nj
            do i = 1 , ni
              im1 = i-1
              jm1 = j-1
              if (im1 == 0) im1 = ni
              if (jm1 == 0) jm1 = nj
              pd(i,j,k) = d_rfour*(px(i,j,k)+px(im1,j,k) + &
                                   px(i,jm1,k)+px(im1,jm1,k))
            end do
          end do
        end do
      else
        nj1 = nj - 1
        do k = 1 , nk
          do j = 2 , nj1
            do i = 1 , ni
              im1 = i-1
              if (im1 == 0) im1 = ni
              pd(i,j,k) = d_rfour*(px(i,j,k)+px(im1,j,k) + &
                                   px(i,j-1,k)+px(im1,j-1,k))
            end do
          end do
          do i = 1 , ni
            im1 = i-1
            if (im1 == 0) im1 = ni
            pd(i,1,k) = d_half*(px(i,1,k)+px(im1,1,k))
            pd(i,nj,k) = d_half*(px(i,nj1,k)+px(im1,nj1,k))
          end do
        end do
      end if
    else
      ni1 = ni - 1
      nj1 = nj - 1
      do k = 1 , nk
        do j = 2 , nj1
          do i = 2 , ni1
            pd(i,j,k) = d_rfour*(px(i,j,k)+px(i-1,j,k) + &
                                 px(i,j-1,k)+px(i-1,j-1,k))
          end do
        end do
        do i = 2 , ni1
          pd(i,1,k) = d_half*(px(i,1,k)+px(i-1,1,k))
          pd(i,nj,k) = d_half*(px(i,nj1,k)+px(i-1,nj1,k))
        end do
        do j = 2 , nj1
          pd(1,j,k) = d_half*(px(1,j,k)+px(1,j-1,k))
          pd(ni,j,k) = d_half*(px(ni1,j,k)+px(ni1,j-1,k))
        end do
        pd(1,1,k) = px(1,1,k)
        pd(1,nj,k) = px(1,nj1,k)
        pd(ni,1,k) = px(ni1,1,k)
        pd(ni,nj,k) = px(ni1,nj1,k)
      end do
    end if
  end subroutine crs2dot_3d

  subroutine dot2crs(px,pd,ni,nj,iband,icrm)
    implicit none
    integer(ik4) , intent(in) :: ni , nj
    real(rkx) , intent(in) , dimension(ni,nj) :: pd
    real(rkx) , intent(out) , dimension(ni,nj) :: px
    integer(ik4) , intent(in) :: iband , icrm
    integer(ik4) :: i , j , ip1 , jp1
    if ( iband == 1 ) then
      if ( icrm == 1 ) then
        do j = 1 , nj
          jp1 = j + 1
          if ( jp1 > nj ) jp1 = 1
          do i = 1 , ni
            ip1 = i + 1
            if ( ip1 > ni ) ip1 = 1
            px(i,j) = ( pd(i  ,j  ) + &
                        pd(ip1,j  ) + &
                        pd(i  ,jp1) + &
                        pd(ip1,jp1) ) * 0.25_rkx
          end do
        end do
      else
        do j = 1 , nj - 1
          do i = 1 , ni
            ip1 = i + 1
            if ( ip1 > ni ) ip1 = 1
            px(i,j) = ( pd(i  ,j  ) + &
                        pd(ip1,j  ) + &
                        pd(i  ,j+1) + &
                        pd(ip1,j+1) ) * 0.25_rkx
          end do
        end do
      end if
    else
      do j = 1 , nj - 1
        do i = 1 , ni - 1
          px(i,j) = ( pd(i  ,j  ) + &
                      pd(i+1,j  ) + &
                      pd(i  ,j+1) + &
                      pd(i+1,j+1) ) * 0.25_rkx
        end do
      end do
    end if
  end subroutine dot2crs
  !
  !-----------------------------------------------------------------------
  !
  subroutine top2btm(x)
    implicit none
    real(rkx) , pointer , intent(inout) , dimension(:,:,:) :: x
    integer(ik4) :: i1 , i2 , j1 , j2 , k1 , k2
    integer(ik4) :: i , j , k , kr
    real(rkx) , dimension(size(x,3)) :: work

    i1 = lbound(x,1)
    i2 = ubound(x,1)
    j1 = lbound(x,2)
    j2 = ubound(x,2)
    k1 = lbound(x,3)
    k2 = ubound(x,3)
    do j = j1 , j2
      do i = i1 , i2
        do k = k1 , k2
          work(k) = x(i,j,k)
        end do
        do k = k1 , k2
          kr = k2 - k + 1
          x(i,j,k) = work(kr)
        end do
      end do
    end do
  end subroutine top2btm
  !
  !-----------------------------------------------------------------------
  !
  subroutine btm2top(x)
    implicit none
    real(rkx) , pointer , intent(inout) , dimension(:,:,:) :: x
    call top2btm(x)
  end subroutine btm2top
  !
  !-----------------------------------------------------------------------
  !
  subroutine relax(chi,ff,imx,jmx,ie,je,ds)
    implicit none
    integer(ik4) , intent(in) :: imx , jmx , ie , je
    real(rkx) , intent(out) , dimension(imx,jmx) :: chi
    real(rkx) , intent(inout) , dimension(imx,jmx) :: ff
    real(rkx) , intent(in) :: ds
    real(rkx) , parameter :: smallres = 1.0e-6_rkx
    integer(ik4) , parameter :: mm = 20000
    real(rkx) , parameter :: alpha = 1.8_rkx
    real(rkx) , parameter :: alphaov4 = alpha * d_rfour
    integer(ik4) :: i , j , iter
    real(rkx) , dimension(jmx) :: chimx
    real(rkx) , dimension(jmx) :: rdmax
    real(rkx) , dimension(imx,jmx) :: rd
    real(rkx) :: epx , fac
    logical :: converged
    converged = .false.
    fac = d_two * ds * ds
    rd(:,:) = d_zero
    chi(:,:) = d_zero
    do j = 1 , je + 1
      do i = 1 , ie + 1
        ff(i,j) = fac * ff(i,j)
      end do
    end do
    iter_loop : do iter = 1 , mm
      chimx(:) = d_zero
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
      rdmax(:) = d_zero
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
    real(rkx) , intent(inout) , dimension(ix,jx) :: f
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
  subroutine meandivf(u,v,psd,dm,sigf,dsg,imx,jmx,kxs,ds,imxm,jmxm)
    implicit none
    integer(ik4) , intent(in) :: imx , jmx , kxs , imxm , jmxm
    real(rkx) , intent(inout) , dimension(imx,jmx,kxs) :: u , v
    real(rkx) , intent(in) , dimension(imx,jmx) :: psd
    real(rkx) , intent(in) , dimension(imx,jmx) :: dm
    real(rkx) , intent(in) , dimension(kxs) :: sigf
    real(rkx) , intent(in) , dimension(kxs-1) :: dsg
    real(rkx) , intent(in) :: ds
    integer(ik4) :: i , j , k
    real(rkx) , dimension(imx,jmx) :: chi , div
    real(rkx) , dimension(imx,jmx) :: dudx , dvdy
    real(rkx) , dimension(imx,jmx) :: udiverg , vdiverg
    real(rkx) , dimension(imx,jmx) :: uslb , vslb
    real(rkx) , dimension(kxs) :: weight
    real(rkx) :: oneov2ds

    oneov2ds = d_one / (d_two * ds)
    do k = 1 , kxs
      weight(k) = d_two * (d_one - sigf(k))
    end do
    !
    ! Integrate p* v/m, compute div, to dot point, to (x,y) format
    !
    do j = 1 , jmx
      do i = 1 , imx
        uslb(i,j) = d_zero
        vslb(i,j) = d_zero
      end do
    end do
    do j = 1 , jmx
      do i = 1 , imx
        do k = 1 , kxs-1
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
        dudx(i,j) = uslb(i+1,j+1) - uslb(i,j+1) + uslb(i+1,j) - uslb(i,j)
        dvdy(i,j) = vslb(i+1,j+1) - vslb(i+1,j) + vslb(i,j+1) - vslb(i,j)
      end do
    end do
    div(:,:) = d_zero
    do j = 1 , jmxm
      do i = 1 , imxm
        div(i,j) = oneov2ds * (dudx(i,j) + dvdy(i,j))
      end do
    end do
    !
    ! Iteratively solve laplacian from good first guess.
    !
    call relax(chi,div,imx,jmx,imx-2,jmx-2,ds)
    !
    ! Get divergent component of wind, 2d field on dot points.
    !
    do j = 2 , jmxm
      do i = 2 , imxm
         vdiverg(i,j) = (chi(i,j) - chi(i,j-1) + &
                         chi(i-1,j) - chi(i-1,j-1)) * oneov2ds
      end do
    end do
    do j = 2 , jmxm
      do i = 2 , imxm
        udiverg(i,j) = (chi(i,j) - chi(i-1,j) + &
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
    do j = 1 , jmx
      do i = 1 , imx
        do k = 1 , kxs
          u(i,j,k) = u(i,j,k) - weight(k) * udiverg(i,j)
          v(i,j,k) = v(i,j,k) - weight(k) * vdiverg(i,j)
        end do
      end do
    end do
  end subroutine meandivf
  !
  !-----------------------------------------------------------------------
  !
  subroutine meandiv(u,v,psd,dm,sigh,dsg,imx,jmx,kxs,ds,imxm,jmxm)
    implicit none
    integer(ik4) , intent(in) :: imx , jmx , kxs , imxm , jmxm
    real(rkx) , intent(inout) , dimension(imx,jmx,kxs) :: u , v
    real(rkx) , intent(in) , dimension(imx,jmx) :: psd
    real(rkx) , intent(in) , dimension(imx,jmx) :: dm
    real(rkx) , intent(in) , dimension(kxs) :: sigh , dsg
    real(rkx) , intent(in) :: ds
    integer(ik4) :: i , j , k
    real(rkx) , dimension(imx,jmx) :: chi , div
    real(rkx) , dimension(imx,jmx) :: dudx , dvdy
    real(rkx) , dimension(imx,jmx) :: udiverg , vdiverg
    real(rkx) , dimension(imx,jmx) :: uslb , vslb
    real(rkx) , dimension(kxs) :: weight
    real(rkx) :: oneov2ds

    oneov2ds = d_one / (d_two * ds)
    do k = 1 , kxs
      weight(k) = d_two * (d_one - sigh(k))
    end do
    !
    ! Integrate p* v/m, compute div, to dot point, to (x,y) format
    !
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
        dudx(i,j) = uslb(i+1,j+1) - uslb(i,j+1) + uslb(i+1,j) - uslb(i,j)
        dvdy(i,j) = vslb(i+1,j+1) - vslb(i+1,j) + vslb(i,j+1) - vslb(i,j)
      end do
    end do
    div(:,:) = d_zero
    do j = 1 , jmxm
      do i = 1 , imxm
        div(i,j) = oneov2ds * (dudx(i,j) + dvdy(i,j))
      end do
    end do
    !
    ! Iteratively solve laplacian from good first guess.
    !
    call relax(chi,div,imx,jmx,imx-2,jmx-2,ds)
    !
    ! Get divergent component of wind, 2d field on dot points.
    !
    do j = 2 , jmxm
      do i = 2 , imxm
         vdiverg(i,j) = (chi(i,j) - chi(i,j-1) + &
                         chi(i-1,j) - chi(i-1,j-1)) * oneov2ds
      end do
    end do
    do j = 2 , jmxm
      do i = 2 , imxm
        udiverg(i,j) = (chi(i,j) - chi(i-1,j) + &
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
