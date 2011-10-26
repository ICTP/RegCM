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
   
module mod_advection
!
! Horizontal and vertical advection.
!
  use mod_runparams
  use mod_memutil
  use mod_service

  private
   
  public :: init_advection, hadv , vadv

  real(8) , pointer , dimension(:,:,:) :: u    ! U wind * ps
  real(8) , pointer , dimension(:,:,:) :: v    ! V wind * ps
  real(8) , pointer , dimension(:,:) :: ps     ! Surface pressure
  real(8) , pointer , dimension(:,:) :: mapfx  ! Map factor Cross
  real(8) , pointer , dimension(:,:) :: mapfd  ! Map factor Dot
  real(8) , pointer , dimension(:,:,:) :: vsv  ! Vertical Sigma Velocity
  integer , pointer , dimension(:,:) :: kpbl   ! Top of PBL
!
! working space used to store the interlated values in vadv.
!
  real(8) , pointer , dimension(:,:) :: fg

  real(8) , parameter :: c287 = 0.287D+00
!
! relaxed upstream scheme factors
!
  real(8) , parameter :: fact1 = 0.60D0
!hy
! real(8) , parameter :: fact1 = 0.75D0
!hy
  real(8) , parameter :: fact2 = d_one - fact1
  real(8) , parameter :: falow = 1.0D-15
!
  contains

    subroutine init_advection(dom,sps,atm,vertvel,kpbltop)
      use mod_atm_interface , only : atmstate , domain , surfpstate
      implicit none
      type(domain) , intent(in) :: dom
      type(surfpstate), intent(in) :: sps
      type(atmstate) , intent(in) :: atm
      real(8) , pointer , dimension(:,:,:) :: vertvel
      integer , pointer , dimension(:,:) :: kpbltop

      call assignpnt(atm%u,u)
      call assignpnt(atm%v,v)
      call assignpnt(sps%ps,ps)
      call assignpnt(dom%msfx,mapfx)
      call assignpnt(dom%msfd,mapfd)
      call assignpnt(vertvel,vsv)
      call assignpnt(kpbltop,kpbl)

      call getmem2d(fg,lbound(atm%t,1),ubound(atm%t,1), &
                       lbound(atm%t,1),ubound(atm%t,2),'mod_advection:fg')
    end subroutine init_advection
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!  HADV                                                               c
!                                                                     c
!     This subroutines computes the horizontal flux-divergence terms. c
!     second-order difference is used.                                c
!                                                                     c
!     ldot   : cross/dot variable flagg                               c
!                                                                     c
!     ften   : is the tendency for variable 'f'.                      c
!                                                                     c
!     f      : is p*f.                                                c
!                                                                     c
!     jstart : is the j'th slice of f anf ften to start               c
!                                                                     c
!     jsstop : is the j'th slice of f anf ften to stop                c
!                                                                     c
!     ind = 1 : for t and qv.                                         c
!         = 2 : for qc and qr.                                        c
!         = 3 : for u and v.                                          c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine hadv(ldot,ften,f,jstart,jstop,ind)
!
      implicit none
!
      logical , intent(in) :: ldot ! Cross/dot flag
      integer , intent (in) :: ind , jstart , jstop
      real(8) , pointer , intent (in) , dimension(:,:,:) :: f
      real(8) , pointer , intent (inout), dimension(:,:,:) :: ften
!
      real(8) :: dxx , fx1 , fx2 , fy1 , fy2
      real(8) :: ucmona , ucmonb , ucmonc , vcmona , vcmonb , vcmonc
      integer :: i , j , k
      integer :: istart , istopd , istopx , kstart , kstop
      integer :: idx , idxm1 , idxp1 , jm1 , jdm1 , jp1 , jdp1
!
      character (len=64) :: subroutine_name='hadv'
      integer :: idindx=0
!
      call time_begin(subroutine_name,idindx)
!
      istart = lbound(ften,1) + 1
      istopd = ubound(ften,1) - 1
      istopx = ubound(ften,1) - 2
      kstart = lbound(ften,2)
      kstop  = ubound(ften,2)

      do j = jstart , jstop
        jm1 = j - 1
        jp1 = j + 1
!
!
!       ua, va : are p*u and p*v.
!       msfd   : is the map scale factor at dot points.
!
!
        if ( ldot ) then
          if ( ind == 3 ) then
!
!           for u and v:
!
            dxx = dx16
#ifdef BAND
            jdp1 = jp1
            jdm1 = jm1
#else
            jdp1 = j + 1
            jdm1 = j - 1
            if ( myid == 0 ) jdm1 = max0(jdm1,2)
            if ( myid == nproc-1 ) jdp1 = min0(jdp1,jendl-1)
#endif
!
            do k = kstart , kstop
              do i = istart , istopd
                idx = i
                idxp1 = i + 1
                idxp1 = min0(idxp1,istopd)
                idxm1 = i - 1
                idxm1 = max0(idxm1,istart)
                ucmona = u(idxp1,k,j)+d_two*u(idx,k,j)+u(idxm1,k,j)
                vcmona = v(idx,k,jdp1)+d_two*v(idx,k,j)+v(idx,k,jdm1)
                ucmonb = u(idxp1,k,jdp1) + d_two*u(idx,k,jdp1) + &
                         u(idxm1,k,jdp1) + ucmona
                vcmonb = v(idxp1,k,jdp1) + d_two*v(idxp1,k,j) +  &
                         v(idxp1,k,jdm1) + vcmona
                ucmonc = u(idxp1,k,jdm1) + d_two*u(idx,k,jdm1) + &
                         u(idxm1,k,jdm1) + ucmona
                vcmonc = v(idxm1,k,jdp1) + d_two*v(idxm1,k,j) +  &
                         v(idxm1,k,jdm1) + vcmona
                ften(i,k,j) = ften(i,k,j) -                  &
                            ((f(i,k,jp1)+f(i,k,j))*ucmonb -  &
                             (f(i,k,j)+f(i,k,jm1))*ucmonc +  &
                             (f(i+1,k,j)+f(i,k,j))*vcmonb -  &
                             (f(i,k,j)+f(i-1,k,j))*vcmonc) / &
                             (dxx*mapfd(i,j)*mapfd(i,j))
              end do
            end do
          else
            call fatal(__FILE__,__LINE__, &
                       'The advection scheme you required is not available.')
          end if

        else  ! This part is for cross point defined variables

          if ( ind == 1 ) then
            dxx = dx4
!
!           for t and qv:
!
            do k = kstart , kstop
              do i = istart , istopx
                ften(i,k,j) = ften(i,k,j) -                             &
                    ((u(i+1,k,jp1)+u(i,k,jp1))*(f(i,k,jp1)+f(i,k,j)) -  &
                     (u(i+1,k,j)+u(i,k,j)) *   (f(i,k,j)+f(i,k,jm1)) +  &
                     (v(i+1,k,jp1)+v(i+1,k,j))*(f(i+1,k,j)+f(i,k,j)) -  &
                     (v(i,k,jp1)+v(i,k,j)) *   (f(i-1,k,j)+f(i,k,j))) / &
                     (dxx*mapfx(i,j)*mapfx(i,j))
              end do
            end do
!
          else if ( ind == 2 ) then
            dxx = dx
!
!           for qc and qr:
!           up-wind values of qc and qr are used.
!
            do k = kstart , kstop
              do i = istart , istopx
                ucmonb = d_half*(u(i+1,k,jp1)+u(i,k,jp1))
                ucmona = d_half*(u(i+1,k,j)+u(i,k,j))
                if ( ucmonb >= d_zero ) then
                  fx2 = fact1*f(i,k,j) + fact2*f(i,k,jp1)
                else
                  fx2 = fact1*f(i,k,jp1) + fact2*f(i,k,j)
                end if
                if ( ucmona >= d_zero ) then
                  fx1 = fact1*f(i,k,jm1) + fact2*f(i,k,j)
                else
                  fx1 = fact1*f(i,k,j) + fact2*f(i,k,jm1)
                end if
                vcmonb = d_half*(v(i+1,k,jp1)+v(i+1,k,j))
                vcmona = d_half*(v(i,k,jp1)+v(i,k,j))
                if ( vcmonb >= d_zero ) then
                  fy2 = fact1*f(i,k,j) + fact2*f(i+1,k,j)
                else
                  fy2 = fact1*f(i+1,k,j) + fact2*f(i,k,j)
                end if
                if ( vcmona >= d_zero ) then
                  fy1 = fact1*f(i-1,k,j) + fact2*f(i,k,j)
                else
                  fy1 = fact1*f(i,k,j) + fact2*f(i-1,k,j)
                end if
                ften(i,k,j) = ften(i,k,j) -                          &
                     (ucmonb*fx2-ucmona*fx1+vcmonb*fy2-vcmona*fy1) / &
                     (dxx*mapfx(i,j)*mapfx(i,j))
              end do
            end do
          else
            call fatal(__FILE__,__LINE__, &
                       'The advection scheme you required is not available.')
          end if
        end if
      end do
!
      call time_end(subroutine_name,idindx)
!
    end subroutine hadv
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
! VADV                                                                c
!                                                                     c
!     This subroutine computes the vertical flux-divergence terms.    c
!                                                                     c
!     ften   : is the tendency of variable 'f'.                       c
!                                                                     c
!     f      : is p*f.                                                c
!                                                                     c
!     jstart : is the j'th slice of f anf ften to start               c
!                                                                     c
!     jsstop : is the j'th slice of f anf ften to stop                c
!                                                                     c
!     ind = 1 : for t.                                                c
!           2 : for qv.                                               c
!           3 : for qc and qr.                                        c
!           4 : for u and v.                                          c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine vadv(ften,f,jstart,jstop,ind)
!
      implicit none
!
      integer , intent(in) :: ind , jstart , jstop
      real(8) , pointer , intent (in) , dimension(:,:,:) :: f
      real(8) , pointer , intent (inout), dimension(:,:,:) :: ften
!
      real(8) :: f1 , f2 , slope
      integer :: i , j , k
      integer :: istart , istopd , istopx , kstart , kstop
!
      character (len=64) :: subroutine_name='vadv'
      integer :: idindx=0
!
!----------------------------------------------------------------------
!
      call time_begin(subroutine_name,idindx)
      istart = lbound(f,1) + 1
      istopd = ubound(f,1) - 1
      istopx = ubound(f,1) - 2
      kstart = lbound(f,2)
      kstop  = ubound(f,2)
!
      do j = jstart , jstop
!
        if ( ind == 1 ) then
!
!         vertical advection terms : interpolate t to full sigma levels
!
          do i = istart , istopx
            fg(i,kstart) = d_zero
          end do
          do k = kstart+1 , kstop
            do i = istart , istopx
              fg(i,k) = twt(k,1)*f(i,k,j) *                              &
                   ((ps(i,j)*sigma(k)+ptop)/(ps(i,j)*a(k)+ptop))**c287 + &
                        twt(k,2)*f(i,k-1,j) *                            &
                   ((ps(i,j)*sigma(k)+ptop)/(ps(i,j)*a(k-1)+ptop))**c287
            end do
          end do
          do i = istart , istopx
            ften(i,kstart,j) = ften(i,kstart,j) - &
                     vsv(i,kstart+1,j)*fg(i,kstart+1)/dsigma(kstart)
          end do
          do k = kstart+1 , kstop-1
            do i = istart , istopx
              ften(i,k,j) = ften(i,k,j) - &
                   (vsv(i,k+1,j)*fg(i,k+1)-vsv(i,k,j)*fg(i,k))/dsigma(k)
            end do
          end do
          do i = istart , istopx
            ften(i,kstop,j) = ften(i,kstop,j) + &
                     vsv(i,kstop,j)*fg(i,kstop)/dsigma(kstop)
          end do
!
        else if ( ind == 2 ) then
!
!         vertical advection term : interpolate qv to full sigma levels
!
          do i = istart , istopx
            fg(i,kstart) = d_zero
          end do
          do k = kstart+1 , kstop
            do i = istart , istopx
              if ( f(i,k,j) > falow .and. f(i,k-1,j) > falow ) then
                fg(i,k) = f(i,k,j)*(f(i,k-1,j)/f(i,k,j))**qcon(k)
              else
                fg(i,k) = d_zero
              end if
            end do
          end do
          do i = istart , istopx
            ften(i,kstart,j) = ften(i,kstart,j) - &
                  vsv(i,kstart+1,j)*fg(i,kstart+1)/dsigma(kstart)
          end do
          do k = kstart+1 , kstop-1
            do i = istart , istopx
              ften(i,k,j) = ften(i,k,j) - &
                       (vsv(i,k+1,j)*fg(i,k+1)-vsv(i,k,j)*fg(i,k))/dsigma(k)
            end do
          end do
          do i = istart , istopx
            ften(i,kstop,j) = ften(i,kstop,j) + &
                      vsv(i,kstop,j)*fg(i,kstop)/dsigma(kstop)
          end do
!
        else if ( ind == 3 ) then
!
!         vertical advection terms for qc and qr:
!
          do i = istart , istopx
            if ( vsv(i,kstart+1,j) >= d_zero ) then
              f2 = f(i,kstart,j)
            else
              f2 = f(i,kstart+1,j)
            end if
            ften(i,kstart,j) = ften(i,kstart,j) - &
                         vsv(i,kstart+1,j)*f2/dsigma(kstart)
          end do
          do k = kstart+1 , kstop-1
            do i = istart , istopx
              if ( vsv(i,k+1,j) >= d_zero ) then
                f2 = f(i,k,j)
              else
                f2 = f(i,k+1,j)
              end if
              if ( vsv(i,k,j) >= d_zero ) then
                f1 = f(i,k-1,j)
              else
                f1 = f(i,k,j)
              end if
              ften(i,k,j) = ften(i,k,j) - &
                    (vsv(i,k+1,j)*f2-vsv(i,k,j)*f1)/dsigma(k)
            end do
          end do
          do i = istart , istopx
            if ( vsv(i,kstop,j) >= d_zero ) then
              f1 = f(i,kstop-1,j)
            else
              f1 = f(i,kstop,j)
            end if
            ften(i,kstop,j) = ften(i,kstop,j) + vsv(i,kstop,j)*f1/dsigma(kstop)
          end do
!
        else if ( ind == 4 ) then
!
!         vertical advection terms : interpolate ua or va to full sigma levels
!
          do i = istart , istopd
            fg(i,kstart) = d_zero
          end do
          do k = kstart+1 , kstop
            do i = istart , istopd
              fg(i,k) = d_half*(f(i,k,j)+f(i,k-1,j))/mapfd(i,j)
            end do
          end do
          do i = istart , istopd
            ften(i,kstart,j) = ften(i,kstart,j) -                    &
                        (vsv(i-1,kstart+1,j-1)+vsv(i,kstart+1,j-1) + &
                         vsv(i,kstart+1,j)+vsv(i-1,kstart+1,j))    * &
                         fg(i,kstart+1)/(d_four*dsigma(kstart))
          end do
          do k = kstart+1 , kstop-1
            do i = istart , istopd
              ften(i,k,j) = ften(i,k,j) -                               &
                          ((vsv(i,k+1,j-1)+vsv(i-1,k+1,j-1)+            &
                            vsv(i,k+1,j)  +vsv(i-1,k+1,j))*fg(i,k+1) -  &
                           (vsv(i,k,j-1)  +vsv(i-1,k,j-1)+              &
                            vsv(i,k,j)    +vsv(i-1,k,j))*fg(i,k)) /     &
                          (d_four*dsigma(k))
            end do
          end do
          do i = istart , istopd
            ften(i,kstop,j) = ften(i,kstop,j) +                 &
                         (vsv(i,kstop,j-1)+vsv(i-1,kstop,j-1) + &
                          vsv(i,kstop,j)+vsv(i-1,kstop,j)) *    &
                          fg(i,kstop)/(d_four*dsigma(kstop))
          end do
!
        else if ( ind == 5 ) then
       
          do k = kstart+1 , kstop
            do i = istart , istopx
              fg(i,k) = twt(k,1)*f(i,k,j) + twt(k,2)*f(i,k-1,j)
            end do
          end do
          do i = istart , istopx
            ften(i,kstart,j) = ften(i,kstart,j) - &
                     vsv(i,kstart+1,j)*fg(i,kstart+1)/dsigma(kstart)
          end do
          do k = kstart+1 , kstop-1
            do i = istart , istopx
              ften(i,k,j) = ften(i,k,j) - &
                      (vsv(i,k+1,j)*fg(i,k+1)-vsv(i,k,j)*fg(i,k))/dsigma(k)
            end do
          end do
          do i = istart , istopx
            ften(i,kstop,j) = ften(i,kstop,j) + &
                    vsv(i,kstop,j)*fg(i,kstop)/dsigma(kstop)
          end do

        else if ( ind == 6 ) then

          do k = kstart+1 , kstop
            do i = istart , istopd
              fg(i,k)= twt(k,1)*f(i,k,j) + twt(k,2)*f(i,k-1,j)
            end do
          end do

          do i = istart , istopd
            if ( kpbl(j,i).gt.kstop ) then
              call fatal(__FILE__,__LINE__,'kpbl is greater than kstop')
            end if
            if ( kpbl(j,i).ge.4 ) then
              ! Calculate slope of scalar in layer above ambiguous layer
              k = kpbl(j,i)-2
              if ( (f(i,k+1,j)-f(i,k,j)) > d_zero .and. &
                   (f(i,k,j)-f(i,k-1,j)) > d_zero ) then
                slope = min((f(i,k+1,j)-f(i,k,j))/(a(k+1)-a(k)), &
                            (f(i,k,j)-f(i,k-1,j))/(a(k)-a(k-1)))
              else if ( (f(i,k+1,j)-f(i,k,j)) < d_zero .and. &
                        (f(i,k,j)-f(i,k-1,j)) < d_zero ) then
                slope = max((f(i,k+1,j)-f(i,k,j))/(a(k+1)-a(k)), &
                            (f(i,k,j)-f(i,k-1,j))/(a(k)-a(k-1)))
              else
                slope = d_zero
              end if
              ! Now replace the values of scalar at top and bottom of ambiguous
              ! layer as long as inversion is actually in the ambiguous layer
              k = kpbl(j,i)
              fg(i,k-1) = f(i,k-2,j) + slope*(sigma(k-1)-a(k-2))
              if (abs(f(i,k-2,j) + slope*(a(k-1)-a(k-2))-f(i,k,j)) > &
                  abs(f(i,k-1,j)-f(i,k,j)) ) then
                fg(i,k) = f(i,k,j)
              else
                fg(i,k) = f(i,k-2,j) + slope*(sigma(k)-a(k-2))
              end if
            end if
          end do

          do i = istart , istopd
            ften(i,kstart,j) = ften(i,kstart,j) - &
                   vsv(i,kstart+1,j)*fg(i,kstart+1)/dsigma(kstart)
          end do
          do k = kstart+1 , kstop-1
            do i = istart , istopd
              ften(i,k,j) = ften(i,k,j) - &
                  (vsv(i,k+1,j)*fg(i,k+1)-vsv(i,k,j)*fg(i,k))/dsigma(k)
            end do
          end do
          do i = istart , istopd
            ften(i,kstop,j) = ften(i,kstop,j) + &
                   vsv(i,kstop,j)*fg(i,kstop)/dsigma(kstop)
          end do

        end if
!
      end do
!
      call time_end(subroutine_name,idindx)

    end subroutine vadv
!
end module mod_advection
