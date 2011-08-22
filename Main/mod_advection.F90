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

    subroutine init_advection(dom,sps,atm,vertvel)
      use mod_atm_interface , only : atmstate , domain , surfpstate
      implicit none
      type(domain) , intent(in) :: dom
      type(surfpstate), intent(in) :: sps
      type(atmstate) , intent(in) :: atm
      real(8) , pointer , dimension(:,:,:) :: vertvel
      u     => atm%u
      v     => atm%v
      ps    => sps%ps
      mapfx => dom%msfx
      mapfd => dom%msfd
      vsv   => vertvel
      call getmem2d(fg,1,iy,1,kz,'mod_advection:fg')
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
!     j      : is the j'th slice of f anf ften.                       c
!                                                                     c
!     ind = 1 : for t and qv.                                         c
!         = 2 : for qc and qr.                                        c
!         = 3 : for u and v.                                          c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
    subroutine hadv(ldot,ften,f,j,ind)
!
    implicit none
!
    logical , intent(in) :: ldot ! Cross/dot flag
    integer , intent (in) :: ind , j
    real(8) , intent (in) , dimension(iy,kz,0:jxp+1) :: f
    real(8) , intent (inout), dimension(iy,kz,jxp) :: ften
!
    integer :: jm1 , jp1
    real(8) :: dxx , fx1 , fx2 , fy1 , fy2
    real(8) :: ucmona , ucmonb , ucmonc , vcmona , vcmonb , vcmonc
    integer :: idx , idxm1 , idxp1 , jdm1 , jdp1
    integer :: i , k
    character (len=64) :: subroutine_name='hadv'
    integer :: idindx=0
!
    call time_begin(subroutine_name,idindx)
!
    jm1 = j - 1
    jp1 = j + 1
!
!----------------------------------------------------------------------
!
! ua, va : are p*u and p*v.
! msfx   : is the map scale factor at cross points.
!
!----------------------------------------------------------------------
!
    if ( ldot ) then
      if ( ind == 3 ) then
        dxx = dx16
!
!-----for u and v:
!
#ifdef BAND
        do k = 1 , kz
          do i = 2 , iym1
            idx = i
            idxp1 = i + 1
            idxp1 = min0(idxp1,iym1)
            idxm1 = i - 1
            idxm1 = max0(idxm1,2)
            ucmona = u(idxp1,k,j)+d_two*u(idx,k,j)+u(idxm1,k,j)
            vcmona = v(idx,k,jp1)+d_two*v(idx,k,j)+v(idx,k,jm1)
            ucmonb = u(idxp1,k,jp1)+d_two*u(idx,k,jp1)+u(idxm1,k,jp1) + ucmona
            vcmonb = v(idxp1,k,jp1)+d_two*v(idxp1,k,j)+v(idxp1,k,jm1) + vcmona
            ucmonc = u(idxp1,k,jm1)+d_two*u(idx,k,jm1)+u(idxm1,k,jm1) + ucmona
            vcmonc = v(idxm1,k,jp1)+d_two*v(idxm1,k,j)+v(idxm1,k,jm1) + vcmona
            ften(i,k,j) = ften(i,k,j) -                  &
                        ((f(i,k,jp1)+f(i,k,j))*ucmonb -  &
                         (f(i,k,j)+f(i,k,jm1))*ucmonc +  &
                         (f(i+1,k,j)+f(i,k,j))*vcmonb -  &
                         (f(i,k,j)+f(i-1,k,j))*vcmonc) / &
                         (dxx*mapfd(i,j)*mapfd(i,j))
          end do
        end do
!
!----------------------------------------------------------------------
#else
        jdp1 = j + 1
        jdm1 = j - 1
        if ( myid == 0 ) jdm1 = max0(jdm1,2)
        if ( myid == nproc-1 ) jdp1 = min0(jdp1,jendl-1)
!
        do k = 1 , kz
          do i = 2 , iym1
            idx = i
            idxp1 = i + 1
            idxp1 = min0(idxp1,iym1)
            idxm1 = i - 1
            idxm1 = max0(idxm1,2)
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
#endif
!
      else
        call fatal(__FILE__,__LINE__, &
             'The advection scheme you required is not available.')
      end if
    else
      if ( ind == 1 ) then
        dxx = dx4
!
!-----for t and qv:
!
        do k = 1 , kz
          do i = 2 , iym2
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
!-----for qc and qr:
!       up-wind values of qc and qr are used.
!
        do k = 1 , kz
          do i = 2 , iym2
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
    call time_end(subroutine_name,idindx)
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
!     fa     : is p*f.                                                c
!                                                                     c
!     j      : j'th slice of variable fa.                             c
!                                                                     c
!     ind = 1 : for t.                                                c
!           2 : for qv.                                               c
!           3 : for qc and qr.                                        c
!           4 : for u and v.                                          c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
subroutine vadv(ften,f,j,ind,kpbl1d)
!
    implicit none
!
    integer , intent(in) :: ind , j
    integer , intent(in) , dimension(iy) :: kpbl1d
    real(8) , intent (in) , dimension(iy,kz,-1:jxp+2) :: f
    real(8) , intent (inout), dimension(iy,kz,jxp) :: ften
!
    real(8) :: f1 , f2 , slope
    integer :: i , k
!
    character (len=64) :: subroutine_name='vadv'
    integer :: idindx=0
!
#ifdef BAND
    integer :: jm1
!----------------------------------------------------------------------
!----------------------------------------------------------------------
    call time_begin(subroutine_name,idindx)
    jm1 = j-1
!
    if ( ind == 1 ) then
!
!-----vertical advection terms for:
!.....interpolate ta to full sigma levels:
!
      do i = 2 , iym2
        fg(i,1) = d_zero
      end do
      do k = 2 , kz
        do i = 2 , iym2
          fg(i,k) = twt(k,1)*f(i,k,j)   * ((ps(i,j)*sigma(k)+r8pt)/  &
                     (ps(i,j)*a(k)+r8pt))**c287 +                   &
                    twt(k,2)*f(i,k-1,j) * ((ps(i,j)*sigma(k)+r8pt) / &
                     (ps(i,j)*a(k-1)+r8pt))**c287
        end do
      end do
!......k = 1
      do i = 2 , iym2
        ften(i,1,j) = ften(i,1,j) - vsv(i,2,j)*fg(i,2)/dsigma(1)
      end do
!......k = 2,kzm1
      do k = 2 , kzm1
        do i = 2 , iym2
          ften(i,k,j) = ften(i,k,j)- &
                (vsv(i,k+1,j)*fg(i,k+1)-vsv(i,k,j)*fg(i,k))/dsigma(k)
        end do
      end do
!,.....k = kz
      do i = 2 , iym2
        ften(i,kz,j) = ften(i,kz,j) + vsv(i,kz,j)*fg(i,kz)/dsigma(kz)
      end do
!
    else if ( ind == 2 ) then
!
!-----vertical advection term for qv:
!.....interpolate qv to full sigma levels:
!
      do i = 2 , iym2
        fg(i,1) = d_zero
      end do
      do k = 2 , kz
        do i = 2 , iym2
! modif !!
          if ( f(i,k,j) > falow .and. f(i,k-1,j) > falow ) then
            fg(i,k) = f(i,k,j)*(f(i,k-1,j)/f(i,k,j))**qcon(k)
          else
            fg(i,k) = d_zero
          end if
        end do
      end do
!......k = 1
      do i = 2 , iym2
        ften(i,1,j) = ften(i,1,j) - vsv(i,2,j)*fg(i,2)/dsigma(1)
      end do
!......k = 2,kzm1
      do k = 2 , kzm1
        do i = 2 , iym2
          ften(i,k,j) = ften(i,k,j)- &
                (vsv(i,k+1,j)*fg(i,k+1)-vsv(i,k,j)*fg(i,k))/dsigma(k)
        end do
      end do
!,.....k = kz
      do i = 2 , iym2
        ften(i,kz,j) = ften(i,kz,j) + vsv(i,kz,j)*fg(i,kz)/dsigma(kz)
      end do
!
    else if ( ind == 3 ) then
!
!-----vertical advection terms for qc and qr:
!
!......k = 1
      do i = 2 , iym2
        if ( vsv(i,2,j) >= d_zero ) then
          f2 = f(i,1,j)
        else
          f2 = f(i,2,j)
        end if
        ften(i,1,j) = ften(i,1,j) - vsv(i,2,j)*f2/dsigma(1)
      end do
!......k = 2,kzm1
      do k = 2 , kzm1
        do i = 2 , iym2
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
          ften(i,k,j) = ften(i,k,j)- &
                 (vsv(i,k+1,j)*f2-vsv(i,k,j)*f1)/dsigma(k)
        end do
      end do
!......k = kz
      do i = 2 , iym2
        if ( vsv(i,kz,j) >= d_zero ) then
          f1 = f(i,kzm1,j)
        else
          f1 = f(i,kz,j)
        end if
        ften(i,kz,j) = ften(i,kz,j) + vsv(i,kz,j)*f1/dsigma(kz)
      end do
!
    else if ( ind == 4 ) then
!
!-----vertical advection terms for u and v:
!.....interpolate ua or va to full sigma levels:
!
      do i = 2 , iym1
        fg(i,1) = d_zero
      end do
      do k = 2 , kz
        do i = 2 , iym1
          fg(i,k) = d_half*(f(i,k,j)+f(i,k-1,j))/mapfd(i,j)
        end do
      end do
!......k = 1
      do i = 2 , iym1
        ften(i,1,j) = ften(i,1,j) -                                 &
                   (vsv(i-1,2,jm1)+vsv(i,2,jm1)+vsv(i,2,j) + &
                    vsv(i-1,2,j))*fg(i,2)/(d_four*dsigma(1))
      end do
!......k = 2,kzm1
      do k = 2 , kzm1
        do i = 2 , iym1
          ften(i,k,j) = ften(i,k,j) -                                     &
                      ((vsv(i,k+1,jm1)+vsv(i-1,k+1,jm1)+            &
                        vsv(i,k+1,j)  +vsv(i-1,k+1,j))*fg(i,k+1) -  &
                       (vsv(i,k,jm1)  +vsv(i-1,k,jm1)+              &
                        vsv(i,k,j)    +vsv(i-1,k,j))*fg(i,k)) /     &
                      (d_four*dsigma(k))
        end do
      end do
!......k = kz
      do i = 2 , iym1
        ften(i,kz,j) = ften(i,kz,j) +                                   &
                    (vsv(i,kz,jm1)+vsv(i-1,kz,jm1)+vsv(i,kz,j) + &
                     vsv(i-1,kz,j))*fg(i,kz)/(d_four*dsigma(kz))
      end do
!
   
   
    else if ( ind == 5 ) then
   
      do k = 2 , kz
        do i = 2 , iym2
          fg(i,k) = twt(k,1)*f(i,k,j) + twt(k,2)*f(i,k-1,j)
        end do
      end do
   
!......k = 1
      do i = 2 , iym2
        ften(i,1,j) = ften(i,1,j) - vsv(i,2,j)*fg(i,2)/dsigma(1)
      end do
!......k = 2,kzm1
      do k = 2 , kzm1
        do i = 2 , iym2
          ften(i,k,j) = ften(i,k,j)- &
             (vsv(i,k+1,j)*fg(i,k+1)-vsv(i,k,j)*fg(i,k))/dsigma(k)
        end do
      end do
!,.....k = kz
      do i = 2 , iym2
        ften(i,kz,j) = ften(i,kz,j) + vsv(i,kz,j)*fg(i,kz)/dsigma(kz)
      end do
   
    end if
!
#else
!----------------------------------------------------------------------
    call time_begin(subroutine_name,idindx)
!
    if ( ind == 1 ) then
!
!-----vertical advection terms for:
!.....interpolate ta to full sigma levels:
!
      do i = 2 , iym2
        fg(i,1) = d_zero
      end do
      do k = 2 , kz
        do i = 2 , iym2
          fg(i,k) = twt(k,1)*f(i,k,j) *                 &
                    ((ps(i,j)*sigma(k)+r8pt)/     &
                     (ps(i,j)*a(k)+r8pt))**c287 + &
                    twt(k,2)*f(i,k-1,j) *               &
                    ((ps(i,j)*sigma(k)+r8pt)/     &
                     (ps(i,j)*a(k-1)+r8pt))**c287
        end do
      end do
!......k = 1
      do i = 2 , iym2
        ften(i,1,j) = ften(i,1,j) - vsv(i,2,j)*fg(i,2)/dsigma(1)
      end do
!......k = 2,kzm1
      do k = 2 , kzm1
        do i = 2 , iym2
          ften(i,k,j) = ften(i,k,j) - &
               (vsv(i,k+1,j)*fg(i,k+1)-vsv(i,k,j)*fg(i,k))/dsigma(k)
        end do
      end do
!,.....k = kz
      do i = 2 , iym2
        ften(i,kz,j) = ften(i,kz,j) + vsv(i,kz,j)*fg(i,kz)/dsigma(kz)
      end do
!
    else if ( ind == 2 ) then
!
!-----vertical advection term for qv:
!.....interpolate qv to full sigma levels:
!
      do i = 2 , iym2
        fg(i,1) = d_zero
      end do
      do k = 2 , kz
        do i = 2 , iym2
! modif !!
          if ( f(i,k,j) > falow .and. f(i,k-1,j) > falow ) then
            fg(i,k) = f(i,k,j)*(f(i,k-1,j)/f(i,k,j))**qcon(k)
          else
            fg(i,k) = d_zero
          end if
        end do
      end do
!......k = 1
      do i = 2 , iym2
        ften(i,1,j) = ften(i,1,j) - vsv(i,2,j)*fg(i,2)/dsigma(1)
      end do
!......k = 2,kzm1
      do k = 2 , kzm1
        do i = 2 , iym2
          ften(i,k,j) = ften(i,k,j) - &
                   (vsv(i,k+1,j)*fg(i,k+1)-vsv(i,k,j)*fg(i,k))/dsigma(k)
        end do
      end do
!,.....k = kz
      do i = 2 , iym2
        ften(i,kz,j) = ften(i,kz,j) + vsv(i,kz,j)*fg(i,kz)/dsigma(kz)
      end do
!
    else if ( ind == 3 ) then
!
!-----vertical advection terms for qc and qr:
!
!......k = 1
      do i = 2 , iym2
        if ( vsv(i,2,j) >= d_zero ) then
          f2 = f(i,1,j)
        else
          f2 = f(i,2,j)
        end if
        ften(i,1,j) = ften(i,1,j) - vsv(i,2,j)*f2/dsigma(1)
      end do
!......k = 2,kzm1
      do k = 2 , kzm1
        do i = 2 , iym2
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
!......k = kz
      do i = 2 , iym2
        if ( vsv(i,kz,j) >= d_zero ) then
          f1 = f(i,kzm1,j)
        else
          f1 = f(i,kz,j)
        end if
        ften(i,kz,j) = ften(i,kz,j) + vsv(i,kz,j)*f1/dsigma(kz)
      end do
!
    else if ( ind == 4 ) then
!
!-----vertical advection terms for u and v:
!.....interpolate ua or va to full sigma levels:
!
      do i = 2 , iym1
        fg(i,1) = d_zero
      end do
      do k = 2 , kz
        do i = 2 , iym1
          fg(i,k) = d_half*(f(i,k,j)+f(i,k-1,j))/mapfd(i,j)
        end do
      end do
!......k = 1
      do i = 2 , iym1
        ften(i,1,j) = ften(i,1,j) -                                  &
                    (vsv(i-1,2,j-1)+vsv(i,2,j-1)+vsv(i,2,j) + &
                     vsv(i-1,2,j))*fg(i,2)/(d_four*dsigma(1))
      end do
!......k = 2,kzm1
      do k = 2 , kzm1
        do i = 2 , iym1
          ften(i,k,j) = ften(i,k,j) -                                     &
                      ((vsv(i,k+1,j-1)+vsv(i-1,k+1,j-1)+            &
                        vsv(i,k+1,j)  +vsv(i-1,k+1,j))*fg(i,k+1) -  &
                       (vsv(i,k,j-1)  +vsv(i-1,k,j-1)+              &
                        vsv(i,k,j)    +vsv(i-1,k,j))*fg(i,k)) /     &
                      (d_four*dsigma(k))
        end do
      end do
!......k = kz
      do i = 2 , iym1
        ften(i,kz,j) = ften(i,kz,j) +                                    &
                     (vsv(i,kz,j-1)+vsv(i-1,kz,j-1)+vsv(i,kz,j) + &
                      vsv(i-1,kz,j))*fg(i,kz)/(d_four*dsigma(kz))
      end do
!
    else if ( ind == 5 ) then
   
      do k = 2 , kz
        do i = 2 , iym2
          fg(i,k) = twt(k,1)*f(i,k,j) + twt(k,2)*f(i,k-1,j)
        end do
      end do
   
!......k = 1
      do i = 2 , iym2
        ften(i,1,j) = ften(i,1,j) - vsv(i,2,j)*fg(i,2)/dsigma(1)
      end do
!......k = 2,kzm1
      do k = 2 , kzm1
        do i = 2 , iym2
          ften(i,k,j) = ften(i,k,j) - &
                  (vsv(i,k+1,j)*fg(i,k+1)-vsv(i,k,j)*fg(i,k))/dsigma(k)
        end do
      end do
!,.....k = kz
      do i = 2 , iym2
        ften(i,kz,j) = ften(i,kz,j) + vsv(i,kz,j)*fg(i,kz)/dsigma(kz)
      end do

    else if ( ind == 6 ) then

      do k = 2 , kz
        do i = 2 , iym1
          fg(i,k)= twt(k,1)*f(i,k,j) + twt(k,2)*f(i,k-1,j)
        end do
      end do

      do i = 2 , iym1
        if ( kpbl1d(i).gt.kz ) then
          call fatal(__FILE__,__LINE__,'kpbl1d is greater than KZ')
        end if
        if ( kpbl1d(i).ge.4 ) then
          ! Calculate slope of scalar in layer above ambiguous layer
          k = kpbl1d(i)-2
          if ( (f(i,k+1,j)-f(i,k,j)) > d_zero .and.   &
               (f(i,k,j)-f(i,k-1,j)) > d_zero ) then
            slope = min((f(i,k+1,j)-f(i,k,j))/(a(k+1)-a(k)),    &
                        (f(i,k,j)-f(i,k-1,j))/(a(k)-a(k-1)))
          else if ( (f(i,k+1,j)-f(i,k,j)) < d_zero .and.   &
                    (f(i,k,j)-f(i,k-1,j)) < d_zero ) then
            slope = max((f(i,k+1,j)-f(i,k,j))/(a(k+1)-a(k)),    &
                        (f(i,k,j)-f(i,k-1,j))/(a(k)-a(k-1)))
          else
            slope = d_zero
          end if
          ! Now replace the values of scalar at top and bottom of ambiguous
          ! layer as long as inversion is actually in the ambiguous layer
          k = kpbl1d(i)
          fg(i,k-1) = f(i,k-2,j) + slope*(sigma(k-1)-a(k-2))
          if (abs(f(i,k-2,j) + slope*(a(k-1)-a(k-2))-f(i,k,j)) > &
              abs(f(i,k-1,j)-f(i,k,j)) ) then
            fg(i,k) = f(i,k,j)
          else
            fg(i,k) = f(i,k-2,j) + slope*(sigma(k)-a(k-2))
          end if
        end if
      end do

      do i = 2 , iym1
        ften(i,1,j) = ften(i,1,j)-vsv(i,2,j)*fg(i,2)/dsigma(1)
      end do
      do k = 2 , kzm1
        do i = 2 , iym1
          ften(i,k,j) = ften(i,k,j)-(vsv(i,k+1,j)*fg(i,k+1)-vsv(i,k,j)*   &
                      fg(i,k))/dsigma(k)
        end do
      end do
      do i = 2 , iym1
        ften(i,kz,j) = ften(i,kz,j)+vsv(i,kz,j)*fg(i,kz)/dsigma(kz)
      end do

    end if
!
#endif
    call time_end(subroutine_name,idindx)
  end subroutine vadv
!
end module mod_advection
