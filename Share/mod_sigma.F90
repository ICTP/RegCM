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

module mod_sigma

  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_memutil
  use mod_message
  use mod_stdio

  implicit none

  private

  real(rk8) , pointer , dimension(:) :: sigma_coordinate
  real(rk8) , pointer , dimension(:) :: sigma_delta
  real(rk8) , pointer , dimension(:) :: half_sigma_coordinate

  public :: sigma_coordinate
  public :: sigma_delta
  public :: half_sigma_coordinate
  public :: init_sigma

  contains
    !
    ! For the RegCM the sigma coordinate is terrain following
    !
    subroutine init_sigma(nk,dmax,dmin)
      implicit none
      integer(ik4) , intent(in) :: nk
      real(rk8) , intent(in) , optional :: dmax , dmin
      real(rk8) , allocatable , dimension(:) :: alph
      real(rk8) :: dsmax , dsmin
      real(rk8) :: jumpsize , apara , bpara , func , funcprev
      integer(ik4) :: k , icount , ierr
      integer(ik4) , parameter :: maxiter = 1000000
      real(rk8) , parameter :: conv_crit = 0.00001D0

      call getmem1d(sigma_coordinate,1,nk+1, &
        'init_sigma:sigma_coordinate')
      call getmem1d(sigma_delta,1,nk, &
        'init_sigma:sigma_delta')
      call getmem1d(half_sigma_coordinate,1,nk, &
        'init_sigma:half_sigma_coordinate')
      !
      ! Setup hardcoded sigma levels
      !
      if ( nk == 14 ) then      ! RegCM2
        sigma_coordinate(1) = 0.0D0
        sigma_coordinate(2) = 0.04D0
        sigma_coordinate(3) = 0.10D0
        sigma_coordinate(4) = 0.17D0
        sigma_coordinate(5) = 0.25D0
        sigma_coordinate(6) = 0.35D0
        sigma_coordinate(7) = 0.46D0
        sigma_coordinate(8) = 0.56D0
        sigma_coordinate(9) = 0.67D0
        sigma_coordinate(10) = 0.77D0
        sigma_coordinate(11) = 0.86D0
        sigma_coordinate(12) = 0.93D0
        sigma_coordinate(13) = 0.97D0
        sigma_coordinate(14) = 0.99D0
        sigma_coordinate(15) = 1.0D0
      else if ( nk == 18 ) then ! RegCM3, default
        sigma_coordinate(1) = 0.0D0
        sigma_coordinate(2) = 0.05D0
        sigma_coordinate(3) = 0.10D0
        sigma_coordinate(4) = 0.16D0
        sigma_coordinate(5) = 0.23D0
        sigma_coordinate(6) = 0.31D0
        sigma_coordinate(7) = 0.39D0
        sigma_coordinate(8) = 0.47D0
        sigma_coordinate(9) = 0.55D0
        sigma_coordinate(10) = 0.63D0
        sigma_coordinate(11) = 0.71D0
        sigma_coordinate(12) = 0.78D0
        sigma_coordinate(13) = 0.84D0
        sigma_coordinate(14) = 0.89D0
        sigma_coordinate(15) = 0.93D0
        sigma_coordinate(16) = 0.96D0
        sigma_coordinate(17) = 0.98D0
        sigma_coordinate(18) = 0.99D0
        sigma_coordinate(19) = 1.0D0
      else if ( nk == 23 ) then ! MM5V3
        sigma_coordinate(1) = 0.0D0
        sigma_coordinate(2) = 0.05D0
        sigma_coordinate(3) = 0.1D0
        sigma_coordinate(4) = 0.15D0
        sigma_coordinate(5) = 0.2D0
        sigma_coordinate(6) = 0.25D0
        sigma_coordinate(7) = 0.3D0
        sigma_coordinate(8) = 0.35D0
        sigma_coordinate(9) = 0.4D0
        sigma_coordinate(10) = 0.45D0
        sigma_coordinate(11) = 0.5D0
        sigma_coordinate(12) = 0.55D0
        sigma_coordinate(13) = 0.6D0
        sigma_coordinate(14) = 0.65D0
        sigma_coordinate(15) = 0.7D0
        sigma_coordinate(16) = 0.75D0
        sigma_coordinate(17) = 0.8D0
        sigma_coordinate(18) = 0.85D0
        sigma_coordinate(19) = 0.89D0
        sigma_coordinate(20) = 0.93D0
        sigma_coordinate(21) = 0.96D0
        sigma_coordinate(22) = 0.98D0
        sigma_coordinate(23) = 0.99D0
        sigma_coordinate(24) = 1.0D0
      else ! Compute (or try to do so)
        dsmax = 0.05D0
        dsmin = 0.0025D0
        if ( present(dmax) ) dsmax = dmax
        if ( present(dmin) ) dsmin = dmin
        allocate(alph(nk), stat=ierr)
        call checkalloc(ierr,__FILE__,__LINE__,'init_sigma:alph')
        write (stdout,*) 'Creating a custom set of sigma levels : '
        if ( (dsmax*dble(nk)) < d_one ) then
          write (stderr,*) 'dsmax must be greater than ', d_one/dble(nk)
          write (stderr,*) 'or kz must be less than ', d_one/dsmax
          call die('init_sigma','Maximum resolution, dsmax, is too low.',1)
        end if
        if ( (dsmin*dble(nk)) >= d_one ) then
          write (stderr,*) 'dsmin must be less than ', d_one/dble(nk)
          write (stderr,*) 'or kz must be greater than ', d_one/dsmax
          call die('init_sigma','Minimum resolution, dsmin, is too large.',1)
        end if
        ! Do a function minimization to determine the a,b coefficients for the
        ! following equation:
        !
        !    sigma_delta(i) = dsmax*a^(i-1)*b^(0.5*(i-2)*(i-1))
        !
        ! which is derived from the recursive relation:
        !
        !    sigma_delta(i) = a(i)*sigma_delta(i-1) with a(i) = b*a(i-1)
        !
        ! The function provides level spacings between dsmin and dsmax such
        ! that level thicknesses vary approximately exponentially up to the
        ! TOA and gives more levels toward the surface.
        !
        ! Set the initial conditions for the function minimization
        !
        jumpsize = 0.0015D0
        bpara = 0.99573D0
        apara = ((dsmin/dsmax)**(d_one/dble(nk-1))) * &
                 (bpara**(-d_half*dble(nk-2)))
        alph(1) = apara/bpara
        sigma_delta(1) = dsmax
        do k = 2 , nk
          alph(k) = bpara*alph(k-1)
          sigma_delta(k) = alph(k)*sigma_delta(k-1)
        end do
        func = sum(sigma_delta)-d_one
        ! Loop through the minimization until the convergence
        ! criterion is satisfied
        do icount = 1 , maxiter
          funcprev = func
          bpara = bpara + jumpsize
          if ( bpara < d_zero ) bpara = 1D-10
          apara = ((dsmin/dsmax)**(d_one/dble(nk-1))) * &
                   (bpara**(-d_half*dble(nk-2)))
          alph(1) = apara/bpara
          sigma_delta(1) = dsmax
          do k = 2 , nk
            alph(k) = bpara*alph(k-1)
            sigma_delta(k) = alph(k)*sigma_delta(k-1)
          end do
          func = sum(sigma_delta)-d_one
          ! If we overshot 0, then reduce the jump size and reverse course
          if ( func*funcprev < d_zero ) then
            jumpsize = -jumpsize/d_two
          else if ( abs(func) > abs(funcprev) ) then
            jumpsize = -jumpsize
          end if
          ! If we converged, then jump out of the loop
          if ( abs(func) < conv_crit ) then
            write (stdout,*) 'Convergence reached.'
            write (stdout,*) '#', apara, bpara, icount
            write (stdout,*) '#', dsmax, dsmin, sum(sigma_delta)
            exit
          end if
          if ( icount == maxiter-1 ) then
            write (stderr,*) 'Failed to converge.'
            write (stderr,*) 'b,a,jumpsize,func,funcprev:', bpara, apara, &
                           jumpsize,func,funcprev
            call die('init_sigma','Error setting up custom sigma levels')
          end if
        end do
        sigma_coordinate(1) = 0.0D0
        do k = 1 , nk-1
          sigma_coordinate(k+1) = sigma_coordinate(k)+sigma_delta(k)
        end do
        sigma_coordinate(nk+1) = 1.0D0
        deallocate(alph)
      end if
      do k = 1 , nk
        sigma_delta(k) = sigma_coordinate(k+1) - sigma_coordinate(k)
        half_sigma_coordinate(k) = (sigma_coordinate(k) + &
                                    sigma_coordinate(k+1)) * d_half
      end do
    end subroutine init_sigma

end module mod_sigma

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
