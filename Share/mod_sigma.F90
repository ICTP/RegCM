!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    Use of this source code is governed by an MIT-style license that can
!    be found in the LICENSE file or at
!
!         https://opensource.org/licenses/MIT.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
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

  real(rkx), pointer, contiguous, dimension(:) :: sigma_coordinate
  real(rkx), pointer, contiguous, dimension(:) :: half_sigma_coordinate
  real(rkx), pointer, contiguous, dimension(:) :: sigma_delta

  real(rkx) :: ptop
  logical :: is_pstar = .true.
  real(rkx), pointer, contiguous, dimension(:,:) :: ps => null()

  real(rkx), pointer, contiguous, dimension(:,:,:) :: pprime => null()

  public :: sigma_coordinate
  public :: sigma_delta
  public :: half_sigma_coordinate
  public :: init_sigma

  public :: init_hydrostatic
  public :: hydrostatic_pressure_full_sigma
  public :: hydrostatic_pressure_half_sigma
  public :: hydrostatic_deltap_full_sigma

  public :: init_non_hydrostatic
  public :: non_hydrostatic_pressure_full_sigma
  public :: non_hydrostatic_pressure_half_sigma
  public :: non_hydrostatic_deltap_full_sigma

  !
  ! Notation:
  !            surface_pressure is the measured surface pressure
  !            pstar = surface_pressure - ptop
  !
  ! Hydrostatic:
  !
  !     pressure = sigma * pstar + ptop =
  !                surface_pressure * sigma + ( 1.0 - sigma ) * ptop
  !
  contains
    !
    ! For the RegCM the sigma coordinate is terrain following
    !
    subroutine init_sigma(nk,dmax,dmin,lzeta)
      implicit none
      integer(ik4), intent(in) :: nk
      real(rkx), intent(in), optional :: dmax, dmin
      logical, intent(in), optional :: lzeta
      real(rkx), allocatable, dimension(:) :: alph
      real(rkx) :: dsmax, dsmin
      real(rkx) :: jumpsize, apara, bpara, func, funcprev
      integer(ik4) :: k, icount, ierr
      integer(ik4), parameter :: maxiter = 1000000
      real(rkx), parameter :: conv_crit = 0.00001_rkx
      logical :: lpress, lcompute

      lpress = .true.
      if ( present(lzeta) ) then
        lpress = .not. lzeta
      end if
      lcompute = .not. lpress .or. ( lpress .and. &
              .not. (nk == 2 .or. nk == 14 .or. nk == 18 .or. &
                     nk == 23 .or. nk == 41) )
      call getmem1d(sigma_coordinate,1,nk+1, &
        'init_sigma:sigma_coordinate')
      call getmem1d(sigma_delta,1,nk, &
        'init_sigma:sigma_delta')
      call getmem1d(half_sigma_coordinate,1,nk, &
        'init_sigma:half_sigma_coordinate')
      !
      ! Setup hardcoded sigma levels
      !
      if ( .not. lcompute ) then
        write (stdout,*) 'Using pre-defined set of sigma levels.'
        if ( nk == 2 ) then            ! For lower only for clmsa
          sigma_coordinate(1) = 0.0_rkx
          sigma_coordinate(2) = 0.99_rkx
          sigma_coordinate(3) = 1.0_rkx
        else if ( nk == 14 ) then      ! RegCM2
          sigma_coordinate(1) = 0.0_rkx
          sigma_coordinate(2) = 0.04_rkx
          sigma_coordinate(3) = 0.10_rkx
          sigma_coordinate(4) = 0.17_rkx
          sigma_coordinate(5) = 0.25_rkx
          sigma_coordinate(6) = 0.35_rkx
          sigma_coordinate(7) = 0.46_rkx
          sigma_coordinate(8) = 0.56_rkx
          sigma_coordinate(9) = 0.67_rkx
          sigma_coordinate(10) = 0.77_rkx
          sigma_coordinate(11) = 0.86_rkx
          sigma_coordinate(12) = 0.93_rkx
          sigma_coordinate(13) = 0.97_rkx
          sigma_coordinate(14) = 0.99_rkx
          sigma_coordinate(15) = 1.0_rkx
        else if ( nk == 18 ) then ! RegCM3, default
          sigma_coordinate(1) = 0.0_rkx
          sigma_coordinate(2) = 0.05_rkx
          sigma_coordinate(3) = 0.10_rkx
          sigma_coordinate(4) = 0.16_rkx
          sigma_coordinate(5) = 0.23_rkx
          sigma_coordinate(6) = 0.31_rkx
          sigma_coordinate(7) = 0.39_rkx
          sigma_coordinate(8) = 0.47_rkx
          sigma_coordinate(9) = 0.55_rkx
          sigma_coordinate(10) = 0.63_rkx
          sigma_coordinate(11) = 0.71_rkx
          sigma_coordinate(12) = 0.78_rkx
          sigma_coordinate(13) = 0.84_rkx
          sigma_coordinate(14) = 0.89_rkx
          sigma_coordinate(15) = 0.93_rkx
          sigma_coordinate(16) = 0.96_rkx
          sigma_coordinate(17) = 0.98_rkx
          sigma_coordinate(18) = 0.99_rkx
          sigma_coordinate(19) = 1.0_rkx
        else if ( nk == 23 ) then ! MM5V3
          sigma_coordinate(1) = 0.0_rkx
          sigma_coordinate(2) = 0.05_rkx
          sigma_coordinate(3) = 0.1_rkx
          sigma_coordinate(4) = 0.15_rkx
          sigma_coordinate(5) = 0.2_rkx
          sigma_coordinate(6) = 0.25_rkx
          sigma_coordinate(7) = 0.3_rkx
          sigma_coordinate(8) = 0.35_rkx
          sigma_coordinate(9) = 0.4_rkx
          sigma_coordinate(10) = 0.45_rkx
          sigma_coordinate(11) = 0.5_rkx
          sigma_coordinate(12) = 0.55_rkx
          sigma_coordinate(13) = 0.6_rkx
          sigma_coordinate(14) = 0.65_rkx
          sigma_coordinate(15) = 0.7_rkx
          sigma_coordinate(16) = 0.75_rkx
          sigma_coordinate(17) = 0.8_rkx
          sigma_coordinate(18) = 0.85_rkx
          sigma_coordinate(19) = 0.89_rkx
          sigma_coordinate(20) = 0.93_rkx
          sigma_coordinate(21) = 0.96_rkx
          sigma_coordinate(22) = 0.98_rkx
          sigma_coordinate(23) = 0.99_rkx
          sigma_coordinate(24) = 1.0_rkx
        else if ( nk == 41 ) then ! RegCM VHR
          sigma_coordinate(1) = 0.0000_rkx
          sigma_coordinate(2) = 0.0500_rkx
          sigma_coordinate(3) = 0.0978_rkx
          sigma_coordinate(4) = 0.1436_rkx
          sigma_coordinate(5) = 0.1875_rkx
          sigma_coordinate(6) = 0.2295_rkx
          sigma_coordinate(7) = 0.2697_rkx
          sigma_coordinate(8) = 0.3082_rkx
          sigma_coordinate(9) = 0.3451_rkx
          sigma_coordinate(10) = 0.3804_rkx
          sigma_coordinate(11) = 0.4143_rkx
          sigma_coordinate(12) = 0.4468_rkx
          sigma_coordinate(13) = 0.4779_rkx
          sigma_coordinate(14) = 0.5078_rkx
          sigma_coordinate(15) = 0.5364_rkx
          sigma_coordinate(16) = 0.5639_rkx
          sigma_coordinate(17) = 0.5903_rkx
          sigma_coordinate(18) = 0.6156_rkx
          sigma_coordinate(19) = 0.6399_rkx
          sigma_coordinate(20) = 0.6632_rkx
          sigma_coordinate(21) = 0.6856_rkx
          sigma_coordinate(22) = 0.7071_rkx
          sigma_coordinate(23) = 0.7277_rkx
          sigma_coordinate(24) = 0.7476_rkx
          sigma_coordinate(25) = 0.7667_rkx
          sigma_coordinate(26) = 0.7850_rkx
          sigma_coordinate(27) = 0.8027_rkx
          sigma_coordinate(28) = 0.8196_rkx
          sigma_coordinate(29) = 0.8359_rkx
          sigma_coordinate(30) = 0.8516_rkx
          sigma_coordinate(31) = 0.8667_rkx
          sigma_coordinate(32) = 0.8812_rkx
          sigma_coordinate(33) = 0.8952_rkx
          sigma_coordinate(34) = 0.9087_rkx
          sigma_coordinate(35) = 0.9216_rkx
          sigma_coordinate(36) = 0.9341_rkx
          sigma_coordinate(37) = 0.9461_rkx
          sigma_coordinate(38) = 0.9577_rkx
          sigma_coordinate(39) = 0.9689_rkx
          sigma_coordinate(40) = 0.9796_rkx
          sigma_coordinate(41) = 0.9900_rkx
          sigma_coordinate(42) = 1.0000_rkx
        end if
      else ! Compute (or try to do so)
        if ( lpress ) then
          dsmax = 0.05_rkx
          dsmin = 0.01_rkx
        else
          dsmax = 0.005_rkx
          dsmin = 0.01_rkx
        end if
        if ( present(dmax) ) dsmax = dmax
        if ( present(dmin) ) dsmin = dmin
        allocate(alph(nk), stat=ierr)
        call checkalloc(ierr,__FILE__,__LINE__,'init_sigma:alph')
        write (stdout,*) 'Creating a set of sigma levels... '
        if ( (dsmax*real(nk,rkx)) < d_one ) then
          write (stderr,*) 'dsmax must be greater than ', d_one/real(nk,rkx)
          write (stderr,*) 'or kz must be less than ', d_one/dsmax
          call die('init_sigma','Maximum resolution, dsmax, is too low.',1)
        end if
        if ( (dsmin*real(nk,rkx)) >= d_one ) then
          write (stderr,*) 'dsmin must be less than ', d_one/real(nk,rkx)
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
        jumpsize = 0.0015_rkx
        bpara = 0.99573_rkx
        apara = ((dsmin/dsmax)**(d_one/real(nk-1,rkx))) * &
                 (bpara**(-d_half*real(nk-2,rkx)))
        alph(1) = apara/bpara
        sigma_delta(1) = dsmax
        do k = 2, nk
          alph(k) = bpara*alph(k-1)
          sigma_delta(k) = alph(k)*sigma_delta(k-1)
        end do
        func = sum(sigma_delta)-d_one
        ! Loop through the minimization until the convergence
        ! criterion is satisfied
        do icount = 1, maxiter
          funcprev = func
          bpara = bpara + jumpsize
          if ( bpara < d_zero ) bpara = 1e-10_rkx
          apara = ((dsmin/dsmax)**(d_one/real(nk-1,rkx))) * &
                   (bpara**(-d_half*real(nk-2,rkx)))
          alph(1) = apara/bpara
          sigma_delta(1) = dsmax
          do k = 2, nk
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
        sigma_coordinate(1) = 0.0_rkx
        do k = 1, nk-1
          sigma_coordinate(k+1) = sigma_coordinate(k)+sigma_delta(k)
        end do
        sigma_coordinate(nk+1) = 1.0_rkx
        deallocate(alph)
      end if
      do k = 1, nk
        sigma_delta(k) = sigma_coordinate(k+1) - sigma_coordinate(k)
        half_sigma_coordinate(k) = (sigma_coordinate(k) + &
                                    sigma_coordinate(k+1)) * d_half
      end do
    end subroutine init_sigma

    subroutine init_hydrostatic(ptin,psin,lpstar)
      implicit none
      real(rkx), intent(in) :: ptin
      real(rkx), pointer, contiguous, dimension(:,:) :: psin
      logical, optional, intent(in) :: lpstar
      if ( present(lpstar) ) is_pstar = lpstar
      ptop = ptin
      ps => psin
    end subroutine init_hydrostatic

    real(rkx) elemental function pstar(surface_pressure)
      implicit none
      real(rkx), intent(in) :: surface_pressure
      pstar = surface_pressure - ptop
    end function pstar

    pure real(rkx) function hydrostatic_pressure_full_sigma(j,i,k) result(p)
      implicit none
      integer(ik4), intent(in) :: j, i, k
      if ( is_pstar ) then
        p = ps(j,i) * sigma_coordinate(k) + ptop
      else
        p = pstar(ps(j,i)) * sigma_coordinate(k) + ptop
      end if
    end function hydrostatic_pressure_full_sigma

    pure real(rkx) function hydrostatic_pressure_half_sigma(j,i,k) result(p)
      implicit none
      integer(ik4), intent(in) :: j, i, k
      if ( is_pstar ) then
        p = ps(j,i) * half_sigma_coordinate(k) + ptop
      else
        p = pstar(ps(j,i)) * half_sigma_coordinate(k) + ptop
      end if
    end function hydrostatic_pressure_half_sigma

    pure real(rkx) function hydrostatic_deltap_full_sigma(j,i,k) result(p)
      implicit none
      integer(ik4), intent(in) :: j, i, k
      if ( is_pstar ) then
        p = ps(j,i) * sigma_delta(k)
      else
        p = pstar(ps(j,i)) * sigma_delta(k)
      end if
    end function hydrostatic_deltap_full_sigma

    subroutine init_non_hydrostatic(ptin,psin,ppin,lpstar)
      implicit none
      real(rkx), intent(in) :: ptin
      real(rkx), pointer, contiguous, dimension(:,:) :: psin
      real(rkx), pointer, contiguous, dimension(:,:,:) :: ppin
      logical, optional, intent(in) :: lpstar
      if ( present(lpstar) ) is_pstar = lpstar
      ptop = ptin
      ps => psin
      pprime => ppin
    end subroutine init_non_hydrostatic

    pure real(rkx) function non_hydrostatic_pressure_full_sigma(j,i,k) result(p)
      implicit none
      integer(ik4), intent(in) :: j, i, k
      if ( k == 1 ) then
         p = ptop
      else
        if ( is_pstar ) then
          p = ps(j,i) * sigma_coordinate(k) + ptop + &
                    d_half * (pprime(j,i,k-1) - pprime(j,i,k))
        else
          p = pstar(ps(j,i)) * sigma_coordinate(k) + ptop + &
                    d_half * (pprime(j,i,k-1) - pprime(j,i,k))
        end if
      end if
    end function non_hydrostatic_pressure_full_sigma

    pure real(rkx) function non_hydrostatic_pressure_half_sigma(j,i,k) result(p)
      implicit none
      integer(ik4), intent(in) :: j, i, k
      if ( is_pstar ) then
        p = ps(j,i) * half_sigma_coordinate(k) + ptop + pprime(j,i,k)
      else
        p = pstar(ps(j,i)) * half_sigma_coordinate(k) + ptop + pprime(j,i,k)
      end if
    end function non_hydrostatic_pressure_half_sigma

    pure real(rkx) function non_hydrostatic_deltap_full_sigma(j,i,k) result(p)
      implicit none
      integer(ik4), intent(in) :: j, i, k
      if ( is_pstar ) then
        p = ps(j,i) * sigma_delta(k) + pprime(j,i,k)
      else
        p = pstar(ps(j,i)) * sigma_delta(k) + pprime(j,i,k)
      end if
    end function non_hydrostatic_deltap_full_sigma

end module mod_sigma

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
