      module parkind

      use mod_realkinds , only : rkx

      implicit none
      save

!------------------------------------------------------------------
! rrtmg kinds
! Define integer and real kinds for various types.
!
! Initial version: MJIacono, AER, jun2006
! Revised: MJIacono, AER, aug2008
!------------------------------------------------------------------

!
!     integer kinds
!     -------------
!
      integer, parameter :: kind_ib = selected_int_kind(13)  ! 8 byte integer
      integer, parameter :: kind_im = selected_int_kind(6)   ! 4 byte integer
      integer, parameter :: kind_in = kind(1)                ! native integer

!
!     real kinds
!     ----------
!
      !integer, parameter :: kind_rb = selected_real_kind(12) ! 8 byte real
      integer, parameter :: kind_rb = rkx
      integer, parameter :: kind_rm = selected_real_kind(6)  ! 4 byte real
      integer, parameter :: kind_rn = kind(1.0)              ! native real
      !
      ! Epsilon for numerical consistency
      !
#ifdef SINGLE_PRECISION_REAL
      real(kind_rb) , parameter :: almostzero = 1.e-10_kind_rb
#else
      real(kind_rb) , parameter :: almostzero = 1.e-20_kind_rb
#endif

      end module parkind
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
