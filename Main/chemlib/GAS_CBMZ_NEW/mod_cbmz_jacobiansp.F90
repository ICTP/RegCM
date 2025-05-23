! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! Sparse Jacobian Data Structures File
!
! Generated by KPP-2.2.3 symbolic chemistry Kinetics PreProcessor
!       (http://www.cs.vt.edu/~asandu/Software/KPP)
! KPP is distributed under GPL, the general public licence
!       (http://www.gnu.org/copyleft/gpl.html)
! (C) 1995-1997, V. Damian & A. Sandu, CGRER, Univ. Iowa
! (C) 1997-2005, A. Sandu, Michigan Tech, Virginia Tech
!     With important contributions from:
!        M. Damian, Villanova University, USA
!        R. Sander, Max-Planck Institute for Chemistry, Mainz, Germany
!
! File                 : mod_cbmz_JacobianSP.f90
! Time                 : Mon Nov 25 13:41:15 2013
! Working directory    : /scratch/ashalaby/kpp-2.2.3/compare/CBMZ
! Equation file        : mod_cbmz.kpp
! Output root filename : mod_cbmz
!
! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

module mod_cbmz_jacobiansp

  public
  save

  ! sparse jacobian data

  integer, parameter, dimension(360) :: lu_irow_0 = (/ &
       1,  1,  1,  1,  1,  1,  1,  2,  2,  2,  3,  3, &
       3,  3,  4,  4,  4,  4,  4,  4,  4,  4,  5,  5, &
       5,  6,  6,  6,  7,  7,  7,  8,  8,  9,  9,  9, &
      10, 10, 10, 10, 11, 11, 11, 11, 12, 12, 13, 13, &
      13, 13, 13, 14, 14, 14, 14, 14, 15, 15, 15, 15, &
      15, 16, 16, 16, 17, 17, 17, 17, 18, 18, 18, 18, &
      19, 19, 19, 19, 19, 20, 20, 20, 20, 20, 21, 21, &
      21, 21, 21, 22, 22, 22, 22, 23, 23, 23, 23, 24, &
      24, 24, 24, 25, 25, 25, 26, 26, 26, 26, 26, 27, &
      27, 27, 27, 27, 27, 27, 28, 28, 28, 28, 28, 28, &
      28, 28, 28, 28, 28, 28, 28, 28, 28, 29, 29, 29, &
      29, 29, 29, 29, 29, 29, 29, 29, 29, 29, 30, 30, &
      30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, &
      30, 30, 31, 31, 31, 31, 31, 31, 31, 32, 32, 32, &
      32, 32, 33, 33, 33, 33, 33, 34, 34, 34, 34, 34, &
      35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 36, 36, &
      36, 36, 36, 36, 36, 36, 36, 36, 36, 36, 37, 37, &
      37, 37, 38, 38, 38, 38, 39, 39, 39, 39, 39, 39, &
      39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, 39, &
      39, 40, 40, 40, 40, 41, 41, 41, 41, 41, 41, 41, &
      41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, 41, &
      42, 42, 42, 42, 42, 42, 42, 42, 42, 42, 42, 42, &
      42, 43, 43, 43, 43, 43, 43, 43, 43, 43, 43, 43, &
      43, 43, 44, 44, 44, 44, 44, 44, 44, 44, 44, 45, &
      45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, 45, &
      45, 45, 45, 45, 45, 45, 45, 46, 46, 46, 46, 46, &
      46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 46, 46, &
      47, 47, 47, 47, 47, 47, 47, 47, 47, 48, 48, 48, &
      48, 48, 48, 48, 48, 48, 48, 48, 48, 49, 49, 49, &
      49, 49, 49, 49, 49, 49, 49, 49, 49, 49, 50, 50 /)
  integer, parameter, dimension(240) :: lu_irow_1 = (/ &
      50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, 50, &
      50, 50, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, &
      51, 51, 51, 51, 51, 51, 52, 52, 52, 52, 52, 52, &
      52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, &
      52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, &
      52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 52, 53, &
      53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, 53, &
      53, 53, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, &
      54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, &
      54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, &
      54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, 54, &
      55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, &
      55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, 55, &
      55, 55, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, &
      56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, &
      56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, 56, &
      57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, &
      57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, 57, &
      57, 57, 57, 57, 58, 58, 58, 58, 58, 58, 58, 58, &
      58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58, 58 /)
  integer, parameter, dimension(600) :: lu_irow = (/&
    lu_irow_0, lu_irow_1 /)

  integer, parameter, dimension(360) :: lu_icol_0 = (/ &
       1, 25, 37, 40, 53, 55, 58,  2, 11, 54,  3, 25, &
      40, 53,  4, 37, 38, 40, 47, 52, 53, 58,  5, 16, &
      54,  6, 41, 52,  7, 56, 58,  8, 54,  9, 35, 53, &
      10, 35, 52, 54, 11, 16, 54, 57, 12, 54, 13, 37, &
      40, 53, 54, 14, 37, 40, 53, 54, 15, 26, 54, 56, &
      57, 16, 54, 57, 17, 52, 54, 56, 18, 40, 53, 54, &
       8, 12, 19, 54, 55, 20, 37, 40, 53, 54, 21, 52, &
      54, 55, 56, 22, 46, 52, 54, 23, 43, 52, 54, 24, &
      35, 56, 57, 25, 53, 54,  8, 12, 26, 54, 57,  9, &
      27, 35, 53, 55, 56, 57, 25, 28, 31, 34, 37, 38, &
      39, 40, 42, 45, 47, 53, 54, 55, 57, 16, 24, 26, &
      29, 35, 39, 42, 45, 47, 52, 54, 56, 57, 12, 30, &
      32, 33, 37, 40, 44, 47, 49, 50, 51, 52, 53, 54, &
      55, 57, 19, 26, 31, 53, 54, 55, 57, 32, 38, 52, &
      55, 57, 33, 38, 52, 54, 55, 34, 47, 52, 54, 55, &
       9, 18, 24, 35, 40, 52, 53, 54, 56, 57, 34, 36, &
      37, 47, 49, 50, 51, 52, 53, 54, 55, 57, 37, 53, &
      54, 57, 38, 53, 54, 57, 20, 22, 25, 31, 33, 34, &
      37, 38, 39, 40, 44, 46, 47, 48, 52, 53, 54, 55, &
      57, 40, 53, 54, 57,  8, 12, 25, 26, 31, 37, 38, &
      40, 41, 42, 47, 49, 50, 51, 52, 53, 54, 55, 57, &
      12, 31, 34, 37, 40, 42, 47, 49, 52, 53, 54, 55, &
      57, 14, 23, 37, 40, 43, 49, 50, 51, 52, 53, 54, &
      55, 57, 37, 40, 44, 51, 52, 53, 54, 55, 57, 23, &
      25, 31, 32, 34, 37, 38, 40, 43, 44, 45, 47, 49, &
      50, 51, 52, 53, 54, 55, 57, 13, 22, 36, 37, 40, &
      45, 46, 47, 49, 50, 51, 52, 53, 54, 55, 57, 58, &
      32, 33, 38, 47, 52, 53, 54, 55, 57, 36, 37, 47, &
      48, 49, 50, 51, 52, 53, 54, 55, 57, 33, 34, 38, &
      47, 48, 49, 50, 51, 52, 53, 54, 55, 57, 30, 32 /)
  integer, parameter, dimension(240) :: lu_icol_1 = (/ &
      33, 37, 38, 40, 44, 47, 49, 50, 51, 52, 53, 54, &
      55, 57, 15, 19, 26, 32, 33, 38, 44, 47, 50, 51, &
      52, 53, 54, 55, 56, 57,  8, 10, 11, 12, 16, 17, &
      18, 19, 20, 22, 23, 25, 26, 28, 31, 32, 33, 34, &
      35, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47, &
      48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, 25, &
      27, 31, 35, 37, 38, 40, 47, 52, 53, 54, 55, 56, &
      57, 58,  8,  9, 10, 11, 12, 13, 14, 16, 17, 18, &
      20, 21, 22, 23, 25, 26, 28, 29, 30, 31, 32, 33, &
      34, 35, 36, 37, 38, 39, 40, 42, 43, 44, 45, 46, &
      47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, &
      19, 21, 27, 32, 33, 34, 35, 38, 40, 41, 42, 43, &
      44, 46, 47, 48, 49, 50, 51, 52, 53, 54, 55, 56, &
      57, 58,  7, 15, 17, 19, 21, 24, 26, 27, 29, 32, &
      33, 34, 35, 38, 39, 40, 41, 42, 43, 44, 45, 46, &
      47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58, &
      16, 24, 26, 27, 29, 35, 37, 38, 39, 40, 41, 42, &
      43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54, &
      55, 56, 57, 58,  7, 31, 36, 37, 38, 40, 42, 45, &
      47, 48, 49, 50, 51, 52, 53, 54, 55, 56, 57, 58 /)
  integer, parameter, dimension(600) :: lu_icol = (/&
    lu_icol_0, lu_icol_1 /)

  integer, parameter, dimension(59) :: lu_crow = (/ &
       1,  8, 11, 15, 23, 26, 29, 32, 34, 37, 41, 45, &
      47, 52, 57, 62, 65, 69, 73, 78, 83, 88, 92, 96, &
     100,103,108,115,130,143,159,166,171,176,181,191, &
     203,207,211,230,234,253,266,279,288,308,325,334, &
     346,359,375,391,432,447,493,519,553,581,601 /)

  integer, parameter, dimension(59) :: lu_diag = (/ &
       1,  8, 11, 15, 23, 26, 29, 32, 34, 37, 41, 45, &
      47, 52, 57, 62, 65, 69, 75, 78, 83, 88, 92, 96, &
     100,105,109,116,133,144,161,166,171,176,184,192, &
     203,207,219,230,242,258,270,281,298,314,328,337, &
     351,368,384,425,441,488,515,550,579,600,601 /)

end module mod_cbmz_jacobiansp

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
