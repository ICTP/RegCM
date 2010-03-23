!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of RegCM model.
!
!    RegCM model is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    RegCM model is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with RegCM model.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      subroutine mtrxclm
!
!=======================================================================
!  based on: clm version 3.5
!=======================================================================
! Matrix CLM was written by Jason Bell for coupling the clm3.0 land
! surface model to RegCM3.
!! 2008.5.8    A Tawfik Revised to work with RegCM3 and clm3.5
!
! mtrxclm must be called outside of the jslc loop. CLM is passed full
! arrays of variables and returns likewise. This is done to avoid
! altering the functionality of CLM as much as possible.
!=======================================================================
!=======================================================================
!  si version  - water fluxes are generally calculated in kg/m**2/s.
!  1000 kg/m**2/s = 1 m/s  - converted to energy units for display:
!                            1 kg/m**2/s = 2.5 x 10.e6  watts/m**2.
!  note also  1 kg/m**2/s = 1 mm/m**2/s so fluxes can be thought of
!                            in mm/s.
!=======================================================================
! include clm modules
      use atmdrvMod       , only : rcmdrv
      use clm_comp        , only : clm_run1, clm_run2
      use shr_kind_mod    , only : r8 => shr_kind_r8
!
      use mod_message
      use mod_clm
      implicit none

!---------------------------------------------------------------------
 
#ifdef MPP1
      call interfclm(1)
#else
      write (*,*) 'CLM not implemented for serial run'
      call fatal(__FILE__,__LINE__, 'CLM not implemented')
!      call interfclm_ser(1)
#endif

      call rcmdrv()
      call clm_run1(r2cdoalb,r2ceccen,r2cobliqr,r2clambm0,r2cmvelpp)
      call clm_run2(r2ceccen,r2cobliqr,r2clambm0,r2cmvelpp)
 
#ifdef MPP1
      call interfclm(2)
#else
      write (*,*) 'CLM not implemented for serial run'
      call fatal(__FILE__,__LINE__, 'CLM not implemented')
!      call interfclm_ser(2)
#endif

      end subroutine mtrxclm
