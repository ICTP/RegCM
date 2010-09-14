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

      module mod_runparams

      use mod_dynparam

      implicit none
 
      integer :: ifrabe , klake , nbatst
!
      real(8) :: dt , dt0 , dt2 , dtbat , dtlake , dtmin
      real(8) :: dx , dx2 , dx4 , dx8 , dx16 , dxsq
      real(8) :: abatm , abemh
      real(8) :: c200 , c203 
      real(8) :: fnudge , gnudge
      real(8) :: xkhmax , xkhz

      integer :: ibltyp , iboudy , ichem , icnt , icup , idirect ,      &
               & iemiss , igcc , iocnflx , ipgf , ipptls , kbats ,      &
               & kchem , lakemod , nradisp , ntrad , ntsave , nttape ,  &
               & idcsst , iseaice

      logical :: ifrest , rfstrt
 
      real(8) :: bdytim , prttim , radfrq , savtim , taptim , tbdybe
      integer :: ndbgfrq , nsavfrq , ntapfrq

      integer :: ispgd , ispgx , k700 , kchi , kclo , kcmd , kt , ncld
!
      real(8) :: r8pt
      real(8) :: akht1 , akht2
      real(8) :: high_nudge , medium_nudge , low_nudge

      real(8) , allocatable , dimension(:) :: dtau
      real(8) , allocatable , dimension(:) :: a , anudg , dsigma , qcon
      real(8) , allocatable , dimension(:) :: sigma
      real(8) , allocatable , dimension(:,:) :: twt
      real(8) , allocatable , dimension(:) :: wgtd
      real(8) , allocatable , dimension(:) :: wgtx

      character(len=3) :: scenario

      contains

      subroutine allocate_mod_runparams
        implicit none
        allocate(a(kz))
        allocate(anudg(kz))
        allocate(dsigma(kz))
        allocate(qcon(kz))
        allocate(sigma(kzp1))
        allocate(twt(kz,2))
        allocate(wgtd(nspgd))
        allocate(wgtx(nspgx))
        allocate(dtau(nsplit))
      end subroutine allocate_mod_runparams
!
      end module mod_runparams
