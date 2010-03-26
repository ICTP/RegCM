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
!
      module mod_postproc_param
      implicit none
!
! PARAMETER definitions
!
      real(4) , parameter :: dtbc  = 6.00
      real(4) , parameter :: dtout = 6.00
      real(4) , parameter :: dtbat = 3.00
      real(4) , parameter :: dtrad = 6.00
      real(4) , parameter :: dtche = 6.00
      real(4) , parameter :: dtsub = 3.00

      integer , parameter :: nhrbc  = 24.001/dtbc
      integer , parameter :: nhrout = 24.001/dtout
      integer , parameter :: nhrbat = 24.001/dtbat
      integer , parameter :: nhrsub = 24.001/dtsub
      integer , parameter :: nhrrad = 24.001/dtrad
      integer , parameter :: nhrche = 24.001/dtche

      integer , parameter :: npl = 11

      real(4) , dimension(npl) :: plev
      real(4) , dimension(npl) :: plevr

      data    plev/1000.,925.,850.,700.,500.,400.,300.,250.,200.,       &
            &      150.,100./
      data    plevr/100.,150.,200.,250.,300.,400.,500.,700.,850.,       &
            &       925., 1000./

      integer , parameter :: nfmax = 999

      character(70) , parameter :: plist = 'postproc.in'
      character(70) , parameter :: ulist = 'user.in'

      end module mod_postproc_param
