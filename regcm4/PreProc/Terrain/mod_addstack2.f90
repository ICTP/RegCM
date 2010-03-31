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

      module mod_addstack2
      use mod_regcm_param , only : iy , jx , nsg , iysg , jxsg
      implicit none
      real(4) :: clong
      real(4) , dimension(iy,jx) :: corc , hscr1 , htsavc , sumc ,      &
                                  & wtmaxc
      real(4) , dimension(iysg,jxsg) :: corc_s , hscr1_s , htsavc_s ,   &
                                  & sumc_s , wtmaxc_s
      real(4) , dimension(iy,jx,2) :: itex , land
      real(4) , dimension(iysg,jxsg,2) :: itex_s , land_s
      integer , dimension(iy,jx) :: nsc
      integer , dimension(iysg,jxsg) :: nsc_s
      end module mod_addstack2
