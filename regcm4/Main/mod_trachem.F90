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

      module mod_trachem

      use mod_regcm_param

      implicit none
!
      character(5) , dimension(ntr) :: chtrname
!
      real(8) , dimension(ntr,2) :: chtrdpv
      real(8) , dimension(nbin,2) :: chtrsize , dustbsiz
      real(8) , dimension(ntr) :: chtrsol
!
      integer :: ichcumtra , ichdrdepo , ichremcvc , ichremlsc ,        &
               & ichsursrc
#ifdef MPP1
      integer , dimension(iy,jxp) :: icumbot , icumdwd , icumtop
#else
      integer , dimension(iy,jx) :: icumbot , icumdwd , icumtop
#endif
!
      integer :: ibchb , ibchl , iochb , iochl , iso2 , iso4 , mixtype
      integer , dimension(nbin) :: idust
!
      real(8) , dimension(iy,2) :: mflx
!
#ifdef MPP1
      real(8) , dimension(iym1,kz,jxp) :: aerasp , aerext , aerssa
      real(8) , dimension(iym1,jxp) :: aersrrf , aertarf
#else
      real(8) , dimension(iym1,kz,jxm1) :: aerasp , aerext , aerssa
      real(8) , dimension(iym1,jxm1) :: aersrrf , aertarf
#endif
!
#ifdef MPP1
      real(8) , dimension(iy,jxp,ntr) :: cemtr , cemtrac , remdrd
      real(8) , dimension(iy,kz) :: rembc , remrat
      real(8) , dimension(iy,kz,jxp,ntr) :: remcvc , remlsc , rxsaq1 ,  &
           & rxsaq2 , rxsg
#else
      real(8) , dimension(iy,jx,ntr) :: cemtr , cemtrac , remdrd
      real(8) , dimension(iy,kz) :: rembc , remrat
      real(8) , dimension(iy,kz,jx,ntr) :: remcvc , remlsc , rxsaq1 ,   &
           & rxsaq2 , rxsg
#endif

      end module mod_trachem
