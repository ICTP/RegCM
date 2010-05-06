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
      integer , allocatable, dimension(:,:) :: icumbot , icumdwd , icumtop
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
      real(8) , allocatable, dimension(:,:,:) :: aerasp , aerext , aerssa
      real(8) , allocatable, dimension(:,:) :: aersrrf , aertarf
#else
      real(8) , dimension(iym1,kz,jxm1) :: aerasp , aerext , aerssa
      real(8) , dimension(iym1,jxm1) :: aersrrf , aertarf
#endif
!
#ifdef MPP1
      real(8) , allocatable, dimension(:,:,:) :: cemtr , cemtrac , remdrd
      real(8) , dimension(iy,kz) :: rembc , remrat
      real(8) , allocatable, dimension(:,:,:,:) :: remcvc , remlsc , rxsaq1 ,  &
           & rxsaq2 , rxsg
#else
      real(8) , dimension(iy,jx,ntr) :: cemtr , cemtrac , remdrd
      real(8) , dimension(iy,kz) :: rembc , remrat
      real(8) , dimension(iy,kz,jx,ntr) :: remcvc , remlsc , rxsaq1 ,   &
           & rxsaq2 , rxsg
#endif

contains
	subroutine allocate_mod_trachem

#ifdef MPP1
	allocate(icumbot(iy,jxp) )
	allocate(icumdwd(iy,jxp) )
	allocate(icumtop(iy,jxp) ) 

        allocate(aerasp(iym1,kz,jxp) )
        allocate(aerext(iym1,kz,jxp) )
        allocate(aerssa(iym1,kz,jxp) )

        allocate(aersrrf(iym1,jxp) )
        allocate(aertarf(iym1,jxp) )

        allocate(cemtr(iy,jxp,ntr))
        allocate(cemtrac(iy,jxp,ntr))
        allocate(remdrd(iy,jxp,ntr))

        allocate(remcvc(iy,kz,jxp,ntr))
        allocate(remlsc(iy,kz,jxp,ntr))
        allocate(rxsaq1(iy,kz,jxp,ntr))
        allocate(rxsaq2(iy,kz,jxp,ntr))
        allocate(rxsg(iy,kz,jxp,ntr))

#else

#endif 

        end subroutine allocate_mod_trachem

      end module mod_trachem
