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

      use mod_dynparam

      implicit none
!
      integer , parameter :: maxntr = 20
      integer , parameter :: maxnbin = 20

      character(5) , allocatable , dimension(:) :: chtrname
!
      real(8) , allocatable , dimension(:,:) :: chtrdpv
      real(8) , allocatable , dimension(:,:) :: chtrsize , dustbsiz
      real(8) , allocatable , dimension(:) :: chtrsol
!
      integer :: ichcumtra , ichdrdepo , ichremcvc , ichremlsc ,        &
               & ichsursrc
      integer , allocatable, dimension(:,:) :: icumbot , icumdwd ,      &
               & icumtop
!
      integer :: ibchb , ibchl , iochb , iochl , iso2 , iso4 , mixtype
      integer , allocatable , dimension(:) :: idust
!
      real(8) , allocatable , dimension(:,:) :: mflx
!
      real(8) , allocatable , dimension(:,:,:) :: aerasp , aerext ,     &
                                & aerssa
      real(8) , allocatable , dimension(:,:) :: aersrrf , aertarf
!
      real(8) , allocatable , dimension(:,:,:) :: cemtr , cemtrac ,     &
                        & remdrd
      real(8) , allocatable , dimension(:,:) :: rembc , remrat
      real(8) , allocatable , dimension(:,:,:,:) :: remcvc , remlsc ,   &
                        & rxsaq1 , rxsaq2 , rxsg

contains
        subroutine allocate_mod_trachem
        use mod_message , only : say , aline
        implicit none
        if ( ntr>maxntr ) then
          write (aline , *) 'In mod_trachem, resetting ntr to maxntr ', &
                 maxntr
          call say
          ntr = maxntr
        end if
        if ( nbin>maxnbin ) then
          write (aline , *) 'In mod_trachem, resetting nbin to maxbin ',&
                 maxnbin
          call say
          nbin = maxnbin
        end if
#ifdef MPP1
        allocate(icumbot(iy,jxp))
        allocate(icumdwd(iy,jxp))
        allocate(icumtop(iy,jxp)) 
        allocate(aerasp(iym1,kz,jxp))
        allocate(aerext(iym1,kz,jxp))
        allocate(aerssa(iym1,kz,jxp))
        allocate(aersrrf(iym1,jxp))
        allocate(aertarf(iym1,jxp))
        allocate(cemtr(iy,jxp,ntr))
        allocate(cemtrac(iy,jxp,ntr))
        allocate(remdrd(iy,jxp,ntr))
        allocate(remcvc(iy,kz,jxp,ntr))
        allocate(remlsc(iy,kz,jxp,ntr))
        allocate(rxsaq1(iy,kz,jxp,ntr))
        allocate(rxsaq2(iy,kz,jxp,ntr))
        allocate(rxsg(iy,kz,jxp,ntr))
#else
        allocate(icumbot((iy,jx))
        allocate(icumdwd((iy,jx))
        allocate(icumtop((iy,jx))
        allocate(aerasp(iym1,kz,jxm1))
        allocate(aerext(iym1,kz,jxm1))
        allocate(aerssa(iym1,kz,jxm1))
        allocate(aersrrf(iym1,jxm1))
        allocate(aertarf(iym1,jxm1))
        allocate(cemtr(iy,jx,ntr))
        allocate(cemtrac(iy,jx,ntr))
        allocate(remdrd(iy,jx,ntr))
        allocate(remcvc(iy,kz,jx,ntr))
        allocate(remlsc(iy,kz,jx,ntr))
        allocate(rxsaq1(iy,kz,jx,ntr))
        allocate(rxsaq2(iy,kz,jx,ntr))
        allocate(rxsg(iy,kz,jx,ntr))
#endif 
        allocate(chtrname(ntr))
        allocate(chtrdpv(ntr,2))
        allocate(chtrsize(nbin,2))
        allocate(dustbsiz(nbin,2))
        allocate(chtrsol(ntr))
        allocate(idust(nbin))
        allocate(mflx(iy,2))
        allocate(rembc(iy,kz))
        allocate(remrat(iy,kz))

        end subroutine allocate_mod_trachem

      end module mod_trachem
