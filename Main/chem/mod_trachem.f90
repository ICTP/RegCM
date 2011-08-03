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

  use m_realkinds
  use mod_dynparam
  use mod_memutil
  use mod_message
  use mod_chem
!
  private
!
  integer , public , parameter :: maxntr = 20
  integer , public , parameter :: maxnbin = 20
!
  character(5) , public , allocatable , dimension(:) :: chtrname
!
  real(dp) , public , pointer , dimension(:,:) :: chtrdpv
  real(dp) , public , pointer , dimension(:,:) :: chtrsize , dustbsiz
  real(dp) , public , pointer , dimension(:) :: chtrsol
!
  integer , public :: ichcumtra , ichdrdepo , ichremcvc , ichremlsc , ichsursrc
  integer , public , pointer, dimension(:,:) :: icumbot , icumdwd , icumtop
!
  integer , public :: ibchb , ibchl , iochb , iochl , iso2 , iso4 , mixtype
  integer , public , pointer , dimension(:) :: idust
!
  real(dp) , public , pointer , dimension(:,:) :: mflx
!
  real(dp) , public , pointer , dimension(:,:,:) :: aerasp , aerext , aerssa
  real(dp) , public , pointer , dimension(:,:) :: aersrrf , aertarf
  real(dp) , public , pointer , dimension(:,:) :: aertalwrf , aersrlwrf 
!
  real(dp) , public , pointer , dimension(:,:,:) :: cemtr , cemtrac , remdrd
  real(dp) , public , pointer , dimension(:,:) :: rembc , remrat
  real(dp) , public , pointer , dimension(:,:,:,:) :: remcvc , remlsc , &
                                            rxsaq1 , rxsaq2 , rxsg

  public :: allocate_mod_trachem

  contains

  subroutine allocate_mod_trachem
    implicit none
    if ( ntr>maxntr ) then
      write (aline , *) 'In mod_trachem, resetting ntr to maxntr ', maxntr
      call say(myid)
      ntr = maxntr
    end if
    if ( nbin>maxnbin ) then
      write (aline , *) 'In mod_trachem, resetting nbin to maxnbin ', maxnbin
      call say(myid)
      nbin = maxnbin
    end if

    call getmem2d(icumbot,1,iy,1,jxp,'trachem:icumbot')
    call getmem2d(icumdwd,1,iy,1,jxp,'trachem:icumdwd')
    call getmem2d(icumtop,1,iy,1,jxp,'trachem:icumtop')
    call getmem3d(aerasp,1,iym1,1,kz,1,jxp,'trachem:aerasp')
    call getmem3d(aerext,1,iym1,1,kz,1,jxp,'trachem:aerext')
    call getmem3d(aerssa,1,iym1,1,kz,1,jxp,'trachem:aerssa')
    call getmem2d(aersrrf,1,iym1,1,jxp,'trachem:aersrrf')
    call getmem2d(aertalwrf,1,iym1,1,jxp,'trachem:aertalwrf')
    call getmem2d(aersrlwrf,1,iym1,1,jxp,'trachem:aersrlwrf')
    call getmem2d(aertarf,1,iym1,1,jxp,'trachem:aertarf')

    if ( lch ) then
      call getmem3d(cemtr,1,iy,1,jxp,1,ntr,'trachem:cemtr')
      call getmem3d(cemtrac,1,iy,1,jxp,1,ntr,'trachem:cemtrac')
      call getmem3d(remdrd,1,iy,1,jxp,1,ntr,'trachem:remdrd')
      call getmem4d(remcvc,1,iy,1,kz,1,jxp,1,ntr,'trachem:remcvc')
      call getmem4d(remlsc,1,iy,1,kz,1,jxp,1,ntr,'trachem:remlsc')
      call getmem4d(rxsaq1,1,iy,1,kz,1,jxp,1,ntr,'trachem:rxsaq1')
      call getmem4d(rxsaq2,1,iy,1,kz,1,jxp,1,ntr,'trachem:rxsaq2')
      call getmem4d(rxsg,1,iy,1,kz,1,jxp,1,ntr,'trachem:rxsg')

      allocate(chtrname(ntr))
      chtrname = ' '

      call getmem2d(chtrdpv,1,ntr,1,2,'trachem:chtrdpv')
      call getmem2d(chtrsize,1,nbin,1,2,'trachem:chtrdpv')
      call getmem1d(chtrsol,1,nbin,'trachem:chtrsol')
      call getmem2d(dustbsiz,1,nbin,1,2,'trachem:dustbsiz')
      call getmem1d(idust,1,nbin,'trachem:idust')
    end if
    call getmem2d(mflx,1,iy,1,2,'trachem:mflx')
    call getmem2d(rembc,1,iy,1,kz,'trachem:rembc')
    call getmem2d(remrat,1,iy,1,kz,'trachem:remrat')
  end subroutine allocate_mod_trachem

end module mod_trachem
