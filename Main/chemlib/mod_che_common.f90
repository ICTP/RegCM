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

module mod_che_common

  use mod_realkinds
  use mod_dynparam
  use mod_memutil
  use mod_mpmessage
  use mod_che_param
  use mod_che_species

  public
!
  integer , parameter :: maxntr =40 
  integer , parameter :: maxnbin = 4
  integer , parameter :: maxnssl = 2
!
  integer       :: iaerosol


  real(dp) , pointer , dimension(:,:,:,:) :: chi
  real(dp) , pointer , dimension(:,:,:,:) :: chic , chiten, chemten
!
  real(dp) , pointer , dimension(:,:,:,:) :: chemsrc
  real(dp) , pointer , dimension(:,:,:,:) :: chia , chib
  real(dp) , pointer , dimension(:,:,:) :: srclp2
  real(dp) , pointer , dimension(:,:,:) :: ddsfc , dtrace , wdcvc , &
                                           wdlsc , wxaq , wxsg,drydepv
  real(dp) , pointer , dimension(:,:,:) :: taucld
  real(dp) , pointer , dimension(:,:,:,:) :: chemall
  integer , pointer , dimension(:,:) :: kcumtop , kcumbot, cveg2d
!
  character(len=5) , pointer , dimension(:) :: chtrname
!
  real(dp) , pointer , dimension(:,:) :: chtrdpv

  ! evap of l-s precip (see mod_precip.f90; [kg_h2o/kg_air/s)
  ! cum h2o vapor tendency for cum precip (kg_h2o/kg_air/s)
  real(dp) , pointer , dimension(:,:) :: chevap , checum

  real(dp) , pointer , dimension(:,:) :: chtrsize , dustbsiz
  real(dp) , pointer , dimension(:,:) :: ssltbsiz
  real(dp) , pointer , dimension(:) :: chtrsol
!
  integer , public :: ichcumtra , ichdrdepo , ichremcvc , ichremlsc , ichsursrc
!
  integer , pointer , dimension(:) :: isslt , icarb , idust
!
  real(dp) , pointer , dimension(:,:,:) :: cemtr , cemtrac , remdrd
 
  real(dp) , pointer , dimension(:,:,:,:) :: remcvc , remlsc , &
                                             rxsaq1 , rxsaq2 , rxsg

!atmospheric variable interface for chemistry 
  real(dp) :: ccalday, crdxsq
  real(dp) , pointer , dimension(:,:,:,:) ::chib3d
  real(dp) , pointer , dimension(:,:,:) :: ctb3d,cqvb3d,cubx3d,cvbx3d,crhob3d,cqcb3d,cfcc,cza,cdzq,ccldfra,crembc,cremrat
  real(dp) , pointer , dimension(:,:) ::  cpsb,ctg,clndcat,cht,cssw2da, &
                                          cvegfrac,csol2d,csdeltk2d,csdelqk2d,ctwt,cuvdrag

 
  real(dp) , pointer , dimension(:) :: hlev, cdsigma,canudg
  real(dp) , pointer , dimension(:,:) :: czen

 real(dp) :: chptop






 contains

  subroutine allocate_mod_che_common(ichem)
    implicit none

    integer , intent(in) :: ichem

    if ( ichem == 1 ) lch = .true.

    if ( ntr>maxntr ) then
      write (aline , *) 'In mod_che_common, resetting ntr to maxntr ', maxntr
      call say
      ntr = maxntr
    end if
    if ( nbin>maxnbin ) then
      write (aline , *) 'In mod_che_common, resetting nbin to maxnbin ', maxnbin
      call say
      nbin = maxnbin
    end if


    if ( lch ) then
      call getmem3d(cemtr,1,iy,1,jxp,1,ntr,'mod_che_common:cemtr')
      call getmem3d(cemtrac,1,iy,1,jxp,1,ntr,'mod_che_common:cemtrac')
      call getmem3d(remdrd,1,iy,1,jxp,1,ntr,'mod_che_common:remdrd')
      call getmem4d(remcvc,1,iy,1,kz,1,jxp,1,ntr,'mod_che_common:remcvc')
      call getmem4d(remlsc,1,iy,1,kz,1,jxp,1,ntr,'mod_che_common:remlsc')
      call getmem4d(rxsaq1,1,iy,1,kz,1,jxp,1,ntr,'mod_che_common:rxsaq1')
      call getmem4d(rxsaq2,1,iy,1,kz,1,jxp,1,ntr,'mod_che_common:rxsaq2')
      call getmem4d(rxsg,1,iy,1,kz,1,jxp,1,ntr,'mod_che_common:rxsg')

      allocate(chtrname(ntr))
      chtrname = ' '

      call getmem1d(chtrsol,1,ntr,'mod_che_common:chtrsol')
      call getmem2d(chtrdpv,1,ntr,1,2,'mod_che_common:chtrdpv')
      call getmem1d(idust,1,nbin,'mod_che_common:idust')
      call getmem1d(isslt,1,sbin,'mod_che_common:isslt')
      call getmem1d(icarb,1,5,'mod_che_common:icarb')
      call getmem2d(chtrsize,1,nbin,1,2,'mod_che_common:chtrsize')
      call getmem2d(dustbsiz,1,nbin,1,2,'mod_che_common:dustbsiz')
      call getmem2d(ssltbsiz,1,sbin,1,2,'mod_che_common:ssltbsiz')

      call getmem4d(chemall,1,iy,1,kz,1,jxp,1,totsp,'mod_che_common:chemall')
      call getmem4d(chi,1,iy,1,kz,0,jxp+1,1,ntr,'mod_che_common:chi')
      call getmem4d(chic,1,iy,1,kz,1,jxp,1,ntr,'mod_che_common:chic')
      call getmem4d(chiten,1,iy,1,kz,1,jxp,1,ntr,'mod_che_common:chiten')
      call getmem4d(chemten,1,iy,1,kz,1,jxp,1,ntr,'mod_che_common:chiten')
      call getmem4d(chemsrc,1,iy,1,jxp,1,mpy,1,ntr,'mod_che_common:chemsrc')
      call getmem4d(chia,1,iy,1,kz,-1,jxp+2,1,ntr,'mod_che_common:chia')
      call getmem4d(chib,1,iy,1,kz,-1,jxp+2,1,ntr,'mod_che_common:chib')
      call getmem3d(taucld,1,iy,1,kz,1,jxp,'mod_che_common:taucld')
      call getmem3d(srclp2,1,iy,1,jxp,1,ntr,'mod_che_common:srclp2')
      call getmem3d(ddsfc,1,iy,1,jxp,1,ntr,'mod_che_common:ddsfc')
      call getmem3d(dtrace,1,iy,1,jxp,1,ntr,'mod_che_common:dtrace')
      call getmem3d(drydepv,1,iy,1,jxp,1,ntr,'mod_che_common:drydepv')
      call getmem3d(wdcvc,1,iy,1,jxp,1,ntr,'mod_che_common:wdcvc')
      call getmem3d(wdlsc,1,iy,1,jxp,1,ntr,'mod_che_common:wdlsc')
      call getmem3d(wxaq,1,iy,1,jxp,1,ntr,'mod_che_common:wxaq')
      call getmem3d(wxsg,1,iy,1,jxp,1,ntr,'mod_che_common:wxsg')
      call getmem2d(chevap,1,iy,1,kz,'mod_che_common:chevap')
      call getmem2d(checum,1,iy,1,kz,'mod_che_common:checum')
    end if

  end subroutine allocate_mod_che_common
!
end module mod_che_common
