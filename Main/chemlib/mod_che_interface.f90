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

module mod_che_interface

  use m_realkinds
  use mod_dynparam
  use mod_memutil
  use mod_chem

  public
!
  real(dp) , pointer , dimension(:,:,:,:) :: chi
  real(dp) , pointer , dimension(:,:,:,:) :: chic , chiten
!
  real(dp) , pointer , dimension(:,:,:,:) :: chemsrc
  real(dp) , pointer , dimension(:,:,:,:) :: chia , chib
  real(dp) , pointer , dimension(:,:,:) :: srclp2
  real(dp) , pointer , dimension(:,:,:) :: ddsfc , dtrace , wdcvc , &
                                           wdlsc , wxaq , wxsg
  contains 

  subroutine allocate_mod_che_interface
    implicit none
    if ( lch ) then
      call getmem4d(chi,1,iy,1,kz,0,jxp+1,1,ntr,'cvaria:chi')
      call getmem4d(chic,1,iy,1,kz,1,jxp,1,ntr,'cvaria:chic')
      call getmem4d(chiten,1,iy,1,kz,1,jxp,1,ntr,'cvaria:chiten')
      call getmem4d(chemsrc,1,iy,1,jxp,1,mpy,1,ntr,'mod_che_interface:chemsrc')
      call getmem4d(chia,1,iy,1,kz,-1,jxp+2,1,ntr,'mod_che_interface:chia')
      call getmem4d(chib,1,iy,1,kz,-1,jxp+2,1,ntr,'mod_che_interface:chib')
      call getmem3d(srclp2,1,iy,1,jxp,1,ntr,'mod_che_interface:srclp2')
      call getmem3d(ddsfc,1,iy,1,jxp,1,ntr,'mod_che_interface:ddsfc')
      call getmem3d(dtrace,1,iy,1,jxp,1,ntr,'mod_che_interface:dtrace')
      call getmem3d(wdcvc,1,iy,1,jxp,1,ntr,'mod_che_interface:wdcvc')
      call getmem3d(wdlsc,1,iy,1,jxp,1,ntr,'mod_che_interface:wdlsc')
      call getmem3d(wxaq,1,iy,1,jxp,1,ntr,'mod_che_interface:wxaq')
      call getmem3d(wxsg,1,iy,1,jxp,1,ntr,'mod_che_interface:wxsg')
    end if

  end subroutine allocate_mod_che_interface
!
end module mod_che_interface
