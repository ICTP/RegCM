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

module mod_che_mppio

  use m_realkinds
  use mod_dynparam
  use mod_memutil
  use mod_message
  use mod_chem
!
  real(dp) , pointer , dimension(:,:,:) :: cemtrac_io , cemtr_io , &
                                          wxaq_io , wxsg_io
  real(dp) , pointer , dimension(:,:,:,:) :: rxsaq1_io , rxsaq2_io , rxsg_io

  real(dp) , private , pointer , dimension(:,:,:) :: chebat

  real(dp) , pointer , dimension(:,:,:,:) :: remcvc_io
  real(dp) , pointer , dimension(:,:,:,:) :: remlsc_io
  real(dp) , pointer , dimension(:,:,:) :: remdrd_io

  real(dp) , pointer , dimension(:,:,:,:) :: chemsrc_io
  real(dp) , pointer , dimension(:,:,:) :: ddsfc_io , dtrace_io , &
                                          wdcvc_io , wdlsc_io

  real(dp) , pointer , dimension(:,:) :: ssw2da_io , sdeltk2d_io ,  &
                         sdelqk2d_io , sfracv2d_io , sfracb2d_io , &
                         sfracs2d_io , svegfrac2d_io

  real(dp) , pointer , dimension(:,:,:,:) :: chia_io , chib_io

  real(dp) , pointer , dimension(:,:,:) :: chem0
  real(dp) , pointer , dimension(:,:,:) :: chem_0
  real(dp) , pointer , dimension(:,:,:,:) :: src0
  real(dp) , pointer , dimension(:,:,:,:) :: src_0
  real(dp) , pointer , dimension(:,:,:) :: src1
  real(dp) , pointer , dimension(:,:,:) :: src_1
!
  real(dp) , pointer , dimension(:,:,:) :: dustsotex_io
!

!---------- DATA init section--------------------------------------------

  contains 
!
!     This routines allocate all the arrays contained in the module
!
  subroutine allocate_mod_che_mppio(lband)
    implicit none
    logical , intent(in) :: lband

    if ( lch ) then

      call getmem3d(chem0,1,iy,1,ntr*kz+kz*3+ntr*7+5, &
                          1,jxp,'mod_che_mppio:chem0')
      call getmem4d(src0,1,iy,1,mpy,1,ntr,1,jxp,'mod_che_mppio:src0')
      call getmem3d(src1,1,iy,1,nats,1,jxp,'mod_che_mppio:src1')

      if (myid == 0) then
        call getmem3d(chem_0,1,iy,1,ntr*kz+kz*3+ntr*7+5, &
                             1,jx,'mod_che_mppio:chem_0')
        call getmem4d(src_0,1,iy,1,mpy,1,ntr,1,jx,'mod_che_mppio:src_0')
        call getmem3d(src_1,1,iy,1,nats,1,jx,'mod_che_mppio:src_1')
      end if

      if (lband) then
        call getmem3d(chebat,1,iym1,1,jx,1,7,'mod_che_mppio:chebat')
      else
        call getmem3d(chebat,1,iym1,1,jxm1,1,7,'mod_che_mppio:chebat')
      end if
      ssw2da_io     => chebat(:,:,1)
      sdelqk2d_io   => chebat(:,:,2)
      sdeltk2d_io   => chebat(:,:,3)
      sfracb2d_io   => chebat(:,:,4)
      sfracs2d_io   => chebat(:,:,5)
      sfracv2d_io   => chebat(:,:,6)
      svegfrac2d_io => chebat(:,:,7)

      call getmem3d(cemtrac_io,1,iy,1,jx,1,ntr,'mod_che_mppio:cemtrac_io')
      call getmem3d(cemtr_io,1,iy,1,jx,1,ntr,'mod_che_mppio:cemtr_io')
      call getmem3d(wxaq_io,1,iy,1,jx,1,ntr,'mod_che_mppio:wxaq_io')
      call getmem3d(wxsg_io,1,iy,1,jx,1,ntr,'mod_che_mppio:wxsg_io')
      call getmem4d(rxsaq1_io,1,iy,1,kz,1,jx,1,ntr,'mod_che_mppio:rxsaq1_io')
      call getmem4d(rxsaq2_io,1,iy,1,kz,1,jx,1,ntr,'mod_che_mppio:rxsaq2_io')
      call getmem4d(rxsg_io,1,iy,1,kz,1,jx,1,ntr,'mod_che_mppio:rxsg_io')
      call getmem4d(remcvc_io,1,iy,1,kz,1,jx,1,ntr,'mod_che_mppio:remcvc_io')
      call getmem4d(remlsc_io,1,iy,1,kz,1,jx,1,ntr,'mod_che_mppio:remlsc_io')
      call getmem3d(remdrd_io,1,iy,1,jx,1,ntr,'mod_che_mppio:remdrd_io')
      call getmem4d(chemsrc_io,1,iy,1,jx,1,mpy,1,ntr,'mod_che_mppio:chemsrc_io')
      call getmem3d(ddsfc_io,1,iy,1,jx,1,ntr,'mod_che_mppio:ddsfc_io')
      call getmem3d(dtrace_io,1,iy,1,jx,1,ntr,'mod_che_mppio:dtrace_io')
      call getmem3d(wdcvc_io,1,iy,1,jx,1,ntr,'mod_che_mppio:wdcvc_io')
      call getmem3d(wdlsc_io,1,iy,1,jx,1,ntr,'mod_che_mppio:wdlsc_io')
      call getmem4d(chia_io,1,iy,1,kz,1,jx,1,ntr,'mod_che_mppio:chia_io')
      call getmem4d(chib_io,1,iy,1,kz,1,jx,1,ntr,'mod_che_mppio:chib_io')
      call getmem3d(dustsotex_io,1,iy,1,jx,1,nats,'mod_che_mppio:dustsotex_io')
    end if

  end subroutine allocate_mod_che_mppio
!
end module mod_che_mppio
