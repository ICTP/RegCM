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
!
  use mod_realkinds
  use mod_dynparam
  use mod_memutil
  use mod_mpmessage
  use mod_che_param
  use mod_che_species
!
  public

  real(dp) , pointer , dimension(:,:,:) :: cemtrac_io , cemtr_io , &
                                           wxaq_io , wxsg_io
  real(dp) , pointer , dimension(:,:,:,:) :: rxsaq1_io , rxsaq2_io , rxsg_io

  real(dp) , pointer , dimension(:,:,:,:) :: remcvc_io
  real(dp) , pointer , dimension(:,:,:,:) :: remlsc_io
  real(dp) , pointer , dimension(:,:,:) :: remdrd_io


  real(dp) , pointer , dimension(:,:,:,:) :: chemsrc_io
  real(dp) , pointer , dimension(:,:,:) :: ddsfc_io , dtrace_io , &
                                           wdcvc_io , wdlsc_io, drydepv_io


  real(8) , pointer , dimension(:,:,:) :: aerasp_io ,           &
                              aerext_io , aerssa_io
  real(8) , pointer , dimension(:,:) :: aersrrf_io , aertarf_io,&
                            aertalwrf_io , aersrlwrf_io


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

  real(dp) , pointer , dimension(:,:,:) :: savch0, savch_0

!
! Boundary conditions arrays
!
  real(dp)  , pointer , dimension(:,:,:,:) :: chebdy_io


  real(dp) , pointer , dimension(:,:,:) :: dustsotex_io
!
  real(dp), pointer, dimension (:,:) :: cpsa_io
!---------- DATA init section--------------------------------------------

  contains 
!
! This routines allocate all the arrays contained in the module
!
  subroutine allocate_mod_che_mppio(lband)
    implicit none
    logical , intent(in) :: lband
    integer :: mmj

    if ( lband ) then
      mmj = jx
    else
      mmj = jxm1
    end if
    if ( lch ) then

      call getmem3d(chem0,1,iy,1,ntr*kz+kz*3+ntr*8+5,1,jxp,'che_mppio:chem0')
      call getmem4d(src0,1,iy,1,mpy,1,ntr,1,jxp,'che_mppio:src0')
      call getmem3d(src1,1,iy,1,nats,1,jxp,'che_mppio:src1')
      call getmem3d(savch0,1,iy,1,kz*25,1,jxp,'che_mppio:savch0')
  
      if (myid == 0) then


        call getmem3d(chem_0,1,iy,1,ntr*kz+kz*3+ntr*8+5,1,jx,'che_mppio:chem_0')
        call getmem4d(src_0,1,iy,1,mpy,1,ntr,1,jx,'che_mppio:src_0')
        call getmem3d(src_1,1,iy,1,nats,1,jx,'che_mppio:src_1')
        call getmem3d(savch_0,1,iy,1,kz*25,1,jx,'che_mppio:savch_0')

        call getmem4d(chebdy_io,1,iy,1,kz,1,jx,1,50,'che_mppio:chebdy_io')


        call getmem2d(cpsa_io,1,iy,1,jx,'che_mppio:cpsa_io')
        call getmem2d(ssw2da_io,1,iym1,1,mmj,'che_mppio:ssw2da_io')
        call getmem2d(sdelqk2d_io,1,iym1,1,mmj,'che_mppio:sdelqk2d_io')
        call getmem2d(sdeltk2d_io,1,iym1,1,mmj,'che_mppio:sdeltk2d_io')
        call getmem2d(sfracb2d_io,1,iym1,1,mmj,'che_mppio:sfracb2d_io')
        call getmem2d(sfracs2d_io,1,iym1,1,mmj,'che_mppio:sfracs2d_io')
        call getmem2d(sfracv2d_io,1,iym1,1,mmj,'che_mppio:sfracv2d_io')
        call getmem2d(svegfrac2d_io,1,iym1,1,mmj,'che_mppio:svegfrac2d_io')

        call getmem3d(cemtrac_io,1,iy,1,jx,1,ntr,'che_mppio:cemtrac_io')
        call getmem3d(cemtr_io,1,iy,1,jx,1,ntr,'che_mppio:cemtr_io')
        call getmem3d(wxaq_io,1,iy,1,jx,1,ntr,'che_mppio:wxaq_io')
        call getmem3d(wxsg_io,1,iy,1,jx,1,ntr,'che_mppio:wxsg_io')
        call getmem4d(rxsaq1_io,1,iy,1,kz,1,jx,1,ntr,'che_mppio:rxsaq1_io')
        call getmem4d(rxsaq2_io,1,iy,1,kz,1,jx,1,ntr,'che_mppio:rxsaq2_io')
        call getmem4d(rxsg_io,1,iy,1,kz,1,jx,1,ntr,'che_mppio:rxsg_io')
        call getmem4d(remcvc_io,1,iy,1,kz,1,jx,1,ntr,'che_mppio:remcvc_io')
        call getmem4d(remlsc_io,1,iy,1,kz,1,jx,1,ntr,'che_mppio:remlsc_io')
        call getmem3d(remdrd_io,1,iy,1,jx,1,ntr,'che_mppio:remdrd_io')
        call getmem3d(drydepv_io,1,iy,1,jx,1,ntr,'che_mppio:drydepv_io')
        call getmem4d(chemsrc_io,1,iy,1,jx,1,mpy,1,ntr,'che_mppio:chemsrc_io')
        call getmem3d(ddsfc_io,1,iy,1,jx,1,ntr,'che_mppio:ddsfc_io')
        call getmem3d(dtrace_io,1,iy,1,jx,1,ntr,'che_mppio:dtrace_io')
        call getmem3d(wdcvc_io,1,iy,1,jx,1,ntr,'che_mppio:wdcvc_io')
        call getmem3d(wdlsc_io,1,iy,1,jx,1,ntr,'che_mppio:wdlsc_io')
        call getmem4d(chia_io,1,iy,1,kz,1,jx,1,ntr,'che_mppio:chia_io')
        call getmem4d(chib_io,1,iy,1,kz,1,jx,1,ntr,'che_mppio:chib_io')
        call getmem3d(dustsotex_io,1,iy,1,jx,1,nats,'che_mppio:dustsotex_io')

        call getmem3d(aerasp_io,1,iym1,1,kz,1,mmj,'che_mppio:aerasp_io')
        call getmem3d(aerext_io,1,iym1,1,kz,1,mmj,'che_mppio:aerext_io')
        call getmem3d(aerssa_io,1,iym1,1,kz,1,mmj,'che_mppio:aerssa_io')
        call getmem2d(aersrrf_io,1,iym1,1,mmj,'che_mppio:aersrrf_io')
        call getmem2d(aertarf_io,1,iym1,1,mmj,'che_mppio:aertarf_io')
        call getmem2d(aertalwrf_io,1,iym1,1,mmj,'che_mppio:aertalwrf_io')
        call getmem2d(aersrlwrf_io,1,iym1,1,mmj,'che_mppio:aersrlwrf_io')


      end if
    end if

  end subroutine allocate_mod_che_mppio
!
end module mod_che_mppio
