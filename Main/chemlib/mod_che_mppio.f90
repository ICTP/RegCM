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
  use mod_che_common
  use mod_che_species
!
  public

  real(dp) , pointer , dimension(:,:,:) :: cemtrac_io , wxaq_io , wxsg_io
  real(dp) , pointer , dimension(:,:,:) :: remdrd_io
  real(dp) , pointer , dimension(:,:,:,:) :: remlsc_io , remcvc_io
  real(dp) , pointer , dimension(:,:,:) :: ddsfc_io , dtrace_io , &
                                           wdcvc_io , wdlsc_io, drydepv_io
  real(8) , pointer , dimension(:,:,:) :: aerasp_io , aerext_io , aerssa_io
  real(8) , pointer , dimension(:,:) :: aersrrf_io , aertarf_io , &
                                        aertalwrf_io , aersrlwrf_io , aeraod_io
  real(dp) , pointer , dimension(:,:) :: ssw2da_io , sdeltk2d_io ,   &
                                         sdelqk2d_io , sfracv2d_io , &
                                         sfracb2d_io , sfracs2d_io , &
                                         svegfrac2d_io
  real(dp) , pointer , dimension(:,:,:) :: chemsrc_io
  real(dp) , pointer , dimension(:,:,:,:) :: chia_io , chib_io, chemdiag_io
!
! Boundary conditions arrays
!
  real(dp) , pointer , dimension(:,:,:,:) :: chebdy_in , chebdy_io0 , &
                                             chebdy_io1 , oxcl_io
  real(dp) , pointer , dimension(:,:,:) :: dustsotex_io
!
  real(dp), pointer, dimension (:,:) :: cpsb_io

!---------- DATA init section--------------------------------------------

  contains 
    !
    ! This routines allocate all the arrays contained in the module
    !
    subroutine allocate_mod_che_mppio(ilcband)
      implicit none
      logical , intent(in) :: ilcband

      lcband = ilcband

      if ( lch ) then
        if (myid == 0) then
          call getmem4d(chebdy_io0,jdot1,jdot2,idot1,idot2, &
                        1,kz,1,ntr,'che_mppio:chebdy_io')
          call getmem4d(chebdy_io1,jdot1,jdot2,idot1,idot2, &
                        1,kz,1,ntr,'che_mppio:chebdy_io')
          call getmem4d(chebdy_in,jdot1,jdot2,idot1,idot2, &
                        1,kz,1,50,'che_mppio:chebdy_in')
          if ( ioxclim == 1 ) then 
            call getmem4d(oxcl_io,jdot1,jdot2,idot1,idot2, &
                          1,kz,1,5,'che_mppio:oxcl_io')
          end if
          call getmem3d(chemsrc_io,jdot1,jdot2,idot1,idot2, &
                        1,ntr,'che_mppio:chemsrc_io')
  
          call getmem3d(dtrace_io,jcross1,jcross2,icross1,icross2,1,ntr, &
                        'che_mppio:dtrace_io')
          call getmem3d(wdlsc_io,jcross1,jcross2,icross1,icross2,1,ntr, &
                        'che_mppio:wdlsc_io')
          call getmem3d(wdcvc_io,jcross1,jcross2,icross1,icross2,1,ntr, &
                        'che_mppio:wdcvc_io')
          call getmem3d(ddsfc_io,jcross1,jcross2,icross1,icross2,1,ntr, &
                        'che_mppio:ddsfc_io')
          call getmem3d(wxsg_io,jcross1,jcross2,icross1,icross2,1,ntr, &
                        'che_mppio:wxsg_io')
          call getmem3d(wxaq_io,jcross1,jcross2,icross1,icross2,1,ntr, &
                        'che_mppio:wxaq_io')
          call getmem3d(cemtrac_io,jcross1,jcross2,icross1,icross2,1,ntr, &
                        'che_mppio:cemtrac_io')
          call getmem3d(drydepv_io,jcross1,jcross2,icross1,icross2,1,ntr, &
                        'che_mppio:drydepv_io')

          call getmem4d(remlsc_io,jcross1,jcross2,icross1,icross2, &
                        1,kz,1,ntr,'che_mppio:remlsc_io')
          call getmem4d(remcvc_io,jcross1,jcross2,icross1,icross2, &
                        1,kz,1,ntr,'che_mppio:remcvc_io')
          call getmem3d(remdrd_io,jcross1,jcross2,icross1,icross2, &
                        1,ntr,'che_mppio:remdrd_io')
          call getmem2d(ssw2da_io,jcross1,jcross2,icross1,icross2, &
                        'che_mppio:ssw2da_io')
          call getmem2d(sdelqk2d_io,jcross1,jcross2,icross1,icross2, &
                        'che_mppio:sdelqk2d_io')
          call getmem2d(sdeltk2d_io,jcross1,jcross2,icross1,icross2, &
                        'che_mppio:sdeltk2d_io')
          call getmem2d(sfracb2d_io,jcross1,jcross2,icross1,icross2, &
                        'che_mppio:sfracb2d_io')
          call getmem2d(sfracs2d_io,jcross1,jcross2,icross1,icross2, &
                        'che_mppio:sfracs2d_io')
          call getmem2d(sfracv2d_io,jcross1,jcross2,icross1,icross2, &
                        'che_mppio:sfracv2d_io')
          call getmem2d(svegfrac2d_io,jcross1,jcross2,icross1,icross2, &
                        'che_mppio:svegfrac2d_io')

          call getmem2d(cpsb_io,jcross1,jcross2,icross1,icross2, &
                        'che_mppio:cpsb_io')
          call getmem4d(chia_io,jcross1,jcross2,icross1,icross2, &
                        1,kz,1,ntr,'che_mppio:chia_io')
          call getmem4d(chib_io,jcross1,jcross2,icross1,icross2, &
                        1,kz,1,ntr,'che_mppio:chib_io')

          call getmem4d(chemdiag_io,jcross1,jcross2,icross1,icross2, &
                        1,kz,1,ntr,'che_mppio:chemdiag_io')


          call getmem3d(dustsotex_io,jdot1,jdot2,idot1,idot2, &
                        1,nats,'che_mppio:dustsotex_io')

          call getmem3d(aerasp_io,jcross1,jcross2,icross1,icross2, &
                        1,kz,'che_mppio:aerasp_io')
          call getmem3d(aerext_io,jcross1,jcross2,icross1,icross2, &
                        1,kz,'che_mppio:aerext_io')
          call getmem3d(aerssa_io,jcross1,jcross2,icross1,icross2, &
                        1,kz,'che_mppio:aerssa_io')
          call getmem2d(aersrrf_io,jcross1,jcross2,icross1,icross2, &
                        'che_mppio:aersrrf_io')
          call getmem2d(aertarf_io,jcross1,jcross2,icross1,icross2, &
                          'che_mppio:aertarf_io')
          call getmem2d(aertalwrf_io,jcross1,jcross2,icross1,icross2, &
                        'che_mppio:aertalwrf_io')
          call getmem2d(aersrlwrf_io,jcross1,jcross2,icross1,icross2, &
                        'che_mppio:aersrlwrf_io')
          call getmem2d(aeraod_io,jcross1,jcross2,icross1,icross2, &
                        'che_mppio:aeraod_io')
        end if
      end if
    end subroutine allocate_mod_che_mppio
!
end module mod_che_mppio
