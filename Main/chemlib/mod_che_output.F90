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

module mod_che_output

  use mod_constants
  use mod_mpmessage
  use mod_service 
  use mod_dynparam
  use mod_che_common
  use mod_che_param
  use mod_che_mppio
  use mod_che_ncio
  use mod_che_indices

  private

  public :: output_chem

  contains

    subroutine output_chem(idatex,aerext,aerssa,aerasp,aertarf, &
                           aersrrf,aertalwrf,aersrlwrf)
      implicit none
      type(rcm_time_and_date) , intent(in) :: idatex
      ! note : aerosol optical properties variables are defined in RAD module
      ! but are outputed if CHE is active !
      ! that is why they need to pass as argument here (an other option
      ! would be to define an interface)!! 
      real(dp) , pointer , dimension(:,:,:)  :: aerext , aerssa , &
                                              aerasp
      real(dp) , pointer , dimension(:,:)  :: aertarf , aersrrf , &
                                              aertalwrf , aersrlwrf            
 
      call grid_collect(chia,chia_io,jce1,jce2,ice1,ice2,1,kz,1,ntr)
      call grid_collect(cpsb,cpsb_io,jce1,jce2,ice1,ice2)
      call grid_collect(aerext,aerext_io,jci1,jci2,ici1,ici2,1,kz)
      call grid_collect(aerssa,aerssa_io,jci1,jci2,ici1,ici2,1,kz)
      call grid_collect(aerasp,aerasp_io,jci1,jci2,ici1,ici2,1,kz)
      call grid_collect(aertarf,aertarf_io,jci1,jci2,ici1,ici2)
      call grid_collect(aersrrf,aersrrf_io,jci1,jci2,ici1,ici2)
      call grid_collect(aertalwrf,aertalwrf_io,jci1,jci2,ici1,ici2)
      call grid_collect(aersrlwrf,aersrlwrf_io,jci1,jci2,ici1,ici2)
      call grid_collect(dtrace,dtrace_io,jce1,jce2,ice1,ice2,1,ntr)
      call grid_collect(wdlsc,wdlsc_io,jce1,jce2,ice1,ice2,1,ntr)
      call grid_collect(wdcvc,wdcvc_io,jce1,jce2,ice1,ice2,1,ntr)
      call grid_collect(ddsfc,ddsfc_io,jce1,jce2,ice1,ice2,1,ntr)
      call grid_collect(cemtrac,cemtrac_io,jce1,jce2,ice1,ice2,1,ntr)
      call grid_collect(drydepv,drydepv_io,jce1,jce2,ice1,ice2,1,ntr)

      if (ichdiag > 0) then   
        call grid_collect(chemdiag,chemdiag_io,jce1,jce2,ice1,ice2,1,kz,1,ntr)
        call grid_collect(cadvhdiag,cadvhdiag_io,jce1,jce2,ice1,ice2,1,kz,1,ntr)
        call grid_collect(cadvvdiag,cadvvdiag_io,jce1,jce2,ice1,ice2,1,kz,1,ntr)
        call grid_collect(cdifhdiag,cdifhdiag_io,jce1,jce2,ice1,ice2,1,kz,1,ntr)
        call grid_collect(cconvdiag,cconvdiag_io,jce1,jce2,ice1,ice2,1,kz,1,ntr)
        call grid_collect(cbdydiag,cbdydiag_io,jce1,jce2,ice1,ice2,1,kz,1,ntr)
        call grid_collect(ctbldiag,ctbldiag_io,jce1,jce2,ice1,ice2,1,kz,1,ntr)
        call grid_collect(remcvc,remcvc_io,jce1,jce2,ice1,ice2,1,kz,1,ntr)
        call grid_collect(remlsc,remlsc_io,jce1,jce2,ice1,ice2,1,kz,1,ntr)
        call grid_collect(cseddpdiag,cseddpdiag_io, &
                          jce1,jce2,ice1,ice2,1,kz,1,ntr)
      end if

      call outche2(idatex) 

      ! put back to zero accumulated variables

 
      remdrd(:,:,:) = d_zero
      drydepv(:,:,:) = d_zero
      cemtrac(:,:,:) = d_zero

      remlsc(:,:,:,:) = d_zero
      remcvc(:,:,:,:) = d_zero
      rxsg(:,:,:,:) = d_zero
      rxsaq1(:,:,:,:) = d_zero
      rxsaq2(:,:,:,:) = d_zero

      wxaq(:,:,:) = d_zero
      wxsg(:,:,:) = d_zero
      ddsfc(:,:,:) = d_zero
      wdcvc(:,:,:) = d_zero
      wdlsc(:,:,:) = d_zero

      if ( ichdiag == 1 ) then 
        chemdiag(:,:,:,:) = d_zero
        cadvhdiag(:,:,:,:) = d_zero
        cadvvdiag(:,:,:,:) = d_zero
        cdifhdiag(:,:,:,:) = d_zero
        cconvdiag(:,:,:,:) = d_zero
        cbdydiag(:,:,:,:) = d_zero
        ctbldiag (:,:,:,:) = d_zero 
        cseddpdiag(:,:,:,:) = d_zero      
       end if

      aertarf(:,:) = d_zero
      aersrrf(:,:) = d_zero
      aertalwrf(:,:) = d_zero
      aersrlwrf(:,:) = d_zero
!
    end subroutine output_chem
!----------------------------------------------------------------
!================================================================

    subroutine outche2(idatex)
      implicit none
      type(rcm_time_and_date) , intent(in) :: idatex

      if ( myid == 0 ) then 
        aeraod_io = sum(aerext_io,3)
        call writerec_che2(chia_io,dtrace_io,wdlsc_io,wdcvc_io,ddsfc_io,   &
                           cemtrac_io,drydepv_io,chemdiag_io,cadvhdiag_io, &
                           cadvvdiag_io,cdifhdiag_io,cconvdiag_io,         &
                           cbdydiag_io,ctbldiag_io,cseddpdiag_io,remlsc_io,remcvc_io,    &
                           aerext_io,aerssa_io,aerasp_io,aeraod_io,        &
                           aertarf_io,aersrrf_io,aertalwrf_io,aersrlwrf_io,&
                           cpsb_io,idatex)
        write (*,*) 'CHE variables written at ' , tochar(idatex) 
        if ( iaerosol > 0 ) then
          write (*,*) 'OPT variables written at ' , tochar(idatex)
        end if
      end if
    end subroutine outche2

end module mod_che_output
