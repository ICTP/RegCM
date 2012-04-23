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
      use mpi
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
 
      call deco1_gather(chia,chia_io,jcross1,jcross2,icross1,icross2,1,kz,1,ntr)
      call deco1_gather(cpsb,cpsb_io,jcross1,jcross2,icross1,icross2)
      call deco1_gather(aerext,aerext_io,jcross1,jcross2,icross1,icross2,1,kz)
      call deco1_gather(aerssa,aerssa_io,jcross1,jcross2,icross1,icross2,1,kz)
      call deco1_gather(aerasp,aerasp_io,jcross1,jcross2,icross1,icross2,1,kz)
      call deco1_gather(aertarf,aertarf_io,jcross1,jcross2,icross1,icross2)
      call deco1_gather(aersrrf,aersrrf_io,jcross1,jcross2,icross1,icross2)
      call deco1_gather(aertalwrf,aertalwrf_io,jcross1,jcross2,icross1,icross2)
      call deco1_gather(aersrlwrf,aersrlwrf_io,jcross1,jcross2,icross1,icross2)
      call deco1_gather(dtrace,dtrace_io,jcross1,jcross2,icross1,icross2,1,ntr)
      call deco1_gather(wdlsc,wdlsc_io,jcross1,jcross2,icross1,icross2,1,ntr)
      call deco1_gather(wdcvc,wdcvc_io,jcross1,jcross2,icross1,icross2,1,ntr)
      call deco1_gather(ddsfc,ddsfc_io,jcross1,jcross2,icross1,icross2,1,ntr)
      call deco1_gather(wxsg,ddsfc_io,jcross1,jcross2,icross1,icross2,1,ntr)
      call deco1_gather(wxaq,ddsfc_io,jcross1,jcross2,icross1,icross2,1,ntr)
      call deco1_gather(cemtrac,cemtrac_io,jcross1,jcross2, &
                        icross1,icross2,1,ntr)
      call deco1_gather(drydepv,drydepv_io,jcross1,jcross2, &
                        icross1,icross2,1,ntr)

      call outche2(idatex) 
              
      ! put back to zero accumulated variables

      remlsc(:,:,:,:) = d_zero
      remcvc(:,:,:,:) = d_zero
      rxsg(:,:,:,:) = d_zero
      rxsaq1(:,:,:,:) = d_zero
      rxsaq2(:,:,:,:) = d_zero
      cemtr(:,:,:) = d_zero
      remdrd(:,:,:) = d_zero

      drydepv(:,:,:) = d_zero
      cemtrac(:,:,:) = d_zero
      wxaq(:,:,:) = d_zero
      wxsg(:,:,:) = d_zero
      ddsfc(:,:,:) = d_zero
      wdcvc(:,:,:) = d_zero
      wdlsc(:,:,:) = d_zero
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

      if (myid ==0) then 
      call writerec_che2(chia_io,wdlsc_io,wdcvc_io,ddsfc_io,cemtrac_io,   &
                         drydepv_io,aerext_io,aerssa_io,aerasp_io,        &
                         aertarf_io,aersrrf_io,aertalwrf_io,aersrlwrf_io, &
                         cpsb_io,idatex)

      write (*,*) 'CHE variables written at ' , tochar(idatex) 
      if ( iaerosol > 0 ) then
        write (*,*) 'OPT variables written at ' , tochar(idatex)
      end if

     end if
    end subroutine outche2

end module mod_che_output
