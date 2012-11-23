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

  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_mpmessage
  use mod_mppparam
  use mod_runparams
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

    subroutine output_chem(idatex)
      implicit none
      type(rcm_time_and_date) , intent(in) :: idatex
 
      call grid_collect(chia,chia_io,jce1,jce2,ice1,ice2,1,kz,1,ntr)
      call grid_collect(cpsb,cpsb_io,jce1,jce2,ice1,ice2)
      call grid_collect(dtrace,dtrace_io,jce1,jce2,ice1,ice2,1,ntr)
      call grid_collect(wdlsc,wdlsc_io,jce1,jce2,ice1,ice2,1,ntr)
      call grid_collect(wdcvc,wdcvc_io,jce1,jce2,ice1,ice2,1,ntr)
      call grid_collect(remdrd,remdrd_io,jce1,jce2,ice1,ice2,1,ntr)
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
        if ( ibltyp == 2 .or. ibltyp == 99 ) then
          call grid_collect(cchifxuw,ccuwdiag_io,jci1,jci2,ici1,ici2,1,ntr)
        end if
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
!
    end subroutine output_chem
!----------------------------------------------------------------
!================================================================

    subroutine outche2(idatex)
      implicit none
      type(rcm_time_and_date) , intent(in) :: idatex

      if ( myid == iocpu ) then 
        call writerec_che2(chia_io,dtrace_io,wdlsc_io,wdcvc_io,remdrd_io,  &
                           cemtrac_io,drydepv_io,chemdiag_io,cadvhdiag_io, &
                           cadvvdiag_io,cdifhdiag_io,cconvdiag_io,         &
                           cbdydiag_io,ctbldiag_io,cseddpdiag_io,          &
                           ccuwdiag_io,remlsc_io,remcvc_io,cpsb_io,idatex)
      end if
    end subroutine outche2

end module mod_che_output
