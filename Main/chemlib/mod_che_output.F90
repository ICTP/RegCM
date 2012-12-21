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
  use mod_runparams
  use mod_dynparam
  use mod_che_param
  use mod_che_common
  use mod_outvars

  private

  public :: fill_chem_outvars

  real(rk8) , parameter :: cfd = 1.D6 * 86400.D0

  contains

    subroutine fill_chem_outvars(itr)
      implicit none
      integer(ik4) , intent(in) :: itr
      integer(ik4) :: k

      if ( associated(che_wdrflx_out) ) then
        che_wdrflx_out = wdlsc(jci1:jci2,ici1:ici2,itr)*cfd
      end if
      wdlsc(:,:,itr) = d_zero
      if ( associated(che_wdcflx_out) ) then
        che_wdcflx_out = wdcvc(jci1:jci2,ici1:ici2,itr)*cfd
      end if
      wdcvc(:,:,itr) = d_zero
      if ( associated(che_ddflx_out) ) then
        che_ddflx_out = remdrd(jci1:jci2,ici1:ici2,itr)*cfd
      end if
      remdrd(:,:,itr) = d_zero
      if ( associated(che_emflx_out) ) then
        che_emflx_out = cemtrac(jci1:jci2,ici1:ici2,itr)*cfd
      end if
      cemtrac(:,:,itr) = d_zero
      if ( associated(che_ddvel_out) ) then
        che_ddvel_out = drydepv(jci1:jci2,ici1:ici2,itr)*cfdout
      end if
      drydepv(:,:,itr) = d_zero
      if ( associated(che_mixrat_out) ) then
        do k = 1 , kz
          che_mixrat_out(:,:,k) = chia(jci1:jci2,ici1:ici2,k,itr) / &
                           cpsb(jci1:jci2,ici1:ici2)
        end do
      end if
      if ( associated(che_burden_out) ) then
        che_burden_out = dtrace(jci1:jci2,ici1:ici2,itr)*1.0D6
      end if
      if ( ichdiag > 0 ) then
        if ( associated(che_cheten_out) ) then
          do k = 1 , kz
            che_cheten_out(:,:,k) = chemdiag(jci1:jci2,ici1:ici2,k,itr) / &
                             cpsb(jci1:jci2,ici1:ici2)
          end do
        end if
        chemdiag(:,:,:,itr) = d_zero
        if ( associated(che_advhten_out) ) then
          do k = 1 , kz
            che_advhten_out(:,:,k) = cadvhdiag(jci1:jci2,ici1:ici2,k,itr) / &
                             cpsb(jci1:jci2,ici1:ici2)
          end do
        end if
        cadvhdiag(:,:,:,itr) = d_zero
        if ( associated(che_advvten_out) ) then
          do k = 1 , kz
            che_advvten_out(:,:,k) = cadvvdiag(jci1:jci2,ici1:ici2,k,itr) / &
                             cpsb(jci1:jci2,ici1:ici2)
          end do
        end if
        cadvvdiag(:,:,:,itr) = d_zero
        if ( associated(che_difhten_out) ) then
          do k = 1 , kz
            che_difhten_out(:,:,k) = cdifhdiag(jci1:jci2,ici1:ici2,k,itr) / &
                             cpsb(jci1:jci2,ici1:ici2)
          end do
        end if
        cdifhdiag(:,:,:,itr) = d_zero
        if ( associated(che_cuten_out) ) then
          do k = 1 , kz
            che_cuten_out(:,:,k) = cconvdiag(jci1:jci2,ici1:ici2,k,itr) / &
                             cpsb(jci1:jci2,ici1:ici2)
          end do
        end if
        cconvdiag(:,:,:,itr) = d_zero
        if ( associated(che_tuten_out) ) then
          do k = 1 , kz
            che_tuten_out(:,:,k) = ctbldiag(jci1:jci2,ici1:ici2,k,itr) / &
                             cpsb(jci1:jci2,ici1:ici2)
          end do
        end if
        ctbldiag (:,:,:,itr) = d_zero 
        if ( associated(che_raiten_out) ) then
          do k = 1 , kz
            che_raiten_out(:,:,k) = remcvc(jci1:jci2,ici1:ici2,k,itr) / &
                             cpsb(jci1:jci2,ici1:ici2)
          end do
        end if
        if ( associated(che_wasten_out) ) then
          do k = 1 , kz
            che_wasten_out(:,:,k) = remlsc(jci1:jci2,ici1:ici2,k,itr) / &
                             cpsb(jci1:jci2,ici1:ici2)
          end do
        end if
        if ( associated(che_bdyten_out) ) then
          do k = 1 , kz
            che_bdyten_out(:,:,k) = cbdydiag(jci1:jci2,ici1:ici2,k,itr) / &
                             cpsb(jci1:jci2,ici1:ici2)
          end do
        end if
        cbdydiag(:,:,:,itr) = d_zero
        if ( associated(che_sedten_out) ) then
          do k = 1 , kz
            che_sedten_out(:,:,k) = cseddpdiag(jci1:jci2,ici1:ici2,k,itr) / &
                             cpsb(jci1:jci2,ici1:ici2)
          end do
        end if
        cseddpdiag(:,:,:,itr) = d_zero      
        if ( associated(che_pblten_out) ) then
          che_pblten_out = cchifxuw(jci1:jci2,ici1:ici2,itr)
        end if
       if ( associated(che_emten_out) ) then
          do k = 1 , kz
            che_emten_out(:,:,k) = cemisdiag(jci1:jci2,ici1:ici2,k,itr) / &
                             cpsb(jci1:jci2,ici1:ici2)
          end do
        end if
       cemisdiag(:,:,:,itr) = d_zero
      endif


      remlsc(:,:,:,itr) = d_zero
      remcvc(:,:,:,itr) = d_zero
      rxsg(:,:,:,itr) = d_zero
      rxsaq1(:,:,:,itr) = d_zero
      rxsaq2(:,:,:,itr) = d_zero
      wxaq(:,:,itr) = d_zero
      wxsg(:,:,itr) = d_zero

    end subroutine fill_chem_outvars

end module mod_che_output
