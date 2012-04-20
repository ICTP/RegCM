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
      integer :: i , j , k , n , ierr
 
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

          do j = 1 , jxp
            do n = 1 , ntr
              do i = 1 , iy
                chem0(i,(ntr+3)*kz+ntr*2+n,j) = wdcvc(i,j,n)
                chem0(i,(ntr+3)*kz+ntr*3+n,j) = ddsfc(i,j,n)
                chem0(i,(ntr+3)*kz+ntr*4+n,j) = wxsg(i,j,n)
                chem0(i,(ntr+3)*kz+ntr*5+n,j) = wxaq(i,j,n)
                chem0(i,(ntr+3)*kz+ntr*6+n,j) = cemtrac(i,j,n)
                chem0(i,(ntr+3)*kz+ntr*7+n,j) = drydepv(i,j,n)

              end do
            end do
          end do

         call mpi_gather(chem0,iy*((ntr+3)*kz+ntr*8+5)*jxp,            &
                        & mpi_real8,chem_0,iy*((ntr+3)*kz+ntr*8+5)*jxp, &
                        & mpi_real8,0,mpi_comm_world,ierr)


          if ( myid.eq.0 ) then
            do j = 1 , jx
              do n = 1 , ntr
                do i = 1 , iy
                  wdcvc_io(i,j,n) = chem_0(i,(ntr+3)*kz+ntr*2+n,j)
                  ddsfc_io(i,j,n) = chem_0(i,(ntr+3)*kz+ntr*3+n,j)
                  wxsg_io(i,j,n) = chem_0(i,(ntr+3)*kz+ntr*4+n,j)
                  wxaq_io(i,j,n) = chem_0(i,(ntr+3)*kz+ntr*5+n,j)
                  cemtrac_io(i,j,n) = chem_0(i,(ntr+3)*kz+ntr*6+n,j)
                  drydepv_io(i,j,n) = chem_0(i,(ntr+3)*kz+ntr*7+n,j)

                end do
              end do
            end do

            call outche2(idatex) 
              
            remlsc_io  = 0.0
            remcvc_io  = 0.0
            rxsg_io    = 0.0
            rxsaq1_io  = 0.0
            rxsaq2_io  = 0.0
            cemtr_io   = 0.0
            remdrd_io  = 0.0
            wdcvc_io   = 0.0
            ddsfc_io   = 0.0
            wxsg_io    = 0.0
            wxaq_io    = 0.0
            cemtrac_io = 0.0
            drydepv_io = 0.0
            aertarf_io = 0.0
            aersrrf_io = 0.0
            aersrlwrf_io=0.0
            aertalwrf_io=0.0
          end if

! put back to zero accumulated variables
          do n = 1 , ntr
            do j = 1 , jxp
              do k = 1 , kz
                do i = 1 , iy
                  remlsc(i,k,j,n) = 0.
                  remcvc(i,k,j,n) = 0.
                  rxsg(i,k,j,n) = 0.
                  rxsaq1(i,k,j,n) = 0.
                  rxsaq2(i,k,j,n) = 0.
                end do
              end do
            end do
          end do
          do n = 1 , ntr
            do j = 1 , jxp
              do i = 1 , iy
                cemtr(i,j,n) = 0.
                remdrd(i,j,n) = 0.
                wdcvc(i,j,n) = 0.
                ddsfc(i,j,n) = 0.
                wxsg(i,j,n) = 0.
                wxaq(i,j,n) = 0.
                cemtrac(i,j,n) = 0.
                drydepv(i,j,n) =0.
              end do
            end do
          end do

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
      integer :: ni , nj , nk , is , ie , js , je

       if ( lcband ) then
         ni = iym2
         nj = jx
         nk = kz
         is = 2
         js = 1
         ie = iym1
         je = jx
       else
         ni = iym2
         nj = jxm2
         nk = kz
         is = 2
         js = 2
         ie = iym1
         je = jxm1
       end if

       call writerec_che2(nj, ni, je, ie, nk, ntr, chia_io,     &
                wdlsc_io, wdcvc_io, ddsfc_io, cemtrac_io,    &
                drydepv_io,  aerext_io, aerssa_io, aerasp_io, &
                aertarf_io, aersrrf_io, &
                aertalwrf_io, aersrlwrf_io, cpsb_io, idatex)

        write (*,*) 'CHE variables written at ' , idatex 
        if (iaerosol > 0)  write (*,*) 'OPT variables written at ' , idatex 
      end subroutine outche2

     end module mod_che_output
