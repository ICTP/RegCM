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
!
  

  private
  public :: output_chem
  contains

!!$!================================================================
  subroutine output_chem (idatex,aerext, aerssa, aerasp, aertarf,aersrrf,aertalwrf, aersrlwrf )
 
 use mpi
 implicit none
 type(rcm_time_and_date) , intent(in) :: idatex
!(jxp, iym1,kz)
! note : aerosol optical properties variables are defined in RAD module but are outputed if CHE is active !
! that is why they need to pass as argument here (an other option would be to define an interface)!! 


 real(dp), dimension(:,:,:), intent(inout) ::   aerext, aerssa, aerasp
 real(dp), dimension(:,:) , intent(inout) :: aertarf,aersrrf,aertalwrf, aersrlwrf            

 integer :: i,j,k,n,ierr
 
 print*, 'FAB, outputchem',myid,  maxval(aerext) 
           
          do j = 1 , jxp
            do n = 1 , ntr
              do k = 1 , kz
                do i = 1 , iy
                  chem0(i,(n-1)*kz+k,j) = chia(j,i,k,n)
                end do
              end do
            end do
          end do
          do j = 1 , jendx
            do k = 1 , kz
              do i = 1 , iym1
                chem0(i,ntr*kz+k,j) = aerext(j,i,k)
                chem0(i,ntr*kz+kz+k,j) = aerssa(j,i,k)
                chem0(i,ntr*kz+kz*2+k,j) = aerasp(j,i,k)
              end do
            end do
          end do
          do j = 1 , jxp
            do n = 1 , ntr
              do i = 1 , iy
                chem0(i,(ntr+3)*kz+n,j) = dtrace(i,j,n)
                chem0(i,(ntr+3)*kz+ntr+n,j) = wdlsc(i,j,n)
                chem0(i,(ntr+3)*kz+ntr*2+n,j) = wdcvc(i,j,n)
                chem0(i,(ntr+3)*kz+ntr*3+n,j) = ddsfc(i,j,n)
                chem0(i,(ntr+3)*kz+ntr*4+n,j) = wxsg(i,j,n)
                chem0(i,(ntr+3)*kz+ntr*5+n,j) = wxaq(i,j,n)
                chem0(i,(ntr+3)*kz+ntr*6+n,j) = cemtrac(i,j,n)
                chem0(i,(ntr+3)*kz+ntr*7+n,j) = drydepv(i,j,n)

              end do
            end do
          end do
          do j = 1 , jendx
            do i = 1 , iym1
              chem0(i,(ntr+3)*kz+ntr*8+1,j) = aertarf(j,i)
              chem0(i,(ntr+3)*kz+ntr*8+2,j) = aersrrf(j,i)
              chem0(i,(ntr+3)*kz+ntr*8+3,j) = aertalwrf(j,i)
              chem0(i,(ntr+3)*kz+ntr*8+4,j) = aersrlwrf(j,i)             

            end do
          end do
          do j = 1 , jxp
            do i = 1 , iy
              chem0(i,(ntr+3)*kz+ntr*8+5,j) = cpsb(j,i)
            end do
          end do
 

         call mpi_gather(chem0,iy*((ntr+3)*kz+ntr*8+5)*jxp,            &
                        & mpi_real8,chem_0,iy*((ntr+3)*kz+ntr*8+5)*jxp, &
                        & mpi_real8,0,mpi_comm_world,ierr)


          if ( myid.eq.0 ) then
            do j = 1 , jx
              do n = 1 , ntr
                do k = 1 , kz
                  do i = 1 , iy
                    chia_io(i,k,j,n) = chem_0(i,(n-1)*kz+k,j)
                  end do
                end do
              end do
            end do
#ifdef BAND
            do j = 1 , jx
              do k = 1 , kz
                do i = 1 , iym1
                  aerext_io(i,k,j) = chem_0(i,ntr*kz+k,j)
                  aerssa_io(i,k,j) = chem_0(i,ntr*kz+kz+k,j)
                  aerasp_io(i,k,j) = chem_0(i,ntr*kz+kz*2+k,j)
#else
            do j = 1 , jxm1
              do k = 1 , kz
                do i = 1 , iym1
                  aerext_io(i,k,j) = chem_0(i,ntr*kz+k,j+1)
                  aerssa_io(i,k,j) = chem_0(i,ntr*kz+kz+k,j+1)
                  aerasp_io(i,k,j) = chem_0(i,ntr*kz+kz*2+k,j+1)
#endif
                end do
              end do
            end do
            do j = 1 , jx
              do n = 1 , ntr
                do i = 1 , iy
                  dtrace_io(i,j,n) = chem_0(i,(ntr+3)*kz+n,j)
                  wdlsc_io(i,j,n) = chem_0(i,(ntr+3)*kz+ntr+n,j)
                  wdcvc_io(i,j,n) = chem_0(i,(ntr+3)*kz+ntr*2+n,j)
                  ddsfc_io(i,j,n) = chem_0(i,(ntr+3)*kz+ntr*3+n,j)
                  wxsg_io(i,j,n) = chem_0(i,(ntr+3)*kz+ntr*4+n,j)
                  wxaq_io(i,j,n) = chem_0(i,(ntr+3)*kz+ntr*5+n,j)
                  cemtrac_io(i,j,n) = chem_0(i,(ntr+3)*kz+ntr*6+n,j)
                  drydepv_io(i,j,n) = chem_0(i,(ntr+3)*kz+ntr*7+n,j)

                end do
              end do
            end do
#ifdef BAND
            do j = 1 , jx
              do i = 1 , iym1
                aertarf_io(i,j) = chem_0(i,(ntr+3)*kz+ntr*8+1,j)
                aersrrf_io(i,j) = chem_0(i,(ntr+3)*kz+ntr*8+2,j)
                aertalwrf_io(i,j) = chem_0(i,(ntr+3)*kz+ntr*8+3,j)
                aersrlwrf_io(i,j) = chem_0(i,(ntr+3)*kz+ntr*8+4,j)
#else
            do j = 1 , jxm1
              do i = 1 , iym1
                aertarf_io(i,j) = chem_0(i,(ntr+3)*kz+ntr*8+1,j+1)
                aersrrf_io(i,j) = chem_0(i,(ntr+3)*kz+ntr*8+2,j+1)
                aertalwrf_io(i,j) = chem_0(i,(ntr+3)*kz+ntr*8+3,j+1)
                aersrlwrf_io(i,j) = chem_0(i,(ntr+3)*kz+ntr*8+4,j+1)
#endif
              end do
            end do


            do j = 1 , jx
              do i = 1 , iy
                cpsa_io(i,j) = chem_0(i,(ntr+3)*kz+ntr*8+5,j)
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
            wdlsc_io   = 0.0
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
                wdlsc(i,j,n) = 0.
                wdcvc(i,j,n) = 0.
                ddsfc(i,j,n) = 0.
                wxsg(i,j,n) = 0.
                wxaq(i,j,n) = 0.
                cemtrac(i,j,n) = 0.
                drydepv(i,j,n) =0.
              end do
            end do
          end do
          do j = 1 , jxp
            do i = 1 , iym1
              aertarf(j,i) = 0.
              aersrrf(j,i) = 0.
              aertalwrf(j,i) = 0.              
              aersrlwrf(j,i) = 0.
            end do
          end do

!


  end subroutine output_chem
!----------------------------------------------------------------
!================================================================

  subroutine outche2(idatex)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!     this subroutine writes the model chem                           c
!                                                                     c
!     iutl : is the output unit number for large-domain variables.    c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none
      type(rcm_time_and_date) , intent(in) :: idatex
!
      integer :: ni , itr , nj , nk , is , ie , js , je


!      character (len=50) :: subroutine_name='outche'
!      integer :: idindx=0
!
!      call time_begin(subroutine_name,idindx)
#ifdef BAND
      ni = iym2
      nj = jx
      nk = kz
      is = 2
      js = 1
      ie = iym1
      je = jx
#else
      ni = iym2
      nj = jxm2
      nk = kz
      is = 2
      js = 2
      ie = iym1
      je = jxm1
#endif

!!$
 call writerec_che2(nj, ni, je, ie, nk, ntr, chia_io,     &
                wdlsc_io, wdcvc_io, ddsfc_io, cemtrac_io,    &
                drydepv_io,  aerext_io, aerssa_io, aerasp_io,aertarf_io, aersrrf_io, &
                aertalwrf_io, aersrlwrf_io, cpsa_io, idatex)

      write (*,*) 'CHE variables written at ' , idatex 
      if (iaerosol > 0)  write (*,*) 'OPT variables written at ' , idatex 


      end subroutine outche2



!!$
!!$      subroutine outopt(idatex)
!!$
!!$!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!$!                                                                     c
!!$!     this subroutine writes the model chem (standard aerosol )       c
!!$!                                                                     c
!!$!     iutl : is the output unit number for large-domain variables.    c
!!$!                                                                     c
!!$!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!!$
!!$      implicit none
!!$!
!!$     type(rcm_time_and_date) , intent(in) :: idatex  
!!$      integer :: ni , itr , nj , nk , is , ie , js , je
!!$
!!$!      character (len=50) :: subroutine_name='outche'
!!$!      integer :: idindx=0
!!$!
!!$!      call time_begin(subroutine_name,idindx)
!!$#ifdef BAND
!!$      ni = iym2
!!$      nj = jx
!!$      nk = kz
!!$      itr = ntr
!!$      is = 2
!!$      js = 1
!!$      ie = iym1
!!$      je = jx
!!$#else
!!$      ni = iym2
!!$      nj = jxm2
!!$      nk = kz
!!$      itr = ntr
!!$      is = 2
!!$      js = 2
!!$      ie = iym1
!!$      je = jxm1
!!$#endif
!!$

!!$     call writerec_opt(nj, ni, je, ie, nk,       &
!!$                aerext_io, aerssa_io, aerasp_io,aertarf_io, aersrrf_io, &
!!$                aertalwrf_io, aersrlwrf_io, cpsa_io, idatex)
                               
!!$
!!$      write (*,*) 'OPT variables written at ' , idatex 
!!$
!!$      end subroutine outopt
 
     end module mod_che_output
