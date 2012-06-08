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

module mod_che_emission
  !
  ! Chemical and aerosol surface emission 
  !
  use netcdf
  use mod_constants
  use mod_mpmessage
  use mod_service 
  use mod_dynparam
  use mod_che_common
  use mod_che_param
  use mod_che_mppio
  use mod_che_dust
  use mod_che_ncio
  use mod_che_indices
  !
  private
  !
  public :: chem_emission , emis_tend 
  !
contains
  !
  ! SURFACE EMIOSSION flux
  !
  subroutine chem_emission(lyear,lmonth,lday,lhour)
    implicit none
    integer , intent(in) :: lyear , lmonth , lday , lhour
    !
    integer , save :: currmonth
    character (len=64) :: subroutine_name = 'chem_emission'
    integer :: idindx = 0
    !
    call time_begin(subroutine_name,idindx)
    !
    ! read the monthly aerosol emission files
    ! FAB: remember for now, we have 1 emission file containing all monthly
    ! emission for the whole simulation period
    ! change that in the future. Also lmonth is not really necessary here,
    ! but KEEP THIS DIMENSION FOR HIGHER TEMPORAL RESOLUTION INVENTORIES 
    !
    if (chemsimtype(1:4)=='DUST') return
    if ( lmonth == currmonth ) then 
       return 
    else
       currmonth = lmonth
    end if
    if ( myid == 0 ) then
       write(*,*)'READ CHEM EMISSION for ', &
            lyear*1000000+lmonth*10000+lday*100+lhour
       ! Also lmonth is not really necessary here, but KEEP THIS DIMENSION
       ! FOR HIGHER TEMPORAL RESOLUTION INVENTORIES
       call read_emission(lyear,lmonth,chemsrc_io)
    end if
    call deco1_scatter(chemsrc_io,chemsrc,jcross1, &
         jcross2,icross1,icross2,1,ntr)
    call time_end(subroutine_name,idindx)
  end subroutine chem_emission
  !
  !
  ! Calculation of emission tendency
  !
#ifdef CLM
#if (defined VOC)
  subroutine emis_tend(ktau,j,lmonth,c2r_voc,bvoc_trmask)
#else
    subroutine emis_tend(ktau,j,lmonth)
#endif  
#else
      !FAB: no dirunal evol for now:  subroutine emis_tend(ktau,j,lmonth,xlat,coszrs,declin,dsigma)
      subroutine emis_tend(ktau,j,lmonth,declin)
#endif

        implicit none

        integer , intent(in) :: j , lmonth
        integer(8) , intent(in) :: ktau
        !    real(dp) , intent(in) , dimension(iy) :: coszrs

        real(dp) , intent(in) ::declin 
#ifdef CLM
#if (defined VOC)
        real(dp) , intent(in) , dimension(:,:,:) :: c2r_voc
        integer , intent(in) , dimension(:) :: bvoc_trmask
#endif
#endif
        integer :: i , itr

        ! real(dp) , intent(in) :: declin
        ! integer :: jj ! Full grid j-component
        ! integer :: ib , k
        real(dp) :: daylen , fact , maxelev , dayhr , isosrc , amp

        ! calculate the tendency linked to emissions from emission fluxes
        ! In the future split these calculations in corresponding module  ??

#ifdef CLM
#if (defined VOC)
        ! jj = j+(jxp*myid)
#endif
#endif

        ! Modify chemsrc for species that need a dirunal cycle

        if (iisop > 0) then 


           do i = ici1 , ici2

           daylen = d_two*acos(-tan(declin)*tan(cxlat(j,i)*degrad))*raddeg
           daylen = daylen*24.0D0/360.0D0
!!$          ! Maximum sun elevation
           maxelev = halfpi - ((cxlat(j,i)*degrad)-declin)
           fact = (halfpi-acos(czen(j,i)))/(d_two*maxelev)
           amp = 12.0D0*mathpi/daylen

           tmpsrc(j,i,iisop) =   chemsrc(j,i,iisop)
           chemsrc(j,i,iisop) =  (amp)*chemsrc(j,i,iisop) * &
                sin(mathpi*fact)*egrav/(cdsigma(kz)*1.0D3)

           end do
         
        end if


        !overwrite chemsrc for biogenic in case of CLM/BVOC option  
#ifdef CLM
#if (defined VOC)
        ! test if CLM BVOC is activated and overwrite chemsrc.
        ! Below included in order to include CLM-MEGAN biogenic emission
        ! into the gas phase chemistry scheme
        ! bvoc_trmask is used to identify which tracer is also a biogenic
        ! emission species contained in MEGAN.  If it is included then
        ! then include MEGAN emissions into the source value
        ! NOTE:  ibvoc=1 means used MEGAN emissions.  ibvoc is forced to
        ! zero when using BATS

        ! jj = j+(jxp*myid) FAB A REVOIRR!
        print*, 'bio emission CLM/BVOC option to be fixed'
        stop 
        do itr =1,ntr
           if ( bvoc_trmask(itr) /= 0 ) then
              if ( ktau == 0 ) c2r_voc(jj,:,bvoc_trmask(itr)) = d_zero
              chemsrc(j,:,itr) = c2r_voc(jj,:,bvoc_trmask(itr))/d_1
           end if
        end do
#endif
#endif


        ! add the source term to tracer tendency

        if ( ichdrdepo /= 2 ) then  

           do itr = 1 , ntr
              do i = ici1 , ici2
                 if ( chtrname(itr)(1:4).ne.'DUST' .or. &
                      chtrname(itr)(1:4).ne.'SSALT' ) then 

                    chiten(j,i,kz,itr) = chiten(j,i,kz,itr) + &
                         chemsrc(j,i,itr)*egrav/(cdsigma(kz)*1.0D3)

                    ! diagnostic for source, cumul
                    cemtr(j,i,itr) = cemtr(j,i,itr) + &
                         chemsrc(j,i,itr)*dtche/d_two
                 end if
              end do
           end do

        elseif ( ichdrdepo ==2) then

           do itr = 1 , ntr
              do i = ici1 , ici2
                 if ( chtrname(itr)(1:4).ne.'DUST' .or. &
                      chtrname(itr)(1:4).ne.'SSALT' ) then 

                    !then emission is injected in the PBL scheme
                    cchifxuw(j,i,itr) = cchifxuw(j,i,itr) +  chemsrc(j,i,itr)
                    ! diagnostic for source, cumul
                    cemtr(j,i,itr) = cemtr(j,i,itr) + &
                         chemsrc(j,i,itr)*dtche/d_two
                  end if 
                 end do
              end do
           end if

! put back isop source to its nominal value 
          
         if (iisop > 0) then 
           chemsrc(j,:,iisop) = tmpsrc(j,:,iisop)
         end if 


         end subroutine emis_tend
         !
       end module mod_che_emission
