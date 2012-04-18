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
  subroutine chem_emission(lmonth)
    implicit none
    integer , intent(in) :: lmonth
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
    if ( lmonth == currmonth ) then 
      return 
    else
      currmonth = lmonth
    end if
    if ( myid == 0 ) then
      write(*,*)'READ CHEM EMISSION for month', lmonth
      chemsrc_io(:,:,:,:) = d_zero
      ! Also lmonth is not really necessary here, but KEEP THIS DIMENSION
      ! FOR HIGHER TEMPORAL RESOLUTION INVENTORIES 
      call read_emission(lmonth,chemsrc_io)
    end if
    call deco1_scatter(chemsrc_io,chemsrc,jcross1, &
                       jcross2,icross1,icross2,1,mpy,1,ntr)
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
subroutine emis_tend(ktau,j,lmonth)
#endif

    implicit none

    integer , intent(in) :: j , lmonth
    integer(8) , intent(in) :: ktau
!    real(dp) , intent(in) , dimension(iy) :: coszrs

!    real(dp) , intent(in) , dimension(iy,jxp) :: xlat
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
    ! real(dp) :: daylen , fact , maxelev , dayhr , isosrc , amp
    
    ! calculate the tendency linked to emissions from emission fluxes
    ! In the future split these calculations in corresponding module  ??

#ifdef CLM
#if (defined VOC)
    ! jj = j+(jxp*myid)
#endif
#endif

    ! 1 General case. In the future: add injection heights.
    do itr = 1 , ntr
      do i = ici1 , ici2
        if ( chtrname(itr)(1:4).ne.'DUST' .or. &
             chtrname(itr)(1:4).ne.'SSALT' ) then 
!!$          daylen = d_two*acos(-tan(declin)*tan(xlat(i,j)*degrad))*raddeg
!!$          daylen = daylen*24.0D0/360.0D0
!!$          ! Maximum sun elevation
!!$          maxelev = halfpi - ((xlat(i,j)*degrad)-declin)
!!$          fact = (halfpi-acos(coszrs(i)))/(d_two*maxelev)
!!$          amp = 12.0D0*mathpi/daylen
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
           if ( bvoc_trmask(itr) /= 0 ) then
             if ( ktau == 0 ) c2r_voc(jj,i,bvoc_trmask(itr)) = d_zero
             chemsrc(j,i,lmonth,itr) = c2r_voc(jj,i,bvoc_trmask(itr))/d_1
           end if
#endif
#endif
           ! update emission tendency according to chemsrc value
!!$           if ( chtrname(itr) == 'ISOP' ) then
!!$             chiten(i,kz,j,itr) = chiten(i,kz,j,itr) + &
!!$                           (amp)*chemsrc(j,i,lmonth,itr) * &
!!$                          sin(mathpi*fact)*egrav/(cdsigma(kz)*1.0D3)
!!$             ! diagnostic for source, cumul
!!$             cemtr(i,j,itr) = cemtr(i,j,itr) + (amp)*chemsrc(j,i,lmonth,itr) * &
!!$                           sin(mathpi*fact)*dtche/d_two
!!$           else
 
  
         if ( ichdrdepo == 1) then  
             chiten(j,i,kz,itr) = chiten(j,i,kz,itr) + &
                           chemsrc(j,i,lmonth,itr)*egrav/(cdsigma(kz)*1.0D3)
         elseif ( ichdrdepo ==2) then
!then emission is injected in the PBL scheme
              cchifxuw(j,i,itr) = cchifxuw(j,i,itr) +  chemsrc(j,i,lmonth,itr)
         end if


             ! diagnostic for source, cumul
             cemtr(i,j,itr) = cemtr(i,j,itr) + &
                              chemsrc(j,i,lmonth,itr)*dtche/d_two
!!$           end if
         end if
       end do
     end do
   end subroutine emis_tend
!
end module mod_che_emission
