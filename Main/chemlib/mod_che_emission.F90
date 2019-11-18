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
  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_mpmessage
  use mod_service
  use mod_dynparam
#ifdef CLM45
  use mod_runparams , only : rcm_megan_enabled
#endif
  use mod_che_common
  use mod_che_param
  use mod_che_mppio
  use mod_che_dust
  use mod_che_ncio
  use mod_che_indices
  use netcdf
  !
  implicit none

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
    integer(ik4) , intent(in) :: lyear , lmonth , lday , lhour
    integer(ik4) , save :: curry = -1
    integer(ik4) , save :: currm = -1
    integer(ik4) , save :: currd = -1
    integer(ik4) , save :: ifreq = -1
    integer(ik4) , save :: currybb = -1
    integer(ik4) , save :: currmbb = -1
    integer(ik4) , save :: currdbb = -1
    integer(ik4) , save :: ifreqbb = -1

#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'chem_emission'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    !
    ! read the aerosol emission files
    !
    if ( chemsimtype(1:2) == 'DU' .or. &
         chemsimtype(1:4) == 'SSLT' .or. &
         chemsimtype(1:4) == 'DUSS' .or. &
         chemsimtype(1:4) == 'MINE' ) then
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
      return
    end if

    if ( ifreq == ifrqmon ) then
      if ( curry == lyear .and. currm == lmonth ) then
        if ( myid == italk ) then
          write(stdout,*) &
            'EMISSION for  ',lyear*1000000+lmonth*10000+100*lday,' ready', &
            ' from ',curry*1000000+currm*10000+100*currd
        end if
#ifdef DEBUG
        call time_end(subroutine_name,idindx)
#endif
!       return
! FAB use eceptionally a goto here to allow reading of BB emission
        goto 99
      end if
    else if ( ifreq == ifrqday ) then
      if ( curry == lyear .and. currm == lmonth .and. currd == lday ) goto 99
    end if
    curry = lyear
    currm = lmonth
    currd = lday
    if ( myid == italk ) then
      write(stdout,*) 'READ CHEM EMISSION for ', &
               lyear*1000000+lmonth*10000+lday
    end if
    ! Also lmonth is not really necessary here, but KEEP THIS DIMENSION
    ! FOR HIGHER TEMPORAL RESOLUTION INVENTORIES
    call read_emission(ifreq,lyear,lmonth,lday,lhour,chemsrcan)

    ! allow to increase OC emission to account for SOA via namelist
    if ( iochb > 0 ) then
      chemsrcan(:,:,iochb) = rocemfac * chemsrcan(:,:,iochb)
    end if

    ! Handle biomass burning emissions , possibly at a different frequency
    !
    99 continue

    if ( ismoke == 1 ) then
      if ( ifreqbb == ifrqmon ) then
        if ( currybb == lyear .and. currmbb == lmonth ) then
          if ( myid == italk ) then
            write(stdout,*) &
              'BB EMISSION for  ',lyear*1000000+lmonth*10000+100*lday,&
              ' ready from ',currybb*1000000+currmbb*10000+100*currdbb
          end if
#ifdef DEBUG
          call time_end(subroutine_name,idindx)
#endif
          return
        end if
      else if ( ifreqbb == ifrqday ) then
        if ( currybb == lyear .and. &
             currmbb == lmonth .and. &
             currdbb == lday ) then
#ifdef DEBUG
          call time_end(subroutine_name,idindx)
#endif
          return
        end if
      end if
      currybb = lyear
      currmbb = lmonth
      currdbb = lday
      if ( myid == italk ) then
        write(stdout,*) 'READ CHEM BIO. BURN. EMISSION for ', &
                 lyear*1000000+lmonth*10000+lday
      end if
      call read_bioburn_emission(ifreqbb,lyear,lmonth,lday,lhour,chemsrcbb)
      ! define the smoke tracer emission
      if ( ism1 > 0 ) then
        chemsrcbb(:,:,ism1) = chemsrcbb(:,:,ibchb) + &
                              rocemfac * chemsrcbb(:,:,iochb)
        chemsrcbb(:,:,iochb) = d_zero
        chemsrcbb(:,:,ibchb) = d_zero
      end if
      ! finally add biomass burning emissions to anthropogenic emissions
      chemsrc(:,:,:) = chemsrcbb(:,:,:) + chemsrcan(:,:,:)
    else
      chemsrc(:,:,:) = chemsrcan(:,:,:)
    end if

#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine chem_emission
  !
  ! Calculation of emission tendency
  !
  subroutine emis_tend(i,declin)
    implicit none
    integer(ik4) , intent(in) :: i
    real(rk8) , intent(in) :: declin

    integer(ik4)  :: j , itr
    real(rkx) :: daylen , fact , maxelev , amp , dayhr
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'emis_tend'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    ! gather anthropogenic and biogenic emissions and
    ! calculate the emission tendency from emission fluxes

#ifdef CLM45
    ! overwrite or add chemsrc for biogenic in case of CLM45/BVOC option
    ! test if CLM BVOC is activated and add bio component to chemsrc.
    ! isoprene ( overwrite biogenic emission from inventory  / rq
    ! this way  we  negelect the anthropog isoprene emissions )
    if ( iisop > 0 ) then
      if ( rcm_megan_enabled ) then
        do j = jci1, jci2
          chemsrc(j,i,iisop) = cvoc_em_clm(j,i,iisop)
        end do
      else
        do j = jci1 , jci2
          dayhr = -tan(declin)*tan(cxlat(j,i)*degrad)
          if ( dayhr < -1 .or. dayhr > 1 ) then
            tmpsrc(j,i,iisop)  = chemsrc(j,i,iisop)
            chemsrc(j,i,iisop) = 0.0
          else
            daylen = d_two*acos(-tan(declin)*tan(cxlat(j,i)*degrad))*raddeg
            daylen = d_two*acos(dayhr)*raddeg
            daylen = daylen*24.0_rkx/360.0_rkx
            ! Maximum sun elevation
            maxelev = halfpi - ((cxlat(j,i)*degrad)-declin)
            fact = (halfpi-acos(czen(j,i)))/(d_two*maxelev)
            amp = 12.0_rkx*mathpi/daylen
            tmpsrc(j,i,iisop)  = chemsrc(j,i,iisop)
            chemsrc(j,i,iisop) = amp*chemsrc(j,i,iisop) * sin(mathpi*fact)
          end if
        end do
      end if
    end if

    ! other BVOC : to do !!

    ! methane : add the biogenic part from CLM to anthropogenic
    ! nb the CN/LCH4 CLM flag have to be activated
#ifdef LCH4
    if ( ich4 > 0 ) then
      do j = jci1, jci2
        chemsrc(j,i,ich4) = chemsrc(j,i,ich4) + cvoc_em_clm(j,i,ich4)
      end do
    end if
#endif
#else
    ! Modify chemsrc for isoprene from inventory to account for
    ! simplified dirunal cycle
    if ( iisop > 0 ) then
      do j = jci1 , jci2
        dayhr = -tan(declin)*tan(cxlat(j,i)*degrad)
        if ( dayhr < -1 .or. dayhr > 1 ) then
          tmpsrc(j,i,iisop)  = chemsrc(j,i,iisop)
          chemsrc(j,i,iisop) = 0.0
        else
          daylen = d_two*acos(-tan(declin)*tan(cxlat(j,i)*degrad))*raddeg
          daylen = d_two*acos(dayhr)*raddeg
          daylen = daylen*24.0_rkx/360.0_rkx
          ! Maximum sun elevation
          maxelev = halfpi - ((cxlat(j,i)*degrad)-declin)
          fact = (halfpi-acos(czen(j,i)))/(d_two*maxelev)
          amp = 12.0_rkx*mathpi/daylen
          tmpsrc(j,i,iisop)  = chemsrc(j,i,iisop)
          chemsrc(j,i,iisop) = amp*chemsrc(j,i,iisop) * sin(mathpi*fact)
        end if
      end do
    end if
#endif
    !
    ! add the source term to tracer tendency
    !
    if ( ichdrdepo /= 2 ) then
      do itr = 1 , ntr
        if ( chtrname(itr)(1:2) == 'DU' .or. &
             chtrname(itr)(1:4) == 'SSLT' .or. &
             chtrname(itr)(1:6) == 'POLLEN'.or. &
             any ( mine_name == chtrname(itr)(1:4) )) cycle
        if ( idynamic == 3 ) then
          do j = jci1 , jci2
            chiten(j,i,kz,itr) = chiten(j,i,kz,itr) + &
                  chemsrc(j,i,itr)/(cdzq(j,i,kz)*crhob3d(j,i,kz))
          end do
        else
          do j = jci1 , jci2
            chiten(j,i,kz,itr) = chiten(j,i,kz,itr) + &
                chemsrc(j,i,itr)/(cdzq(j,i,kz)*crhob3d(j,i,kz))*cpsb(j,i)
          end do
        end if
        do j = jci1 , jci2
          ! diagnostic for source, cumul
          cemtrac(j,i,itr) = cemtrac(j,i,itr) + chemsrc(j,i,itr)*cfdout
        end do
        if ( ichdiag > 0 ) then
          do j = jci1 , jci2
            cemisdiag(j,i,kz,itr) = cemisdiag(j,i,kz,itr) + &
                chemsrc(j,i,itr)/(cdzq(j,i,kz)*crhob3d(j,i,kz)) * cfdout
          end do
        end if
      end do
    else if ( ichdrdepo == 2 ) then
      do itr = 1 , ntr
        if ( chtrname(itr)(1:2) == 'DU' .or. &
             chtrname(itr)(1:4) == 'SSLT' .or. &
             chtrname(itr)(1:6) == 'POLLEN'.or. &
             any ( mine_name == chtrname(itr)(1:4) )) cycle
        if ( ibltyp /= 2 ) then
          ! if PBL scheme is not UW then calculate emission tendency
          if ( idynamic == 3 ) then
            do j = jci1 , jci2
              chiten(j,i,kz,itr) = chiten(j,i,kz,itr) + &
                   chemsrc(j,i,itr)/(cdzq(j,i,kz)*crhob3d(j,i,kz))
            end do
          else
            do j = jci1 , jci2
              chiten(j,i,kz,itr) = chiten(j,i,kz,itr) + &
                   chemsrc(j,i,itr)/(cdzq(j,i,kz)*crhob3d(j,i,kz))*cpsb(j,i)
            end do
          end if
        else
          ! otherwise emission is injected in the PBL scheme ( together
          ! with dry deposition) for tend calculation
          do j = jci1 , jci2
            chifxuw(j,i,itr) = chifxuw(j,i,itr) + chemsrc(j,i,itr)
          end do
        end if
        ! diagnostic for source, cumul
        do j = jci1 , jci2
          cemtrac(j,i,itr) = cemtrac(j,i,itr) + chemsrc(j,i,itr)*cfdout
        end do
        if ( ichdiag > 0 ) then
          if ( ibltyp /= 2 ) then
            do j = jci1 , jci2
            ! in this case we calculate emission tendency diagnostic, but
            ! this term will also be included in BL tendency diagnostic
            ! if UW scheme is used.
              cemisdiag(j,i,kz,itr) = cemisdiag(j,i,kz,itr) + &
                 chemsrc(j,i,itr)/ ( cdzq(j,i,kz)*crhob3d(j,i,kz)) * cfdout
            end do
          end if
        end if
      end do
    end if

    ! put back chemsrc to its mean value after diurnal cycle
#ifndef CLM45
    if ( iisop > 0 ) then
      chemsrc(:,i,iisop) = tmpsrc(:,i,iisop)
    end if
#else
    if ( .not. rcm_megan_enabled ) then
      if ( iisop > 0 ) then
        chemsrc(:,i,iisop) = tmpsrc(:,i,iisop)
      end if
    end if
#endif

#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine emis_tend
!
end module mod_che_emission
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
