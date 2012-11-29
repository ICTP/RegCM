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
!
module mod_che_ncio
!
  use mod_intkinds
  use mod_realkinds
  use mod_nchelper
  use mod_dynparam
  use mod_mpmessage
  use mod_che_indices
  use mod_che_common
  use mod_runparams
  use mod_domain
  use netcdf
!
  private
!
  public :: read_texture , read_emission , recc
  public :: prepare_chem_out , init_mod_che_ncio
  public :: open_chbc , close_chbc , chbc_search , read_chbc

  public :: chbc_ivar , n_chbcvar , n_aebcvar

  integer(ik4) :: istatus
  integer(ik4) :: recc

  integer(ik4) , parameter :: n_chevar = 20
  integer(ik4) , parameter :: n_oxbcvar = 5
  integer(ik4) , parameter :: n_optvar = 10
  integer(ik4) , parameter :: n_chbcvar = 25
  integer(ik4) :: n_aebcvar
  integer(ik4) :: ichin , iaein , ioxin

  character(len=8) , dimension(n_chbcvar) :: chbcname
  character(len=8) , dimension(n_oxbcvar) :: oxbcname
  character(len=8) , target , dimension(4) :: aedust
  character(len=8) , target , dimension(2) :: aesslt
  character(len=8) , target , dimension(4) :: aecarb
  character(len=8) , target , dimension(2) :: aesulf
  character(len=8) , target , dimension(6) :: aesuca
  character(len=8) , target , dimension(12) :: aeaero

  character(len=8) , pointer , dimension(:) :: aebcname
  integer(ik4) , dimension(n_chbcvar) :: chbc_ivar
  integer(ik4) , dimension(n_oxbcvar) :: oxbc_ivar
  integer(ik4) , dimension(:) , pointer :: aebc_ivar
  
  type(rcm_time_and_date) , dimension(:) , allocatable :: chbc_idate
  integer(ik4) :: ibcrec , ibcnrec

  integer(ik4) , public , parameter :: ifrqmon = 1
  integer(ik4) , public , parameter :: ifrqday = 2
  integer(ik4) , public , parameter :: ifrqhrs = 3

  character(256) :: dname , icbcname

  data ichin   /-1/
  data iaein   /-1/
  data ibcrec  / 1/
  data ibcnrec / 0/
  data ioxin  /-1/

  data chbcname /'O3      ','NO      ','NO2     ','HNO3    ', &
                 'N2O5    ','H2O2    ','CH4     ','CO      ', &
                 'CH2O    ','CH3OH   ','C2H5OH  ','C2H4    ', &
                 'C2H6    ','CH3CHO  ','CH3COCH3','BIGENE  ', &
                 'BIGALK  ','C3H6    ','C3H8    ','ISOP    ', &
                 'TOLUENE ','PAN     ','SO2     ','SO4     ', &
                 'DMS     '/
  data oxbcname /'OH      ','HO2     ','O3      ', 'NO3    ','H2O2   ' /
  data aedust / 'DUST01' , 'DUST02' , 'DUST03', 'DUST04' /
  data aesslt / 'SSLT01' , 'SSLT02' /
  data aecarb / 'BC_HL' , 'BC_HB' , 'OC_HL' , 'OC_HB' /
  data aesulf / 'SO2' , 'SO4' /
  data aesuca / 'BC_HL' , 'BC_HB' , 'OC_HL' , 'OC_HB' , 'SO2' , 'SO4' /
  data aeaero / 'BC_HL' , 'BC_HB' , 'OC_HL' , 'OC_HB' , 'SO2' , 'SO4' , &
                'SSLT01' , 'SSLT02', 'DUST01', 'DUST02', 'DUST03' , &
                'DUST04' /

  contains

    subroutine init_mod_che_ncio(chemsymtype)
      implicit none
      character(len=8) , intent(in) :: chemsymtype

      dname = trim(dirter)//pthsep//trim(domname)//'_DOMAIN000.nc'
      n_aebcvar = 0
      select case ( chemsymtype )
        case ( 'DUST' )
          n_aebcvar = 4
          aebcname => aedust
        case ( 'SSLT' )
          n_aebcvar = 2
          aebcname => aesslt
        case ( 'CARB' )
          n_aebcvar = 4
          aebcname => aecarb
        case ( 'SULF' )
          n_aebcvar = 2
          aebcname => aesulf
        case ( 'SUCA' )
          n_aebcvar = 6
          aebcname => aesuca
        case ( 'AERO' )
          n_aebcvar = 12
          aebcname => aeaero
      end select
      if ( n_aebcvar > 0 ) then
        call getmem1d(aebc_ivar,1,n_aebcvar,'ncio:aebc_ivar')
      end if
    end subroutine init_mod_che_ncio

    subroutine read_texture(nats,texture)
      implicit none
      integer(ik4) , intent(in) :: nats
      real(rk8) , pointer , dimension(:,:,:) , intent(out) :: texture

      integer(ik4) :: ivarid
      integer(ik4) :: i , j , n
      integer(ik4) :: idmin
      integer(ik4) , dimension(3) :: istart , icount
      character(256) :: dname
      real(rk4), dimension(jx,iy) ::  toto

      dname = trim(dirter)//pthsep//trim(domname)//'_DOMAIN000.nc'

      istatus = nf90_open(dname, nf90_nowrite, idmin)
      call check_ok(__FILE__,__LINE__, &
                  'Error Opening Domain file '//trim(dname),'DOMAIN FILE OPEN')
      istatus = nf90_inq_varid(idmin, 'texture_fraction', ivarid)
      call check_ok(__FILE__,__LINE__,'Variable texture_fraction miss', &
                  'DOMAIN FILE')
      istart(2) = 1
      istart(1) = 1
      icount(3) = 1
      icount(2) = iy
      icount(1) = jx
      do n = 1 , nats
        istart(3) = n
        istatus = nf90_get_var(idmin, ivarid, toto, istart, icount)
        call check_ok(__FILE__,__LINE__,'Variable texture_frac read error', &
                      'DOMAIN FILE')
        do i = 1 , iy
          do j = 1 , jx
            texture(j,i,n) = dble(toto(j,i))*0.01D0
            if ( texture(j,i,n) < d_zero ) texture(j,i,n) = d_zero
          end do
        end do
      end do
      istatus = nf90_close(idmin)
      call check_ok(__FILE__,__LINE__,'Domain file close error','DOMAIN FILE')
    end subroutine read_texture

    subroutine read_emission(ifreq,lyear,lmonth,lday,lhour,echemsrc)
      implicit none
      integer(ik4) , intent(in) :: lyear , lmonth , lday , lhour
      integer(ik4) , intent(out) :: ifreq
      real(rk8) , pointer , dimension(:,:,:) , intent(out) :: echemsrc
      character(256) :: aername
      integer(ik4) :: n,ncid , itvar, idimid, chmnrec,sdim
      character(64) ::chemi_timeunits , chemi_timecal
      real(rk8) , dimension(:) , allocatable :: emtimeval
      integer(ik4) , dimension(4) :: istart , icount
      integer(ik4) :: year , month , day , hour 
      type(rcm_time_and_date) :: tchdate

      ! FAB: remember for now, we have 1 emission file containing all monthly
      ! emission for the whole simulation period 
      ! change that in the future.,

      aername = trim(dirglob)//pthsep//trim(domname)//'_CHEMISS.nc'
      print *, 'Opening ch. emission file ', trim(aername)

      istatus = nf90_open(aername, nf90_nowrite, ncid)
      call check_ok(__FILE__,__LINE__, &
                    'Error Opening chem emissiom file '//trim(aername), &
                    'CHE EMISS FILE OPEN ERROR')

      istatus = nf90_inq_dimid(ncid, 'time', idimid)
      call check_ok(__FILE__,__LINE__,'Dimension time miss', 'CHEMI FILE')
      istatus = nf90_inquire_dimension(ncid, idimid, len=chmnrec)
      call check_ok(__FILE__,__LINE__,'Dimension time read error', 'CHEMI FILE')

      allocate (emtimeval(chmnrec))

      istatus = nf90_inq_varid(ncid, 'time', itvar)
      call check_ok(__FILE__,__LINE__,'variable time miss', 'CHEMISS FILE')
      istatus = nf90_get_att(ncid, itvar, 'units', chemi_timeunits)
      call check_ok(__FILE__,__LINE__,'variable time units miss', &
                    'CHEMISS FILE')
      if ( chemi_timeunits(1:6) == 'months' ) then
        ifreq = ifrqmon
      else if ( chemi_timeunits(1:4) == 'days' ) then
        ifreq = ifrqday
      else if ( chemi_timeunits(1:5) == 'hours' ) then
        ifreq = ifrqhrs
      else
        call fatal(__FILE__,__LINE__,'NO CODED FREQUENCY IN CHEMISS FILE')
      end if
        
      istatus = nf90_get_att(ncid, itvar, 'calendar', chemi_timecal)
      if ( istatus /= nf90_noerr ) then
        chemi_timecal = 'gregorian'
      end if
      istatus = nf90_get_var(ncid, itvar, emtimeval)
      call check_ok(__FILE__,__LINE__,'variable time read error', 'ICBC FILE')
    
      recc = 0
      looprec: &
      do n = 1 , chmnrec
        tchdate = timeval2date(emtimeval(n),chemi_timeunits,chemi_timecal)
        call split_idate(tchdate,year,month,day,hour)
        select case (ifreq)
          case (ifrqmon)
            if ( year == lyear .and. month == lmonth ) then
              recc = n
              exit looprec
            end if
          case (ifrqday)
            if ( year == lyear .and. month == lmonth .and. day == lday ) then
              recc = n
              exit looprec
            end if
          case (ifrqhrs)
            if ( year == lyear .and. month == lmonth .and. &
                  day == lday  .and. hour  == lhour ) then
              recc = n
              exit looprec
            end if
        end select
      end do looprec
       
      if ( recc == 0 ) then
        print *,'chem emission : time record not found emission file, STOP ! '
        call fatal(__FILE__,__LINE__,'IO ERROR in CHEM EMISSION')
      end if  

      !*** intialized in start_chem
      !*** Advice record counter
      istart = 0
      icount = 0 
      istatus = nf90_inq_dimid(ncid, 'lev', idimid)

      if(istatus /= nf90_noerr) then
        ! no lev diemsion in emission variables
        istart(1) = 1
        istart(2) = 1
        istart(3) = recc
        icount(1) = jx
        icount(2) = iy
        icount(3) = 1
        sdim = 3
      else
        istart(1) = 1
        istart(2) = 1
        istart(3) = 1
        istart(4) = recc
        icount(1) = jx
        icount(2) = iy
        icount(3) = 1
        icount(4) = 1
        sdim=4
      end if
      ! CO emission                  
      if ( ico /= 0 ) then
        call rvar(ncid,istart,icount,ico,echemsrc,'CO_flux',.false.,sdim)
        print*, 'FAB emis testco','ico', maxval(echemsrc)
      end if

      ! NO emission                  
      if ( ino /= 0 ) then
        call rvar(ncid,istart,icount,ino,echemsrc, &
                  'NO_flux',.false.,sdim)
      end if

      ! HCHO emission                  
      if ( ihcho /= 0 ) then
         call rvar(ncid,istart,icount,ihcho,echemsrc, &
           'HCHO_flux',.false.,sdim)
      end if
      ! ACET emission                  
      if ( iacet /= 0 ) then
        call rvar(ncid,istart,icount,iacet,echemsrc, &
            'ACET_flux',.false.,sdim)
      end if
      ! SO2 emission
      if ( iso2 /= 0 ) then
        call rvar(ncid,istart,icount,iso2,echemsrc, &
                  'SO2_flux',.false.,sdim)
      end if
      !NH3
      if ( iNH3 /= 0 ) then
        call rvar(ncid,istart,icount,inh3,echemsrc, &
                 'NH3_flux',.false.,sdim)
      end if
      ! CH4
      if ( ich4 /= 0 ) then
        call rvar(ncid,istart,icount,ich4,echemsrc, &
                  'CH4_flux',.false.,sdim)
      end if
      ! Ethane
      if ( ic2h6 /= 0 ) then
        call rvar(ncid,istart,icount,ic2h6,echemsrc, &
                  'C2H6_flux',.false.,sdim)
      end if
      ! PAR
      if ( ipar /= 0 ) then
       call rvar(ncid,istart,icount,ipar,echemsrc, &
               'PAR_flux',.false.,sdim)
      end if
      ! Ethene
      if ( iethe /= 0 ) then
        call rvar(ncid,istart,icount,iethe,echemsrc, &
                  'ETHE_flux',.false.,sdim)
      end if
      ! Termenal Alkene
      if ( iolt /= 0 ) then
        call rvar(ncid,istart,icount,iolt,echemsrc, &
                  'OLT_flux',.false.,sdim)
      end if

      ! Internal Alkene
!     if ( ioli /= 0 ) then
!       call rvar(ncid,istart,icount,ioli,echemsrc,'OLI_flux',.false.)
!     end if

      ! Isoprene
      if ( iisop /= 0 ) then
        call rvar(ncid,istart,icount,iisop,echemsrc, &
                  'ISOP_BIO_flux',.false.,sdim)
        ! here use io3(never emuitted) to temporarily read anthropo
        ! isoprene and add to biogenic. Should be refined 
        call rvar(ncid,istart,icount,io3,echemsrc,'ISO_flux',.false.,sdim)
        echemsrc(:,:,iisop) =  echemsrc(:,:,iisop) + echemsrc(:,:,io3)
        echemsrc(:,:,io3) = d_zero
      end if

      ! Toluene
      if ( itolue /= 0 ) then
        call rvar(ncid,istart,icount,itolue,echemsrc, &
                  'TOL_flux',.false.,sdim)
      end if

      ! Xylene
      if ( ixyl /= 0 ) then
        call rvar(ncid,istart,icount,ixyl,echemsrc, &
                 'XYL_flux',.false.,sdim)
      end if
      ! Acetaldehyde
      if ( iald2 /= 0 ) then
        call rvar(ncid,istart,icount,iald2,echemsrc,'ALD2_flux',.false.,sdim)
      end if
      ! Methanol + Ethanol
      if ( imoh /= 0 ) then
        call rvar(ncid,istart,icount,imoh,echemsrc, &
                 'MOH_flux',.false.,sdim)
      end if           
      !acids
      if ( ircooh /= 0 ) then
        call rvar(ncid,istart,icount,ircooh,echemsrc, &
                  'RCOOH_flux',.false.,sdim)
      end if

!!$   ! DMS
!!$   if ( idms /= 0 ) then
!!$     ! call rvar(ncid,istart,icount,idms,echemsrc,'o_DMS',.false.)
!!$   end if

      ! OC and BC anthropogenic + biomass burning
      if ( ibchb /= 0 ) then
        call rvar(ncid,istart,icount,ibchb,echemsrc, &
                  'BC_flux',.false.,sdim)
      end if
      if ( iochb /= 0 ) then
        call rvar(ncid,istart,icount,iochb,echemsrc, &
                  'OC_flux',.false.,sdim)
      end if

      if (ipollen /=0 ) then 
        call  rvar(ncid,istart,icount,ipollen,echemsrc, &
                  'POLLEN',.false.,sdim)
      end if

      where (echemsrc(:,:,:) < d_zero ) echemsrc(:,:,:) = d_zero

      istatus = nf90_close(ncid)
      call check_ok(__FILE__,__LINE__, &
                    'Error Closing Chem emission file '//trim(aername), &
                    'CH EMISS FILE CLOSE ERROR')
      deallocate (emtimeval)
    end subroutine read_emission

    subroutine rvar(ncid,istart,icount,ind,echemsrc,cna,lh,sdim,cnb,cnc,cnd)
      implicit none

      integer(ik4) , intent(in) :: ncid,sdim
      integer(ik4) , dimension(4) , intent(in) :: istart , icount
      real(rk8) , pointer , dimension(:,:,:) , intent(out) :: echemsrc
      logical , intent(in) :: lh
      character(len=*) , intent(in) :: cna
      character(len=*) , intent(in) , optional :: cnb
      character(len=*) , intent(in) , optional :: cnc
      character(len=*) , intent(in) , optional :: cnd
      integer(ik4) :: ivarid 
      real(rk4) , dimension(jx,iy) :: toto
      integer(ik4) :: i , j , ind

      istatus = nf90_inq_varid(ncid, cna, ivarid)
      call check_ok(__FILE__,__LINE__, &
                    'Variable '//cna//' miss','CHEM_EMISS FILE')
      istatus = nf90_get_var(ncid,ivarid,toto,istart(1:sdim),icount(1:sdim))
      call check_ok(__FILE__,__LINE__, &
                    'Variable '//cna//' read err','CHEM_EMISS FILE')
      if ( lh ) then  ! half of lumped Aromatics
        do i = 1 , iy
          do j = 1 , jx
            echemsrc(j,i,ind) = d_half*toto(j,i)
          end do
        end do
      else
        do i = 1 , iy
          do j = 1 , jx
            echemsrc(j,i,ind) = toto(j,i)
          end do
        end do
      end if
      if ( present(cnb) ) then
        istatus = nf90_inq_varid(ncid, cnb, ivarid)
        call check_ok(__FILE__,__LINE__, &
                      'Variable '//cnb//' miss','CHEM_EMISS FILE')
        istatus = nf90_get_var(ncid,ivarid,toto,istart(1:sdim),icount(1:sdim))
        call check_ok(__FILE__,__LINE__, &
                      'Variable '//cnb//' read err','CHEM_EMISS FILE')
        do i = 1 , iy
          do j = 1 , jx
            echemsrc(j,i,ind) = toto(j,i) + echemsrc(j,i,ind)
          end do
        end do
      end if
      if ( present(cnc) ) then
        istatus = nf90_inq_varid(ncid, cnc, ivarid)
        call check_ok(__FILE__,__LINE__, &
                      'Variable '//cnc//' miss','CHEM_EMISS FILE')
        istatus = nf90_get_var(ncid,ivarid,toto,istart(1:sdim),icount(1:sdim))
        call check_ok(__FILE__,__LINE__, &
                      'Variable '//cnc//' read err','CHEM_EMISS FILE')
        do i = 1 , iy
          do j = 1 , jx
            echemsrc(j,i,ind) = toto(j,i) + echemsrc(j,i,ind)
          end do
        end do
      end if
      if ( present(cnd) ) then
        istatus = nf90_inq_varid(ncid, cnd, ivarid)
        call check_ok(__FILE__,__LINE__, &
                      'Variable '//cnd//' miss','CHEM_EMISS FILE')
        istatus = nf90_get_var(ncid,ivarid,toto,istart(1:sdim),icount(1:sdim))
        call check_ok(__FILE__,__LINE__, &
                      'Variable '//cnd//' read err','CHEM_EMISS FILE')
        do i = 1 , iy
          do j = 1 , jx
            echemsrc(j,i,ind) = toto(j,i) + echemsrc(j,i,ind)
          end do
        end do
      end if
    end subroutine rvar

!============================================================================

    integer(ik4) function chbc_search(idate)
      implicit none
      type(rcm_time_and_date) , intent(in) :: idate
      type(rcm_time_interval) :: tdif
      character(len=32) :: appdat1, appdat2
      if (idate > chbc_idate(ibcnrec) .or. idate < chbc_idate(1)) then
        chbc_search = -1
      else
        tdif = idate - chbc_idate(1)
        ibcrec = (idnint(tohours(tdif))/ibdyfrq)+1 
        if ( ibcrec < 1 .or. ibcrec > ibcnrec ) then
          appdat1 = tochar(idate)
          write (6,*) 'Record is not found in CHBC file for ',appdat1
          appdat1 = tochar(chbc_idate(1))
          appdat2 = tochar(chbc_idate(ibcnrec))
          write (6,*) 'Range is : ', appdat1, '-', appdat2
          call fatal(__FILE__,__LINE__,'CHBC READ')
        end if
        chbc_search = ibcrec
      end if 
    end function chbc_search

!============================================

    subroutine open_chbc(idate)
      implicit none
      type(rcm_time_and_date) , intent(in) :: idate
      character(10) :: ctime
      integer(ik4) :: ibcid , idimid , itvar , i , chkdiff
      real(rk8) , dimension(:) , allocatable :: icbc_nctime
      character(64) :: icbc_timeunits , icbc_timecal

      call close_chbc
      write (ctime, '(i10)') toint10(idate)
      if ( igaschem == 1 ) then
        icbcname = trim(dirglob)//pthsep//trim(domname)//'_CHBC.'//ctime//'.nc'
        istatus = nf90_open(icbcname, nf90_nowrite, ichin)
        call check_ok(__FILE__,__LINE__, &
              'Error Opening ICBC file '//trim(icbcname),'CHBC FILE OPEN')
        call check_dims(ichin)
        ibcid = ichin
      end if
      if ( iaerosol == 1 ) then
        icbcname = trim(dirglob)//pthsep//trim(domname)//'_AEBC.'//ctime//'.nc'
        istatus = nf90_open(icbcname, nf90_nowrite, iaein)
        call check_ok(__FILE__,__LINE__, &
              'Error Opening ICBC file '//trim(icbcname),'AEBC FILE OPEN')
        call check_dims(iaein)
        ibcid = iaein
      end if
      if ( ioxclim == 1 ) then
        icbcname = trim(dirglob)//pthsep//trim(domname)//'_OXCL.'//ctime//'.nc'
        istatus = nf90_open(icbcname, nf90_nowrite, ioxin)
        call check_ok(__FILE__,__LINE__, &
              'Error Opening ICBC file '//trim(icbcname),'OXBC FILE OPEN')
        call check_dims(iaein)
        ibcid = ioxin
      end if
      ibcrec = 1
      ibcnrec = 0
      istatus = nf90_inq_dimid(ibcid, 'time', idimid)
      call check_ok(__FILE__,__LINE__,'Dimension time miss', 'ICBC FILE')
      istatus = nf90_inquire_dimension(ibcid, idimid, len=ibcnrec)
      call check_ok(__FILE__,__LINE__,'Dimension time read error', 'ICBC FILE')
      if ( ibcnrec < 1 ) then
        write (6,*) 'Time var in ICBC has zero dim.'
        call fatal(__FILE__,__LINE__,'ICBC READ')
      end if
      istatus = nf90_inq_varid(ibcid, 'time', itvar)
      call check_ok(__FILE__,__LINE__,'variable time miss', 'ICBC FILE')
      istatus = nf90_get_att(ibcid, itvar, 'units', icbc_timeunits)
      call check_ok(__FILE__,__LINE__,'variable time units miss','ICBC FILE')
      istatus = nf90_get_att(ibcid, itvar, 'calendar', icbc_timecal)
      call check_ok(__FILE__,__LINE__,'variable time calendar miss','ICBC FILE')
      allocate(icbc_nctime(ibcnrec), stat=istatus)
      if ( istatus /= 0 ) then
        write(6,*) 'Memory allocation error in ICBC for time real values'
        call fatal(__FILE__,__LINE__,'ICBC READ')
      end if
      allocate(chbc_idate(ibcnrec), stat=istatus)
      if ( istatus /= 0 ) then
        write(6,*) 'Memory allocation error in ICBC for time array'
        call fatal(__FILE__,__LINE__,'ICBC READ')
      end if
      istatus = nf90_get_var(ibcid, itvar, icbc_nctime)
      call check_ok(__FILE__,__LINE__,'variable time read error', 'ICBC FILE')
      do i = 1 , ibcnrec
        chbc_idate(i) = timeval2date(icbc_nctime(i), &
                                     icbc_timeunits,icbc_timecal)
      end do
      if ( ibcnrec > 1 ) then
        chkdiff = idnint(icbc_nctime(2) - icbc_nctime(1))
        if (chkdiff /= ibdyfrq) then
          write (6,*) 'Time var in ICBC inconsistency.'
          write (6,*) 'Expecting ibdyfrq = ', ibdyfrq
          write (6,*) 'Found     ibdyfrq = ', chkdiff
          call fatal(__FILE__,__LINE__,'ICBC READ')
        end if
      end if
      deallocate(icbc_nctime)
      if ( igaschem == 1 ) then
        do i = 1 , n_chbcvar
          istatus = nf90_inq_varid(ichin, trim(chbcname(i)), chbc_ivar(i))
          call check_ok(__FILE__,__LINE__, &
               'variable '//trim(chbcname(i))//' missing','CHBC FILE ERROR')
        end do
      end if
      if ( iaerosol == 1 ) then
        do i = 1 , n_aebcvar
          istatus = nf90_inq_varid(iaein, trim(aebcname(i)), aebc_ivar(i))
          call check_ok(__FILE__,__LINE__, &
               'variable '//trim(aebcname(i))//' missing','AEBC FILE ERROR')
        end do
      end if
      if ( ioxclim == 1 ) then
        do i = 1 , n_oxbcvar
          istatus = nf90_inq_varid(ioxin, trim(oxbcname(i)), oxbc_ivar(i))
          call check_ok(__FILE__,__LINE__, &
               'variable '//trim(oxbcname(i))//' missing','OXBC FILE ERROR')
        end do
      end if
    end subroutine open_chbc

    subroutine read_chbc(chebdio)
      implicit none
      real(rk8) , dimension (:,:,:,:), intent(out) :: chebdio 
      integer(ik4) , dimension(4) :: istart , icount
      real(rk4) , dimension(jx,iy,kz) :: xread
      integer(ik4) :: i , j , k, n , iafter
      istart(4) = ibcrec
      istart(3) = 1
      istart(2) = 1
      istart(1) = 1
      icount(4) = 1
      icount(3) = kz
      icount(2) = iy
      icount(1) = jx
      iafter = 0
      if ( igaschem == 1 ) then
        do n = 1 , n_chbcvar
          istatus = nf90_get_var(ichin, chbc_ivar(n), xread, istart, icount)
          call check_ok(__FILE__,__LINE__, &
               'variable '//trim(chbcname(n))//' read error','CHBC FILE ERROR')
          do k = 1 , kz
            do j = 1 , jx
              do i = 1 , iy
                chebdio(j,i,k,n) = xread(j,i,k)
              end do
            end do
          end do
          iafter = iafter + 1
        end do
      end if
      if ( iaerosol == 1 ) then
        do n = 1 , n_aebcvar
          istatus = nf90_get_var(iaein, aebc_ivar(n), xread, istart, icount)
          call check_ok(__FILE__,__LINE__, &
               'variable '//trim(aebcname(n))//' read error','AEBC FILE ERROR')
          do k = 1 , kz
            do j = 1 , jx
              do i = 1 , iy
                chebdio(j,i,k,iafter+1) = xread(j,i,k)
              end do
            end do
          end do
        iafter = iafter + 1
        end do
      end if
      if ( ioxclim == 1 ) then
        do n = 1 , n_oxbcvar
          istatus = nf90_get_var(ioxin, oxbc_ivar(n), xread, istart, icount)
          call check_ok(__FILE__,__LINE__, &
               'variable '//trim(oxbcname(n))//' read error','OXBC FILE ERROR')
          do k = 1 , kz
            do j = 1 , jx
              do i = 1 , iy
                chebdio(j,i,k,iafter+n) = xread(j,i,k)
              end do
            end do
          end do
        end do
      end if

     where (chebdio < d_zero) chebdio = d_zero

    end subroutine read_chbc

    subroutine close_chbc
      implicit none
      if ( ichin >= 0 ) then
        istatus = nf90_close(ichin)
        call check_ok(__FILE__,__LINE__, 'Error Close file', 'CHBC FILE ERROR')
        ichin = -1
      end if
      if ( iaein >= 0 ) then
        istatus = nf90_close(iaein)
        call check_ok(__FILE__,__LINE__, 'Error Close file', 'AEBC FILE ERROR')
        iaein = -1
      end if
      if ( allocated(chbc_idate) ) deallocate(chbc_idate)
    end subroutine close_chbc

    subroutine check_ok(f,l,m1,mf)
      implicit none
      character(*) , intent(in) :: f, m1 , mf
      integer(ik4) , intent(in) :: l
      if (istatus /= nf90_noerr) then
        write (6,*) trim(m1)
        write (6,*) nf90_strerror(istatus)
        call fatal(f,l,trim(mf))
      end if
    end subroutine check_ok

end module mod_che_ncio
