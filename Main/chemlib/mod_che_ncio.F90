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
  public :: prepare_chem_out , init_mod_che_ncio , writerec_che2
  public :: open_chbc , close_chbc , chbc_search , read_chbc

  public :: chbc_ivar,n_chbcvar,n_aebcvar

  integer(ik4) :: istatus
  integer(ik4) :: recc

  integer(ik4) , parameter :: n_chevar = 20
   integer(ik4) , parameter :: n_oxbcvar = 5
  integer(ik4) , parameter :: n_optvar = 10
  integer(ik4) , parameter :: n_chbcvar = 25
  integer(ik4) :: n_aebcvar
  integer(ik4) :: ichin  , iaein, ioxin

  integer(ik4) , dimension(:) , pointer :: ncche     
  integer(ik4) , dimension(n_chevar) :: ichevar
  integer(ik4) , dimension(n_optvar) ::ioptvar 

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
  type(rcm_time_and_date) , save :: icherefdate
  integer(ik4) , dimension(9) :: idims 
  integer(ik4) ::idmin , icherec , ioptrec
  integer(ik4) :: ibcrec , ibcnrec
  real(rk8) :: tpd , cfd
  real(rk8) :: rpt

  integer , public , parameter :: ifrqmon = 1
  integer , public , parameter :: ifrqday = 2
  integer , public , parameter :: ifrqhrs = 3

  integer(ik4) :: o_is
  integer(ik4) :: o_ie
  integer(ik4) :: o_js
  integer(ik4) :: o_je
  integer(ik4) :: o_ni
  integer(ik4) :: o_nj
  integer(ik4) :: o_nz
  logical :: lwrap  
  character(256) :: dname , icbcname
  real(rk4) , dimension(:,:) , pointer :: ioxlat
  real(rk4) , dimension(:,:) , pointer :: ioxlon
  real(rk4) , dimension(:,:) , pointer :: iotopo
  real(rk4) , dimension(:,:) , pointer :: iomask
  real(rk4) , dimension(:,:,:) , pointer :: dumio

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

    subroutine init_mod_che_ncio(chemsymtype,lcband)
      implicit none
      character(len=8) , intent(in) :: chemsymtype
      logical , intent(in) :: lcband
      dname = trim(dirter)//pthsep//trim(domname)//'_DOMAIN000.nc'
      if (lcband) then
        o_is = 2
        o_ie = iy-2
        o_js = 1
        o_je = jx
        o_ni = iy-3
        o_nj = jx
        o_nz = kz
        lwrap = .true.
      else
        o_is = 2
        o_ie = iy-2
        o_js = 2
        o_je = jx-2
        o_ni = iy-3
        o_nj = jx-3
        o_nz = kz
        lwrap = .false.
      end if

      call getmem2d(ioxlat,1,o_nj,1,o_ni,'ncio:ioxlat')
      call getmem2d(ioxlon,1,o_nj,1,o_ni,'ncio:ioxlon')
      call getmem2d(iotopo,1,o_nj,1,o_ni,'ncio:iotopo')
      call getmem2d(iomask,1,o_nj,1,o_ni,'ncio:iomask')
      call getmem3d(dumio,1,o_nj,1,o_ni,1,o_nz,'ncio:dumio')

      ioxlat(:,:) = real(mddom_io%xlat(o_js:o_je,o_is:o_ie))
      ioxlon(:,:) = real(mddom_io%xlon(o_js:o_je,o_is:o_ie))
      iotopo(:,:) = real(mddom_io%ht(o_js:o_je,o_is:o_ie))
      iomask(:,:) = real(mddom_io%mask(o_js:o_je,o_is:o_ie))
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

!------------------------------------------------------------------------------      
!       IROUTINE: prepare_chem_out
!       SUBROUTINE INTERFACE:

    subroutine prepare_chem_out(idate, ifrest)
    
!   !DESCRIPTION:
!   prepare the dimensions variables, write global attributes
!   define the chemistry variables

      implicit none
      type(rcm_time_and_date) , intent(in) :: idate
      logical, intent(in) :: ifrest
      integer  :: itr
      character(8) ::chevarnam
      character(32) :: fbname
      character(16) :: fterr
      character(256) :: ofname
      real(rk4) :: hptop
      real(rk4) , dimension(iysg) :: yiy
      real(rk4) , dimension(jxsg) :: xjx
      integer(ik4) :: ncid
      integer(ik4) , dimension(3) :: izvar
      integer(ik4) , dimension(2) :: ivvar
      integer(ik4) , dimension(5) :: illtpvar
      integer(ik4) :: itvar , i , j
      integer(ik4) :: noutf
      character(len=129) :: cdum

      integer(ik4) , dimension(9) :: tyx
      integer(ik4) , dimension(9) :: tzyx
          
      character(len=36) :: ctime

      ctime = tochar(idate)

      ! total number of output  = ntr + 1 for optical properties file
      noutf = ntr
      if ( iaerosol == 1 ) noutf = ntr + 1

      if ( .not. associated(ncche) ) then
        call getmem1d(ncche,1,noutf,'che_ncio:ncche')
      end if
      do itr = 1, noutf
        if ( ncche(itr) >= 0 ) then
          istatus = nf90_close(ncche(itr))
          call check_ok(__FILE__,__LINE__,'Error close chem file output','CHE')
          ncche(itr) = -1
        end if
      end do
      ! tracer loop , since we are generating one output per
      ! tracer + 1 output for OPT
      do itr = 1, noutf
        if ( itr < noutf ) then  
          ncid = ncche(itr)     
          chevarnam =  chtrname(itr)
        else if ( itr == noutf) then
          chevarnam = 'OPT'
        end if 
        icherefdate = idate
        icherec = 1

        write (fterr, '(a3,a)') chevarnam, ' FILE'
        write (fbname,'(a,a,i10)') trim(chevarnam), '.', toint10(idate)
        ofname = trim(dirout)//pthsep//trim(domname)// &
                   '_'//trim(fbname)//'.nc'

        write (aline, *) 'Opening new output file ', trim(ofname)
        call say
#ifdef DEBUG
        write (aline,*) 'NOUTF = ', noutf, ', NTR = ', ntr
        call say
#endif

        call createfile_withname(ofname,ncid)
        call add_common_global_params(ncid,'Model (Chemistry '// &
                        trim(chevarnam)//')',.false.)
        ncche(itr) = ncid
!
!       ADD RUN PARAMETERS
!
        istatus = nf90_put_att(ncid,nf90_global,'model_IPCC_scenario',scenario)
        call check_ok(__FILE__,__LINE__,'Error add scenario', fterr)

        call cdumlbcs(cdum)
        istatus = nf90_put_att(ncid, nf90_global,  &
               'model_boundary_conditions' , trim(cdum))
        call check_ok(__FILE__,__LINE__,'Error add lbcs', fterr)
        call cdumcums(cdum)
        istatus = nf90_put_att(ncid, nf90_global,  &
                'model_cumulous_convection_scheme' , trim(cdum))
        call check_ok(__FILE__,__LINE__,'Error add icup', fterr)
        if ( icup == 2 .or. icup == 99 .or. icup == 98 ) then
          call cdumcumcl(cdum)
          istatus = nf90_put_att(ncid, nf90_global,  &
                'model_convective_closure_assumption' , trim(cdum))
          call check_ok(__FILE__,__LINE__,'Error add igcc', fterr)
        end if
        call cdumpbl(cdum)
        istatus = nf90_put_att(ncid, nf90_global,  &
                'model_boundary_layer_scheme' , trim(cdum))
        call check_ok(__FILE__,__LINE__,'Error add ibltyp', fterr)
        call cdummoist(cdum)
        istatus = nf90_put_att(ncid, nf90_global,  &
                'model_moist_physics_scheme' , trim(cdum))
        call check_ok(__FILE__,__LINE__,'Error add ipptls', fterr)
        call cdumocnflx(cdum)
        istatus = nf90_put_att(ncid, nf90_global,  &
                'model_ocean_flux_scheme' , trim(cdum))
        call check_ok(__FILE__,__LINE__,'Error add iocnflx', fterr)
        call cdumpgfs(cdum)
        istatus = nf90_put_att(ncid, nf90_global,  &
                'model_pressure_gradient_force_scheme' , trim(cdum))
        call check_ok(__FILE__,__LINE__,'Error add ipgf', fterr)
        call cdumemiss(cdum)
        istatus = nf90_put_att(ncid, nf90_global,  &
                'model_use_emission_factor' , trim(cdum))
        call check_ok(__FILE__,__LINE__,'Error add iemiss', fterr)
        call cdumlakes(cdum)
        istatus = nf90_put_att(ncid, nf90_global,  &
                'model_use_lake_model' , trim(cdum))
        call check_ok(__FILE__,__LINE__,'Error add lakemod', fterr)
        call cdumchems(cdum)
        istatus = nf90_put_att(ncid, nf90_global,  &
                'model_chemistry' , trim(cdum))
        call check_ok(__FILE__,__LINE__,'Error add ichem', fterr)
        call cdumdcsst(cdum)
        istatus = nf90_put_att(ncid, nf90_global,  &
                'model_diurnal_cycle_sst' , trim(cdum))
        call check_ok(__FILE__,__LINE__,'Error add dcsst', fterr)
        call cdumseaice(cdum)
        istatus = nf90_put_att(ncid, nf90_global,  &
                'model_seaice_effect' , trim(cdum))
        call check_ok(__FILE__,__LINE__,'Error add seaice', fterr)
        call cdumdesseas(cdum)
        istatus = nf90_put_att(ncid, nf90_global,  &
                'model_seasonal_desert_albedo_effect' , trim(cdum))
        call check_ok(__FILE__,__LINE__,'Error add desseas', fterr)
        istatus = nf90_put_att(ncid, nf90_global,  &
                 'model_simulation_initial_start' , tochar(globidate1))
        call check_ok(__FILE__,__LINE__,'Error add globidate1', fterr)

        istatus = nf90_put_att(ncid, nf90_global,  &
                'model_simulation_start' , tochar(idate1))
        call check_ok(__FILE__,__LINE__,'Error add idate1', fterr)
        istatus = nf90_put_att(ncid, nf90_global,  &
                'model_simulation_expected_end' , tochar(idate2))
        call check_ok(__FILE__,__LINE__,'Error add idate2', fterr)

        if ( ifrest ) then
          cdum = 'Yes'
        else
          cdum = 'No'
        end if
        istatus = nf90_put_att(ncid, nf90_global,  &
                'model_simulation_is_a_restart' , trim(cdum))
        call check_ok(__FILE__,__LINE__,'Error add ifrest', fterr)

        istatus = nf90_put_att(ncid, nf90_global,  &
                 'model_timestep_in_seconds' , dt)
        call check_ok(__FILE__,__LINE__,'Error add dt', fterr)
        istatus = nf90_put_att(ncid, nf90_global,  &
                'model_timestep_in_minutes_solar_rad_calc' , dtrad)
        call check_ok(__FILE__,__LINE__,'Error add dtrad', fterr)
        istatus = nf90_put_att(ncid, nf90_global,  &
               'model_timestep_in_seconds_bats_calc' , dtsrf)
        call check_ok(__FILE__,__LINE__,'Error add dtsrf', fterr)
        istatus = nf90_put_att(ncid, nf90_global,  &
                'model_timestep_in_hours_radiation_calc' , dtabem)
        call check_ok(__FILE__,__LINE__,'Error add dtabem', fterr)
        istatus = nf90_put_att(ncid, nf90_global,  &
                'model_timestep_in_hours_boundary_input' , ibdyfrq)
        call check_ok(__FILE__,__LINE__,'Error add dtabem', fterr)
!
!       End of Global Attributes
!
!       ADD DIMENSIONS
!
        istatus = nf90_def_dim(ncid, 'iy', o_ni, idims(2))
        call check_ok(__FILE__,__LINE__,'Error create dim iy', fterr)
        istatus = nf90_def_dim(ncid, 'jx', o_nj, idims(1))
        call check_ok(__FILE__,__LINE__,'Error create dim jx', fterr)
        istatus = nf90_def_dim(ncid, 'time', nf90_unlimited, idims(3))
        call check_ok(__FILE__,__LINE__,'Error create dim time', fterr)
        istatus = nf90_def_dim(ncid, 'kz', kz, idims(4))
        call check_ok(__FILE__,__LINE__,'Error create dim kz', fterr)
!
!       OUT TYPE DEPENDENT DIMENSIONS
!
        istatus = nf90_def_var(ncid, 'sigma', nf90_float, idims(4), izvar(1))
        call check_ok(__FILE__,__LINE__,'Error add var sigma', fterr)
        istatus = nf90_put_att(ncid, izvar(1), 'standard_name', &
                         'atmosphere_sigma_coordinate')
        call check_ok(__FILE__,__LINE__,'Error add sigma standard_name', fterr)
        istatus = nf90_put_att(ncid, izvar(1), 'long_name', &
                         'Sigma at model layers')
        call check_ok(__FILE__,__LINE__,'Error add sigma long_name', fterr)
        istatus = nf90_put_att(ncid, izvar(1), 'units', '1')
        call check_ok(__FILE__,__LINE__,'Error add sigma units', fterr)
        istatus = nf90_put_att(ncid, izvar(1), 'axis', 'Z')
        call check_ok(__FILE__,__LINE__,'Error add sigma axis', fterr)
        istatus = nf90_put_att(ncid, izvar(1), 'positive', 'down')
        call check_ok(__FILE__,__LINE__,'Error add sigma positive', fterr)
        istatus = nf90_put_att(ncid, izvar(1), 'formula_terms',  &
                         'sigma: sigma ps: ps ptop: ptop')
        call check_ok(__FILE__,__LINE__,'Error add sigma formula_terms', fterr)
   
        istatus = nf90_def_var(ncid, 'ptop', nf90_float, varid=izvar(2))
        call check_ok(__FILE__,__LINE__,'Error add var ptop', fterr)
        istatus = nf90_put_att(ncid, izvar(2), 'standard_name', 'air_pressure')
        call check_ok(__FILE__,__LINE__,'Error add ptop standard_name', fterr)
        istatus = nf90_put_att(ncid, izvar(2), 'long_name', &
                               'Pressure at model top')
        call check_ok(__FILE__,__LINE__,'Error add ptop long_name', fterr)
        istatus = nf90_put_att(ncid, izvar(2), 'units', 'hPa')
        call check_ok(__FILE__,__LINE__,'Error add ptop units', fterr)
        istatus = nf90_def_var(ncid, 'iy', nf90_float, idims(2), ivvar(1))
        call check_ok(__FILE__,__LINE__,'Error add var iy', fterr)
        istatus = nf90_put_att(ncid, ivvar(1), 'standard_name', &
                           'projection_y_coordinate')
        call check_ok(__FILE__,__LINE__,'Error add iy standard_name', fterr)
        istatus = nf90_put_att(ncid, ivvar(1), 'long_name', &
                           'y-coordinate in Cartesian system')
        call check_ok(__FILE__,__LINE__,'Error add iy long_name', fterr)
        istatus = nf90_put_att(ncid, ivvar(1), 'units', 'km')
        call check_ok(__FILE__,__LINE__,'Error add iy units', fterr)
        istatus = nf90_def_var(ncid, 'jx', nf90_float, idims(1), ivvar(2))
        call check_ok(__FILE__,__LINE__,'Error add var jx', fterr)
        istatus = nf90_put_att(ncid, ivvar(2), 'standard_name', &
                         'projection_x_coordinate')
        call check_ok(__FILE__,__LINE__,'Error add jx standard_name', fterr)
        istatus = nf90_put_att(ncid, ivvar(2), 'long_name', &
                         'x-coordinate in Cartesian system')
        call check_ok(__FILE__,__LINE__,'Error add jx long_name', fterr)
        istatus = nf90_put_att(ncid, ivvar(2), 'units', 'km')
        call check_ok(__FILE__,__LINE__,'Error add jx units', fterr)
        istatus = nf90_def_var(ncid,'xlat',nf90_float,idims(1:2),illtpvar(1))
        call check_ok(__FILE__,__LINE__,'Error add var xlat', fterr)
        istatus = nf90_put_att(ncid, illtpvar(1), 'standard_name', 'latitude')
        call check_ok(__FILE__,__LINE__,'Error add xlat standard_name', fterr)
        istatus = nf90_put_att(ncid, illtpvar(1), 'long_name', &
                         'Latitude at cross points')
        call check_ok(__FILE__,__LINE__,'Error add xlat long_name', fterr)
        istatus = nf90_put_att(ncid, illtpvar(1), 'units', 'degrees_north')
        call check_ok(__FILE__,__LINE__,'Error add xlat units', fterr)
        istatus = nf90_def_var(ncid,'xlon',nf90_float,idims(1:2),illtpvar(2))
        call check_ok(__FILE__,__LINE__,'Error add var xlon', fterr)
        istatus = nf90_put_att(ncid, illtpvar(2), 'standard_name', 'longitude')
        call check_ok(__FILE__,__LINE__,'Error add xlon standard_name', fterr)
        istatus = nf90_put_att(ncid, illtpvar(2), 'long_name', &
                         'Longitude at cross points')
        call check_ok(__FILE__,__LINE__,'Error add xlon long_name', fterr)
        istatus = nf90_put_att(ncid, illtpvar(2), 'units', 'degrees_east')
        call check_ok(__FILE__,__LINE__,'Error add xlon units', fterr)
        istatus = nf90_def_var(ncid,'topo',nf90_float,idims(1:2),illtpvar(3))
        call check_ok(__FILE__,__LINE__,'Error add var topo', fterr)
        istatus = nf90_put_att(ncid, illtpvar(3), 'standard_name', &
                         'surface_altitude')
        call check_ok(__FILE__,__LINE__,'Error add topo standard_name', fterr)
        istatus = nf90_put_att(ncid, illtpvar(3), 'long_name',     &
                         'Domain surface elevation')
        call check_ok(__FILE__,__LINE__,'Error add topo long_name', fterr)
        istatus = nf90_put_att(ncid, illtpvar(3), 'units', 'm')
        call check_ok(__FILE__,__LINE__,'Error add topo units', fterr)
        istatus = nf90_put_att(ncid, illtpvar(3), 'coordinates', 'xlat xlon')
        call check_ok(__FILE__,__LINE__,'Error add topo coord', fterr)
        istatus = nf90_put_att(ncid, illtpvar(3), 'grid_mapping', 'rcm_map')
        call check_ok(__FILE__,__LINE__,'Error add topo grid_mapping', fterr)
        istatus = nf90_def_var(ncid, 'mask',nf90_float,idims(1:2),illtpvar(4))
        call check_ok(__FILE__,__LINE__,'Error add var mask', fterr)
        istatus = nf90_put_att(ncid, illtpvar(4), 'standard_name', 'landmask')
        call check_ok(__FILE__,__LINE__,'Error add mask standard_name', fterr)
        istatus = nf90_put_att(ncid, illtpvar(4), 'long_name',     &
                         'Domain land/ocean mask')
        call check_ok(__FILE__,__LINE__,'Error add mask long_name', fterr)
        istatus = nf90_put_att(ncid, illtpvar(4), 'units', '1')
        call check_ok(__FILE__,__LINE__,'Error add mask units', fterr)
        istatus = nf90_put_att(ncid, illtpvar(4), 'coordinates', 'xlat xlon')
        call check_ok(__FILE__,__LINE__,'Error add mask coord', fterr)
        istatus = nf90_put_att(ncid, illtpvar(4), 'grid_mapping', 'rcm_map')
        call check_ok(__FILE__,__LINE__,'Error add mask grid_mapping', fterr)
        istatus = nf90_def_var(ncid, 'time', nf90_double, idims(3:3), itvar)
        call check_ok(__FILE__,__LINE__,'Error add var time', fterr)
        istatus = nf90_put_att(ncid, itvar, 'standard_name', 'time')
        call check_ok(__FILE__,__LINE__,'Error add time standard_name', fterr)
        istatus = nf90_put_att(ncid, itvar, 'long_name', 'time')
        call check_ok(__FILE__,__LINE__,'Error add time long_name', fterr)
        istatus = nf90_put_att(ncid, itvar, 'calendar', calstr(idate%calendar))
        call check_ok(__FILE__,__LINE__,'Error add time calendar', fterr)
        istatus = nf90_put_att(ncid, itvar, 'units', 'hours since '//ctime)
        call check_ok(__FILE__,__LINE__,'Error add time units', fterr)
   
        istatus = nf90_def_var(ncid, 'ps', nf90_float, idims(1:3), illtpvar(5))
        call check_ok(__FILE__,__LINE__,'Error add var ps', fterr)
        istatus = nf90_put_att(ncid, illtpvar(5), 'standard_name', &
                         'surface_air_pressure')
        call check_ok(__FILE__,__LINE__,'Error add ps standard_name', fterr)
        istatus = nf90_put_att(ncid, illtpvar(5),'long_name','Surface pressure')
        call check_ok(__FILE__,__LINE__,'Error add ps long_name', fterr)
        istatus = nf90_put_att(ncid, illtpvar(5), 'units', 'hPa')
        call check_ok(__FILE__,__LINE__,'Error add ps units', fterr)
        istatus = nf90_put_att(ncid, illtpvar(5), 'coordinates', 'xlat xlon')
        call check_ok(__FILE__,__LINE__,'Error add ps coord', fterr)
        istatus = nf90_put_att(ncid, illtpvar(5), 'grid_mapping', 'rcm_map')
        call check_ok(__FILE__,__LINE__,'Error add ps grid_mapping', fterr)
        istatus = nf90_put_att(ncid, illtpvar(5), 'cell_methods', 'time: point')
        call check_ok(__FILE__,__LINE__,'Error add ps cell_methods', fterr)

        tyx = (/idims(1),idims(2),idims(3),-1,-1,-1,-1,-1,-1/)
        tzyx = (/idims(1),idims(2),idims(4),idims(3),-1,-1,-1,-1,-1/)

        if (  itr < noutf  ) then      
          ichevar = -1
          ichevar(1) = itvar
          ichevar(2) = illtpvar(5)
          call ch_addvara(ncid,chevarnam,chevarnam, &
                    'atmosphere_mixing_ratio_of_tracer', &
                    'Tracers mixing ratios','kg kg-1', &
                    tzyx,.false.,ichevar(3))
          call ch_addvara(ncid,chevarnam,'wetdep_rainout_flx', &
                    'wet_deposition_from_large_scale_precip', &
                    'Wet deposition LS','mg/m2/d', &
                    tyx,.false.,ichevar(4))
          call ch_addvara(ncid,chevarnam,'wetdep_washout_flx', &
                    'wet_deposition_from_convective_precip', &
                    'Wet deposition CONV','mg/m2/d', &
                    tyx,.false.,ichevar(5))
          call ch_addvara(ncid,chevarnam,'drydep_flx', &
                    'dry_deposition', &
                    'Dry deposition rate','mg/m2/d', &
                    tyx,.false.,ichevar(6))
          call ch_addvara(ncid,chevarnam,'emiss_flx', &
                    'surface_emission_rate', &
                    'Emission rate','mg/m2/d', &
                    tyx,.false.,ichevar(7))
          call ch_addvara(ncid,chevarnam,'drydep_vel', &
                    'dry deposition velocity', &
                    'dr. dep. vel','m.s-1', &
                    tyx,.false.,ichevar(8))
          call ch_addvara(ncid,chevarnam,'trac_burden', &
                    'tracer_burden', &
                    'trac bud ','mg.m-2', &
                    tyx,.false.,ichevar(9))

          if ( ichdiag == 1 ) then
            call ch_addvara(ncid,chevarnam,'chem_tend', &
                    'chemical_prod_loss', &
                    'chem tendency','kg kg-1 s-1', &
                    tzyx,.false.,ichevar(10))
            call ch_addvara(ncid,chevarnam,'advh_tend', &
                    'advect. hor. tend', &
                    'chem tendency','kg kg-1 s-1', &
                    tzyx,.false.,ichevar(11))
            call ch_addvara(ncid,chevarnam,'advv_tend', &
                    'advec. ver. tend', &
                    'chem tendency','kg kg-1 s-1', &
                    tzyx,.false.,ichevar(12))
            call ch_addvara(ncid,chevarnam,'difh_tend', &
                    'diff. hor. tend', &
                    'chem tendency','kg kg-1 s-1', &
                    tzyx,.false.,ichevar(13))
            call ch_addvara(ncid,chevarnam,'conv_tend', &
                    'convec_tend', &
                    'chem tendency','kg kg-1 s-1', &
                    tzyx,.false.,ichevar(14))
            call ch_addvara(ncid,chevarnam,'tubl_tend', &
                    'turb. vert. tend', &
                    'chem tendency','kg kg-1 s-1', &
                    tzyx,.false.,ichevar(15))
            call ch_addvara(ncid,chevarnam,'rainout_tend', &
                    'large + conv scal rainout removal tend', &
                    'chem tendency','kg kg-1 s-1', &
                    tzyx,.false.,ichevar(16))
            call ch_addvara(ncid,chevarnam,'washout_tend', &
                    'large + conv scale washout removal tend', &
                    'chem tendency','kg kg-1 s-1', &
                    tzyx,.false.,ichevar(17))
            call ch_addvara(ncid,chevarnam,'bdyc_tend', &
                    'bdy. cond. tend', &
                    'chem tendency','kg kg-1 s-1', &
                    tzyx,.false.,ichevar(18))
            call ch_addvara(ncid,chevarnam,'sedddp_tend', &
                    'sedim. ddep. tend', &
                    'chem tendency','kg kg-1 s-1', &
                    tzyx,.false.,ichevar(19))
            if ( ibltyp == 2 .or. ibltyp == 99 ) then
              call ch_addvara(ncid,chevarnam,'chifxuw', &
                      'UWpbl. net flux', &
                      'chem tendency','kg kg-1 s-1', &
                      tyx,.false.,ichevar(20))
            end if
          end if

         else if ( itr == noutf ) then 
           ioptvar = -1
           ioptvar(1) = itvar
           ioptvar(2) = illtpvar(4)
           call ch_addvara(ncid,chevarnam,'aext8', &
                    'aerosol_optical_depth', &
                    'aer mix. aod.','1',tzyx,.false.,ioptvar(3))
           call ch_addvara(ncid,chevarnam,'assa8', &
                    'aerosol_single_scattering_albedo', &
                    'aer mix. sin. scat. alb','1',tzyx,.false.,ioptvar(4))
           call ch_addvara(ncid,chevarnam,'agfu8', &
                    'aerosol_asymmetry_parameter', &
                    'aer mix. sin. scat. alb','1',tzyx,.false.,ioptvar(5))
           call ch_addvara(ncid,chevarnam,'acstoarf', &
                    'toa_instantaneous_shortwave_radiative_forcing', &
                    'TOArad SW forcing av.','W m-2', &
                    tyx,.false.,ioptvar(6))
           call ch_addvara(ncid,chevarnam,'acstsrrf', &
                    'surface_shortwave_radiative_forcing', &
                    'SRFrad SW forcing av.','W m-2', &
                    tyx,.false.,ioptvar(7))
           call ch_addvara(ncid,chevarnam,'acstalrf', &
                    'toa_longwave_radiative_forcing', &
                    'TOArad LW forcing av.','W m-2', &
                    tyx,.false.,ioptvar(8))
           call ch_addvara(ncid,chevarnam,'acssrlrf', &
                    'surface_longwave_radiative_forcing', &
                    'SRFrad LW forcing av.','W m-2', &
                    tyx,.false.,ioptvar(9))
           call ch_addvara(ncid,chevarnam,'aod', &
                    'aerosol optical depth vis band', &
                    'aer aod vis band',' ', &
                    tyx,.false.,ioptvar(10))
         end if 

        istatus = nf90_enddef(ncid)
        call check_ok(__FILE__,__LINE__, &
           'Error End Definitions NetCDF output',fterr)

        ! write variables which are not time dependant in the file

        istatus = nf90_put_var(ncid, izvar(1), hsigma)
        call check_ok(__FILE__,__LINE__,'Error var sigma write', fterr)
        hptop = real(ptop*d_10)
        istatus = nf90_put_var(ncid, izvar(2), hptop)
        call check_ok(__FILE__,__LINE__,'Error var ptop write', fterr)
 
        yiy(1) = -real((dble(o_ni-1)/2.0D0)*ds)
        xjx(1) = -real((dble(o_nj-1)/2.0D0)*ds)
        do i = 2 , o_ni
          yiy(i) = yiy(i-1)+real(ds)
        end do
        do j = 2 , o_nj
          xjx(j) = xjx(j-1)+real(ds)
        end do
        istatus = nf90_put_var(ncid, ivvar(1), yiy(1:o_ni))
        call check_ok(__FILE__,__LINE__,'Error var iy write', fterr)
        istatus = nf90_put_var(ncid, ivvar(2), xjx(1:o_nj))
        call check_ok(__FILE__,__LINE__,'Error var jx write', fterr)
        istatus = nf90_put_var(ncid, illtpvar(1), ioxlat)
        call check_ok(__FILE__,__LINE__,'Error var xlat write', fterr)
        istatus = nf90_put_var(ncid, illtpvar(2), ioxlon)
        call check_ok(__FILE__,__LINE__,'Error var xlon write', fterr)
        istatus = nf90_put_var(ncid, illtpvar(3), iotopo)
        call check_ok(__FILE__,__LINE__,'Error var topo write', fterr)
        istatus = nf90_put_var(ncid, illtpvar(4), iomask)
        call check_ok(__FILE__,__LINE__,'Error var mask write', fterr)
  
        istatus = nf90_sync(ncid)
        call check_ok(__FILE__,__LINE__,'Error sync file', fterr)

        ncche(itr) = ncid

        ! end of tracer loop       
      end do      
    end subroutine prepare_chem_out

!=====================================================================

    subroutine ch_addvara(ncid,ctype,vname,vst,vln,vuni,idims,lmiss,ivar)
      implicit none
      integer(ik4) , intent(in) :: ncid
      character(3) , intent(in) :: ctype
      character(len=*) , intent(in) :: vname
      character(len=*) , intent(in) :: vst , vln , vuni
      integer(ik4) , dimension(5) , intent(in) :: idims
      logical , intent(in) :: lmiss
      integer(ik4) , intent(out) :: ivar
      character(64) :: cdum
      real(rk4) , parameter :: fillv = +1E+20
      integer(ik4) :: i , ndims
       
      ndims = 0
      do i = 1 , 5
        if ( idims(i) > 0 ) ndims = ndims+1
      end do

      cdum = vname
      istatus = nf90_def_var(ncid, cdum, nf90_float, &
                             idims(1:ndims), ivar)
      call check_ok(__FILE__,__LINE__,'Error adding variable '//vname, &
                    ctype//' FILE ERROR')
#ifdef NETCDF4_HDF5
      istatus = nf90_def_var_deflate(ncid, ivar, 1, 1, 9)
      call check_ok(__FILE__,__LINE__, &
           'Error setting deflate on variable '//vname,ctype//' FILE ERROR')
#endif
      cdum = vst
      istatus = nf90_put_att(ncid, ivar, 'standard_name', cdum)
      call check_ok(__FILE__,__LINE__, &
            'Error adding '//vname//' standard_name',ctype//' FILE ERROR')
      cdum = vln
      istatus = nf90_put_att(ncid, ivar, 'long_name', cdum)
      call check_ok(__FILE__,__LINE__,'Error adding '//vname//' long_name', &
                     ctype//' FILE ERROR')
      cdum = vuni
      istatus = nf90_put_att(ncid, ivar, 'units', cdum)
      call check_ok(__FILE__,__LINE__,'Error adding '//vname//' units', &
                     ctype//' FILE ERROR')
      istatus = nf90_put_att(ncid, ivar, 'coordinates','xlat xlon')
      call check_ok(__FILE__,__LINE__,'Error adding '//vname//' coordinates', &
                     ctype//' FILE ERROR')
      istatus = nf90_put_att(ncid, ivar, 'grid_mapping','rcm_map')
      call check_ok(__FILE__,__LINE__,'Error adding '//vname//' grid_mapping', &
                    ctype//' FILE ERROR')
      if (lmiss) then
        istatus = nf90_put_att(ncid, ivar, '_FillValue',fillv)
        call check_ok(__FILE__,__LINE__,'Error adding '//vname//' coordinates',&
                      ctype//' FILE ERROR')
      end if
    end subroutine ch_addvara

!============================================================================

    subroutine writerec_che2(chia,dtrace,wdlsc,wdcvc,remdrd,cemtrac,drydepv, &
                             chemdiag,cadvhdiag,cadvvdiag,cdifhdiag,         &
                             cconvdiag,cbdydiag,ctbldiag,cseddpdiag,         &
                             ccuwdiag,remlsc,remcvc,ext,ssa,asp,aod,tarf,    &
                             ssrf,talwrf,srlwrf,ps,idate)
      implicit none
          
      type(rcm_time_and_date) , intent(in) :: idate
      real(rk8) , pointer , dimension(:,:,:,:) , intent(in) :: chia , &
             chemdiag , cadvhdiag , cadvvdiag , cdifhdiag , cconvdiag , &
             cbdydiag , ctbldiag , remlsc , remcvc  ,cseddpdiag
      real(rk8) , pointer , dimension(:,:,:) , intent(in) :: ccuwdiag
      real(rk8) , pointer , dimension(:,:) , intent(in) :: ps
      real(rk8) , pointer , dimension(:,:,:) , intent(in) :: wdlsc , wdcvc , &
                        remdrd , cemtrac , drydepv , dtrace
      real(rk8) , pointer , dimension(:,:,:) , intent(in) :: ext , ssa , asp 
      real(rk8) , pointer , dimension(:,:) , intent(in) :: tarf , ssrf , aod , &
                                                          talwrf , srlwrf
      integer(ik4) :: n , k , noutf
      integer(ik4) , dimension(5) :: istart , icount
      real(rk8) , dimension(1) :: nctime
      type(rcm_time_interval) :: tdif
      character(len=36) :: ctime
      real(rk8) :: cfd2

      noutf = ntr

      if ( iaerosol == 1 ) noutf = ntr + 1 
 
      do n = 1 , noutf        
        istart(1) = icherec
        icount(1) = 1
        ctime = tochar(idate)
        tdif = idate-icherefdate
        nctime(1) = tohours(tdif)
          
        istatus = nf90_put_var(ncche(n), ichevar(1), nctime, &
                               istart(1:1), icount(1:1))
        call check_ok(__FILE__,__LINE__, &
                      'Error writing itime '//ctime, 'CHE FILE ERROR')
        dumio(:,:,1) = real((ps(o_js:o_je,o_is:o_ie)+rpt)*10.0D0)
        istart(3) = icherec
        istart(2) = 1
        istart(1) = 1
        icount(3) = 1
        icount(2) = o_ni
        icount(1) = o_nj
        istatus = nf90_put_var(ncche(n), ichevar(2), &
                               dumio(:,:,1), istart(1:3), icount(1:3))
        call check_ok(__FILE__,__LINE__, &
                      'Error writing ps at '//ctime, 'CHE FILE ERROR')

        if ( n < noutf ) then 
          istart(4) = icherec
          istart(3) = 1
          istart(2) = 1
          istart(1) = 1
          icount(5) = 1
          icount(4) = 1
          icount(3) = o_nz
          icount(2) = o_ni
          icount(1) = o_nj

          !*** tracer concentration
          do k = 1 , kz
            dumio(:,:,k) = real(chia(o_js:o_je,o_is:o_ie,k,n) / &
                                ps(o_js:o_je,o_is:o_ie))
          end do
          istatus = nf90_put_var(ncche(n), ichevar(3), &
                                 dumio, istart, icount)
          call check_ok(__FILE__,__LINE__, &
               'Error writing '//chtrname(n)//' at '//ctime,'CHE FILE ERROR')
          
          !*** output 2-D gas chem fields
          istart(4) = 1
          istart(3) = icherec
          istart(2) = 1
          istart(1) = 1
          icount(4) = 1
          icount(3) = 1
          icount(2) = o_ni
          icount(1) = o_nj

          ! accumulated quantities between two output steps are converted
          ! to deposition/emission mean rate (mg /m2/per day)  
          cfd = 1.D6 * 86400.D0 ! from kg/m2/s to mg/m2/day
          cfd2 = dtche / (chemfrq *3600.0D0)

          !*** wet deposition from large-scale precip
          dumio(:,:,1) = real(wdlsc(o_js:o_je,o_is:o_ie,n)* cfd)
          istatus = nf90_put_var(ncche(n), ichevar(4), &
                                 dumio(:,:,1), istart(1:3), icount(1:3))
          call check_ok(__FILE__,__LINE__, &
               'Error writing wet dep LS at '//ctime,'CHE FILE ERROR')

          !*** wet deposition from convective precip
          dumio(:,:,1) = real(wdcvc(o_js:o_je,o_is:o_ie,n)* cfd)
          istatus = nf90_put_var(ncche(n), ichevar(5), &
                                 dumio(:,:,1), istart(1:3), icount(1:3))
          call check_ok(__FILE__,__LINE__, &
               'Error writing wet dep CONV at '//ctime,'CHE FILE ERROR')

          !*** dry deposition
          dumio(:,:,1) = real(remdrd(o_js:o_je,o_is:o_ie,n)*cfd)
          istatus = nf90_put_var(ncche(n), ichevar(6), &
                                 dumio(:,:,1), istart(1:3), icount(1:3))
          call check_ok(__FILE__,__LINE__, &
               'Error writing dry dep '//ctime, 'CHE FILE ERROR')

          !*** emission rates
          dumio(:,:,1) = real(cemtrac(o_js:o_je,o_is:o_ie,n)*cfd)
          istatus = nf90_put_var(ncche(n), ichevar(7), &
                                 dumio(:,:,1), istart(1:3), icount(1:3))
          call check_ok(__FILE__,__LINE__, &
               'Error writing emission rate '//ctime, 'CHE FILE ERROR')

          !*** dry dep vel 
          dumio(:,:,1) = real(drydepv(o_js:o_je,o_is:o_ie,n)*cfd2)
          istatus = nf90_put_var(ncche(n), ichevar(8), &
                                 dumio(:,:,1), istart(1:3), icount(1:3))
          call check_ok(__FILE__,__LINE__, &
               'Error writing dr.dep.vel '//ctime, 'CHE FILE ERROR')

          !*** tracer burden (instantaneous in mg/m2) 
          dumio(:,:,1) = real(dtrace(o_js:o_je,o_is:o_ie,n)*1.0D6)
          istatus = nf90_put_var(ncche(n), ichevar(9), &
                                 dumio(:,:,1), istart(1:3), icount(1:3))
          call check_ok(__FILE__,__LINE__, &
               'Error writing trac burden '//ctime, 'CHE FILE ERROR')

          if ( ichdiag == 1 ) then 
            !*** 3D tracer diagnostic : chemical productio/loss
            istart(4) = icherec
            istart(3) = 1
            istart(2) = 1
            istart(1) = 1
            icount(5) = 1
            icount(4) = 1
            icount(3) = o_nz
            icount(2) = o_ni
            icount(1) = o_nj
            do k = 1 , kz
              dumio(:,:,k) = real(chemdiag(o_js:o_je,o_is:o_ie,k,n) / &
                                  ps(o_js:o_je,o_is:o_ie))
            end do
            istatus = nf90_put_var(ncche(n), ichevar(10), &
                                   dumio, istart, icount)
            call check_ok(__FILE__,__LINE__, &
                         'Error writing chemdiag at '//ctime,'CHE FILE ERROR')
            do k = 1 , kz
              dumio(:,:,k) = real(cadvhdiag(o_js:o_je,o_is:o_ie,k,n) / &
                                  ps(o_js:o_je,o_is:o_ie))
            end do
            istatus = nf90_put_var(ncche(n), ichevar(11), &
                                   dumio, istart, icount)
            call check_ok(__FILE__,__LINE__, &
                         'Error writing cadvhdiag at '//ctime,'CHE FILE ERROR')
            do k = 1 , kz
              dumio(:,:,k) = real(cadvvdiag(o_js:o_je,o_is:o_ie,k,n) / &
                                  ps(o_js:o_je,o_is:o_ie))
            end do
            istatus = nf90_put_var(ncche(n), ichevar(12), &
                                   dumio, istart, icount)
            call check_ok(__FILE__,__LINE__, &
                         'Error writing cadvvdiag at '//ctime,'CHE FILE ERROR')
            do k = 1 , kz
              dumio(:,:,k) = real(cdifhdiag(o_js:o_je,o_is:o_ie,k,n) / &
                                  ps(o_js:o_je,o_is:o_ie))
            end do
            istatus = nf90_put_var(ncche(n), ichevar(13), &
                                   dumio, istart, icount)
            call check_ok(__FILE__,__LINE__, &
                         'Error writing cdifhdiag at '//ctime,'CHE FILE ERROR')
            do k = 1 , kz
              dumio(:,:,k) = real(cconvdiag(o_js:o_je,o_is:o_ie,k,n) / &
                                  ps(o_js:o_je,o_is:o_ie))
            end do
            istatus = nf90_put_var(ncche(n), ichevar(14), &
                                   dumio, istart, icount)
            call check_ok(__FILE__,__LINE__, &
                          'Error writing cconvdiag at '//ctime,'CHE FILE ERROR')
            do k = 1 , kz
              dumio(:,:,k) = real(ctbldiag(o_js:o_je,o_is:o_ie,k,n) / &
                                  ps(o_js:o_je,o_is:o_ie))
            end do
            istatus = nf90_put_var(ncche(n), ichevar(15), &
                                   dumio, istart, icount)
            call check_ok(__FILE__,__LINE__, &
                         'Error writing ctbldiag at '//ctime,'CHE FILE ERROR')
            do k = 1 , kz
              dumio(:,:,k) = real(remlsc(o_js:o_je,o_is:o_ie,k,n) / &
                                  ps(o_js:o_je,o_is:o_ie))
            end do
            istatus = nf90_put_var(ncche(n), ichevar(16), &
                                   dumio, istart, icount)
            call check_ok(__FILE__,__LINE__, &
                         'Error writing remlsc at '//ctime,'CHE FILE ERROR')
            do k = 1 , kz
              dumio(:,:,k) = real(remcvc(o_js:o_je,o_is:o_ie,k,n) / &
                                  ps(o_js:o_je,o_is:o_ie))
            end do
            istatus = nf90_put_var(ncche(n), ichevar(17), &
                                   dumio, istart, icount)
            call check_ok(__FILE__,__LINE__, &
                          'Error writing remcvc at '//ctime,'CHE FILE ERROR')
            do k = 1 , kz
              dumio(:,:,k) = real(cbdydiag(o_js:o_je,o_is:o_ie,k,n) / &
                                  ps(o_js:o_je,o_is:o_ie))
            end do
            istatus = nf90_put_var(ncche(n), ichevar(18), &
                                   dumio, istart, icount)
            call check_ok(__FILE__,__LINE__, &
                         'Error writing cbdydiag at '//ctime,'CHE FILE ERROR')

            do k = 1 , kz
              dumio(:,:,k) = real(cseddpdiag(o_js:o_je,o_is:o_ie,k,n) / &
                                  ps(o_js:o_je,o_is:o_ie))
            end do
            istatus = nf90_put_var(ncche(n), ichevar(19), &
                                   dumio, istart, icount)
            call check_ok(__FILE__,__LINE__, &
                         'Error writing cbdydiag at '//ctime,'CHE FILE ERROR')
            if ( ibltyp == 2 .or. ibltyp == 99 ) then
              !TAO: Re-set istart/icount, since ccuwdiag is a 2D variable
              !and the previous istart/icount setting was for 3D variables
              istart(4) = 1
              istart(3) = icherec
              istart(2) = 1
              istart(1) = 1
              icount(4) = 1
              icount(3) = 1
              icount(2) = o_ni
              icount(1) = o_nj
              dumio(:,:,1) = real(ccuwdiag(o_js:o_je,o_is:o_ie,n))
              istatus = nf90_put_var(ncche(n), ichevar(20), &
                                     dumio(:,:,1), istart, icount)
              call check_ok(__FILE__,__LINE__, &
                           'Error writing cbdydiag at '//ctime,'CHE FILE ERROR')
            end if
          end if 

          !closing
          istatus = nf90_sync(ncche(n))
          call check_ok(__FILE__,__LINE__, &
                        'Error sync at '//ctime, 'CHE FILE ERROR')

        else if ( n == noutf ) then 

          istart(4) = icherec
          istart(3) = 1
          istart(2) = 1
          istart(1) = 1
          icount(5) = 1
          icount(4) = 1
          icount(3) = o_nz
          icount(2) = o_ni
          icount(1) = o_nj

          !*** extinction
          do k = 1 , kz
            dumio(:,:,k) = real(ext(o_js:o_je,o_is:o_ie,kz-k+1))
          end do
          istatus = nf90_put_var(ncche(n), ioptvar(3), &
                                 dumio, istart(1:4), icount(1:4))
          call check_ok(__FILE__,__LINE__,'Error writing EXT at '//ctime,&
                         'CHE FILE ERROR')

          !*** SSAE
          do k = 1 , kz
            dumio(:,:,k) = real(ssa(o_js:o_je,o_is:o_ie,kz-k+1))
          end do
          istatus = nf90_put_var(ncche(n), ioptvar(4), &
                                 dumio, istart(1:4), icount(1:4))
          call check_ok(__FILE__,__LINE__,'Error writing SSAE at '//ctime,&
                        'CHE FILE ERROR')

          !*** ASP
          do k = 1 , kz
            dumio(:,:,k) = real(asp(o_js:o_je,o_is:o_ie,kz-k+1))
          end do
          istatus = nf90_put_var(ncche(n), ioptvar(5), &
                                 dumio, istart(1:4), icount(1:4))
          call check_ok(__FILE__,__LINE__,'Error writing ASP at '//ctime,&
                        'CHE FILE ERROR')

          !*** output 2-D Aerosl Radiative forcing : NB these variables
          ! are already accumulated and averagesd in aerout

          istart(3) = icherec
          istart(2) = 1
          istart(1) = 1
          icount(3) = 1
          icount(2) = o_ni
          icount(1) = o_nj
          dumio(:,:,1) = real(tarf(o_js:o_je,o_is:o_ie))
          istatus = nf90_put_var(ncche(n), ioptvar(6) , &
                                 dumio(:,:,1), istart(1:3), icount(1:3))
          call check_ok(__FILE__,__LINE__,'Error writing aertarf at '//ctime, &
                       'OPT FILE ERROR')
          dumio(:,:,1) = real(ssrf(o_js:o_je,o_is:o_ie))
          istatus = nf90_put_var(ncche(n),  ioptvar(7), &
                                 dumio(:,:,1), istart(1:3), icount(1:3))
          call check_ok(__FILE__,__LINE__,'Error writing aersrrf at '//ctime, &
                       'OPT FILE ERROR')
          dumio(:,:,1) = real(talwrf(o_js:o_je,o_is:o_ie))
          istatus = nf90_put_var(ncche(n),ioptvar(8), &
                                 dumio(:,:,1), istart(1:3), icount(1:3))
          call check_ok(__FILE__,__LINE__, &
                     'Error writing aertalwrf at '//ctime,'OPT FILE ERROR')
          dumio(:,:,1) = real(srlwrf(o_js:o_je,o_is:o_ie))
          istatus = nf90_put_var(ncche(n), ioptvar(9), &
                                 dumio(:,:,1), istart(1:3), icount(1:3))
          call check_ok(__FILE__,__LINE__, &
                     'Error writing aersrlwrf at '//ctime,'OPT FILE ERROR')
          dumio(:,:,1) = real(aod(o_js:o_je,o_is:o_ie))
          istatus = nf90_put_var(ncche(n), ioptvar(10), &
                                 dumio(:,:,1), istart(1:3), icount(1:3))
          call check_ok(__FILE__,__LINE__, &
                     'Error writing aod at '//ctime,'OPT FILE ERROR')
 
          istatus = nf90_sync(ncche(n))
          call check_ok(__FILE__,__LINE__, &
                     'Error sync at '//ctime, 'OPT FILE ERROR')
        end if 
      end do     !main species looop
      icherec = icherec + 1
    end subroutine writerec_che2

!===========================================

    integer function chbc_search(idate)
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
