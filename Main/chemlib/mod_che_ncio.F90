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

module mod_che_ncio

  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_date
  use mod_nchelper
  use mod_dynparam
  use mod_mpmessage
  use mod_che_indices
  use mod_che_common
  use mod_runparams
  use mod_domain
  use netcdf

  implicit none

  private

  public :: read_texture , read_emission , read_miner , recc , reccbb
  public :: read_bioburn_emission
  public :: init_mod_che_ncio
  public :: open_chbc , close_chbc , chbc_search , read_chbc , read_bionem
  public :: read_dust_param
  public :: chbc_ivar , n_chbcvar , n_aebcvar, chbcname, aeaero, aedu12

  integer(ik4) :: istatus
  integer(ik4) :: recc , reccbb

  integer(ik4) , parameter :: n_chevar = 20
  integer(ik4) , parameter :: n_oxbcvar = 5
  integer(ik4) , parameter :: n_optvar = 10
!ashalaby  integer(ik4) , parameter :: n_chbcvar = 25
  integer(ik4) , parameter :: n_chbcvar = 33
  integer(ik4) :: n_aebcvar
  integer(ik4) :: ichin , iaein , ioxin

  character(len=8) , dimension(n_chbcvar) :: chbcname
  character(len=8) , dimension(n_oxbcvar) :: oxbcname
  character(len=8) , target , dimension(4) :: aedust
  character(len=8) , target , dimension(12) :: aedu12
  character(len=8) , target , dimension(2) :: aesslt
  character(len=8) , target , dimension(6) :: aecarb
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
  real(rkx) , dimension(:,:) , pointer :: rspace2
  real(rkx) , dimension(:,:) , pointer :: rspace2_loc
  real(rkx) , dimension(:,:,:) , pointer :: rspace3
  real(rkx) , dimension(:,:,:) , pointer :: rspace3_loc

  character(256) :: icbcname

  data ichin   /-1/
  data iaein   /-1/
  data ibcrec  / 1/
  data ibcnrec / 0/
  data ioxin   /-1/

!  data chbcname /'O3      ','NO      ','NO2     ','HNO3    ', &
!                 'N2O5    ','H2O2    ','CH4     ','CO      ', &
!                 'CH2O    ','CH3OH   ','C2H5OH  ','C2H4    ', &
!                 'C2H6    ','CH3CHO  ','CH3COCH3','BIGENE  ', &
!                 'BIGALK  ','C3H6    ','C3H8    ','ISOP    ', &
!                 'TOLUENE ','PAN     ','SO2     ','SO4     ', &
!                 'DMS     '/
  data chbcname / 'O3      ','NO      ','NO2     ','HNO3    ', &
                  'HNO4    ','N2O5    ','H2O2    ','CH4     ', &
                  'CO      ','SO2     ','H2SO4   ','DMS     ', &
                  'PAR     ','C2H6    ','ETH     ','OLET    ', &
                  'OLEI    ','TOL     ','XYL     ','ISOP    ', &
                  'CRES    ','OPEN    ','ISOPN   ','ISOPRD  ', &
                  'ONIT    ','MGLY    ','AONE    ','PAN     ', &
                  'CH3OOH  ','ETHOOH  ','ALD2    ','HCHO    ', &
                  'CH3OH   '/
  data oxbcname /'OH      ','HO2     ','O3      ', 'NO3    ','H2O2   ' /
  data aedust / 'DUST01' , 'DUST02' , 'DUST03', 'DUST04' /
  data aedu12 / 'DUST01', 'DUST02', 'DUST03', 'DUST04',  &
                'DUST05', 'DUST06', 'DUST07', 'DUST08',  &
                'DUST09', 'DUST10', 'DUST11', 'DUST12' /
  data aesslt / 'SSLT01' , 'SSLT02' /
  data aecarb / 'BC_HB' , 'BC_HL' , 'OC_HB' , 'OC_HL' , 'SM1' , 'SM2' /
  data aesulf / 'SO2' , 'SO4' /
  data aesuca / 'BC_HB' , 'BC_HL' , 'OC_HB' , 'OC_HL' , 'SO2' , 'SO4' /
  data aeaero / 'BC_HB' , 'BC_HL' , 'OC_HB' , 'OC_HL' , 'SO2' , 'SO4' , &
                'SSLT01' , 'SSLT02', 'DUST01', 'DUST02', 'DUST03' , 'DUST04' /

  contains

    subroutine init_mod_che_ncio(chemsymtype)
      implicit none
      character(len=8) , intent(in) :: chemsymtype

      n_aebcvar = 0
      select case ( chemsymtype )
        case ( 'DUST' )
          n_aebcvar = 4
          aebcname => aedust
        case ( 'DU12' )
          n_aebcvar = 12
          aebcname => aedu12
        case ( 'SSLT' )
          n_aebcvar = 2
          aebcname => aesslt
        case ( 'CARB' )
          n_aebcvar = 4
          aebcname => aecarb
        case ( 'SULF' )
          n_aebcvar = 2
          aebcname => aesulf
        case ( 'SUCA' , 'SUCE' )
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

    subroutine read_texture(nats,rtex)
      implicit none
      integer(ik4) , intent(in) :: nats
      real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: rtex
      integer(ik4) :: idmin
      integer(ik4) , dimension(3) :: istart , icount
      character(len=256) :: dname
      real(rkx) , pointer , dimension(:,:,:) ::  rspace

      dname = trim(dirter)//pthsep//trim(domname)//'_DOMAIN000.nc'

      if ( do_parallel_netcdf_in ) then
        call openfile_withname(dname,idmin)
        istart(1) = jde1
        istart(2) = ide1
        istart(3) = 1
        icount(1) = jde2-jde1+1
        icount(2) = ide2-ide1+1
        icount(3) = nats
        allocate(rspace(jde1:jde2,ide1:ide2,nats))
        call read_var3d_static(idmin,'texture_fraction',rspace, &
                               istart=istart,icount=icount)
        rtex(jci1:jci2,ici1:ici2,:) = &
          max(rspace(jci1:jci2,ici1:ici2,:)*0.01_rkx,d_zero)
        call closefile(idmin)
        deallocate(rspace)
      else
        if ( myid == iocpu ) then
          call openfile_withname(dname,idmin)
          istart(1) = 1
          istart(2) = 1
          istart(3) = 1
          icount(1) = jx
          icount(2) = iy
          icount(3) = nats
          allocate(rspace(jx,iy,nats))
          call read_var3d_static(idmin,'texture_fraction',rspace, &
            istart=istart,icount=icount)
          rspace = max(rspace*0.01_rkx,d_zero)
          call grid_distribute(rspace,rtex,jci1,jci2,ici1,ici2,1,nats)
          call closefile(idmin)
          deallocate(rspace)
        else
          call grid_distribute(rspace,rtex,jci1,jci2,ici1,ici2,1,nats)
        end if
      end if
    end subroutine read_texture

   subroutine read_dust_param(erodfc, aez0)
!read dust emission relevant parameters
!for now : erod_dsfc = source function / erodibility mask
!e.g. see  Zender et al., Laurent et al.
!place holder for other relevant geographical data afecting dust ( e.g. non erodibe zo)
      implicit none

      real(rkx) , pointer , dimension(:,:) , intent(inout) :: erodfc
      real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: aez0
      integer(ik4) :: idmin
      integer(ik4) , dimension(2) :: istart , icount
      integer(ik4) , dimension(3) :: istart3 , icount3
      character(len=256) :: dname
      real(rkx) , pointer , dimension(:,:) ::  rspace
      real(rkx), pointer , dimension(:,:,:) :: rspace3

      dname = trim(dirter)//pthsep//trim(domname)//'_DUSTPARAM.nc'

      if ( do_parallel_netcdf_in ) then
        call openfile_withname(dname,idmin)
        istart(1) = jde1
        istart(2) = ide1
        icount(1) = jde2-jde1+1
        icount(2) = ide2-ide1+1
        allocate(rspace(jde1:jde2,ide1:ide2))
        call read_var2d_static(idmin,'erodfc',rspace, &
          istart=istart,icount=icount)
        erodfc(jci1:jci2,ici1:ici2) = &
          max(rspace(jci1:jci2,ici1:ici2),d_zero)
        call closefile(idmin)
        deallocate(rspace)
      else
        if ( myid == iocpu ) then
          call openfile_withname(dname,idmin)
          istart(1) = 1
          istart(2) = 1
          icount(1) = jx
          icount(2) = iy
          allocate(rspace(jx,iy))
          call read_var2d_static(idmin,'erodfc',rspace, &
            istart=istart,icount=icount)
          rspace = max(rspace,d_zero)
          call grid_distribute(rspace,erodfc,jci1,jci2,ici1,ici2)

          istart3(1) = 1
          istart3(2) = 1
          istart3(3) = 1
          icount3(1) = jx
          icount3(2) = iy
          icount3(3) = 12
          allocate(rspace3(jx,iy,12))
          call read_var3d_static(idmin,'z0',rspace3, &
            istart=istart3,icount=icount3)
          rspace3 = max(rspace3,d_zero)
          call grid_distribute(rspace3,aez0,jci1,jci2,ici1,ici2,1,12)

          call closefile(idmin)
          deallocate(rspace)
          deallocate(rspace3)
        else
          call grid_distribute(rspace,erodfc,jci1,jci2,ici1,ici2)
          call grid_distribute(rspace3,aez0,jci1,jci2,ici1,ici2,1,12)
        end if
      end if
    end subroutine read_dust_param

    subroutine read_bionem(nfert,nmanure,soilph)
      implicit none

      real(rkx) , pointer , dimension(:,:) , intent(inout) :: nfert
      real(rkx) , pointer , dimension(:,:) , intent(inout) :: nmanure
      real(rkx) , pointer , dimension(:,:) , intent(inout) :: soilph
      integer(ik4) :: idmin
      integer(ik4) , dimension(2) :: istart , icount
      character(len=256) :: dname
      real(rkx) , pointer , dimension(:,:) ::  rspace

      dname = trim(dirter)//pthsep//trim(domname)//'_BIONEM.nc'

      if ( do_parallel_netcdf_in ) then
        call openfile_withname(dname,idmin)
        istart(1) = jde1
        istart(2) = ide1
        icount(1) = jde2-jde1+1
        icount(2) = ide2-ide1+1
        allocate(rspace(jde1:jde2,ide1:ide2))
        call read_var2d_static(idmin,'nfert',rspace, &
          istart=istart,icount=icount)
        nfert(jci1:jci2,ici1:ici2) = &
          max(rspace(jci1:jci2,ici1:ici2),d_zero)

        call read_var2d_static(idmin,'nmanure',rspace, &
          istart=istart,icount=icount)
           nmanure(jci1:jci2,ici1:ici2) = &
          max(rspace(jci1:jci2,ici1:ici2),d_zero)

         call read_var2d_static(idmin,'SoilpH',rspace, &
          istart=istart,icount=icount)
           soilph(jci1:jci2,ici1:ici2) = &
          max(rspace(jci1:jci2,ici1:ici2),d_zero)

        call closefile(idmin)
        deallocate(rspace)
      else
        if ( myid == iocpu ) then
          call openfile_withname(dname,idmin)
          istart(1) = 1
          istart(2) = 1
          icount(1) = jx
          icount(2) = iy
          allocate(rspace(jx,iy))
          call read_var2d_static(idmin,'nfert',rspace, &
            istart=istart,icount=icount)
          rspace = max(rspace,d_zero)
          call grid_distribute(rspace,nfert,jci1,jci2,ici1,ici2)

          call read_var2d_static(idmin,'nmanure',rspace, &
            istart=istart,icount=icount)
          rspace = max(rspace,d_zero)
          call grid_distribute(rspace,nmanure,jci1,jci2,ici1,ici2)

           call read_var2d_static(idmin,'SoilpH',rspace, &
            istart=istart,icount=icount)
          rspace = max(rspace,d_zero)
          call grid_distribute(rspace,soilph,jci1,jci2,ici1,ici2)
          call closefile(idmin)
          deallocate(rspace)
        else
          call grid_distribute(rspace,nfert,jci1,jci2,ici1,ici2)
          call grid_distribute(rspace,nmanure,jci1,jci2,ici1,ici2)
          call grid_distribute(rspace,soilph,jci1,jci2,ici1,ici2)
        end if
      end if
    end subroutine read_bionem

    subroutine read_miner(nmine,cminer,sminer)
      implicit none
      integer(ik4), intent(in) :: nmine
      real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: cminer , sminer

      integer(ik4) :: idmin
      integer(ik4) , dimension(2) :: istart , icount
      character(len=256) :: dname
      real(rkx) , pointer , dimension(:,:) ::  rspace
      real(rkx) , pointer , dimension(:,:,:) :: tmp

      dname = trim(dirter)//pthsep//trim(domname)//'_MINER.nc'

      if ( do_parallel_netcdf_in ) then
        ! THIS PART IS NOT TESTED
        call openfile_withname(dname,idmin)
        istart(1) = jde1
        istart(2) = ide1
        icount(1) = jde2-jde1+1
        icount(2) = ide2-ide1+1
        call fatal(__FILE__,__LINE__, &
                       'NOT IMPLEMENTED PARALLEL')
!        allocate(rspace(icount(1),icount(2)))
!        call read_var2d_static(idmin,'CIRON',rspace, &
!                               istart=istart,icount=icount)
!        cminer(jci1:jci2,ici1:ici2,isciron) = &
!            max(rspace(jci1:jci2,ici1:ici2),d_zero)*0.01_rkx
!        call read_var2d_static(idmin,'CHMT',rspace, &
!                               istart=istart,icount=icount)
!        cminer(jci1:jci2,ici1:ici2,ischmt) = &
!               max(rspace(jci1:jci2,ici1:ici2),d_zero)*0.01_rkx
!        call closefile(idmin)
      else
        if ( myid == iocpu ) then
          call openfile_withname(dname,idmin)
          istart(1) = 1
          istart(2) = 1
          icount(1) = jx
          icount(2) = iy
          allocate(rspace(icount(1),icount(2)))
          allocate(tmp(size(rspace,1),size(rspace,2),nmine))
          ! read all miner variable in file ( percent converted to fraction)
          tmp = 0_rkx   !*****************jj**************
          call read_var2d_static(idmin,'CIRON',rspace, &
                                 istart=istart,icount=icount)
          rspace = max(rspace,d_zero)
          ! replace by the miner index
          tmp(:,:,iiron) = rspace(:,:)*0.01_rkx
          call read_var2d_static(idmin,'CHMT',rspace, &
                                 istart=istart,icount=icount)
          rspace = max(rspace,d_zero)
          tmp(:,:,ihmt) = rspace(:,:)*0.01_rkx

          call read_var2d_static(idmin,'CCAL',rspace, &
                                 istart=istart,icount=icount)
          rspace = max(rspace,d_zero)
          tmp(:,:,icalc) = rspace(:,:)*0.01_rkx
          !********************************added --jj*******************
          call read_var2d_static(idmin,'CGTH',rspace, &
                                 istart=istart,icount=icount)
          rspace = max(rspace,d_zero)
          tmp(:,:,igth) = rspace(:,:)*0.01_rkx
          call read_var2d_static(idmin,'CCHL',rspace, &
                                 istart=istart,icount=icount)
          rspace = max(rspace,d_zero)
          tmp(:,:,ichl) = rspace(:,:)*0.01_rkx
          call read_var2d_static(idmin,'CFLD',rspace, &
                                 istart=istart,icount=icount)
          rspace = max(rspace,d_zero)
          tmp(:,:,ifld) = rspace(:,:)*0.01_rkx
          call read_var2d_static(idmin,'CQTZ',rspace, &
                                 istart=istart,icount=icount)
          rspace = max(rspace,d_zero)
          tmp(:,:,iqtz) = rspace(:,:)*0.01_rkx
          call read_var2d_static(idmin,'CSMEC',rspace, &
                                 istart=istart,icount=icount)
          rspace = max(rspace,d_zero)
          tmp(:,:,ismec) = rspace(:,:)*0.01_rkx
          call read_var2d_static(idmin,'CVMC',rspace, &
                                 istart=istart,icount=icount)
          rspace = max(rspace,d_zero)
          tmp(:,:,ivmc) = rspace(:,:)*0.01_rkx
          call read_var2d_static(idmin,'CKAL',rspace, &
                                 istart=istart,icount=icount)
          rspace = max(rspace,d_zero)
          tmp(:,:,ikal) = rspace(:,:)*0.01_rkx
          call read_var2d_static(idmin,'CILT',rspace, &
                                 istart=istart,icount=icount)
          rspace = max(rspace,d_zero)
          tmp(:,:,iilt) = rspace(:,:)*0.01_rkx
          !***********************************************************
          call grid_distribute(tmp,cminer,jci1,jci2,ici1,ici2,1,nmine)
          !! now the silt fraction
          tmp = 0_rkx
          call read_var2d_static(idmin,'SCAL',rspace, &
                                 istart=istart,icount=icount)
          rspace = max(rspace,d_zero)
          tmp(:,:,icalc) = rspace(:,:)*0.01_rkx
          !****************8added --jjoshi****************************
          call read_var2d_static(idmin,'SGTH',rspace, &
                                 istart=istart,icount=icount)
          rspace = max(rspace,d_zero)
          tmp(:,:,igth) = rspace(:,:)*0.01_rkx
          call read_var2d_static(idmin,'SCHL',rspace, &
                                 istart=istart,icount=icount)
          rspace = max(rspace,d_zero)
          tmp(:,:,ichl) = rspace(:,:)*0.01_rkx

          call read_var2d_static(idmin,'SFLD',rspace, &
                                 istart=istart,icount=icount)
          rspace = max(rspace,d_zero)
          tmp(:,:,ifld) = rspace(:,:)*0.01_rkx
          call read_var2d_static(idmin,'SQTZ',rspace, &
                                 istart=istart,icount=icount)
          rspace = max(rspace,d_zero)
          tmp(:,:,iqtz) = rspace(:,:)*0.01_rkx
          call read_var2d_static(idmin,'SGYP',rspace, &
                                 istart=istart,icount=icount)
          rspace = max(rspace,d_zero)
          tmp(:,:,igyp) = rspace(:,:)*0.01_rkx
          call read_var2d_static(idmin,'SMICA',rspace, &
                                 istart=istart,icount=icount)
          rspace = max(rspace,d_zero)
          tmp(:,:,imica) = rspace(:,:)*0.01_rkx
          !************************************************************
          call grid_distribute(tmp,sminer,jci1,jci2,ici1,ici2,1,nmine)
          call closefile(idmin)
        else
          call grid_distribute(tmp,cminer,jci1,jci2,ici1,ici2,1,nmine)
          call grid_distribute(tmp,sminer,jci1,jci2,ici1,ici2,1,nmine)
        end if
      end if
    end subroutine read_miner

    subroutine read_emission(ifreq,lyear,lmonth,lday,lhour,echemsrc)
      implicit none
      integer(ik4) , intent(in) :: lyear , lmonth , lday , lhour
      integer(ik4) , intent(out) :: ifreq
      real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: echemsrc
      character(256) :: aername
      integer(ik4) :: n,ncid , itvar, idimid, chmnrec,sdim
      character(64) ::chemi_timeunits , chemi_timecal
      real(rkx) , dimension(:) , allocatable :: emtimeval
      integer(ik4) , dimension(4) :: istart , icount
      integer(ik4) :: year , month , day , hour
      type(rcm_time_and_date) :: tchdate

      ! FAB: remember for now, we have 1 emission file containing all monthly
      ! emission for the whole simulation period
      ! change that in the future.,

      aername = trim(dirglob)//pthsep//trim(domname)//'_CHEMISS.nc'
      if ( myid == italk ) then
        write(stdout,*) 'Opening ch. emission file ', trim(aername)
      end if

      if ( do_parallel_netcdf_in .or. myid == iocpu ) then
        call openfile_withname(aername,ncid)
        istatus = nf90_inq_dimid(ncid, 'time', idimid)
        call check_ok(__FILE__,__LINE__, &
                      'Dimension time miss', 'CHEMI FILE')
        istatus = nf90_inquire_dimension(ncid, idimid, len=chmnrec)
        call check_ok(__FILE__,__LINE__, &
                'Dimension time read error', 'CHEMI FILE')
        allocate (emtimeval(chmnrec))
        istatus = nf90_inq_varid(ncid, 'time', itvar)
        call check_ok(__FILE__,__LINE__, &
                      'variable time miss', 'CHEMISS FILE')
        istatus = nf90_get_att(ncid, itvar, 'units', chemi_timeunits)
        call check_ok(__FILE__,__LINE__, &
                      'variable time units miss', &
                      'CHEMISS FILE')
        if ( chemi_timeunits(1:5) == 'month' ) then
          ifreq = ifrqmon
        else if ( chemi_timeunits(1:3) == 'day' ) then
          ifreq = ifrqday
        else if ( chemi_timeunits(1:4) == 'hour' ) then
          ifreq = ifrqhrs
        else
          call fatal(__FILE__,__LINE__, &
                        'NO CODED FREQUENCY IN CHEMISS FILE')
        end if
        istatus = nf90_get_att(ncid, itvar, 'calendar', chemi_timecal)
        if ( istatus /= nf90_noerr ) then
          chemi_timecal = 'gregorian'
        end if
        istatus = nf90_get_var(ncid, itvar, emtimeval)
        call check_ok(__FILE__,__LINE__, &
                      'variable time read error', 'ICBC FILE')
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
          write(stderr,*) 'chem emission : time record not found emission file'
          call fatal(__FILE__,__LINE__, &
                        'IO ERROR in CHEM EMISSION')
        else
          if ( myid == italk ) then
            write(stdout,*) 'CHE_EMISS: Reading record ', recc
          end if
        end if

        !*** intialized in start_chem
        !*** Advice record counter
        istart = 0
        icount = 0

        if ( do_parallel_netcdf_in ) then
          istart(1) = jde1
          istart(2) = ide1
          icount(1) = jde2-jde1+1
          icount(2) = ide2-ide1+1
          allocate(rspace2(jde1:jde2,ide1:ide2))
        else
          istart(1) = 1
          istart(2) = 1
          icount(1) = jx
          icount(2) = iy
          allocate(rspace2_loc(jci1:jci2,ici1:ici2))
          allocate(rspace2(jx,iy))
        end if
        istatus = nf90_inq_dimid(ncid, 'lev', idimid)
        if ( istatus /= nf90_noerr ) then
          ! no lev diemsion in emission variables
          istart(3) = recc
          icount(3) = 1
          sdim = 3
        else
          istart(3) = 1
          istart(4) = recc
          icount(3) = 1
          icount(4) = 1
          sdim=4
        end if
      else
        allocate(rspace2_loc(jci1:jci2,ici1:ici2))
      end if

! ALD2
      if ( iald2 /= 0 ) then
        call rvar(ncid,istart,icount,iald2,echemsrc,'ALD2_flux',.false.,sdim)
      end if
! AONE
      if ( iaone /= 0 ) then
        call rvar(ncid,istart,icount,iaone,echemsrc,'AONE_flux',.false.,sdim)
      end if
!BCHB
      if ( ibchb /= 0 ) then
        call rvar(ncid,istart,icount,ibchb,echemsrc,'BC_flux',.false.,sdim)
      end if

!C2H6
      if ( ic2h6 /= 0 ) then
        call rvar(ncid,istart,icount,ic2h6,echemsrc,'C2H6_flux',.false.,sdim)
      end if
!CH3OH
      if ( ich3oh /= 0 ) then
        call rvar(ncid,istart,icount,ich3oh,echemsrc,'CH3OH_flux',.false.,sdim)
      end if
!CH4
      if ( ich4 /= 0 ) then
        call rvar(ncid,istart,icount,ich4,echemsrc,'CH4_flux',.false.,sdim)
      end if
!CO
      if ( ico /= 0 ) then
        call rvar(ncid,istart,icount,ico,echemsrc,'CO_flux',.false.,sdim)
      end if
!ETH
      if ( ieth /= 0 ) then
        call rvar(ncid,istart,icount,ieth,echemsrc,'ETH_flux',.false.,sdim)
      end if
!HCHO
      if ( ihcho /= 0 ) then
        call rvar(ncid,istart,icount,ihcho,echemsrc,'HCHO_flux',.false.,sdim)
      end if
!NH3
      if ( inh3 /= 0 ) then
        call rvar(ncid,istart,icount,inh3,echemsrc,'NH3_flux',.false.,sdim)
      end if
!NO
      if ( ino /= 0 ) then
        call rvar(ncid,istart,icount,ino,echemsrc,'NOx_flux',.false.,sdim)
      end if
!OCHB
      if ( iochb /= 0 ) then
        call rvar(ncid,istart,icount,iochb,echemsrc,'OC_flux',.false.,sdim)
      end if
!SMOKE
!      if ( ism1 /= 0 ) then
!        call rvar(ncid,istart,icount,ism1,echemsrc,'SM_flux',.false.,sdim)
!      end if
!OLET
      if ( iolet /= 0 ) then
        call rvar(ncid,istart,icount,iolet,echemsrc,'OLET_flux',.false.,sdim)
      end if
!OLEI
      if ( iolei /= 0 ) then
        call rvar(ncid,istart,icount,iolei,echemsrc,'OLEI_flux',.false.,sdim)
      end if
!PAR
      if ( ipar/= 0 ) then
        call rvar(ncid,istart,icount,ipar,echemsrc,'PAR_flux',.false.,sdim)
      end if
!RCOOH
      if ( ircooh/= 0 ) then
        call rvar(ncid,istart,icount,ircooh,echemsrc,'RCOOH_flux',.false.,sdim)
      end if
!SO2
      if ( iso2/= 0 ) then
        call rvar(ncid,istart,icount,iso2,echemsrc,'SO2_flux',.false.,sdim)
      end if
!TOL
      if ( itol/= 0 ) then
        call rvar(ncid,istart,icount,itol,echemsrc,'TOL_flux',.false.,sdim)
      end if
!XYL
      if ( ixyl/= 0 ) then
        call rvar(ncid,istart,icount,ixyl,echemsrc,'XYL_flux',.false.,sdim)
      end if
!ISOP
      if ( iisop/= 0 ) then
        call rvar(ncid,istart,icount,iisop,echemsrc,'ISOP_BIO_flux',.false.,sdim)
      end if

      ! NO2 emission
      if ( ino2 /= 0 ) then
!       call rvar(ncid,istart,icount,ino2,echemsrc, &
!                 'NO_flux',.false.,sdim)
!        echemsrc(:,:,ino2) = 0.1_rkx * echemsrc(:,:,ino)
!        echemsrc(:,:,ino)  = 0.9_rkx * echemsrc(:,:,ino)
      end if

      if (ipollen /=0 ) then
        call  rvar(ncid,istart,icount,ipollen,echemsrc, &
                  'POLLEN',.false.,sdim)
      end if

      where (echemsrc(:,:,:) < d_zero ) echemsrc(:,:,:) = d_zero

      if ( do_parallel_netcdf_in ) then
        call closefile(ncid)
        deallocate (emtimeval)
        deallocate(rspace2)
      else
        if ( myid == iocpu ) then
          call closefile(ncid)
          deallocate (emtimeval)
          deallocate(rspace2)
          deallocate(rspace2_loc)
        else
          deallocate(rspace2_loc)
        end if
        call bcast(ifreq)
      end if
      if ( myid == italk ) then
        write(stdout,*) 'Emission read in success.'
      end if
    end subroutine read_emission

    subroutine read_bioburn_emission(ifreq,lyear,lmonth,lday,lhour,echemsrc)
      implicit none
      integer(ik4) , intent(in) :: lyear , lmonth , lday , lhour
      integer(ik4) , intent(out) :: ifreq
      real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: echemsrc
      character(256) :: aername
      integer(ik4) :: n,ncid , itvar, idimid, chmnrec,sdim
      character(64) ::chemi_timeunits , chemi_timecal
      real(rkx) , dimension(:) , allocatable :: emtimeval
      integer(ik4) , dimension(4) :: istart , icount
      integer(ik4) :: year , month , day , hour
      type(rcm_time_and_date) :: tchdate


      aername = trim(dirglob)//pthsep//trim(domname)//'_CHEMISS_BB.nc'
      if ( myid == italk ) then
        write(stdout,*) 'Opening ch. biomass burning emission file ', trim(aername)
      end if

      if ( do_parallel_netcdf_in .or. myid == iocpu ) then
        call openfile_withname(aername,ncid)
        istatus = nf90_inq_dimid(ncid, 'time', idimid)
        call check_ok(__FILE__,__LINE__, &
                      'Dimension time miss', 'CHEMI BB FILE')
        istatus = nf90_inquire_dimension(ncid, idimid, len=chmnrec)
        call check_ok(__FILE__,__LINE__, &
                'Dimension time read error', 'CHEMI BB FILE')
        allocate (emtimeval(chmnrec))
        istatus = nf90_inq_varid(ncid, 'time', itvar)
        call check_ok(__FILE__,__LINE__, &
                      'variable time miss', 'CHEMISS BB FILE')
        istatus = nf90_get_att(ncid, itvar, 'units', chemi_timeunits)
        call check_ok(__FILE__,__LINE__, &
                      'variable time units miss', &
                      'CHEMISS BB FILE')
        if ( chemi_timeunits(1:6) == 'months' ) then
          ifreq = ifrqmon
        else if ( chemi_timeunits(1:4) == 'days' ) then
          ifreq = ifrqday
        else if ( chemi_timeunits(1:5) == 'hours' ) then
          ifreq = ifrqhrs
        else
          call fatal(__FILE__,__LINE__, &
                        'NO CODED FREQUENCY IN CHEMISS FILE')
        end if
        istatus = nf90_get_att(ncid, itvar, 'calendar', chemi_timecal)
        if ( istatus /= nf90_noerr ) then
          chemi_timecal = 'gregorian'
        end if
        istatus = nf90_get_var(ncid, itvar, emtimeval)
        call check_ok(__FILE__,__LINE__, &
                      'variable time read error', 'chemiss FILE')

        reccbb = 0
        looprec: &
        do n = 1 , chmnrec
          tchdate = timeval2date(emtimeval(n),chemi_timeunits,chemi_timecal)
          call split_idate(tchdate,year,month,day,hour)
          select case (ifreq)
            case (ifrqmon)
              if ( year == lyear .and. month == lmonth ) then
                reccbb = n
                exit looprec
              end if
            case (ifrqday)
              if ( year == lyear .and. month == lmonth .and. day == lday ) then
                reccbb = n
                exit looprec
              end if
            case (ifrqhrs)
              if ( year == lyear .and. month == lmonth .and. &
                    day == lday  .and. hour  == lhour ) then
                reccbb = n
                exit looprec
              end if
          end select
        end do looprec

        if ( reccbb == 0 ) then
          write(stderr,*) &
             'chem BB emission : time record not found emission file'
          call fatal(__FILE__,__LINE__, &
                        'IO ERROR in CHEM BB EMISSION')
        else
          if ( myid == italk ) then
            write(stdout,*) 'CHE_EMISS BB: Reading record ', reccbb
          end if
        end if

        !*** intialized in start_chem
        !*** Advice record counter
        istart = 0
        icount = 0

        if ( do_parallel_netcdf_in ) then
          istart(1) = jde1
          istart(2) = ide1
          icount(1) = jde2-jde1+1
          icount(2) = ide2-ide1+1
          allocate(rspace2(jde1:jde2,ide1:ide2))
        else
          istart(1) = 1
          istart(2) = 1
          icount(1) = jx
          icount(2) = iy
          allocate(rspace2_loc(jci1:jci2,ici1:ici2))
          allocate(rspace2(jx,iy))
        end if
        istatus = nf90_inq_dimid(ncid, 'lev', idimid)
        if ( istatus /= nf90_noerr ) then
          ! no lev diemsion in emission variables
          istart(3) = reccbb
          icount(3) = 1
          sdim = 3
        else
          istart(3) = 1
          istart(4) = reccbb
          icount(3) = 1
          icount(4) = 1
          sdim=4
        end if
      else
        allocate(rspace2_loc(jci1:jci2,ici1:ici2))
      end if

! ALD2
      if ( iald2 /= 0 ) then
        call rvar(ncid,istart,icount,iald2,echemsrc,'ALD2_flux',.false.,sdim)
      end if
! AONE
      if ( iaone /= 0 ) then
        call rvar(ncid,istart,icount,iaone,echemsrc,'AONE_flux',.false.,sdim)
      end if
!C2H6
      if ( ic2h6 /= 0 ) then
        call rvar(ncid,istart,icount,ic2h6,echemsrc,'C2H6_flux',.false.,sdim)
      end if
!CH3OH
      if ( ich3oh /= 0 ) then
        call rvar(ncid,istart,icount,ich3oh,echemsrc,'CH3OH_flux',.false.,sdim)
      end if
!CH4
      if ( ich4 /= 0 ) then
        call rvar(ncid,istart,icount,ich4,echemsrc,'CH4_flux',.false.,sdim)
      end if
!CO
      if ( ico /= 0 ) then
        call rvar(ncid,istart,icount,ico,echemsrc,'CO_flux',.false.,sdim)
      end if
!ETH
      if ( ieth /= 0 ) then
        call rvar(ncid,istart,icount,ieth,echemsrc,'ETH_flux',.false.,sdim)
      end if
!HCHO
      if ( ihcho /= 0 ) then
        call rvar(ncid,istart,icount,ihcho,echemsrc,'HCHO_flux',.false.,sdim)
      end if
!NH3
      if ( inh3 /= 0 ) then
        call rvar(ncid,istart,icount,inh3,echemsrc,'NH3_flux',.false.,sdim)
      end if
!NO
      if ( ino /= 0 ) then
        call rvar(ncid,istart,icount,ino,echemsrc,'NOx_flux',.false.,sdim)
      end if
!SMOKE
!      if ( ism1 /= 0 ) then
!        call rvar(ncid,istart,icount,ism1,echemsrc,'SM_flux',.false.,sdim)
!      end if
!OLET
      if ( iolet /= 0 ) then
        call rvar(ncid,istart,icount,iolet,echemsrc,'OLET_flux',.false.,sdim)
      end if
!OLEI
      if ( iolei /= 0 ) then
        call rvar(ncid,istart,icount,iolei,echemsrc,'OLEI_flux',.false.,sdim)
      end if
!PAR
      if ( ipar/= 0 ) then
        call rvar(ncid,istart,icount,ipar,echemsrc,'PAR_flux',.false.,sdim)
      end if
!RCOOH
      if ( ircooh/= 0 ) then
        call rvar(ncid,istart,icount,ircooh,echemsrc,'RCOOH_flux',.false.,sdim)
      end if
!SO2
      if ( iso2/= 0 ) then
        call rvar(ncid,istart,icount,iso2,echemsrc,'SO2_flux',.false.,sdim)
      end if
!TOL
      if ( itol/= 0 ) then
        call rvar(ncid,istart,icount,itol,echemsrc,'TOL_flux',.false.,sdim)
      end if
!XYL
      if ( ixyl/= 0 ) then
        call rvar(ncid,istart,icount,ixyl,echemsrc,'XYL_flux',.false.,sdim)
      end if

      ! Carbonaceous species / will be lumped to smoke tracers
      if ( ibchb/= 0 .and. ism1/=0) then
        call rvar(ncid,istart,icount,ibchb,echemsrc,'BC_flux',.false.,sdim)
      end if
      if ( iochb/= 0 .and. ism1/=0) then
        call rvar(ncid,istart,icount,iochb,echemsrc,'OC_flux',.false.,sdim)
      end if

      where (echemsrc(:,:,:) < d_zero ) echemsrc(:,:,:) = d_zero

      if ( do_parallel_netcdf_in ) then
        call closefile(ncid)
        deallocate (emtimeval)
        deallocate(rspace2)
      else
        if ( myid == iocpu ) then
          call closefile(ncid)
          deallocate (emtimeval)
          deallocate(rspace2)
          deallocate(rspace2_loc)
        else
          deallocate(rspace2_loc)
        end if
        call bcast(ifreq)
      end if
      if ( myid == italk ) then
        write(stdout,*) 'Emission read in success.'
      end if
    end subroutine read_bioburn_emission

    subroutine rvar(ncid,istart,icount,ind,echemsrc,cna,lh,sdim,cnb,cnc,cnd)
      implicit none

      integer(ik4) , intent(in) :: ncid,sdim
      integer(ik4) , dimension(4) , intent(in) :: istart , icount
      real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: echemsrc
      logical , intent(in) :: lh
      character(len=*) , intent(in) :: cna
      character(len=*) , intent(in) , optional :: cnb
      character(len=*) , intent(in) , optional :: cnc
      character(len=*) , intent(in) , optional :: cnd
      integer(ik4) :: ivarid
      integer(ik4) :: i , j , ind

      if ( do_parallel_netcdf_in ) then
        istatus = nf90_inq_varid(ncid, cna, ivarid)
        call check_ok(__FILE__,__LINE__, &
                      'Variable '//cna//' miss','CHEM_EMISS FILE')
        istatus = nf90_get_var(ncid,ivarid,rspace2, &
                istart(1:sdim),icount(1:sdim))
        call check_ok(__FILE__,__LINE__, &
                      'Variable '//cna//' read err','CHEM_EMISS FILE')
        if ( lh ) then  ! half of lumped Aromatics
          do i = ici1 , ici2
            do j = jci1 , jci2
              echemsrc(j,i,ind) = d_half*rspace2(j,i)
            end do
          end do
        else
          do i = ici1 , ici2
            do j = jci1 , jci2
              echemsrc(j,i,ind) = rspace2(j,i)
            end do
          end do
        end if
        if ( present(cnb) ) then
          istatus = nf90_inq_varid(ncid, cnb, ivarid)
          call check_ok(__FILE__,__LINE__, &
                        'Variable '//cnb//' miss','CHEM_EMISS FILE')
          istatus = nf90_get_var(ncid,ivarid,rspace2, &
            istart(1:sdim),icount(1:sdim))
          call check_ok(__FILE__,__LINE__, &
                        'Variable '//cnb//' read err','CHEM_EMISS FILE')
          do i = ici1 , ici2
            do j = jci1 , jci2
              echemsrc(j,i,ind) = rspace2(j,i) + echemsrc(j,i,ind)
            end do
          end do
        end if
        if ( present(cnc) ) then
          istatus = nf90_inq_varid(ncid, cnc, ivarid)
          call check_ok(__FILE__,__LINE__, &
                        'Variable '//cnc//' miss','CHEM_EMISS FILE')
          istatus = nf90_get_var(ncid,ivarid,rspace2, &
            istart(1:sdim),icount(1:sdim))
          call check_ok(__FILE__,__LINE__, &
                        'Variable '//cnc//' read err','CHEM_EMISS FILE')
          do i = ici1 , ici2
            do j = jci1 , jci2
              echemsrc(j,i,ind) = rspace2(j,i) + echemsrc(j,i,ind)
            end do
          end do
        end if
        if ( present(cnd) ) then
          istatus = nf90_inq_varid(ncid, cnd, ivarid)
          call check_ok(__FILE__,__LINE__, &
                      'Variable '//cnd//' miss','CHEM_EMISS FILE')
          istatus = nf90_get_var(ncid,ivarid,rspace2, &
            istart(1:sdim),icount(1:sdim))
          call check_ok(__FILE__,__LINE__, &
                        'Variable '//cnd//' read err','CHEM_EMISS FILE')
          do i = ici1 , ici2
            do j = jci1 , jci2
              echemsrc(j,i,ind) = rspace2(j,i) + echemsrc(j,i,ind)
            end do
          end do
        end if
      else
        if ( myid == iocpu ) then
          istatus = nf90_inq_varid(ncid, cna, ivarid)
          call check_ok(__FILE__,__LINE__, &
                        'Variable '//cna//' miss','CHEM_EMISS FILE')
          istatus = nf90_get_var(ncid,ivarid,rspace2, &
                  istart(1:sdim),icount(1:sdim))
          call check_ok(__FILE__,__LINE__, &
                        'Variable '//cna//' read err','CHEM_EMISS FILE')
          call grid_distribute(rspace2,rspace2_loc,jci1,jci2,ici1,ici2)
          if ( lh ) then  ! half of lumped Aromatics
            do i = ici1 , ici2
              do j = jci1 , jci2
                echemsrc(j,i,ind) = d_half*rspace2_loc(j,i)
              end do
            end do
          else
            do i = ici1 , ici2
              do j = jci1 , jci2
                echemsrc(j,i,ind) = rspace2_loc(j,i)
              end do
            end do
          end if
          if ( present(cnb) ) then
            istatus = nf90_inq_varid(ncid, cnb, ivarid)
            call check_ok(__FILE__,__LINE__, &
                          'Variable '//cnb//' miss','CHEM_EMISS FILE')
            istatus = nf90_get_var(ncid,ivarid,rspace2, &
              istart(1:sdim),icount(1:sdim))
            call check_ok(__FILE__,__LINE__, &
                          'Variable '//cnb//' read err','CHEM_EMISS FILE')
            call grid_distribute(rspace2,rspace2_loc,jci1,jci2,ici1,ici2)
            do i = ici1 , ici2
              do j = jci1 , jci2
                echemsrc(j,i,ind) = rspace2_loc(j,i) + echemsrc(j,i,ind)
              end do
            end do
          end if
          if ( present(cnc) ) then
            istatus = nf90_inq_varid(ncid, cnc, ivarid)
            call check_ok(__FILE__,__LINE__, &
                          'Variable '//cnc//' miss','CHEM_EMISS FILE')
            istatus = nf90_get_var(ncid,ivarid,rspace2, &
              istart(1:sdim),icount(1:sdim))
            call check_ok(__FILE__,__LINE__, &
                          'Variable '//cnc//' read err','CHEM_EMISS FILE')
            call grid_distribute(rspace2,rspace2_loc,jci1,jci2,ici1,ici2)
            do i = ici1 , ici2
              do j = jci1 , jci2
                echemsrc(j,i,ind) = rspace2_loc(j,i) + echemsrc(j,i,ind)
              end do
            end do
          end if
          if ( present(cnd) ) then
            istatus = nf90_inq_varid(ncid, cnd, ivarid)
            call check_ok(__FILE__,__LINE__, &
                         'Variable '//cnd//' miss','CHEM_EMISS FILE')
            istatus = nf90_get_var(ncid,ivarid,rspace2, &
              istart(1:sdim),icount(1:sdim))
            call check_ok(__FILE__,__LINE__, &
                          'Variable '//cnd//' read err','CHEM_EMISS FILE')
            call grid_distribute(rspace2,rspace2_loc,jci1,jci2,ici1,ici2)
            do i = ici1 , ici2
              do j = jci1 , jci2
                echemsrc(j,i,ind) = rspace2_loc(j,i) + echemsrc(j,i,ind)
              end do
            end do
          end if
        else
          call grid_distribute(rspace2,rspace2_loc,jci1,jci2,ici1,ici2)
          if ( lh ) then  ! half of lumped Aromatics
            do i = ici1 , ici2
              do j = jci1 , jci2
                echemsrc(j,i,ind) = d_half*rspace2_loc(j,i)
              end do
            end do
          else
            do i = ici1 , ici2
              do j = jci1 , jci2
                echemsrc(j,i,ind) = rspace2_loc(j,i)
              end do
            end do
          end if
          if ( present(cnb) ) then
            call grid_distribute(rspace2,rspace2_loc,jci1,jci2,ici1,ici2)
            do i = ici1 , ici2
              do j = jci1 , jci2
                echemsrc(j,i,ind) = rspace2_loc(j,i) + echemsrc(j,i,ind)
              end do
            end do
          end if
          if ( present(cnc) ) then
            call grid_distribute(rspace2,rspace2_loc,jci1,jci2,ici1,ici2)
            do i = ici1 , ici2
              do j = jci1 , jci2
                echemsrc(j,i,ind) = rspace2_loc(j,i) + echemsrc(j,i,ind)
              end do
            end do
          end if
          if ( present(cnd) ) then
            call grid_distribute(rspace2,rspace2_loc,jci1,jci2,ici1,ici2)
            do i = ici1 , ici2
              do j = jci1 , jci2
                echemsrc(j,i,ind) = rspace2_loc(j,i) + echemsrc(j,i,ind)
              end do
            end do
          end if
        end if
      end if
    end subroutine rvar

    integer(ik4) function chbc_search(idate)
      implicit none
      type(rcm_time_and_date) , intent(in) :: idate
      type(rcm_time_interval) :: tdif
      character(len=32) :: appdat1, appdat2
      if ( .not. do_parallel_netcdf_in .and. myid /= iocpu ) then
        chbc_search = 1
        return
      end if
      if (idate > chbc_idate(ibcnrec) .or. idate < chbc_idate(1)) then
        chbc_search = -1
      else
        tdif = idate - chbc_idate(1)
        ibcrec = (nint(tohours(tdif))/ibdyfrq)+1
        if ( ibcrec < 1 .or. ibcrec > ibcnrec ) then
          appdat1 = tochar(idate)
          write (stderr,*) 'Record is not found in CHBC file for ',appdat1
          appdat1 = tochar(chbc_idate(1))
          appdat2 = tochar(chbc_idate(ibcnrec))
          write (stderr,*) 'Range is : ', appdat1, '-', appdat2
          call fatal(__FILE__,__LINE__, &
                        'CHBC READ')
        end if
        chbc_search = ibcrec
      end if
    end function chbc_search

    subroutine open_chbc(idate)
      implicit none
      type(rcm_time_and_date) , intent(in) :: idate
      character(len=11) :: ctime
      integer(ik4) :: ibcid , idimid , itvar , i , chkdiff
      real(rkx) , dimension(:) , allocatable :: icbc_nctime
      character(len=64) :: icbc_timeunits , icbc_timecal

      if ( .not. do_parallel_netcdf_in .and. myid /= iocpu ) then
        allocate(rspace3_loc(jce1:jce2,ice1:ice2,kz))
        return
      end if

      call close_chbc
      write (ctime, '(a)') tochar10(idate)
      if ( igaschem == 1 ) then
        icbcname = trim(dirglob)//pthsep//trim(domname)// &
                   '_CHBC.'//trim(ctime)//'.nc'
        call openfile_withname(icbcname,ichin)
        call check_dims(ichin)
        ibcid = ichin
      end if
      if ( iaerosol == 1 ) then
        icbcname = trim(dirglob)//pthsep//trim(domname)// &
                   '_AEBC.'//trim(ctime)//'.nc'
        call openfile_withname(icbcname,iaein)
        call check_dims(iaein)
        ibcid = iaein
      end if
      if ( ioxclim == 1 ) then
        icbcname = trim(dirglob)//pthsep//trim(domname)// &
                   '_OXBC.'//trim(ctime)//'.nc'
        call openfile_withname(icbcname,ioxin)
        call check_dims(iaein)
        ibcid = ioxin
      end if
      ibcrec = 1
      ibcnrec = 0
      istatus = nf90_inq_dimid(ibcid, 'time', idimid)
      call check_ok(__FILE__,__LINE__, &
                    'Dimension time miss', 'ICBC FILE')
      istatus = nf90_inquire_dimension(ibcid, idimid, len=ibcnrec)
      call check_ok(__FILE__,__LINE__, &
                    'Dimension time read error', 'ICBC FILE')
      if ( ibcnrec < 1 ) then
        write (stderr,*) 'Time var in ICBC has zero dim.'
        call fatal(__FILE__,__LINE__, &
                      'ICBC READ')
      end if
      istatus = nf90_inq_varid(ibcid, 'time', itvar)
      call check_ok(__FILE__,__LINE__, &
                    'variable time miss', 'ICBC FILE')
      istatus = nf90_get_att(ibcid, itvar, 'units', icbc_timeunits)
      call check_ok(__FILE__,__LINE__, &
                    'variable time units miss','ICBC FILE')
      istatus = nf90_get_att(ibcid, itvar, 'calendar', icbc_timecal)
      call check_ok(__FILE__,__LINE__, &
                    'variable time calendar miss','ICBC FILE')
      allocate(icbc_nctime(ibcnrec), stat=istatus)
      if ( istatus /= 0 ) then
        write(stderr,*) 'Memory allocation error in ICBC for time real values'
        call fatal(__FILE__,__LINE__, &
                      'ICBC READ')
      end if
      allocate(chbc_idate(ibcnrec), stat=istatus)
      if ( istatus /= 0 ) then
        write(stderr,*) 'Memory allocation error in ICBC for time array'
        call fatal(__FILE__,__LINE__, &
                      'ICBC READ')
      end if
      istatus = nf90_get_var(ibcid, itvar, icbc_nctime)
      call check_ok(__FILE__,__LINE__, &
                    'variable time read error', 'ICBC FILE')
      do i = 1 , ibcnrec
        chbc_idate(i) = timeval2date(icbc_nctime(i), &
                                     icbc_timeunits,icbc_timecal)
      end do
      if ( ibcnrec > 1 ) then
        chkdiff = nint(icbc_nctime(2) - icbc_nctime(1))
        if (chkdiff /= ibdyfrq) then
          write (stderr,*) 'Time var in ICBC inconsistency.'
          write (stderr,*) 'Expecting ibdyfrq = ', ibdyfrq
          write (stderr,*) 'Found     ibdyfrq = ', chkdiff
          call fatal(__FILE__,__LINE__, &
                        'ICBC READ')
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
      if ( do_parallel_netcdf_in ) then
        allocate(rspace3(jde1:jde2,ide1:ide2,kz))
      else
        allocate(rspace3(jx,iy,kz))
        allocate(rspace3_loc(jce1:jce2,ice1:ice2,kz))
      end if
    end subroutine open_chbc

    subroutine read_chbc(chebdio)
      implicit none
      real(rkx) , dimension (:,:,:,:) , pointer , intent(inout) :: chebdio
      integer(ik4) , dimension(4) :: istart , icount
      integer(ik4) :: i , j , k, n , iafter

      if ( do_parallel_netcdf_in ) then
        istart(1) = jde1
        istart(2) = ide1
        istart(3) = 1
        istart(4) = ibcrec
        icount(1) = jde2-jde1+1
        icount(2) = ide2-ide1+1
        icount(3) = kz
        icount(4) = 1
        iafter = 0
        if ( igaschem == 1 ) then
          do n = 1 , n_chbcvar
            istatus = nf90_get_var(ichin, chbc_ivar(n), rspace3, istart, icount)
            call check_ok(__FILE__,__LINE__, &
                'variable '//trim(chbcname(n))//' read error','CHBC FILE ERROR')
            do k = 1 , kz
              do i = ice1 , ice2
                do j = jce1 , jce2
                  chebdio(j,i,k,n) = rspace3(j,i,k)
                end do
              end do
            end do
            iafter = iafter + 1
          end do
        end if
        if ( iaerosol == 1 ) then
          do n = 1 , n_aebcvar
            istatus = nf90_get_var(iaein, aebc_ivar(n), rspace3, istart, icount)
            call check_ok(__FILE__,__LINE__, &
                'variable '//trim(aebcname(n))//' read error','AEBC FILE ERROR')
            do k = 1 , kz
              do i = ice1 , ice2
                do j = jce1 , jce2
                  chebdio(j,i,k,iafter+1) = rspace3(j,i,k)
                end do
              end do
            end do
            iafter = iafter + 1
          end do
        end if
        if ( ioxclim == 1 ) then
          do n = 1 , n_oxbcvar
            istatus = nf90_get_var(ioxin, oxbc_ivar(n), rspace3, istart, icount)
            call check_ok(__FILE__,__LINE__, &
                'variable '//trim(oxbcname(n))//' read error','OXBC FILE ERROR')
            do k = 1 , kz
              do i = ice1 , ice2
                do j = jce1 , jce2
                  chebdio(j,i,k,iafter+n) = rspace3(j,i,k)
                end do
              end do
            end do
          end do
        end if
      else
        if ( myid == iocpu ) then
          istart(1) = 1
          istart(2) = 1
          istart(3) = 1
          istart(4) = ibcrec
          icount(1) = jx
          icount(2) = iy
          icount(3) = kz
          icount(4) = 1
          iafter = 0
          if ( igaschem == 1 ) then
            do n = 1 , n_chbcvar
              istatus = nf90_get_var(ichin,chbc_ivar(n),rspace3,istart,icount)
              call check_ok(__FILE__,__LINE__, &
                'variable '//trim(chbcname(n))//' read error','CHBC FILE ERROR')
              call grid_distribute(rspace3,rspace3_loc,jce1,jce2,ice1,ice2,1,kz)
              do k = 1 , kz
                do i = ice1 , ice2
                  do j = jce1 , jce2
                    chebdio(j,i,k,n) = rspace3_loc(j,i,k)
                  end do
                end do
              end do
              iafter = iafter + 1
            end do
          end if
          if ( iaerosol == 1 ) then
            do n = 1 , n_aebcvar
              istatus = nf90_get_var(iaein,aebc_ivar(n),rspace3,istart,icount)
              call check_ok(__FILE__,__LINE__, &
                'variable '//trim(aebcname(n))//' read error','AEBC FILE ERROR')
              call grid_distribute(rspace3,rspace3_loc,jce1,jce2,ice1,ice2,1,kz)
              do k = 1 , kz
                do i = ice1 , ice2
                  do j = jce1 , jce2
                    chebdio(j,i,k,iafter+1) = rspace3_loc(j,i,k)
                  end do
                end do
              end do
              iafter = iafter + 1
            end do
          end if
          if ( ioxclim == 1 ) then
            do n = 1 , n_oxbcvar
              istatus = nf90_get_var(ioxin,oxbc_ivar(n),rspace3,istart,icount)
              call check_ok(__FILE__,__LINE__, &
                'variable '//trim(oxbcname(n))//' read error','OXBC FILE ERROR')
              call grid_distribute(rspace3,rspace3_loc,jce1,jce2,ice1,ice2,1,kz)
              do k = 1 , kz
                do i = ice1 , ice2
                  do j = jce1 , jce2
                    chebdio(j,i,k,iafter+n) = rspace3_loc(j,i,k)
                  end do
                end do
              end do
            end do
          end if
        else
          iafter = 0
          if ( igaschem == 1 ) then
            do n = 1 , n_chbcvar
              call grid_distribute(rspace3,rspace3_loc,jce1,jce2,ice1,ice2,1,kz)
              do k = 1 , kz
                do i = ice1 , ice2
                  do j = jce1 , jce2
                    chebdio(j,i,k,n) = rspace3_loc(j,i,k)
                  end do
                end do
              end do
              iafter = iafter + 1
            end do
          end if
          if ( iaerosol == 1 ) then
            do n = 1 , n_aebcvar
              call grid_distribute(rspace3,rspace3_loc,jce1,jce2,ice1,ice2,1,kz)
              do k = 1 , kz
                do i = ice1 , ice2
                  do j = jce1 , jce2
                    chebdio(j,i,k,iafter+1) = rspace3_loc(j,i,k)
                  end do
                end do
              end do
              iafter = iafter + 1
            end do
          end if
          if ( ioxclim == 1 ) then
            do n = 1 , n_oxbcvar
              call grid_distribute(rspace3,rspace3_loc,jce1,jce2,ice1,ice2,1,kz)
              do k = 1 , kz
                do i = ice1 , ice2
                  do j = jce1 , jce2
                    chebdio(j,i,k,iafter+n) = rspace3_loc(j,i,k)
                  end do
                end do
              end do
            end do
          end if
        end if
      end if
      where (chebdio < d_zero) chebdio = d_zero
    end subroutine read_chbc

    subroutine close_chbc
      implicit none
      if ( ichin >= 0 ) then
        call closefile(ichin)
        ichin = -1
      end if
      if ( iaein >= 0 ) then
        call closefile(iaein)
        iaein = -1
      end if
      if ( ioxin >= 0 ) then
        call closefile(ioxin)
        ioxin = -1
      end if
      if ( allocated(chbc_idate) )   deallocate(chbc_idate)
      if ( associated(rspace3) )     deallocate(rspace3)
      if ( associated(rspace3_loc) ) deallocate(rspace3_loc)
    end subroutine close_chbc

    subroutine check_ok(f,l,m1,mf)
      implicit none
      character(*) , intent(in) :: f, m1 , mf
      integer(ik4) , intent(in) :: l
      if (istatus /= nf90_noerr) then
        write (stderr,*) trim(m1)
        write (stderr,*) nf90_strerror(istatus)
        call fatal(f,l,trim(mf))
      end if
    end subroutine check_ok

end module mod_che_ncio
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
