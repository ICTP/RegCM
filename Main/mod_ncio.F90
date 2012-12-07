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
module mod_ncio
!
  use netcdf
  use mod_runparams
  use mod_dynparam
  use mod_ensemble
  use mod_mpmessage
  use mod_mppparam
  use mod_memutil
  use mod_nchelper
  use mod_domain
!
  private

  public :: ivarname_lookup
  public :: read_domain_info , read_subdomain_info
  public :: open_icbc , icbc_search , read_icbc , close_icbc
!
  integer(ik4) :: ibcin
  integer(ik4) :: istatus
  integer(ik4) :: ibcrec , ibcnrec
  type(rcm_time_and_date) , dimension(:) , allocatable :: icbc_idate
  integer(ik4) , dimension(7) :: icbc_ivar

  data ibcin   /-1/
  data ibcrec  / 1/
  data ibcnrec / 0/

  real(rk8) , dimension(:,:) , allocatable :: rspace2
  real(rk8) , dimension(:,:,:) , allocatable :: rspace3

  contains

  subroutine read_domain_info(ht,lnd,mask,xlat,xlon,dlat,dlon, &
                              msfx,msfd,coriol,snowam,hlake)
    implicit none
    real(rk8) , pointer , dimension(:,:) , intent(out) :: ht
    real(rk8) , pointer , dimension(:,:) , intent(out) :: lnd
    real(rk8) , pointer , dimension(:,:) , intent(out) :: mask
    real(rk8) , pointer , dimension(:,:) , intent(out) :: xlat
    real(rk8) , pointer , dimension(:,:) , intent(out) :: xlon
    real(rk8) , pointer , dimension(:,:) , intent(out) :: dlat
    real(rk8) , pointer , dimension(:,:) , intent(out) :: dlon
    real(rk8) , pointer , dimension(:,:) , intent(out) :: msfx
    real(rk8) , pointer , dimension(:,:) , intent(out) :: msfd
    real(rk8) , pointer , dimension(:,:) , intent(out) :: coriol
    real(rk8) , pointer , dimension(:,:) , intent(out) :: snowam
    real(rk8) , pointer , dimension(:,:) , intent(out) :: hlake
    character(len=256) :: dname
    integer(ik4) :: idmin
    integer(ik4) , dimension(2) :: istart , icount
    real(rk8) , dimension(:,:) , allocatable :: rspace

    dname = trim(dirter)//pthsep//trim(domname)//'_DOMAIN000.nc'
    if ( myid == italk ) then
      write(stdout,*) 'Reading Domain file : ', trim(dname)
    end if
    istart(1) = global_dot_jstart
    istart(2) = global_dot_istart
    icount(1) = global_dot_jend-global_dot_jstart+1
    icount(2) = global_dot_iend-global_dot_istart+1
    allocate(rspace(icount(1),icount(2)))
    call openfile_withname(dname,idmin)
    call check_domain(idmin)
    call read_var1d_static(idmin,'sigma',sigma)
    call read_var2d_static(idmin,'xlat',rspace,istart=istart,icount=icount)
    xlat(jde1:jde2,ide1:ide2) = rspace
    call read_var2d_static(idmin,'xlon',rspace,istart=istart,icount=icount)
    xlon(jde1:jde2,ide1:ide2) = rspace
    call read_var2d_static(idmin,'dlat',rspace,istart=istart,icount=icount)
    dlat(jde1:jde2,ide1:ide2) = rspace
    call read_var2d_static(idmin,'dlon',rspace,istart=istart,icount=icount)
    dlon(jde1:jde2,ide1:ide2) = rspace
    call read_var2d_static(idmin,'topo',rspace,istart=istart,icount=icount)
    if ( ensemble_run ) then
      write(stdout,*) 'Appling perturbation to input dataset:'
      if ( lperturb_topo ) then
        write(stdout,'(a,f7.2,a)') 'Topo with value ', &
          perturb_frac_topo*d_100,'%'
        call randify(rspace,perturb_frac_topo,icount(1),icount(2))
      end if
    end if
    ht(jde1:jde2,ide1:ide2) = rspace
    call read_var2d_static(idmin,'mask',rspace,istart=istart,icount=icount)
    mask(jde1:jde2,ide1:ide2) = rspace
    call read_var2d_static(idmin,'landuse',rspace,istart=istart,icount=icount)
    lnd(jde1:jde2,ide1:ide2) = rspace
    call read_var2d_static(idmin,'xmap',rspace,istart=istart,icount=icount)
    msfx(jde1:jde2,ide1:ide2) = rspace
    call read_var2d_static(idmin,'dmap',rspace,istart=istart,icount=icount)
    msfd(jde1:jde2,ide1:ide2) = rspace
    call read_var2d_static(idmin,'coriol',rspace,istart=istart,icount=icount)
    coriol(jde1:jde2,ide1:ide2) = rspace
    rspace = d_zero
    call read_var2d_static(idmin,'snowam',rspace,.true.,istart,icount)
    snowam(jde1:jde2,ide1:ide2) = rspace
    if ( lakemod == 1 ) then
      call read_var2d_static(idmin,'dhlake',rspace,.true.,istart,icount)
      hlake(jde1:jde2,ide1:ide2) = rspace
    end if
    call closefile(idmin)
    deallocate(rspace)
  end subroutine read_domain_info

  subroutine read_subdomain_info(ht1,lnd1,mask1,xlat1,xlon1,hlake1)
    implicit none
    real(rk8) , pointer , dimension(:,:,:) , intent(out) :: ht1
    real(rk8) , pointer , dimension(:,:,:) , intent(out) :: lnd1
    real(rk8) , pointer , dimension(:,:,:) , intent(out) :: mask1
    real(rk8) , pointer , dimension(:,:,:) , intent(out) :: xlat1
    real(rk8) , pointer , dimension(:,:,:) , intent(out) :: xlon1
    real(rk8) , pointer , dimension(:,:,:) , intent(out) :: hlake1
    character(len=256) :: dname
    integer(ik4) :: idmin
    integer(ik4) , dimension(2) :: istart , icount
    real(rk8) , dimension(:,:) , allocatable :: rspace
    character(len=3) :: sbstring

    write (sbstring,'(i0.3)') nsg
    dname = trim(dirter)//pthsep//trim(domname)//'_DOMAIN'//sbstring//'.nc'
    if ( myid == italk ) then
      write(stdout,*) 'Reading Sub-Domain file : ', trim(dname)
    end if
    istart(1) = (global_dot_jstart-1)*nsg+1
    istart(2) = (global_dot_istart-1)*nsg+1
    icount(1) = (global_dot_jend-global_dot_jstart+1)*nsg
    icount(2) = (global_dot_iend-global_dot_istart+1)*nsg
    allocate(rspace(icount(1),icount(2)))
    call openfile_withname(dname,idmin)
    call read_var2d_static(idmin,'xlat',rspace,istart=istart,icount=icount)
    call input_reorder(rspace,xlat1,jde1,jde2,ide1,ide2)
    call read_var2d_static(idmin,'xlon',rspace,istart=istart,icount=icount)
    call input_reorder(rspace,xlon1,jde1,jde2,ide1,ide2)
    call read_var2d_static(idmin,'topo',rspace,istart=istart,icount=icount)
    if ( ensemble_run ) then
      write(stdout,*) 'Appling perturbation to input dataset:'
      if ( lperturb_topo ) then
        write(stdout,'(a,f7.2,a)') 'Topo with value ', &
          perturb_frac_topo*d_100,'%'
        call randify(rspace,perturb_frac_topo,icount(1),icount(2))
      end if
    end if
    call input_reorder(rspace,ht1,jde1,jde2,ide1,ide2)
    call read_var2d_static(idmin,'mask',rspace,istart=istart,icount=icount)
    call input_reorder(rspace,mask1,jde1,jde2,ide1,ide2)
    call read_var2d_static(idmin,'landuse',rspace,istart=istart,icount=icount)
    call input_reorder(rspace,lnd1,jde1,jde2,ide1,ide2)
    if ( lakemod == 1 ) then
      call read_var2d_static(idmin,'dhlake',rspace,.true.,istart,icount)
      call input_reorder(rspace,hlake1,jde1,jde2,ide1,ide2)
    end if
    call closefile(idmin)
    deallocate(rspace)
  end subroutine read_subdomain_info

  integer function icbc_search(idate)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate
    type(rcm_time_interval) :: tdif
    character(len=32) :: appdat1, appdat2
    if (idate > icbc_idate(ibcnrec) .or. idate < icbc_idate(1)) then
      icbc_search = -1
    else
      tdif = idate-icbc_idate(1)
      ibcrec = (idnint(tohours(tdif))/ibdyfrq)+1
      if ( ibcrec < 1 .or. ibcrec > ibcnrec ) then
        appdat1 = tochar(idate)
        write (stderr,*) 'Record is not found in ICBC file for ',appdat1
        appdat1 = tochar(icbc_idate(1))
        appdat2 = tochar(icbc_idate(ibcnrec))
        write (stderr,*) 'Range is : ', appdat1, '-', appdat2
        call fatal(__FILE__,__LINE__,'ICBC READ')
      end if
      icbc_search = ibcrec
    end if 
  end function icbc_search

  subroutine open_icbc(idate)
    type(rcm_time_and_date) , intent(in) :: idate
    character(10) :: ctime
    integer(ik4) :: idimid , itvar , i , chkdiff , nnj , nni
    real(rk8) , dimension(:) , allocatable :: icbc_nctime
    character(64) :: icbc_timeunits , icbc_timecal
    character(len=256) :: icbcname

    call close_icbc
    write (ctime, '(i10)') toint10(idate)
    icbcname = trim(dirglob)//pthsep//trim(domname)//'_ICBC.'//ctime//'.nc'
    call openfile_withname(icbcname,ibcin)
    ibcrec = 1
    ibcnrec = 0
    call check_domain(ibcin,.true.)
    istatus = nf90_inq_dimid(ibcin, 'time', idimid)
    call check_ok(__FILE__,__LINE__,'Dimension time miss', 'ICBC FILE')
    istatus = nf90_inquire_dimension(ibcin, idimid, len=ibcnrec)
    call check_ok(__FILE__,__LINE__,'Dimension time read error', 'ICBC FILE')
    if ( ibcnrec < 1 ) then
      write (stderr,*) 'Time var in ICBC has zero dim.'
      call fatal(__FILE__,__LINE__,'ICBC READ')
    end if
    istatus = nf90_inq_varid(ibcin, 'time', itvar)
    call check_ok(__FILE__,__LINE__,'variable time miss', 'ICBC FILE')
    istatus = nf90_get_att(ibcin, itvar, 'units', icbc_timeunits)
    call check_ok(__FILE__,__LINE__,'variable time units miss','ICBC FILE')
    istatus = nf90_get_att(ibcin, itvar, 'calendar', icbc_timecal)
    call check_ok(__FILE__,__LINE__,'variable time calendar miss','ICBC FILE')
    allocate(icbc_nctime(ibcnrec), stat=istatus)
    if ( istatus /= 0 ) then
      write(stderr,*) 'Memory allocation error in ICBC for time real values'
      call fatal(__FILE__,__LINE__,'ICBC READ')
    end if
    allocate(icbc_idate(ibcnrec), stat=istatus)
    if ( istatus /= 0 ) then
      write(stderr,*) 'Memory allocation error in ICBC for time array'
      call fatal(__FILE__,__LINE__,'ICBC READ')
    end if
    istatus = nf90_get_var(ibcin, itvar, icbc_nctime)
    call check_ok(__FILE__,__LINE__,'variable time read error', 'ICBC FILE')
    do i = 1 , ibcnrec
      icbc_idate(i) = timeval2date(icbc_nctime(i), icbc_timeunits, icbc_timecal)
    end do
    if ( ibcnrec > 1 ) then
      chkdiff = idnint(icbc_nctime(2) - icbc_nctime(1))
      if (chkdiff /= ibdyfrq) then
        write (stderr,*) 'Time var in ICBC inconsistency.'
        write (stderr,*) 'Expecting ibdyfrq = ', ibdyfrq
        write (stderr,*) 'Found     ibdyfrq = ', chkdiff
        call fatal(__FILE__,__LINE__,'ICBC READ')
      end if
    end if
    deallocate(icbc_nctime)
    istatus = nf90_inq_varid(ibcin, 'ps', icbc_ivar(1))
    call check_ok(__FILE__,__LINE__,'variable ps miss', 'ICBC FILE')
    istatus = nf90_inq_varid(ibcin, 'ts', icbc_ivar(2))
    call check_ok(__FILE__,__LINE__,'variable ts miss', 'ICBC FILE')
    istatus = nf90_inq_varid(ibcin, 'u', icbc_ivar(3))
    call check_ok(__FILE__,__LINE__,'variable u miss', 'ICBC FILE')
    istatus = nf90_inq_varid(ibcin, 'v', icbc_ivar(4))
    call check_ok(__FILE__,__LINE__,'variable v miss', 'ICBC FILE')
    istatus = nf90_inq_varid(ibcin, 't', icbc_ivar(5))
    call check_ok(__FILE__,__LINE__,'variable t miss', 'ICBC FILE')
    istatus = nf90_inq_varid(ibcin, 'qv', icbc_ivar(6))
    call check_ok(__FILE__,__LINE__,'variable qv miss', 'ICBC FILE')

    nnj = global_dot_jend-global_dot_jstart+1
    nni = global_dot_iend-global_dot_istart+1
    allocate(rspace2(nnj,nni))
    allocate(rspace3(nnj,nni,kz))
  end subroutine open_icbc

  subroutine read_icbc(ps,ts,u,v,t,qv)
    implicit none
    real(rk8) , pointer , dimension(:,:) , intent(out) :: ps
    real(rk8) , pointer , dimension(:,:) , intent(out) :: ts
    real(rk8) , pointer , dimension(:,:,:) , intent(out) :: u
    real(rk8) , pointer , dimension(:,:,:) , intent(out) :: v
    real(rk8) , pointer , dimension(:,:,:) , intent(out) :: t
    real(rk8) , pointer , dimension(:,:,:) , intent(out) :: qv
    integer(ik4) , dimension(4) :: istart , icount

    istart(1) = global_dot_jstart
    istart(2) = global_dot_istart
    istart(3) = ibcrec
    icount(1) = global_dot_jend-global_dot_jstart+1
    icount(2) = global_dot_iend-global_dot_istart+1
    icount(3) = 1
    istatus = nf90_get_var(ibcin,icbc_ivar(1),rspace2,istart(1:3),icount(1:3))
    call check_ok(__FILE__,__LINE__,'variable ps read error', 'ICBC FILE')
    if ( ensemble_run ) then
      if ( lperturb_ps ) then
        write(stdout,'(a,f7.2,a)') 'PS with value ',perturb_frac_ps*d_100,'%'
        call randify(rspace2,perturb_frac_ps,icount(1),icount(2))
      end if
    end if
    ps(jce1:jce2,ice1:ice2) = rspace2(jce1:jce2,ice1:ice2)
    istatus = nf90_get_var(ibcin,icbc_ivar(2),rspace2,istart(1:3),icount(1:3))
    call check_ok(__FILE__,__LINE__,'variable ts read error', 'ICBC FILE')
    if ( ensemble_run ) then
      write(stdout,*) 'Appling perturbation to input dataset:'
      if ( lperturb_ts ) then
        write(stdout,'(a,f7.2,a)') 'TS with value ',perturb_frac_ts*d_100,'%'
        call randify(rspace2,perturb_frac_ts,icount(1),icount(2))
      end if
    end if
    ts(jce1:jce2,ice1:ice2) = rspace2(jce1:jce2,ice1:ice2)
    istart(1) = global_dot_jstart
    istart(2) = global_dot_istart
    istart(3) = 1
    istart(4) = ibcrec
    icount(1) = global_dot_jend-global_dot_jstart+1
    icount(2) = global_dot_iend-global_dot_istart+1
    icount(3) = kz
    icount(4) = 1
    istatus = nf90_get_var(ibcin,icbc_ivar(3),rspace3,istart,icount)
    call check_ok(__FILE__,__LINE__,'variable u read error', 'ICBC FILE')
    if ( ensemble_run ) then
      if ( lperturb_u ) then
        write(stdout,'(a,f7.2,a)') 'U  with value ',perturb_frac_u*d_100,'%'
        call randify(rspace3,perturb_frac_u,icount(1),icount(2),icount(3))
      end if
    end if
    u(jde1:jde2,ide1:ide2,1:kz) = rspace3
    istatus = nf90_get_var(ibcin,icbc_ivar(4),rspace3,istart,icount)
    call check_ok(__FILE__,__LINE__,'variable v read error', 'ICBC FILE')
    if ( ensemble_run ) then
      if ( lperturb_v ) then
        write(stdout,'(a,f7.2,a)') 'V  with value ',perturb_frac_v*d_100,'%'
        call randify(rspace3,perturb_frac_v,icount(1),icount(2),icount(3))
      end if
    end if
    v(jde1:jde2,ide1:ide2,1:kz) = rspace3
    istatus = nf90_get_var(ibcin,icbc_ivar(5),rspace3,istart,icount)
    call check_ok(__FILE__,__LINE__,'variable t read error', 'ICBC FILE')
    if ( ensemble_run ) then
      if ( lperturb_t ) then
        write(stdout,'(a,f7.2,a)') 'T  with value ',perturb_frac_t*d_100,'%'
        call randify(rspace3,perturb_frac_t,icount(1),icount(2),icount(3))
      end if
    end if
    t(jce1:jce2,ice1:ice2,1:kz) = rspace3(jce1:jce2,ice1:ice2,1:kz)
    istatus = nf90_get_var(ibcin,icbc_ivar(6),rspace3,istart,icount)
    call check_ok(__FILE__,__LINE__,'variable qx read error', 'ICBC FILE')
    if ( ensemble_run ) then
      if ( lperturb_q ) then
        write(stdout,'(a,f7.2,a)') 'Q  with value ',perturb_frac_q*d_100,'%'
        call randify(rspace3,perturb_frac_q,icount(1),icount(2),icount(3))
      end if
    end if
    qv(jce1:jce2,ice1:ice2,1:kz) = rspace3(jce1:jce2,ice1:ice2,1:kz)
  end subroutine read_icbc

  subroutine close_icbc
    implicit none
    if (ibcin >= 0) then
      istatus = nf90_close(ibcin)
      call check_ok(__FILE__,__LINE__,'Error Close ICBC file','ICBC FILE')
      if ( allocated(icbc_idate) ) deallocate(icbc_idate)
      ibcin = -1
    end if
    if ( allocated(rspace2) ) deallocate(rspace2)
    if ( allocated(rspace3) ) deallocate(rspace3)
  end subroutine close_icbc

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

end module mod_ncio
