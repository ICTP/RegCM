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

  contains

  subroutine read_domain_info_par
    implicit none
    character(len=256) :: dname
    integer(ik4) :: idmin

    if ( myid == iocpu ) then
      dname = trim(dirter)//pthsep//trim(domname)//'_DOMAIN000.nc'
      write(stdout,*) 'Reading Domain file : ', trim(dname)
      call openfile_withname(dname,idmin)
      call read_domain(idmin)
      call closefile(idmin)
      if ( ensemble_run ) then
        write(stdout,*) 'Appling perturbation to input dataset:'
        if ( lperturb_topo ) then
          write(stdout,'(a,f7.2,a)') 'Topo with value ', &
            perturb_frac_topo*d_100,'%'
          call randify(mddom_io%ht,perturb_frac_topo,jx,iy)
        end if
      end if
    end if
  end subroutine read_domain_info_par

  subroutine read_domain_info
    implicit none
    character(len=256) :: dname
    integer(ik4) :: idmin

    if ( myid == iocpu ) then
      dname = trim(dirter)//pthsep//trim(domname)//'_DOMAIN000.nc'
      write(stdout,*) 'Reading Domain file : ', trim(dname)
      call openfile_withname(dname,idmin)
      call read_domain(idmin)
      call closefile(idmin)
      if ( ensemble_run ) then
        write(stdout,*) 'Appling perturbation to input dataset:'
        if ( lperturb_topo ) then
          write(stdout,'(a,f7.2,a)') 'Topo with value ', &
            perturb_frac_topo*d_100,'%'
          call randify(mddom_io%ht,perturb_frac_topo,jx,iy)
        end if
      end if
    end if
  end subroutine read_domain_info

  subroutine read_subdomain_info(ht1,lnd1,mask1,xlat1,xlon1,hlake1)
    implicit none
    real(rk8) , pointer , dimension(:,:,:) , intent(out) :: ht1
    real(rk8) , pointer , dimension(:,:,:) , intent(out) :: lnd1
    real(rk8) , pointer , dimension(:,:,:) , intent(out) :: mask1
    real(rk8) , pointer , dimension(:,:,:) , intent(out) :: xlat1
    real(rk8) , pointer , dimension(:,:,:) , intent(out) :: xlon1
    real(rk8) , pointer , dimension(:,:,:) , intent(out) :: hlake1
    integer(ik4) :: isdmin , ivarid
    character(len=3) :: sbstring
    character(len=256) :: sdname
    real(rk4) , dimension(:,:) , allocatable :: sp2d1

    if ( nsg < 2 ) return

    write (sbstring,'(i0.3)') nsg
    sdname = trim(dirter)//pthsep//trim(domname)//'_DOMAIN'//sbstring//'.nc'
    write(stdout,*) 'Reading Sub-Domain file : ', trim(sdname)

    allocate(sp2d1(jxsg,iysg))
    call openfile_withname(sdname,isdmin)
    istatus = nf90_inq_varid(isdmin, 'topo', ivarid)
    call check_ok(__FILE__,__LINE__,'Variable topo miss', 'SUBDOMAIN FILE')
    istatus = nf90_get_var(isdmin, ivarid, sp2d1)
    call check_ok(__FILE__,__LINE__,'Variable topo read error','SUBDOMAIN FILE')
    call input_reorder(sp2d1,ht1)
    istatus = nf90_inq_varid(isdmin, 'landuse', ivarid)
    call check_ok(__FILE__,__LINE__,'Variable landuse miss','SUBDOMAIN FILE')
    istatus = nf90_get_var(isdmin, ivarid, sp2d1)
    call check_ok(__FILE__,__LINE__,'Variable landuse read error', &
                  'SUBDOMAIN FILE')
    call input_reorder(sp2d1,lnd1)
    istatus = nf90_inq_varid(isdmin, 'mask', ivarid)
    call check_ok(__FILE__,__LINE__,'Variable mask miss','SUBDOMAIN FILE')
    istatus = nf90_get_var(isdmin, ivarid, sp2d1)
    call check_ok(__FILE__,__LINE__,'Variable mask read error', &
                  'SUBDOMAIN FILE')
    call input_reorder(sp2d1,mask1)
    istatus = nf90_inq_varid(isdmin, 'xlat', ivarid)
    call check_ok(__FILE__,__LINE__,'Variable xlat miss','SUBDOMAIN FILE')
    istatus = nf90_get_var(isdmin, ivarid, sp2d1)
    call check_ok(__FILE__,__LINE__,'Variable xlat read error','SUBDOMAIN FILE')
    call input_reorder(sp2d1,xlat1)
    istatus = nf90_inq_varid(isdmin, 'xlon', ivarid)
    call check_ok(__FILE__,__LINE__,'Variable xlon miss','SUBDOMAIN FILE')
    istatus = nf90_get_var(isdmin, ivarid, sp2d1)
    call check_ok(__FILE__,__LINE__,'Variable xlon read error','SUBDOMAIN FILE')
    call input_reorder(sp2d1,xlon1)
    istatus = nf90_inq_varid(isdmin, 'mask', ivarid)
    call check_ok(__FILE__,__LINE__,'Variable mask miss','SUBDOMAIN FILE')
    istatus = nf90_get_var(isdmin, ivarid, sp2d1)
    call input_reorder(sp2d1,mask1)
    call check_ok(__FILE__,__LINE__,'Variable mask read error','SUBDOMAIN FILE')
    if ( lakemod == 1 ) then
      istatus = nf90_inq_varid(isdmin, 'dhlake', ivarid)
      call check_ok(__FILE__,__LINE__,'Variable dhlake miss','SUBDOMAIN FILE')
      istatus = nf90_get_var(isdmin, ivarid, sp2d1)
      call check_ok(__FILE__,__LINE__,'Variable dhlake read error', &
                    'SUBDOMAIN FILE')
      call input_reorder(sp2d1,hlake1)
    end if
    istatus = nf90_close(isdmin)
    call check_ok(__FILE__,__LINE__,'SubDomain file close error', &
                 'SUBDOMAIN FILE')
    deallocate(sp2d1)
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
        write (6,*) 'Record is not found in ICBC file for ',appdat1
        appdat1 = tochar(icbc_idate(1))
        appdat2 = tochar(icbc_idate(ibcnrec))
        write (6,*) 'Range is : ', appdat1, '-', appdat2
        call fatal(__FILE__,__LINE__,'ICBC READ')
      end if
      icbc_search = ibcrec
    end if 
  end function icbc_search

  subroutine open_icbc(idate)
    type(rcm_time_and_date) , intent(in) :: idate
    character(10) :: ctime
    integer(ik4) :: idimid , itvar , i , chkdiff
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
      write (6,*) 'Time var in ICBC has zero dim.'
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
      write(6,*) 'Memory allocation error in ICBC for time real values'
      call fatal(__FILE__,__LINE__,'ICBC READ')
    end if
    allocate(icbc_idate(ibcnrec), stat=istatus)
    if ( istatus /= 0 ) then
      write(6,*) 'Memory allocation error in ICBC for time array'
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
        write (6,*) 'Time var in ICBC inconsistency.'
        write (6,*) 'Expecting ibdyfrq = ', ibdyfrq
        write (6,*) 'Found     ibdyfrq = ', chkdiff
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
  end subroutine open_icbc

  subroutine read_icbc(ps,ts,u,v,t,qv)
    implicit none
    real(rk8) , pointer , dimension(:,:,:) , intent(out) :: u
    real(rk8) , pointer , dimension(:,:,:) , intent(out) :: v
    real(rk8) , pointer , dimension(:,:,:) , intent(out) :: t
    real(rk8) , pointer , dimension(:,:,:) , intent(out) :: qv
    real(rk8) , pointer , dimension(:,:) , intent(out) :: ps
    real(rk8) , pointer , dimension(:,:) , intent(out) :: ts

    integer(ik4) , dimension(4) :: istart , icount

    istart(3) = ibcrec
    istart(2) = 1
    istart(1) = 1
    icount(3) = 1
    icount(2) = iy
    icount(1) = jx
    istatus = nf90_get_var(ibcin,icbc_ivar(1),ps,istart(1:3),icount(1:3))
    call check_ok(__FILE__,__LINE__,'variable ps read error', 'ICBC FILE')
    istatus = nf90_get_var(ibcin,icbc_ivar(2),ts,istart(1:3),icount(1:3))
    call check_ok(__FILE__,__LINE__,'variable ts read error', 'ICBC FILE')
    istart(4) = ibcrec
    istart(3) = 1
    istart(2) = 1
    istart(1) = 1
    icount(4) = 1
    icount(3) = kz
    icount(2) = iy
    icount(1) = jx
    istatus = nf90_get_var(ibcin,icbc_ivar(3),u,istart,icount)
    call check_ok(__FILE__,__LINE__,'variable u read error', 'ICBC FILE')
    istatus = nf90_get_var(ibcin,icbc_ivar(4),v,istart,icount)
    call check_ok(__FILE__,__LINE__,'variable v read error', 'ICBC FILE')
    istatus = nf90_get_var(ibcin,icbc_ivar(5),t,istart,icount)
    call check_ok(__FILE__,__LINE__,'variable t read error', 'ICBC FILE')
    istatus = nf90_get_var(ibcin,icbc_ivar(6),qv,istart,icount)
    call check_ok(__FILE__,__LINE__,'variable qx read error', 'ICBC FILE')
    if ( ensemble_run ) then
      write(stdout,*) 'Appling perturbation to input dataset:'
      if ( lperturb_ts ) then
        write(stdout,'(a,f7.2,a)') 'TS with value ',perturb_frac_ts*d_100,'%'
        call randify(ts,perturb_frac_ts,jx,iy)
      end if
      if ( lperturb_ps ) then
        write(stdout,'(a,f7.2,a)') 'PS with value ',perturb_frac_ps*d_100,'%'
        call randify(ps,perturb_frac_ps,jx,iy)
      end if
      if ( lperturb_t ) then
        write(stdout,'(a,f7.2,a)') 'T  with value ',perturb_frac_t*d_100,'%'
        call randify(t,perturb_frac_t,jx,iy,kz)
      end if
      if ( lperturb_q ) then
        write(stdout,'(a,f7.2,a)') 'Q  with value ',perturb_frac_q*d_100,'%'
        call randify(qv,perturb_frac_q,jx,iy,kz)
      end if
      if ( lperturb_u ) then
        write(stdout,'(a,f7.2,a)') 'U  with value ',perturb_frac_u*d_100,'%'
        call randify(u,perturb_frac_u,jx,iy,kz)
      end if
      if ( lperturb_v ) then
        write(stdout,'(a,f7.2,a)') 'V  with value ',perturb_frac_v*d_100,'%'
        call randify(v,perturb_frac_v,jx,iy,kz)
      end if
    end if
  end subroutine read_icbc

  subroutine close_icbc
    implicit none
    if (ibcin >= 0) then
      istatus = nf90_close(ibcin)
      call check_ok(__FILE__,__LINE__,'Error Close ICBC file','ICBC FILE')
      if ( allocated(icbc_idate) ) deallocate(icbc_idate)
      ibcin = -1
    end if
  end subroutine close_icbc

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

end module mod_ncio
