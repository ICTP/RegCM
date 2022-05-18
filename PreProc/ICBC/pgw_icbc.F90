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

#ifdef PNETCDF
subroutine myabort
  use mod_stdio
  use mpi
  implicit none
  integer :: ierr
  write(stderr,*) ' Execution terminated because of runtime error'
  call mpi_abort(mpi_comm_self,1,ierr)
end subroutine myabort
#else
subroutine myabort
  implicit none
  stop ' Execution terminated because of runtime error'
end subroutine myabort
#endif

program pgw_icbc
  use mod_intkinds
  use mod_realkinds
  use mod_dynparam , only : npgwlev
  use mod_message
  use mod_header
  use mod_stdio
  use mod_memutil
  use mod_date
  use mod_date
  use netcdf
  use mod_vertint
  use mod_vectutil
  use mod_constants
#ifdef PNETCDF
  use mpi
#endif

  implicit none

  integer(ik4) :: ierr , idynamic
  integer(ik4) :: jx , iy , kz , icheck , jcheck
  integer(ik4) :: pgwin , icbcin
  integer(ik4) :: idimid , ivarid
  integer(ik4) , dimension(7) :: pgw_ivar
  integer(ik4) , dimension(8) :: icbc_ivar
  character(len=256) :: prgname , icbcfile , pgwfile
  real(rkx) , pointer , dimension(:) :: plev , sigma , ak , bk
  real(rkx) , pointer , dimension(:) :: sigmar
  real(rkx) , pointer , dimension(:) :: dsigma , sigmaf
  real(rkx) , pointer , dimension(:,:) :: ps , pd , ts , topo , xmap , ps0
  real(rkx) , pointer , dimension(:,:,:) :: z0 , p0 , t0 , p
  real(rkx) , pointer , dimension(:,:,:) :: u , v , t , q , pp , ww
  real(rkx) , pointer , dimension(:,:) :: bps , ps1 , ps2 , ps3
  real(rkx) , pointer , dimension(:,:) :: bts , ts1 , ts2 , ts3
  real(rkx) , pointer , dimension(:,:,:) :: bu , u1 , u2 , u3
  real(rkx) , pointer , dimension(:,:,:) :: bv , v1 , v2 , v3
  real(rkx) , pointer , dimension(:,:,:) :: bt , t1 , t2 , t3
  real(rkx) , pointer , dimension(:,:,:) :: bq , q1 , q2 , q3
  real(rkx) , pointer , dimension(:,:,:) :: bias
  integer(ik4) :: n , nsteps
  integer(ik4) :: i , j , k
  type(rcm_time_and_date) :: idate
  character(len=64) :: timeunit , timecal
  real(rkx) :: ds , pss , pst , ptop
  real(rkx) , dimension(1) :: nctime
  real(rkx) :: w1 , w2 , w3 , fm , hm
  integer(ik4) :: year , im2 , day , hour , im1 , im3

#ifdef PNETCDF
  call mpi_init(ierr)
#endif

  call header('pgw_icbc')
  !
  ! Read input global namelist
  !
  call get_command_argument(0,value=prgname)
  call get_command_argument(1,value=icbcfile)
  call get_command_argument(2,value=pgwfile)

  ierr = nf90_open(icbcfile, nf90_write, icbcin)
  if ( ierr /= nf90_noerr ) then
    write ( stderr, * ) 'Cannot open icbc file: ', nf90_strerror(ierr)
    call usage( )
  end if
  ierr = nf90_open(pgwfile, nf90_nowrite, pgwin)
  if ( ierr /= nf90_noerr ) then
    write ( stderr, * ) 'Cannot open pgwbc file: ', nf90_strerror(ierr)
    write ( stderr, * ) nf90_strerror(ierr)
    call usage( )
  end if
!
  call memory_init

  idynamic = 1
  ierr = nf90_inq_varid(icbcin, 'a', ivarid)
  if ( ierr == nf90_noerr ) then
    idynamic = 3
  end if
  ierr = nf90_inq_varid(icbcin, 'pp', ivarid)
  if ( ierr == nf90_noerr ) then
    idynamic = 2
  end if

  ierr = nf90_inq_varid(pgwin, 'ps', pgw_ivar(1))
  call check_ok(__FILE__,__LINE__,'variable ps miss', 'PGWBC FILE')
  ierr = nf90_inq_varid(pgwin, 'ts', pgw_ivar(2))
  call check_ok(__FILE__,__LINE__,'variable ts miss', 'PGWBC FILE')
  ierr = nf90_inq_varid(pgwin, 'u', pgw_ivar(3))
  call check_ok(__FILE__,__LINE__,'variable u miss', 'PGWBC FILE')
  ierr = nf90_inq_varid(pgwin, 'v', pgw_ivar(4))
  call check_ok(__FILE__,__LINE__,'variable v miss', 'PGWBC FILE')
  ierr = nf90_inq_varid(pgwin, 't', pgw_ivar(5))
  call check_ok(__FILE__,__LINE__,'variable t miss', 'PGWBC FILE')
  ierr = nf90_inq_varid(pgwin, 'qv', pgw_ivar(6))
  call check_ok(__FILE__,__LINE__,'variable qv miss', 'PGWBC FILE')
  ierr = nf90_inq_varid(pgwin, 'z', pgw_ivar(7))
  call check_ok(__FILE__,__LINE__,'variable z miss', 'PGWBC FILE')

  ierr = nf90_inq_varid(icbcin, 'ps', icbc_ivar(1))
  call check_ok(__FILE__,__LINE__,'variable ps miss', 'ICBC FILE')
  ierr = nf90_inq_varid(icbcin, 'ts', icbc_ivar(2))
  call check_ok(__FILE__,__LINE__,'variable ts miss', 'ICBC FILE')
  ierr = nf90_inq_varid(icbcin, 'u', icbc_ivar(3))
  call check_ok(__FILE__,__LINE__,'variable u miss', 'ICBC FILE')
  ierr = nf90_inq_varid(icbcin, 'v', icbc_ivar(4))
  call check_ok(__FILE__,__LINE__,'variable v miss', 'ICBC FILE')
  ierr = nf90_inq_varid(icbcin, 't', icbc_ivar(5))
  call check_ok(__FILE__,__LINE__,'variable t miss', 'ICBC FILE')
  ierr = nf90_inq_varid(icbcin, 'qv', icbc_ivar(6))
  call check_ok(__FILE__,__LINE__,'variable qv miss', 'ICBC FILE')
  if ( idynamic == 2 ) then
    ierr = nf90_inq_varid(icbcin, 'pp', icbc_ivar(7))
    call check_ok(__FILE__,__LINE__,'variable pp miss', 'ICBC FILE')
    ierr = nf90_inq_varid(icbcin, 'w', icbc_ivar(8))
    call check_ok(__FILE__,__LINE__,'variable w miss', 'ICBC FILE')
  end if

  ierr = nf90_inq_dimid(icbcin, 'iy', idimid)
  call check_ok(__FILE__,__LINE__,'dimension iy miss', 'ICBC FILE')
  ierr = nf90_inquire_dimension(icbcin, idimid, len=iy)
  call check_ok(__FILE__,__LINE__,'dimension iy error', 'ICBC FILE')
  ierr = nf90_inq_dimid(icbcin, 'jx', idimid)
  call check_ok(__FILE__,__LINE__,'dimension jx miss', 'ICBC FILE')
  ierr = nf90_inquire_dimension(icbcin, idimid, len=jx)
  call check_ok(__FILE__,__LINE__,'dimension jx error', 'ICBC FILE')
  ierr = nf90_inq_dimid(icbcin, 'kz', idimid)
  call check_ok(__FILE__,__LINE__,'dimension kz miss', 'ICBC FILE')
  ierr = nf90_inquire_dimension(icbcin, idimid, len=kz)
  call check_ok(__FILE__,__LINE__,'dimension kz error', 'ICBC FILE')
  ierr = nf90_inq_dimid(icbcin, 'time', idimid)
  call check_ok(__FILE__,__LINE__,'dimension time miss', 'ICBC FILE')
  ierr = nf90_inquire_dimension(icbcin, idimid, len=nsteps)
  call check_ok(__FILE__,__LINE__,'dimension time error', 'ICBC FILE')

  ierr = nf90_inq_dimid(pgwin, 'iy', idimid)
  call check_ok(__FILE__,__LINE__,'dimension iy miss', 'PGWBC FILE')
  ierr = nf90_inquire_dimension(pgwin, idimid, len=icheck)
  call check_ok(__FILE__,__LINE__,'dimension iy error', 'PGWBC FILE')
  ierr = nf90_inq_dimid(pgwin, 'jx', idimid)
  call check_ok(__FILE__,__LINE__,'dimension jx miss', 'PGWBC FILE')
  ierr = nf90_inquire_dimension(pgwin, idimid, len=jcheck)
  call check_ok(__FILE__,__LINE__,'dimension jx error', 'PGWBC FILE')

  if ( iy /= icheck .or. jx /= jcheck ) then
    write(stderr,*) 'Dimensions do not match between ICBC and PGWBC files'
    write(stderr,*) 'Please check.'
    call die('pgw_icbc','Dimension mismatch',2)
  end if

  call getmem1d(plev,1,npgwlev,'pgw_icbc:plev')
  call getmem2d(bps,1,jx,1,iy,'pgw_icbc:bps')
  call getmem2d(ps1,1,jx,1,iy,'pgw_icbc:ps1')
  call getmem2d(ps2,1,jx,1,iy,'pgw_icbc:ps2')
  call getmem2d(ps3,1,jx,1,iy,'pgw_icbc:ps3')
  call getmem2d(bts,1,jx,1,iy,'pgw_icbc:bts')
  call getmem2d(ts1,1,jx,1,iy,'pgw_icbc:ts1')
  call getmem2d(ts2,1,jx,1,iy,'pgw_icbc:ts2')
  call getmem2d(ts3,1,jx,1,iy,'pgw_icbc:ts3')
  call getmem3d(bu,1,jx,1,iy,1,npgwlev,'pgw_icbc:bu')
  call getmem3d(u1,1,jx,1,iy,1,npgwlev,'pgw_icbc:u1')
  call getmem3d(u2,1,jx,1,iy,1,npgwlev,'pgw_icbc:u2')
  call getmem3d(u3,1,jx,1,iy,1,npgwlev,'pgw_icbc:u3')
  call getmem3d(bv,1,jx,1,iy,1,npgwlev,'pgw_icbc:bv')
  call getmem3d(v1,1,jx,1,iy,1,npgwlev,'pgw_icbc:v1')
  call getmem3d(v2,1,jx,1,iy,1,npgwlev,'pgw_icbc:v2')
  call getmem3d(v3,1,jx,1,iy,1,npgwlev,'pgw_icbc:v3')
  call getmem3d(bt,1,jx,1,iy,1,npgwlev,'pgw_icbc:bt')
  call getmem3d(t1,1,jx,1,iy,1,npgwlev,'pgw_icbc:t1')
  call getmem3d(t2,1,jx,1,iy,1,npgwlev,'pgw_icbc:t2')
  call getmem3d(t3,1,jx,1,iy,1,npgwlev,'pgw_icbc:t3')
  call getmem3d(bq,1,jx,1,iy,1,npgwlev,'pgw_icbc:bq')
  call getmem3d(q1,1,jx,1,iy,1,npgwlev,'pgw_icbc:q1')
  call getmem3d(q2,1,jx,1,iy,1,npgwlev,'pgw_icbc:q2')
  call getmem3d(q3,1,jx,1,iy,1,npgwlev,'pgw_icbc:q3')
  call getmem3d(bias,1,jx,1,iy,1,kz,'pgw_icbc:bias')

  call getmem1d(sigma,1,kz,'pgw_icbc:sigma')
  call getmem2d(topo,1,jx,1,iy,'pgw_icbc:topo')
  call getmem2d(ps,1,jx,1,iy,'pgw_icbc:ps')
  call getmem2d(ts,1,jx,1,iy,'pgw_icbc:ts')
  call getmem3d(u,1,jx,1,iy,1,kz,'pgw_icbc:u')
  call getmem3d(v,1,jx,1,iy,1,kz,'pgw_icbc:v')
  call getmem3d(t,1,jx,1,iy,1,kz,'pgw_icbc:t')
  call getmem3d(q,1,jx,1,iy,1,kz,'pgw_icbc:q')
  if ( idynamic == 1 ) then
    call getmem1d(sigmar,1,npgwlev,'pgw_icbc:sigmar')
    call getmem2d(pd,1,jx,1,iy,'pgw_icbc:pd')
  end if
  if ( idynamic == 2 ) then
    call getmem1d(sigmaf,1,kz+1,'pgw_icbc:sigmaf')
    call getmem1d(dsigma,1,kz,'pgw_icbc:dsigma')
    call getmem2d(xmap,1,jx,1,iy,'pgw_icbc:xmap')
    call getmem2d(ps0,1,jx,1,iy,'pgw_icbc:ps0')
    call getmem2d(pd,1,jx,1,iy,'pgw_icbc:pd')
    call getmem3d(p,1,jx,1,iy,1,kz,'pgw_icbc:p')
    call getmem3d(p0,1,jx,1,iy,1,kz,'pgw_icbc:p0')
    call getmem3d(t0,1,jx,1,iy,1,kz,'pgw_icbc:t0')
    call getmem3d(pp,1,jx,1,iy,1,kz,'pgw_icbc:pp')
    call getmem3d(ww,1,jx,1,iy,1,kz,'pgw_icbc:ww')
  end if
  if ( idynamic == 3 ) then
    call getmem3d(z0,1,jx,1,iy,1,kz,'pgw_icbc:z0')
    call getmem1d(ak,1,kz,'pgw_icbc:ak')
    call getmem1d(bk,1,kz,'pgw_icbc:bk')
  end if

  ierr = nf90_inq_varid(pgwin, 'plev', ivarid)
  call check_ok(__FILE__,__LINE__,'variable plev miss', 'PGWBC FILE')
  ierr = nf90_get_var(pgwin,ivarid,plev)
  call check_ok(__FILE__,__LINE__,'variable plev read error', 'PGWBC FILE')
  ierr = nf90_inq_varid(icbcin, 'kz', ivarid)
  call check_ok(__FILE__,__LINE__,'variable kz miss', 'ICBC FILE')
  ierr = nf90_get_var(icbcin,ivarid,sigma)
  call check_ok(__FILE__,__LINE__,'variable kz read error', 'ICBC FILE')
  ierr = nf90_inq_varid(icbcin, 'topo', ivarid)
  call check_ok(__FILE__,__LINE__,'variable topo miss', 'ICBC FILE')
  ierr = nf90_get_var(icbcin,ivarid,topo)
  call check_ok(__FILE__,__LINE__,'variable topo read error', 'ICBC FILE')
  if ( idynamic == 3 ) then
    ierr = nf90_inq_varid(icbcin, 'a', ivarid)
    call check_ok(__FILE__,__LINE__,'variable a miss', 'ICBC FILE')
    ierr = nf90_get_var(icbcin,ivarid,ak)
    call check_ok(__FILE__,__LINE__,'variable a read error', 'ICBC FILE')
    ierr = nf90_inq_varid(icbcin, 'b', ivarid)
    call check_ok(__FILE__,__LINE__,'variable b miss', 'ICBC FILE')
    ierr = nf90_get_var(icbcin,ivarid,bk)
    call check_ok(__FILE__,__LINE__,'variable b read error', 'ICBC FILE')
    do k = 1 , kz
      do i = 1 , iy
        do j = 1 , jx
          z0(j,i,k) = topo(j,i)*(bk(k)-1.0_rkx) + ak(k)
        end do
      end do
    end do
  else if ( idynamic == 2 ) then
    ierr = nf90_get_att(icbcin, nf90_global, 'grid_size_in_meters', ds)
    call check_ok(__FILE__,__LINE__,'attribute grid_size_in_meters miss', &
                  'ICBC FILE')
    ierr = nf90_inq_varid(icbcin, 'xmap', ivarid)
    call check_ok(__FILE__,__LINE__,'variable xmap miss', 'ICBC FILE')
    ierr = nf90_get_var(icbcin,ivarid,xmap)
    call check_ok(__FILE__,__LINE__,'variable xmap read error', 'ICBC FILE')
    ierr = nf90_inq_varid(icbcin, 'ps0', ivarid)
    call check_ok(__FILE__,__LINE__,'variable ps0 miss', 'ICBC FILE')
    ierr = nf90_get_var(icbcin,ivarid,ps0)
    call check_ok(__FILE__,__LINE__,'variable ps0 read error', 'ICBC FILE')
    ierr = nf90_inq_varid(icbcin, 'p0', ivarid)
    call check_ok(__FILE__,__LINE__,'variable p0 miss', 'ICBC FILE')
    ierr = nf90_get_var(icbcin,ivarid,p0)
    call check_ok(__FILE__,__LINE__,'variable p0 read error', 'ICBC FILE')
    ierr = nf90_inq_varid(icbcin, 't0', ivarid)
    call check_ok(__FILE__,__LINE__,'variable t0 miss', 'ICBC FILE')
    ierr = nf90_get_var(icbcin,ivarid,t0)
    call check_ok(__FILE__,__LINE__,'variable t0 read error', 'ICBC FILE')
    sigmaf(1) = 0.0_rkx
    sigmaf(kz+1) = 1.0_rkx
    do k = 2 , kz
      sigmaf(k) = (sigma(k-1)+sigma(k))*0.5_rkx
    end do
    do k = 1 , kz
      dsigma(k) = sigmaf(k+1)-sigmaf(k)
    end do
    ierr = nf90_inq_varid(icbcin, 'ptop', ivarid)
    call check_ok(__FILE__,__LINE__,'variable ptop miss', 'ICBC FILE')
    ierr = nf90_get_var(icbcin,ivarid,ptop)
    call check_ok(__FILE__,__LINE__,'variable ptop read error', 'ICBC FILE')
    ps0 = ps0 + ptop*100.0_rkx
  else
    ierr = nf90_inq_varid(icbcin, 'ptop', ivarid)
    call check_ok(__FILE__,__LINE__,'variable ptop miss', 'ICBC FILE')
    ierr = nf90_get_var(icbcin,ivarid,ptop)
    call check_ok(__FILE__,__LINE__,'variable ptop read error', 'ICBC FILE')
  end if
  if ( idynamic == 1 ) then
    pss = (plev(1)-plev(npgwlev))/100.0_rkx
    pst = plev(npgwlev)/100.0_rkx
    do k = 1 , npgwlev
      sigmar(k) = (plev(k)-plev(npgwlev))/(plev(1)-plev(npgwlev))
    end do
  else if ( idynamic == 2 ) then
    plev = plev/100.0_rkx
  end if

  ierr = nf90_inq_varid(icbcin, 'time', ivarid)
  call check_ok(__FILE__,__LINE__,'variable time miss', 'ICBC FILE')
  ierr = nf90_get_att(icbcin, ivarid, 'units', timeunit)
  call check_ok(__FILE__,__LINE__,'variable time units missing', 'ICBC FILE')
  ierr = nf90_get_att(icbcin, ivarid, 'calendar', timecal)
  call check_ok(__FILE__,__LINE__,'variable time calendar missing', 'ICBC FILE')
  ierr = nf90_get_var(icbcin,ivarid,nctime,[1],[1])
  call check_ok(__FILE__,__LINE__,'variable time read error', 'ICBC FILE')
  idate = timeval2date(nctime(1),timeunit,timecal)
  call split_idate(idate,year,im2,day,hour)
  im1 = im2 - 1
  if ( im1 == 0 ) im1 = 12
  im3 = im2 + 1
  if ( im3 == 13 ) im3 = 1

  call read_pgw(im1,ps1,ts1,u1,v1,t1,q1)
  call read_pgw(im2,ps2,ts2,u2,v2,t2,q2)
  call read_pgw(im3,ps3,ts3,u3,v3,t3,q3)

  fm = real(nsteps,rkx)
  hm = fm/2.0_rkx

  write(stdout,* ) 'Initial time is ', tochar(idate)

  do n = 1 , nsteps

    write(stdout,* ) 'Processing time record ', n

    if ( n < int(hm) ) then
      w1 = (hm-n)/fm
      w2 = (hm+n)/fm
      w3 = 0.0_rkx
    else
      w1 = 0.0_rkx
      w2 = (fm-n)/fm
      w3 = (n-hm)/fm
    end if
    bps = w1*ps1 + w2*ps2 + w3*ps3
    bts = w1*ts1 + w2*ts2 + w3*ts3
    bu = w1*u1 + w2*u2 + w3*u3
    bv = w1*v1 + w2*v2 + w3*v3
    bt = w1*t1 + w2*t2 + w3*t3
    bq = w1*q1 + w2*q2 + w3*q3

    call read_icbc(n,ps,ts,u,v,t,q,pp,ww)

    ps = (ps + bps/100.0_rkx)
    ts = ts + bts

    ! interpolate to model levels
    if ( idynamic == 1 ) then
      ps = ps - ptop
      !!!!!!!!!!!!!!!!! ATTENTION !!!!!!!!!!!!!!!!
      !!!!!!!!!!!! ASSUME NON BAND RUN !!!!!!!!!!!
      !!!!!!!!!!!!!!!!! ATTENTION !!!!!!!!!!!!!!!!
      call crs2dot(pd,ps,jx,iy,0,0)
      call intv1(bias,bu,pd,sigma,pss,sigmar,ptop,pst,jx,iy,kz,npgwlev,1)
      u = u + bias
      call intv1(bias,bv,pd,sigma,pss,sigmar,ptop,pst,jx,iy,kz,npgwlev,1)
      v = v + bias
      call intv1(bias,bt,ps,sigma,pss,sigmar,ptop,pst,jx,iy,kz,npgwlev,1)
      t = t + bias
      call intv1(bias,bq,ps,sigma,pss,sigmar,ptop,pst,jx,iy,kz,npgwlev,1)
      q = q + bias
      ps = ps + ptop
    else if ( idynamic == 2 ) then
      p = (p0 + pp)/100.0_rkx
      call intlinreg(bias,bu,ps,plev,1,jx,1,iy,npgwlev,p,kz)
      u = u + bias
      call intlinreg(bias,bv,ps,plev,1,jx,1,iy,npgwlev,p,kz)
      v = v + bias
      call intlinreg(bias,bt,ps,plev,1,jx,1,iy,npgwlev,p,kz)
      t = t + bias
      call intlinreg(bias,bq,ps,plev,1,jx,1,iy,npgwlev,p,kz)
      q = q + bias
      ps = (ps - ptop)/10.0_rkx
      call crs2dot(pd,ps,jx,iy,0,0)
      p = p * 100.0_rkx
      call compute_w(ds,pd,ps,xmap,p,u,v,t,ww)
      ps = (ps * 10.0_rkx + ptop) * 100.0_rkx
      call compute_pp(sigmaf,t,q,ps,ps0,p0,t0,pp)
      ps = ps / 100.0_rkx
    else ! if ( idynamic == 3 ) then

    end if

    call write_icbc(n,ps,ts,u,v,t,q,pp,ww)
  end do

  ierr = nf90_close(pgwin)
  call check_ok(__FILE__,__LINE__,'Error Close PGWBC file','PGWBC FILE')
  ierr = nf90_close(icbcin)
  call check_ok(__FILE__,__LINE__,'Error Close ICBC file','ICBC FILE')

  call memory_destroy

  call finaltime(0)
  write(stdout,*) 'Successfully completed PGWBC'

#ifdef PNETCDF
  call mpi_finalize(ierr)
#endif

  contains

  subroutine usage( )
    write ( stderr, * ) 'Usage : '
    write ( stderr, * ) '          ', trim(prgname), ' icbcfile.nc pgwfile.nc'
    write ( stderr, * ) ' '
    write ( stderr, * ) 'Check arguments.'
    call die('pgw_icbc','Check arguments',1)
  end subroutine usage

  subroutine check_ok(f,l,m1,mf)
    implicit none
    character(*) , intent(in) :: f, m1 , mf
    integer(ik4) , intent(in) :: l
    if (ierr /= nf90_noerr) then
      write (stderr,*) trim(m1)
      write (stderr,*) nf90_strerror(ierr)
      call die(f,trim(mf),l)
    end if
  end subroutine check_ok

  subroutine read_pgw(irec,ps,ts,u,v,t,q)
    implicit none
    integer(ik4) , intent(in) :: irec
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: ps
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: ts
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: u
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: v
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: t
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: q
    !real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: z
    integer(ik4) , dimension(4) :: istart , icount

    istart(1) = 1
    istart(2) = 1
    istart(3) = irec
    icount(1) = jx
    icount(2) = iy
    icount(3) = 1
    ierr = nf90_get_var(pgwin,pgw_ivar(1),ps,istart(1:3),icount(1:3))
    call check_ok(__FILE__,__LINE__,'variable ps read error', 'PGWBC FILE')
    ierr = nf90_get_var(pgwin,pgw_ivar(2),ts,istart(1:3),icount(1:3))
    call check_ok(__FILE__,__LINE__,'variable ts read error', 'PGWBC FILE')
    istart(1) = 1
    istart(2) = 1
    istart(3) = 1
    istart(4) = irec
    icount(1) = jx
    icount(2) = iy
    icount(3) = npgwlev
    icount(4) = 1
    ierr = nf90_get_var(pgwin,pgw_ivar(3),u,istart,icount)
    call check_ok(__FILE__,__LINE__,'variable u read error', 'PGWBC FILE')
    ierr = nf90_get_var(pgwin,pgw_ivar(4),v,istart,icount)
    call check_ok(__FILE__,__LINE__,'variable v read error', 'PGWBC FILE')
    ierr = nf90_get_var(pgwin,pgw_ivar(5),t,istart,icount)
    call check_ok(__FILE__,__LINE__,'variable t read error', 'PGWBC FILE')
    ierr = nf90_get_var(pgwin,pgw_ivar(6),q,istart,icount)
    call check_ok(__FILE__,__LINE__,'variable q read error', 'PGWBC FILE')
    !ierr = nf90_get_var(pgwin,pgw_ivar(7),z,istart,icount)
    !call check_ok(__FILE__,__LINE__,'variable z read error', 'PGWBC FILE')
  end subroutine read_pgw

  subroutine read_icbc(irec,ps,ts,u,v,t,q,pp,ww)
    implicit none
    integer(ik4) , intent(in) :: irec
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: ps
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: ts
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: u
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: v
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: t
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: q
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: pp
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: ww
    integer(ik4) , dimension(4) :: istart , icount

    istart(1) = 1
    istart(2) = 1
    istart(3) = irec
    icount(1) = jx
    icount(2) = iy
    icount(3) = 1
    ierr = nf90_get_var(icbcin,icbc_ivar(1),ps,istart(1:3),icount(1:3))
    call check_ok(__FILE__,__LINE__,'variable ps read error', 'ICBC FILE')
    ierr = nf90_get_var(icbcin,icbc_ivar(2),ts,istart(1:3),icount(1:3))
    call check_ok(__FILE__,__LINE__,'variable ts read error', 'ICBC FILE')
    istart(1) = 1
    istart(2) = 1
    istart(3) = 1
    istart(4) = irec
    icount(1) = jx
    icount(2) = iy
    icount(3) = kz
    icount(4) = 1
    ierr = nf90_get_var(icbcin,icbc_ivar(3),u,istart,icount)
    call check_ok(__FILE__,__LINE__,'variable u read error', 'ICBC FILE')
    ierr = nf90_get_var(icbcin,icbc_ivar(4),v,istart,icount)
    call check_ok(__FILE__,__LINE__,'variable v read error', 'ICBC FILE')
    ierr = nf90_get_var(icbcin,icbc_ivar(5),t,istart,icount)
    call check_ok(__FILE__,__LINE__,'variable t read error', 'ICBC FILE')
    ierr = nf90_get_var(icbcin,icbc_ivar(6),q,istart,icount)
    call check_ok(__FILE__,__LINE__,'variable q read error', 'ICBC FILE')
    if ( idynamic == 2 ) then
      ierr = nf90_get_var(icbcin,icbc_ivar(7),pp,istart,icount)
      call check_ok(__FILE__,__LINE__,'variable pp read error', 'ICBC FILE')
      ierr = nf90_get_var(icbcin,icbc_ivar(8),ww,istart,icount)
      call check_ok(__FILE__,__LINE__,'variable ww read error', 'ICBC FILE')
    end if
  end subroutine read_icbc

  subroutine write_icbc(irec,ps,ts,u,v,t,q,pp,ww)
    implicit none
    integer(ik4) , intent(in) :: irec
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: ps
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: ts
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: u
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: v
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: t
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: q
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: pp
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: ww
    integer(ik4) , dimension(4) :: istart , icount

    istart(1) = 1
    istart(2) = 1
    istart(3) = irec
    icount(1) = jx
    icount(2) = iy
    icount(3) = 1
    ierr = nf90_put_var(icbcin,icbc_ivar(1),ps,istart(1:3),icount(1:3))
    call check_ok(__FILE__,__LINE__,'variable ps read error', 'ICBC FILE')
    ierr = nf90_put_var(icbcin,icbc_ivar(2),ts,istart(1:3),icount(1:3))
    call check_ok(__FILE__,__LINE__,'variable ts read error', 'ICBC FILE')
    istart(1) = 1
    istart(2) = 1
    istart(3) = 1
    istart(4) = irec
    icount(1) = jx
    icount(2) = iy
    icount(3) = kz
    icount(4) = 1
    ierr = nf90_put_var(icbcin,icbc_ivar(3),u,istart,icount)
    call check_ok(__FILE__,__LINE__,'variable u read error', 'ICBC FILE')
    ierr = nf90_put_var(icbcin,icbc_ivar(4),v,istart,icount)
    call check_ok(__FILE__,__LINE__,'variable v read error', 'ICBC FILE')
    ierr = nf90_put_var(icbcin,icbc_ivar(5),t,istart,icount)
    call check_ok(__FILE__,__LINE__,'variable t read error', 'ICBC FILE')
    ierr = nf90_put_var(icbcin,icbc_ivar(6),q,istart,icount)
    call check_ok(__FILE__,__LINE__,'variable q read error', 'ICBC FILE')
    if ( idynamic == 2 ) then
      ierr = nf90_put_var(icbcin,icbc_ivar(7),pp,istart,icount)
      call check_ok(__FILE__,__LINE__,'variable pp read error', 'ICBC FILE')
      ierr = nf90_put_var(icbcin,icbc_ivar(8),ww,istart,icount)
      call check_ok(__FILE__,__LINE__,'variable ww read error', 'ICBC FILE')
    end if
  end subroutine write_icbc

  subroutine compute_w(dx,pd,ps,xm,p,u,v,t,w)
    implicit none
    real(rkx) , intent(in) :: dx
    real(rkx) , dimension(:,:) , pointer , intent(in) :: ps , pd , xm
    real(rkx) , dimension(:,:,:) , pointer , intent(in) :: p , u , v , t
    real(rkx) , dimension(:,:,:) , pointer , intent(inout) :: w
    integer(ik4) :: i , j , k , jp , ip , jm , im , km1 , kp1
    real(rkx) :: ua , ub , va , vb , ubar , vbar , rho , dx2
    real(rkx) , dimension(jx,iy) :: dummy , dummy1
    real(rkx) , dimension(kz) :: mdv
    real(rkx) , dimension(kz+1) :: qdt
    real(rkx) , dimension(jx,iy,kz+1) :: omega

    dx2 = d_two * dx
    dummy = (xm * xm) / dx2
    dummy1 = xm / dx2
    qdt(kz+1) = 0.0_rkx

    do i = 1 , iy
      ip = min(i+1,iy)
      im = max(i-1,1)
      do j = 1 , jx
        jp = min(j+1,jx)
        jm = max(j-1,1)
        do k = 1 , kz
          ua = u(j ,i ,k) * pd(j,i)  + &
               u(j ,ip,k) * pd(j,ip)
          ub = u(jp, i,k) * pd(jp,i) + &
               u(jp,ip,k) * pd(jp,ip)
          va = v(j ,i ,k) * pd(j,i) + &
               v(jp,i ,k) * pd(jp,i)
          vb = v(j ,ip,k) * pd(j,ip) + &
               v(jp,ip,k) * pd(jp,ip)
          mdv(k) = (ub-ua + vb-va) * dummy(j,i) / ps(j,i)
        end do
        do k = kz , 1 , -1
          qdt(k) = qdt(k+1) + mdv(k) * dsigma(k)
        end do
        do k = kz+1 , 1 , -1
          km1 = max(k-1,1)
          kp1 = min(k+1,kz)
          ubar = 0.125_rkx * (u(j ,i ,km1) + u(j ,ip,km1) + &
                              u(jp,i ,km1) + u(jp,ip,km1) + &
                              u(j ,i ,kp1) + u(j ,ip,kp1) + &
                              u(jp,i ,kp1) + u(jp,ip,kp1))
          vbar = 0.125_rkx * (v(j ,i ,km1) + v(j ,ip,km1) + &
                              v(jp,i ,km1) + v(jp,ip,km1) + &
                              v(j ,i ,kp1) + v(j ,ip,kp1) + &
                              v(jp,i ,kp1) + v(jp,ip,kp1))
          omega(j,i,k) = ps(j,i) * qdt(k) + sigmaf(k) * &
                       ((ps(jp,i) - ps(jm,i)) * ubar + &
                        (ps(j,ip) - ps(j,im)) * vbar) * dummy1(j,i)
        end do
      end do
    end do
    call smtdsmt(omega,1,iy,1,jx,1,kz+1)
    do k = 2 , kz+1
      do i = 1 , iy
        do j = 1 , jx
          rho = p(j,i,k-1) / rgas / t(j,i,k-1)
          w(j,i,k-1) = - 1000.0_rkx * omega(j,i,k)/rho * regrav
        end do
      end do
    end do
  end subroutine compute_w

  subroutine smtdsmt(slab,i1,i2,j1,j2,k1,k2)
    implicit none
    integer(ik4) , intent(in) :: i1 , i2 , j1 , j2 , k1 , k2
    real(rkx) , intent(inout) , dimension(j1:j2,i1:i2,k1:k2) :: slab
    real(rkx) :: aplus , asv , cell
    integer(ik4) :: i , is , ie , j , js , je , k , kp , np
    real(rkx) , dimension(2) :: xnu
    integer(ik4) , parameter :: npass = 16
    !
    ! purpose: spatially smooth data in slab to dampen short
    ! wavelength components
    !
    ie = i2-1
    je = j2-1
    is = i1+1
    js = j1+1
    xnu(1) =  0.50_rkx
    xnu(2) = -0.52_rkx
    do k = k1 , k2
      do np = 1 , npass
        do kp = 1 , 2
          ! first smooth in the ni direction
          do i = i1 , i2
            asv = slab(j1,i,k)
            do j = js , je
              cell = slab(j,i,k)
              aplus = slab(j+1,i,k)
              slab(j,i,k) = cell + xnu(kp)*( (asv+aplus)/d_two - cell)
              asv = cell
            end do
          end do
          ! smooth in the nj direction
          do j = j1 , j2
            asv = slab(j,i1,k)
            do i = is , ie
              cell = slab(j,i,k)
              aplus = slab(j,i+1,k)
              slab(j,i,k) = cell + xnu(kp)*((asv+aplus)/d_two - cell)
              asv = cell
            end do
          end do
        end do
        slab(j1,:,k) = slab(j1+1,:,k)
        slab(:,i1,k) = slab(:,i1+1,k)
        slab(j2-1,:,k) = slab(j2-2,:,k)
        slab(:,i2-1,k) = slab(:,i2-2,k)
        slab(j2,:,k) = slab(j2-1,:,k)
        slab(:,i2,k) = slab(:,i2-1,k)
      end do
    end do
  end subroutine smtdsmt

  subroutine compute_pp(sigma,t,q,ps,ps0,p0,t0,pp)
    implicit none
    real(rkx) , pointer , intent(in) , dimension(:) :: sigma
    real(rkx) , pointer , intent(in) , dimension(:,:) :: ps , ps0
    real(rkx) , pointer , intent(in) , dimension(:,:,:) :: p0 , t0
    real(rkx) , pointer , intent(in) , dimension(:,:,:) :: t , q
    real(rkx) , pointer , intent(inout) , dimension(:,:,:) :: pp
    real(rkx) :: aa , bb , cc , delp0
    real(rkx) :: psp , tk , tkp1 , tvk , tvkp1 , tvpot , wtl , wtu
    integer :: i , j , k
    do i = 1 , iy
      do j = 1 , jx
        psp = ps(j,i) - ps0(j,i)
        delp0 = ps0(j,i) - p0(j,i,kz)
        tvk = t(j,i,kz) * (1.0_rkx + ep1*q(j,i,kz))
        tk = t(j,i,kz)
        tvpot = (tvk - t0(j,i,kz)) / tk
        pp(j,i,kz) = (tvpot*delp0 + psp) / (d_one + delp0/p0(j,i,kz))
        do k = kz-1 , 1 , -1
          tvkp1 = t(j,i,k+1) * (1.0_rkx + ep1*q(j,i,k+1))
          tvk = t(j,i,k) * (1.0_rkx + ep1*q(j,i,k))
          tkp1 = t(j,i,k+1)
          tk = t(j,i,k)
          wtl = (sigma(k+1) - sigma(k)) / (sigma(k+2) - sigma(k))
          wtu = d_one - wtl
          aa = egrav / (p0(j,i,k+1) - p0(j,i,k))
          bb = egrav * wtl / p0(j,i,k+1) * t0(j,i,k+1) / tkp1
          cc = egrav * wtu / p0(j,i,k  ) * t0(j,i,k  ) / tk
          tvpot = wtl * ((tvkp1 - t0(j,i,k+1)) / tkp1) + &
                  wtu * ((tvk   - t0(j,i,k  )) / tk  )
          pp(j,i,k) = (egrav * tvpot + pp(j,i,k+1) * (aa - bb)) / (aa + cc)
        end do
      end do
    end do
  end subroutine compute_pp

end program pgw_icbc
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
