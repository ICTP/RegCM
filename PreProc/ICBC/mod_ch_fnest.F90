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

module mod_ch_fnest
  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_constants
  use mod_stdio
  use mod_memutil
  use mod_grid
  use mod_wrtoxd
  use mod_kdinterp
  use mod_date
  use mod_nchelper
  use mod_posix
  use mod_message
  use mod_vertint
  use netcdf

  implicit none

  private

  integer(ik4) :: fchem , nrec
  integer(ik4) , dimension(:) , allocatable :: ncid

  ! Remember to increase this if we get more chemistry variables !
  integer(ik4) , parameter :: maxchvar = 64
  character(len=8) , dimension(maxchvar) :: chnames
  character(len=256) :: cbase

  integer(ik4) :: iy_in , jx_in , kz_in , np

  real(rkx) , dimension(:,:,:) , pointer :: mxc , mxcp , mxcp4
  real(rkx) , dimension(:,:,:,:) , pointer :: mxc4_1
  real(rkx) , dimension(:,:,:) , pointer :: pp3d , p3d
  real(rkx) , dimension(:,:) , pointer :: xlat_in , xlon_in , ht_in
  real(rkx) , dimension(:,:) , pointer :: p0_in , pstar0 , ps , xps , xps3
  real(rkx) , dimension(:) , pointer :: sigma_in , plev , sigmar
  real(rkx) :: pss
  integer(ik4) :: oidyn
  character(len=6) :: iproj_in
  real(rkx) :: clat_in , clon_in , ds_in
  real(rkx) :: plat_in , plon_in , ptop_in , xcone_in
  type(rcm_time_and_date) , dimension(:) , pointer :: itimes
  real(rkx) , dimension(:) , pointer :: xtimes
  character(len=64) :: timeunits , timecal
  integer(ik4) :: ncicbc , ivarps , irec
  type(rcm_time_and_date) :: iodate
  type(h_interpolator) :: hint
  logical , dimension(3) :: mapping

  data mapping /.false.,.false.,.false./
  data ncicbc /-1/

  public :: init_fnest , get_fnest , close_fnest

  contains

  subroutine init_fnest(idate,cdir,cname,dochem,dooxcl,doaero)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate
    character(len=*) , intent(in) :: cdir , cname
    logical , intent(in) :: dochem , dooxcl , doaero
    type(direntry) , pointer , dimension(:) :: listf => null()
    type(rcm_time_and_date) :: imf
    character(len=11) :: cdate
    character(len=256) :: fname , icbcfilename
    integer(ik4) :: fnum , nf , is , ie , ip , i , j , k
    real(rkx) , dimension(2) :: trlat
    real(rkx) :: xsign
    integer(ik4) :: istatus , idimid , ivarid
    integer(ik4) , dimension(3) :: istart , icount

    if ( dochem ) mapping(1) = .true.
    if ( dooxcl ) mapping(2) = .true.
    if ( doaero ) mapping(3) = .true.

    iodate = idate

    call dirlist(cdir,listf)
    if ( .not. associated(listf) ) then
      call die('chem_icbc','Cannot read from coarse '//trim(cdir),1)
    end if
    fnum = size(listf)
    fchem = 0
    imf = monfirst(idate)
    write(cdate,'(a)') tochar10(imf)
    write (icbcfilename,'(a,a,a,a,a,a)') trim(dirglob), pthsep, &
           trim(domname), '_ICBC.', trim(tochar10(idate)), '.nc'
    do nf = 1 , fnum
      if ( index(listf(nf)%ename,trim(cname)) /= 0 .and. &
           index(listf(nf)%ename,'.nc') /= 0 .and. &
           index(listf(nf)%ename,trim(cdate)) /= 0 .and. &
           index(listf(nf)%ename,'ATM') == 0 .and. &
           index(listf(nf)%ename,'SRF') == 0 .and. &
           index(listf(nf)%ename,'STS') == 0 .and. &
           index(listf(nf)%ename,'SAV') == 0 .and. &
           index(listf(nf)%ename,'OPT') == 0 .and. &
           index(listf(nf)%ename,'RAD') == 0 .and. &
           index(listf(nf)%ename,'SUB') == 0 .and. &
           index(listf(nf)%ename,'LAK') == 0 ) then
        fchem = fchem + 1
        is = len_trim(cname) + 2
        ie = index(listf(nf)%ename,trim(cdate)) - 2
        j = 1
        chnames(fchem)(:) = ' '
        do i = is , ie
          chnames(fchem)(j:j) = listf(nf)%ename(i:i)
          j = j + 1
        end do
      end if
    end do
    if ( fchem == 0 ) then
      call die('chem_icbc','No chemistry coarse files available.',1)
    end if
    allocate(ncid(fchem))
    write(cbase,'(a,a1,a,a1)') trim(cdir),pthsep,trim(cname),'_'
    do nf = 1 , fchem
      write(fname,'(a,a,a1,a,a)') trim(cbase), trim(chnames(nf)), &
                  '.', trim(cdate), '.nc'
      istatus = nf90_open(fname,nf90_nowrite,ncid(nf))
      call checkncerr(istatus,__FILE__,__LINE__, &
              'Error open file chemical '//trim(fname))
    end do
    istatus = nf90_inq_dimid(ncid(1),'jx',idimid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error find dim jx')
    istatus = nf90_inquire_dimension(ncid(1),idimid,len=jx_in)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error inquire dim jx')
    istatus = nf90_inq_dimid(ncid(1),'iy',idimid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error find dim iy')
    istatus = nf90_inquire_dimension(ncid(1),idimid,len=iy_in)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error inquire dim iy')
    istatus = nf90_inq_dimid(ncid(1),'kz',idimid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error find dim kz')
    istatus = nf90_inquire_dimension(ncid(1),idimid,len=kz_in)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error inquire dim kz')
    istatus = nf90_inq_dimid(ncid(1), 'time', idimid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Dimension time missing')
    istatus = nf90_inquire_dimension(ncid(1), idimid, len=nrec)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Dimension time read error')
    istatus = nf90_inq_varid(ncid(1), 'time', ivarid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'variable time missing')
    istatus = nf90_get_att(ncid(1), ivarid, 'units', timeunits)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'variable time units missing')
    istatus = nf90_get_att(ncid(1), ivarid, 'calendar', timecal)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'variable time calendar missing')
    call getmem1d(itimes,1,nrec,'mod:nest:itimes')
    call getmem1d(xtimes,1,nrec,'mod:nest:xtimes')
    istatus = nf90_get_var(ncid(1), ivarid, xtimes)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'variable time read error')
    do i = 1 , nrec
      itimes(i) = timeval2date(xtimes(i), timeunits, timecal)
    end do

    istatus = nf90_open(icbcfilename,nf90_nowrite, ncicbc)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error open ICBC file '//trim(icbcfilename))
    istatus = nf90_inq_varid(ncicbc,'ps',ivarps)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error find var ps in icbc file '//trim(icbcfilename))
    irec = 1

    call getmem2d(xlat_in,1,jx_in,1,iy_in,'init_fnest:xlat_in')
    call getmem2d(xlon_in,1,jx_in,1,iy_in,'init_fnest:xlon_in')
    call getmem2d(ht_in,1,jx_in,1,iy_in,'init_fnest:xlon_in')
    call getmem1d(sigma_in,1,kz_in,'init_fnest:sigma_in')
    call getmem3d(mxc,1,jx_in,1,iy_in,1,kz_in,'init_fnest:mxc')

    istatus = nf90_inq_varid(ncid(1), 'kz', ivarid)
    if ( istatus /= nf90_noerr ) then
      istatus = nf90_inq_varid(ncid(1), 'sigma', ivarid)
    end if
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'variable sigma or kz error')
    istatus = nf90_get_var(ncid(1), ivarid, sigma_in)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'variable sigma or kz read error')
    istatus = nf90_inq_varid(ncid(1), 'xlat', ivarid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'variable xlat error')
    istatus = nf90_get_var(ncid(1), ivarid, xlat_in)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'variable xlat read error')
    istatus = nf90_inq_varid(ncid(1), 'xlon', ivarid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'variable xlon error')
    istatus = nf90_get_var(ncid(1), ivarid, xlon_in)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'variable xlon read error')
    istatus = nf90_inq_varid(ncid(1), 'topo', ivarid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'variable topo error')
    istatus = nf90_get_var(ncid(1), ivarid, ht_in)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'variable topo read error')
    istatus = nf90_inq_varid(ncid(1), 'ptop', ivarid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'variable ptop error')
    istatus = nf90_get_var(ncid(1), ivarid, ptop_in)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'variable ptop read error')
    istatus = nf90_get_att(ncid(1), nf90_global, &
                      'grid_size_in_meters', ds_in)
    ds_in = ds_in * sqrt(d_two) * d_r1000
    istatus = nf90_get_att(ncid(1), nf90_global,'projection', iproj_in)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'attribure iproj read error')
    istatus = nf90_get_att(ncid(1), nf90_global, &
                      'latitude_of_projection_origin', clat_in)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'attribure clat read error')
    istatus = nf90_get_att(ncid(1), nf90_global, &
                      'longitude_of_projection_origin', clon_in)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'attribure clat read error')
    istatus = nf90_get_att(ncid(1), nf90_global, 'dynamical_core', oidyn)
    if ( istatus /= nf90_noerr ) then
      oidyn = 1 ! Assume non-hydrostatic
    end if
    if ( iproj_in == 'LAMCON' ) then
      istatus = nf90_get_att(ncid(1), nf90_global, 'standard_parallel', trlat)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'attribure truelat read error')
      if ( clat_in < 0. ) then
        xsign = -1.0_rkx       ! SOUTH HEMESPHERE
      else
        xsign = 1.0_rkx        ! NORTH HEMESPHERE
      end if
      if ( abs(trlat(1)-trlat(2)) > 1.E-1 ) then
        xcone_in = real((log10(cos(trlat(1)*degrad))                       &
                    -log10(cos(trlat(2)*degrad))) /                        &
                    (log10(tan((45.0_rk8-xsign*trlat(1)/2.0_rk8)*degrad))  &
                    -log10(tan((45.0_rk8-xsign*trlat(2)/2.0_rk8)*degrad))),rkx)
      else
        xcone_in = real(xsign*sin(real(trlat(1),rk8)*degrad),rkx)
      end if
    else if ( iproj_in == 'POLSTR' ) then
      xcone_in = 1.0_rkx
    else if ( iproj_in == 'NORMER' ) then
      xcone_in = 0.0_rkx
    else
      istatus = nf90_get_att(ncid(1), nf90_global, &
                      'grid_north_pole_latitude', plat_in)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'attribure plat read error')
      istatus = nf90_get_att(ncid(1), nf90_global, &
                      'grid_north_pole_longitude', plon_in)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'attribure plon read error')
      xcone_in = 0.0_rkx
    end if

    if ( ptop_in > ptop*10.0_rkx ) then
      write(stderr,*) 'WARNING : top pressure higher than PTOP detected.'
      write(stderr,*) 'WARNING : Extrapolation will be performed.'
    end if

    ptop_in = ptop_in * d_100

    np = kz_in - 1
    call getmem1d(plev,1,np,'mod_nest:plev')
    call getmem1d(sigmar,1,np,'mod_nest:sigmar')

    call getmem2d(p0_in,1,jx_in,1,iy_in,'mod_nest:p0_in')
    call getmem2d(pstar0,1,jx_in,1,iy_in,'mod_nest:pstar0')
    call getmem2d(ps,1,jx_in,1,iy_in,'mod_nest:ps')

    if ( oidyn == 2 ) then
      call getmem3d(pp3d,1,jx_in,1,iy_in,1,kz_in,'mod_nest:pp3d')
      call getmem3d(p3d,1,jx_in,1,iy_in,1,kz_in,'mod_nest:p3d')
      istatus = nf90_inq_varid(ncid(1), 'p0', ivarid)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'variable p0 error')
      istatus = nf90_get_var(ncid(1), ivarid, p0_in)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'variable p0 read error')
      pstar0 = p0_in - ptop_in
    else
      istatus = nf90_inq_varid(ncid(1), 'ps', ivarid)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'variable ps error')
      istart(:) = 1
      icount(1) = jx_in
      icount(2) = iy_in
      icount(3) = 1
      istatus = nf90_get_var(ncid(1), ivarid, p0_in, istart, icount)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'variable ps read error')
      pstar0 = p0_in - ptop_in
    end if

    do ip = 1 , np/2
      plev(ip) = d_half * (minval(pstar0*sigma_in(ip+1)) + &
                           maxval(pstar0*sigma_in(ip))) + ptop_in
    end do
    do ip = np/2+1 , np
      plev(ip) = maxval(pstar0*sigma_in(ip+1)) + ptop_in
    end do

    call h_interpolator_create(hint,xlat_in,xlon_in,xlat,xlon,ds_in)

    call getmem3d(mxcp,1,jx_in,1,iy_in,1,np,'init_fnest:mxcp')
    call getmem3d(mxcp4,1,jx,1,iy,1,np,'init_fnest:mxcp4')
    call getmem4d(mxc4_1,1,jx,1,iy,1,kz,1,fchem,'init_fnest:mxc4_1')
    call getmem2d(xps,1,jx,1,iy,'mod_nest:xps')
    call getmem2d(xps3,1,jx,1,iy,'mod_nest:xps3')

    do k = 1 , np
      sigmar(k) = plev(k)/plev(np)
    end do
    pss = plev(np)

  end subroutine init_fnest

  subroutine get_fnest(idate)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate
    character(len=11) :: cdate
    character(len=256) :: fname , icbcfilename
    type(rcm_time_and_date) :: imf
    integer(ik4) :: nf , i , j , k , l , crec , k0
    integer(ik4) :: istatus , idimid , ivarid
    integer(ik4) , dimension(4) :: istart , icount
    real(rkx) :: prcm , pmpi , pmpj , wt1 , wt2
    character(len=8) :: specname

    if ( idate > itimes(nrec) ) then
      imf = monfirst(idate)
      write(cdate,'(a)') tochar10(imf)
      do nf = 1 , fchem
        istatus = nf90_close(ncid(nf))
        write(fname,'(a,a,a1,a,a)') trim(cbase), trim(chnames(nf)), &
                  '.', trim(cdate), '.nc'
        istatus = nf90_open(fname,nf90_nowrite,ncid(nf))
        call checkncerr(istatus,__FILE__,__LINE__, &
                'Error open file chemical '//trim(fname))
      end do
      istatus = nf90_inq_dimid(ncid(1), 'time', idimid)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Dimension time missing')
      istatus = nf90_inquire_dimension(ncid(1), idimid, len=nrec)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Dimension time read error')
      istatus = nf90_inq_varid(ncid(1), 'time', ivarid)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'variable time missing')
      istatus = nf90_get_att(ncid(1), ivarid, 'units', timeunits)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'variable time units missing')
      istatus = nf90_get_att(ncid(1), ivarid, 'calendar', timecal)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'variable time calendar missing')
      call getmem1d(itimes,1,nrec,'mod:nest:itimes')
      call getmem1d(xtimes,1,nrec,'mod:nest:xtimes')
      istatus = nf90_get_var(ncid(1), ivarid, xtimes)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'variable time read error')
      do i = 1 , nrec
        itimes(i) = timeval2date(xtimes(i), timeunits, timecal)
      end do
    else if ( idate < itimes(1) ) then
      imf = prevmon(idate)
      write(cdate,'(a)') tochar10(imf)
      do nf = 1 , fchem
        istatus = nf90_close(ncid(nf))
        write(fname,'(a,a,a1,a,a)') trim(cbase), trim(chnames(nf)), &
                  '.', trim(cdate), '.nc'
        istatus = nf90_open(fname,nf90_nowrite,ncid(nf))
        call checkncerr(istatus,__FILE__,__LINE__, &
                'Error open file chemical '//trim(fname))
      end do
      istatus = nf90_inq_dimid(ncid(1), 'time', idimid)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Dimension time missing')
      istatus = nf90_inquire_dimension(ncid(1), idimid, len=nrec)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Dimension time read error')
      istatus = nf90_inq_varid(ncid(1), 'time', ivarid)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'variable time missing')
      istatus = nf90_get_att(ncid(1), ivarid, 'units', timeunits)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'variable time units missing')
      istatus = nf90_get_att(ncid(1), ivarid, 'calendar', timecal)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'variable time calendar missing')
      call getmem1d(itimes,1,nrec,'mod:nest:itimes')
      call getmem1d(xtimes,1,nrec,'mod:nest:xtimes')
      istatus = nf90_get_var(ncid(1), ivarid, xtimes)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'variable time read error')
      do i = 1 , nrec
        itimes(i) = timeval2date(xtimes(i), timeunits, timecal)
      end do
    end if

    if (.not. lsamemonth(idate, iodate) ) then
      istatus = nf90_close(ncicbc)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error close ICBC file')
      write (icbcfilename,'(a,a,a,a,a,a)') trim(dirglob), pthsep, &
             trim(domname), '_ICBC.', trim(tochar10(idate)), '.nc'
      istatus = nf90_open(icbcfilename,nf90_nowrite, ncicbc)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error open ICBC file '//trim(icbcfilename))
      istatus = nf90_inq_varid(ncicbc,'ps',ivarps)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find var ps in icbc file '//trim(icbcfilename))
      iodate = idate
      irec = 1
    end if

    crec = -1
    do i = 1 , nrec
      if (idate == itimes(i)) then
        crec = i
        exit
      end if
    end do

    if ( crec < 0 ) then
      call die('chem_icbc','Cannot find timestep '//tochar(idate),1)
    end if

    istart(3) = crec
    istart(2) = 1
    istart(1) = 1
    icount(3) = 1
    icount(2) = iy_in
    icount(1) = jx_in
    istatus = nf90_inq_varid(ncid(1), 'ps', ivarid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'variable ps missing')
    istatus = nf90_get_var(ncid(1), ivarid, ps, istart(1:3), icount(1:3))
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'variable ps read error')
    call h_interpolate_cont(hint,ps,xps)

    istart(1) = 1
    istart(2) = 1
    istart(3) = irec
    icount(1) = jx
    icount(2) = iy
    icount(3) = 1
    istatus = nf90_get_var(ncicbc,ivarps,xps3,istart(1:3),icount(1:3))
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read icbc var ps')
    irec = irec + 1

    istart(4) = crec
    istart(3) = 1
    istart(2) = 1
    istart(1) = 1
    icount(4) = 1
    icount(3) = kz_in
    icount(2) = iy_in
    icount(1) = jx_in

    if ( oidyn == 2 ) then
      istatus = nf90_inq_varid(ncid(1), 'ppa', ivarid)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'variable ppa missing')
      istatus = nf90_get_var(ncid(1), ivarid, pp3d, istart, icount)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'variable ppa read error')
      do k = 1 , kz_in
        do i = 1 , iy_in
          do j = 1 , jx_in
            p3d(j,i,k) = pstar0(j,i) * sigma_in(k) + &
                         ptop_in + pp3d(j,i,k)
          end do
        end do
      end do
    end if
    do nf = 1 , fchem
      istatus = nf90_inq_varid(ncid(nf), 'mixrat', ivarid)
      istatus = nf90_get_var(ncid(nf), ivarid, mxc, istart, icount)
      call checkncerr(istatus,__FILE__,__LINE__, &
             'variable mixrat for '//trim(chnames(nf))//' read error')
      if ( oidyn == 1 ) then
        call intlin(mxcp,mxc,ps,sigma_in,ptop_in,jx_in,iy_in,kz_in,plev,np)
      else
        call intlin(mxcp,mxc,p3d,jx_in,iy_in,kz_in,plev,np)
      end if
      call h_interpolate_cont(hint,mxcp,mxcp4)
      do i = 1 , iy
        do j = 1 , jx
          do l = 1 , kz
            prcm = ((xps3(j,i)*0.1_rkx-ptop)*sigmah(l)+ptop)*1000.0_rkx
            k0 = -1
            do k = np , 1 , -1
              pmpi = plev(k)
              k0 = k
              if (prcm > pmpi) exit
            end do
            if ( k0 == np ) then
              pmpj = plev(np-1)
              pmpi = plev(np)
              mxc4_1(j,i,l,nf) = max(d_zero,mxcp4(j,i,np) + &
                 (mxcp4(j,i,np-1) - mxcp4(j,i,np)) * (prcm-pmpi)/(pmpi-pmpj))
            else if (k0 >= 1) then
              pmpj = plev(k0)
              pmpi = plev(k0+1)
              wt1 = (prcm-pmpj)/(pmpi-pmpj)
              wt2 = 1.0_rkx - wt1
              mxc4_1(j,i,l,nf) = mxcp4(j,i,k0+1)*wt1 + mxcp4(j,i,k0)*wt2
            end if
          end do
        end do
      end do
    end do
    !
    ! Now we have to map....
    !
    if ( mapping(1) ) then
      do i = 1 , ncbmz
        chv4(:,:,:,i) = mxc4_1(:,:,:,findex(cbmzspec(i)))
      end do
      call write_ch_icbc(idate)
    end if
    if ( mapping(2) ) then
      do i = 1 , noxsp
        oxv4(:,:,:,i) = mxc4_1(:,:,:,findex(oxspec(i)))
      end do
      call write_ox_icbc(idate)
    end if
    if ( mapping(3) ) then
      do i = 1 , naesp
        if ( aespec(i) == 'SOA' ) cycle
        if ( aespec(i) == 'SSLT03' ) cycle
        if ( aespec(i) == 'SSLT04' ) cycle
        if ( aespec(i)(1:3) == 'DST' ) then
          specname = 'DUST'//aespec(i)(4:5)
        else if ( aespec(i)(1:3) == 'D12' ) then
          specname = 'DUST'//aespec(i)(4:5)
        else if ( aespec(i)(1:3) == 'OC1' ) then
          specname = 'OC_HB'
        else if ( aespec(i)(1:3) == 'OC2' ) then
          specname = 'OC_HL'
        else if ( aespec(i)(1:3) == 'CB1' ) then
          specname = 'BC_HB'
        else if ( aespec(i)(1:3) == 'CB2' ) then
          specname = 'BC_HL'
        else
          specname = aespec(i)
        end if
        aev4(:,:,:,i) = mxc4_1(:,:,:,findex(specname))
      end do
      call write_ae_icbc(idate)
    end if
  end subroutine get_fnest

  integer(ik4) function findex(sname) result(i)
    implicit none
    character(len=*) , intent(in) :: sname
    integer(ik4) :: nf
    i = -1
    do nf = 1 , fchem
      if ( trim(chnames(nf)) == trim(sname) ) then
        i = nf
        return
      end if
    end do
    if ( i == -1 ) then
      call die('ch_fnest','Cannot find specie : '//trim(sname)//' in coarse',1)
    end if
  end function findex

  subroutine close_fnest
    use netcdf
    implicit none
    integer(ik4) :: istatus , fid
    if ( allocated(ncid) ) then
      do fid = 1 , size(ncid)
        if ( ncid(fid) > 0 ) then
          istatus = nf90_close(ncid(fid))
          call checkncerr(istatus,__FILE__,__LINE__, &
                          'Error close fnest file')
        end if
      end do
      deallocate(ncid)
    end if
    if ( ncicbc > 0 ) then
      istatus = nf90_close(ncicbc)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error close icbc file')
    end if
  end subroutine close_fnest

end module mod_ch_fnest
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
