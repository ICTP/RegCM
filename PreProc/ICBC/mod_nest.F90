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

module mod_nest

  use mod_intkinds
  use mod_realkinds
  use mod_stdio
  use mod_dynparam
  use mod_constants
  use mod_date
  use mod_grid
  use mod_write
  use mod_interp
  use mod_vertint
  use mod_hgt
  use mod_humid
  use mod_mksst
  use mod_uvrot
  use mod_vectutil
  use mod_message
  use mod_memutil
  use mod_nchelper

  private

  public :: get_nest , headernest

  character(len=256) , public :: coarsedir , coarsedom

  integer(ik4) , parameter :: np = 15

  integer(ik4) :: nrec

  real(rk8) , pointer , dimension(:,:,:) :: b3
  real(rk8) , pointer , dimension(:,:,:) :: d3
  real(rk8) , pointer , dimension(:,:,:) :: z1

  real(rk8) , pointer , dimension(:,:,:) :: b2
  real(rk8) , pointer , dimension(:,:,:) :: d2

  real(rk8) , pointer , dimension(:,:,:) :: q , t
  real(rk8) , pointer , dimension(:,:,:) :: u , v
  real(rk8) , pointer , dimension(:,:) :: ps
  real(rk8) , pointer , dimension(:,:) :: ht_in
  real(rk8) , pointer , dimension(:,:) :: xlat_in , xlon_in

  real(rk8) , pointer , dimension(:,:,:) :: h3 , q3 , t3
  real(rk8) , pointer , dimension(:,:,:) :: u3 , v3
 
  real(rk8) , pointer , dimension(:,:,:) :: hp , qp , tp
  real(rk8) , pointer , dimension(:,:,:) :: up , vp

  real(rk8) , dimension(np) :: plev , sigmar
  real(rk8) , pointer , dimension(:) :: sig

  integer(ik4) :: iy_in , jx_in , kz_in
  character(len=6) :: iproj_in
  real(rk8) :: clat_in , clon_in , plat_in , plon_in , ptop_in , xcone_in
!
  character(len=14) :: fillin
  character(len=256) :: inpfile
!
  integer(ik4) :: ncinp
  type(rcm_time_and_date) , dimension(:) , pointer :: itimes
  real(rk8) , dimension(:) , pointer :: xtimes
  character(len=64) :: timeunits , timecal
!
  contains

  subroutine get_nest(idate)
    use netcdf
    implicit none

    type(rcm_time_and_date) , intent(in) :: idate

    integer(ik4) :: i , istatus , ivarid , idimid , irec
    integer(ik4) , dimension(4) :: istart , icount
    type(rcm_time_and_date) :: imf
    logical :: lspch

    if (.not. associated(b2)) then
      call die('get_nest','Called get_nest before headernest !',1)
    end if

    if ( idate > itimes(nrec) ) then
      istatus = nf90_close(ncinp)
      call checkncerr(istatus,__FILE__,__LINE__, 'Error close')
      imf = monfirst(idate)
      write (fillin,'(a,i10)') 'ATM.', toint10(imf)
      if ( coarsedir(1:5) == '     ' ) then
        inpfile = trim(inpglob)//pthsep//'RegCM'//pthsep
      else
        inpfile = trim(coarsedir)//pthsep
      end if
      if ( coarsedom(1:5) == '     ' ) then
        inpfile = trim(inpfile)//fillin//'.nc'
      else
        inpfile = trim(inpfile)//trim(coarsedom)//'_'//fillin//'.nc'
      end if
      istatus = nf90_open(inpfile, nf90_nowrite, ncinp)
      call checkncerr(istatus,__FILE__,__LINE__,'Error opening '//trim(inpfile))
      istatus = nf90_inq_dimid(ncinp, 'time', idimid)
      call checkncerr(istatus,__FILE__,__LINE__,'Dimension time missing')
      istatus = nf90_inquire_dimension(ncinp, idimid, len=nrec)
      call checkncerr(istatus,__FILE__,__LINE__,'Dimension time read error')
      istatus = nf90_inq_varid(ncinp, 'time', ivarid)
      call checkncerr(istatus,__FILE__,__LINE__,'variable time missing')
      call getmem1d(itimes,1,nrec,'mod_nest:itimes')
      call getmem1d(xtimes,1,nrec,'mod_nest:xtimes')
      istatus = nf90_get_var(ncinp, ivarid, xtimes)
      call checkncerr(istatus,__FILE__,__LINE__,'variable time read error')
      do i = 1 , nrec
        itimes(i) = timeval2date(xtimes(i), timeunits,timecal)
      end do
    end if

    irec = -1
    do i = 1 , nrec
      if (idate == itimes(i)) then
        irec = i
        exit
      end if
    end do
    if (irec < 0) then
      write (stderr,*) 'Error : time ', tochar(idate), ' not in file'
      call die('get_nest')
    end if

    istart(4) = irec
    istart(3) = 1
    istart(2) = 1
    istart(1) = 1
    icount(4) = 1
    icount(3) = kz_in
    icount(2) = iy_in
    icount(1) = jx_in
    istatus = nf90_inq_varid(ncinp, 'ua', ivarid)
    if ( istatus /= nf90_noerr ) then
      istatus = nf90_inq_varid(ncinp, 'u', ivarid)
      call checkncerr(istatus,__FILE__,__LINE__,'variable u/ua missing')
    end if
    istatus = nf90_get_var(ncinp, ivarid, u, istart, icount)
    call checkncerr(istatus,__FILE__,__LINE__,'variable u/ua read error')
    istatus = nf90_inq_varid(ncinp, 'va', ivarid)
    if ( istatus /= nf90_noerr ) then
      istatus = nf90_inq_varid(ncinp, 'v', ivarid)
      call checkncerr(istatus,__FILE__,__LINE__,'variable v/va missing')
    end if
    istatus = nf90_get_var(ncinp, ivarid, v, istart, icount)
    call checkncerr(istatus,__FILE__,__LINE__,'variable v/va read error')
    istatus = nf90_inq_varid(ncinp, 'ta', ivarid)
    if ( istatus /= nf90_noerr ) then
      istatus = nf90_inq_varid(ncinp, 't', ivarid)
      call checkncerr(istatus,__FILE__,__LINE__,'variable t/ta missing')
    end if
    istatus = nf90_get_var(ncinp, ivarid, t, istart, icount)
    call checkncerr(istatus,__FILE__,__LINE__,'variable t/ta read error')
    lspch = .true.
    istatus = nf90_inq_varid(ncinp, 'qas', ivarid)
    if ( istatus /= nf90_noerr ) then
      lspch = .false.
      istatus = nf90_inq_varid(ncinp, 'qv', ivarid)
      call checkncerr(istatus,__FILE__,__LINE__,'variable qv/qas missing')
    end if
    istatus = nf90_get_var(ncinp, ivarid, q, istart, icount)
    call checkncerr(istatus,__FILE__,__LINE__,'variable qv read error')
    ! Transform specific humidity in mixing ratio
    if ( lspch ) then
      q = q/(d_one-q)
    end if
    istart(3) = irec
    istart(2) = 1
    istart(1) = 1
    icount(3) = 1
    icount(2) = iy_in
    icount(1) = jx_in
    istatus = nf90_inq_varid(ncinp, 'ps', ivarid)
    call checkncerr(istatus,__FILE__,__LINE__,'variable ps missing')
    istatus = nf90_get_var(ncinp, ivarid, ps, istart(1:3), icount(1:3))
    call checkncerr(istatus,__FILE__,__LINE__,'variable ps read error')

    write (stdout,*) 'READ IN fields at DATE:' , tochar(idate)
    !
    ! Calculate Heights on sigma surfaces.
    !
    call htsig_o(t,z1,ps,ht_in,sig,ptop_in,jx_in,iy_in,kz_in)
    !
    ! Interpolate H,U,V,T,Q
    !
    ! 1. For Heights
    call height_o(hp,z1,t,ps,ht_in,sig,ptop_in,jx_in,iy_in,kz_in,plev,np)
    ! 2. For Zonal and Meridional Winds
    call intlin_o(up,u,ps,sig,ptop_in,jx_in,iy_in,kz_in,plev,np)
    call intlin_o(vp,v,ps,sig,ptop_in,jx_in,iy_in,kz_in,plev,np)
    ! 3. For Temperatures
    call intlog_o(tp,t,ps,sig,ptop_in,jx_in,iy_in,kz_in,plev,np)
    ! 4. For Moisture qva
    call humid1_o(t,q,ps,sig,ptop_in,jx_in,iy_in,kz_in)
    call intlin_o(qp,q,ps,sig,ptop_in,jx_in,iy_in,kz_in,plev,np)
    call uvrot4nx(up,vp,xlon_in,xlat_in,clon_in,clat_in, &
                  xcone_in,jx_in,iy_in,np,plon_in,plat_in,iproj_in)
    !
    ! HORIZONTAL INTERPOLATION OF BOTH THE SCALAR AND VECTOR FIELDS
    !
    call cressmcr(b3,b2,xlon,xlat,xlon_in,xlat_in,jx,iy,jx_in,iy_in,np,3)
    call cressmdt(d3,d2,dlon,dlat,xlon_in,xlat_in,jx,iy,jx_in,iy_in,np,2)
    !
    ! ROTATE U-V FIELDS AFTER HORIZONTAL INTERPOLATION
    !
    call uvrot4(u3,v3,dlon,dlat,clon,clat,xcone,jx,iy,np,plon,plat,iproj)
    !
    ! X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
    ! V E R T I C A L   I N T E R P O L A T I O N
    ! X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
    !
    call top2btm(t3,jx,iy,np)
    call top2btm(q3,jx,iy,np)
    call top2btm(h3,jx,iy,np)
    call top2btm(u3,jx,iy,np)
    call top2btm(v3,jx,iy,np)
    !
    ! NEW CALCULATION OF P* ON RegCM TOPOGRAPHY.
    !
    call intgtb(pa,za,tlayer,topogm,t3,h3,sigmar,jx,iy,np)
   
    call intpsn(ps4,topogm,pa,za,tlayer,ptop,jx,iy)
    if(i_band == 1) then
       call p1p2_band(b3pd,ps4,jx,iy)
    else
       call p1p2(b3pd,ps4,jx,iy)
    endif
    !
    ! F0    DETERMINE SURFACE TEMPS ON RegCM TOPOGRAPHY.
    ! INTERPOLATION FROM PRESSURE LEVELS AS IN INTV2
    !
    call intv3(ts4,t3,ps4,sigmar,ptop,jx,iy,np)
    call readsst(ts4,idate)
    !
    ! F3    INTERPOLATE U, V, T, AND Q.
    !
    call intv1(u4,u3,b3pd,sigma2,sigmar,ptop,jx,iy,kz,np)
    call intv1(v4,v3,b3pd,sigma2,sigmar,ptop,jx,iy,kz,np)
    call intv2(t4,t3,ps4,sigma2,sigmar,ptop,jx,iy,kz,np)
    call intv1(q4,q3,ps4,sigma2,sigmar,ptop,jx,iy,kz,np)
    call humid2(t4,q4,ps4,ptop,sigma2,jx,iy,kz)
    !
    ! F4    DETERMINE H
    !
    call hydrost(h4,t4,topogm,ps4,ptop,sigma2,jx,iy,kz)
  end subroutine get_nest
!
  subroutine headernest
    use netcdf
    implicit none
!
    real(rk8) :: xsign
    integer(ik4) :: i , k , istatus , idimid , ivarid
    type(rcm_time_and_date) :: imf
    real(rk8) , dimension(2) :: trlat
    real(rk8) , dimension(:) , allocatable :: sigfix
!
    plev(1) = 50.
    plev(2) = 70.
    plev(3) = 100.
    plev(4) = 150.
    plev(5) = 200.
    plev(6) = 250.
    plev(7) = 300.
    plev(8) = 400.
    plev(9) = 500.
    plev(10) = 600.
    plev(11) = 700.
    plev(12) = 775.
    plev(13) = 850.
    plev(14) = 925.
    plev(15) = 1000.
   
    do k = 1 , np
      sigmar(k) = plev(k)*0.001
    end do
   
    imf = monfirst(globidate1)
    write (fillin,'(a,i10)') 'ATM.', toint10(imf)

    if ( coarsedir(1:5) == '     ' ) then
      inpfile = trim(inpglob)//pthsep//'RegCM'//pthsep
    else
      inpfile = trim(coarsedir)//pthsep
    end if
    if ( coarsedom(1:5) == '     ' ) then
      inpfile = trim(inpfile)//fillin//'.nc'
    else
      inpfile = trim(inpfile)//trim(coarsedom)//'_'//fillin//'.nc'
    end if

    istatus = nf90_open(inpfile, nf90_nowrite, ncinp)
    call checkncerr(istatus,__FILE__,__LINE__, 'Error opening '//trim(inpfile))

    istatus = nf90_inq_dimid(ncinp, 'iy', idimid)
    call checkncerr(istatus,__FILE__,__LINE__,'Dimension iy missing')
    istatus = nf90_inquire_dimension(ncinp, idimid, len=iy_in)
    call checkncerr(istatus,__FILE__,__LINE__,'Dimension iy read error')
    istatus = nf90_inq_dimid(ncinp, 'jx', idimid)
    call checkncerr(istatus,__FILE__,__LINE__,'Dimension jx missing')
    istatus = nf90_inquire_dimension(ncinp, idimid, len=jx_in)
    call checkncerr(istatus,__FILE__,__LINE__,'Dimension jx read error')
    istatus = nf90_inq_dimid(ncinp, 'kz', idimid)
    call checkncerr(istatus,__FILE__,__LINE__,'Dimension kz missing')
    istatus = nf90_inquire_dimension(ncinp, idimid, len=kz_in)
    call checkncerr(istatus,__FILE__,__LINE__,'Dimension kz read error')
    istatus = nf90_inq_dimid(ncinp, 'time', idimid)
    call checkncerr(istatus,__FILE__,__LINE__,'Dimension time missing')
    istatus = nf90_inquire_dimension(ncinp, idimid, len=nrec)
    call checkncerr(istatus,__FILE__,__LINE__,'Dimension time read error')
    istatus = nf90_inq_varid(ncinp, 'time', ivarid)
    call checkncerr(istatus,__FILE__,__LINE__,'variable time missing')
    istatus = nf90_get_att(ncinp, ivarid, 'units', timeunits)
    call checkncerr(istatus,__FILE__,__LINE__,'variable time units missing')
    istatus = nf90_get_att(ncinp, ivarid, 'calendar', timecal)
    call checkncerr(istatus,__FILE__,__LINE__,'variable time calendar missing')
    call getmem1d(itimes,1,nrec,'mod:nest:itimes')
    call getmem1d(xtimes,1,nrec,'mod:nest:xtimes')
    istatus = nf90_get_var(ncinp, ivarid, xtimes)
    call checkncerr(istatus,__FILE__,__LINE__,'variable time read error')
    do i = 1 , nrec
      itimes(i) = timeval2date(xtimes(i), timeunits, timecal)
    end do

    ! Reserve space for I/O

    call getmem1d(sig,1,kz_in,'mod_nest:sig')
    call getmem3d(b2,1,jx_in,1,iy_in,1,np*3,'mod_nest:b2')
    call getmem3d(d2,1,jx_in,1,iy_in,1,np*2,'mod_nest:d2')
    call getmem3d(q,1,jx_in,1,iy_in,1,kz_in,'mod_nest:q')
    call getmem3d(t,1,jx_in,1,iy_in,1,kz_in,'mod_nest:t')
    call getmem3d(u,1,jx_in,1,iy_in,1,kz_in,'mod_nest:u')
    call getmem3d(v,1,jx_in,1,iy_in,1,kz_in,'mod_nest:v')
    call getmem3d(z1,1,jx_in,1,iy_in,1,kz_in,'mod_nest:z1')
    call getmem2d(ps,1,jx_in,1,iy_in,'mod_nest:ps')
    call getmem2d(xlat_in,1,jx_in,1,iy_in,'mod_nest:xlat_in')
    call getmem2d(xlon_in,1,jx_in,1,iy_in,'mod_nest:xlon_in')
    call getmem2d(ht_in,1,jx_in,1,iy_in,'mod_nest:ht_in')
    call getmem3d(b3,1,iy,1,jx,1,np*3,'mod_nest:b3')
    call getmem3d(d3,1,iy,1,jx,1,np*2,'mod_nest:d3')

    istatus = nf90_inq_varid(ncinp, 'sigma', ivarid) 
    call checkncerr(istatus,__FILE__,__LINE__,'variable sigma error')
    istatus = nf90_get_var(ncinp, ivarid, sig)
    call checkncerr(istatus,__FILE__,__LINE__,'variable sigma read error')
    if ( sig(1) < dlowval ) then
      ! Fix V. 4.3.x bug in sigma levels
      allocate(sigfix(kz_in+1))
      sigfix(1:kz_in) = sig(:)
      sigfix(kz_in+1) = d_one
      do k = 1 , kz_in
        sig(k) = d_half*(sigfix(k)+sigfix(k+1))
      end do
      deallocate(sigfix)
    end if
    istatus = nf90_inq_varid(ncinp, 'xlat', ivarid) 
    call checkncerr(istatus,__FILE__,__LINE__,'variable xlat error')
    istatus = nf90_get_var(ncinp, ivarid, xlat_in)
    call checkncerr(istatus,__FILE__,__LINE__,'variable xlat read error')
    istatus = nf90_inq_varid(ncinp, 'xlon', ivarid) 
    call checkncerr(istatus,__FILE__,__LINE__,'variable xlon error')
    istatus = nf90_get_var(ncinp, ivarid, xlon_in)
    call checkncerr(istatus,__FILE__,__LINE__,'variable xlon read error')
    istatus = nf90_inq_varid(ncinp, 'topo', ivarid) 
    call checkncerr(istatus,__FILE__,__LINE__,'variable topo error')
    istatus = nf90_get_var(ncinp, ivarid, ht_in)
    call checkncerr(istatus,__FILE__,__LINE__,'variable topo read error')
    istatus = nf90_inq_varid(ncinp, 'ptop', ivarid) 
    call checkncerr(istatus,__FILE__,__LINE__,'variable ptop error')
    istatus = nf90_get_var(ncinp, ivarid, ptop_in)
    call checkncerr(istatus,__FILE__,__LINE__,'variable ptop read error')

    istatus = nf90_get_att(ncinp, nf90_global,'projection', iproj_in)
    call checkncerr(istatus,__FILE__,__LINE__,'attribure iproj read error')
    istatus = nf90_get_att(ncinp, nf90_global, &
                      'latitude_of_projection_origin', clat_in)
    call checkncerr(istatus,__FILE__,__LINE__,'attribure clat read error')
    istatus = nf90_get_att(ncinp, nf90_global, &
                      'longitude_of_projection_origin', clon_in)
    call checkncerr(istatus,__FILE__,__LINE__,'attribure clat read error')

    if ( iproj_in == 'LAMCON' ) then
      istatus = nf90_get_att(ncinp, nf90_global, 'standard_parallel', trlat)
      call checkncerr(istatus,__FILE__,__LINE__,'attribure truelat read error')
      if ( clat_in < 0. ) then
        xsign = -1.0D0       ! SOUTH HEMESPHERE
      else
        xsign = 1.0D0        ! NORTH HEMESPHERE
      end if
      if ( abs(trlat(1)-trlat(2)) > 1.E-1 ) then
        xcone_in = (log10(cos(trlat(1)*degrad))                    &
                    -log10(cos(trlat(2)*degrad))) /                &
                    (log10(tan((45.0-xsign*trlat(1)/2.0)*degrad))  &
                    -log10(tan((45.0-xsign*trlat(2)/2.0)*degrad)))
      else
        xcone_in = xsign*sin(dble(trlat(1))*degrad)
      end if
    else if ( iproj_in == 'POLSTR' ) then
      xcone_in = 1.0D0
    else if ( iproj_in == 'NORMER' ) then
      xcone_in = 0.0D0
    else
      istatus = nf90_get_att(ncinp, nf90_global, &
                      'grid_north_pole_latitude', plat_in)
      call checkncerr(istatus,__FILE__,__LINE__,'attribure plat read error')
      istatus = nf90_get_att(ncinp, nf90_global, &
                      'grid_north_pole_longitude', plon_in)
      call checkncerr(istatus,__FILE__,__LINE__,'attribure plon read error')
      xcone_in = 0.0D0
    end if
   
    ! Set up pointers
   
    tp => b2(:,:,1:np)
    qp => b2(:,:,np+1:2*np)
    hp => b2(:,:,2*np+1:3*np)
    up => d2(:,:,1:np)
    vp => d2(:,:,np+1:2*np)
    t3 => b3(:,:,1:np)
    q3 => b3(:,:,np+1:2*np)
    h3 => b3(:,:,2*np+1:3*np)
    u3 => d3(:,:,1:np)
    v3 => d3(:,:,np+1:2*np)

  end subroutine headernest
!
end module mod_nest
