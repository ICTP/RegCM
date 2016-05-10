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
  use mod_nhinterp
  use mod_hgt
  use mod_humid
  use mod_uvrot
  use mod_vectutil
  use mod_message
  use mod_memutil
  use mod_nchelper

  private

  public :: get_nest , headernest

  character(len=256) , public :: coarsedir , coarsedom

  integer(ik4) :: nrec

  real(rkx) , pointer , dimension(:,:,:) :: b3
  real(rkx) , pointer , dimension(:,:,:) :: d3
  real(rkx) , pointer , dimension(:,:,:) :: z1

  real(rkx) , pointer , dimension(:,:,:) :: b2
  real(rkx) , pointer , dimension(:,:,:) :: d2

  real(rkx) , pointer , dimension(:,:,:) :: q , t
  real(rkx) , pointer , dimension(:,:,:) :: u , v
  real(rkx) , pointer , dimension(:,:,:) :: pp3d , p3d , t0_in
  real(rkx) , pointer , dimension(:,:) :: ts
  real(rkx) , pointer , dimension(:,:) :: ps , xts , p0_in , pstar0
  real(rkx) , pointer , dimension(:,:) :: ht_in
  real(rkx) , pointer , dimension(:,:) :: xlat_in , xlon_in

  real(rkx) , pointer , dimension(:,:,:) :: h3 , q3 , t3
  real(rkx) , pointer , dimension(:,:,:) :: u3 , v3

  real(rkx) , pointer , dimension(:,:,:) :: hp , qp , tp
  real(rkx) , pointer , dimension(:,:,:) :: up , vp

  real(rkx) , pointer , dimension(:) :: plev , sigmar
  real(rkx) :: pss
  real(rkx) , pointer , dimension(:) :: sigma_in

  integer(ik4) :: iy_in , jx_in , kz_in
  integer(ik4) :: np

  integer(ik4) :: oidyn
  character(len=6) :: iproj_in
  real(rkx) :: clat_in , clon_in , plat_in , plon_in , ptop_in , xcone_in

  character(len=14) :: fillin
  character(len=256) :: inpfile

  integer(ik4) :: ncinp
  type(rcm_time_and_date) , dimension(:) , pointer :: itimes
  real(rkx) , dimension(:) , pointer :: xtimes
  character(len=64) :: timeunits , timecal

  contains

  subroutine get_nest(idate)
    use netcdf
    implicit none

    type(rcm_time_and_date) , intent(in) :: idate

    integer(ik4) :: i , j , k , istatus , ivarid , idimid , irec
    integer(ik4) , dimension(4) :: istart , icount
    type(rcm_time_and_date) :: imf
    logical :: lspch
    real(rkx) :: ptoppa

    ptoppa = ptop * d_1000

    if (.not. associated(b2)) then
      call die('get_nest','Called get_nest before headernest !',1)
    end if

    if ( idate > itimes(nrec) ) then
      istatus = nf90_close(ncinp)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error close')
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
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error opening '//trim(inpfile))
      istatus = nf90_inq_dimid(ncinp, 'time', idimid)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Dimension time missing')
      istatus = nf90_inquire_dimension(ncinp, idimid, len=nrec)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Dimension time read error')
      istatus = nf90_inq_varid(ncinp, 'time', ivarid)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'variable time missing')
      call getmem1d(itimes,1,nrec,'mod_nest:itimes')
      call getmem1d(xtimes,1,nrec,'mod_nest:xtimes')
      istatus = nf90_get_var(ncinp, ivarid, xtimes)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'variable time read error')
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
      ! Try previous month !
      istatus = nf90_close(ncinp)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error close '//trim(inpfile))
      imf = prevmon(idate)
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
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error opening '//trim(inpfile))
      istatus = nf90_inq_dimid(ncinp, 'time', idimid)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Dimension time missing')
      istatus = nf90_inquire_dimension(ncinp, idimid, len=nrec)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Dimension time read error')
      istatus = nf90_inq_varid(ncinp, 'time', ivarid)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'variable time missing')
      call getmem1d(itimes,1,nrec,'mod_nest:itimes')
      call getmem1d(xtimes,1,nrec,'mod_nest:xtimes')
      istatus = nf90_get_var(ncinp, ivarid, xtimes)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'variable time read error')
      do i = 1 , nrec
        itimes(i) = timeval2date(xtimes(i), timeunits,timecal)
      end do
      irec = -1
      do i = 1 , nrec
        if (idate == itimes(i)) then
          irec = i
          exit
        end if
      end do
      if ( irec < 0 ) then
        write (stderr,*) 'Error : time ', tochar(idate), ' not in file'
        call die('get_nest')
      end if
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
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'variable u/ua missing')
    end if
    istatus = nf90_get_var(ncinp, ivarid, u, istart, icount)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'variable u/ua read error')
    istatus = nf90_inq_varid(ncinp, 'va', ivarid)
    if ( istatus /= nf90_noerr ) then
      istatus = nf90_inq_varid(ncinp, 'v', ivarid)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'variable v/va missing')
    end if
    istatus = nf90_get_var(ncinp, ivarid, v, istart, icount)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'variable v/va read error')
    istatus = nf90_inq_varid(ncinp, 'ta', ivarid)
    if ( istatus /= nf90_noerr ) then
      istatus = nf90_inq_varid(ncinp, 't', ivarid)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'variable t/ta missing')
    end if
    istatus = nf90_get_var(ncinp, ivarid, t, istart, icount)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'variable t/ta read error')
    lspch = .true.
    istatus = nf90_inq_varid(ncinp, 'qas', ivarid)
    if ( istatus /= nf90_noerr ) then
      lspch = .false.
      istatus = nf90_inq_varid(ncinp, 'qv', ivarid)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'variable qv/qas missing')
    end if
    istatus = nf90_get_var(ncinp, ivarid, q, istart, icount)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'variable qv read error')
    ! Transform specific humidity in mixing ratio
    if ( lspch ) then
      call sph2mxr(q,jx_in,iy_in,kz_in)
    end if
    if ( oidyn == 2 ) then
      istatus = nf90_inq_varid(ncinp, 'ppa', ivarid)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'variable ppa missing')
      istatus = nf90_get_var(ncinp, ivarid, pp3d, istart, icount)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'variable ppa read error')
      do k = 1 , kz_in
        do i = 1 , iy_in
          do j = 1 , jx_in
            p3d(j,i,k) = pstar0(j,i) * sigma_in(k) + &
                         ptop_in*d_100 + pp3d(j,i,k)
          end do
        end do
      end do
    end if
    istart(3) = irec
    istart(2) = 1
    istart(1) = 1
    icount(3) = 1
    icount(2) = iy_in
    icount(1) = jx_in
    istatus = nf90_inq_varid(ncinp, 'ps', ivarid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'variable ps missing')
    istatus = nf90_get_var(ncinp, ivarid, ps, istart(1:3), icount(1:3))
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'variable ps read error')
    istatus = nf90_inq_varid(ncinp, 'ts', ivarid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'variable ts missing')
    istatus = nf90_get_var(ncinp, ivarid, xts, istart(1:3), icount(1:3))
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'variable ts read error')

    write (stdout,*) 'READ IN fields at DATE:' , tochar(idate)
    !
    ! Calculate Heights on sigma surfaces.
    !
    ! Interpolate H,U,V,T,Q
    !
    ! In this module, all pressures are in Pascal.
    !
    if ( oidyn == 1 ) then
      call htsig_o(t,z1,ps,ht_in,sigma_in,ptop_in,jx_in,iy_in,kz_in)
      call mxr2rh(t,q,ps,sigma_in,ptop_in,jx_in,iy_in,kz_in)
      call intlog(tp,t,ps,sigma_in,ptop_in,jx_in,iy_in,kz_in,plev,np)
      call height_o(hp,z1,t,ps,ht_in,sigma_in,ptop_in,jx_in,iy_in,kz_in,plev,np)
      call intlin(up,u,ps,sigma_in,ptop_in,jx_in,iy_in,kz_in,plev,np)
      call intlin(vp,v,ps,sigma_in,ptop_in,jx_in,iy_in,kz_in,plev,np)
      call intlin(qp,q,ps,sigma_in,ptop_in,jx_in,iy_in,kz_in,plev,np)
    else
      call nonhydrost(z1,t0_in,p0_in,ptop_in,ht_in,sigma_in,jx_in,iy_in,kz_in)
      call mxr2rh(t,q,p3d,jx_in,iy_in,kz_in)
      call intlog(tp,t,ps,p3d,jx_in,iy_in,kz_in,plev,np)
      call height_o(hp,z1,t,ps,ht_in,p3d,jx_in,iy_in,kz_in,plev,np)
      call intlin(up,u,ps,p3d,jx_in,iy_in,kz_in,plev,np)
      call intlin(vp,v,ps,p3d,jx_in,iy_in,kz_in,plev,np)
      call intlin(qp,q,ps,p3d,jx_in,iy_in,kz_in,plev,np)
    end if

    call uvrot4nx(up,vp,xlon_in,xlat_in,clon_in,clat_in, &
                  xcone_in,jx_in,iy_in,np,plon_in,plat_in,iproj_in)
    !
    ! Horizontal interpolation of both the scalar and vector fields
    !
    call cressmcr(b3,b2,xlon,xlat,xlon_in,xlat_in,jx,iy,jx_in,iy_in,np,3)
    call cressmdt(d3,d2,dlon,dlat,xlon_in,xlat_in,jx,iy,jx_in,iy_in,np,2)
    call cressmcr(ts,xts,xlon,xlat,xlon_in,xlat_in,jx,iy,jx_in,iy_in)
    !
    ! Rotate U-V fields after horizontal interpolation
    !
    call uvrot4(u3,v3,dlon,dlat,clon,clat,xcone,jx,iy,np,plon,plat,iproj)
    !
    ! Vertical interpolation
    !
    call top2btm(t3,jx,iy,np)
    call top2btm(q3,jx,iy,np)
    call top2btm(h3,jx,iy,np)
    call top2btm(u3,jx,iy,np)
    call top2btm(v3,jx,iy,np)
    !
    ! New calculation of P* on RegCM topography.
    !
    call intgtb(pa,za,tlayer,topogm,t3,h3,pss,sigmar,jx,iy,np)
    call intpsn(ps4,topogm,pa,za,tlayer,ptoppa,jx,iy)
    call crs2dot(pd4,ps4,jx,iy,i_band)
    !
    ! Determine surface temps on RegCM topography.
    ! Interpolation from pressure levels
    !
    call intv3(ts4,t3,ps4,pss,sigmar,ptoppa,jx,iy,np)
    !
    ! Overwrite SST using TS from ATM file.
    !
    where ( mask == 0 )
      ts4(:,:) = ts(:,:)
    end where
    !
    ! Interpolate U, V, T, and Q.
    !
    call intv1(u4,u3,pd4,sigmah,pss,sigmar,ptoppa,jx,iy,kz,np)
    call intv1(v4,v3,pd4,sigmah,pss,sigmar,ptoppa,jx,iy,kz,np)
    call intv2(t4,t3,ps4,sigmah,pss,sigmar,ptoppa,jx,iy,kz,np)
    call intv1(q4,q3,ps4,sigmah,pss,sigmar,ptoppa,jx,iy,kz,np)
    !
    ! Put surface pressures in cb now to be conforming to other modules.
    !
    ps4 = ps4 * d_r1000
    pd4 = pd4 * d_r1000
    call rh2mxr(t4,q4,ps4,ptop,sigmah,jx,iy,kz)
  end subroutine get_nest

  subroutine headernest
    use netcdf
    implicit none

    real(rkx) :: xsign
    integer(ik4) :: i , j , k , ip , istatus , idimid , ivarid
    type(rcm_time_and_date) :: imf
    real(rkx) , dimension(2) :: trlat
    real(rkx) , dimension(:) , allocatable :: sigfix
    real(rkx) :: tlp , pr0_in

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
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error opening '//trim(inpfile))

    istatus = nf90_inq_dimid(ncinp, 'iy', idimid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Dimension iy missing')
    istatus = nf90_inquire_dimension(ncinp, idimid, len=iy_in)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Dimension iy read error')
    istatus = nf90_inq_dimid(ncinp, 'jx', idimid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Dimension jx missing')
    istatus = nf90_inquire_dimension(ncinp, idimid, len=jx_in)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Dimension jx read error')
    istatus = nf90_inq_dimid(ncinp, 'kz', idimid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Dimension kz missing')
    istatus = nf90_inquire_dimension(ncinp, idimid, len=kz_in)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Dimension kz read error')
    istatus = nf90_inq_dimid(ncinp, 'time', idimid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Dimension time missing')
    istatus = nf90_inquire_dimension(ncinp, idimid, len=nrec)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Dimension time read error')
    istatus = nf90_inq_varid(ncinp, 'time', ivarid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'variable time missing')
    istatus = nf90_get_att(ncinp, ivarid, 'units', timeunits)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'variable time units missing')
    istatus = nf90_get_att(ncinp, ivarid, 'calendar', timecal)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'variable time calendar missing')
    call getmem1d(itimes,1,nrec,'mod:nest:itimes')
    call getmem1d(xtimes,1,nrec,'mod:nest:xtimes')
    istatus = nf90_get_var(ncinp, ivarid, xtimes)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'variable time read error')
    do i = 1 , nrec
      itimes(i) = timeval2date(xtimes(i), timeunits, timecal)
    end do

    ! Reserve space for I/O

    call getmem1d(sigma_in,1,kz_in,'mod_nest:sigma_in')
    call getmem3d(q,1,jx_in,1,iy_in,1,kz_in,'mod_nest:q')
    call getmem3d(t,1,jx_in,1,iy_in,1,kz_in,'mod_nest:t')
    call getmem3d(u,1,jx_in,1,iy_in,1,kz_in,'mod_nest:u')
    call getmem3d(v,1,jx_in,1,iy_in,1,kz_in,'mod_nest:v')
    call getmem3d(z1,1,jx_in,1,iy_in,1,kz_in,'mod_nest:z1')
    call getmem2d(ps,1,jx_in,1,iy_in,'mod_nest:ps')
    call getmem2d(xts,1,jx_in,1,iy_in,'mod_nest:xts')
    call getmem2d(xlat_in,1,jx_in,1,iy_in,'mod_nest:xlat_in')
    call getmem2d(xlon_in,1,jx_in,1,iy_in,'mod_nest:xlon_in')
    call getmem2d(ht_in,1,jx_in,1,iy_in,'mod_nest:ht_in')

    istatus = nf90_inq_varid(ncinp, 'sigma', ivarid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'variable sigma error')
    istatus = nf90_get_var(ncinp, ivarid, sigma_in)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'variable sigma read error')
    if ( sigma_in(1) < dlowval ) then
      ! Fix V. 4.3.x bug in sigma levels
      allocate(sigfix(kz_in+1))
      sigfix(1:kz_in) = sigma_in(:)
      sigfix(kz_in+1) = d_one
      do k = 1 , kz_in
        sigma_in(k) = d_half*(sigfix(k)+sigfix(k+1))
      end do
      deallocate(sigfix)
    end if
    istatus = nf90_inq_varid(ncinp, 'xlat', ivarid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'variable xlat error')
    istatus = nf90_get_var(ncinp, ivarid, xlat_in)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'variable xlat read error')
    istatus = nf90_inq_varid(ncinp, 'xlon', ivarid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'variable xlon error')
    istatus = nf90_get_var(ncinp, ivarid, xlon_in)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'variable xlon read error')
    istatus = nf90_inq_varid(ncinp, 'topo', ivarid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'variable topo error')
    istatus = nf90_get_var(ncinp, ivarid, ht_in)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'variable topo read error')
    istatus = nf90_inq_varid(ncinp, 'ptop', ivarid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'variable ptop error')
    istatus = nf90_get_var(ncinp, ivarid, ptop_in)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'variable ptop read error')

    istatus = nf90_get_att(ncinp, nf90_global,'projection', iproj_in)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'attribure iproj read error')
    istatus = nf90_get_att(ncinp, nf90_global, &
                      'latitude_of_projection_origin', clat_in)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'attribure clat read error')
    istatus = nf90_get_att(ncinp, nf90_global, &
                      'longitude_of_projection_origin', clon_in)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'attribure clat read error')
    istatus = nf90_get_att(ncinp, nf90_global, 'dynamical_core', oidyn)
    if ( istatus /= nf90_noerr ) then
      oidyn = 1 ! Assume non-hydrostatic
    end if

    if ( iproj_in == 'LAMCON' ) then
      istatus = nf90_get_att(ncinp, nf90_global, 'standard_parallel', trlat)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'attribure truelat read error')
      if ( clat_in < 0. ) then
        xsign = -1.0_rkx       ! SOUTH HEMESPHERE
      else
        xsign = 1.0_rkx        ! NORTH HEMESPHERE
      end if
      if ( abs(trlat(1)-trlat(2)) > 1.E-1 ) then
        xcone_in = (log10(cos(trlat(1)*degrad))                    &
                    -log10(cos(trlat(2)*degrad))) /                &
                    (log10(tan((45.0-xsign*trlat(1)/2.0)*degrad))  &
                    -log10(tan((45.0-xsign*trlat(2)/2.0)*degrad)))
      else
        xcone_in = xsign*sin(real(trlat(1),rkx)*degrad)
      end if
    else if ( iproj_in == 'POLSTR' ) then
      xcone_in = 1.0_rkx
    else if ( iproj_in == 'NORMER' ) then
      xcone_in = 0.0_rkx
    else
      istatus = nf90_get_att(ncinp, nf90_global, &
                      'grid_north_pole_latitude', plat_in)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'attribure plat read error')
      istatus = nf90_get_att(ncinp, nf90_global, &
                      'grid_north_pole_longitude', plon_in)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'attribure plon read error')
      xcone_in = 0.0_rkx
    end if

    if ( ptop_in > ptop*10.0_rkx ) then
      write(stderr,*) 'WARNING : top pressure higher than PTOP detected.'
      write(stderr,*) 'WARNING : Extrapolation will be performed.'
    end if

    np = kz_in - 1
    call getmem1d(plev,1,np,'mod_nest:plev')
    call getmem1d(sigmar,1,np,'mod_nest:sigmar')

    call getmem2d(p0_in,1,jx_in,1,iy_in,'mod_nest:p0_in')
    call getmem2d(pstar0,1,jx_in,1,iy_in,'mod_nest:pstar0')

    if ( oidyn == 2 ) then
      istatus = nf90_get_att(ncinp, nf90_global, 'logp_lapse_rate', tlp)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'attribure logp_lapse_rate read error')
      call getmem3d(pp3d,1,jx_in,1,iy_in,1,kz_in,'mod_nest:pp3d')
      call getmem3d(p3d,1,jx_in,1,iy_in,1,kz_in,'mod_nest:p3d')
      call getmem3d(t0_in,1,jx_in,1,iy_in,1,kz_in,'mod_nest:t0_in')
      istatus = nf90_inq_varid(ncinp, 'p0', ivarid)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'variable p0 error')
      istatus = nf90_get_var(ncinp, ivarid, p0_in)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'variable p0 read error')
      pstar0 = p0_in - ptop_in * d_100
      do k = 1 , kz_in
        do i = 1 , iy_in
          do j = 1 , jx_in
            pr0_in = pstar0(j,i) * sigma_in(k) + ptop_in * d_100
            t0_in(j,i,k) = temppres(pr0_in)
          end do
        end do
      end do
    else
      istatus = nf90_inq_varid(ncinp, 'ps', ivarid)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'variable ps error')
      istatus = nf90_get_var(ncinp, ivarid, p0_in)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'variable ps read error')
      pstar0 = p0_in - ptop_in * d_100
    end if

    do ip = 1 , np/2
      plev(ip) = d_half * (minval(pstar0*sigma_in(ip+1)) + &
                           maxval(pstar0*sigma_in(ip))) + ptop_in * d_100
    end do
    do ip = np/2+1 , np
      plev(ip) = maxval(pstar0*sigma_in(ip+1)) + ptop_in * d_100
    end do

    ! Set up pointers

    call getmem3d(b2,1,jx_in,1,iy_in,1,np*3,'mod_nest:b2')
    call getmem3d(d2,1,jx_in,1,iy_in,1,np*2,'mod_nest:d2')
    call getmem3d(b3,1,jx,1,iy,1,np*3,'mod_nest:b3')
    call getmem3d(d3,1,jx,1,iy,1,np*2,'mod_nest:d3')
    call getmem2d(ts,1,jx,1,iy,'mod_nest:ts')
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

    do k = 1 , np
      sigmar(k) = plev(k)/plev(np)
    end do
    pss = plev(np)
  end subroutine headernest

end module mod_nest
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
