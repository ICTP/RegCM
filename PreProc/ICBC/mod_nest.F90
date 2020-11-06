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
  use mod_kdinterp
  use mod_vertint
  use mod_stdatm
  use mod_hgt
  use mod_humid
  use mod_projections
  use mod_vectutil
  use mod_message
  use mod_memutil
  use mod_nchelper
  use mod_domain , only : read_reference_surface_temp

  private

  public :: get_nest , init_nest , conclude_nest

  character(len=256) , public :: coarsedir , coarsedom

  integer(ik4) :: nrec

  type(regcm_projection) :: pj
  type(anyprojparams) :: pjparam

  real(rkx) , pointer , dimension(:,:,:) :: q_in , t_in , p_in , pr0_in
  real(rkx) , pointer , dimension(:,:,:) :: u_in , v_in , z_in
  real(rkx) , pointer , dimension(:,:,:) :: ppa_in , t0_in , qc_in , qi_in
  real(rkx) , pointer , dimension(:,:) :: ht_in , ps_in , ts_in , p0_in
  real(rkx) , pointer , dimension(:,:) :: xlat_in , xlon_in
  real(rkx) , pointer , dimension(:) :: ak_in , bk_in
  real(rkx) , pointer , dimension(:,:) :: pstar_in
  real(rkx) , pointer , dimension(:,:) :: ts , ps , zs
  real(rkx) , pointer , dimension(:,:) :: topou , topov

  real(rkx) , pointer , dimension(:,:,:) :: t3 , q3 , z3 , zud3 , zvd3
  real(rkx) , pointer , dimension(:,:,:) :: u3 , v3 , p3 , pd3 , qc3 , qi3
  real(rkx) , pointer , dimension(:,:,:) :: u3v , v3u

  real(rkx) , pointer , dimension(:,:,:) :: p_out , pd_out

  real(rkx) :: ts0_in , bsp_in
  real(rkx) , pointer , dimension(:) :: sigma_in

  integer(ik4) :: iy_in , jx_in , kz_in

  integer(ik4) :: oidyn
  real(rkx) :: ptop_in , ptop_out

  logical :: uvrotate = .false.
  logical :: has_qc = .false.
  logical :: has_qi = .false.

  character(len=14) :: fillin
  character(len=256) :: inpfile

  integer(ik4) :: ncinp
  type(rcm_time_and_date) , dimension(:) , pointer :: itimes
  real(rkx) , dimension(:) , pointer :: xtimes
  character(len=64) :: timeunits , timecal

  type(h_interpolator) :: cross_hint , udot_hint , vdot_hint

  contains

  subroutine init_nest
    use netcdf
    implicit none

    integer(ik4) :: i , j , k , istatus , idimid , ivarid
    type(rcm_time_and_date) :: imf
    real(rkx) , dimension(2) :: trlat
    character(len=6) :: iproj_in
    real(rkx) :: clat_in , clon_in , plat_in , plon_in , ds_in
    real(rkx) , dimension(:) , allocatable :: sigfix
    real(rkx) :: tlp , alnp
    character(len=16) :: charatt

    imf = monfirst(globidate1)
    write (fillin,'(a,a)') 'ATM.', trim(tochar10(imf))

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
    istatus = nf90_inq_varid(ncinp, 'clw', ivarid)
    if ( istatus == nf90_noerr ) then
      has_qc = .true.
    end if
    istatus = nf90_inq_varid(ncinp, 'cli', ivarid)
    if ( istatus == nf90_noerr ) then
      has_qi = .true.
    end if
    do i = 1 , nrec
      itimes(i) = timeval2date(xtimes(i), timeunits, timecal)
    end do

    istatus = nf90_get_att(ncinp, nf90_global, &
                    'wind_rotated_eastward_northward', charatt)
    if ( istatus == nf90_noerr ) then
      uvrotate = .true.
    end if

    istatus = nf90_get_att(ncinp, nf90_global, 'dynamical_core', oidyn)
    if ( istatus /= nf90_noerr ) then
      oidyn = 1 ! Assume non-hydrostatic
    end if

    ! Reserve space for I/O

    call getmem1d(sigma_in,1,kz_in,'mod_nest:sigma_in')
    call getmem3d(q_in,1,jx_in,1,iy_in,1,kz_in,'mod_nest:q_in')
    if ( has_qc ) then
      call getmem3d(qc_in,1,jx_in,1,iy_in,1,kz_in,'mod_nest:qc_in')
    end if
    if ( has_qi ) then
      call getmem3d(qi_in,1,jx_in,1,iy_in,1,kz_in,'mod_nest:qi_in')
    end if
    call getmem3d(t_in,1,jx_in,1,iy_in,1,kz_in,'mod_nest:t_in')
    call getmem3d(u_in,1,jx_in,1,iy_in,1,kz_in,'mod_nest:u_in')
    call getmem3d(v_in,1,jx_in,1,iy_in,1,kz_in,'mod_nest:v_in')
    call getmem3d(z_in,1,jx_in,1,iy_in,1,kz_in,'mod_nest:z_in')
    call getmem3d(p_in,1,jx_in,1,iy_in,1,kz_in,'mod_nest:p_in')
    if ( oidyn == 2 ) then
      call getmem3d(pr0_in,1,jx_in,1,iy_in,1,kz_in,'mod_nest:pr0_in')
    end if
    call getmem2d(ts_in,1,jx_in,1,iy_in,'mod_nest:ts_in')
    call getmem2d(ps_in,1,jx_in,1,iy_in,'mod_nest:ps_in')
    call getmem2d(xlat_in,1,jx_in,1,iy_in,'mod_nest:xlat_in')
    call getmem2d(xlon_in,1,jx_in,1,iy_in,'mod_nest:xlon_in')
    call getmem2d(ht_in,1,jx_in,1,iy_in,'mod_nest:ht_in')

    istatus = nf90_inq_varid(ncinp, 'kz', ivarid)
    if ( istatus /= nf90_noerr ) then
      istatus = nf90_inq_varid(ncinp, 'sigma', ivarid)
    end if
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'variable sigma or kz error')
    istatus = nf90_get_var(ncinp, ivarid, sigma_in)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'variable sigma or kz read error')
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
    if ( oidyn /= 3 ) then
      istatus = nf90_inq_varid(ncinp, 'ptop', ivarid)
      call checkncerr(istatus,__FILE__,__LINE__, &
                    'variable ptop error')
      istatus = nf90_get_var(ncinp, ivarid, ptop_in)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'variable ptop read error')
      ptop_in = ptop_in * d_100
      ptop_out = ptop * d_1000
    else
      ptop_in = d_zero
      ptop_out = d_zero
    end if

    if ( .not. uvrotate ) then
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
      istatus = nf90_get_att(ncinp, nf90_global, &
                        'grid_size_in_meters', ds_in)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'attribure ds read error')
      plat_in = clat_in
      plon_in = clon_in
      trlat(1) = clat_in
      trlat(2) = clat_in
      if ( iproj_in == 'LAMCON' ) then
        istatus = nf90_get_att(ncinp, nf90_global, 'standard_parallel', trlat)
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'attribure truelat read error')
      else if ( iproj_in == 'ROTMER' .or. iproj_in == 'ROTLLR' ) then
        istatus = nf90_get_att(ncinp, nf90_global, &
                        'grid_north_pole_latitude', plat_in)
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'attribure plat read error')
        istatus = nf90_get_att(ncinp, nf90_global, &
                        'grid_north_pole_longitude', plon_in)
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'attribure plon read error')
      end if

      pjparam%pcode = iproj_in
      pjparam%ds = ds_in
      pjparam%clat = clat_in
      pjparam%clon = clon_in
      pjparam%plat = plat_in
      pjparam%plon = plon_in
      pjparam%trlat1 = trlat(1)
      pjparam%trlat2 = trlat(2)
      pjparam%nlon = jx_in
      pjparam%nlat = iy_in
      pjparam%rotparam = .true.
      pjparam%staggerx = .false.
      pjparam%staggery = .false.
      call pj%initialize(pjparam)
    end if

    call getmem2d(pstar_in,1,jx_in,1,iy_in,'mod_nest:pstar_in')
    if ( oidyn == 2 ) then
      call getmem2d(p0_in,1,jx_in,1,iy_in,'mod_nest:p0_in')
    end if

    if ( oidyn == 2 ) then
      istatus = nf90_get_att(ncinp, nf90_global, 'logp_lapse_rate', tlp)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'attribure logp_lapse_rate read error')
      istatus = nf90_get_att(ncinp, nf90_global, &
           'base_state_surface_temperature', ts0_in)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'attribure logp_lapse_rate read error')
      istatus = nf90_get_att(ncinp, nf90_global, &
           'base_state_pressure', bsp_in)
      if ( istatus /= nf90_noerr ) then
        write(stdout,*) 'Assuming base state pressure to be ', stdp
        bsp_in = stdp
      end if
      call getmem3d(ppa_in,1,jx_in,1,iy_in,1,kz_in,'mod_nest:ppa_in')
      call getmem3d(t0_in,1,jx_in,1,iy_in,1,kz_in,'mod_nest:t0_in')
      istatus = nf90_inq_varid(ncinp, 'p0', ivarid)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'variable p0 error')
      istatus = nf90_get_var(ncinp, ivarid, p0_in)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'variable p0 read error')
      pstar_in = p0_in - ptop_in
      do k = 1 , kz_in
        do i = 1 , iy_in
          do j = 1 , jx_in
            pr0_in(j,i,k) = pstar_in(j,i) * sigma_in(k) + ptop_in
            t0_in(j,i,k) = max(ts0_in +  &
                      tlp * log(pr0_in(j,i,k) / bsp_in),tiso)
            alnp = log(pr0_in(j,i,k)/p0_in(j,i))
            z_in(j,i,k) = max(-(d_half*rovg*tlp*alnp*alnp + &
                              rovg*ts0_in*alnp),d_zero)
          end do
        end do
      end do
    else if ( oidyn == 3 ) then
      call getmem1d(ak_in,1,kz_in,'mod_nest:ak_in')
      call getmem1d(bk_in,1,kz_in,'mod_nest:bk_in')
      istatus = nf90_inq_varid(ncinp, 'a', ivarid)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'variable a error')
      istatus = nf90_get_var(ncinp, ivarid, ak_in)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'variable a read error')
      istatus = nf90_inq_varid(ncinp, 'b', ivarid)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'variable b error')
      istatus = nf90_get_var(ncinp, ivarid, bk_in)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'variable b read error')
    end if

    call h_interpolator_create(cross_hint,xlat_in,xlon_in,xlat,xlon)
    if ( idynamic == 3 ) then
      call h_interpolator_create(udot_hint,xlat_in,xlon_in,ulat,ulon)
      call h_interpolator_create(vdot_hint,xlat_in,xlon_in,vlat,vlon)
    else
      call h_interpolator_create(udot_hint,xlat_in,xlon_in,dlat,dlon)
    end if

    ! Fields on output horizontal grid
    call getmem3d(t3,1,jx,1,iy,1,kz_in,'mod_nest:t3')
    call getmem3d(u3,1,jx,1,iy,1,kz_in,'mod_nest:u3')
    call getmem3d(v3,1,jx,1,iy,1,kz_in,'mod_nest:v3')
    call getmem3d(q3,1,jx,1,iy,1,kz_in,'mod_nest:q3')
    call getmem3d(z3,1,jx,1,iy,1,kz_in,'mod_nest:z3')
    call getmem3d(p3,1,jx,1,iy,1,kz_in,'mod_nest:p3')
    if ( has_qc ) then
      call getmem3d(qc3,1,jx,1,iy,1,kz_in,'mod_nest:qc3')
    end if
    if ( has_qi ) then
      call getmem3d(qi3,1,jx,1,iy,1,kz_in,'mod_nest:qi3')
    end if
    if ( idynamic == 3 ) then
      call getmem2d(topou,1,jx,1,iy,'mod_nest:topou')
      call getmem2d(topov,1,jx,1,iy,'mod_nest:topov')
      call getmem3d(zud3,1,jx,1,iy,1,kz_in,'mod_nest:zud3')
      call getmem3d(zvd3,1,jx,1,iy,1,kz_in,'mod_nest:zvd3')
      call getmem3d(u3v,1,jx,1,iy,1,kz_in,'mod_nest:u3v')
      call getmem3d(v3u,1,jx,1,iy,1,kz_in,'mod_nest:v3u')
    else
      call getmem3d(pd3,1,jx,1,iy,1,kz_in,'mod_nest:pd3')
    end if
    call getmem2d(ts,1,jx,1,iy,'mod_nest:ts')
    call getmem2d(ps,1,jx,1,iy,'mod_nest:ps')
    call getmem2d(zs,1,jx,1,iy,'mod_nest:zs')

    if ( idynamic /= 3 ) then
      call getmem3d(p_out,1,jx,1,iy,1,kz,'mod_nest:p_out')
      call getmem3d(pd_out,1,jx,1,iy,1,kz,'mod_nest:pd_out')
    end if

    if ( oidyn == 3 ) then
      ! Compute vertical model level of input MOLOCH grid
      do k = 1 , kz_in
        do i = 1 , iy_in
          do j = 1 , jx_in
            z_in(j,i,k) = ak_in(k) + bk_in(k) * ht_in(j,i)
          end do
        end do
      end do
      call h_interpolate_cont(cross_hint,z_in,z3)
      call top2btm(z3)
    end if
    if ( idynamic == 2 ) then
      do k = 1 , kz
        do i = 1 , iy
          do j = 1 , jx
            p_out(j,i,k) = ps0(j,i)*sigmah(k) + ptop_out
          end do
        end do
      end do
      call crs2dot(pd_out,p_out,jx,iy,kz,i_band,i_crm)
    end if
    if ( oidyn == 2 ) then
      call h_interpolate_cont(cross_hint,z_in,z3)
      call top2btm(z3)
    end if

    call h_interpolate_cont(cross_hint,ht_in,zs)

    if ( idynamic == 3 ) then
      call top2btm(z0)
      call ucrs2dot(zud4,z0,jx,iy,kz,i_band)
      call vcrs2dot(zvd4,z0,jx,iy,kz,i_crm)
      call ucrs2dot(topou,topogm,jx,iy,i_band)
      call vcrs2dot(topov,topogm,jx,iy,i_crm)
    end if

  end subroutine init_nest

  subroutine get_nest(idate)
    use netcdf
    implicit none

    type(rcm_time_and_date) , intent(in) :: idate

    integer(ik4) :: i , j , k , istatus , ivarid , idimid , irec
    integer(ik4) , dimension(4) :: istart , icount
    type(rcm_time_and_date) :: imf
    logical :: lspch

    if ( idate > itimes(nrec) ) then
      istatus = nf90_close(ncinp)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error close')
      imf = monfirst(idate)
      write (fillin,'(a,a)') 'ATM.', trim(tochar10(imf))
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
    else if ( idate < itimes(1) ) then
      ! Try previous month !
      istatus = nf90_close(ncinp)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error close '//trim(inpfile))
      imf = prevmon(idate)
      write (fillin,'(a,a)') 'ATM.', trim(tochar10(imf))
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
    if ( irec < 0 ) then
      write (stderr,*) 'Error : time ', trim(tochar(idate)), ' not in file'
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
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'variable u/ua missing')
    end if
    istatus = nf90_get_var(ncinp, ivarid, u_in, istart, icount)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'variable u/ua read error')
    istatus = nf90_inq_varid(ncinp, 'va', ivarid)
    if ( istatus /= nf90_noerr ) then
      istatus = nf90_inq_varid(ncinp, 'v', ivarid)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'variable v/va missing')
    end if
    istatus = nf90_get_var(ncinp, ivarid, v_in, istart, icount)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'variable v/va read error')
    istatus = nf90_inq_varid(ncinp, 'ta', ivarid)
    if ( istatus /= nf90_noerr ) then
      istatus = nf90_inq_varid(ncinp, 't', ivarid)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'variable t/ta missing')
    end if
    istatus = nf90_get_var(ncinp, ivarid, t_in, istart, icount)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'variable t/ta read error')
    lspch = .true.
    istatus = nf90_inq_varid(ncinp, 'hus', ivarid)
    if ( istatus /= nf90_noerr ) then
      istatus = nf90_inq_varid(ncinp, 'qas', ivarid)
      if ( istatus /= nf90_noerr ) then
        lspch = .false.
        istatus = nf90_inq_varid(ncinp, 'qv', ivarid)
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'variable qv/qas/hus missing')
      end if
    end if
    istatus = nf90_get_var(ncinp, ivarid, q_in, istart, icount)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'variable qv read error')
    ! Transform specific humidity in mixing ratio
    if ( lspch ) then
      call sph2mxr(q_in,jx_in,iy_in,kz_in)
    end if
    q_in = max(minqq,q_in)
    if ( has_qc ) then
      istatus = nf90_inq_varid(ncinp, 'clw', ivarid)
      call checkncerr(istatus,__FILE__,__LINE__, &
                        'variable clw missing')
      istatus = nf90_get_var(ncinp, ivarid, qc_in, istart, icount)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'variable clw read error')
    end if
    if ( has_qi ) then
      istatus = nf90_inq_varid(ncinp, 'cli', ivarid)
      call checkncerr(istatus,__FILE__,__LINE__, &
                        'variable cli missing')
      istatus = nf90_get_var(ncinp, ivarid, qi_in, istart, icount)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'variable cli read error')
    end if
    if ( oidyn == 2 ) then
      istatus = nf90_inq_varid(ncinp, 'ppa', ivarid)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'variable ppa missing')
      istatus = nf90_get_var(ncinp, ivarid, ppa_in, istart, icount)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'variable ppa read error')
      do k = 1 , kz_in
        do i = 1 , iy_in
          do j = 1 , jx_in
            p_in(j,i,k) = pr0_in(j,i,k) + ppa_in(j,i,k)
          end do
        end do
      end do
    else if ( oidyn == 3 ) then
      istatus = nf90_inq_varid(ncinp, 'pai', ivarid)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'variable pai missing')
      istatus = nf90_get_var(ncinp, ivarid, p_in, istart, icount)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'variable pai read error')
      do k = 1 , kz_in
        do i = 1 , iy_in
          do j = 1 , jx_in
            p_in(j,i,k) = (p_in(j,i,k)**cpovr) * p00
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
    istatus = nf90_get_var(ncinp, ivarid, ps_in, istart(1:3), icount(1:3))
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'variable ps read error')
    istatus = nf90_inq_varid(ncinp, 'ts', ivarid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'variable ts missing')
    istatus = nf90_get_var(ncinp, ivarid, ts_in, istart(1:3), icount(1:3))
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'variable ts read error')
    if ( oidyn == 1 ) then
      pstar_in = ps_in - ptop_in
      do k = 1 , kz_in
        do i = 1 , iy_in
          do j = 1 , jx_in
            p_in(j,i,k) = pstar_in(j,i) * sigma_in(k) + ptop_in
          end do
        end do
      end do
    end if

    write (stdout,*) 'READ IN fields at DATE:' , trim(tochar(idate))
    !
    ! Calculate Heights on sigma surfaces.
    !
    ! Interpolate H,U,V,T,Q
    !
    ! In this module, all pressures are in Pascal.
    !
    if ( .not. uvrotate ) then
      call pj%wind_antirotate(u_in,v_in)
    end if

    if ( oidyn == 1 ) then
      call hydrost(z_in,t_in,ht_in,ps_in,ptop_in,sigma_in,jx_in,iy_in,kz_in)
      call h_interpolate_cont(cross_hint,z_in,z3)
      call top2btm(z3)
    end if
    !
    ! Horizontal interpolation of both the scalar and vector fields
    !
    call h_interpolate_cont(cross_hint,t_in,t3)
    call h_interpolate_cont(cross_hint,q_in,q3)
    call h_interpolate_cont(cross_hint,p_in,p3)
    if ( idynamic == 3 ) then
      call h_interpolate_cont(udot_hint,u_in,u3)
      call h_interpolate_cont(udot_hint,v_in,v3u)
      call h_interpolate_cont(vdot_hint,u_in,u3v)
      call h_interpolate_cont(vdot_hint,v_in,v3)
    else
      call h_interpolate_cont(udot_hint,u_in,u3)
      call h_interpolate_cont(udot_hint,v_in,v3)
    end if
    call h_interpolate_cont(cross_hint,ts_in,ts)
    call h_interpolate_cont(cross_hint,ps_in,ps)
    if ( has_qc ) then
      call h_interpolate_cont(cross_hint,qc_in,qc3)
    end if
    if ( has_qi ) then
      call h_interpolate_cont(cross_hint,qi_in,qi3)
    end if
    !
    ! Rotate U-V fields after horizontal interpolation
    !
    if ( idynamic == 3 ) then
      call pju%wind_rotate(u3,v3u)
      call pjv%wind_rotate(u3v,v3)
    else
      call pjd%wind_rotate(u3,v3)
    end if
    !
    ! Vertical interpolation
    !
!$OMP SECTIONS
!$OMP SECTION
    call top2btm(t3)
!$OMP SECTION
    call top2btm(q3)
!$OMP SECTION
    call top2btm(u3)
!$OMP SECTION
    call top2btm(v3)
!$OMP SECTION
    call top2btm(p3)
!$OMP SECTION
    if ( has_qc ) then
      call top2btm(qc3)
      where ( qc3 < 1.01e-8_rkx ) qc3 = 0.0_rkx
    end if
!$OMP SECTION
    if ( has_qi ) then
      call top2btm(qi3)
      where ( qi3 < 1.01e-8_rkx ) qi3 = 0.0_rkx
    end if
!$OMP END SECTIONS
    !
    ! New calculation of P* on RegCM topography.
    !
    call intpsn(ps4,topogm,ps,zs,ts,ptop_out,jx,iy)

    if ( idynamic /= 3 ) then
      if ( idynamic == 1 ) then
        do k = 1 , kz
          do i = 1 , iy
            do j = 1 , jx
              p_out(j,i,k) = ps4(j,i) * sigmah(k) + ptop_out
            end do
          end do
        end do
        call crs2dot(pd_out,p_out,jx,iy,kz,i_band,i_crm)
      end if
      call crs2dot(pd4,ps4,jx,iy,i_band,i_crm)
      call crs2dot(pd3,p3,jx,iy,kz_in,i_band,i_crm)
      ps = ps4 + ptop
      call intp3(ts4,t3,p3,ps,jx,iy,kz_in,0.0_rkx,0.05_rkx,0.05_rkx)
    else
      call ucrs2dot(zud3,z3,jx,iy,kz_in,i_band)
      call vcrs2dot(zvd3,z3,jx,iy,kz_in,i_crm)
      call intz3(ts4,t3,z3,topogm,jx,iy,kz_in,0.0_rkx,0.05_rkx,0.05_rkx)
    end if

    where ( mask == 0 )
      ts4(:,:) = ts(:,:)
    end where
    !
    ! Interpolate U, V, T, and Q.
    !
    if ( idynamic == 3 ) then
!$OMP SECTIONS
!$OMP SECTION
      call intz1(u4,u3,zud4,zud3,topou,jx,iy,kz,kz_in,0.6_rkx,0.2_rkx,0.2_rkx)
      call top2btm(u4)
!$OMP SECTION
      call intz1(v4,v3,zvd4,zvd3,topov,jx,iy,kz,kz_in,0.6_rkx,0.2_rkx,0.2_rkx)
      call top2btm(v4)
!$OMP SECTION
      call intz1(t4,t3,z0,z3,topogm,jx,iy,kz,kz_in,0.6_rkx,0.85_rkx,0.5_rkx)
      call top2btm(t4)
!$OMP SECTION
      call intz1(q4,q3,z0,z3,topogm,jx,iy,kz,kz_in,0.7_rkx,0.7_rkx,0.4_rkx)
      call top2btm(q4)
!$OMP SECTION
      if ( has_qc ) then
        call intz1(qc4,qc3,z0,z3,topogm,jx,iy,kz,kz_in,0.7_rkx,0.7_rkx,0.4_rkx)
        call top2btm(qc4)
      end if
      if ( has_qi ) then
        call intz1(qi4,qi3,z0,z3,topogm,jx,iy,kz,kz_in,0.7_rkx,0.7_rkx,0.4_rkx)
        call top2btm(qi4)
      end if
!$OMP END SECTIONS
    else
!$OMP SECTIONS
!$OMP SECTION
      call intp1(u4,u3,pd_out,pd3,jx,iy,kz,kz_in,0.6_rkx,0.2_rkx,0.2_rkx)
!$OMP SECTION
      call intp1(v4,v3,pd_out,pd3,jx,iy,kz,kz_in,0.6_rkx,0.2_rkx,0.2_rkx)
!$OMP SECTION
      call intp1(t4,t3,p_out,p3,jx,iy,kz,kz_in,0.6_rkx,0.85_rkx,0.5_rkx)
!$OMP SECTION
      call intp1(q4,q3,p_out,p3,jx,iy,kz,kz_in,0.7_rkx,0.7_rkx,0.4_rkx)
!$OMP SECTION
      if ( has_qc ) then
        call intp1(qc4,qc3,p_out,p3,jx,iy,kz,kz_in,0.7_rkx,0.7_rkx,0.4_rkx)
      end if
      if ( has_qi ) then
        call intp1(qi4,qi3,p_out,p3,jx,iy,kz,kz_in,0.7_rkx,0.7_rkx,0.4_rkx)
      end if
!$OMP END SECTIONS
    end if
    !
    ! Put surface pressures as ps and in cb now to be conforming to
    ! other modules.
    !
    ps4 = ps4 * d_r1000
    if ( idynamic /= 3 ) then
      pd4 = pd4 * d_r1000
    end if
  end subroutine get_nest

  subroutine conclude_nest
    implicit none
    call h_interpolator_destroy(cross_hint)
    call h_interpolator_destroy(udot_hint)
    if ( idynamic == 3 ) then
      call h_interpolator_destroy(vdot_hint)
    end if
    call pj%destruct( )
  end subroutine conclude_nest

end module mod_nest
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
