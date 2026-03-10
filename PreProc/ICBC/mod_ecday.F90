!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    Use of this source code is governed by an MIT-style license that can
!    be found in the LICENSE file or at
!
!         https://opensource.org/licenses/MIT.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

module mod_ecday

  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_stdio
  use mod_grid
  use mod_write
  use mod_date
  use mod_constants
  use mod_kdinterp
  use mod_vertint
  use mod_hgt
  use mod_humid
  use mod_mksst
  use mod_projections
  use mod_vectutil
  use mod_message
  use mod_memutil
  use mod_nchelper
  implicit none

  private

  integer(ik4) :: klev
  integer(ik4) :: jlat
  integer(ik4) :: ilon

  real(rkx), pointer, contiguous, dimension(:) :: glat
  real(rkx), pointer, contiguous, dimension(:) :: glon
  real(rkx), pointer, contiguous, dimension(:,:) :: psvar
  real(rkx), pointer, contiguous, dimension(:) :: sigmar, plevs
  real(rkx) :: pss, pst

  real(rkx), pointer, contiguous, dimension(:,:,:) :: b2
  real(rkx), pointer, contiguous, dimension(:,:,:) :: d2
  real(rkx), pointer, contiguous, dimension(:,:,:) :: b3
  real(rkx), pointer, contiguous, dimension(:,:,:) :: d3
  real(rkx), pointer, contiguous, dimension(:,:,:) :: d3u
  real(rkx), pointer, contiguous, dimension(:,:,:) :: d3v
  real(rkx), pointer, contiguous, dimension(:,:,:) :: h3v, h3u
  real(rkx), pointer, contiguous, dimension(:,:) :: topou, topov
  real(rkx), pointer, contiguous, dimension(:,:,:) :: u3, v3
  real(rkx), pointer, contiguous, dimension(:,:,:) :: u3v, v3u
  real(rkx), pointer, contiguous, dimension(:,:,:) :: h3, q3, t3
  real(rkx), pointer, contiguous, dimension(:,:,:) :: uvar, vvar
  real(rkx), pointer, contiguous, dimension(:,:,:) :: hvar, qvar, tvar

  integer(ik4) :: year, month, day, hour
  integer(ik4) :: itcfs = 0

  public :: get_ecday, init_ecday, conclude_ecday

  type(h_interpolator) :: cross_hint, udot_hint, vdot_hint

  contains

  subroutine init_ecday
    use netcdf
    implicit none

    integer(ik4) :: k, year, month, day, hour
    integer(ik4) :: istatus, inet, iddim, idv
    character(len=256) :: inpfile

    call split_idate(globidate1, year, month, day, hour)
    write (inpfile,'(a,i0.4,a)') trim(inpglob)//pthsep// &
           dattyp//pthsep//'ta'//pthsep//'ta_',year,'.nc'
    istatus = nf90_open(inpfile,nf90_nowrite,inet)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error opening '//trim(inpfile))
    istatus = nf90_inq_dimid(inet,'lon',iddim)
    if ( istatus /= nf90_noerr ) then
      istatus = nf90_inq_dimid(inet,'longitude',iddim)
    end if
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error find dim lon')
    istatus = nf90_inquire_dimension(inet,iddim, len=ilon)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error inquire dim lon')
    istatus = nf90_inq_dimid(inet,'lat',iddim)
    if ( istatus /= nf90_noerr ) then
      istatus = nf90_inq_dimid(inet,'latitude',iddim)
    end if
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error find dim lat')
    istatus = nf90_inquire_dimension(inet,iddim, len=jlat)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error inquire dim lat')
    istatus = nf90_inq_dimid(inet,'plev',iddim)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error find dim plev')
    istatus = nf90_inquire_dimension(inet,iddim, len=klev)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error inquire dim plev')

    call getmem(glon,1,ilon,'mod_ecday:glon')
    call getmem(glat,1,jlat,'mod_ecday:glat')
    call getmem(sigmar,1,klev,'mod_ecday:sigmar')
    call getmem(plevs,1,klev,'mod_ecday:plevs')

    istatus = nf90_inq_varid(inet,'plev',idv)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error find var plev')
    istatus = nf90_get_var(inet,idv,plevs)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read plev')
    plevs = plevs/100.0_rkx ! Put them in mb.
    do k = 1, klev
      sigmar(k) = (plevs(k)-plevs(klev))/(plevs(1)-plevs(klev))
    end do
    pss = (plevs(1)-plevs(klev)) / d_10 ! centibars
    pst = plevs(klev) / d_10 ! centibars
    !
    ! INITIAL GLOBAL GRID-POINT LONGITUDE & LATITUDE
    !
    istatus = nf90_inq_varid(inet,'lon',idv)
    if ( istatus /= nf90_noerr ) then
      istatus = nf90_inq_varid(inet,'longitude',idv)
    end if
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error find var lon')
    istatus = nf90_get_var(inet,idv,glon)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read lon')
    istatus = nf90_inq_varid(inet,'lat',idv)
    if ( istatus /= nf90_noerr ) then
      istatus = nf90_inq_varid(inet,'latitude',idv)
    end if
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error find var lon')
    istatus = nf90_get_var(inet,idv,glat)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read lat')
    istatus = nf90_close(inet)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Cannot close file')

    write(stdout,*) 'Read static data'

    call h_interpolator_create(cross_hint,glat,glon,xlat,xlon)
    if ( idynamic == 3 ) then
      call h_interpolator_create(udot_hint,glat,glon,ulat,ulon)
      call h_interpolator_create(vdot_hint,glat,glon,vlat,vlon)
    else
      call h_interpolator_create(udot_hint,glat,glon,dlat,dlon)
    end if

    call getmem(psvar,1,ilon,1,jlat,'mod_ecday:psvar')
    call getmem(b2,1,ilon,1,jlat,1,klev*3,'mod_ecday:b3')
    call getmem(d2,1,ilon,1,jlat,1,klev*2,'mod_ecday:d3')
    call getmem(b3,1,jx,1,iy,1,klev*3,'mod_ecday:b3')
    if ( idynamic == 3 ) then
      call getmem(d3u,1,jx,1,iy,1,klev*2,'mod_ecday:d3u')
      call getmem(d3v,1,jx,1,iy,1,klev*2,'mod_ecday:d3v')
      call getmem(h3u,1,jx,1,iy,1,klev,'mod_ecday:h3u')
      call getmem(h3v,1,jx,1,iy,1,klev,'mod_ecday:h3v')
      call getmem(topou,1,jx,1,iy,'mod_ecday:topou')
      call getmem(topov,1,jx,1,iy,'mod_ecday:topov')
    else
      call getmem(d3,1,jx,1,iy,1,klev*2,'mod_ecday:d3')
    end if

    ! Set up pointers

    if ( idynamic == 3 ) then
      u3 => d3u(:,:,1:klev)
      v3u => d3u(:,:,klev+1:2*klev)
      u3v => d3v(:,:,1:klev)
      v3 => d3v(:,:,klev+1:2*klev)
    else
      u3 => d3(:,:,1:klev)
      v3 => d3(:,:,klev+1:2*klev)
    end if
    t3 => b3(:,:,1:klev)
    h3 => b3(:,:,klev+1:2*klev)
    q3 => b3(:,:,2*klev+1:3*klev)
    uvar => d2(:,:,1:klev)
    vvar => d2(:,:,klev+1:2*klev)
    tvar => b2(:,:,1:klev)
    hvar => b2(:,:,klev+1:2*klev)
    qvar => b2(:,:,2*klev+1:3*klev)
    if ( idynamic == 3 ) then
      call ucrs2dot(zud4,z0,jx,iy,kz,i_band)
      call vcrs2dot(zvd4,z0,jx,iy,kz,i_crm)
      call ucrs2dot(topou,topogm,jx,iy,i_band)
      call vcrs2dot(topov,topogm,jx,iy,i_crm)
    end if
  end subroutine init_ecday

  subroutine get_ecday(idate)
    implicit none
    type(rcm_time_and_date), intent(in) :: idate

    call split_idate(idate,year,month,day,hour)
    call eceday(idate)

    write (stdout,*) 'READ IN fields at DATE:', tochar(idate)
    !
    ! Horizontal interpolation of both the scalar and vector fields
    !
    call h_interpolate_cont(cross_hint,b2,b3)
    if ( idynamic == 3 ) then
      call h_interpolate_cont(udot_hint,d2,d3u)
      call h_interpolate_cont(vdot_hint,d2,d3v)
    else
      call h_interpolate_cont(udot_hint,d2,d3)
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
    ! New calculation of P* on rcm topography.
    !
    if ( idynamic == 3 ) then
      call ucrs2dot(h3u,h3,jx,iy,klev,i_band)
      call vcrs2dot(h3v,h3,jx,iy,klev,i_crm)
      call intzps(ps4,topogm,t3,h3,pss,sigmar,pst, &
                  xlat,yeardayfrac(idate),dayspy,jx,iy,klev)
      call intz3(ts4,t3,h3,topogm,0.6_rkx,0.5_rkx,0.85_rkx,jx,iy,klev)
    else
      call intgtb(pa,za,tlayer,topogm,t3,h3,pss,sigmar,pst,jx,iy,klev)
      call intpsn(ps4,topogm,pa,za,tlayer,ptop,jx,iy)
      call crs2dot(pd4,ps4,jx,iy,i_band,i_crm)
      call intv3(ts4,t3,ps4,pss,sigmar,ptop,pst,jx,iy,klev)
    end if

    call readsst(ts4,idate)
    !
    ! Interpolate U, V, T, and Q.
    !
    if ( idynamic == 3 ) then
!$OMP SECTIONS
!$OMP SECTION
      call intz1(u4,u3,zud4,h3u,topou,jx,iy,kz,klev,0.6_rkx,0.2_rkx,0.2_rkx)
!$OMP SECTION
      call intz1(v4,v3,zvd4,h3v,topov,jx,iy,kz,klev,0.6_rkx,0.2_rkx,0.2_rkx)
!$OMP SECTION
      call intz1(t4,t3,z0,h3,topogm,jx,iy,kz,klev,0.6_rkx,0.5_rkx,0.85_rkx)
!$OMP SECTION
      call intz1(q4,q3,z0,h3,topogm,jx,iy,kz,klev,0.7_rkx,0.4_rkx,0.7_rkx)
!$OMP END SECTIONS
    else
!$OMP SECTIONS
!$OMP SECTION
      call intv1(u4,u3,pd4,sigmah,pss,sigmar,ptop,pst,jx,iy,kz,klev,1)
!$OMP SECTION
      call intv1(v4,v3,pd4,sigmah,pss,sigmar,ptop,pst,jx,iy,kz,klev,1)
!$OMP SECTION
      call intv2(t4,t3,ps4,sigmah,pss,sigmar,ptop,pst,jx,iy,kz,klev)
!$OMP SECTION
      call intv1(q4,q3,ps4,sigmah,pss,sigmar,ptop,pst,jx,iy,kz,klev,1)
!$OMP END SECTIONS
    end if
  end subroutine get_ecday

  subroutine eceday(idate)
    use netcdf
    implicit none
    type(rcm_time_and_date), intent (in) :: idate
    integer(ik4) :: i, ilev, inet, it, j, kkrec, k, istatus
    character(len=256), save :: pathaddname
    character(len=5), dimension(6) :: varname
    integer(ik4), dimension(4) :: icount, istart
    integer(ik4), dimension(6), save :: inet5, ivar5
    data varname/'ta', 'zg', 'hus', 'ua', 'va', 'psl'/
    do kkrec = 1, 6
      if ( idate == globidate1 .or. lfdoyear(idate) ) then
        write (pathaddname,'(a,i0.4,a)') trim(inpglob)//pthsep// &
           dattyp//pthsep//trim(varname(kkrec))//pthsep// &
           trim(varname(kkrec))//'_',year,'.nc'
        istatus = nf90_open(pathaddname,nf90_nowrite,inet5(kkrec))
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error opening '//trim(pathaddname))
        istatus = nf90_inq_varid(inet5(kkrec),varname(kkrec),ivar5(kkrec))
        call checkncerr(istatus,__FILE__,__LINE__, &
             'Variable '//varname(kkrec)//' error in file'//trim(pathaddname))
        write (stdout,*) inet5(kkrec), trim(pathaddname)
      end if
      it = dayofyear(idate)
      istart(1) = 1
      icount(1) = ilon
      istart(2) = 1
      icount(2) = jlat
      istart(3) = 1
      icount(3) = klev
      istart(4) = it
      icount(4) = 1
      inet = inet5(kkrec)
      if ( kkrec == 1 ) then
        istatus = nf90_get_var(inet,ivar5(kkrec),tvar,istart,icount)
      else if ( kkrec == 2 ) then
        istatus = nf90_get_var(inet,ivar5(kkrec),hvar,istart,icount)
      else if ( kkrec == 3 ) then
        istatus = nf90_get_var(inet,ivar5(kkrec),qvar,istart,icount)
      else if ( kkrec == 3 ) then
        istatus = nf90_get_var(inet,ivar5(kkrec),uvar,istart,icount)
      else if ( kkrec == 5 ) then
        istatus = nf90_get_var(inet,ivar5(kkrec),vvar,istart,icount)
      else if ( kkrec == 6 ) then
        istart(3) = it
        icount(3) = 1
        istatus = nf90_get_var(inet,ivar5(kkrec),psvar,istart(1:3),icount(1:3))
      end if
      call checkncerr(istatus,__FILE__,__LINE__, &
                  'Variable '//varname(kkrec)// &
                  'read error in file '//trim(pathaddname))
    end do
  end subroutine eceday

  subroutine conclude_ecday
    implicit none
    call h_interpolator_destroy(cross_hint)
    call h_interpolator_destroy(udot_hint)
    if ( idynamic == 3 ) then
      call h_interpolator_destroy(vdot_hint)
    end if
  end subroutine conclude_ecday

end module mod_ecday
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
