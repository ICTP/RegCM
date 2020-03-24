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

module mod_ifs

  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_stdio
  use mod_memutil
  use mod_grid
  use mod_write
  use mod_vertint
  use mod_earth
  use mod_hgt
  use mod_humid
  use mod_projections
  use mod_vectutil
  use mod_message
  use mod_nchelper
  use mod_kdinterp
  use netcdf

  private

  integer(ik4) :: jlat , ilon , klev , timlen

  real(rkx) , pointer , dimension(:,:,:) :: b3
  real(rkx) , pointer , dimension(:,:,:) :: d3
  real(rkx) , pointer , dimension(:,:,:) :: d3u , d3v

  real(rkx) , pointer :: u3(:,:,:) , v3(:,:,:)
  real(rkx) , pointer :: u3v(:,:,:) , v3u(:,:,:)
  real(rkx) , pointer :: h3(:,:,:) , q3(:,:,:) , t3(:,:,:)
  real(rkx) , pointer :: h3u(:,:,:) , h3v(:,:,:)
  real(rkx) , pointer :: uvar(:,:,:) , vvar(:,:,:)
  real(rkx) , pointer :: hvar(:,:,:) , qvar(:,:,:) , tvar(:,:,:)
  real(rkx) , pointer :: topou(:,:) , topov(:,:)
  real(rkx) , pointer :: xps(:,:) , xts(:,:)
  real(rkx) , pointer :: ps(:,:) , ts(:,:)

  real(rkx) , pointer , dimension(:,:,:) :: b2
  real(rkx) , pointer , dimension(:,:,:) :: d2
  real(rkx) , pointer , dimension(:) :: glat
  real(rkx) , pointer , dimension(:) :: ghelp
  real(rkx) , pointer , dimension(:) :: glon
  real(rkx) , pointer , dimension(:) :: slev
  real(rkx) , pointer , dimension(:) :: sigma1 , sigmar
  real(rkx) :: pss

  integer(ik4) :: ncin
  integer(ik4) , dimension(8) :: ivar5

  type(global_domain) :: gdomain
  type(h_interpolator) :: cross_hint , udot_hint , vdot_hint

  public :: init_ifs , get_ifs , conclude_ifs

  contains

  subroutine init_ifs
    implicit none
    integer(ik4) :: k , kr
    integer(ik4) :: year , month , monthp1 , day , hour
    character(len=256) :: pathaddname
    integer(ik4) :: istatus , ncid , ivarid , idimid
    character(len=64) :: inname

    call split_idate(globidate1,year,month,day,hour)
    write(inname,'(a,i0.4,i0.2,i0.2,i0.2,a)') &
      'IFS_',year,month,day,hour,'+000.nc'
    pathaddname = trim(inpglob)//pthsep//'IFS'//pthsep//inname
    istatus = nf90_open(pathaddname,nf90_nowrite,ncid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error open file '//trim(pathaddname))
    istatus = nf90_inq_dimid(ncid,'lat',idimid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Missing lat dimension in file '//trim(pathaddname))
    istatus = nf90_inquire_dimension(ncid,idimid,len=jlat)
    call checkncerr(istatus,__FILE__,__LINE__, &
          'Error reading lat dimelen in file '//trim(pathaddname))
    istatus = nf90_inq_dimid(ncid,'lon',idimid)
    call checkncerr(istatus,__FILE__,__LINE__, &
          'Missing long dimension in file '//trim(pathaddname))
    istatus = nf90_inquire_dimension(ncid,idimid,len=ilon)
    call checkncerr(istatus,__FILE__,__LINE__, &
          'Error reading lon dimelen in file '//trim(pathaddname))
    istatus = nf90_inq_dimid(ncid,'lev',idimid)
    call checkncerr(istatus,__FILE__,__LINE__, &
          'Missing lev dimension in file '//trim(pathaddname))
    istatus = nf90_inquire_dimension(ncid,idimid,len=klev)
    call checkncerr(istatus,__FILE__,__LINE__, &
          'Error reading levelist dimelen in file '//trim(pathaddname))
    !
    ! Allocate working space
    !
    call getmem1d(slev,1,klev,'mod_ifs:slev')
    call getmem1d(glat,1,jlat,'mod_ifs:glat')
    call getmem1d(glon,1,ilon,'mod_ifs:glon')
    call getmem1d(ghelp,1,max(jlat,ilon),'mod_ifs:ghelp')
    call getmem1d(sigma1,1,klev,'mod_ifs:sigma1')
    call getmem1d(sigmar,1,klev,'mod_ifs:sigmar')
    call getmem3d(b3,1,jx,1,iy,1,klev*3,'mod_ifs:b3')
    if ( idynamic == 3 ) then
      call getmem3d(d3u,1,jx,1,iy,1,klev*2,'mod_ifs:d3u')
      call getmem3d(d3v,1,jx,1,iy,1,klev*2,'mod_ifs:d3v')
      call getmem3d(h3u,1,jx,1,iy,1,klev,'mod_ifs:h3u')
      call getmem3d(h3v,1,jx,1,iy,1,klev,'mod_ifs:h3v')
      call getmem2d(topou,1,jx,1,iy,'mod_ifs:topou')
      call getmem2d(topov,1,jx,1,iy,'mod_ifs:topov')
    else
      call getmem3d(d3,1,jx,1,iy,1,klev*2,'mod_ifs:d3')
    end if
    call getmem2d(ps,1,jx,1,iy,'mod_ifs:ps')
    call getmem2d(ts,1,jx,1,iy,'mod_ifs:ts')

    istatus = nf90_inq_varid(ncid,'lat',ivarid)
    call checkncerr(istatus,__FILE__,__LINE__, &
          'Missing lat variable in file '//trim(pathaddname))
    istatus = nf90_get_var(ncid,ivarid,glat)
    call checkncerr(istatus,__FILE__,__LINE__, &
          'Error reading lat variable in file '//trim(pathaddname))
    istatus = nf90_inq_varid(ncid,'lon',ivarid)
    call checkncerr(istatus,__FILE__,__LINE__, &
          'Missing longitude variable in file '//trim(pathaddname))
    istatus = nf90_get_var(ncid,ivarid,glon)
    call checkncerr(istatus,__FILE__,__LINE__, &
          'Error reading lon variable in file '//trim(pathaddname))
    istatus = nf90_inq_varid(ncid,'lev',ivarid)
    call checkncerr(istatus,__FILE__,__LINE__, &
          'Missing lev variable in file '//trim(pathaddname))
    istatus = nf90_get_var(ncid,ivarid,slev)
    call checkncerr(istatus,__FILE__,__LINE__, &
          'Error reading lev variable in file '//trim(pathaddname))
    istatus = nf90_close(ncid)
    call checkncerr(istatus,__FILE__,__LINE__, &
          'Error close file '//trim(pathaddname))
    !
    ! Find window to read
    !
    call get_window(glat,glon,xlat,xlon,i_band,gdomain)

    ghelp(1:jlat) = glat
    jlat = gdomain%nj
    call getmem1d(glat,1,jlat,'mod_ifs:glat')
    glat = ghelp(gdomain%jgstart:gdomain%jgstop)
    ghelp(1:ilon) = glon
    ilon = sum(gdomain%ni)
    call getmem1d(glon,1,ilon,'mod_ifs:glon')
    glon(1:gdomain%ni(1)) = ghelp(gdomain%igstart(1):gdomain%igstop(1))
    if ( gdomain%ntiles == 2 ) then
      glon(gdomain%ni(1)+1:ilon) = ghelp(gdomain%igstart(2):gdomain%igstop(2))
    end if

    print *, glon
    print *, glat

    call h_interpolator_create(cross_hint,glat,glon,xlat,xlon)
    if ( idynamic == 3 ) then
      call h_interpolator_create(udot_hint,glat,glon,ulat,ulon)
      call h_interpolator_create(vdot_hint,glat,glon,vlat,vlon)
    else
      call h_interpolator_create(udot_hint,glat,glon,dlat,dlon)
    end if

    call getmem3d(b2,1,ilon,1,jlat,1,klev*3,'mod_ifs:b2')
    call getmem3d(d2,1,ilon,1,jlat,1,klev*2,'mod_ifs:d2')
    call getmem2d(xps,1,ilon,1,jlat,'mod_ifs:xps')
    call getmem2d(xts,1,ilon,1,jlat,'mod_ifs:xts')
    !
    ! Set up pointers
    !
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
  end subroutine init_ifs

  subroutine get_ifs(idate)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate
    !
    ! Read data at idate
    !
    call ifs6hour(dattyp,idate,globidate1)
    write (stdout,*) 'READ IN fields at DATE:' , tochar(idate)
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
    call h_interpolate_cont(cross_hint,xps,ps)
    call h_interpolate_cont(cross_hint,xts,ts)
    !
    ! Rotate u-v fields after horizontal interpolation
    !
    if ( idynamic == 3 ) then
      call pju%wind_rotate(u3,v3u)
      call pju%wind_rotate(u3v,v3)
    else
      call pjd%wind_rotate(u3,v3)
    end if
    !
    ! Invert vertical order top -> bottom for RegCM convention
    !
!$OMP SECTIONS
!$OMP SECTION
    call top2btm(t3)
!$OMP SECTION
    call top2btm(q3)
!$OMP SECTION
    call top2btm(h3)
!$OMP SECTION
    call top2btm(u3)
!$OMP SECTION
    call top2btm(v3)
!$OMP END SECTIONS
    !
    ! Vertical interpolation
    ! New calculation of p* on rcm topography.
    !
    if ( idynamic == 3 ) then
      call ucrs2dot(h3u,h3,jx,iy,klev,i_band)
      call vcrs2dot(h3v,h3,jx,iy,klev,i_crm)
      call intzps(ps4,topogm,t3,h3,pss,sigmar,xlat,julianday(idate),jx,iy,klev)
    else
      call intgtb(pa,za,tlayer,topogm,t3,h3,pss,sigmar,jx,iy,klev)
      call intpsn(ps4,topogm,pa,za,tlayer,ptop,jx,iy)
      call crs2dot(pd4,ps4,jx,iy,i_band,i_crm)
    end if

    ts4 = ts

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
      call intz1(t4,t3,z0,h3,topogm,jx,iy,kz,klev,0.6_rkx,0.85_rkx,0.5_rkx)
!$OMP SECTION
      call intz1(q4,q3,z0,h3,topogm,jx,iy,kz,klev,0.7_rkx,0.7_rkx,0.4_rkx)
!$OMP END SECTIONS
    else
!$OMP SECTIONS
!$OMP SECTION
      call intv1(u4,u3,pd4,sigmah,pss,sigmar,ptop,jx,iy,kz,klev,1)
!$OMP SECTION
      call intv1(v4,v3,pd4,sigmah,pss,sigmar,ptop,jx,iy,kz,klev,1)
!$OMP SECTION
      call intv2(t4,t3,ps4,sigmah,pss,sigmar,ptop,jx,iy,kz,klev)
!$OMP SECTION
      call intv1(q4,q3,ps4,sigmah,pss,sigmar,ptop,jx,iy,kz,klev,1)
!$OMP END SECTIONS
    end if
    q4 = d_10**q4
  end subroutine get_ifs

  subroutine ifs6hour(dattyp,idate,idate0)
    implicit none
    character(len=5) , intent(in) :: dattyp
    type(rcm_time_and_date) , intent(in) :: idate , idate0
    integer(ik4) :: i , it , j , kkrec , istatus , ivar
    integer(ik4) :: timid
    character(len=64) :: inname
    character(len=256) :: pathaddname
    character(len=1) , dimension(7) :: varname
    character(len=64) :: cunit , ccal
    real(rkx) :: xadd , xscale
    integer(ik4) :: year , month , day , hour , monthp1
    integer(ik4) , save :: lastmonth , lastyear
    type(rcm_time_and_date) :: xdate
    type(rcm_time_interval) :: tdif
    !
    ! This is the latitude, longitude dimension of the grid to be read.
    ! This corresponds to the lat and lon dimension variables in the
    ! netCDF file.
    !
    data varname /'t' , 'z' , 'u' , 'v' , 'q' , 'lnsp' , 'skt' /

    call split_idate(idate,year,month,day,hour)
    write(inname,'(a,i0.4,i0.2,i0.2,i0.2,a)') &
           'IFS_',year,month,day,hour,'+000.nc'
    pathaddname = trim(inpglob)//pthsep//'IFS'//pthsep//inname
    istatus = nf90_open(pathaddname,nf90_nowrite,ncin)
    call checkncerr(istatus,__FILE__,__LINE__, &
                     'Error open file '//trim(pathaddname))
    do kkrec = 1 , size(varname)
      istatus = nf90_inq_varid(ncin,varname(kkrec), ivar5(kkrec))
      call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error find var '//varname(kkrec))
    end do

    it = 1

    do kkrec = 1 , size(varname)
      ivar = ivar5(kkrec)
      if ( kkrec == 1 ) then
        call getwork3(kkrec,tvar)
      else if ( kkrec == 2 ) then
        call getwork3(kkrec,hvar)
        hvar = hvar/9.80616_rk4
      else if ( kkrec == 3 ) then
        call getwork3(kkrec,qvar)
        call sph2mxr(qvar,ilon,jlat,klev)
        qvar = log10(qvar)
      else if ( kkrec == 4 ) then
        call getwork3(kkrec,uvar)
      else if ( kkrec == 5 ) then
        call getwork3(kkrec,vvar)
      else if ( kkrec == 6 ) then
        call getwork2(kkrec,xps)
        xps = exp(xps)
      else if ( kkrec == 7 ) then
        call getwork2(kkrec,xts)
      end if
    end do

    contains

      subroutine getwork3(irec,var)
        implicit none
        real(rkx) , pointer , intent(inout) , dimension(:,:,:) :: var
        integer(ik4) , intent(in) :: irec
        integer(ik4) :: itile , iti , itf
        integer(ik4) , dimension(4) :: icount , istart
        istart(3) = 1
        icount(3) = klev
        istart(4) = it
        icount(4) = 1
        iti = 1
        do itile = 1 , gdomain%ntiles
          istart(1) = gdomain%igstart(itile)
          icount(1) = gdomain%ni(itile)
          ! Latitudes are reversed in original file
          istart(2) = gdomain%jgstart
          icount(2) = gdomain%nj
          itf = iti + gdomain%ni(itile) - 1
          istatus = nf90_get_var(ncin,ivar,var(iti:itf,:,:),istart,icount)
          call checkncerr(istatus,__FILE__,__LINE__, &
                          'Error read var '//varname(irec))
          iti = iti + gdomain%ni(itile)
        end do
      end subroutine getwork3

      subroutine getwork2(irec,var)
        implicit none
        real(rkx) , pointer , intent(inout) , dimension(:,:) :: var
        integer(ik4) , intent(in) :: irec
        integer(ik4) :: itile , iti , itf
        integer(ik4) , dimension(3) :: icount , istart
        istart(3) = it
        icount(3) = 1
        iti = 1
        do itile = 1 , gdomain%ntiles
          istart(1) = gdomain%igstart(itile)
          icount(1) = gdomain%ni(itile)
          ! Latitudes are reversed in original file
          istart(2) = gdomain%jgstart
          icount(2) = gdomain%nj
          itf = iti + gdomain%ni(itile) - 1
          istatus = nf90_get_var(ncin,ivar,var(iti:itf,:),istart,icount)
          call checkncerr(istatus,__FILE__,__LINE__, &
                          'Error read var '//varname(irec))
          iti = iti + gdomain%ni(itile)
        end do
      end subroutine getwork2

    end subroutine ifs6hour

  subroutine conclude_ifs
    implicit none
    call h_interpolator_destroy(cross_hint)
    call h_interpolator_destroy(udot_hint)
    if ( idynamic == 3 ) then
      call h_interpolator_destroy(vdot_hint)
    end if
  end subroutine conclude_ifs

end module mod_ifs
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
