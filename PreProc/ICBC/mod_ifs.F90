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

  integer(ik4) :: jlat , ilon , klev , hynlev

  real(rkx) , pointer , dimension(:,:,:) :: b3
  real(rkx) , pointer , dimension(:,:,:) :: d3
  real(rkx) , pointer , dimension(:,:,:) :: d3u , d3v

  real(rkx) , pointer :: u3(:,:,:) , v3(:,:,:)
  real(rkx) , pointer :: u3v(:,:,:) , v3u(:,:,:)
  real(rkx) , pointer :: p3(:,:,:) , pd3(:,:,:)
  real(rkx) , pointer :: q3(:,:,:) , t3(:,:,:) , z3(:,:,:)
  real(rkx) , pointer :: z3u(:,:,:) , z3v(:,:,:)
  real(rkx) , pointer :: uvar(:,:,:) , vvar(:,:,:) , pvar(:,:,:)
  real(rkx) , pointer :: qvar(:,:,:) , tvar(:,:,:)
  real(rkx) , pointer :: p_out(:,:,:) , pd_out(:,:,:)
  real(rkx) , pointer :: topou(:,:) , topov(:,:)
  real(rkx) , pointer :: xps(:,:) , xts(:,:) , xzs(:,:)
  real(rkx) , pointer :: ps(:,:) , ts(:,:) , zs(:,:)

  real(rkx) , pointer , dimension(:,:,:) :: b2
  real(rkx) , pointer , dimension(:,:,:) :: d2
  real(rkx) , pointer , dimension(:) :: glat
  real(rkx) , pointer , dimension(:) :: ghelp
  real(rkx) , pointer , dimension(:) :: glon
  integer(ik4) , pointer , dimension(:) :: slev
  real(rkx) , pointer , dimension(:) :: hyam , hybm
  real(rkx) :: pss

  integer(ik4) :: ncin
  integer(ik4) , parameter :: nrvar = 7
  character(len=4) , dimension(nrvar) , parameter :: varname = &
           ['t   ' , 'q   ' , 'u   ' , 'v   ' , 'lnsp' , 'st  ' , 'z   ']
  integer(ik4) , dimension(nrvar) :: ivar5

  type(global_domain) :: gdomain
  type(h_interpolator) :: cross_hint , udot_hint , vdot_hint

  public :: init_ifs , get_ifs , conclude_ifs

  contains

  subroutine init_ifs
    implicit none
    integer(ik4) :: year , month , day , hour
    character(len=256) :: pathaddname
    integer(ik4) :: istatus , ncid , ivarid , idimid
    integer(ik4) :: i , j , k
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
    istatus = nf90_inq_dimid(ncid,'nhym',idimid)
    call checkncerr(istatus,__FILE__,__LINE__, &
          'Missing nhym dimension in file '//trim(pathaddname))
    istatus = nf90_inquire_dimension(ncid,idimid,len=hynlev)
    call checkncerr(istatus,__FILE__,__LINE__, &
          'Error reading levelist dimelen in file '//trim(pathaddname))
    !
    ! Allocate working space
    !
    call getmem1d(slev,1,klev,'mod_ifs:slev')
    call getmem1d(hyam,1,hynlev,'mod_ifs:hyam')
    call getmem1d(hybm,1,hynlev,'mod_ifs:hybm')
    call getmem1d(glat,1,jlat,'mod_ifs:glat')
    call getmem1d(glon,1,ilon,'mod_ifs:glon')
    call getmem1d(ghelp,1,max(jlat,ilon),'mod_ifs:ghelp')
    call getmem3d(b3,1,jx,1,iy,1,klev*3,'mod_ifs:b3')
    if ( idynamic == 3 ) then
      call getmem3d(d3u,1,jx,1,iy,1,klev*2,'mod_ifs:d3u')
      call getmem3d(d3v,1,jx,1,iy,1,klev*2,'mod_ifs:d3v')
      call getmem3d(z3,1,jx,1,iy,1,klev,'mod_ifs:z3')
      call getmem3d(z3u,1,jx,1,iy,1,klev,'mod_ifs:z3u')
      call getmem3d(z3v,1,jx,1,iy,1,klev,'mod_ifs:z3v')
      call getmem2d(topou,1,jx,1,iy,'mod_ifs:topou')
      call getmem2d(topov,1,jx,1,iy,'mod_ifs:topov')
    else
      call getmem3d(d3,1,jx,1,iy,1,klev*2,'mod_ifs:d3')
    end if
    call getmem2d(ps,1,jx,1,iy,'mod_ifs:ps')
    call getmem2d(ts,1,jx,1,iy,'mod_ifs:ts')
    call getmem2d(zs,1,jx,1,iy,'mod_ifs:zs')

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
    istatus = nf90_inq_varid(ncid,'hyam',ivarid)
    call checkncerr(istatus,__FILE__,__LINE__, &
          'Missing hyam variable in file '//trim(pathaddname))
    istatus = nf90_get_var(ncid,ivarid,hyam)
    call checkncerr(istatus,__FILE__,__LINE__, &
          'Error reading hyam variable in file '//trim(pathaddname))
    istatus = nf90_inq_varid(ncid,'hybm',ivarid)
    call checkncerr(istatus,__FILE__,__LINE__, &
          'Missing hybm variable in file '//trim(pathaddname))
    istatus = nf90_get_var(ncid,ivarid,hybm)
    call checkncerr(istatus,__FILE__,__LINE__, &
          'Error reading hybm variable in file '//trim(pathaddname))
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
    call getmem2d(xzs,1,ilon,1,jlat,'mod_ifs:xzs')
    if ( idynamic /= 3 ) then
      call getmem3d(pd3,1,jx,1,iy,1,klev,'mod_nest:pd3')
      call getmem3d(p_out,1,jx,1,iy,1,kz,'mod_nest:p_out')
      call getmem3d(pd_out,1,jx,1,iy,1,kz,'mod_nest:pd_out')
    end if

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
    q3 => b3(:,:,1*klev+1:2*klev)
    p3 => b3(:,:,2*klev+1:3*klev)
    uvar => d2(:,:,1:klev)
    vvar => d2(:,:,klev+1:2*klev)
    tvar => b2(:,:,1:klev)
    qvar => b2(:,:,1*klev+1:2*klev)
    pvar => b2(:,:,2*klev+1:3*klev)
    if ( idynamic == 3 ) then
      call ucrs2dot(zud4,z0,jx,iy,kz,i_band)
      call vcrs2dot(zvd4,z0,jx,iy,kz,i_crm)
      call ucrs2dot(topou,topogm,jx,iy,i_band)
      call vcrs2dot(topov,topogm,jx,iy,i_crm)
    else if ( idynamic == 2 ) then
      do k = 1 , kz
        do i = 1 , iy
          do j = 1 , jx
            p_out(j,i,k) = ps0(j,i)*d_r1000*sigmah(k) + ptop
          end do
        end do
      end do
      call crs2dot(pd_out,p_out,jx,iy,kz,i_band,i_crm)
    end if
  end subroutine init_ifs

  subroutine get_ifs(idate)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate
    integer(ik4) :: i , j , k
    !
    ! Read data at idate
    !
    call ifs6hour(idate,globidate1)
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
    call h_interpolate_cont(cross_hint,xzs,zs)
    ps = ps * d_r1000
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
    ! Vertical interpolation
    ! New calculation of p* on rcm topography.
    !
    call intpsn(ps4,topogm,ps,zs,ts,ptop,jx,iy)

    if ( idynamic == 3 ) then
      z3(:,:,klev) = zs(:,:) + log(ps(:,:)/p3(:,:,klev))*rovg* &
                  t3(:,:,klev)*(d_one+ep1*q3(:,:,klev))
      do k = klev-1, 1 , -1
        z3(:,:,k) = z3(:,:,k+1) + log(p3(:,:,k+1)/p3(:,:,k))*rovg* &
                d_half * (t3(:,:,k+1)*(d_one+ep1*q3(:,:,k+1)) + &
                          t3(:,:,k)*(d_one+ep1*q3(:,:,k)))
      end do
      call ucrs2dot(z3u,z3,jx,iy,klev,i_band)
      call vcrs2dot(z3v,z3,jx,iy,klev,i_crm)
      call intz3(ts4,t3,z3,topogm,jx,iy,klev,0.0_rkx,0.05_rkx,0.05_rkx)
    else
      if ( idynamic == 1 ) then
        do k = 1 , kz
          do i = 1 , iy
            do j = 1 , jx
              p_out(j,i,k) = ps4(j,i) * sigmah(k) + ptop
            end do
          end do
        end do
        call crs2dot(pd_out,p_out,jx,iy,kz,i_band,i_crm)
      end if
      call crs2dot(pd4,ps4,jx,iy,i_band,i_crm)
      call crs2dot(pd3,p3,jx,iy,klev,i_band,i_crm)
      ps = ps4 + ptop
      call intp3(ts4,t3,p3,ps,jx,iy,klev,0.0_rkx,0.05_rkx,0.05_rkx)
    end if

    where ( mask == 0 )
      ts4 = ts
    end where
    !
    ! Interpolate U, V, T, and Q.
    !
    if ( idynamic == 3 ) then
!$OMP SECTIONS
!$OMP SECTION
      call intz1(u4,u3,zud4,z3u,topou,jx,iy,kz,klev,0.6_rkx,0.2_rkx,0.2_rkx)
!$OMP SECTION
      call intz1(v4,v3,zvd4,z3v,topov,jx,iy,kz,klev,0.6_rkx,0.2_rkx,0.2_rkx)
!$OMP SECTION
      call intz1(t4,t3,z0,z3,topogm,jx,iy,kz,klev,0.6_rkx,0.5_rkx,0.85_rkx)
!$OMP SECTION
      call intz1(q4,q3,z0,z3,topogm,jx,iy,kz,klev,0.7_rkx,0.4_rkx,0.7_rkx)
!$OMP END SECTIONS
    else
!$OMP SECTIONS
!$OMP SECTION
      call intp1(u4,u3,pd_out,pd3,jx,iy,kz,klev,0.6_rkx,0.2_rkx,0.2_rkx)
!$OMP SECTION
      call intp1(v4,v3,pd_out,pd3,jx,iy,kz,klev,0.6_rkx,0.2_rkx,0.2_rkx)
!$OMP SECTION
      call intp1(t4,t3,p_out,p3,jx,iy,kz,klev,0.6_rkx,0.85_rkx,0.5_rkx)
!$OMP SECTION
      call intp1(q4,q3,p_out,p3,jx,iy,kz,klev,0.7_rkx,0.7_rkx,0.4_rkx)
!$OMP END SECTIONS
    end if
  end subroutine get_ifs

  subroutine ifs6hour(idate,idate0)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate , idate0
    integer(ik4) :: k , it , kkrec , istatus , ivar
    character(len=64) :: inname
    character(len=256) :: pathaddname
    integer(ik4) :: year , month , day , hour
    !
    ! This is the latitude, longitude dimension of the grid to be read.
    ! This corresponds to the lat and lon dimension variables in the
    ! netCDF file.
    !

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
        call getwork3(kkrec,qvar)
        call sph2mxr(qvar,ilon,jlat,klev)
      else if ( kkrec == 3 ) then
        call getwork3(kkrec,uvar)
      else if ( kkrec == 4 ) then
        call getwork3(kkrec,vvar)
      else if ( kkrec == 5 ) then
        call getwork2(kkrec,xps)
        xps = exp(xps)
      else if ( kkrec == 6 ) then
        call getwork2(kkrec,xts)
      else if ( kkrec == 7 ) then
        call getwork2(kkrec,xzs)
        xzs = xzs / 9.80616_rkx
      end if
    end do

    do k = 1 , klev
      pvar(:,:,k) = d_r1000*(hyam(slev(k)) + xps(:,:) * hybm(slev(k))) ! cb
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
