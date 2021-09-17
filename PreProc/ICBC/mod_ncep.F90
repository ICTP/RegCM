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

module mod_ncep

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

  private

  integer(ik4) :: klev
  integer(ik4) :: jlat
  integer(ik4) :: ilon

  real(rkx) , pointer , dimension(:) :: glat
  real(rkx) , pointer , dimension(:) :: glon
  real(rkx) , pointer , dimension(:) :: sigmar , plevs
  real(rkx) :: pss , pst

  real(rkx) , pointer , dimension(:,:,:) :: b2
  real(rkx) , pointer , dimension(:,:,:) :: d2
  real(rkx) , pointer , dimension(:,:,:) :: b3
  real(rkx) , pointer , dimension(:,:,:) :: d3
  real(rkx) , pointer , dimension(:,:,:) :: d3u
  real(rkx) , pointer , dimension(:,:,:) :: d3v
  real(rkx) , pointer , dimension(:,:,:) :: h3v , h3u
  real(rkx) , pointer , dimension(:,:) :: topou , topov
  !
  ! The data are packed into short integers (INTEGER*2).  The array
  ! work will be used to hold the packed integers.
  !
  integer(2) , pointer , dimension(:,:,:) :: work

  real(rkx) , pointer :: u3(:,:,:) , v3(:,:,:)
  real(rkx) , pointer :: u3v(:,:,:) , v3u(:,:,:)
  real(rkx) , pointer :: h3(:,:,:) , q3(:,:,:) , t3(:,:,:)
  real(rkx) , pointer :: uvar(:,:,:) , vvar(:,:,:)
  real(rkx) , pointer :: hvar(:,:,:) , rhvar(:,:,:) , tvar(:,:,:)

  integer(ik4) :: year , month , day , hour
  integer(ik4) :: itcfs = 0

  public :: get_ncep , init_ncep , conclude_ncep

  type(h_interpolator) :: cross_hint , udot_hint , vdot_hint


  contains

  subroutine init_ncep
    use netcdf
    implicit none

    integer(ik4) :: k , year , month , day , hour
    integer(ik4) :: istatus , inet , iddim , idv
    character(len=256) :: inpfile

    call split_idate(globidate1, year, month, day, hour)
    if ( dattyp(1:3) == 'CFS' ) then
      write(inpfile,'(a,i0.4,i0.2,i0.2,i0.2,a,i0.4,i0.2,i0.2,i0.2,a)') &
        trim(inpglob)//'/CFS/',year,month,day,hour, &
        '/'//dattyp(4:5)//'/PLEV/air.',year,month,day,hour,'.nc'
    else
      write (inpfile,'(a,i0.4,a,i0.4,a)') trim(inpglob)//'/'// &
           dattyp//'/',year , '/air.' , year,'.nc'
    end if
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
    istatus = nf90_inq_dimid(inet,'level',iddim)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error find dim level')
    istatus = nf90_inquire_dimension(inet,iddim, len=klev)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error inquire dim level')

    call getmem1d(glon,1,ilon,'mod_ncep:glon')
    call getmem1d(glat,1,jlat,'mod_ncep:glat')
    call getmem1d(sigmar,1,klev,'mod_ncep:sigmar')
    call getmem1d(plevs,1,klev,'mod_ncep:plevs')

    istatus = nf90_inq_varid(inet,'level',idv)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error find var level')
    istatus = nf90_get_var(inet,idv,plevs)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read level')
    do k = 1 , klev
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

    call getmem3d(work,1,ilon,1,jlat,1,klev,'mod_ncep:work')
    call getmem3d(b2,1,ilon,1,jlat,1,klev*3,'mod_ncep:b3')
    call getmem3d(d2,1,ilon,1,jlat,1,klev*2,'mod_ncep:d3')
    call getmem3d(b3,1,jx,1,iy,1,klev*3,'mod_ncep:b3')
    if ( idynamic == 3 ) then
      call getmem3d(d3u,1,jx,1,iy,1,klev*2,'mod_ncep:d3u')
      call getmem3d(d3v,1,jx,1,iy,1,klev*2,'mod_ncep:d3v')
      call getmem3d(h3u,1,jx,1,iy,1,klev,'mod_era5:h3u')
      call getmem3d(h3v,1,jx,1,iy,1,klev,'mod_era5:h3v')
      call getmem2d(topou,1,jx,1,iy,'mod_era5:topou')
      call getmem2d(topov,1,jx,1,iy,'mod_era5:topov')
    else
      call getmem3d(d3,1,jx,1,iy,1,klev*2,'mod_ncep:d3')
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
    rhvar => b2(:,:,2*klev+1:3*klev)
    if ( idynamic == 3 ) then
      call ucrs2dot(zud4,z0,jx,iy,kz,i_band)
      call vcrs2dot(zvd4,z0,jx,iy,kz,i_crm)
      call ucrs2dot(topou,topogm,jx,iy,i_band)
      call vcrs2dot(topov,topogm,jx,iy,i_crm)
    end if
  end subroutine init_ncep

  subroutine get_ncep(idate)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate

    call split_idate(idate,year,month,day,hour)
    if ( dattyp(1:3) == 'CFS' ) then
      call cfs6hour
    else
      call cdc6hour(idate)
    end if

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
                  xlat,yeardayfrac(idate),jx,iy,klev)
      call intz3(ts4,t3,h3,topogm,jx,iy,klev,0.6_rkx,0.5_rkx,0.85_rkx)
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
    call rh2mxr(t4,q4,ps4,ptop,sigmah,jx,iy,kz)
  end subroutine get_ncep

  subroutine cfs6hour
    use netcdf
    implicit none
    integer(ik4) :: i , j , k , inet , it , kkrec , istatus
    character(len=256) , save :: pathaddname
    character(len=5) , dimension(5) :: varname
    real(rkx) :: xadd , xscale
    integer(ik4) , dimension(4) :: icount , istart
    integer(ik4) , dimension(5) , save :: inet5 , ivar5
    real(rkx) , dimension(5) , save :: xoff , xscl
    integer(2) , dimension(5) , save :: xfil
    data varname/'air' , 'hgt' , 'rhum' , 'uwnd' , 'vwnd'/
    !
    if ( itcfs == 0 ) then
      do kkrec = 1 , 5
        write(pathaddname,'(a,i0.4,i0.2,i0.2,i0.2,a,i0.4,i0.2,i0.2,i0.2,a)') &
          trim(inpglob)//'/CFS/',year,month,day,hour, &
          '/'//dattyp(4:5)//'/PLEV/'//trim(varname(kkrec))// &
          '.',year,month,day,hour,'.nc'
        istatus = nf90_open(pathaddname,nf90_nowrite,inet5(kkrec))
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error opening '//trim(pathaddname))
        istatus = nf90_inq_varid(inet5(kkrec),varname(kkrec),ivar5(kkrec))
        call checkncerr(istatus,__FILE__,__LINE__, &
             'Variable '//varname(kkrec)//' error in file'//trim(pathaddname))
        istatus = nf90_get_att(inet5(kkrec),ivar5(kkrec), &
                              'scale_factor',xscl(kkrec))
        call checkncerr(istatus,__FILE__,__LINE__, &
             'Variable '//varname(kkrec)// &
             ':scale_factor in file'//trim(pathaddname))
        istatus = nf90_get_att(inet5(kkrec),ivar5(kkrec), &
                               'add_offset',xoff(kkrec))
        call checkncerr(istatus,__FILE__,__LINE__, &
              'Variable '//varname(kkrec)// &
              ':add_offset in file'//trim(pathaddname))
        write (stdout,*) inet5(kkrec) , trim(pathaddname) , &
                         xscl(kkrec) , xoff(kkrec)
        istatus = nf90_get_att(inet5(kkrec),ivar5(kkrec), &
                               '_FillValue',xfil(kkrec))
        call checkncerr(istatus,__FILE__,__LINE__, &
              'Variable '//varname(kkrec)// &
              ':_FillValue in file'//trim(pathaddname))
        write (stdout,*) inet5(kkrec) , trim(pathaddname) , &
                         xscl(kkrec) , xoff(kkrec)
        itcfs = 1
      end do
    end if

    it = itcfs
    itcfs = itcfs + 1

    do kkrec = 1 , 5
      xscale = xscl(kkrec)
      xadd = xoff(kkrec)
      istart(1) = 1
      icount(1) = ilon
      istart(2) = 1
      icount(2) = jlat
      istart(3) = 1
      icount(3) = klev
      istart(4) = it
      icount(4) = 1
      inet = inet5(kkrec)
      istatus = nf90_get_var(inet,ivar5(kkrec),work,istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Variable '//varname(kkrec)// &
                      'read error in file '//trim(pathaddname))
      if ( kkrec == 1 ) then
        do k = 1 , klev
          do j = 1 , jlat
            do i = 1 , ilon
              tvar(i,j,k) = real(work(i,j,k),rkx)*xscale+xadd
            end do
          end do
        end do
      else if ( kkrec == 2 ) then
        do k = 1 , klev
          do j = 1 , jlat
            do i = 1 , ilon
              hvar(i,j,k) = real(work(i,j,k),rkx)*xscale+xadd
            end do
          end do
        end do
      else if ( kkrec == 3 ) then
        do k = 1 , klev
          do j = 1 , jlat
            do i = 1 , ilon
              if ( work(i,j,k) /= xfil(kkrec) ) then
                rhvar(i,j,k) = (real(work(i,j,k),rkx)*xscale+xadd)*0.01_rkx
              else
                rhvar(i,j,k) = 0.01_rkx
              end if
            end do
          end do
        end do
      else if ( kkrec == 4 ) then
        do k = 1 , klev
          do j = 1 , jlat
            do i = 1 , ilon
              uvar(i,j,k) = real(work(i,j,k),rkx)*xscale+xadd
            end do
          end do
        end do
      else if ( kkrec == 5 ) then
        do k = 1 , klev
          do j = 1 , jlat
            do i = 1 , ilon
              vvar(i,j,k) = real(work(i,j,k),rkx)*xscale+xadd
            end do
          end do
        end do
      end if
    end do
  end subroutine cfs6hour

  subroutine cdc6hour(idate)
    use netcdf
    implicit none
    type(rcm_time_and_date) , intent (in) :: idate
    integer(ik4) :: i , ilev , inet , it , j , kkrec , k , nlev , istatus
    character(len=256) , save :: pathaddname
    character(len=5) , dimension(5) :: varname
    real(rkx) :: xadd , xscale
    integer(ik4) , dimension(4) :: icount , istart
    integer(ik4) , dimension(5) , save :: inet5 , ivar5
    real(rkx) , dimension(5) , save :: xoff , xscl
    data varname/'air' , 'hgt' , 'rhum' , 'uwnd' , 'vwnd'/
    !
    xadd = 0.0_rkx
    xscale = 1.0_rkx
    nlev = 0
    do kkrec = 1 , 5
      nlev = klev
      if ( dattyp == 'NNRP1' ) then
        if ( kkrec == 3 ) nlev = 8 ! Relative humidity has less levels
      end if
      if ( idate == globidate1 .or. &
           (lfdoyear(idate) .and. lmidnight(idate))) then
        write(pathaddname,'(a,i0.4,a,i0.4,a)') trim(inpglob)//'/'// &
            dattyp//'/',year , '/'//trim(varname(kkrec))//'.' , year,'.nc'
        istatus = nf90_open(pathaddname,nf90_nowrite,inet5(kkrec))
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error opening '//trim(pathaddname))
        istatus = nf90_inq_varid(inet5(kkrec),varname(kkrec),ivar5(kkrec))
        call checkncerr(istatus,__FILE__,__LINE__, &
             'Variable '//varname(kkrec)//' error in file'//trim(pathaddname))
        istatus = nf90_get_att(inet5(kkrec),ivar5(kkrec), &
                              'scale_factor',xscl(kkrec))
        call checkncerr(istatus,__FILE__,__LINE__, &
             'Variable '//varname(kkrec)// &
             ':scale_factor in file'//trim(pathaddname))
        istatus = nf90_get_att(inet5(kkrec),ivar5(kkrec), &
                               'add_offset',xoff(kkrec))
        call checkncerr(istatus,__FILE__,__LINE__, &
              'Variable '//varname(kkrec)// &
              ':add_offset in file'//trim(pathaddname))
        write (stdout,*) inet5(kkrec) , trim(pathaddname) , &
                         xscl(kkrec) , xoff(kkrec)
      end if

      it = (day-1)*4 + hour/6 + 1
      if ( month == 2 ) it = it + 31*4
      if ( month == 3 ) it = it + 59*4
      if ( month == 4 ) it = it + 90*4
      if ( month == 5 ) it = it + 120*4
      if ( month == 6 ) it = it + 151*4
      if ( month == 7 ) it = it + 181*4
      if ( month == 8 ) it = it + 212*4
      if ( month == 9 ) it = it + 243*4
      if ( month == 10 ) it = it + 273*4
      if ( month == 11 ) it = it + 304*4
      if ( month == 12 ) it = it + 334*4
      if ( mod(day,4) == 0 .and. month > 2 ) it = it + 4
      if ( mod(day,100) == 0 .and. month > 2 ) it = it - 4
      if ( mod(day,400) == 0 .and. month > 2 ) it = it + 4

      istart(1) = 1
      icount(1) = ilon
      istart(2) = 1
      icount(2) = jlat
      istart(3) = 1
      icount(3) = nlev
      istart(4) = it
      icount(4) = 1
      inet = inet5(kkrec)
      istatus = nf90_get_var(inet,ivar5(kkrec),work,istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Variable '//varname(kkrec)// &
                      'read error in file '//trim(pathaddname))
      xscale = xscl(kkrec)
      xadd = xoff(kkrec)
      do ilev = 1 , nlev
        if ( kkrec == 1 ) then
          do j = 1 , jlat
            do i = 1 , ilon
              tvar(i,j,ilev) = real(work(i,j,ilev),rkx)*xscale+xadd
            end do
          end do
        else if ( kkrec == 2 ) then
          do j = 1 , jlat
            do i = 1 , ilon
              hvar(i,j,ilev) = real(work(i,j,ilev),rkx)*xscale+xadd
            end do
          end do
        else if ( kkrec == 3 ) then
          do j = 1 , jlat
            do i = 1 , ilon
              rhvar(i,j,ilev) = min((real(work(i,j,ilev),rkx)* &
                            xscale+xadd)*0.01_rkx,1._rkx)
            end do
          end do
        else if ( kkrec == 4 ) then
          do j = 1 , jlat
            do i = 1 , ilon
              uvar(i,j,ilev) = real(work(i,j,ilev),rkx)*xscale+xadd
            end do
          end do
        else if ( kkrec == 5 ) then
          do j = 1 , jlat
            do i = 1 , ilon
              vvar(i,j,ilev) = real(work(i,j,ilev),rkx)*xscale+xadd
            end do
          end do
        end if
      end do
      if ( dattyp == 'NNRP1' ) then
        if ( kkrec == 3 ) then
          do k = nlev+1, klev
            do j = 1 , jlat
              do i = 1 , ilon
                rhvar(i,j,k) = rhvar(i,j,k-1)
              end do
            end do
          end do
        end if
      end if
    end do
  end subroutine cdc6hour

  subroutine conclude_ncep
    implicit none
    call h_interpolator_destroy(cross_hint)
    call h_interpolator_destroy(udot_hint)
    if ( idynamic == 3 ) then
      call h_interpolator_destroy(vdot_hint)
    end if
  end subroutine conclude_ncep

end module mod_ncep
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
