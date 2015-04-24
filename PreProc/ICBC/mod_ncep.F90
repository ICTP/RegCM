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

  integer(ik4) :: klev
  integer(ik4) :: jlat
  integer(ik4) :: ilon

  real(rk8) , pointer , dimension(:) :: glat , glat1
  real(rk8) , pointer , dimension(:) :: glon
  real(rk8) , pointer , dimension(:) :: sigmar , sigma1

  real(rk8) , pointer , dimension(:,:,:) :: b2
  real(rk8) , pointer , dimension(:,:,:) :: d2
  real(rk8) , pointer , dimension(:,:,:) :: b3
  real(rk8) , pointer , dimension(:,:,:) :: d3
  !
  ! The data are packed into short integers (INTEGER*2).  The array
  ! work will be used to hold the packed integers.
  !
  integer(2) , pointer , dimension(:,:,:) :: work

  real(rk8) , pointer :: u3(:,:,:) , v3(:,:,:)
  real(rk8) , pointer :: h3(:,:,:) , q3(:,:,:) , t3(:,:,:)
  real(rk8) , pointer :: uvar(:,:,:) , vvar(:,:,:)
  real(rk8) , pointer :: hvar(:,:,:) , rhvar(:,:,:) , tvar(:,:,:)

  integer(ik4) :: year , month , day , hour
  integer(ik4) :: itcfs = 0

  public :: getncep , headernc

  contains

  subroutine getncep(idate)
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
    call bilinx2(b3,b2,xlon,xlat,glon,glat,ilon,jlat,jx,iy,klev*3)
    call bilinx2(d3,d2,dlon,dlat,glon,glat,ilon,jlat,jx,iy,klev*2)
    !
    ! Rotate U-V fields after horizontal interpolation
    !
    call uvrot4(u3,v3,dlon,dlat,clon,clat,xcone,jx,iy,klev,plon,plat,iproj)
    !
    ! Vertical interpolation
    !
    ! New calculation of P* on rcm topography.
    !
    call intgtb(pa,za,tlayer,topogm,t3,h3,sigmar,jx,iy,klev)
    call intpsn(ps4,topogm,pa,za,tlayer,ptop,jx,iy)
    call crs2dot(pd4,ps4,jx,iy,i_band)
    !
    ! Determine surface temps on rcm topography.
    ! Interpolation from pressure levels
    !
    call intv3(ts4,t3,ps4,sigmar,ptop,jx,iy,klev)

    call readsst(ts4,idate)
    !
    ! Interpolate U, V, T, and Q.
    !
    call intv1(u4,u3,pd4,sigmah,sigmar,ptop,jx,iy,kz,klev)
    call intv1(v4,v3,pd4,sigmah,sigmar,ptop,jx,iy,kz,klev)
    call intv2(t4,t3,ps4,sigmah,sigmar,ptop,jx,iy,kz,klev)
    call intv1(q4,q3,ps4,sigmah,sigmar,ptop,jx,iy,kz,klev)
    call humid2(t4,q4,ps4,ptop,sigmah,jx,iy,kz)
  end subroutine getncep

  subroutine cfs6hour
    use netcdf
    implicit none
    integer(ik4) :: i , j , k , inet , it , kkrec , istatus
    character(len=256) :: pathaddname
    character(len=5) , dimension(5) :: varname
    real(rk8) :: xadd , xscale
    integer(ik4) , dimension(4) :: icount , istart
    integer(ik4) , dimension(5) , save :: inet5 , ivar5
    real(rk8) , dimension(5) , save :: xoff , xscl
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
      call checkncerr(istatus,__FILE__,__LINE__,'Variable '//varname(kkrec)// &
                    'read error in file'//trim(pathaddname))
      if ( kkrec == 1 ) then
        do k = 1 , klev
          do j = 1 , jlat
            do i = 1 , ilon
              tvar(i,j,k) = dble(work(i,j,k))*xscale+xadd
            end do
          end do
        end do
      else if ( kkrec == 2 ) then
        do k = 1 , klev
          do j = 1 , jlat
            do i = 1 , ilon
              hvar(i,j,k) = dble(work(i,j,k))*xscale+xadd
            end do
          end do
        end do
      else if ( kkrec == 3 ) then
        do k = 1 , klev
          do j = 1 , jlat
            do i = 1 , ilon
              rhvar(i,j,k) = dmin1((dble(work(i,j,k))* &
                            xscale+xadd)*0.01D0,1.D0)
            end do
          end do
        end do
      else if ( kkrec == 4 ) then
        do k = 1 , klev
          do j = 1 , jlat
            do i = 1 , ilon
              uvar(i,j,k) = dble(work(i,j,k))*xscale+xadd
            end do
          end do
        end do
      else if ( kkrec == 5 ) then
        do k = 1 , klev
          do j = 1 , jlat
            do i = 1 , ilon
              vvar(i,j,k) = dble(work(i,j,k))*xscale+xadd
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
    character(len=256) :: pathaddname
    character(len=5) , dimension(5) :: varname
    real(rk8) :: xadd , xscale
    integer(ik4) , dimension(4) :: icount , istart
    integer(ik4) , dimension(5) , save :: inet5 , ivar5
    real(rk8) , dimension(5) , save :: xoff , xscl
    data varname/'air' , 'hgt' , 'rhum' , 'uwnd' , 'vwnd'/
    !
    xadd = 0.0D0
    xscale = 1.0D0
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
      call checkncerr(istatus,__FILE__,__LINE__,'Variable '//varname(kkrec)// &
                    'read error in file'//trim(pathaddname))
      xscale = xscl(kkrec)
      xadd = xoff(kkrec)
      do ilev = 1 , nlev
        if ( kkrec == 1 ) then
          do j = 1 , jlat
            do i = 1 , ilon
              tvar(i,jlat+1-j,ilev) = dble(work(i,j,ilev))*xscale+xadd
            end do
          end do
        else if ( kkrec == 2 ) then
          do j = 1 , jlat
            do i = 1 , ilon
              hvar(i,jlat+1-j,ilev) = dble(work(i,j,ilev))*xscale+xadd
            end do
          end do
        else if ( kkrec == 3 ) then
          do j = 1 , jlat
            do i = 1 , ilon
              rhvar(i,jlat+1-j,ilev) = dmin1((dble(work(i,j,ilev))* &
                            xscale+xadd)*0.01D0,1.D0)
            end do
          end do
        else if ( kkrec == 4 ) then
          do j = 1 , jlat
            do i = 1 , ilon
              uvar(i,jlat+1-j,ilev) = dble(work(i,j,ilev))*xscale+xadd
            end do
          end do
        else if ( kkrec == 5 ) then
          do j = 1 , jlat
            do i = 1 , ilon
              vvar(i,jlat+1-j,ilev) = dble(work(i,j,ilev))*xscale+xadd
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

  subroutine headernc
    use netcdf
    implicit none

    integer(ik4) :: j , k , year , month , day , hour
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
    call checkncerr(istatus,__FILE__,__LINE__,'Error opening '//trim(inpfile))
    istatus = nf90_inq_dimid(inet,'lon',iddim)
    if ( istatus /= nf90_noerr ) then
      istatus = nf90_inq_dimid(inet,'longitude',iddim)
    end if
    call checkncerr(istatus,__FILE__,__LINE__,'Error find dim lon')
    istatus = nf90_inquire_dimension(inet,iddim, len=ilon)
    call checkncerr(istatus,__FILE__,__LINE__,'Error inquire dim lon')
    istatus = nf90_inq_dimid(inet,'lat',iddim)
    if ( istatus /= nf90_noerr ) then
      istatus = nf90_inq_dimid(inet,'latitude',iddim)
    end if
    call checkncerr(istatus,__FILE__,__LINE__,'Error find dim lat')
    istatus = nf90_inquire_dimension(inet,iddim, len=jlat)
    call checkncerr(istatus,__FILE__,__LINE__,'Error inquire dim lat')
    istatus = nf90_inq_dimid(inet,'level',iddim)
    call checkncerr(istatus,__FILE__,__LINE__,'Error find dim level')
    istatus = nf90_inquire_dimension(inet,iddim, len=klev)
    call checkncerr(istatus,__FILE__,__LINE__,'Error inquire dim level')

    call getmem1d(glon,1,ilon,'mod_ncep:glon')
    call getmem1d(glat,1,jlat,'mod_ncep:glat')
    call getmem1d(glat1,1,jlat,'mod_ncep:glat')
    call getmem1d(sigmar,1,klev,'mod_ncep:sigmar')
    call getmem1d(sigma1,1,klev,'mod_ncep:sigma1')

    istatus = nf90_inq_varid(inet,'level',idv)
    call checkncerr(istatus,__FILE__,__LINE__,'Error find var level')
    istatus = nf90_get_var(inet,idv,sigma1)
    call checkncerr(istatus,__FILE__,__LINE__,'Error read level')
    ! Invert levels
    do k = 1 , klev
      sigmar(k) = sigma1(klev-k+1)/1000.0D0
    end do
    !
    ! INITIAL GLOBAL GRID-POINT LONGITUDE & LATITUDE
    !
    istatus = nf90_inq_varid(inet,'lon',idv)
    if ( istatus /= nf90_noerr ) then
      istatus = nf90_inq_varid(inet,'longitude',idv)
    end if
    call checkncerr(istatus,__FILE__,__LINE__,'Error find var lon')
    istatus = nf90_get_var(inet,idv,glon)
    call checkncerr(istatus,__FILE__,__LINE__,'Error read lon')
    istatus = nf90_inq_varid(inet,'lat',idv)
    if ( istatus /= nf90_noerr ) then
      istatus = nf90_inq_varid(inet,'latitude',idv)
    end if
    call checkncerr(istatus,__FILE__,__LINE__,'Error find var lon')
    istatus = nf90_get_var(inet,idv,glat1)
    call checkncerr(istatus,__FILE__,__LINE__,'Error read lat')
    istatus = nf90_close(inet)
    call checkncerr(istatus,__FILE__,__LINE__,'Cannot close file')

    write(stdout,*) 'Read static data'

    if ( dattyp(1:3) /= 'CFS' ) then
      do j = 1 , jlat
        glat(j) = glat1(jlat-j+1)
      end do
    else
      glat(:) = glat1(:)
    end if

    call getmem3d(work,1,ilon,1,jlat,1,klev,'mod_ncep:work')
    call getmem3d(b2,1,ilon,1,jlat,1,klev*3,'mod_ncep:b3')
    call getmem3d(d2,1,ilon,1,jlat,1,klev*2,'mod_ncep:d3')
    call getmem3d(b3,1,jx,1,iy,1,klev*3,'mod_ncep:b3')
    call getmem3d(d3,1,jx,1,iy,1,klev*2,'mod_ncep:d3')

    ! Set up pointers

    u3 => d3(:,:,1:klev)
    v3 => d3(:,:,klev+1:2*klev)
    t3 => b3(:,:,1:klev)
    h3 => b3(:,:,klev+1:2*klev)
    q3 => b3(:,:,2*klev+1:3*klev)
    uvar => d2(:,:,1:klev)
    vvar => d2(:,:,klev+1:2*klev)
    tvar => b2(:,:,1:klev)
    hvar => b2(:,:,klev+1:2*klev)
    rhvar => b2(:,:,2*klev+1:3*klev)
  end subroutine headernc

end module mod_ncep
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
