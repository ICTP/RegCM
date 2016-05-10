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

module mod_ae_icbc

  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_memutil
  use mod_grid
  use mod_wrtoxd
  use mod_interp
  use mod_date
  use mod_nchelper
  use netcdf

  private
!
  integer(ik4) :: aeilon , aejlat , aeilev , aeitime
  real(rkx) , pointer , dimension(:) :: aet42lon
  real(rkx) , pointer , dimension(:) :: aet42lat
  real(rkx) , pointer , dimension(:) :: aet42hyam , aet42hybm
!
! Oxidant climatology variables
!
  real(rkx) :: p0 , r4pt
  real(rkx) , pointer , dimension(:,:) :: xps
  real(rkx) , pointer , dimension(:,:,:,:,:) :: aev2
  real(rkx) , pointer , dimension(:,:,:) :: xps2
  real(rkx) , pointer , dimension(:,:,:) :: xinp
  real(rkx) , pointer , dimension(:,:) :: paeid_3
  real(rkx) , pointer , dimension(:,:,:,:) :: aev3
  integer(ik4) :: iyear
  character(len=8) , dimension(4) :: scendir

  real(rkx) :: prcm , pmpi , pmpj
  integer(ik4) :: ncid , istatus , iscen

  public :: header_ae_icbc , get_ae_icbc , close_ae_icbc

  data ncid /-1/
  data scendir / 'RF', 'RCP26', 'RCP45', 'RCP85' /

  contains

  subroutine header_ae_icbc(idate)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate
    integer(ik4) :: ivarid , istatus , dimid , is
    integer(ik4) :: nyear , month , nday , nhour
    character(len=256) :: aefilename

    call split_idate(idate,nyear,month,nday,nhour)

    r4pt = real(ptop)

    iyear = nyear/10*10
    select case ( dattyp(4:5) )
      case ( '26' )
        iscen = 2
      case ( '45' )
        iscen = 3
      case ( '85' )
        iscen = 4
      case default
        iscen = 1
    end select
    write(aefilename,'(a,i0.4,a,i0.4,a)') &
       trim(inpglob)//pthsep//'AERGLOB'//pthsep// &
       trim(scendir(iscen))//pthsep//'aero_1.9x2.5_L26_', &
       iyear, '-', iyear+9, '.nc'

    write(stdout,*) 'Opening ',trim(aefilename)
    istatus = nf90_open(aefilename,nf90_nowrite, ncid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error open aeid file')

    istatus = nf90_inq_dimid(ncid,'lon',dimid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error find dim lon')
    istatus = nf90_inquire_dimension(ncid,dimid,len=aeilon)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error inquire lon')
    istatus = nf90_inq_dimid(ncid,'lat',dimid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error find dim lat')
    istatus = nf90_inquire_dimension(ncid,dimid,len=aejlat)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error inquire lat')
    istatus = nf90_inq_dimid(ncid,'lev',dimid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error find dim lev')
    istatus = nf90_inquire_dimension(ncid,dimid,len=aeilev)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error inquire lev')
    istatus = nf90_inq_dimid(ncid,'time',dimid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error find dim time')
    istatus = nf90_inquire_dimension(ncid,dimid,len=aeitime)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error inquire time')

    call getmem1d(aet42lon,1,aeilon,'mod_ae_icbc:aeilon')
    call getmem1d(aet42lat,1,aejlat,'mod_ae_icbc:aejlat')
    call getmem1d(aet42hyam,1,aeilev,'mod_ae_icbc:aet42hyam')
    call getmem1d(aet42hybm,1,aeilev,'mod_ae_icbc:aet42hybm')

    istatus = nf90_inq_varid(ncid,'lon',ivarid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error find var lon')
    istatus = nf90_get_var(ncid,ivarid,aet42lon)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read var lon')
    istatus = nf90_inq_varid(ncid,'lat',ivarid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error find var lat')
    istatus = nf90_get_var(ncid,ivarid,aet42lat)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read var lat')
    istatus = nf90_inq_varid(ncid,'hyam',ivarid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error find var hyam')
    istatus = nf90_get_var(ncid,ivarid,aet42hyam)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read var hyam')
    istatus = nf90_inq_varid(ncid,'hybm',ivarid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error find var hybm')
    istatus = nf90_get_var(ncid,ivarid,aet42hybm)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read var hybm')
    istatus = nf90_inq_varid(ncid,'P0',ivarid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error find var P0')
    istatus = nf90_get_var(ncid,ivarid,p0)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read var P0')
    p0 = p0*0.01

    call getmem2d(paeid_3,1,jx,1,iy,'mod_ae_icbc:paeid_3')
    call getmem2d(xps,1,aeilon,1,aejlat,'mod_ae_icbc:xps')
    call getmem3d(xps2,1,aeilon,1,aejlat,1,aeitime,'mod_ae_icbc:xps2')
    call getmem3d(xinp,1,aeilon,1,aejlat,1,aeilev,'mod_ae_icbc:xinp')
    call getmem4d(aev3,1,jx,1,iy,1,aeilev,1,naesp,'mod_ae_icbc:aev3')
    call getmem5d(aev2,1,aeilon,1,aejlat,1,aeilev,1,aeitime, &
                  1,naesp,'mod_ae_icbc:aev2')

    istatus = nf90_inq_varid(ncid,'PS',ivarid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error find var PS')
    istatus = nf90_get_var(ncid,ivarid,xps2)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read var PS')
    xps2 = xps2*0.01
    do is = 1 , naesp
      istatus = nf90_inq_varid(ncid,aespec(is),ivarid)
      if ( istatus == nf90_noerr ) then
        istatus = nf90_get_var(ncid,ivarid,aev2(:,:,:,:,is))
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error read var '//aespec(is))
      else
        aev2(:,:,:,:,is) = d_zero
      end if
    end do
  end subroutine header_ae_icbc

  subroutine get_ae_icbc(idate)
    implicit none

    type(rcm_time_and_date) , intent(in) :: idate
    integer(ik4) :: i , l , is , j , k , k0
    real(rkx) :: wt1 , wt2
    type(rcm_time_and_date) :: d1 , d2
    type(rcm_time_interval) :: t1 , tt
    integer(ik4) :: m1 , m2
    integer(ik4) :: ivarid , istatus
    integer(ik4) :: nyear , month , nday , nhour
    character(len=256) :: aefilename

    call split_idate(idate,nyear,month,nday,nhour)

    if ( nyear > iyear + 9 ) then
      iyear = nyear/10*10
      write(aefilename,'(a,i0.4,a,i0.4,a)') &
         trim(inpglob)//pthsep//'AERGLOB'//pthsep// &
         trim(scendir(iscen))//pthsep//'aero_1.9x2.5_L26_', &
         iyear, '-', iyear+9, '.nc'
      write(stdout,*) 'Opening ',trim(aefilename)
      istatus = nf90_open(aefilename,nf90_nowrite, ncid)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error open aeid file')
      istatus = nf90_inq_varid(ncid,'PS',ivarid)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find var PS')
      istatus = nf90_get_var(ncid,ivarid,xps2)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read var PS')
      xps2 = xps2*0.01
      do is = 1 , naesp
        istatus = nf90_inq_varid(ncid,aespec(is),ivarid)
        if ( istatus == nf90_noerr ) then
          istatus = nf90_get_var(ncid,ivarid,aev2(:,:,:,:,is))
          call checkncerr(istatus,__FILE__,__LINE__, &
                          'Error read var '//aespec(is))
        else
          aev2(:,:,:,:,is) = d_zero
        end if
      end do
    end if

    d1 = monfirst(idate)
    d2 = nextmon(d1)
    m1 = getmonth(d1)
    m2 = getmonth(d2)
    t1 = idate-d1
    tt = d2-d1
    wt1 = real(tohours(t1)/tohours(tt))
    wt2 = 1.0 - wt1

    do is = 1 , naesp
      do l = 1 , aeilev
        do i = 1 , aejlat
          do j = 1 , aeilon
            xinp(j,i,l) = aev2(j,i,l,m1,is)*wt2+aev2(j,i,l,m2,is)*wt1
          end do
        end do
      end do
      call bilinx2(aev3(:,:,:,is),xinp,xlon,xlat,aet42lon,aet42lat, &
                   aeilon,aejlat,jx,iy,aeilev)
    end do

    do i = 1 , aejlat
      do j = 1 , aeilon
        xps(j,i) = xps2(j,i,m1)*wt1+xps2(j,i,m2)*wt2
      end do
    end do

    call bilinx2(paeid_3,xps,xlon,xlat,aet42lon,aet42lat, &
                 aeilon,aejlat,jx,iy)

    do i = 1 , iy
      do j = 1 , jx
        do l = 1 , kz
          prcm=((paeid_3(j,i)*0.1-r4pt)*sigmah(l)+r4pt)*10.
          k0 = -1
          do k = aeilev , 1 , -1
            pmpi = aet42hyam(k)*p0+paeid_3(j,i)*aet42hybm(k)
            k0 = k
            if (prcm > pmpi) exit
          end do
          if (k0 == aeilev) then
            pmpj = aet42hyam(aeilev-1)*p0+paeid_3(j,i)*aet42hybm(aeilev-1)
            pmpi = aet42hyam(aeilev  )*p0+paeid_3(j,i)*aet42hybm(aeilev  )
            do is = 1 , naesp
              aev4(j,i,l,is) = aev3(j,i,aeilev,is) + &
                 (aev3(j,i,aeilev-1,is) - aev3(j,i,aeilev,is)) * &
                 (prcm-pmpi)/(pmpi-pmpj)
            end do
          else if (k0 >= 1) then
            pmpj = aet42hyam(k0  )*p0+paeid_3(j,i)*aet42hybm(k0  )
            pmpi = aet42hyam(k0+1)*p0+paeid_3(j,i)*aet42hybm(k0+1)
            wt1 = (prcm-pmpj)/(pmpi-pmpj)
            wt2 = 1.0 - wt1
            do is = 1 , naesp
              aev4(j,i,l,is) = aev3(j,i,k0+1,is)*wt1 + aev3(j,i,k0,is)*wt2
            end do
          end if
        end do
      end do
    end do

    call write_ae_icbc(idate)

  end subroutine get_ae_icbc

  subroutine close_ae_icbc
    use netcdf
    implicit none
    if ( ncid > 0 ) then
      istatus = nf90_close(ncid)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error close aeid file')
    end if
  end subroutine close_ae_icbc

end module mod_ae_icbc
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
