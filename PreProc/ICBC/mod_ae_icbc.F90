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

  use netcdf
  use mod_dynparam
  use mod_memutil
  use mod_grid
  use mod_wrtoxd
  use mod_interp
  use mod_date
  use mod_realkinds
  use mod_nchelper

  private
!
  integer , parameter :: aeilon = 144 , aejlat = 96 , aeilev = 26 , aeitime = 12
  real(sp) , dimension(aeilon) :: aet42lon
  real(sp) , dimension(aejlat) :: aet42lat
  real(sp) , dimension(aeilev) :: aet42hyam , aet42hybm
!
! Oxidant climatology variables
!
  real(sp) :: p0
  real(sp) , dimension(aeilon,aejlat) :: xps
  real(sp) , dimension(aeilon,aejlat) :: paeid_2
  real(sp) , pointer , dimension(:,:,:,:,:) :: aev2
  real(sp) , dimension(aeilon,aejlat,aeitime) :: xps2
  real(sp) , pointer , dimension(:,:) :: paeid_3
  real(sp) , pointer , dimension(:,:,:,:) :: aev3
  integer :: iyear
  character(len=8) , dimension(4) :: scendir

  real(sp) :: prcm , pmpi , pmpj
  integer :: ncid , istatus , iscen

  public :: header_ae_icbc , get_ae_icbc , close_ae_icbc

  data ncid /-1/
  data scendir / 'RF', 'RCP26', 'RCP45', 'RCP85' /

  contains

  subroutine header_ae_icbc(idate)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate
    integer :: ivarid , istatus , is
    integer :: nyear , month , nday , nhour
    character(len=256) :: aefilename

    call split_idate(idate,nyear,month,nday,nhour)
    call getmem2d(paeid_3,1,jx,1,iy,'mod_ae_icbc:paeid_3')
    call getmem4d(aev3,1,jx,1,iy,1,aeilev,1,naesp,'mod_ae_icbc:aev3')
    call getmem5d(aev2,1,aeilon,1,aejlat,1,aeilev,1,aeitime, &
                  1,naesp,'mod_ae_icbc:aev2')

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
    istatus = nf90_open(aefilename,nf90_nowrite, ncid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error open aeid file')

    istatus = nf90_inq_varid(ncid,'lon',ivarid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error find var lon')
    istatus = nf90_get_var(ncid,ivarid,aet42lon)
    call checkncerr(istatus,__FILE__,__LINE__,'Error read var lon')
    istatus = nf90_inq_varid(ncid,'lat',ivarid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error find var lat')
    istatus = nf90_get_var(ncid,ivarid,aet42lat)
    call checkncerr(istatus,__FILE__,__LINE__,'Error read var lat')
    istatus = nf90_inq_varid(ncid,'hyam',ivarid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error find var hyam')
    istatus = nf90_get_var(ncid,ivarid,aet42hyam)
    call checkncerr(istatus,__FILE__,__LINE__,'Error read var hyam')
    istatus = nf90_inq_varid(ncid,'hybm',ivarid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error find var hybm')
    istatus = nf90_get_var(ncid,ivarid,aet42hybm)
    call checkncerr(istatus,__FILE__,__LINE__,'Error read var hybm')
    istatus = nf90_inq_varid(ncid,'P0',ivarid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error find var P0')
    istatus = nf90_get_var(ncid,ivarid,p0)
    call checkncerr(istatus,__FILE__,__LINE__,'Error read var P0')
    istatus = nf90_inq_varid(ncid,'PS',ivarid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error find var PS')
    istatus = nf90_get_var(ncid,ivarid,xps2)
    call checkncerr(istatus,__FILE__,__LINE__,'Error read var PS')
    do is = 1 , naesp
      istatus = nf90_inq_varid(ncid,aespec(is),ivarid)
      call checkncerr(istatus,__FILE__,__LINE__,'Error find var '//aespec(is))
      istatus = nf90_get_var(ncid,ivarid,aev2(:,:,:,:,is))
      call checkncerr(istatus,__FILE__,__LINE__,'Error read var '//aespec(is))
    end do
  end subroutine header_ae_icbc

  subroutine get_ae_icbc(idate)
    implicit none
!
    type(rcm_time_and_date) , intent(in) :: idate
    integer :: i , l , is , j , k , k0
    real(sp) , dimension(aeilon,aejlat,aeilev) :: xinp
    real(sp) :: wt1 , wt2 , r4pt
    type(rcm_time_and_date) :: d1 , d2
    type(rcm_time_interval) :: t1 , tt
    integer :: m1 , m2
    integer :: ivarid , istatus
    integer :: nyear , month , nday , nhour
    character(len=256) :: aefilename

    call split_idate(idate,nyear,month,nday,nhour)

    if ( nyear > iyear + 9 ) then
      iyear = nyear/10*10
      write(aefilename,'(a,i0.4,a,i0.4,a)') &
         trim(inpglob)//pthsep//'AERGLOB'//pthsep// &
         trim(scendir(iscen))//pthsep//'aero_1.9x2.5_L26_', &
         iyear, '-', iyear+9, '.nc'
      istatus = nf90_open(aefilename,nf90_nowrite, ncid)
      call checkncerr(istatus,__FILE__,__LINE__,'Error open aeid file')
      istatus = nf90_inq_varid(ncid,'PS',ivarid)
      call checkncerr(istatus,__FILE__,__LINE__,'Error find var PS')
      istatus = nf90_get_var(ncid,ivarid,xps2)
      call checkncerr(istatus,__FILE__,__LINE__,'Error read var PS')
      do is = 1 , naesp
        istatus = nf90_inq_varid(ncid,aespec(is),ivarid)
        call checkncerr(istatus,__FILE__,__LINE__,'Error find var '//aespec(is))
        istatus = nf90_get_var(ncid,ivarid,aev2(:,:,:,:,is))
        call checkncerr(istatus,__FILE__,__LINE__,'Error read var '//aespec(is))
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
                   aeilon,aejlat,iy,jx,aeilev) 
    end do

    do i = 1 , aejlat
      do j = 1 , aeilon
        xps(j,i) = xps2(j,i,m1)*wt1+xps2(j,i,m2)*wt2
      end do
    end do

    paeid_2 = xps*0.01
    p0 = p0*0.01
    r4pt = real(ptop)

    call bilinx2(paeid_3,paeid_2,xlon,xlat,aet42lon,aet42lat, &
                 aeilon,aejlat,iy,jx,1)

    do i = 1 , iy 
      do j = 1 , jx
        do l = 1 , kz
          prcm=((paeid_3(j,i)*0.1-r4pt)*sigma2(l)+r4pt)*10.
          k0 = -1
          do k = aeilev , 1 , -1
            pmpi = paeid_3(j,i)*aet42hybm(k)+aet42hyam(k)*p0
            k0 = k
            if (prcm > pmpi) exit
          end do
          if (k0 == aeilev) then        
            pmpj = paeid_3(j,i)*aet42hybm(aeilev-1)+aet42hyam(aeilev-1)*p0
            pmpi = paeid_3(j,i)*aet42hybm(aeilev)+aet42hyam(aeilev)*p0

            do is = 1 , naesp
              aev4(j,i,l,is) = aev3(j,i,aeilev,is) + &
                 (aev3(j,i,aeilev,is) - aev3(j,i,aeilev-1,is)) * &
                 (prcm-pmpi)/(pmpi-pmpj)
            end do
          else if (k0 >= 1) then
            pmpj = paeid_3(j,i)*aet42hybm(k0)+aet42hyam(k0)*p0
            pmpi = paeid_3(j,i)*aet42hybm(k0+1)+aet42hyam(k0+1)*p0
            do is = 1 , naesp
              aev4(j,i,l,is) = (aev3(j,i,k0+1,is)*(prcm-pmpj) + &
                                aev3(j,i,k0,is)*(prcm-pmpi))/(pmpi-pmpj)
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
      call checkncerr(istatus,__FILE__,__LINE__,'Error close aeid file')
    end if
  end subroutine close_ae_icbc

end module mod_ae_icbc
