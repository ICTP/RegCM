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

module mod_ch_icbc

  use mod_intkinds
  use mod_realkinds
  use mod_stdio
  use mod_dynparam
  use mod_grid
  use mod_wrtoxd
  use mod_interp
  use mod_date
  use mod_memutil
  use mod_message
  use mod_nchelper
  use netcdf

  private
!
  integer(ik4) :: chilon , chjlat , chilev

  real(rk8) , pointer , dimension(:) :: cht42lon
  real(rk8) , pointer , dimension(:) :: cht42lat
  real(rk8) , pointer , dimension(:) :: cht42hyam , cht42hybm
!
! Oxidant climatology variables
!
  real(rk8) :: p0
  real(rk8) , pointer , dimension(:,:) :: pchem_3
  real(rk8) , pointer , dimension(:,:,:,:) :: chv3
  real(rk8) , pointer , dimension(:,:) :: xps
  real(rk8) , pointer , dimension(:,:,:,:) :: xinp
  real(rk8) , pointer , dimension(:,:,:,:) :: chv4_1
  real(rk8) , pointer , dimension(:,:,:,:) :: chv4_2

  real(rk8) :: prcm , pmpi , pmpj
  real(rk8) :: r4pt
  integer(ik4) :: ism
  type (rcm_time_and_date) , save :: iref1 , iref2
  
  public :: header_ch_icbc , get_ch_icbc , close_ch_icbc

  contains

  subroutine header_ch_icbc(idate)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate
    type (rcm_time_and_date) :: imonmidd
    integer(ik4) :: ivarid , idimid
    integer(ik4) :: nyear , month , nday , nhour
    character(len=256) :: chfilename
    integer(ik4) :: ncid , istatus
    integer(ik4) :: im1 , im2

    call split_idate(idate,nyear,month,nday,nhour)
    imonmidd = monmiddle(idate)
    im1 = month
    im2 = month
    if ( idate >= imonmidd ) then
      im2 = inextmon(im2)
      iref1 = imonmidd
      iref2 = monmiddle(nextmon(idate))
    else
      im1 = iprevmon(im1)
      iref1 = monmiddle(prevmon(idate))
      iref2 = imonmidd
    end if
    ism = im1

    write(chfilename,'(a,i0.2,a)') &
       trim(inpglob)//pthsep//'OXIGLOB'//pthsep// &
       'mz4_avg_2000-2007_',im1,'.nc'
    istatus = nf90_open(chfilename,nf90_nowrite,ncid)
    call checkncerr(istatus,__FILE__,__LINE__, &
       'Error open file chemical '//trim(chfilename))
    istatus = nf90_inq_dimid(ncid,'lon',idimid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error find dim lon')
    istatus = nf90_inquire_dimension(ncid,idimid,len=chilon)
    call checkncerr(istatus,__FILE__,__LINE__,'Error inquire dim lon')
    istatus = nf90_inq_dimid(ncid,'lat',idimid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error find dim lat')
    istatus = nf90_inquire_dimension(ncid,idimid,len=chjlat)
    call checkncerr(istatus,__FILE__,__LINE__,'Error inquire dim lat')
    istatus = nf90_inq_dimid(ncid,'lev',idimid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error find dim lev')
    istatus = nf90_inquire_dimension(ncid,idimid,len=chilev)
    call checkncerr(istatus,__FILE__,__LINE__,'Error inquire dim lev')

    call getmem2d(pchem_3,1,jx,1,iy,'mod_ch_icbc:pchem_3_1')
    call getmem4d(chv3,1,jx,1,iy,1,chilev,1,nchsp,'mod_ch_icbc:chv3')
    call getmem4d(chv4_1,1,jx,1,iy,1,kz,1,nchsp,'mod_ch_icbc:chv4_1')
    call getmem4d(chv4_2,1,jx,1,iy,1,kz,1,nchsp,'mod_ch_icbc:chv4_2')

    call getmem1d(cht42lon,1,chilon,'mod_ch_icbc:cht42lon')
    call getmem1d(cht42lat,1,chjlat,'mod_ch_icbc:cht42lat')
    call getmem1d(cht42hyam,1,chilev,'mod_ch_icbc:cht42hyam')
    call getmem1d(cht42hybm,1,chilev,'mod_ch_icbc:cht42hybm')
    call getmem2d(xps,1,chilon,1,chjlat,'mod_ch_icbc:xps1')
    call getmem4d(xinp,1,chilon,1,chjlat,1,chilev,1,nchsp,'mod_ch_icbc:xinp')

    istatus = nf90_inq_varid(ncid,'lon',ivarid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error find var lon')
    istatus = nf90_get_var(ncid,ivarid,cht42lon)
    call checkncerr(istatus,__FILE__,__LINE__,'Error read var lon')
    istatus = nf90_inq_varid(ncid,'lat',ivarid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error find var lat')
    istatus = nf90_get_var(ncid,ivarid,cht42lat)
    call checkncerr(istatus,__FILE__,__LINE__,'Error read var lat')
    istatus = nf90_inq_varid(ncid,'hyam',ivarid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error find var hyam')
    istatus = nf90_get_var(ncid,ivarid,cht42hyam)
    call checkncerr(istatus,__FILE__,__LINE__,'Error read var hyam')
    istatus = nf90_inq_varid(ncid,'hybm',ivarid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error find var hybm')
    istatus = nf90_get_var(ncid,ivarid,cht42hybm)
    call checkncerr(istatus,__FILE__,__LINE__,'Error read var hybm')
    istatus = nf90_inq_varid(ncid,'P0',ivarid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error find var P0')
    istatus = nf90_get_var(ncid,ivarid,p0)
    call checkncerr(istatus,__FILE__,__LINE__,'Error read var P0')
    istatus = nf90_close(ncid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error close file chemical')

    p0 = p0*0.01
    r4pt = real(ptop)
    write(stdout,*) 'Static read OK.'

    call read2m(im1,im2)

  end subroutine header_ch_icbc

  subroutine get_ch_icbc(idate)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate
    integer(ik4) :: nyear , month , nday , nhour
    logical :: doread
    type (rcm_time_and_date) :: imonmidd
    type (rcm_time_interval) :: tdif
    real(rk8) :: xfac1 , xfac2 , odist
    integer(ik4) :: im1 , im2, n

    call split_idate(idate,nyear,month,nday,nhour)
    imonmidd = monmiddle(idate)
    im1 = month
    im2 = month
    if ( idate >= imonmidd ) then
      im2 = inextmon(im2)
      iref1 = imonmidd
      iref2 = monmiddle(nextmon(idate))
    else
      im1 = iprevmon(im1)
      iref1 = monmiddle(prevmon(idate))
      iref2 = imonmidd
    end if
    doread = .false.
    if ( ism /= im1 ) then
      ism = im1
      doread = .true.
    end if

    if ( doread ) then
      call read2m(im1,im2)
    end if

    tdif = idate-iref1
    xfac1 = tohours(tdif)
    tdif = iref2-idate
    xfac2 = tohours(tdif)
    odist = xfac1 + xfac2
    xfac1 = xfac1/odist
    xfac2 = d_one-xfac1
! rq: pppv(mozart) to pppm 
    chv4 = (chv4_1*xfac2+chv4_2*xfac1)

   do n=1, nchsp
    chv4(:,:,:,n)  = chv4(:,:,:,n) * mw(n)/amd  
   end do

    call write_ch_icbc(idate)

  end subroutine get_ch_icbc

  subroutine read2m(im1,im2)
    implicit none
    integer(ik4) , intent(in) :: im1 , im2
    integer(ik4) :: i , is , j , k , l , k0
    character(len=256) :: chfilename
    real(rk8) :: wt1 , wt2
    integer(ik4) :: ncid , istatus , ivarid

    write(chfilename,'(a,i0.2,a)') &
       trim(inpglob)//pthsep//'OXIGLOB'//pthsep// &
       'mz4_avg_2000-2007_',im1,'.nc'
    istatus = nf90_open(chfilename,nf90_nowrite,ncid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error open file chemical')
    write(stdout, *) trim(chfilename)
    istatus = nf90_inq_varid(ncid,'PS',ivarid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error find var PS')
    istatus = nf90_get_var(ncid,ivarid,xps)
    call checkncerr(istatus,__FILE__,__LINE__,'Error read var PS')
    xps = xps*0.01
    do is = 1 , nchsp
      istatus = nf90_inq_varid(ncid,trim(chspec(is))//'_VMR_avrg',ivarid)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find var '//trim(chspec(is)))
      istatus = nf90_get_var(ncid,ivarid,xinp(:,:,:,is))
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read var '//trim(chspec(is)))
      call bilinx2(chv3(:,:,:,is),xinp(:,:,:,is),xlon,xlat,cht42lon,cht42lat, &
                   chilon,chjlat,iy,jx,chilev) 
    end do
    call bilinx2(pchem_3,xps,xlon,xlat,cht42lon,cht42lat, &
                 chilon,chjlat,iy,jx,1)
    do i = 1 , iy 
      do j = 1 , jx
        do l = 1 , kz
          prcm=((pchem_3(j,i)*0.1-r4pt)*sigma2(l)+r4pt)*10.
          k0 = -1
          do k = chilev , 1 , -1
            pmpi = cht42hyam(k)*p0+pchem_3(j,i)*cht42hybm(k)
            k0 = k
            if (prcm > pmpi) exit
          end do
          if (k0 == chilev) then
            pmpj = cht42hyam(chilev-1)*p0+pchem_3(j,i)*cht42hybm(chilev-1)
            pmpi = cht42hyam(chilev  )*p0+pchem_3(j,i)*cht42hybm(chilev  )
            do is = 1 , nchsp
              chv4_1(j,i,l,is) = chv3(j,i,chilev,is) + &
                 (chv3(j,i,chilev-1,is) - chv3(j,i,chilev,is)) * &
                 (prcm-pmpi)/(pmpi-pmpj)
            end do
          else if (k0 >= 1) then
            pmpj = cht42hyam(k0  )*p0+pchem_3(j,i)*cht42hybm(k0  )
            pmpi = cht42hyam(k0+1)*p0+pchem_3(j,i)*cht42hybm(k0+1)
            wt1 = (prcm-pmpj)/(pmpi-pmpj)
            wt2 = 1.0 - wt1
            do is = 1 , nchsp
              chv4_1(j,i,l,is) = chv3(j,i,k0+1,is)*wt1 + chv3(j,i,k0,is)*wt2
            end do
          end if
        end do
      end do
    end do
    istatus = nf90_close(ncid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error close file chemical')
    write(chfilename,'(a,i0.2,a)') &
       trim(inpglob)//pthsep//'OXIGLOB'//pthsep// &
       'mz4_avg_2000-2007_',im2,'.nc'
    istatus = nf90_open(chfilename,nf90_nowrite,ncid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error open file chemical')
    write(stdout, *) trim(chfilename)
    istatus = nf90_inq_varid(ncid,'PS',ivarid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error find var PS')
    istatus = nf90_get_var(ncid,ivarid,xps)
    call checkncerr(istatus,__FILE__,__LINE__,'Error read var PS')
    xps = xps*0.01
    do is = 1 , nchsp
      istatus = nf90_inq_varid(ncid,trim(chspec(is))//'_VMR_avrg',ivarid)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find var '//trim(chspec(is)))
      istatus = nf90_get_var(ncid,ivarid,xinp(:,:,:,is))
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read var '//trim(chspec(is)))
      call bilinx2(chv3(:,:,:,is),xinp(:,:,:,is),xlon,xlat,cht42lon,cht42lat, &
                   chilon,chjlat,iy,jx,chilev) 
    end do
    call bilinx2(pchem_3,xps,xlon,xlat,cht42lon,cht42lat, &
                 chilon,chjlat,iy,jx,1)
    do i = 1 , iy 
      do j = 1 , jx
        do l = 1 , kz
          prcm=((pchem_3(j,i)*0.1-r4pt)*sigma2(l)+r4pt)*10.
          k0 = -1
          do k = chilev , 1 , -1
            pmpi = cht42hyam(k)*p0+pchem_3(j,i)*cht42hybm(k)
            k0 = k
            if (prcm > pmpi) exit
          end do
          if (k0 == chilev) then
            pmpj = cht42hyam(chilev-1)*p0+pchem_3(j,i)*cht42hybm(chilev-1)
            pmpi = cht42hyam(chilev  )*p0+pchem_3(j,i)*cht42hybm(chilev  )
            do is = 1 , nchsp
              chv4_2(j,i,l,is) = chv3(j,i,chilev,is) + &
                 (chv3(j,i,chilev-1,is) - chv3(j,i,chilev,is)) * &
                 (prcm-pmpi)/(pmpi-pmpj)
            end do
          else if (k0 >= 1) then
            pmpj = cht42hyam(k0  )*p0+pchem_3(j,i)*cht42hybm(k0  )
            pmpi = cht42hyam(k0+1)*p0+pchem_3(j,i)*cht42hybm(k0+1)
            wt1 = (prcm-pmpj)/(pmpi-pmpj)
            wt2 = 1.0 - wt1
            do is = 1 , nchsp
              chv4_2(j,i,l,is) = chv3(j,i,k0+1,is)*wt1 + chv3(j,i,k0,is)*wt2
            end do
          end if
        end do
      end do
    end do
    istatus = nf90_close(ncid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error close file chemical')
  end subroutine read2m

  subroutine close_ch_icbc
    implicit none
    return
  end subroutine close_ch_icbc

  integer(ik4) function inextmon(im)
    implicit none
    integer(ik4) , intent(in) :: im
    inextmon = im+1
    if ( inextmon == 13 ) inextmon = 1
  end function inextmon

  integer(ik4) function iprevmon(im)
    implicit none
    integer(ik4) , intent(in) :: im
    iprevmon = im-1
    if ( iprevmon == 0 ) iprevmon = 12
  end function iprevmon

end module mod_ch_icbc
