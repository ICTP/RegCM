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
  use mod_ch_param

  private

  integer(ik4) :: chilon , chjlat , chilev

  real(rk8) , pointer , dimension(:) :: cht42lon
  real(rk8) , pointer , dimension(:) :: cht42lat
  real(rk8) , pointer , dimension(:) :: cht42hyam , cht42hybm
  !
  ! Oxidant climatology variables
  !
  real(rk8) , pointer , dimension(:) :: xtimes
  type(rcm_time_and_date) , pointer , dimension(:) :: itimes
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
    integer(ik4) :: ivarid , idimid
    integer(ik4) :: nyear , month , nday , nhour
    character(len=256) :: chfilename
    integer(ik4) :: ncid , istatus

    call split_idate(idate,nyear,month,nday,nhour)

    write(chfilename,'(a,i0.2,a)') &
       trim(inpglob)//pthsep//'OXIGLOB'//pthsep// &
       'mz4_19990401.nc'
    istatus = nf90_open(chfilename,nf90_nowrite,ncid)
    call checkncerr(istatus,__FILE__,__LINE__, &
       'Error open file chemical '//trim(chfilename))
    istatus = nf90_inq_dimid(ncid,'lon',idimid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error find dim lon')
    istatus = nf90_inquire_dimension(ncid,idimid,len=chilon)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error inquire dim lon')
    istatus = nf90_inq_dimid(ncid,'lat',idimid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error find dim lat')
    istatus = nf90_inquire_dimension(ncid,idimid,len=chjlat)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error inquire dim lat')
    istatus = nf90_inq_dimid(ncid,'lev',idimid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error find dim lev')
    istatus = nf90_inquire_dimension(ncid,idimid,len=chilev)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error inquire dim lev')

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
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error find var lon')
    istatus = nf90_get_var(ncid,ivarid,cht42lon)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read var lon')
    istatus = nf90_inq_varid(ncid,'lat',ivarid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error find var lat')
    istatus = nf90_get_var(ncid,ivarid,cht42lat)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read var lat')
    istatus = nf90_inq_varid(ncid,'hyam',ivarid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error find var hyam')
    istatus = nf90_get_var(ncid,ivarid,cht42hyam)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read var hyam')
    istatus = nf90_inq_varid(ncid,'hybm',ivarid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error find var hybm')
    istatus = nf90_get_var(ncid,ivarid,cht42hybm)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read var hybm')
    istatus = nf90_inq_varid(ncid,'P0',ivarid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error find var P0')
    istatus = nf90_get_var(ncid,ivarid,p0)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read var P0')
    istatus = nf90_close(ncid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error close file chemical')

    p0 = 1000.0 !p0*0.01
    r4pt = real(ptop)
    write(stdout,*) 'Static read OK.'

    call readps

  end subroutine header_ch_icbc

  subroutine get_ch_icbc(idate)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate
    character(len=256)     :: chfilename
    integer(ik4) :: year1 , month1 , day1 , hour1
    integer(ik4) :: nyear , month , nday , nhour

    call split_idate(globidate1,year1,month1,day1,hour1)
    call split_idate(idate,nyear,month,nday,nhour)
    call find_data(idate,globidate1,chfilename)
    write(*,*)chfilename
    call readmz4(idate,chfilename)

    ! Lumping of Mozart species to CBMZ species

    chv4(:,:,:,cb_O3)    = chv4_1(:,:,:,mz_O3)*w_o3/amd
    chv4(:,:,:,cb_NO)    = chv4_1(:,:,:,mz_NO)*w_no/amd
    chv4(:,:,:,cb_NO2)   = chv4_1(:,:,:,mz_NO2)*w_no2/amd
    chv4(:,:,:,cb_HNO3)  = chv4_1(:,:,:,mz_HNO3)*w_hno3/amd
    chv4(:,:,:,cb_HNO4)  = chv4_1(:,:,:,mz_HO2NO2)*w_hno4/amd
    chv4(:,:,:,cb_N2O5)  = chv4_1(:,:,:,mz_N2O5)*w_n2o5/amd
    chv4(:,:,:,cb_H2O2)  = chv4_1(:,:,:,mz_H2O2)*w_h2o2/amd
    chv4(:,:,:,cb_CH4)   = chv4_1(:,:,:,mz_CH4)*w_ch4/amd
    chv4(:,:,:,cb_CO)    = chv4_1(:,:,:,mz_CO)*w_co/amd
    chv4(:,:,:,cb_SO2)   = chv4_1(:,:,:,mz_SO2)*w_so2/amd
    chv4(:,:,:,cb_H2SO4) = chv4_1(:,:,:,mz_SO4)*w_h2so4/amd
    chv4(:,:,:,cb_DMS)   = chv4_1(:,:,:,mz_DMS)*w_dms/amd
    chv4(:,:,:,cb_PAR)   = (3*chv4_1(:,:,:,mz_C3H8) + &
                            4*chv4_1(:,:,:,mz_BIGALK) + &
                            chv4_1(:,:,:,mz_C3H6) + &
                            chv4_1(:,:,:,mz_BIGENE))*w_par/amd
    chv4(:,:,:,cb_C2H6)  = chv4_1(:,:,:,mz_C2H6)*w_c2h6/amd
    chv4(:,:,:,cb_ETH)   = chv4_1(:,:,:,mz_C2H4)*w_eth/amd
    chv4(:,:,:,cb_OLET)  = chv4_1(:,:,:,mz_BIGENE)*w_olet/amd
    chv4(:,:,:,cb_OLEI)  = chv4_1(:,:,:,mz_BIGENE)*w_olei/amd
    chv4(:,:,:,cb_TOL)   = chv4_1(:,:,:,mz_TOLUENE)*w_tol/amd
    chv4(:,:,:,cb_XYL)   = chv4_1(:,:,:,mz_TOLUENE)*w_xyl/amd
    chv4(:,:,:,cb_ISOP)  = chv4_1(:,:,:,mz_ISOP)*w_isop/amd
    chv4(:,:,:,cb_CRES)  = chv4_1(:,:,:,mz_CRESOL)*w_cres/amd
    chv4(:,:,:,cb_OPEN)  = chv4_1(:,:,:,mz_BIGALD)*w_open/amd
    chv4(:,:,:,cb_ISOPN) = chv4_1(:,:,:,mz_ISOPNO3)*w_isopn/amd
    chv4(:,:,:,cb_ISOPRD)= (chv4_1(:,:,:,mz_MVK) + &
                            chv4_1(:,:,:,mz_MACR)+ &
                            chv4_1(:,:,:,mz_HYDRALD))*w_isoprd/amd
    chv4(:,:,:,cb_ONIT)  = (chv4_1(:,:,:,mz_ONIT) + &
                            chv4_1(:,:,:,mz_ONITR))*w_onit/amd
    chv4(:,:,:,cb_MGLY)  = chv4_1(:,:,:,mz_CH3COCHO)*w_mgly/amd
    chv4(:,:,:,cb_AONE)  = (chv4_1(:,:,:,mz_CH3COCH3) + &
                            chv4_1(:,:,:,mz_HYAC) + &
                            chv4_1(:,:,:,mz_MEK))*w_aone/amd
    chv4(:,:,:,cb_PAN)   = chv4_1(:,:,:,mz_PAN)*w_pan/amd
    chv4(:,:,:,cb_CH3OOH)= chv4_1(:,:,:,mz_CH3OOH)*w_ch3ooh/amd
    chv4(:,:,:,cb_ETHOOH)= chv4_1(:,:,:,mz_C2H5OOH)*w_ethooh/amd
    chv4(:,:,:,cb_ALD2)  = (chv4_1(:,:,:,mz_CH3CHO) + &
                            chv4_1(:,:,:,mz_GLYALD))*w_ald2/amd
    chv4(:,:,:,cb_HCHO)  = chv4_1(:,:,:,mz_CH2O)*w_hcho/amd
    chv4(:,:,:,cb_CH3OH) = chv4_1(:,:,:,mz_CH3OH)*w_ch3oh/amd

    call write_ch_icbc(idate)
  end subroutine get_ch_icbc

  subroutine readps
    implicit none
    character(len=256) :: chfilename
    integer :: ncid , istatus , ivarid

    write(chfilename,'(a,i0.2,a)') &
       trim(inpglob)//pthsep//'OXIGLOB'//pthsep// &
       'mz4_19990401.nc'
    istatus = nf90_open(chfilename,nf90_nowrite,ncid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error open file chemical')
    istatus = nf90_inq_varid(ncid,'PS',ivarid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error find var PS')
    istatus = nf90_get_var(ncid,ivarid,xps)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read var PS')
    xps = xps*0.01
    call bilinx2(pchem_3,xps,xlon,xlat,cht42lon,cht42lat, &
                 chilon,chjlat,jx,iy)
  end subroutine readps

  subroutine find_data(idate,idate0,chfilename)
    implicit none
    integer, parameter                   ::nfile=127
    type(rcm_time_and_date) , intent(in) :: idate
    type(rcm_time_and_date) , intent(in) :: idate0
    real(rk8) , dimension(nfile,124),save :: timearray

    character(len=256),intent(out) :: chfilename
    character(len=44),dimension(nfile),save:: filename
    integer(ik4) :: i , recc , it
    integer(ik4) :: timid , timlen
    character(len=32) :: chfilemm,chfilemmm1,chfilemmp1
    character(len=64) :: cunit , ccal
    character(len=32) :: datename,datenamem1,datenamep1
    character(len=4)  :: yyyy
    character(len=2)  :: mm
    character(len=7)  :: yyyy_mm,yyyy_mmm1,yyyy_mmp1
    integer(ik4) :: year , month , day , hour
    integer(ik4) :: year1 , month1 , day1 , hour1
    integer(ik4) :: nyear , nmonth , nday , nhour
    integer :: ncid , istatus
    type (rcm_time_interval) :: tdif

    call split_idate(idate,year,month,day,hour)
    datename   = tochar(idate)
    datenamem1 = tochar(prevmon(idate))
    datenamep1 = tochar(nextmon(idate))
    yyyy     = datename(1:4)
    mm       = datename(6:7)
    yyyy_mm  = datename(1:7)
    yyyy_mmm1  = datenamem1(1:7)
    yyyy_mmp1  = datenamep1(1:7)
    chfilemm='MZ4-synoz-NCEPT42.mz4.h0.'//yyyy_mm
    chfilemmm1='MZ4-synoz-NCEPT42.mz4.h0.'//yyyy_mmm1
    chfilemmp1='MZ4-synoz-NCEPT42.mz4.h0.'//yyyy_mmp1
    open(10,file=trim(inpglob)//pthsep//'OXIGLOB'//pthsep//'list')

    call split_idate(idate0,year1,month1,day1,hour1)
    call split_idate(idate,nyear,nmonth,nday,nhour)
    if ( nyear .eq. year1 .and. nmonth .eq. month1 .and. nhour .eq. hour1 ) then
      do i = 1 , nfile
        read(10,*) filename(i)
      end do
      recc=0
      do i = 1 , nfile
        chfilename=trim(trim(inpglob)//pthsep//'OXIGLOB'//pthsep//filename(i))
        istatus = nf90_open(chfilename,nf90_nowrite,ncid)
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error open file chemical')
        istatus = nf90_inq_dimid(ncid,'time',timid)
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error find dim time')

        istatus = nf90_inquire_dimension(ncid,timid,len=timlen)
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error inquire time')

        istatus = nf90_inq_varid(ncid,'time',timid)
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error find var time')

        istatus = nf90_get_att(ncid,timid,'units',cunit)
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error read time units')
        ccal = 'gregorian'
        call getmem1d(xtimes,1,timlen,'mod_ein:xtimes')
        istatus = nf90_get_var(ncid,timid,xtimes)
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error read time')
        cunit="days since 0000-01-01 00:00:00"
        do it = 1 , timlen
          timearray(i,it) = xtimes(it)
        end do
        istatus = nf90_close(ncid)
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error close file')
      end do
    end if
    ccal = 'gregorian'
    cunit = "days since 1950-01-01 00:00:00"
    timlen = 124
    call getmem1d(itimes,1,timlen,'mod_ein:itimes')
    do i = 1 , nfile
      do it = 1 , 124
        itimes(it)=timeval2date(dble(timearray(i,it)),cunit,ccal)
        tdif = itimes(it)-idate
        if(tohours(tdif) == 0) then
          chfilename=trim(trim(inpglob)//pthsep//'OXIGLOB'//pthsep//filename(i))
        end if
      end do
    end do
    rewind(10)
  end subroutine find_data

  subroutine readmz4(idate,chfilename)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate
    character(len=256),intent(in) :: chfilename
    integer(ik4) :: i , is , j , k , l , k0,recc
    integer(ik4) :: timid , timlen
    character(len=64) :: cunit , ccal
    real(rk8) :: wt1 , wt2
    integer(ik4) :: ncid , istatus , ivarid
    integer(ik4) , dimension(4) :: icount , istart
    integer(ik4) :: year , month , day , hour
    integer(ik4) :: year1 , month1 , day1 , hour1
    integer(ik4) :: it
    type(rcm_time_interval) :: tdif

    istatus = nf90_open(chfilename,nf90_nowrite,ncid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error open file chemical')
    istatus = nf90_inq_dimid(ncid,'time',timid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error find dim time')

    istatus = nf90_inquire_dimension(ncid,timid,len=timlen)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error inquire time')

    istatus = nf90_inq_varid(ncid,'time',timid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error find var time')

    istatus = nf90_get_att(ncid,timid,'units',cunit)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read time units')

    ccal = 'gregorian'
    call getmem1d(itimes,1,timlen,'mod_ein:itimes')
    call getmem1d(xtimes,1,timlen,'mod_ein:xtimes')
    istatus = nf90_get_var(ncid,timid,xtimes)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read time')
    cunit="days since 1950-01-01 00:00:00"
    do it = 1 , timlen
      itimes(it) = timeval2date(dble(xtimes(it)),cunit,ccal)

      tdif = itimes(it)-idate
      if ( tohours(tdif) == 0 ) then
        recc = it
      end if
    end do
    write(*,*)chfilename,'  ',recc,'  ',tochar(itimes(recc))

    do k = 1 , 4
      istart(k) = 1
    end do

    icount(1) = chilon
    icount(2) = chjlat
    icount(3) = chilev
    icount(4) = 1
    call split_idate(idate,year,month,day,hour)
    call split_idate(globidate1,year1,month1,day1,hour1)

    do is = 1 , nchsp
      istart(4) = recc
      istatus = nf90_inq_varid(ncid,trim(chspec(is))//'_VMR_inst',ivarid)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find var '//trim(chspec(is)))
      istatus = nf90_get_var(ncid,ivarid,xinp(:,:,:,is),istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read var '//trim(chspec(is)))
      call bilinx2(chv3(:,:,:,is),xinp(:,:,:,is),xlon,xlat,cht42lon,cht42lat, &
                   chilon,chjlat,jx,iy,chilev)
    end do
    do i = 1 , iy
      do j = 1 , jx
        do l = 1 , kz
          prcm=((pchem_3(j,i)*0.1-r4pt)*sigmah(l)+r4pt)*10.
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
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error close file chemical')
  end subroutine readmz4

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
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
