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

module mod_ch_icbc_clim

  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_stdio
  use mod_dynparam
  use mod_grid
  use mod_wrtoxd
  use mod_kdinterp
  use mod_date
  use mod_memutil
  use mod_message
  use mod_nchelper
  use mod_ch_param
  use netcdf

  private

  integer(ik4) :: chilon, chjlat, chilev

  real(rkx), pointer, contiguous, dimension(:) :: cht42lon
  real(rkx), pointer, contiguous, dimension(:) :: cht42lat
  real(rkx), pointer, contiguous, dimension(:) :: cht42hyam, cht42hybm
  !
  ! Oxidant climatology variables
  !
  real(rkx) :: p0
  real(rkx), pointer, contiguous, dimension(:,:) :: pchem_3
  real(rkx), pointer, contiguous, dimension(:,:) :: xps31, xps32
  real(rkx), pointer, contiguous, dimension(:,:,:,:) :: chv31, chv32
  real(rkx), pointer, contiguous, dimension(:,:) :: xps
  real(rkx), pointer, contiguous, dimension(:,:,:) :: xinp
  real(rkx), pointer, contiguous, dimension(:,:,:,:) :: chv4_1
  real(rkx), pointer, contiguous, dimension(:,:,:,:) :: chv4_2
  real(rkx), pointer, contiguous, dimension(:,:,:,:) :: chv4_3

  real(rkx) :: prcm, pmpi, pmpj
  real(rkx) :: xpt
  integer(ik4) :: ism = -1
  type (rcm_time_and_date), save :: iref1, iref2

  integer(ik4) :: ncicbc, ivarps, irec
  public :: init_ch_icbc_clim, get_ch_icbc_clim, close_ch_icbc_clim

  type(rcm_time_and_date) :: iodate
  type(h_interpolator) :: hint

  data ncicbc /-1/

  contains

  subroutine init_ch_icbc_clim(idate)
    implicit none
    type(rcm_time_and_date), intent(in) :: idate
    type(rcm_time_interval) :: tdif
    integer(ik4) :: ivarid, idimid
    character(len=256) :: chfilename, icbcfilename
    integer(ik4) :: ncid, istatus

    iodate = monfirst(idate)
    write (icbcfilename,'(a,a,a,a,a,a)') trim(dirglob), pthsep, &
           trim(domname), '_ICBC.', trim(tochar10(iodate)), '.nc'
    write(stdout,*) 'Opening ',trim(icbcfilename)
    istatus = nf90_open(icbcfilename,nf90_nowrite, ncicbc)
    if ( istatus /= nf90_noerr ) then
      iodate = idate
      write (icbcfilename,'(a,a,a,a,a,a)') trim(dirglob), pthsep, &
             trim(domname), '_ICBC.', trim(tochar10(iodate)), '.nc'
      write(stdout,*) 'Opening ',trim(icbcfilename)
      istatus = nf90_open(icbcfilename,nf90_nowrite, ncicbc)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error open ICBC file '//trim(icbcfilename))
    end if
    istatus = nf90_inq_varid(ncicbc,'ps',ivarps)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error find var ps in icbc file '//trim(icbcfilename))
    tdif = (idate-iodate)
    irec = int(tohours(tdif)/ibdyfrq)+1

    write(chfilename,'(a,i0.2,a)') &
       trim(inpglob)//pthsep//'OXIGLOB'//pthsep//'mz4_19990401.nc'
    istatus = nf90_open(chfilename,nf90_nowrite,ncid)
    call checkncerr(istatus,__FILE__,__LINE__, &
       'Error open file '//trim(chfilename))
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
    call getmem2d(xps31,1,jx,1,iy,'mod_ch_icbc:xps31')
    call getmem2d(xps32,1,jx,1,iy,'mod_ch_icbc:xps32')
    call getmem4d(chv31,1,jx,1,iy,1,chilev,1,nchsp,'mod_ch_icbc:chv31')
    call getmem4d(chv32,1,jx,1,iy,1,chilev,1,nchsp,'mod_ch_icbc:chv32')
    call getmem4d(chv4_1,1,jx,1,iy,1,kz,1,nchsp,'mod_ch_icbc:chv4_1')
    call getmem4d(chv4_2,1,jx,1,iy,1,kz,1,nchsp,'mod_ch_icbc:chv4_2')
    call getmem4d(chv4_3,1,jx,1,iy,1,kz,1,nchsp,'mod_ch_icbc:chv4_2')

    call getmem1d(cht42lon,1,chilon,'mod_ch_icbc:cht42lon')
    call getmem1d(cht42lat,1,chjlat,'mod_ch_icbc:cht42lat')
    call getmem1d(cht42hyam,1,chilev,'mod_ch_icbc:cht42hyam')
    call getmem1d(cht42hybm,1,chilev,'mod_ch_icbc:cht42hybm')
    call getmem2d(xps,1,chilon,1,chjlat,'mod_ch_icbc:xps')
    call getmem3d(xinp,1,chilon,1,chjlat,1,chilev,'mod_ch_icbc:xinp')

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

    call h_interpolator_create(hint,cht42lat,cht42lon,xlat,xlon)

    xpt = real(ptop)
    write(stdout,*) 'Static read OK.'
  end subroutine init_ch_icbc_clim

  subroutine get_ch_icbc_clim(idate)
    implicit none
    type(rcm_time_and_date), intent(in) :: idate
    integer(ik4) :: nyear, month, nday, nhour
    logical :: doread, lfuture
    character(len=256) :: icbcfilename
    integer(ik4), dimension(3) :: istart, icount
    type(rcm_time_and_date) :: imonmidd
    type(rcm_time_interval) :: tdif
    real(rk8) :: xfac1, xfac2, odist
    integer(ik4) :: im1, im2, istatus

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
    if ( nyear > 2020 ) then
      lfuture = .true.
    else
      lfuture = .false.
    end if

    if ( .not. lsamemonth(idate, iodate) ) then
      istatus = nf90_close(ncicbc)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error close ICBC file')
      write (icbcfilename,'(a,a,a,a,a,a)') trim(dirglob), pthsep, &
             trim(domname), '_ICBC.', trim(tochar10(idate)), '.nc'
      write(stdout,*) 'Opening ',trim(icbcfilename)
      istatus = nf90_open(icbcfilename,nf90_nowrite, ncicbc)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error open ICBC file '//trim(icbcfilename))
      istatus = nf90_inq_varid(ncicbc,'ps',ivarps)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find var ps in icbc file '//trim(icbcfilename))
      iodate = idate
      irec = 1
    end if

    istart(1) = 1
    istart(2) = 1
    istart(3) = irec
    icount(1) = jx
    icount(2) = iy
    icount(3) = 1
    istatus = nf90_get_var(ncicbc,ivarps,pchem_3,istart,icount)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read var ps')
    irec = irec + 1

    call read2m(im1,im2,doread,lfuture)

    tdif = idate-iref1
    xfac1 = tohours(tdif)
    tdif = iref2-idate
    xfac2 = tohours(tdif)
    odist = xfac1 + xfac2
    xfac1 = xfac1/odist
    xfac2 = 1.0_rk8-xfac1
    ! rq: pppv(mozart) to pppm
    chv4_3 = (chv4_1*real(xfac2,rkx)+chv4_2*real(xfac1,rkx))
    chv4(:,:,:,cb_O3)    = chv4_3(:,:,:,mz_O3)*w_o3/amd
    chv4(:,:,:,cb_NO)    = chv4_3(:,:,:,mz_NO)*w_no/amd
    chv4(:,:,:,cb_NO2)   = chv4_3(:,:,:,mz_NO2)*w_no2/amd
    chv4(:,:,:,cb_HNO3)  = chv4_3(:,:,:,mz_HNO3)*w_hno3/amd
    chv4(:,:,:,cb_HNO4)  = chv4_3(:,:,:,mz_HO2NO2)*w_hno4/amd
    chv4(:,:,:,cb_N2O5)  = chv4_3(:,:,:,mz_N2O5)*w_n2o5/amd
    chv4(:,:,:,cb_H2O2)  = chv4_3(:,:,:,mz_H2O2)*w_h2o2/amd
    chv4(:,:,:,cb_CH4)   = chv4_3(:,:,:,mz_CH4)*w_ch4/amd
    chv4(:,:,:,cb_CO)    = chv4_3(:,:,:,mz_CO)*w_co/amd
    chv4(:,:,:,cb_SO2)   = chv4_3(:,:,:,mz_SO2)*w_so2/amd
    chv4(:,:,:,cb_H2SO4) = chv4_3(:,:,:,mz_SO4)*w_h2so4/amd
    chv4(:,:,:,cb_DMS)   = chv4_3(:,:,:,mz_DMS)*w_dms/amd
    chv4(:,:,:,cb_PAR)   = (3*chv4_3(:,:,:,mz_C3H8) + &
                            4*chv4_3(:,:,:,mz_BIGALK) + &
                              chv4_3(:,:,:,mz_C3H6) + &
                              chv4_3(:,:,:,mz_BIGENE))*w_par/amd
    chv4(:,:,:,cb_C2H6)  = chv4_3(:,:,:,mz_C2H6)*w_c2h6/amd
    chv4(:,:,:,cb_ETH)   = chv4_3(:,:,:,mz_C2H4)*w_eth/amd
    chv4(:,:,:,cb_OLET)  = chv4_3(:,:,:,mz_BIGENE)*w_olet/amd
    chv4(:,:,:,cb_OLEI)  = chv4_3(:,:,:,mz_BIGENE)*w_olei/amd
    chv4(:,:,:,cb_TOL)   = chv4_3(:,:,:,mz_TOLUENE)*w_tol/amd
    chv4(:,:,:,cb_XYL)   = chv4_3(:,:,:,mz_TOLUENE)*w_xyl/amd
    chv4(:,:,:,cb_ISOP)  = chv4_3(:,:,:,mz_ISOP)*w_isop/amd
    chv4(:,:,:,cb_CRES)  = chv4_3(:,:,:,mz_CRESOL)*w_cres/amd
    chv4(:,:,:,cb_OPEN)  = chv4_3(:,:,:,mz_BIGALD)*w_open/amd
    chv4(:,:,:,cb_ISOPN) = chv4_3(:,:,:,mz_ISOPNO3)*w_isopn/amd
    chv4(:,:,:,cb_ISOPRD)= (chv4_3(:,:,:,mz_MVK) + &
                            chv4_3(:,:,:,mz_MACR) + &
                            chv4_3(:,:,:,mz_HYDRALD))*w_isoprd/amd
    chv4(:,:,:,cb_ONIT)  = (chv4_3(:,:,:,mz_ONIT) + &
                            chv4_3(:,:,:,mz_ONITR))*w_onit/amd
    chv4(:,:,:,cb_MGLY)  = chv4_3(:,:,:,mz_CH3COCHO)*w_mgly/amd
    chv4(:,:,:,cb_AONE)  = (chv4_3(:,:,:,mz_CH3COCH3) + &
                            chv4_3(:,:,:,mz_HYAC) + &
                            chv4_3(:,:,:,mz_MEK))*w_aone/amd
    chv4(:,:,:,cb_PAN)   = chv4_3(:,:,:,mz_PAN)*w_pan/amd
    chv4(:,:,:,cb_CH3OOH)= chv4_3(:,:,:,mz_CH3OOH)*w_ch3ooh/amd
    chv4(:,:,:,cb_ETHOOH)= chv4_3(:,:,:,mz_C2H5OOH)*w_ethooh/amd
    chv4(:,:,:,cb_ALD2)  = (chv4_3(:,:,:,mz_CH3CHO) + &
                            chv4_3(:,:,:,mz_GLYALD))*w_ald2/amd
    chv4(:,:,:,cb_HCHO)  = chv4_3(:,:,:,mz_CH2O)*w_hcho/amd
    chv4(:,:,:,cb_CH3OH) = chv4_3(:,:,:,mz_CH3OH)*w_ch3oh/amd

    call write_ch_icbc(idate)

  end subroutine get_ch_icbc_clim

  subroutine read2m(im1,im2,doread,lfuture)
    implicit none
    integer(ik4), intent(in) :: im1, im2
    logical, intent(in) :: doread, lfuture
    integer(ik4) :: i, is, j, k, l, k0
    character(len=256) :: chfilename
    real(rkx) :: wt1, wt2
    integer(ik4) :: ncid, istatus, ivarid
    real(rkx), parameter :: rglrog = rgas*lrate*regrav

    if ( doread ) then
      if ( lfuture ) then
        write(chfilename,'(a,i0.2,a)') &
           trim(inpglob)//pthsep//'OXIGLOB'//pthsep// &
           'mz4-EMAC_avg_2040-2050_',im1,'.nc'
      else
        write(chfilename,'(a,i0.2,a)') &
           trim(inpglob)//pthsep//'OXIGLOB'//pthsep// &
           'mz4_avg_1999-2009_',im1,'.nc'
      end if
      istatus = nf90_open(chfilename,nf90_nowrite,ncid)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error open file '//trim(chfilename))
      write(stdout, *) 'Opening ', trim(chfilename)
      istatus = nf90_inq_varid(ncid,'PS',ivarid)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find var PS')
      istatus = nf90_get_var(ncid,ivarid,xps)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read var PS')
      call h_interpolate_cont(hint,xps,xps31)
      do is = 1, nchsp
        istatus = nf90_inq_varid(ncid,trim(chspec(is))//'_VMR_inst',ivarid)
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error find var '//trim(chspec(is)))
        istatus = nf90_get_var(ncid,ivarid,xinp)
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error read var '//trim(chspec(is)))
        where ( xinp < d_zero ) xinp = d_zero
        call h_interpolate_cont(hint,xinp,chv31(:,:,:,is))
      end do
      istatus = nf90_close(ncid)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error close file chemical')
    end if
    do i = 1, iy
      do j = 1, jx
        do l = 1, kz
          prcm = ((pchem_3(j,i)*0.1_rkx-xpt)*sigmah(l)+xpt)*1000.0_rkx
          k0 = -1
          do k = chilev, 1, -1
            pmpi = cht42hyam(k)*p0+xps31(j,i)*cht42hybm(k)
            k0 = k
            if (prcm > pmpi) exit
          end do
          if ( prcm < pmpi ) then
            chv4_1(j,i,l,is) = chv31(j,i,1,is)
          else if (k0 == chilev) then
            do is = 1, nchsp
              chv4_1(j,i,l,is) = 0.5_rkx * &
                (chv31(j,i,chilev,is) + chv31(j,i,chilev-1,is)) * &
                  exp(rglrog*log(prcm/pmpi))
            end do
          else if (k0 >= 1) then
            pmpj = cht42hyam(k0  )*p0+xps31(j,i)*cht42hybm(k0  )
            pmpi = cht42hyam(k0+1)*p0+xps31(j,i)*cht42hybm(k0+1)
            wt1 = log(prcm/pmpj)/log(pmpi/pmpj)
            wt2 = 1.0 - wt1
            do is = 1, nchsp
              chv4_1(j,i,l,is) = chv31(j,i,k0+1,is)*wt1 + chv31(j,i,k0,is)*wt2
            end do
          end if
        end do
      end do
    end do
    if ( doread ) then
      if ( lfuture ) then
        write(chfilename,'(a,i0.2,a)') &
           trim(inpglob)//pthsep//'OXIGLOB'//pthsep// &
           'mz4-EMAC_avg_2040-2050_',im2,'.nc'
      else
        write(chfilename,'(a,i0.2,a)') &
           trim(inpglob)//pthsep//'OXIGLOB'//pthsep// &
           'mz4_avg_1999-2009_',im2,'.nc'
      end if
      istatus = nf90_open(chfilename,nf90_nowrite,ncid)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error open file '//trim(chfilename))
      write(stdout, *) 'Opening ', trim(chfilename)
      istatus = nf90_inq_varid(ncid,'PS',ivarid)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find var PS')
      istatus = nf90_get_var(ncid,ivarid,xps)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read var PS')
      call h_interpolate_cont(hint,xps,xps32)
      do is = 1, nchsp
        istatus = nf90_inq_varid(ncid,trim(chspec(is))//'_VMR_inst',ivarid)
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error find var '//trim(chspec(is)))
        istatus = nf90_get_var(ncid,ivarid,xinp)
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error read var '//trim(chspec(is)))
        where ( xinp < d_zero ) xinp = d_zero
        call h_interpolate_cont(hint,xinp,chv32(:,:,:,is))
      end do
      istatus = nf90_close(ncid)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error close file chemical')
    end if
    do i = 1, iy
      do j = 1, jx
        do l = 1, kz
          prcm = ((pchem_3(j,i)*0.1_rkx-xpt)*sigmah(l)+xpt)*1000.0_rkx
          k0 = -1
          do k = chilev, 1, -1
            pmpi = cht42hyam(k)*p0+xps32(j,i)*cht42hybm(k)
            k0 = k
            if (prcm > pmpi) exit
          end do
          if ( prcm < pmpi ) then
            chv4_2(j,i,l,is) = chv32(j,i,1,is)
          else if (k0 == chilev) then
            pmpj = cht42hyam(chilev-1)*p0+xps32(j,i)*cht42hybm(chilev-1)
            pmpi = cht42hyam(chilev  )*p0+xps32(j,i)*cht42hybm(chilev  )
            do is = 1, nchsp
              chv4_2(j,i,l,is) = 0.5_rkx * &
                (chv32(j,i,chilev,is) + chv32(j,i,chilev-1,is)) * &
                  exp(rglrog*log(prcm/pmpi))
            end do
          else if (k0 >= 1) then
            pmpj = cht42hyam(k0  )*p0+xps32(j,i)*cht42hybm(k0  )
            pmpi = cht42hyam(k0+1)*p0+xps32(j,i)*cht42hybm(k0+1)
            wt1 = log(prcm/pmpj)/log(pmpi/pmpj)
            wt2 = 1.0 - wt1
            do is = 1, nchsp
              chv4_2(j,i,l,is) = chv32(j,i,k0+1,is)*wt1 + chv32(j,i,k0,is)*wt2
            end do
          end if
        end do
      end do
    end do
  end subroutine read2m

  subroutine close_ch_icbc_clim
    implicit none
    integer(ik4) :: istatus
    call h_interpolator_destroy(hint)
    if ( ncicbc > 0 ) then
      istatus = nf90_close(ncicbc)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error close icbc file')
    end if
    return
  end subroutine close_ch_icbc_clim

  integer(ik4) function inextmon(im)
    implicit none
    integer(ik4), intent(in) :: im
    inextmon = im+1
    if ( inextmon == 13 ) inextmon = 1
  end function inextmon

  integer(ik4) function iprevmon(im)
    implicit none
    integer(ik4), intent(in) :: im
    iprevmon = im-1
    if ( iprevmon == 0 ) iprevmon = 12
  end function iprevmon

end module mod_ch_icbc_clim
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
