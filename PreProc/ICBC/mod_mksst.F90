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

module mod_mksst

  use mod_intkinds
  use mod_realkinds
  use mod_stdio
  use mod_grid
  use mod_stdio
  use mod_date
  use mod_memutil
  use mod_constants
  use mod_dynparam
  use mod_message
  use mod_nchelper
  use mod_earth
  use netcdf

  private

  logical :: lopen, lhasice
  integer(ik4) :: ncst, ntime
  integer(ik4), dimension(3) :: ivar
  type(rcm_time_and_date), dimension(:), contiguous, pointer :: itime
  real(rkx), dimension(:,:), contiguous, pointer :: work1, work2
  real(rkx), dimension(:,:), contiguous, pointer :: work3, work4
  real(rkx), dimension(:), contiguous, pointer :: xtime
  real(rkx) :: icet

  data lopen/.false./
  character(len=256) :: sstfile

  public :: readsst, closesst

  contains

  subroutine readsst(tsccm, idate)
    implicit none
    real(rkx), dimension(jx,iy), intent(inout) :: tsccm
    type(rcm_time_and_date), intent(in) :: idate
    integer(ik4) :: istatus, idimid, itvar
    integer(ik4), dimension(3) :: istart, icount
    character(len=64) :: timeunits, timecal
    integer(ik4) :: i, j, irec
    integer(ik4) :: iyy, im, id, ih
    type(rcm_time_interval) :: ks1, ks2
    real(rkx) :: wt, a, b

    call split_idate(idate,iyy,im,id,ih)

    select case (ssttyp)
      case ('EIN15','EIN75','EIXXX')
        icet = 271.46_rkx
      case default
        icet = 271.35_rkx
    end select
    if (.not. lopen) then
      sstfile=trim(dirglob)//pthsep//trim(domname)//'_SST.nc'
      istatus = nf90_open(sstfile, nf90_nowrite, ncst)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error Opening '//trim(sstfile))
      istatus = nf90_inq_dimid(ncst, 'time', idimid)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error time dimension '//trim(sstfile))
      istatus = nf90_inquire_dimension(ncst, idimid, len=ntime)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error time dimension '//trim(sstfile))
      if ( ntime == 0 ) then
        write (stderr,*) 'No timesteps for SST data in SST file '// &
                    trim(sstfile)
        call die('readsst')
      end if
      istatus = nf90_inq_varid(ncst, "time", itvar)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error time var '//trim(sstfile))
      istatus = nf90_inq_varid(ncst, "sst", ivar(2))
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error sst var '//trim(sstfile))
      lhasice = .true.
      istatus = nf90_inq_varid(ncst, "ice", ivar(3))
      if ( istatus /= nf90_noerr) then
        lhasice = .false.
      end if
      istatus = nf90_get_att(ncst, itvar, "units", timeunits)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error time var units '//trim(sstfile))
      istatus = nf90_get_att(ncst, itvar, "calendar", timecal)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error time var units '//trim(sstfile))

      call getmem2d(work1,1,jx,1,iy,'mod_mksst:work1')
      call getmem2d(work2,1,jx,1,iy,'mod_mksst:work2')
      if (lhasice) then
        call getmem2d(work3,1,jx,1,iy,'mod_mksst:work3')
        call getmem2d(work4,1,jx,1,iy,'mod_mksst:work4')
      end if
      call getmem1d(xtime,1,ntime,'mod_mksst:xtime')
      call getmem1d(itime,1,ntime,'mod_mksst:itime')

      istatus = nf90_get_var(ncst, itvar, xtime)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error time var read '//trim(sstfile))
      do i = 1, ntime
        itime(i) = timeval2date(xtime(i), timeunits, timecal)
      end do
      lopen = .true.
    end if

    if (idate > itime(ntime) .or. idate < itime(1)) then
      write (stderr,*) 'Cannot find ', tochar(idate), ' in SST file'
      write (stderr,*) 'Range is : ', tochar(itime(1)), ' to ', &
                       tochar(itime(ntime))
      call die('readsst')
    end if

    irec = 1
    do i = 1, ntime
      if ( idate < itime(i) ) then
        exit
      end if
      irec = irec + 1
    end do
    irec = irec - 1

    istart(3) = irec
    istart(2) = 1
    istart(1) = 1
    icount(3) = 1
    icount(2) = iy
    icount(1) = jx
    istatus = nf90_get_var(ncst, ivar(2), work1, istart, icount)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error sst var read '//trim(sstfile))
    if (lhasice) then
      istatus = nf90_get_var(ncst, ivar(3), work3, istart, icount)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error ice var read '//trim(sstfile))
    end if
    if (idate == itime(irec)) then
      do i = 1, iy
        do j = 1, jx
          if ( (landuse(j,i) > 14.9 .and. landuse(j,i) < 15.1) ) then
            if ( work1(j,i) > -900.0_rkx ) then
              tsccm(j,i) = work1(j,i)
              if (lhasice) then
                if ( work3(j,i) > -900.0_rkx) then
                  if ( work3(j,i) > 35._rkx ) then
                    tsccm(j,i) = 271.0_rkx
                  end if
                end if
              end if
            else
              ! Find nearest water points
              a = nearn(j,i,work1)
              if ( a > -900.0_rkx ) tsccm(j,i) = a
            end if
          end if
        end do
      end do
    else
      istart(3) = irec+1
      istart(2) = 1
      istart(1) = 1
      icount(3) = 1
      icount(2) = iy
      icount(1) = jx
      istatus = nf90_get_var(ncst, ivar(2), work2, istart, icount)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error sst var read '//trim(sstfile))
      if (lhasice) then
        istatus = nf90_get_var(ncst, ivar(3), work4, istart, icount)
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error ice var read '//trim(sstfile))
      end if
      ks1 = itime(irec+1)-idate
      ks2 = itime(irec+1)-itime(irec)
      wt = real(tohours(ks1)/tohours(ks2),rkx)
      do i = 1, iy
        do j = 1, jx
          if ( (landuse(j,i) > 14.9_rkx .and. landuse(j,i) < 15.1_rkx) ) then
            if ( (work1(j,i) > -900.0_rkx .and. work2(j,i) > -900.0_rkx) ) then
              tsccm(j,i) = wt*work1(j,i) + (1.0_rkx-wt)*work2(j,i)
              if (lhasice) then
                if ( work3(j,i) > -900.0_rkx .and. &
                     work4(j,i) > -900.0_rkx ) then
                  if ( wt*work3(j,i)+(1-wt)*work4(j,i) > 35.0_rkx ) then
                    tsccm(j,i) = 271.0_rkx
                  endif
                end if
              end if
            else
              a = nearn(j,i,work1)
              b = nearn(j,i,work2)
              ! Find nearest water points
              if ( a > -900_rkx .and. b > -900_rkx) then
                tsccm(j,i) = a*wt + (1.0_rkx-wt)*b
              end if
            end if
          end if
        end do
      end do
    end if
  end subroutine readsst

  real(rkx) function nearn(jp,ip,sst)
    implicit none
    integer(ik4), intent(in) :: jp, ip
    real(rkx), dimension(:,:), intent(in) :: sst
    real(rkx) :: wt, wtsum
    ! real(rkx) :: distsig
    integer(ik4) :: i, j, nr, np, maxn
    if ( all(sst < -900_rkx) ) then
      nearn = -999.0_rkx
      return
    end if
    maxn = min(size(sst,1),size(sst,2))/4
    nr = 1
    np = -1
    nearn = 0.0_rkx
    wtsum = 0.0_rkx
    do while ( np < 0 .and. nr < maxn )
      do i = ip - nr, ip + nr
        do j = jp - nr, jp + nr
          if ( j == jp .and. i == ip ) cycle
          if ( i < 1 .or. i > iy ) cycle
          if ( j < 1 .or. j > jx ) cycle
          if ( sst(j,i) > -900.0_rkx ) then
            wt = 1.0_rkx/sqrt(real(i-ip,rkx)**2+real(j-jp,rkx)**2)
            !distsig = sign(1.0_rkx,xlat(jp,ip)-xlat(j,i))
            !if ( (xlat(jp,ip) - xlat(j,i)) > 2.0_rkx*epsilon(d_zero) ) then
            !  nearn = nearn + (max(icet,sst(j,i) - lrate * &
            !   (max(0.0_rkx,topogm(jp,ip)-topogm(j,i))) - distsig * &
            !   gcdist_simple(xlat(jp,ip),xlon(jp,ip), &
            !                 xlat(j,i),xlon(jp,ip))/100.0_rkx))*wt
            !else
            !  nearn = nearn + max(icet,sst(j,i) - lrate * &
            !   (max(0.0_rkx,topogm(jp,ip)-topogm(j,i))))*wt
            !  nearn = nearn + max(icet,sst(j,i))*wt
            !end if
            nearn = nearn + max(icet,sst(j,i))*wt
            wtsum = wtsum + wt
            np = np + 1
          end if
        end do
      end do
      nr = nr + 1
    end do
    if ( wtsum > 0.0_rkx ) then
      nearn = (nearn / wtsum)
    else
      nearn = -999.0_rkx
    end if
  end function nearn

  subroutine closesst
    implicit none
    integer(ik4) :: istatus
    istatus = nf90_close(ncst)
  end subroutine closesst

end module mod_mksst
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
