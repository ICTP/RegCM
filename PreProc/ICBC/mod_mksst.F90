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

module mod_mksst

  use mod_intkinds
  use mod_realkinds
  use mod_stdio
  use mod_grid
  use mod_memutil
  use mod_constants
  use mod_dynparam
  use mod_message
  use mod_nchelper
  use mod_interp
  use netcdf

  private

  logical :: lopen , lhasice
  integer(ik4) :: ncst , ntime
  integer(ik4) , dimension(3) :: ivar
  type(rcm_time_and_date) , dimension(:) , pointer :: itime
  real(rk8) , dimension(:,:) , pointer :: work1 , work2
  real(rk8) , dimension(:,:) , pointer :: work3 , work4
  real(rk8) , dimension(:) , pointer :: xtime

  data lopen/.false./
  character(len=256) :: sstfile

  public :: readsst , closesst

  contains
!
!-----------------------------------------------------------------------
!
  subroutine readsst(tsccm, idate)
    implicit none
    real(rk8) , dimension(jx,iy) , intent(inout) :: tsccm
    type(rcm_time_and_date) , intent(in) :: idate
    integer(ik4) :: istatus , idimid , itvar
    integer(ik4) , dimension(3) :: istart , icount
    character(len=64) :: timeunits , timecal
    integer(ik4) :: i , j , irec
    integer(ik4) :: iyy , im , id , ih
    type(rcm_time_interval) :: ks1 , ks2
    real(rk8) :: wt , a , b

    call split_idate(idate,iyy,im,id,ih)

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
      do i = 1 , ntime
        itime(i) = timeval2date(xtime(i), timeunits, timecal)
      end do
      lopen = .true.
    end if

    if (idate > itime(ntime) .or. idate < itime(1)) then
      write (stderr,*) 'Cannot find ', tochar(idate), ' in SST file'
      write (stderr,*) 'Range is : ', tochar(itime(1)) , ' to ', &
                       tochar(itime(ntime))
      call die('readsst')
    end if

    irec = 1
    do i = 1 , ntime
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
      do i = 1 , iy
        do j = 1 , jx
          if ( (landuse(j,i) > 13.9 .and. landuse(j,i) < 15.1) ) then
            if ( work1(j,i) > -900.0 ) then
              tsccm(j,i) = work1(j,i)
              if (lhasice) then
                if ( work3(j,i) > -900.0) then
                  if ( work3(j,i) > 35. ) then
                    tsccm(j,i) = 271.0
                  end if
                end if
              end if
            else
              ! If latitude is greater than 70 and month is in winter,
              ! assume ice coverage if missing
              if ( xlat(j,i) > 55.0 .and. ( im > 9 .or. im < 4 ) ) then
                tsccm(j,i) = 271.0
              else if ( xlat(j,i) < -55.0 .and. ( im > 5 .and. im < 9 ) ) then
                tsccm(j,i) = 271.0
              else
                ! Find nearest water points
                a = nearn(j,i,work1)
                if ( a > -900.0 ) tsccm(j,i) = a
              end if
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
      wt = dble(tohours(ks1)/tohours(ks2))
      do i = 1 , iy
        do j = 1 , jx
          if ( (landuse(j,i) > 13.9 .and. landuse(j,i) < 15.1) ) then
            if ( (work1(j,i) > -900.0 .and. work2(j,i) > -900.0) ) then
              tsccm(j,i) = wt*work1(j,i) + (1.0-wt)*work2(j,i)
              if (lhasice) then
                if ( work3(j,i) > -900.0 .and. work4(j,i) > -900.0 ) then
                  if ( wt*work3(j,i)+(1-wt)*work4(j,i) > 35. ) then
                    tsccm(j,i) = 271.0
                  endif
                end if
              end if
            else
              ! If latitude is greater than 70 and month is in winter,
              ! assume ice coverage if missing
              if ( xlat(j,i) > 55.0 .and. ( im > 9 .or. im < 4 ) ) then
                tsccm(j,i) = 271.0
              else if ( xlat(j,i) < -55.0 .and. ( im > 5 .and. im < 9 ) ) then
                tsccm(j,i) = 271.0
              else
                a = nearn(j,i,work1)
                b = nearn(j,i,work2)
                ! Find nearest water points
                if ( a > -900 .and. b > -900) tsccm(j,i) = a*wt + (1.0-wt)*b
              end if
            end if
          end if
        end do
      end do
    end if
  end subroutine readsst

  real(rk8) function nearn(jp,ip,sst)
    implicit none
    integer(ik4) , intent(in) :: jp , ip
    real(rk8) , dimension(:,:) , intent(in) :: sst
    real(rk8) :: wt , wtsum , distsig
    integer(ik4) :: i , j , nr , np
    if ( all(sst < -900) ) then
      nearn = -999.0
      return
    end if
    nr = 1
    np = -1
    nearn = 0.0D0
    wtsum = 0.0D0
    do while ( np < 0 )
      do i = ip - nr , ip + nr
        do j = jp - nr , jp + nr
          if ( j == jp .and. i == ip ) cycle
          if ( i < 1 .or. i > iy ) cycle
          if ( j < 1 .or. j > jx ) cycle
          if ( sst(j,i) > -900.0 ) then
            wt = 1.0D0/sqrt(dble(i-ip)**2+dble(j-jp)**2)
            distsig = sign(1.0D0,xlat(jp,ip)-xlat(j,i))
            if ( (xlat(jp,ip) - xlat(j,i)) > 2.0D0*epsilon(d_zero) ) then
              nearn = nearn + (max(icetemp,sst(j,i) - lrate * &
               (max(0.0D0,topogm(jp,ip)-topogm(j,i))) - distsig * &
               gcdist(xlat(jp,ip),xlon(jp,ip), &
                      xlat(j,i),xlon(jp,ip))/100000.0D0))*wt
            else
              nearn = nearn + max(icetemp,sst(j,i) - lrate * &
               (max(0.0D0,topogm(jp,ip)-topogm(j,i))))*wt
            end if
            wtsum = wtsum + wt
            np = np + 1
          end if
        end do
      end do
      nr = nr + 1
    end do
    nearn = (nearn / wtsum)
  end function nearn

  subroutine closesst
    implicit none
    integer(ik4) :: istatus
    istatus = nf90_close(ncst)
  end subroutine closesst
!
end module mod_mksst
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
