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
!
module mod_che_ncio
!
  use m_realkinds
  use mod_dynparam
  use mod_message
!
  private
!
  public :: read_texture , read_aerosol
!
  integer :: istatus

  contains

  subroutine read_texture(nats,texture)
    use netcdf
    implicit none

    integer , intent(in) :: nats
    real(dp) , dimension(iy,jx,nats) , intent(out) :: texture

    integer :: ivarid
    integer :: i , j , n
    integer :: idmin
    integer , dimension(3) :: istart , icount
    character(256) :: dname
    real(sp), dimension(jx,iy) ::  toto

    dname = trim(dirter)//pthsep//trim(domname)//'_DOMAIN000.nc'
    istatus = nf90_open(dname, nf90_nowrite, idmin)
    call check_ok(__FILE__,__LINE__, &
                  'Error Opening Domain file '//trim(dname),'DOMAIN FILE OPEN')
    istatus = nf90_inq_varid(idmin, 'texture_fraction', ivarid)
    call check_ok(__FILE__,__LINE__,'Variable texture_fraction miss', &
                  'DOMAIN FILE')
    istart(2) = 1
    istart(1) = 1
    icount(3) = 1
    icount(2) = iy
    icount(1) = jx
    do n = 1 , nats
      istart(3) = n
      istatus = nf90_get_var(idmin, ivarid, toto, istart, icount)
      call check_ok(__FILE__,__LINE__,'Variable texture_frac read error', &
                    'DOMAIN FILE')
      do j = 1 , jx
        do i = 1 , iy
          texture(i,j,n) = dble(toto(j,i))*0.01D0
          if (texture(i,j,n)<d_zero) texture(i,j,n)=d_zero
        end do
      end do
    end do
    istatus = nf90_close(idmin)
    call check_ok(__FILE__,__LINE__,'Domain file close error','DOMAIN FILE')
  end subroutine read_texture

  subroutine read_aerosol(chtrname,chemsrc)
    use netcdf
    implicit none

    character(256) :: aername
    character(5) , dimension(ntr) , intent(in) :: chtrname
    real(dp) , dimension(iy,jx,12,ntr) , intent(out) :: chemsrc

    integer :: ncid , ivarid
    real(sp) , dimension(jx,iy) :: toto
    character(5) :: aerctl
    integer , dimension(3) :: istart , icount
    integer :: itr , i , j , m

    aername = trim(dirglob)//pthsep//trim(domname)//'_AERO.nc'
    istatus = nf90_open(aername, nf90_nowrite, ncid)
    call check_ok(__FILE__,__LINE__, &
         'Error Opening Aerosol file '//trim(aername),'AEROSOL FILE OPEN')

    do itr = 1 , ntr
      aerctl = chtrname(itr)
      write (aline, *) itr , aerctl
      call say(myid)
      if ( aerctl(1:4) /= 'DUST') then
        if ( aerctl(1:3) == 'SO2' ) then
          if ( aertyp(4:4) == '1' ) then
            istatus = nf90_inq_varid(ncid, 'so2', ivarid)
            call check_ok(__FILE__,__LINE__, &
                          'Variable so2 miss','AEROSOL FILE')
            istatus = nf90_get_var(ncid, ivarid, toto)
            call check_ok(__FILE__,__LINE__, &
                          'Variable so2 read error','AEROSOL FILE')
            do m = 1 , 12
              do j = 1 , jx
                do i = 1 , iy
                  chemsrc(i,j,m,itr) = dble(toto(j,i))
                end do
              end do
            end do
          end if
          if ( aertyp(5:5) == '1' ) then
            istatus = nf90_inq_varid(ncid, 'so2_monthly', ivarid)
            call check_ok(__FILE__,__LINE__, &
                          'Variable so2_mon miss','AEROSOL FILE')
            istart(1) = 1
            istart(2) = 1
            icount(1) = jx
            icount(2) = iy
            icount(3) = 1
            do m = 1 , 12
              istart(3) = m
              istatus = nf90_get_var(ncid,ivarid,toto,istart,icount)
              call check_ok(__FILE__,__LINE__, &
                            'Variable so2_mon read err','AEROSOL FILE')
              do j = 1 , jx
                do i = 1 , iy
                  chemsrc(i,j,m,itr) = chemsrc(i,j,m,itr) + dble(toto(j,i))
                end do
              end do
            end do
          end if
        else if ( aerctl(1:2) == 'BC' ) then
          if ( aertyp(4:4) == '1' ) then
            istatus = nf90_inq_varid(ncid, 'bc', ivarid)
            call check_ok(__FILE__,__LINE__, &
                          'Variable bc miss','AEROSOL FILE')
            istatus = nf90_get_var(ncid, ivarid, toto)
            call check_ok(__FILE__,__LINE__, &
                          'Variable bc read error','AEROSOL FILE')
            do m = 1 , 12
              do j = 1 , jx
                do i = 1 , iy
                  chemsrc(i,j,m,itr) = dble(toto(j,i))
                end do
              end do
            end do
          end if
          if ( aertyp(5:5) == '1' ) then
            istatus = nf90_inq_varid(ncid, 'bc_monthly', ivarid)
            call check_ok(__FILE__,__LINE__, &
                          'Variable bc_mon miss','AEROSOL FILE')
            istart(1) = 1
            istart(2) = 1
            icount(1) = jx
            icount(2) = iy
            icount(3) = 1
            do m = 1 , 12
              istart(3) = m
              istatus = nf90_get_var(ncid,ivarid,toto,istart,icount)
              call check_ok(__FILE__,__LINE__, &
                            'Variable bc_mon read err','AEROSOL FILE')
              do j = 1 , jx
                do i = 1 , iy
                  chemsrc(i,j,m,itr) = chemsrc(i,j,m,itr) + dble(toto(j,i))
                end do
              end do
            end do
          end if
        else if ( aerctl(1:2) == 'OC' ) then
          if ( aertyp(4:4) == '1' ) then
            istatus = nf90_inq_varid(ncid, 'oc', ivarid)
            call check_ok(__FILE__,__LINE__, &
                          'Variable oc miss','AEROSOL FILE')
            istatus = nf90_get_var(ncid, ivarid, toto)
            call check_ok(__FILE__,__LINE__, &
                          'Variable oc read error','AEROSOL FILE')
            do m = 1 , 12
              do j = 1 , jx
                do i = 1 , iy
                  chemsrc(i,j,m,itr) = dble(toto(j,i))
                end do
              end do
            end do
          end if
          if ( aertyp(5:5) == '1' ) then
            istatus = nf90_inq_varid(ncid, 'oc_monthly', ivarid)
            call check_ok(__FILE__,__LINE__, &
                          'Variable oc_mon miss','AEROSOL FILE')
            istart(1) = 1
            istart(2) = 1
            icount(1) = jx
            icount(2) = iy
            icount(3) = 1
            do m = 1 , 12
              istart(3) = m
              istatus = nf90_get_var(ncid,ivarid,toto,istart,icount)
              call check_ok(__FILE__,__LINE__, &
                            'Variable oc_mon read err','AEROSOL FILE')
              do j = 1 , jx
                do i = 1 , iy
                  chemsrc(i,j,m,itr) = chemsrc(i,j,m,itr) + dble(toto(j,i))
                end do
              end do
            end do
          end if
        end if
      end if
    end do

    istatus = nf90_close(ncid)
    call check_ok(__FILE__,__LINE__, &
                  'Error Close Aerosol file '//trim(aername), &
                  'AEROSOL FILE CLOSE')

  end subroutine read_aerosol

  subroutine check_ok(f,l,m1,mf)
    use netcdf
    implicit none
    character(*) , intent(in) :: f, m1 , mf
    integer , intent(in) :: l
    if (istatus /= nf90_noerr) then
      write (6,*) trim(m1)
      write (6,*) nf90_strerror(istatus)
      call fatal(f,l,trim(mf))
    end if
  end subroutine check_ok

end module mod_che_ncio
