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
  use mod_realkinds
  use mod_dynparam
  use mod_mpmessage
  use mod_che_indices
  use netcdf
!
  private
!
  public :: read_texture , read_aerosol , read_emission
!
  integer :: istatus
  integer :: recc

  contains

  subroutine read_texture(nats,texture)
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
      call say
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

  subroutine read_emission(chtrname,lmonth,echemsrc)
    implicit none

    integer , intent(in) :: lmonth
    character(5) , dimension(ntr) , intent(in) :: chtrname
    real(dp) , dimension(iy,jx,12,ntr) , intent(inout) :: echemsrc

    character(256) :: aername
    character(5) :: aerctl
    integer :: ncid , istatus
    integer , dimension(3) :: istart , icount
    integer :: itr , ivarid
    integer :: currmonth , oldmonth
    save    :: oldmonth

    aername = trim(dirglob)//pthsep//trim(domname)//'_AERO.nc'
    istatus = nf90_open(aername, nf90_nowrite, ncid)
    call check_ok(__FILE__,__LINE__, &
                  'Error Opening Aerosol file '//trim(aername), &
                  'AEROSOL FILE OPEN ERROR')

    !*** intialized in start_chem
    !*** Advice record counter
    recc = recc + 1
    do itr = 1 , ntr
      aerctl = chtrname(itr)
      write (*, *) itr, aerctl, lmonth, recc
    end do

!   if ( oldmonth /= lmonth ) then
!     recc = recc+1
!   else 
!     recc = recc
!   end if

    istart(1) = 1
    istart(2) = 1
    istart(3) = recc
    icount(1) = jx
    icount(2) = iy
    icount(3) = 1

    ! NO emission                  
    if ( ino /= 0 ) then
      call rvar(ncid,istart,icount,ino,lmonth,echemsrc, &
                'a_NO',.false.,'bio_nox')
    end if
    ! CO emission
    if ( ico /= 0 ) then
      call rvar(ncid,istart,icount,ico,lmonth,echemsrc, &
                'a_CO',.false.,'bio_co','o_co')
    end if
    ! HCHO emission                  
    if ( ihcho /= 0 ) then
      call rvar(ncid,istart,icount,ihcho,lmonth,echemsrc,'a_HCHO',.false.)
    end if
    ! ACET emission                  
    if ( iacet /= 0 ) then
      call rvar(ncid,istart,icount,iacet,lmonth,echemsrc, &
                'a_ACET',.false.,'bio_acet')
    end if
    ! SO2 emission
    if ( iso2 /= 0 ) then
      call rvar(ncid,istart,icount,iso2,lmonth,echemsrc, &
                'a_SO2',.false.,'b_so2')
    end if
    ! CH4
    if ( ich4 /= 0 ) then
      call rvar(ncid,istart,icount,ich4,lmonth,echemsrc, &
                'a_ch4',.false.,'bio_ch4')
    end if
    ! Ethane
    if ( ic2h6 /= 0 ) then
      call rvar(ncid,istart,icount,ic2h6,lmonth,echemsrc, &
                'a_ETHANE',.false.,'bio_c2h6','o_c2h6')
    end if
    ! PAR
    if ( ipar /= 0 ) then
      call rvar(ncid,istart,icount,ipar,lmonth,echemsrc, &
                'a_c3h8',.false.,'a_butane','bio_c3h8','o_c3h8')
    end if
    ! Ethene
    if ( iethe /= 0 ) then
      call rvar(ncid,istart,icount,iethe,lmonth,echemsrc, &
                'a_ETHENE',.false.,'bio_c2h4','o_c2h4')
    end if
    ! Termenal Alkene
    if ( iolt /= 0 ) then
      call rvar(ncid,istart,icount,iolt,lmonth,echemsrc, &
                'a_PROPENE',.false.,'bio_c3h6')
    end if
    ! Internal Alkene
    if ( ioli /= 0 ) then
      call rvar(ncid,istart,icount,ioli,lmonth,echemsrc,'a_BIGENE',.true.)
    end if
    ! Isoprene
    if ( iisop /= 0 ) then
      call rvar(ncid,istart,icount,iisop,lmonth,echemsrc,'bio_isop',.false.)
    end if
    ! Toluene
    if ( itolue /= 0 ) then
      istatus = nf90_inq_varid(ncid, 'a_XYLENE', ivarid)
      if ( istatus == nf90_noerr ) then
        call rvar(ncid,istart,icount,itolue,lmonth,echemsrc, &
                  'a_TOLUENE',.true.,'b_TOLUENE')
      end if
    end if
    ! Xylene
    if ( ixyl /= 0 ) then
      istatus = nf90_inq_varid(ncid, 'a_TOLUENE', ivarid)
      if ( istatus == nf90_noerr ) then
        call rvar(ncid,istart,icount,ixyl,lmonth,echemsrc, &
                  'a_XYLENE',.true.)
      else
        call rvar(ncid,istart,icount,ixyl,lmonth,echemsrc, &
                  'a_TOLUENE',.true.,'b_XYLENE')
      end if
    end if
    ! Acetaldehyde
    if ( iald2 /= 0 ) then
      call rvar(ncid,istart,icount,iald2,lmonth,echemsrc,'a_ALD2',.false.)
    end if
    ! Methanol + Ethanol
    if ( imoh /= 0 ) then
      call rvar(ncid,istart,icount,imoh,lmonth,echemsrc, &
                'a_MOH',.false.,'b_ch3oh','bio_ch3oh')
    end if           
    ! DMS
    if ( idms /= 0 ) then
      call rvar(ncid,istart,icount,idms,lmonth,echemsrc,'o_DMS',.false.)
    end if
    ! OC and BC anthropogenic + biomass burning
    if ( ibchb /= 0 ) then
      call rvar(ncid,istart,icount,ibchb,lmonth,echemsrc, &
                'b_BC',.false.,'a_BC')
    end if
    if ( iochb /= 0 ) then
      call rvar(ncid,istart,icount,iochb,lmonth,echemsrc, &
                'b_OC',.false.,'a_OC')
    end if

    oldmonth = lmonth

    istatus = nf90_close(ncid)
    call check_ok(__FILE__,__LINE__, &
                  'Error Closing Aerosol file '//trim(aername), &
                  'AEROSOL FILE CLOSE ERROR')

  end subroutine read_emission

  subroutine rvar(ncid,istart,icount,ind,lmonth,echemsrc,cna,lh,cnb,cnc,cnd)
    integer , intent(in) :: ncid
    integer , dimension(3) , intent(in) :: istart , icount
    integer , intent(in) :: lmonth
    real(dp) , dimension(iy,jx,12,ntr) , intent(inout) :: echemsrc
    logical , intent(in) :: lh
    character(len=*) , intent(in) :: cna
    character(len=*) , intent(in) , optional :: cnb
    character(len=*) , intent(in) , optional :: cnc
    character(len=*) , intent(in) , optional :: cnd
    integer :: ivarid , istatus
    real(sp) , dimension(jx,iy) :: toto
    integer :: i , j

    istatus = nf90_inq_varid(ncid, cna, ivarid)     
    call check_ok(__FILE__,__LINE__, &
                  'Variable '//cna//' miss','AEROSOL FILE')
    istatus = nf90_get_var(ncid,ivarid,toto)
    call check_ok(__FILE__,__LINE__, &
                  'Variable '//cna//' read err','AEROSOL FILE')
    if ( lh ) then  ! half of lumped Aromatics
      do j = 1 , jx
        do i = 1 , iy
          echemsrc(i,j,lmonth,ind) = d_half*toto(j,i)
        end do
      end do
    else
      do j = 1 , jx
        do i = 1 , iy
          echemsrc(i,j,lmonth,ind) = toto(j,i)
        end do
      end do
    end if
    if ( present(cnb) ) then
      istatus = nf90_inq_varid(ncid, cnb, ivarid)
      call check_ok(__FILE__,__LINE__, &
                    'Variable '//cnb//' miss','AEROSOL FILE')
      istatus = nf90_get_var(ncid,ivarid,toto,istart,icount)
      call check_ok(__FILE__,__LINE__, &
                    'Variable '//cnb//' read err','AEROSOL FILE')
      do j = 1 , jx
        do i = 1 , iy
          echemsrc(i,j,lmonth,ind) = toto(j,i) + echemsrc(i,j,lmonth,ind)
        end do
      end do
    end if
    if ( present(cnc) ) then
      istatus = nf90_inq_varid(ncid, cnc, ivarid)
      call check_ok(__FILE__,__LINE__, &
                    'Variable '//cnc//' miss','AEROSOL FILE')
      istatus = nf90_get_var(ncid,ivarid,toto,istart,icount)
      call check_ok(__FILE__,__LINE__, &
                    'Variable '//cnc//' read err','AEROSOL FILE')
      do j = 1 , jx
        do i = 1 , iy
          echemsrc(i,j,lmonth,ind) = toto(j,i) + echemsrc(i,j,lmonth,ind)
        end do
      end do
    end if
    if ( present(cnd) ) then
      istatus = nf90_inq_varid(ncid, cnd, ivarid)
      call check_ok(__FILE__,__LINE__, &
                    'Variable '//cnd//' miss','AEROSOL FILE')
      istatus = nf90_get_var(ncid,ivarid,toto,istart,icount)
      call check_ok(__FILE__,__LINE__, &
                    'Variable '//cnd//' read err','AEROSOL FILE')
      do j = 1 , jx
        do i = 1 , iy
          echemsrc(i,j,lmonth,ind) = toto(j,i) + echemsrc(i,j,lmonth,ind)
        end do
      end do
    end if
  end subroutine rvar

  subroutine check_ok(f,l,m1,mf)
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
