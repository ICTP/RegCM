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
!
! Names of the physics schemes in the namelist (physicsparam) and the
! integer codes used inside the model:
!
!   pbl_scheme       ibltyp   none 0, holtslag 1, uw 2, gfs 3, myj 4,
!                             shin-hong 5
!   cumulus_land,    icup_lnd none 0, shallow -1, grell 2, emanuel 4,
!   cumulus_ocean    icup_ocn tiedtke 5, kain-fritsch 6
!   microphysics     ipptls   none 0, noto 2, wsm5 3, wsm7 4, wdm7 5
!   cloud_fraction   icldfrac sundqvist 0, xu-randall 1, thompson 2,
!                             gultepe-isaac 3, teixeira 4, tompkins 5,
!                             echam5 6, digital 7
!   radiation_scheme irrtm    ccm 0, rrtmg 1
!
! The names are case insensitive, with '-' or '_'. The codes of the
! removed schemes (Kuo and Betts-Miller cumulus, SUBEX microphysics, in
! Main/historical) are not reused.
!
module mod_scheme_names

  use mod_intkinds, only : ik4

  implicit none

  private

  ! Returned for a name that is not a scheme of the kind
  integer(ik4), parameter, public :: unknown_scheme = -9999

  integer(ik4), parameter :: maxlen = 16

  type scheme_entry
    character(len=maxlen) :: name
    integer(ik4) :: code
  end type scheme_entry

  type(scheme_entry), dimension(6), parameter :: pbl = [ &
    scheme_entry('none',0), scheme_entry('holtslag',1), &
    scheme_entry('uw',2), scheme_entry('gfs',3), &
    scheme_entry('myj',4), scheme_entry('shin-hong',5) ]
  type(scheme_entry), dimension(6), parameter :: cumulus = [ &
    scheme_entry('none',0), scheme_entry('shallow',-1), &
    scheme_entry('grell',2), scheme_entry('emanuel',4), &
    scheme_entry('tiedtke',5), scheme_entry('kain-fritsch',6) ]
  type(scheme_entry), dimension(5), parameter :: micro = [ &
    scheme_entry('none',0), scheme_entry('noto',2), &
    scheme_entry('wsm5',3), scheme_entry('wsm7',4), &
    scheme_entry('wdm7',5) ]
  type(scheme_entry), dimension(8), parameter :: cloud = [ &
    scheme_entry('sundqvist',0), scheme_entry('xu-randall',1), &
    scheme_entry('thompson',2), scheme_entry('gultepe-isaac',3), &
    scheme_entry('teixeira',4), scheme_entry('tompkins',5), &
    scheme_entry('echam5',6), scheme_entry('digital',7) ]
  type(scheme_entry), dimension(2), parameter :: rad = [ &
    scheme_entry('ccm',0), scheme_entry('rrtmg',1) ]

  public :: scheme_code, scheme_name, scheme_list

  contains
  !
  ! Code of the scheme name of a kind (pbl, cumulus, microphysics,
  ! cloud_fraction, radiation), unknown_scheme if not a scheme of the kind;
  ! msg explains why.
  !
  integer(ik4) function scheme_code(kind,name,msg) result(code)
    implicit none
    character(len=*), intent(in) :: kind, name
    character(len=*), intent(out) :: msg
    character(len=maxlen) :: n
    n = normalize(name)
    msg = ''
    code = unknown_scheme
    ! aliases
    select case ( trim(n) )
      case ( 'kf', 'kainfritsch' )
        n = 'kain-fritsch'
      case ( 'shinhong' )
        n = 'shin-hong'
      case ( 'xuran', 'xurandall' )
        n = 'xu-randall'
      case ( 'frictionless' )
        n = 'none'
    end select
    select case ( kind )
      case ( 'pbl' )
        code = lookup(pbl)
      case ( 'cumulus' )
        code = lookup(cumulus)
        if ( n == 'kuo' .or. n == 'betts-miller' .or. n == 'bm' ) then
          msg = 'the '//trim(name)//' cumulus scheme has been removed '// &
                '(source in Main/historical)'
          return
        end if
      case ( 'microphysics' )
        code = lookup(micro)
        if ( n == 'subex' ) then
          msg = 'the SUBEX microphysics has been removed (source in '// &
                'Main/historical); its cloud fraction is sundqvist'
          return
        end if
      case ( 'cloud_fraction' )
        code = lookup(cloud)
        if ( n == 'subex' ) then
          msg = 'the subex cloud fraction is now named sundqvist'
          return
        end if
      case ( 'radiation' )
        code = lookup(rad)
      case default
        msg = 'unknown kind of scheme '//trim(kind)
        return
    end select
    if ( code == unknown_scheme ) then
      msg = 'unknown '//trim(kind)//' scheme '''//trim(name)// &
            ''', one of: '//scheme_list(kind)
    end if

    contains

      integer(ik4) function lookup(table) result(c)
        implicit none
        type(scheme_entry), dimension(:), intent(in) :: table
        integer(ik4) :: k
        c = unknown_scheme
        do k = 1, size(table)
          if ( table(k)%name == n ) then
            c = table(k)%code
            exit
          end if
        end do
      end function lookup

  end function scheme_code
  !
  ! Name of the scheme of a kind with the given code
  !
  character(len=maxlen) function scheme_name(kind,code) result(name)
    implicit none
    character(len=*), intent(in) :: kind
    integer(ik4), intent(in) :: code
    select case ( kind )
      case ( 'pbl' )
        name = rlookup(pbl)
      case ( 'cumulus' )
        name = rlookup(cumulus)
      case ( 'microphysics' )
        name = rlookup(micro)
      case ( 'cloud_fraction' )
        name = rlookup(cloud)
      case ( 'radiation' )
        name = rlookup(rad)
      case default
        name = 'unknown'
    end select

    contains

      character(len=maxlen) function rlookup(table) result(nm)
        implicit none
        type(scheme_entry), dimension(:), intent(in) :: table
        integer(ik4) :: k
        nm = 'unknown'
        do k = 1, size(table)
          if ( table(k)%code == code ) then
            nm = table(k)%name
            exit
          end if
        end do
      end function rlookup

  end function scheme_name
  !
  ! The names of the schemes of a kind, comma separated
  !
  function scheme_list(kind) result(list)
    implicit none
    character(len=*), intent(in) :: kind
    character(len=256) :: list
    list = ''
    select case ( kind )
      case ( 'pbl' )
        call append(pbl)
      case ( 'cumulus' )
        call append(cumulus)
      case ( 'microphysics' )
        call append(micro)
      case ( 'cloud_fraction' )
        call append(cloud)
      case ( 'radiation' )
        call append(rad)
    end select

    contains

      subroutine append(table)
        implicit none
        type(scheme_entry), dimension(:), intent(in) :: table
        integer(ik4) :: k
        list = table(1)%name
        do k = 2, size(table)
          list = trim(list)//', '//trim(table(k)%name)
        end do
      end subroutine append

  end function scheme_list

  pure function normalize(name) result(n)
    implicit none
    character(len=*), intent(in) :: name
    character(len=maxlen) :: n
    integer(ik4) :: i, c
    n = adjustl(name)
    do i = 1, len_trim(n)
      c = iachar(n(i:i))
      if ( c >= iachar('A') .and. c <= iachar('Z') ) then
        n(i:i) = achar(c + 32)
      else if ( n(i:i) == '_' ) then
        n(i:i) = '-'
      end if
    end do
  end function normalize

end module mod_scheme_names
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
