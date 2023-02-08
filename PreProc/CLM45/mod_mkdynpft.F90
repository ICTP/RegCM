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
module mod_mkdynpft
#ifdef DYNPFT
  use mod_realkinds
  use mod_intkinds
  use mod_dynparam
  use mod_constants
  use mod_grid
  use mod_stdio
  use mod_message
  use mod_rdldtr

  implicit none

  private

  public :: mkdynpft

  character(len=16) , parameter :: varname = 'PCT_PFT'

  real(rkx) , parameter :: vcutoff = 1.0_rkx

  contains

  subroutine mkdynpft(mask,pft,year)
    implicit none
    real(rkx) , dimension(:,:) , intent(in) :: mask
    real(rkx) , dimension(:,:,:) , intent(out) :: pft
    integer(ik4) , intent(in) :: year
    integer(ik4) :: i , j , n , iy , npft
    integer(ik4) , dimension(1) :: il
    character(len=32) :: p1 , p2
    character(len=4) :: cy
    type(globalfile) :: gfile
    character(len=256) :: inpfile
    character(len=8) :: scenario = 'RCP4.5'

    iy = max(min(year,2100),1950)
    p1 = 'dynamic'
    p2 = '.'
    write(cy,'(i0.4)') iy

    if ( iy > 2005 ) then
      select case (dattyp(4:5))
        case ('RF')
          continue
        case ('26')
          p2 = 'SCENARIO'//pthsep//'RCP2.6'
        case ('45')
          p2 = 'SCENARIO'//pthsep//'RCP4.5'
        case ('60')
          p2 = 'SCENARIO'//pthsep//'RCP6.0'
        case ('85', '15')
          p2 = 'SCENARIO'//pthsep//'RCP8.5'
        case default
          write(stderr,*) 'WARNING : Using CMIP5 '//scenario//'scenario !'
          p2 = 'SCENARIO'//pthsep//scenario
      end select
    end if

    npft = size(pft,3)

    inpfile = trim(inpglob)//pthsep//'CLM45'//pthsep//'surface'// &
              pthsep//trim(p1)//pthsep//trim(p2)//pthsep//&
              'mksrf_landuse_'//cy//'.nc'

    call gfopen(gfile,inpfile,xlat,xlon,ds*nsg,roidem,i_band)
    call gfread(gfile,varname,pft,h_missing_value)
    call gfclose(gfile)
    !
    ! Mask and fill grid
    !
    where (pft < d_zero
      pft = h_missing_value
    end where
    do n = 1 , npft
      call bestaround(pft(:,:,n),h_missing_value)
      do i = 1 , iysg
        do j = 1 , jxsg
          if ( mask(j,i) < 0.5_rkx ) then
            pft(j,i,n) = h_missing_value
          else
            pft(j,i,n) = min(max(0,nint(pft(j,i,n))),100)
          end if
        end do
      end do
    end do
    !
    ! Keep only classes with percentage greater than cutoff
    !
    do i = 1 , iysg
      do j = 1 , jxsg
        if ( mask(j,i) > 0.5_rkx ) then
          il = maxloc(pft(j,i,:))
          do n = 1 , npft
            if ( n == il(1) ) cycle
            if ( pft(j,i,n) < vcutoff ) then
              pft(j,i,il(1)) = pft(j,i,il(1)) + pft(j,i,n)
              pft(j,i,n) = d_zero
            end if
          end do
          do n = 1 , npft
            pft(j,i,n) = min(max(0,nint(pft(j,i,n))),100)
          end do
        end if
      end do
    end do
  end subroutine mkdynpft
#endif

end module mod_mkdynpft
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
