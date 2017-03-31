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
#ifdef CN
  use mod_realkinds
  use mod_intkinds
  use mod_dynparam
  use mod_grid
  use mod_message
  use mod_rdldtr

  implicit none

  private

  public :: mkdynpft

  character(len=16) , parameter :: varname = 'PCT_PFT'

  real(rkx) :: vcutoff = 19.0_rkx

  contains

  subroutine mkdynpft(dynpft,year)
    implicit none
    real(rkx) , dimension(:,:,:) , intent(out) :: dynpft
    integer(ik4) , intent(in) :: year
    integer(ik4) :: i , j , n , iy , npft
    integer(ik4) , dimension(1) :: il
    character(len=32) :: p1 , p2
    character(len=4) :: cy
    type(globalfile) :: gfile
    character(len=256) :: inpfile

    iy = min(max(year,2100),1950)
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
          if ( dattyp /= "EIN15" .and. &
               dattyp(1:4) /= "NNRP" .and. &
               dattyp /= "JRA55" ) then
            call die(__FILE__, &
              'Dynamic landuse only supported for CMIP5',__LINE__)
          end if
      end select
    end if

    npft = size(dynpft,3)

    call gfopen(gfile,inpfile,xlat,xlon,ds*nsg,i_band)
    call gfread(gfile,varname,dynpft,h_missing_value)

    do i = 1 , iysg
      do j = 1 , jxsg
        if ( dynpft(j,i,1) > h_missing_value ) then
          il = maxloc(dynpft(j,i,:))
          do n = 1 , npft
            if ( n == il(1) ) cycle
            if ( dynpft(j,i,n) < vcutoff ) then
              dynpft(j,i,il(1)) = dynpft(j,i,il(1)) + dynpft(j,i,n)
              dynpft(j,i,n) = d_zero
            end if
          end do
          do n = 1 , npft
            dynpft(j,i,n) = min(dynpft(j,i,n),100.0_rkx)
          end do
        end if
      end do
    end do
    call gfclose(gfile)
  end subroutine mkdynpft
#endif

end module mod_mkdynpft
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
