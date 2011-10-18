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

module mod_bats_romsocn
!
! ROMS ocean model
!
  use mod_realkinds
  use mod_dynparam
  use mod_service
  use mod_bats_common
!
  private
!
  public :: romsocndrv
  public :: allocate_mod_bats_romsocn
!
  real(dp), public, pointer, dimension(:,:) :: sst2d
  real(dp), parameter :: MISSING_R8 = 1.0d20 
!
  contains

  subroutine allocate_mod_bats_romsocn()
    implicit none
    call getmem2d(sst2d,1,iy,1,jxp,'roms:sst2d') 
    sst2d = MISSING_R8
  end subroutine allocate_mod_bats_romsocn

  subroutine romsocndrv(j,istart,iend,ktau)
!
    implicit none
!
    integer i, n
    
!
    character (len=64) :: subroutine_name='romsocndrv'
    integer :: idindx=0
!
    integer , intent(in) :: j , istart , iend
    integer(8) , intent(in) :: ktau
!
    call time_begin(subroutine_name,idindx)
!
    do i = istart , iend
      do n = 1 , nnsg
!       Feed back ground temperature (in Kelvin) 
        if (sst2d(i,j) .lt. MISSING_R8) then  
          tg1d(n,i) = sst2d(i,j)
          tgb1d(n,i) = sst2d(i,j)
        end if
      end do
    end do
!
    call time_end(subroutine_name,idindx)
  end subroutine romsocndrv
!
end module mod_bats_romsocn
