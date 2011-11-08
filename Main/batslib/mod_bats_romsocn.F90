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

  subroutine romsocndrv(jstart,jend,istart,iend,ktau)
    implicit none
    integer , intent(in) :: jstart , jend , istart , iend
    integer(8) , intent(in) :: ktau
!
    integer i , j , n
    character (len=64) :: subroutine_name='romsocndrv'
    integer :: idindx=0
!
    call time_begin(subroutine_name,idindx)
!
    do i = istart , iend
      do j = jstart , jend
        do n = 1 , nnsg
          ! Feed back ground temperature (in Kelvin) 
          if (sst2d(i,j) .lt. MISSING_R8) then  
            tgrd(n,j,i) = sst2d(i,j)
            tgbrd(n,j,i) = sst2d(i,j)
          end if
        end do
      end do
    end do
!
    call time_end(subroutine_name,idindx)
  end subroutine romsocndrv
!
  subroutine print_matrix_r8(inp, iskip, jskip, pet, header)
    implicit none
!
    real(dp), intent(in) :: inp(:,:)
    integer, intent(in) ::  iskip, jskip, pet
    character(len=*), intent(in) :: header
!
    integer :: i, j, imin, imax, jmin, jmax
    character(100) :: fmt_123
!
    imin = lbound(inp, dim=1)
    imax = ubound(inp, dim=1)
    jmin = lbound(inp, dim=2)
    jmax = ubound(inp, dim=2)
!
    write(6, fmt="('PET(',I2,') - ',A)") pet, trim(header)
!
    write(fmt_123, fmt="('(/, 5X, ', I3, 'I10)')") (jmax-jmin)+1
    write(6, fmt=trim(fmt_123))  (j, j=jmin, jmax, jskip)
!  
    write(fmt_123, fmt="('(I5, ', I3, 'F10.2)')") jmax
    do i=imin, imax, iskip
      write(6, fmt=trim(fmt_123)) i, (inp(i,j),j=jmin, jmax, jskip)
    end do
!
    return
  end subroutine print_matrix_r8
!
end module mod_bats_romsocn
