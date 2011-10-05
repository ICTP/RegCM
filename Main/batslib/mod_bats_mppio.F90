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

module mod_bats_mppio

  use mod_realkinds
  use mod_dynparam
  use mod_bats_common

  public

  integer , pointer , dimension(:,:,:) :: idep2d_io
  real(dp) , pointer , dimension(:,:,:) :: dhlake1_io
  real(dp) , pointer , dimension(:,:,:) :: eta2d_io
  real(dp) , pointer , dimension(:,:,:) :: hi2d_io
  real(dp) , pointer , dimension(:,:,:) :: aveice2d_io
  real(dp) , pointer , dimension(:,:,:) :: hsnow2d_io
  real(dp) , pointer , dimension(:,:,:) :: evl2d_io
  real(dp) , pointer , dimension(:,:,:,:) :: tlak3d_io

  contains

    subroutine allocate_mod_bats_mppio(lakemod)
      implicit none
      integer , intent(in) :: lakemod
      if ( lakemod == 1 ) then
        if ( myid == 0 ) then
          call getmem3d(dhlake1_io,1,nnsg,1,iy,1,jx,'mod_mppio:dhlake1_io')
          call getmem3d(idep2d_io,1,nnsg,1,iym1,1,jx,'mod_mppio:idep2d_io')
          call getmem3d(eta2d_io,1,nnsg,1,iym1,1,jx,'mod_mppio:eta2d_io')
          call getmem3d(hi2d_io,1,nnsg,1,iym1,1,jx,'mod_mppio:hi2d_io')
          call getmem3d(aveice2d_io,1,nnsg,1,iym1,1,jx,'mod_mppio:aveice2d_io')
          call getmem3d(hsnow2d_io,1,nnsg,1,iym1,1,jx,'mod_mppio:hsnow2d_io')
          call getmem3d(evl2d_io,1,nnsg,1,iym1,1,jx,'mod_mppio:evl2d_io')
          call getmem4d(tlak3d_io,1,ndpmax,1,nnsg, &
                                  1,iym1,1,jx,'mod_mppio:tlak3d_io')
        end if
      end if
    end subroutine allocate_mod_bats_mppio

end module mod_bats_mppio
