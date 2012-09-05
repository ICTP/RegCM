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

  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_bats_common
  use mod_mppparam

  public

  integer(ik4) , pointer , dimension(:,:,:) :: idep_io
  real(rk8) , pointer , dimension(:,:,:) :: dhlake1_io
  real(rk8) , pointer , dimension(:,:,:) :: eta_io
  real(rk8) , pointer , dimension(:,:,:) :: hi_io
  real(rk8) , pointer , dimension(:,:,:) :: aveice_io
  real(rk8) , pointer , dimension(:,:,:) :: hsnow_io
  real(rk8) , pointer , dimension(:,:,:) :: evl_io
  real(rk8) , pointer , dimension(:,:,:,:) :: tlak_io

  contains

    subroutine allocate_mod_bats_mppio(lakemod)
      implicit none
      integer(ik4) , intent(in) :: lakemod
      if ( lakemod == 1 ) then
        if ( myid == iocpu ) then
          call getmem3d(dhlake1_io,1,nnsg,jdot1,jdot2,idot1,idot2, &
                        'bats_mppio:dhlake1_io')
          call getmem3d(idep_io,1,nnsg,jcross1,jcross2,icross1,icross2, &
                        'bats_mppio:idep_io')
          call getmem3d(eta_io,1,nnsg,jcross1,jcross2,icross1,icross2, &
                        'bats_mppio:eta_io')
          call getmem3d(hi_io,1,nnsg,jcross1,jcross2,icross1,icross2, &
                        'bats_mppio:hi_io')
          call getmem3d(aveice_io,1,nnsg,jcross1,jcross2,icross1,icross2, &
                        'bats_mppio:aveice_io')
          call getmem3d(hsnow_io,1,nnsg,jcross1,jcross2,icross1,icross2, &
                        'bats_mppio:hsnow_io')
          call getmem3d(evl_io,1,nnsg,jcross1,jcross2,icross1,icross2, &
                        'bats_mppio:evl_io')
          call getmem4d(tlak_io,1,nnsg,jcross1,jcross2,icross1,icross2, &
                        1,ndpmax,'bats_mppio:tlak_io')
        end if
      end if
    end subroutine allocate_mod_bats_mppio

end module mod_bats_mppio
