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

module mod_maputils

  use mod_constants
  use mod_intkinds
  use mod_realkinds
  use mod_projections

  private

  public :: getcoord , corpar , mappar , anyprojparams , regcm_projection

  contains

  subroutine getcoord(pjpara,lon,lat,pj,jx,iy)
    implicit none
    type(anyprojparams) , intent(in) :: pjpara
    type(regcm_projection) , intent(out) :: pj
    integer(ik4) , intent(in) :: jx , iy
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: lat , lon
    integer :: i , j

    call pj%initialize(pjpara)
    do i = 1 , iy
      do j = 1 , jx
        call pj%ijll(real(j,rkx),real(i,rkx),lat(j,i),lon(j,i))
      end do
    end do
  end subroutine getcoord

  subroutine corpar(lat,coriol)
    implicit none
    real(rkx) , pointer , dimension(:,:) , intent(in) :: lat
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: coriol
    coriol = real(eomeg2*sin(lat*degrad),rkx)
  end subroutine corpar

  subroutine mappar(pj,xlat,xlon,mapf)
    implicit none
    type(regcm_projection) , intent(in) :: pj
    real(rkx) , pointer , dimension(:,:) , intent(in) :: xlat , xlon
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: mapf
    call pj%mapfac(xlat,xlon,mapf)
  end subroutine mappar

end module mod_maputils

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
