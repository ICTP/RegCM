!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of RegCM model.
!
!    RegCM model is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    RegCM model is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with RegCM model.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      module mod_cvaria

      use mod_regcm_param

      implicit none
!
! COMMON /CVARIA/
!
#ifdef MPP1
      real(8) , dimension(ix,kx,jxp) :: diffq , difft , difuu , difuv , &
                                      & omega , qcc , qcten , qvc ,     &
                                      & qvten , tc , td , tten , uc ,   &
                                      & uten , vc , vten , xkc
      real(8) , dimension(ix,jxp) :: pdota , psc , pten
      real(8) , dimension(ix,kx,0:jxp) :: phi
      real(8) , dimension(ix,0:jxp+1) :: psd
      real(8) , dimension(ix,kx,0:jxp+1) :: qc , qv , t , u , v
      real(8) , dimension(ix,kxp1,0:jxp+1) :: qdot
      real(8) , dimension(ix,kx) :: qvcs
#else
      real(8) , dimension(ix,kx,jx) :: diffq , difft , difuu , difuv ,  &
                                     & omega , phi , qc , qcc , qcten , &
                                     & qv , qvc , qvten , t , tc , td , &
                                     & tten , u , uc , uten , v , vc ,  &
                                     & vten , xkc
      real(8) , dimension(ix,jx) :: pdota , psc , psd , pten
      real(8) , dimension(ix,kxp1,jx) :: qdot
      real(8) , dimension(ix,kx) :: qvcs
#endif
!
! COMMON /TRACER/
!
#ifdef MPP1
      real(8) , dimension(ix,kx,0:jxp+1,ntr) :: chi
      real(8) , dimension(ix,kx,jxp,ntr) :: chic , chiten
#else
      real(8) , dimension(ix,kx,jx,ntr) :: chi , chic , chiten
#endif
      end module mod_cvaria
