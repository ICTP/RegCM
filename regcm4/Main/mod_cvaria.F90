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

      module mod_cvaria

      use mod_regcm_param

      implicit none
!
#ifdef MPP1
      real(8) , dimension(iy,kz,jxp) :: diffq , difft , difuu , difuv , &
                                      & omega , qcc , qcten , qvc ,     &
                                      & qvten , tc , td , tten , uc ,   &
                                      & uten , vc , vten , xkc
      real(8) , dimension(iy,jxp) :: pdota , psc , pten
      real(8) , dimension(iy,kz,0:jxp) :: phi
      real(8) , dimension(iy,0:jxp+1) :: psd
      real(8) , dimension(iy,kz,0:jxp+1) :: qc , qv , t , u , v
      real(8) , dimension(iy,kzp1,0:jxp+1) :: qdot
      real(8) , dimension(iy,kz) :: qvcs
#else
      real(8) , dimension(iy,kz,jx) :: diffq , difft , difuu , difuv ,  &
                                     & omega , phi , qc , qcc , qcten , &
                                     & qv , qvc , qvten , t , tc , td , &
                                     & tten , u , uc , uten , v , vc ,  &
                                     & vten , xkc
      real(8) , dimension(iy,jx) :: pdota , psc , psd , pten
      real(8) , dimension(iy,kzp1,jx) :: qdot
      real(8) , dimension(iy,kz) :: qvcs
#endif
!
#ifdef MPP1
      real(8) , dimension(iy,kz,0:jxp+1,ntr) :: chi
      real(8) , dimension(iy,kz,jxp,ntr) :: chic , chiten
#else
      real(8) , dimension(iy,kz,jx,ntr) :: chi , chic , chiten
#endif

      end module mod_cvaria
