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
      real(8) ,allocatable, dimension(:,:,:) :: diffq , difft , difuu , difuv , &
                                      & omega , qcc , qcten , qvc ,     &
                                      & qvten , tc , td , tten , uc ,   &
                                      & uten , vc , vten , xkc
      real(8) ,allocatable, dimension(:,:) :: pdota , psc , pten
      real(8) ,allocatable, dimension(:,:,:) :: phi
      real(8) ,allocatable, dimension(:,:) :: psd
      real(8) ,allocatable, dimension(:,:,:) :: qc , qv , t , u , v
      real(8) ,allocatable, dimension(:,:,:) :: qdot
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
      real(8) ,allocatable, dimension(:,:,:,:) :: chi
      real(8) ,allocatable, dimension(:,:,:,:) :: chic , chiten
#else
      real(8) , dimension(iy,kz,jx,ntr) :: chi , chic , chiten
#endif

contains 
	subroutine allocate_mod_cvaria

#ifdef MPP1
	allocate(diffq(iy,kz,jxp))
	allocate(difft(iy,kz,jxp))
	allocate(difuu(iy,kz,jxp))
	allocate(difuv(iy,kz,jxp))
	allocate(omega(iy,kz,jxp))
	allocate(qcc(iy,kz,jxp))
	allocate(qcten(iy,kz,jxp))
	allocate(qvc(iy,kz,jxp))
	allocate(qvten(iy,kz,jxp))
	allocate(tc(iy,kz,jxp))
	allocate(td(iy,kz,jxp))
	allocate(tten(iy,kz,jxp))
	allocate(uc(iy,kz,jxp))
	allocate(uten(iy,kz,jxp))
	allocate(vc(iy,kz,jxp))
	allocate(vten(iy,kz,jxp))
	allocate(xkc(iy,kz,jxp))


	allocate(pdota(iy,jxp))
	allocate(psc(iy,jxp))
	allocate(pten(iy,jxp))
	allocate(phi(iy,kz,0:jxp))
	allocate(psd(iy,0:jxp+1))

        allocate(qc(iy,kz,0:jxp+1))
        allocate(qv(iy,kz,0:jxp+1))
        allocate(t(iy,kz,0:jxp+1))
        allocate(u(iy,kz,0:jxp+1))
        allocate(v(iy,kz,0:jxp+1))

	allocate(qdot(iy,kzp1,0:jxp+1))
	allocate(chi(iy,kz,0:jxp+1,ntr))
	allocate(chic(iy,kz,jxp,ntr))
	allocate(chiten(iy,kz,jxp,ntr))

#endif
      end  subroutine allocate_mod_cvaria
      end module mod_cvaria
