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
!
! Storage for the prognostic variables at tau+1,
!     decoupled variables, diagnostic variables and
!     working spaces needed in the model.
!
      use mod_dynparam

      implicit none
!
      real(8) , allocatable , dimension(:,:,:) :: diffq , difft ,       &
                                      & difuu , difuv ,                 &
                                      & omega , qcc , qcten , qvc ,     &
                                      & qvten , tc , td , tten , uc ,   &
                                      & uten , vc , vten , xkc
      real(8) , allocatable , dimension(:,:) :: pdota , psc , pten
      real(8) , allocatable , dimension(:,:,:) :: phi
      real(8) , allocatable , dimension(:,:) :: psd
      real(8) , allocatable , dimension(:,:,:) :: qc , qv , t , u , v
      real(8) , allocatable , dimension(:,:,:) :: qdot
      real(8) , allocatable , dimension(:,:) :: qvcs
!
      real(8) , allocatable , dimension(:,:,:,:) :: chi
      real(8) , allocatable , dimension(:,:,:,:) :: chic , chiten

      contains 

        subroutine allocate_mod_cvaria
        implicit none
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
#else
        allocate(diffq(iy,kz,jx))
        allocate(difft(iy,kz,jx))
        allocate(difuu(iy,kz,jx))
        allocate(difuv(iy,kz,jx))
        allocate(omega(iy,kz,jx))
        allocate(qcc(iy,kz,jx))
        allocate(qcten(iy,kz,jx))
        allocate(qvc(iy,kz,jx))
        allocate(qvten(iy,kz,jx))
        allocate(tc(iy,kz,jx))
        allocate(td(iy,kz,jx))
        allocate(tten(iy,kz,jx))
        allocate(uc(iy,kz,jx))
        allocate(uten(iy,kz,jx))
        allocate(vc(iy,kz,jx))
        allocate(vten(iy,kz,jx))
        allocate(xkc(iy,kz,jx))
        allocate(pdota(iy,jx))
        allocate(psc(iy,jx))
        allocate(pten(iy,jx))
        allocate(phi(iy,kz,jx))
        allocate(psd(iy,jx))
        allocate(qc(iy,kz,jx))
        allocate(qv(iy,kz,jx))
        allocate(t(iy,kz,jx))
        allocate(u(iy,kz,jx))
        allocate(v(iy,kz,jx))
        allocate(chi(iy,kz,jx,ntr))
        allocate(chic(iy,kz,jx,ntr))
        allocate(chiten(iy,kz,jx,ntr))
        allocate(qdot(iy,kzp1,jx))
#endif
        allocate(qvcs(iy,kz))
!
        diffq = 0.0D0
        difft = 0.0D0
        difuu = 0.0D0
        difuv = 0.0D0
        omega = 0.0D0
        qcc = 0.0D0
        qcten = 0.0D0
        qvc = 0.0D0
        qvten = 0.0D0
        tc = 0.0D0
        td = 0.0D0
        tten = 0.0D0
        uc = 0.0D0
        uten = 0.0D0
        vc = 0.0D0
        vten = 0.0D0
        xkc = 0.0D0
        pdota = 0.0D0
        psc = 0.0D0
        pten = 0.0D0
        phi = 0.0D0
        psd = 0.0D0
        qc = 0.0D0
        qv = 0.0D0
        t = 0.0D0
        u = 0.0D0
        v = 0.0D0
        chi = 0.0D0
        chic = 0.0D0
        chiten = 0.0D0
        qdot = 0.0D0
        qvcs = 0.0D0
      end  subroutine allocate_mod_cvaria
!
      end module mod_cvaria
