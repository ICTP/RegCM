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
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
! This is the Regional Climatic Model RegCM from ICTP Trieste
!
!  References :
!
!     1. model:
!
!     Anthes, R. A., and T. T. Warner, 1978: Development of
!     hydrodynamic models suitable for air pollution and
!     other mesometeorological studies. Mon. Wea. Rev.,
!     106, 1045-1078.
!
!     2. cumulus parameterization :
!
!     Anthes, R. A., 1977: A cumulus parameterization scheme
!     utilizing a one-dimensional cloud model. Mon. Wea.
!     Rev., 105, 270-286.
!
!     Kuo, Y.-H., 1983: A diagnostic case study of the effects
!     of deep extratropical convection on the large-scale
!     temperature and moisture structure. PH.D. thesis,
!     Department of Meteorology, the Pennsylvania State
!     University, 222 pp.
!
!     Grell,
!
!     3. explicit moisture :
!
!     Hsie, E.-Y., 1983: Frontogenesis in a moist atmosphere.
!     PH.D. thesis, Department of Meteorology, the
!     Pennsylvania State University, 251 pp.
!
!     4. pbl parameterization :
!
!     Holtslag, De Bruijn and Pan - MWR - 8/90
!
!     5. radiation parameterization :
!
!     CCM2 radiation column model, bruce briegleb, jan '92
!
!     CCM3 radiation column model, NCAR/TN-422+PPR, Description
!     of the NCAR CCM,      J. T. Kiehl, J. Hack et al.,
!     introduced by   Filippo Giorgi, N. Keiichi, Yun Qian
!
!     CCM3.6.6 code introduced by Xunqiang Bi, 2000
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

subroutine myabort
  use mod_dynparam , only : mycomm
  use mod_intkinds
  implicit none
  include 'mpif.h'
  integer(ik4) :: ierr
  call mpi_abort(mycomm,1,ierr)
end subroutine myabort

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
! The main program is in either regcm_nonesmf.F90 or regcm_esmf.F90
! Check wichever applies
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
