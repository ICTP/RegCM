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

#ifdef CLM

module mod_clm
  !
  ! Storage and parameters for CLM model v3.5
  !
  use mod_intkinds
  use mod_realkinds
  use mod_dynparam , only : domname , pthsep
  use mod_runparams , only : ichem

  implicit none

  integer(ik4) :: r2comm        ! RegCM MPI communicator
  integer(ik4) :: r2cdtime      ! timestep in seconds
  integer(ik4) :: r2cnsrest     ! 0=initial, 1=restart
  integer(ik4) :: r2cnestep     ! final timestep (or day if negative) number
  integer(ik4) :: r2cnelapse    ! # of timesteps (or days if negative)
                           ! to extend a run
  integer(ik4) :: r2cnstep      ! current timestep (ktau)
!
  integer(ik4) :: r2cstart_ymd  ! starting date for run in yearmmdd format
  integer(ik4) :: r2cstart_tod  ! starting time of day for run in seconds
  integer(ik4) :: r2cstop_ymd   ! stopping date for run in yearmmdd format
  integer(ik4) :: r2cstop_tod   ! stopping time of day for run in seconds
  integer(ik4) :: r2cref_ymd    ! reference date for time coordinate
                           ! in yearmmdd format
  integer(ik4) :: r2cref_tod    ! reference time of day for time
                           ! coordinate in seconds

  character(len=32) :: r2cclndr ! Calendar type ('NO_LEAP' or 'GREGORIAN')

  logical :: r2cwrtdia     ! write output true/false
  integer(ik4) :: r2cirad       ! radiation calculation
  integer(ik4) :: r2cmss_irt    ! NCAR mass store retention time set to 0

  integer(ik4) :: r2clsmlon     ! number of longitude points
  integer(ik4) :: r2clsmlat     ! number of latitude points
  real(rkx) :: r2cdx         ! model resolution (m)
  real(rkx) :: r2carea       ! area of each grid cell (km^2)
  real(rkx) :: r2cedgen      ! N edge of grid
  real(rkx) :: r2cedges      ! S edge of grid
  real(rkx) :: r2cedgee      ! E edge of grid
  real(rkx) :: r2cedgew      ! W edge of grid

  real(rkx) :: r2ceccen      ! orbital eccentricity
  real(rkx) :: r2cobliqr     ! Earths obliquity in rad
  real(rkx) :: r2clambm0     ! Mean long of perihelion at
                           ! vernal equinox (radians)
  real(rkx) :: r2cmvelpp     ! moving vernal equinox long
  real(rkx) :: r2ceccf       ! Earth-sun distance factor (1/r)**2

  logical :: r2cdoalb      ! time for next albedo call

  integer(ik4) :: c2rcnts

  real(rkx) , pointer , dimension(:,:) :: r2ctb
  real(rkx) , pointer , dimension(:,:) :: r2ctb_all
  real(rkx) , pointer , dimension(:,:) :: r2cqb
  real(rkx) , pointer , dimension(:,:) :: r2cqb_all
  real(rkx) , pointer , dimension(:,:) :: r2czga
  real(rkx) , pointer , dimension(:,:) :: r2czga_all
  real(rkx) , pointer , dimension(:,:) :: r2cpsb
  real(rkx) , pointer , dimension(:,:) :: r2cpsb_all
  real(rkx) , pointer , dimension(:,:) :: r2cuxb
  real(rkx) , pointer , dimension(:,:) :: r2cuxb_all
  real(rkx) , pointer , dimension(:,:) :: r2cvxb
  real(rkx) , pointer , dimension(:,:) :: r2cvxb_all
  real(rkx) , pointer , dimension(:,:) :: r2crnc
  real(rkx) , pointer , dimension(:,:) :: r2crnc_all
  real(rkx) , pointer , dimension(:,:) :: r2crnnc
  real(rkx) , pointer , dimension(:,:) :: r2crnnc_all
  real(rkx) , pointer , dimension(:,:) :: r2csols
  real(rkx) , pointer , dimension(:,:) :: r2csols_all
  real(rkx) , pointer , dimension(:,:) :: r2csoll
  real(rkx) , pointer , dimension(:,:) :: r2csoll_all
  real(rkx) , pointer , dimension(:,:) :: r2csolsd
  real(rkx) , pointer , dimension(:,:) :: r2csolsd_all
  real(rkx) , pointer , dimension(:,:) :: r2csolld
  real(rkx) , pointer , dimension(:,:) :: r2csolld_all
  real(rkx) , pointer , dimension(:,:) :: r2cflwd
  real(rkx) , pointer , dimension(:,:) :: r2cflwd_all
  real(rkx) , pointer , dimension(:,:) :: r2ccosz_all
  real(rkx) , pointer , dimension(:,:) :: r2cxlat_all     ! xlat in radians
  real(rkx) , pointer , dimension(:,:) :: r2cxlon_all     ! xlon in radians
  real(rkx) , pointer , dimension(:,:) :: r2cxlatd_all    ! xlat in degrees
  real(rkx) , pointer , dimension(:,:) :: r2cxlond_all    ! xlon in degrees
  real(rkx) , pointer , dimension(:,:) :: r2cxlat
  real(rkx) , pointer , dimension(:,:) :: r2cxlon
  real(rkx) , pointer , dimension(:,:) :: r2cxlatd
  real(rkx) , pointer , dimension(:,:) :: r2cxlond
  real(rkx) , pointer , dimension(:,:) :: c2rtgb
  real(rkx) , pointer , dimension(:,:) :: c2rsenht
  real(rkx) , pointer , dimension(:,:) :: c2rlatht
  real(rkx) , pointer , dimension(:,:) :: c2ralbdirs
  real(rkx) , pointer , dimension(:,:) :: c2ralbdirl
  real(rkx) , pointer , dimension(:,:) :: c2ralbdifs
  real(rkx) , pointer , dimension(:,:) :: c2ralbdifl
  real(rkx) , pointer , dimension(:,:) :: c2rtaux
  real(rkx) , pointer , dimension(:,:) :: c2rtauy
  real(rkx) , pointer , dimension(:,:) :: c2ruvdrag
  real(rkx) , pointer , dimension(:,:) :: c2rlsmask
  real(rkx) , pointer , dimension(:,:) :: c2rtgbb
  real(rkx) , pointer , dimension(:,:) :: c2rsnowc
  real(rkx) , pointer , dimension(:,:) :: c2rtest
  real(rkx) , pointer , dimension(:,:) :: c2r2mt
  real(rkx) , pointer , dimension(:,:) :: c2r2mq
  real(rkx) , pointer , dimension(:,:) :: c2rtlef
  real(rkx) , pointer , dimension(:,:) :: c2ru10
  real(rkx) , pointer , dimension(:,:) :: c2rsm10cm
  real(rkx) , pointer , dimension(:,:) :: c2rsm1m
  real(rkx) , pointer , dimension(:,:) :: c2rsmtot
  real(rkx) , pointer , dimension(:,:) :: c2rinfl
  real(rkx) , pointer , dimension(:,:) :: c2rro_sur
  real(rkx) , pointer , dimension(:,:) :: c2rro_sub
  real(rkx) , pointer , dimension(:,:) :: c2rfracsno
  real(rkx) , pointer , dimension(:,:) :: c2rfvegnosno

  integer(ik4) , pointer , dimension(:,:) :: c2rprocmap
  integer(ik4) , pointer , dimension(:) :: c2rngc
  integer(ik4) , pointer , dimension(:) :: c2rdisps
!
  real(rkx) , pointer , dimension(:,:) :: rs2d
  real(rkx) , pointer , dimension(:,:) :: ra2d

end module mod_clm
!
#endif
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
