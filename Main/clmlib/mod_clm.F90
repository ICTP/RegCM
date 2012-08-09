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
  use mod_dynparam
  use mod_realkinds
  use mod_constants
  use mod_memutil
  use mod_runparams , only : ichem , iocnflx
!
  implicit none

  integer :: imask
  real(dp) :: clmfrq
!
  integer :: r2comm        ! RegCM MPI communicator
  integer :: r2cdtime      ! timestep in seconds
  integer :: r2cnsrest     ! 0=initial, 1=restart
  integer :: r2cnestep     ! final timestep (or day if negative) number
  integer :: r2cnelapse    ! # of timesteps (or days if negative) 
                           ! to extend a run
  integer :: r2cnstep      ! current timestep (ktau)
!
  integer :: r2cstart_ymd  ! starting date for run in yearmmdd format
  integer :: r2cstart_tod  ! starting time of day for run in seconds
  integer :: r2cstop_ymd   ! stopping date for run in yearmmdd format
  integer :: r2cstop_tod   ! stopping time of day for run in seconds
  integer :: r2cref_ymd    ! reference date for time coordinate 
                           ! in yearmmdd format
  integer :: r2cref_tod    ! reference time of day for time 
                           ! coordinate in seconds

  character(len=32) :: r2cclndr ! Calendar type ('NO_LEAP' or 'GREGORIAN')

  logical :: r2cwrtdia     ! write output true/false
  integer :: r2cirad       ! radiation calculation 
  integer :: r2cmss_irt    ! NCAR mass store retention time set to 0

  integer :: r2clsmlon     ! number of longitude points
  integer :: r2clsmlat     ! number of latitude points
  real(dp) :: r2cdx         ! model resolution (m)
  real(dp) :: r2carea       ! area of each grid cell (km^2)
  real(dp) :: r2cedgen      ! N edge of grid
  real(dp) :: r2cedges      ! S edge of grid
  real(dp) :: r2cedgee      ! E edge of grid
  real(dp) :: r2cedgew      ! W edge of grid

  real(dp) :: r2ceccen      ! orbital eccentricity
  real(dp) :: r2cobliqr     ! Earths obliquity in rad
  real(dp) :: r2clambm0     ! Mean long of perihelion at
                           ! vernal equinox (radians)
  real(dp) :: r2cmvelpp     ! moving vernal equinox long
  real(dp) :: r2ceccf       ! Earth-sun distance factor (1/r)**2

  logical :: r2cdoalb      ! time for next albedo call

  integer :: c2rcnts

  real(dp) , pointer , dimension(:,:) :: r2ctb
  real(dp) , pointer , dimension(:,:) :: r2ctb_all
  real(dp) , pointer , dimension(:,:) :: r2cqb
  real(dp) , pointer , dimension(:,:) :: r2cqb_all
  real(dp) , pointer , dimension(:,:) :: r2czga
  real(dp) , pointer , dimension(:,:) :: r2czga_all
  real(dp) , pointer , dimension(:,:) :: r2cpsb
  real(dp) , pointer , dimension(:,:) :: r2cpsb_all
  real(dp) , pointer , dimension(:,:) :: r2cuxb
  real(dp) , pointer , dimension(:,:) :: r2cuxb_all
  real(dp) , pointer , dimension(:,:) :: r2cvxb
  real(dp) , pointer , dimension(:,:) :: r2cvxb_all
  real(dp) , pointer , dimension(:,:) :: r2crnc
  real(dp) , pointer , dimension(:,:) :: r2crnc_all
  real(dp) , pointer , dimension(:,:) :: r2crnnc
  real(dp) , pointer , dimension(:,:) :: r2crnnc_all
  real(dp) , pointer , dimension(:,:) :: r2csols
  real(dp) , pointer , dimension(:,:) :: r2csols_all
  real(dp) , pointer , dimension(:,:) :: r2csoll
  real(dp) , pointer , dimension(:,:) :: r2csoll_all
  real(dp) , pointer , dimension(:,:) :: r2csolsd
  real(dp) , pointer , dimension(:,:) :: r2csolsd_all
  real(dp) , pointer , dimension(:,:) :: r2csolld
  real(dp) , pointer , dimension(:,:) :: r2csolld_all
  real(dp) , pointer , dimension(:,:) :: r2cflwd
  real(dp) , pointer , dimension(:,:) :: r2cflwd_all
  real(dp) , pointer , dimension(:,:) :: r2ccosz_all
  real(dp) , pointer , dimension(:,:) :: r2cxlat_all     ! xlat in radians
  real(dp) , pointer , dimension(:,:) :: r2cxlon_all     ! xlon in radians
  real(dp) , pointer , dimension(:,:) :: r2cxlatd_all    ! xlat in degrees
  real(dp) , pointer , dimension(:,:) :: r2cxlond_all    ! xlon in degrees
  real(dp) , pointer , dimension(:,:) :: r2cxlat
  real(dp) , pointer , dimension(:,:) :: r2cxlon
  real(dp) , pointer , dimension(:,:) :: r2cxlatd
  real(dp) , pointer , dimension(:,:) :: r2cxlond
  real(dp) , pointer , dimension(:,:) :: c2rtgb
  real(dp) , pointer , dimension(:,:) :: c2rsenht
  real(dp) , pointer , dimension(:,:) :: c2rlatht
  real(dp) , pointer , dimension(:,:) :: c2ralbdirs
  real(dp) , pointer , dimension(:,:) :: c2ralbdirl
  real(dp) , pointer , dimension(:,:) :: c2ralbdifs
  real(dp) , pointer , dimension(:,:) :: c2ralbdifl
  real(dp) , pointer , dimension(:,:) :: c2rtaux
  real(dp) , pointer , dimension(:,:) :: c2rtauy
  real(dp) , pointer , dimension(:,:) :: c2ruvdrag
  real(dp) , pointer , dimension(:,:) :: c2rlsmask
  real(dp) , pointer , dimension(:,:) :: c2rtgbb
  real(dp) , pointer , dimension(:,:) :: c2rsnowc
  real(dp) , pointer , dimension(:,:) :: c2rtest
  real(dp) , pointer , dimension(:,:) :: c2r2mt
  real(dp) , pointer , dimension(:,:) :: c2r2mq
  real(dp) , pointer , dimension(:,:) :: c2rtlef
  real(dp) , pointer , dimension(:,:) :: c2ru10
  real(dp) , pointer , dimension(:,:) :: c2rsm10cm
  real(dp) , pointer , dimension(:,:) :: c2rsm1m
  real(dp) , pointer , dimension(:,:) :: c2rsmtot
  real(dp) , pointer , dimension(:,:) :: c2rinfl
  real(dp) , pointer , dimension(:,:) :: c2rro_sur
  real(dp) , pointer , dimension(:,:) :: c2rro_sub
  real(dp) , pointer , dimension(:,:) :: c2rfracsno      
  real(dp) , pointer , dimension(:,:) :: c2rfvegnosno 
  real(dp) , pointer , dimension(:,:) :: voc_em
  real(dp) , pointer , dimension(:,:,:) :: dep_vels

  integer , pointer , dimension(:,:) :: c2rprocmap
  integer , pointer , dimension(:) :: c2rngc
  integer , pointer , dimension(:) :: c2rdisps
!
  ! Direct solar rad incident on surface (<0.7)
  real(dp) , pointer , dimension(:,:) :: sols2d
  ! Direct solar rad incident on surface (>=0.7)
  real(dp) , pointer , dimension(:,:) :: soll2d
  ! Diffuse solar rad incident on surface (<0.7)
  real(dp) , pointer , dimension(:,:) :: solsd2d
  ! Diffuse solar rad incident on surface (>=0.7)
  real(dp) , pointer , dimension(:,:) :: solld2d
  real(dp) , pointer , dimension(:,:) :: aldirs2d
  real(dp) , pointer , dimension(:,:) :: aldirl2d
  real(dp) , pointer , dimension(:,:) :: aldifs2d
  real(dp) , pointer , dimension(:,:) :: aldifl2d
  real(dp) , pointer , dimension(:,:) :: lndcat2d
  real(dp) , pointer , dimension(:,:) :: rs2d
  real(dp) , pointer , dimension(:,:) :: ra2d
  ! 2 meter specific humidity
  real(dp) , pointer , dimension(:,:) :: q2d

  real(dp) , pointer , dimension(:,:) :: htf      ! mddom_io%ht
  real(dp) , pointer , dimension(:,:) :: lndcatf  ! mddom_io%lndcat
!
  contains
!
  subroutine allocate_mod_clm(n_tr,igases)

    implicit none

    integer, intent(in) :: igases
    integer, intent(in) :: n_tr

    call getmem2d(r2ctb,1,jxp,1,iyp,'clm:r2ctb')
    call getmem2d(r2cqb,1,jxp,1,iyp,'clm:r2cqb')
    call getmem2d(r2czga,1,jxp,1,iyp,'clm:r2czga')
    call getmem2d(r2cpsb,1,jxp,1,iyp,'clm:r2cpsb')
    call getmem2d(r2cuxb,1,jxp,1,iyp,'clm:r2cuxb')
    call getmem2d(r2cvxb,1,jxp,1,iyp,'clm:r2cvxb')
    call getmem2d(r2crnc,1,jxp,1,iyp,'clm:r2crnc')
    call getmem2d(r2crnnc,1,jxp,1,iyp,'clm:r2crnnc')
    call getmem2d(r2csols,1,jxp,1,iyp,'clm:r2csols')
    call getmem2d(r2csoll,1,jxp,1,iyp,'clm:r2csoll')
    call getmem2d(r2csolsd,1,jxp,1,iyp,'clm:r2csolsd')
    call getmem2d(r2csolld,1,jxp,1,iyp,'clm:r2csolld')
    call getmem2d(r2cflwd,1,jxp,1,iyp,'clm:r2cflwd')
    call getmem2d(r2cxlat,1,jxp,1,iyp,'clm:r2cxlat')
    call getmem2d(r2cxlon,1,jxp,1,iyp,'clm:r2cxlon')
    call getmem2d(r2cxlatd,1,jxp,1,iyp,'clm:r2cxlatd')
    call getmem2d(r2cxlond,1,jxp,1,iyp,'clm:r2cxlond')

    call getmem2d(r2ctb_all,1,jx,1,iy,'clm:r2ctb_all')
    call getmem2d(r2cqb_all,1,jx,1,iy,'clm:r2cqb_all')
    call getmem2d(r2czga_all,1,jx,1,iy,'clm:r2czga_all')
    call getmem2d(r2cpsb_all,1,jx,1,iy,'clm:r2cpsb_all')
    call getmem2d(r2cuxb_all,1,jx,1,iy,'clm:r2cuxb_all')
    call getmem2d(r2cvxb_all,1,jx,1,iy,'clm:r2cvxb_all')
    call getmem2d(r2crnc_all,1,jx,1,iy,'clm:r2crnc_all')
    call getmem2d(r2crnnc_all,1,jx,1,iy,'clm:r2crnnc_all')
    call getmem2d(r2csols_all,1,jx,1,iy,'clm:r2csols_all')
    call getmem2d(r2csoll_all,1,jx,1,iy,'clm:r2csoll_all')
    call getmem2d(r2csolsd_all,1,jx,1,iy,'clm:r2csolsd_all')
    call getmem2d(r2csolld_all,1,jx,1,iy,'clm:r2csolld_all')
    call getmem2d(r2cflwd_all,1,jx,1,iy,'clm:r2cflwd_all')
    call getmem2d(r2ccosz_all,1,jx,1,iy,'clm:r2ccosz_all')
    call getmem2d(r2cxlat_all,1,jx,1,iy,'clm:r2cxlat_all')
    call getmem2d(r2cxlon_all,1,jx,1,iy,'clm:r2cxlon_all')
    call getmem2d(r2cxlatd_all,1,jx,1,iy,'clm:r2cxlatd_all')
    call getmem2d(r2cxlond_all,1,jx,1,iy,'clm:r2cxlond_all')

    call getmem2d(c2rtgb,1,jx,1,iy,'clm:c2rtgb')
    call getmem2d(c2rsenht,1,jx,1,iy,'clm:c2rsenht')
    call getmem2d(c2rlatht,1,jx,1,iy,'clm:c2rlatht')
    call getmem2d(c2ralbdirs,1,jx,1,iy,'clm:c2ralbdirs')
    call getmem2d(c2ralbdirl,1,jx,1,iy,'clm:c2ralbdirl')
    call getmem2d(c2ralbdifs,1,jx,1,iy,'clm:c2ralbdifs')
    call getmem2d(c2ralbdifl,1,jx,1,iy,'clm:c2ralbdifl')
    call getmem2d(c2rtaux,1,jx,1,iy,'clm:c2rtaux')
    call getmem2d(c2rtauy,1,jx,1,iy,'clm:c2rtauy')
    call getmem2d(c2ruvdrag,1,jx,1,iy,'clm:c2ruvdrag')
    call getmem2d(c2rlsmask,1,jx,1,iy,'clm:c2rlsmask')
    call getmem2d(c2rtgbb,1,jx,1,iy,'clm:c2rtgbb')
    call getmem2d(c2rsnowc,1,jx,1,iy,'clm:c2rsnowc')
    call getmem2d(c2rtest,1,jx,1,iy,'clm:c2rtest')
    call getmem2d(c2r2mt,1,jx,1,iy,'clm:c2r2mt')
    call getmem2d(c2r2mq,1,jx,1,iy,'clm:c2r2mq')
    call getmem2d(c2rtlef,1,jx,1,iy,'clm:c2rtlef')
    call getmem2d(c2ru10,1,jx,1,iy,'clm:c2ru10')
    call getmem2d(c2rsm10cm,1,jx,1,iy,'clm:c2rsm10cm')
    call getmem2d(c2rsm1m,1,jx,1,iy,'clm:c2rsm1m')
    call getmem2d(c2rsmtot,1,jx,1,iy,'clm:c2rsmtot')
    call getmem2d(c2rinfl,1,jx,1,iy,'clm:c2rinfl')
    call getmem2d(c2rro_sur,1,jx,1,iy,'clm:c2rro_sur')
    call getmem2d(c2rro_sub,1,jx,1,iy,'clm:c2rro_sub')
    call getmem2d(c2rfracsno,1,jx,1,iy,'clm:c2rfracsno')
    call getmem2d(c2rfvegnosno,1,jx,1,iy,'clm:c2rfvegnosno')
    call getmem2d(c2rprocmap,1,jx,1,iy,'clm:c2rprocmap')
#if (defined VOC)
    call getmem2d(voc_em,1,jx,1,iy,'clm:voc_em')
#endif
    if ( igases == 1 ) then
      call getmem3d(dep_vels,1,jx,1,iy,1,n_tr,'clm:dep_vels')
    end if

    call getmem1d(c2rngc,1,nproc,'clm:c2rngc')
    call getmem1d(c2rdisps,1,nproc,'clm:c2rdisps')

    call getmem2d(sols2d,jci1,jci2,ici1,ici2,'clm:sols2d')
    call getmem2d(soll2d,jci1,jci2,ici1,ici2,'clm:soll2d')
    call getmem2d(solsd2d,jci1,jci2,ici1,ici2,'clm:solsd2d')
    call getmem2d(solld2d,jci1,jci2,ici1,ici2,'clm:solld2d')
    call getmem2d(aldirs2d,jci1,jci2,ici1,ici2,'clm:aldirs2d')
    call getmem2d(aldirl2d,jci1,jci2,ici1,ici2,'clm:aldirl2d')
    call getmem2d(aldifs2d,jci1,jci2,ici1,ici2,'clm:aldifs2d')
    call getmem2d(aldifl2d,jci1,jci2,ici1,ici2,'clm:aldifl2d')
    call getmem2d(rs2d,jci1,jci2,ici1,ici2,'clm:rs2d')
    call getmem2d(ra2d,jci1,jci2,ici1,ici2,'clm:ra2d')
    call getmem2d(q2d,jci1,jci2,ici1,ici2,'clm:q2d')
    call getmem2d(lndcat2d,jci1,jci2,ici1,ici2,'clm:lndcat2d')
 
  end subroutine allocate_mod_clm
!
end module mod_clm
!
#endif
