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

      use mod_dynparam

      implicit none

#ifdef MPP1

      integer :: imask

      integer :: r2cdtime      ! timestep in seconds
      integer :: r2cnsrest     ! 0=initial, 1=restart
      integer :: r2cnestep     ! final timestep (or day if negative) number
      integer :: r2cnelapse    ! # of timesteps (or days if negative) 
                               ! to extend a run
      integer :: r2cnstep      ! current timestep (ktau)

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
      real(8) :: r2cdx         ! model resolution (m)
      real(8) :: r2carea       ! area of each grid cell (km^2)
      real(8) :: r2cedgen      ! N edge of grid
      real(8) :: r2cedges      ! S edge of grid
      real(8) :: r2cedgee      ! E edge of grid
      real(8) :: r2cedgew      ! W edge of grid

      real(8) :: r2ceccen      ! orbital eccentricity
      real(8) :: r2cobliqr     ! Earths obliquity in rad
      real(8) :: r2clambm0     ! Mean long of perihelion at
                               ! vernal equinox (radians)
      real(8) :: r2cmvelpp     ! moving vernal equinox long
      real(8) :: r2ceccf       ! Earth-sun distance factor (1/r)**2

      logical :: r2cdoalb      ! time for next albedo call

      integer :: c2rcnts

      real(8) , target , allocatable , dimension(:,:,:) :: clmspace
      private :: clmspace

      real(8) , pointer , dimension(:,:) :: r2ctb_all   
      real(8) , pointer , dimension(:,:) :: r2cqb_all    
      real(8) , pointer , dimension(:,:) :: r2czga_all   
      real(8) , pointer , dimension(:,:) :: r2cpsb_all   
      real(8) , pointer , dimension(:,:) :: r2cuxb_all   
      real(8) , pointer , dimension(:,:) :: r2cvxb_all   
      real(8) , pointer , dimension(:,:) :: r2crnc_all   
      real(8) , pointer , dimension(:,:) :: r2crnnc_all  
      real(8) , pointer , dimension(:,:) :: r2csols_all  
      real(8) , pointer , dimension(:,:) :: r2csoll_all  
      real(8) , pointer , dimension(:,:) :: r2csolsd_all 
      real(8) , pointer , dimension(:,:) :: r2csolld_all
      real(8) , pointer , dimension(:,:) :: r2cflwd_all
      real(8) , pointer , dimension(:,:) :: r2ccosz_all
      real(8) , pointer , dimension(:,:) :: r2cxlat_all     ! xlat in radians
      real(8) , pointer , dimension(:,:) :: r2cxlon_all     ! xlon in radians
      real(8) , pointer , dimension(:,:) :: r2cxlatd_all    ! xlat in degrees
      real(8) , pointer , dimension(:,:) :: r2cxlond_all    ! xlon in degrees
      real(8) , pointer , dimension(:,:) :: c2rtgb
      real(8) , pointer , dimension(:,:) :: c2rsenht
      real(8) , pointer , dimension(:,:) :: c2rlatht
      real(8) , pointer , dimension(:,:) :: c2ralbdirs
      real(8) , pointer , dimension(:,:) :: c2ralbdirl
      real(8) , pointer , dimension(:,:) :: c2ralbdifs
      real(8) , pointer , dimension(:,:) :: c2ralbdifl
      real(8) , pointer , dimension(:,:) :: c2rtaux
      real(8) , pointer , dimension(:,:) :: c2rtauy
      real(8) , pointer , dimension(:,:) :: c2ruvdrag
      real(8) , pointer , dimension(:,:) :: c2rlsmask
      real(8) , pointer , dimension(:,:) :: c2rtgbb
      real(8) , pointer , dimension(:,:) :: c2rcosz
      real(8) , pointer , dimension(:,:) :: c2rsnowc
      real(8) , pointer , dimension(:,:) :: c2rtest
      real(8) , pointer , dimension(:,:) :: c2r2mt
      real(8) , pointer , dimension(:,:) :: c2r2mq
      real(8) , pointer , dimension(:,:) :: c2rtlef
      real(8) , pointer , dimension(:,:) :: c2ru10
      real(8) , pointer , dimension(:,:) :: c2rsm10cm
      real(8) , pointer , dimension(:,:) :: c2rsm1m
      real(8) , pointer , dimension(:,:) :: c2rsmtot
      real(8) , pointer , dimension(:,:) :: c2rinfl
      real(8) , pointer , dimension(:,:) :: c2rro_sur
      real(8) , pointer , dimension(:,:) :: c2rro_sub
      real(8) , pointer , dimension(:,:) :: c2rfracsno      
      real(8) , pointer , dimension(:,:) :: c2rfvegnosno 
      integer , allocatable , dimension(:,:) :: c2rprocmap
      integer , allocatable , dimension(:) :: c2rngc
      integer , allocatable , dimension(:) :: c2rdisps

#endif

      contains

      subroutine allocate_mod_clm
      implicit none
!     About the dimension ordering:
!     regcm: ix=lat,jx=lon, arrays are lat by lon
!     clm: i=lon, j=lat, arrays are lon by lat
#ifdef MPP1
      allocate(clmspace(jx,iy,45))
      r2ctb_all     => clmspace(:,:,1)
      r2cqb_all     => clmspace(:,:,2)
      r2czga_all    => clmspace(:,:,3)
      r2cpsb_all    => clmspace(:,:,4)
      r2cuxb_all    => clmspace(:,:,5)
      r2cvxb_all    => clmspace(:,:,6)
      r2crnc_all    => clmspace(:,:,7)
      r2crnnc_all   => clmspace(:,:,8)
      r2csols_all   => clmspace(:,:,9)
      r2csoll_all   => clmspace(:,:,10)
      r2csolsd_all  => clmspace(:,:,11)
      r2csolld_all  => clmspace(:,:,12)
      r2cflwd_all   => clmspace(:,:,13)
      r2ccosz_all   => clmspace(:,:,14)
      r2cxlat_all   => clmspace(:,:,15)
      r2cxlon_all   => clmspace(:,:,16)
      r2cxlatd_all  => clmspace(:,:,17)
      r2cxlond_all  => clmspace(:,:,18)
      c2rtgb        => clmspace(:,:,19)
      c2rsenht      => clmspace(:,:,20)
      c2rlatht      => clmspace(:,:,21)
      c2ralbdirs    => clmspace(:,:,22)
      c2ralbdirl    => clmspace(:,:,23)
      c2ralbdifs    => clmspace(:,:,24)
      c2ralbdifl    => clmspace(:,:,25)
      c2rtaux       => clmspace(:,:,26)
      c2rtauy       => clmspace(:,:,27)
      c2ruvdrag     => clmspace(:,:,28)
      c2rlsmask     => clmspace(:,:,29)
      c2rtgbb       => clmspace(:,:,30)
      c2rcosz       => clmspace(:,:,31)
      c2rsnowc      => clmspace(:,:,32)
      c2rtest       => clmspace(:,:,33)
      c2r2mt        => clmspace(:,:,34)
      c2r2mq        => clmspace(:,:,35)
      c2rtlef       => clmspace(:,:,36)
      c2ru10        => clmspace(:,:,37)
      c2rsm10cm     => clmspace(:,:,38)
      c2rsm1m       => clmspace(:,:,39)
      c2rsmtot      => clmspace(:,:,40)
      c2rinfl       => clmspace(:,:,41)
      c2rro_sur     => clmspace(:,:,42)
      c2rro_sub     => clmspace(:,:,43)
      c2rfracsno    => clmspace(:,:,44)
      c2rfvegnosno  => clmspace(:,:,45)
      allocate(c2rprocmap(jx,iy))
      allocate(c2rngc(nproc))
      allocate(c2rdisps(nproc))
#endif
      end subroutine allocate_mod_clm

      end module mod_clm

#endif
