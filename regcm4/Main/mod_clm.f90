      module mod_clm

      implicit none

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

      real(8) :: r2ctb_all   
      real(8) :: r2cqb_all    
      real(8) :: r2czga_all   
      real(8) :: r2cpsb_all   
      real(8) :: r2cuxb_all   
      real(8) :: r2cvxb_all   
      real(8) :: r2crnc_all   
      real(8) :: r2crnnc_all  
      real(8) :: r2csols_all  
      real(8) :: r2csoll_all  
      real(8) :: r2csolsd_all 
      real(8) :: r2csolld_all
      real(8) :: r2cflwd_all
      real(8) :: r2ccosz_all
      real(8) :: r2cxlat_all     ! xlat in radians
      real(8) :: r2cxlon_all     ! xlon in radians
      real(8) :: r2cxlatd_all    ! xlat in degrees
      real(8) :: r2cxlond_all    ! xlon in degrees
      real(8) :: c2rtgb
      real(8) :: c2rsenht
      real(8) :: c2rlatht
      real(8) :: c2ralbdirs
      real(8) :: c2ralbdirl
      real(8) :: c2ralbdifs
      real(8) :: c2ralbdifl
      real(8) :: c2rtaux
      real(8) :: c2rtauy
      real(8) :: c2ruvdrag
      real(8) :: c2rlsmask
      real(8) :: c2rtgbb
      real(8) :: c2rcosz
      real(8) :: c2rsnowc
      real(8) :: c2rtest
      real(8) :: c2r2mt
      real(8) :: c2r2mq
      real(8) :: c2rtlef
      real(8) :: c2ru10
      real(8) :: c2rsm10cm
      real(8) :: c2rsm1m
      real(8) :: c2rsmtot
      real(8) :: c2rinfl
      real(8) :: c2rro_sur
      real(8) :: c2rro_sub
      real(8) :: c2rfracsno      
      real(8) :: c2rfvegnosno 

      integer :: c2rprocmap
      integer :: c2rngc
      integer :: c2rdisps
      integer :: c2rcnts

      end module mod_clm
