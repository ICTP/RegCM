      module mod_clm

      use mod_dynparam , only : iy , jx , nproc

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

      real(8) , dimension(jx,iy) :: r2ctb_all   
      real(8) , dimension(jx,iy) :: r2cqb_all    
      real(8) , dimension(jx,iy) :: r2czga_all   
      real(8) , dimension(jx,iy) :: r2cpsb_all   
      real(8) , dimension(jx,iy) :: r2cuxb_all   
      real(8) , dimension(jx,iy) :: r2cvxb_all   
      real(8) , dimension(jx,iy) :: r2crnc_all   
      real(8) , dimension(jx,iy) :: r2crnnc_all  
      real(8) , dimension(jx,iy) :: r2csols_all  
      real(8) , dimension(jx,iy) :: r2csoll_all  
      real(8) , dimension(jx,iy) :: r2csolsd_all 
      real(8) , dimension(jx,iy) :: r2csolld_all
      real(8) , dimension(jx,iy) :: r2cflwd_all
      real(8) , dimension(jx,iy) :: r2ccosz_all
      real(8) , dimension(jx,iy) :: r2cxlat_all     ! xlat in radians
      real(8) , dimension(jx,iy) :: r2cxlon_all     ! xlon in radians
      real(8) , dimension(jx,iy) :: r2cxlatd_all    ! xlat in degrees
      real(8) , dimension(jx,iy) :: r2cxlond_all    ! xlon in degrees
      real(8) , dimension(jx,iy) :: c2rtgb
      real(8) , dimension(jx,iy) :: c2rsenht
      real(8) , dimension(jx,iy) :: c2rlatht
      real(8) , dimension(jx,iy) :: c2ralbdirs
      real(8) , dimension(jx,iy) :: c2ralbdirl
      real(8) , dimension(jx,iy) :: c2ralbdifs
      real(8) , dimension(jx,iy) :: c2ralbdifl
      real(8) , dimension(jx,iy) :: c2rtaux
      real(8) , dimension(jx,iy) :: c2rtauy
      real(8) , dimension(jx,iy) :: c2ruvdrag
      real(8) , dimension(jx,iy) :: c2rlsmask
      real(8) , dimension(jx,iy) :: c2rtgbb
      real(8) , dimension(jx,iy) :: c2rcosz
      real(8) , dimension(jx,iy) :: c2rsnowc
      real(8) , dimension(jx,iy) :: c2rtest
      real(8) , dimension(jx,iy) :: c2r2mt
      real(8) , dimension(jx,iy) :: c2r2mq
      real(8) , dimension(jx,iy) :: c2rtlef
      real(8) , dimension(jx,iy) :: c2ru10
      real(8) , dimension(jx,iy) :: c2rsm10cm
      real(8) , dimension(jx,iy) :: c2rsm1m
      real(8) , dimension(jx,iy) :: c2rsmtot
      real(8) , dimension(jx,iy) :: c2rinfl
      real(8) , dimension(jx,iy) :: c2rro_sur
      real(8) , dimension(jx,iy) :: c2rro_sub
      real(8) , dimension(jx,iy) :: c2rfracsno      
      real(8) , dimension(jx,iy) :: c2rfvegnosno 

      integer :: c2rprocmap(jx,iy)
      integer :: c2rngc(nproc)
      integer :: c2rdisps(nproc)
      integer :: c2rcnts

      end module mod_clm
