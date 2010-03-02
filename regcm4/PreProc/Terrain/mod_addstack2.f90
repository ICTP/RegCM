      module mod_addstack2
      use mod_domain , only : iy , jx , nsg
      implicit none
      real(4) :: clong
      real(4) , dimension(iy,jx) :: corc , hscr1 , htsavc , sumc ,      &
                                  & wtmaxc
      real(4) , dimension(iy*nsg,jx*nsg) :: corc_s , hscr1_s ,          &
           & htsavc_s , sumc_s , wtmaxc_s
      real(4) , dimension(iy,jx,2) :: itex , land
      real(4) , dimension(iy*nsg,jx*nsg,2) :: itex_s , land_s
      integer , dimension(iy,jx) :: nsc
      integer , dimension(iy*nsg,jx*nsg) :: nsc_s
      end module mod_addstack2
