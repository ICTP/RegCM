      module aerm

      use regcm_param
      use parrad

      implicit none

#ifdef MPP1
      real(8) , dimension(ix-1,kx,jxp) :: aermm
#else
      real(8) , dimension(ix-1,kx,jx-1) :: aermm
#endif

!     Background aerosol mass mixing ratio
      real(8) , dimension(plon,plevr) :: aermmb
!     Radiation level aerosol mass mixing ratio
      real(8) , dimension(plond,plevr,ntr) :: aermmr

      end module aerm
