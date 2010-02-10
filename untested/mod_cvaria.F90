      module mod_cvaria

      use mod_regcm_param

      implicit none
!
! COMMON /CVARIA/
!
#ifdef MPP1
      real(8) , dimension(ix,kx,jxp) :: diffq , difft , difuu , difuv , &
                                      & omega , qcc , qcten , qvc ,     &
                                      & qvten , tc , td , tten , uc ,   &
                                      & uten , vc , vten , xkc
      real(8) , dimension(ix,jxp) :: pdota , psc , pten
      real(8) , dimension(ix,kx,0:jxp) :: phi
      real(8) , dimension(ix,0:jxp+1) :: psd
      real(8) , dimension(ix,kx,0:jxp+1) :: qc , qv , t , u , v
      real(8) , dimension(ix,kxp1,0:jxp+1) :: qdot
      real(8) , dimension(ix,kx) :: qvcs
#else
      real(8) , dimension(ix,kx,jx) :: diffq , difft , difuu , difuv ,  &
                                     & omega , phi , qc , qcc , qcten , &
                                     & qv , qvc , qvten , t , tc , td , &
                                     & tten , u , uc , uten , v , vc ,  &
                                     & vten , xkc
      real(8) , dimension(ix,jx) :: pdota , psc , psd , pten
      real(8) , dimension(ix,kxp1,jx) :: qdot
      real(8) , dimension(ix,kx) :: qvcs
#endif
!
! COMMON /TRACER/
!
#ifdef MPP1
      real(8) , dimension(ix,kx,0:jxp+1,ntr) :: chi
      real(8) , dimension(ix,kx,jxp,ntr) :: chic , chiten
#else
      real(8) , dimension(ix,kx,jx,ntr) :: chi , chic , chiten
#endif
      end module mod_cvaria
