      module mod_mppio

      use regcm_param

      implicit none

#ifdef MPP1
!
! COMMON /BAT2D1IO/
!
      real(8) , dimension(nnsg,ilx,jxbb) :: col2d_io , dew2d_io ,       &
           & evpa2d_io , gwet2d_io , ircp2d_io , ocld2d_io , rno2d_io , &
           & rnos2d_io , sag2d_io , scv2d_io , sena2d_io , sice2d_io ,  &
           & srw2d_io , ssw2d_io , swt2d_io , taf2d_io , text2d_io ,    &
           & tg2d_io , tgb2d_io , tlef2d_io , veg2d1_io
      real(8) , dimension(nnsg,ix,mjx) :: ht1_io , satbrt1_io
!
! COMMON /BAT2DIO/
!
      real(8) , dimension(ilx,jxbb) :: flw2d_io , flwd2d_io , fsw2d_io ,&
                                     & sabv2d_io , sdelqk2d_io ,        &
                                     & sdeltk2d_io , sfracb2d_io ,      &
                                     & sfracs2d_io , sfracv2d_io ,      &
                                     & sinc2d_io , sol2d_io ,           &
                                     & solvd2d_io , solvs2d_io ,        &
                                     & ssw2da_io , svegfrac2d_io ,      &
                                     & veg2d_io

      
      real(kind=4) , dimension(mjx-2,ix-2,numbat) :: fbat_io
      real(4) , dimension(nnsg,mjx-2,ix-2,numsub) :: fsub_io
#endif

      end module mod_mppio
