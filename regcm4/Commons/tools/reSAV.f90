      program resav
      use mod_regcm_param , only : ix => iy
      use mod_regcm_param , only : jx => jx
      use mod_regcm_param , only : kx => kz
      implicit none
!
! PARAMETER definitions
!
      integer , parameter :: ixb = ix - 1 , jxb = jx - 1 ,              &
                           & ilx = ix - 1 , jlx = jx - 1 ,              &
                           & kxp1 = kx + 1 , nsplit = 2
!
! Local variables
!
      real(8) , dimension(ixb,kx,4,jxb) :: absnxt
      real(8) , dimension(ixb,kx,kx,jxb) :: abstot
      real(8) , dimension(ixb,jxb) :: col2d , dew2d , flw2d ,      &
                & flwd2d , fsw2d , gwet2d , ircp2d , ocld2d , pptc ,    &
                & pptnc , prca2d , prnca2d , sabv2d , sag2d , scv2d ,   &
                & sice2d , sinc2d , sol2d , solvd2d , solvs2d , srw2d , &
                & ssw2d , swt2d , taf2d , text2d , tgb2d , tlef2d ,     &
                & veg2d
      real(8) , dimension(ix,jx,nsplit) :: dstor , hstor
      real(8) , dimension(ixb,kx,jxb) :: emstot
      real(8) , dimension(ix,jx) :: f , hfx , ht , msfd , msfx ,   &
                & psa , psb , qfx , rainc , rainnc , rh0 , satbrt ,     &
                & snowc , tga , tgb , tgbb , uvdrag , xlat , xlong
      real(8) , dimension(ix,jx,kx) :: fcc , qca , qcb , qva ,     &
                & qvb , rsheat , rswat , ta , tb , ua , ub , va , vb
      real(8) , dimension(ilx,jlx,kx) :: heatrt
      integer :: iutl , ktau , ldatez , lday , lhour , lmonth , lyear , &
               & mdate0
      integer*8 :: ntime
      real(8) , dimension(ilx,jlx,kxp1) :: o3prof
      real(4) , dimension(ix,jx) :: ps0 , ts0
      real(4) , dimension(ix,jx,kx) :: qb0 , tb0 , ub0 , vb0
      real(8) :: tdadv , tdini , tqadv , tqeva , tqini , tqrai ,   &
                    & xtime
      real(8) , dimension(jx,kx) :: ui1 , ui2 , uil , uilx , vi1 , &
                & vi2 , vil , vilx
      real(8) , dimension(ix,kx) :: uj1 , uj2 , ujl , ujlx , vj1 , &
                & vj2 , vjl , vjlx
!
      iutl = 14
      read (iutl) mdate0
      read (iutl) ktau , xtime , ldatez , lyear , lmonth , lday ,       &
                & lhour , ntime
      read (iutl) ub0 , vb0 , qb0 , tb0 , ps0 , ts0
      read (iutl) ua
      read (iutl) ub
      read (iutl) va
      read (iutl) vb
      read (iutl) ta
      read (iutl) tb
      read (iutl) qva
      read (iutl) qvb
      read (iutl) qca
      read (iutl) qcb
      read (iutl) psa , psb , satbrt , f
      read (iutl) ht , msfx , msfd , xlat , xlong
      read (iutl) tga , tgb , rainc , rainnc
      read (iutl) hfx , qfx , snowc , uvdrag
      read (iutl) tdini , tdadv , tqini , tqadv , tqeva , tqrai
      read (iutl) absnxt , abstot , emstot
      read (iutl) fcc , rh0
      read (iutl) sol2d , solvd2d , solvs2d , flw2d , flwd2d , fsw2d ,  &
                & sabv2d , sinc2d
      read (iutl) taf2d , tlef2d , tgbb , ssw2d , srw2d , tgb2d ,       &
                & swt2d , scv2d , gwet2d , veg2d , sag2d , sice2d ,     &
                & dew2d , ircp2d , text2d , col2d , ocld2d , rsheat ,   &
                & rswat , heatrt , o3prof
      read (iutl) pptnc , pptc , prca2d , prnca2d
      read (iutl) dstor
      read (iutl) hstor
      read (iutl) uj1 , uj2 , ujlx , ujl
      read (iutl) ui1 , ui2 , uilx , uil
      read (iutl) vj1 , vj2 , vjlx , vjl
      read (iutl) vi1 , vi2 , vilx , vil
!
      end program resav
