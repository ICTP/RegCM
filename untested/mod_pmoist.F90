      module pmoist

      use regcm_param

      implicit none
!
! COMMON /PMOIST/
!
      real(8) :: avt , bvt , caccr , cevap , clfrcv , clfrcvmax ,       &
               & cllwcv , conf , dtauc , edtmax , edtmaxo , edtmaxx ,   &
               & edtmin , edtmino , edtminx , ep1 , ep2 , fcmax , g3pb ,&
               & g4pb , g5pb , gulland , guloce , htmax , htmin ,       &
               & mincld , n0r , pbcmax , ppi , prac , prec1 , prec2 ,   &
               & qck10 , qck1land , qck1oce , qcth , qdcrit , rh0land , &
               & rh0oce , rhmax , rv , shrmax , shrmin , skbmax , svp1 ,&
               & svp2 , svp3 , tc0 , trel , vtc , xlv , xlvocp
#ifdef MPP1
      real(8) , dimension(ix,jxp) :: cbmf2d , cgul , dtauc2d ,          &
                                   & edtmax2d , edtmaxo2d , edtmaxx2d , &
                                   & edtmin2d , edtmino2d , edtminx2d , &
                                   & htmax2d , htmin2d , mincld2d ,     &
                                   & pbcmax2d , qck1 , rh0 , shrmax2d , &
                                   & shrmin2d
      real(8) , dimension(ix,kx,jxp) :: fcc , rsheat , rswat
#else
      real(8) , dimension(ix,jx) :: cbmf2d , cgul , dtauc2d , edtmax2d ,&
                                  & edtmaxo2d , edtmaxx2d , edtmin2d ,  &
                                  & edtmino2d , edtminx2d , htmax2d ,   &
                                  & htmin2d , mincld2d , pbcmax2d ,     &
                                  & qck1 , rh0 , shrmax2d , shrmin2d
      real(8) , dimension(ix,kx,jx) :: fcc , rsheat , rswat
#endif
      real(8) , dimension(kx) :: qwght
      real(8) , dimension(kx,5:kx,1:kx-3) :: twght , vqflx
!
! COMMON /PMOISTINT/
!
#ifdef MPP1
      integer , dimension(jxp) :: icon
      integer , dimension(ix,jxp) :: kbmax2d
#else
      integer , dimension(jx) :: icon
      integer , dimension(ix,jx) :: kbmax2d
#endif
      integer :: kbmax

#ifdef MPP1
!
! COMMON /PMOISTIO/
!
      real(8) , dimension(ix,mjx) :: cbmf2d_io
      real(8) , dimension(ix,kx,mjx) :: fcc_io , rsheat_io , rswat_io
#endif

      end module pmoist
