      module mod_bdycod

      use mod_regcm_param

      implicit none
!
! COMMON /BCVARS/
!
#ifdef MPP1
      real(8) , dimension(ix,0:jxp+1) :: ps0 , ps1
      real(8) , dimension(ix,kx,jxp) :: qb0 , qb1 , so0 , so1 , tb0 ,   &
                                      & tb1 , ub0 , ub1 , vb0 , vb1
      real(8) , dimension(ix,jxp) :: ts0 , ts1
#else
      real(8) , dimension(ix,jx) :: ps0 , ps1 , ts0 , ts1
      real(8) , dimension(ix,kx,jx) :: qb0 , qb1 , so0 , so1 , tb0 ,    &
                                     & tb1 , ub0 , ub1 , vb0 , vb1
#endif
!
! COMMON /BDYCOD/
!
#ifdef MPP1
      real(8) , dimension(ix,0:jxp+1) :: peb , pebt , pwb , pwbt
      real(8) , dimension(nspgx,0:jxp+1) :: pnb , pnbt , psbt , pss
      real(8) , dimension(ix,kx,0:jxp+1) :: qeb , qebt , qwb , qwbt ,   &
           & teb , tebt , twb , twbt , ueb , uebt , uwb , uwbt , veb ,  &
           & vebt , vwb , vwbt
      real(8) , dimension(nspgx,kx,0:jxp+1) :: qnb , qnbt , qsb , qsbt ,&
           & tnb , tnbt , tsb , tsbt
      real(8) , dimension(kx,0:jxp+1) :: ui1 , ui2 , uil , uilx , vi1 , &
           & vi2 , vil , vilx
      real(8) , dimension(ix,kx) :: uj1 , uj2 , ujl , ujlx , vj1 , vj2 ,&
                                  & vjl , vjlx
      real(8) , dimension(nspgd,kx,0:jxp+1) :: unb , unbt , usb , usbt ,&
           & vnb , vnbt , vsb , vsbt
#else
      real(8) , dimension(ix,nspgx) :: peb , pebt , pwb , pwbt
      real(8) , dimension(nspgx,jx) :: pnb , pnbt , psbt , pss
      real(8) , dimension(ix,kx,nspgx) :: qeb , qebt , qwb , qwbt ,     &
           & teb , tebt , twb , twbt
      real(8) , dimension(nspgx,kx,jx) :: qnb , qnbt , qsb , qsbt ,     &
           & tnb , tnbt , tsb , tsbt
      real(8) , dimension(ix,kx,nspgd) :: ueb , uebt , uwb , uwbt ,     &
           & veb , vebt , vwb , vwbt
      real(8) , dimension(kx,jx) :: ui1 , ui2 , uil , uilx , vi1 , vi2 ,&
                                  & vil , vilx
      real(8) , dimension(ix,kx) :: uj1 , uj2 , ujl , ujlx , vj1 , vj2 ,&
                                  & vjl , vjlx
      real(8) , dimension(nspgd,kx,jx) :: unb , unbt , usb , usbt ,     &
           & vnb , vnbt , vsb , vsbt
#endif

#ifdef MPP1
!
! COMMON /BCVARSIO/
!
      real(8) , dimension(ix,jx) :: ps0_io , ps1_io , ts0_io , ts1_io
      real(8) , dimension(ix,kx,jx) :: qb0_io , qb1_io , so0_io ,      &
                                      & so1_io , tb0_io , tb1_io ,      &
                                      & ub0_io , ub1_io , vb0_io ,      &
                                      & vb1_io
!
! COMMON /BDYCODIO/
!
      real(8) , dimension(kx,jx) :: ui1_io , ui2_io , uilx_io ,        &
                                   & uil_io , vi1_io , vi2_io ,         &
                                   & vilx_io , vil_io
#endif

      end module mod_bdycod
