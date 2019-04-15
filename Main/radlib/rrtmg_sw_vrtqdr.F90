!     path:      $Source: /storm/rc1/cvsroot/rc/rrtmg_sw/src/rrtmg_sw_vrtqdr.f90,v $
!     author:    $Author: mike $
!     revision:  $Revision: 1.4 $
!     created:   $Date: 2009/05/22 22:22:22 $
!
      module rrtmg_sw_vrtqdr

!  --------------------------------------------------------------------------
! |                                                                          |
! |  Copyright 2002-2009, Atmospheric & Environmental Research, Inc. (AER).  |
! |  This software may be used, copied, or redistributed as long as it is    |
! |  not sold and this copyright notice is reproduced on each copy made.     |
! |  This model is provided as is without any express or implied warranties. |
! |                       (http://www.rtweb.aer.com/)                        |
! |                                                                          |
!  --------------------------------------------------------------------------

! ------- Modules -------

      use parkind, only: im => kind_im, rb => kind_rb , almostzero
!      use parrrsw, only: ngptsw

      implicit none

      contains

! --------------------------------------------------------------------------
      subroutine vrtqdr_sw(klev, kw, &
                           pref, prefd, ptra, ptrad, &
                           pdbt, prdnd, prup, prupd, ptdbt, &
                           pfd, pfu)
! --------------------------------------------------------------------------

! Purpose: This routine performs the vertical quadrature integration
!
! Interface:  *vrtqdr_sw* is called from *spcvrt_sw* and *spcvmc_sw*
!
! Modifications.
!
! Original: H. Barker
! Revision: Integrated with rrtmg_sw, J.-J. Morcrette, ECMWF, Oct 2002
! Revision: Reformatted for consistency with rrtmg_lw: MJIacono, AER, Jul 2006
!
!-----------------------------------------------------------------------

! ------- Declarations -------

! Input

      integer(kind=im), intent (in) :: klev                   ! number of model layers
      integer(kind=im), intent (in) :: kw                     ! g-point index

      real(kind=rb), intent(in) :: pref(:)                    ! direct beam reflectivity
                                                              !   Dimensions: (nlayers+1)
      real(kind=rb), intent(in) :: prefd(:)                   ! diffuse beam reflectivity
                                                              !   Dimensions: (nlayers+1)
      real(kind=rb), intent(in) :: ptra(:)                    ! direct beam transmissivity
                                                              !   Dimensions: (nlayers+1)
      real(kind=rb), intent(in) :: ptrad(:)                   ! diffuse beam transmissivity
                                                              !   Dimensions: (nlayers+1)

      real(kind=rb), intent(in) :: pdbt(:)
                                                              !   Dimensions: (nlayers+1)
      real(kind=rb), intent(in) :: ptdbt(:)
                                                              !   Dimensions: (nlayers+1)

      real(kind=rb), intent(inout) :: prdnd(:)
                                                              !   Dimensions: (nlayers+1)
      real(kind=rb), intent(inout) :: prup(:)
                                                              !   Dimensions: (nlayers+1)
      real(kind=rb), intent(inout) :: prupd(:)
                                                              !   Dimensions: (nlayers+1)

! Output
      real(kind=rb), intent(out) :: pfd(:,:)                  ! downwelling flux (W/m2)
                                                              !   Dimensions: (nlayers+1,ngptsw)
                                                              ! unadjusted for earth/sun distance or zenith angle
      real(kind=rb), intent(out) :: pfu(:,:)                  ! upwelling flux (W/m2)
                                                              !   Dimensions: (nlayers+1,ngptsw)
                                                              ! unadjusted for earth/sun distance or zenith angle

! Local

      integer(kind=im) :: ikp, ikx, jk

      real(kind=rb) :: zreflect , a1 , a2
      real(kind=rb) :: ztdn(klev+1)

! Definitions
!
! pref(jk)   direct reflectance
! prefd(jk)  diffuse reflectance
! ptra(jk)   direct transmittance
! ptrad(jk)  diffuse transmittance
!
! pdbt(jk)   layer mean direct beam transmittance
! ptdbt(jk)  total direct beam transmittance at levels
!
!-----------------------------------------------------------------------------

! Link lowest layer with surface

      zreflect = 1._rb / (1._rb - prefd(klev+1) * prefd(klev))
      prup(klev) = pref(klev) + (ptrad(klev) * &
                 ((ptra(klev) - pdbt(klev)) * prefd(klev+1) + &
                   pdbt(klev) * pref(klev+1))) * zreflect
      prupd(klev) = prefd(klev) + ptrad(klev) * ptrad(klev) * &
                    prefd(klev+1) * zreflect

! Pass from bottom to top

      do jk = 1,klev-1
         ikp = klev+1-jk
         ikx = ikp-1
         zreflect = 1._rb / (1._rb -prupd(ikp) * prefd(ikx))
         prup(ikx) = pref(ikx) + (ptrad(ikx) * &
                   ((ptra(ikx) - pdbt(ikx)) * prupd(ikp) + &
                     pdbt(ikx) * prup(ikp))) * zreflect
         prupd(ikx) = prefd(ikx) + ptrad(ikx) * ptrad(ikx) * &
                      prupd(ikp) * zreflect
      enddo

! Upper boundary conditions

      ztdn(1) = 1._rb
      prdnd(1) = 0._rb
      ztdn(2) = ptra(1)
      prdnd(2) = prefd(1)

! Pass from top to bottom

      do jk = 2,klev
         ikp = jk+1
         if ( ptra(jk) > 1.0e-5_rb ) then
            ztdn(ikp) = ptdbt(jk) * ptra(jk)
         else
            ztdn(ikp) = 0.0_rb
         end if
         zreflect = 1._rb / (1._rb - prefd(jk) * prdnd(jk))
         if ( zreflect > almostzero .and. ptrad(jk) > almostzero ) then
           a1 = ztdn(jk) - ptdbt(jk)
           if ( ptdbt(jk) > almostzero .and. &
                pref(jk)  > almostzero .and. &
                prdnd(jk) > almostzero ) then
              a2 = ptdbt(jk) * pref(jk) * prdnd(jk)
           else
              a2 = 0.0_rb
           end if
           if ( a1 < almostzero ) a1 = 0.0_rb
           if ( a2 < almostzero ) a2 = 0.0_rb
           ztdn(ikp) = ztdn(ikp) + zreflect * ptrad(jk) * (a1 + a2)
         end if
         ! ztdn(ikp) = ptdbt(jk) * ptra(jk) + &
         !            (ptrad(jk) * ((ztdn(jk) - ptdbt(jk)) + &
         !             ptdbt(jk) * pref(jk) * prdnd(jk))) * zreflect
         prdnd(ikp) = prefd(jk) + ptrad(jk) * ptrad(jk) * &
                      prdnd(jk) * zreflect
         if ( prdnd(ikp) < almostzero ) prdnd(ikp) = 0.0_rb
      enddo

! Up and down-welling fluxes at levels

      do jk = 1,klev+1
         zreflect = 1._rb / (1._rb - prdnd(jk) * prupd(jk))
         if ( zreflect > almostzero ) then
            if ( ptdbt(jk) > almostzero .and. &
                 prup(jk)  > almostzero ) then
               a1 = ptdbt(jk) * prup(jk)
            else
               a1 = 0.0_rb
            end if
            if ( prupd(jk) > almostzero ) then
               a2 = (ztdn(jk) - ptdbt(jk)) * prupd(jk)
            else
               a2 = 0.0_rb
            end if
            if ( a1 < almostzero ) a1 = 0.0_rb
            if ( a2 < almostzero ) a2 = 0.0_rb
            pfu(jk,kw) = (a1 + a2) * zreflect
            a1 = ztdn(jk) - ptdbt(jk)
            if ( ptdbt(jk) > almostzero .and. &
                 prup(jk)  > almostzero .and. &
                 prdnd(jk) > almostzero ) then
               a2 = ptdbt(jk) * prup(jk) * prdnd(jk)
            else
               a2 = 0.0_rb
            end if
            if ( a1 < almostzero ) a1 = 0.0_rb
            if ( a2 < almostzero ) a2 = 0.0_rb
            pfd(jk,kw) = ptdbt(jk) + (a1 + a2) * zreflect
         else
            pfu(jk,kw) = 0.0_rb
            pfd(jk,kw) = ptdbt(jk)
         end if
      enddo

      end subroutine vrtqdr_sw

      end module rrtmg_sw_vrtqdr
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
