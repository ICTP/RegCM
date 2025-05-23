!     path:      $Source$
!     author:    $Author$
!     revision:  $Revision$
!     created:   $Date$
!
      module rrtmg_lw_cldprmc

!----------------------------------------------------------------------------
! Copyright (c) 2002-2020, Atmospheric & Environmental Research, Inc. (AER)
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!  * Redistributions of source code must retain the above copyright
!    notice, this list of conditions and the following disclaimer.
!  * Redistributions in binary form must reproduce the above copyright
!    notice, this list of conditions and the following disclaimer in the
!    documentation and/or other materials provided with the distribution.
!  * Neither the name of Atmospheric & Environmental Research, Inc., nor
!    the names of its contributors may be used to endorse or promote products
!    derived from this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
! ARE DISCLAIMED. IN NO EVENT SHALL ATMOSPHERIC & ENVIRONMENTAL RESEARCH, INC.,
! BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
! CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
! SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
! INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
! CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
! ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
! THE POSSIBILITY OF SUCH DAMAGE.
!                        (http://www.rtweb.aer.com/)
!----------------------------------------------------------------------------

! --------- Modules ----------

      use parkind, only : im => kind_im, rb => kind_rb
      use parrrtm, only : ngptlw, nbndlw
      use rrlw_cld, only: abscld1, absliq0, absliq1, &
                          absice0, absice1, absice2, absice3
      use rrlw_wvn, only: ngb
      use rrlw_vsn, only: hvrclc, hnamclc

      implicit none

      contains

! ------------------------------------------------------------------------------
      subroutine cldprmc(nlayers, inflag, iceflag, liqflag, cldfmc, &
                         ciwpmc, clwpmc, reicmc, relqmc, ncbands, taucmc)
! ------------------------------------------------------------------------------

! Purpose:  Compute the cloud optical depth(s) for each cloudy layer.

! ------- Input -------

      integer(kind=im), intent(in) :: nlayers         ! total number of layers
      integer(kind=im), intent(in) :: inflag          ! see definitions
      integer(kind=im), intent(in) :: iceflag         ! see definitions
      integer(kind=im), intent(in) :: liqflag         ! see definitions

      real(kind=rb), intent(in) :: cldfmc(:,:)        ! cloud fraction [mcica]
                                                      !    Dimensions: (ngptlw,nlayers)
      real(kind=rb), intent(in) :: ciwpmc(:,:)        ! cloud ice water path [mcica]
                                                      !    Dimensions: (ngptlw,nlayers)
      real(kind=rb), intent(in) :: clwpmc(:,:)        ! cloud liquid water path [mcica]
                                                      !    Dimensions: (ngptlw,nlayers)
      real(kind=rb), intent(in) :: relqmc(:)          ! liquid particle effective radius (microns)
                                                      !    Dimensions: (nlayers)
      real(kind=rb), intent(in) :: reicmc(:)          ! ice particle effective radius (microns)
                                                      !    Dimensions: (nlayers)
                                                      ! specific definition of reicmc depends on setting of iceflag:
                                                      ! iceflag = 0: ice effective radius, r_ec, (Ebert and Curry, 1992),
                                                      !              r_ec must be >= 10.0 microns
                                                      ! iceflag = 1: ice effective radius, r_ec, (Ebert and Curry, 1992),
                                                      !              r_ec range is limited to 13.0 to 130.0 microns
                                                      ! iceflag = 2: ice effective radius, r_k, (Key, Streamer Ref. Manual, 1996)
                                                      !              r_k range is limited to 5.0 to 131.0 microns
                                                      ! iceflag = 3: generalized effective size, dge, (Fu, 1996),
                                                      !              dge range is limited to 5.0 to 140.0 microns
                                                      !              [dge = 1.0315 * r_ec]

! ------- Output -------

      integer(kind=im), intent(out) :: ncbands        ! number of cloud spectral bands
      real(kind=rb), intent(inout) :: taucmc(:,:)     ! cloud optical depth [mcica]
                                                      !    Dimensions: (ngptlw,nlayers)

! ------- Local -------

      integer(kind=im) :: lay                         ! Layer index
      integer(kind=im) :: ib                          ! spectral band index
      integer(kind=im) :: ig                          ! g-point interval index
      integer(kind=im) :: indx

      real(kind=rb) :: abscoice(ngptlw)               ! ice absorption coefficients
      real(kind=rb) :: abscoliq(ngptlw)               ! liquid absorption coefficients
      real(kind=rb) :: cwp                            ! cloud water path
      real(kind=rb) :: radice                         ! cloud ice effective size (microns)
      real(kind=rb) :: factor                         !
      real(kind=rb) :: fint                           !
      real(kind=rb) :: radliq                         ! cloud liquid droplet radius (microns)
      real(kind=rb), parameter :: cldmin = 1.e-20_rb  ! minimum value for cloud quantities

! ------- Definitions -------

!     Explanation of the method for each value of INFLAG.  Values of
!     0 or 1 for INFLAG do not distingish being liquid and ice clouds.
!     INFLAG = 2 does distinguish between liquid and ice clouds, and
!     requires further user input to specify the method to be used to
!     compute the aborption due to each.
!     INFLAG = 0:  For each cloudy layer, the cloud fraction and (gray)
!                  optical depth are input.
!     INFLAG = 1:  For each cloudy layer, the cloud fraction and cloud
!                  water path (g/m2) are input.  The (gray) cloud optical
!                  depth is computed as in CCM2.
!     INFLAG = 2:  For each cloudy layer, the cloud fraction, cloud
!                  water path (g/m2), and cloud ice fraction are input.
!       ICEFLAG = 0:  The ice effective radius (microns) is input and the
!                     optical depths due to ice clouds are computed as in CCM3.
!       ICEFLAG = 1:  The ice effective radius (microns) is input and the
!                     optical depths due to ice clouds are computed as in
!                     Ebert and Curry, JGR, 97, 3831-3836 (1992).  The
!                     spectral regions in this work have been matched with
!                     the spectral bands in RRTM to as great an extent
!                     as possible:
!                     E&C 1      IB = 5      RRTM bands 9-16
!                     E&C 2      IB = 4      RRTM bands 6-8
!                     E&C 3      IB = 3      RRTM bands 3-5
!                     E&C 4      IB = 2      RRTM band 2
!                     E&C 5      IB = 1      RRTM band 1
!       ICEFLAG = 2:  The ice effective radius (microns) is input and the
!                     optical properties due to ice clouds are computed from
!                     the optical properties stored in the RT code,
!                     STREAMER v3.0 (Reference: Key. J., Streamer
!                     User's Guide, Cooperative Institute for
!                     Meteorological Satellite Studies, 2001, 96 pp.).
!                     Valid range of values for re are between 5.0 and
!                     131.0 micron.
!       ICEFLAG = 3: The ice generalized effective size (dge) is input
!                    and the optical properties, are calculated as in
!                    Q. Fu, J. Climate, (1998). Q. Fu provided high resolution
!                    tables which were appropriately averaged for the
!                    bands in RRTM_LW.  Linear interpolation is used to
!                    get the coefficients from the stored tables.
!                    Valid range of values for dge are between 5.0 and
!                    140.0 micron.
!       LIQFLAG = 0:  The optical depths due to water clouds are computed as
!                     in CCM3.
!       LIQFLAG = 1:  The water droplet effective radius (microns) is input
!                     and the optical depths due to water clouds are computed
!                     as in Hu and Stamnes, J., Clim., 6, 728-742, (1993).
!                     The values for absorption coefficients appropriate for
!                     the spectral bands in RRTM have been obtained for a
!                     range of effective radii by an averaging procedure
!                     based on the work of J. Pinto (private communication).
!                     Linear interpolation is used to get the absorption
!                     coefficients for the input effective radius.

      integer(kind=im), dimension(nbndlw), parameter :: icb = &
        [1, 2, 3, 3, 3, 4, 4, 4, 5, 5, 5, 5, 5, 5, 5, 5]

      hvrclc = '$Revision$'

      ncbands = 1

! This initialization is done in rrtmg_lw_subcol.F90.
!      do lay = 1, nlayers
!         do ig = 1, ngptlw
!            taucmc(ig,lay) = 0.0_rb
!         enddo
!      enddo

! Main layer loop
      do lay = 1, nlayers

        do ig = 1, ngptlw
          cwp = ciwpmc(ig,lay) + clwpmc(ig,lay)
          if (cldfmc(ig,lay) .ge. cldmin .and. &
             (cwp .ge. cldmin .or. taucmc(ig,lay) .ge. cldmin)) then

! Ice clouds and water clouds combined.
            if (inflag .eq. 0) then
! Cloud optical depth already defined in taucmc, return to main program
               return

            elseif(inflag .eq. 1) then
                stop 'INFLAG = 1 OPTION NOT AVAILABLE WITH MCICA'
!               cwp = ciwpmc(ig,lay) + clwpmc(ig,lay)
!               taucmc(ig,lay) = abscld1 * cwp

! Separate treatement of ice clouds and water clouds.
            elseif(inflag .eq. 2) then
               radice = reicmc(lay)

! Calculation of absorption coefficients due to ice clouds.
               if (ciwpmc(ig,lay) .eq. 0.0_rb) then
                  abscoice(ig) = 0.0_rb

               elseif (iceflag .eq. 0) then
#ifdef DEBUG
                  if (radice .lt. 10.0_rb) stop 'ICE RADIUS TOO SMALL'
#endif
                  abscoice(ig) = absice0(1) + absice0(2)/radice

               elseif (iceflag .eq. 1) then
#ifdef DEBUG
                  if (radice .lt. 13.0_rb .or. radice .gt. 130._rb) &
                      stop 'ICE RADIUS OUT OF BOUNDS'
#endif
                  ncbands = 5
                  ib = icb(ngb(ig))
                  abscoice(ig) = absice1(1,ib) + absice1(2,ib)/radice

! For iceflag=2 option, ice particle effective radius is limited to 5.0 to 131.0 microns

               elseif (iceflag .eq. 2) then
#ifdef DEBUG
                  if (radice .lt. 5.0_rb .or. radice .gt. 131.0_rb) &
                    stop 'ICE RADIUS OUT OF BOUNDS'
#endif
                     ncbands = 16
                     factor = (radice - 2._rb)/3._rb
                     indx = int(factor)
                     if (indx .eq. 43) indx = 42
                     fint = factor - real(indx)
                     ib = ngb(ig)
                     abscoice(ig) = &
                         absice2(indx,ib) + fint * &
                         (absice2(indx+1,ib) - (absice2(indx,ib)))

! For iceflag=3 option, ice particle generalized effective size is limited to 5.0 to 140.0 microns

               elseif (iceflag .eq. 3) then
#ifdef DEBUG
                  if (radice .lt. 5.0_rb .or. radice .gt. 140.0_rb) &
                    stop 'ICE GENERALIZED EFFECTIVE SIZE OUT OF BOUNDS'
#endif
                     ncbands = 16
                     factor = (radice - 2._rb)/3._rb
                     indx = int(factor)
                     if (indx .eq. 46) indx = 45
                     fint = factor - real(indx)
                     ib = ngb(ig)
                     abscoice(ig) = &
                         absice3(indx,ib) + fint * &
                         (absice3(indx+1,ib) - (absice3(indx,ib)))

               endif

! Calculation of absorption coefficients due to water clouds.
               if (clwpmc(ig,lay) .eq. 0.0_rb) then
                  abscoliq(ig) = 0.0_rb

               elseif (liqflag .eq. 0) then
                   abscoliq(ig) = absliq0

               elseif (liqflag .eq. 1) then
                  radliq = relqmc(lay)
#ifdef DEBUG
                  if (radliq .lt. 2.5_rb .or. radliq .gt. 60._rb) &
                    stop 'LIQUID EFFECTIVE RADIUS OUT OF BOUNDS'
#endif
                  indx = int(radliq - 1.5_rb)
                  if (indx .eq. 0) indx = 1
                  if (indx .eq. 58) indx = 57
                  fint = radliq - 1.5_rb - real(indx)
                  ib = ngb(ig)
                  abscoliq(ig) = &
                        absliq1(indx,ib) + fint * &
                        (absliq1(indx+1,ib) - (absliq1(indx,ib)))
               endif

               taucmc(ig,lay) = ciwpmc(ig,lay) * abscoice(ig) + &
                                clwpmc(ig,lay) * abscoliq(ig)

            endif
         endif
         enddo
      enddo

      end subroutine cldprmc

      end module rrtmg_lw_cldprmc
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
