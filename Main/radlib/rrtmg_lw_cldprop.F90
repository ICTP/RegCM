!     path:      $Source$
!     author:    $Author$
!     revision:  $Revision$
!     created:   $Date$
!
      module rrtmg_lw_cldprop

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
      use parrrtm, only : nbndlw
      use rrlw_cld, only: abscld1, absliq0, absliq1, &
                          absice0, absice1, absice2, absice3
      use rrlw_vsn, only: hvrcld, hnamcld

      implicit none

      contains

! ------------------------------------------------------------------------------
      subroutine cldprop(nlayers, inflag, iceflag, liqflag, cldfrac, tauc, &
                         ciwp, clwp, rei, rel, ncbands, taucloud)
! ------------------------------------------------------------------------------

! Purpose:  Compute the cloud optical depth(s) for each cloudy layer.

! ------- Input -------

      integer(kind=im), intent(in) :: nlayers         ! total number of layers
      integer(kind=im), intent(in) :: inflag          ! see definitions
      integer(kind=im), intent(in) :: iceflag         ! see definitions
      integer(kind=im), intent(in) :: liqflag         ! see definitions

      real(kind=rb), intent(in) :: cldfrac(:)         ! cloud fraction
                                                      !    Dimensions: (nlayers)
      real(kind=rb), intent(in) :: ciwp(:)            ! cloud ice water path
                                                      !    Dimensions: (nlayers)
      real(kind=rb), intent(in) :: clwp(:)            ! cloud liquid water path
                                                      !    Dimensions: (nlayers)
      real(kind=rb), intent(in) :: rei(:)             ! cloud ice particle effective size (microns)
                                                      !    Dimensions: (nlayers)
                                                      ! specific definition of rei depends on setting of iceflag:
                                                      ! iceflag = 0: ice effective radius, r_ec, (Ebert and Curry, 1992),
                                                      !              r_ec must be >= 10.0 microns
                                                      ! iceflag = 1: ice effective radius, r_ec, (Ebert and Curry, 1992),
                                                      !              r_ec range is limited to 13.0 to 130.0 microns
                                                      ! iceflag = 2: ice effective radius, r_k, (Key, Streamer Ref. Manual, 1996)
                                                      !              r_k range is limited to 5.0 to 131.0 microns
                                                      ! iceflag = 3: generalized effective size, dge, (Fu, 1996),
                                                      !              dge range is limited to 5.0 to 140.0 microns
                                                      !              [dge = 1.0315 * r_ec]
      real(kind=rb), intent(in) :: rel(:)             ! cloud liquid particle effective radius (microns)
                                                      !    Dimensions: (nlayers)
      real(kind=rb), intent(in) :: tauc(:,:)          ! cloud optical depth
                                                      !    Dimensions: (nbndlw,nlayers)

! ------- Output -------

      integer(kind=im), intent(out) :: ncbands        ! number of cloud spectral bands
      real(kind=rb), intent(out) :: taucloud(:,:)     ! cloud optical depth
                                                      !    Dimensions: (nlayers,nbndlw)

! ------- Local -------

      integer(kind=im) :: lay                         ! Layer index
      integer(kind=im) :: ib                          ! spectral band index
      integer(kind=im) :: index
      integer(kind=im) :: iceind
      integer(kind=im) :: liqind
      integer(kind=im) :: icb(nbndlw,0:2)

      real(kind=rb) :: abscoice(nbndlw)               ! ice absorption coefficients
      real(kind=rb) :: abscoliq(nbndlw)               ! liquid absorption coefficients
      real(kind=rb) :: cwp                            ! cloud water path
      real(kind=rb) :: radliq                         ! cloud liquid droplet radius (microns)
      real(kind=rb) :: radice                         ! cloud ice effective size (microns)
      real(kind=rb) :: factor                         !
      real(kind=rb) :: fint                           !
      real(kind=rb) :: tauctot(nlayers)               ! band integrated cloud optical depth
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

      data icb /1,1,1,1,1,1,1,1,1, 1, 1, 1, 1, 1, 1, 1, &
                1,2,3,3,3,4,4,4,5, 5, 5, 5, 5, 5, 5, 5, &
                1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16/

      hvrcld = '$Revision$'

      ncbands = 1
      liqind = 0
      iceind = 0
      tauctot(:) = 0._rb

      do lay = 1, nlayers
         do ib = 1, nbndlw
            taucloud(lay,ib) = 0.0_rb
            tauctot(lay) = tauctot(lay) + tauc(ib,lay)
         enddo
      enddo

! Main layer loop
      do lay = 1, nlayers
         cwp = ciwp(lay) + clwp(lay)
         if (cldfrac(lay) .ge. cldmin .and. &
            (cwp .ge. cldmin .or. tauctot(lay) .ge. cldmin)) then

! Ice clouds and water clouds combined.
            if (inflag .eq. 0) then
               ncbands = 16
               do ib = 1, ncbands
                  taucloud(lay,ib) = tauc(ib,lay)
               end do

            elseif (inflag .eq. 1) then
               ncbands = 16
               do ib = 1, ncbands
                  taucloud(lay,ib) = abscld1 * cwp
               end do

! Separate treatement of ice clouds and water clouds.
            elseif (inflag .eq. 2) then
               radice = rei(lay)

! Calculation of absorption coefficients due to ice clouds.
               if (ciwp(lay) .eq. 0.0_rb) then
                  abscoice(1) = 0.0_rb
                  iceind = 0

               elseif (iceflag .eq. 0) then
#ifdef DEBUG
                  if (radice .lt. 10.0_rb) stop 'ICE RADIUS TOO SMALL'
#endif
                  abscoice(1) = absice0(1) + absice0(2)/radice
                  iceind = 0

               elseif (iceflag .eq. 1) then
#ifdef DEBUG
                  if (radice .lt. 13.0_rb .or. radice .gt. 130._rb) &
                       stop 'ICE RADIUS OUT OF BOUNDS'
#endif
                  ncbands = 5
                  do ib = 1, ncbands
                     abscoice(ib) = absice1(1,ib) + absice1(2,ib)/radice
                  enddo
                  iceind = 1

! For iceflag=2 option, ice particle effective radius is limited to 5.0 to 131.0 microns

               elseif (iceflag .eq. 2) then
#ifdef DEBUG
                  if (radice .lt. 5.0_rb .or. radice .gt. 131.0_rb) &
                     stop 'ICE RADIUS OUT OF BOUNDS'
#endif
                     ncbands = 16
                     factor = (radice - 2._rb)/3._rb
                     index = int(factor)
                     if (index .eq. 43) index = 42
                     fint = factor - real(index)
                     do ib = 1, ncbands
                        abscoice(ib) = &
                            absice2(index,ib) + fint * &
                            (absice2(index+1,ib) - (absice2(index,ib)))
                     enddo
                     iceind = 2

! For iceflag=3 option, ice particle generalized effective size is limited to 5.0 to 140.0 microns

               elseif (iceflag .eq. 3) then
#ifdef DEBUG
                  if (radice .lt. 5.0_rb .or. radice .gt. 140.0_rb) &
                     stop 'ICE GENERALIZED EFFECTIVE SIZE OUT OF BOUNDS'
#endif
                     ncbands = 16
                     factor = (radice - 2._rb)/3._rb
                     index = int(factor)
                     if (index .eq. 46) index = 45
                     fint = factor - real(index)
                     do ib = 1, ncbands
                        abscoice(ib) = &
                          absice3(index,ib) + fint * &
                          (absice3(index+1,ib) - (absice3(index,ib)))
                     enddo
                     iceind = 2

               endif

! Calculation of absorption coefficients due to water clouds.
               if (clwp(lay) .eq. 0.0_rb) then
                  abscoliq(1) = 0.0_rb
                  liqind = 0
                  if (iceind .eq. 1) iceind = 2

               elseif (liqflag .eq. 0) then
                  abscoliq(1) = absliq0
                  liqind = 0
                  if (iceind .eq. 1) iceind = 2

               elseif (liqflag .eq. 1) then
                  radliq = rel(lay)
#ifdef DEBUG
                  if (radliq .lt. 2.5_rb .or. radliq .gt. 60._rb) &
                    stop 'LIQUID EFFECTIVE RADIUS OUT OF BOUNDS'
#endif
                  index = int(radliq - 1.5_rb)
                  if (index .eq. 0) index = 1
                  if (index .eq. 58) index = 57
                  fint = radliq - 1.5_rb - real(index)
                  ncbands = 16
                  do ib = 1, ncbands
                     abscoliq(ib) = &
                         absliq1(index,ib) + fint * &
                         (absliq1(index+1,ib) - (absliq1(index,ib)))
                  enddo
                  liqind = 2
               endif

               do ib = 1, ncbands
                  taucloud(lay,ib) = &
                            ciwp(lay) * abscoice(icb(ib,iceind)) + &
                            clwp(lay) * abscoliq(icb(ib,liqind))
               enddo
            endif
         endif
      enddo

      end subroutine cldprop

      end module rrtmg_lw_cldprop
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
