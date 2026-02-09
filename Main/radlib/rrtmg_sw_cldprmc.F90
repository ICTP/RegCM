!     path:      $Source$
!     author:    $Author$
!     revision:  $Revision$
!     created:   $Date$

      module rrtmg_sw_cldprmc

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

! ------- Modules -------

      use parkind, only : im => kind_im, rb => kind_rb
      use parrrsw, only : ngptsw, jpband, jpb1, jpb2
      use rrsw_cld, only : extliq1, ssaliq1, asyliq1, &
                           extice2, ssaice2, asyice2, &
                           extice3, ssaice3, asyice3, fdlice3, &
                           abari, bbari, cbari, dbari, ebari, fbari
      use rrsw_wvn, only : wavenum1, wavenum2, ngb
      use rrsw_vsn, only : hvrclc, hnamclc

      implicit none

      contains

! ----------------------------------------------------------------------------
      subroutine cldprmc_sw(nlayers, inflag, iceflag, liqflag, cldfmc, &
                            ciwpmc, clwpmc, reicmc, relqmc, &
                            taormc, taucmc, ssacmc, asmcmc, fsfcmc)
! ----------------------------------------------------------------------------

! Purpose: Compute the cloud optical properties for each cloudy layer
! and g-point interval for use by the McICA method.
! Note: Only inflag = 0 and inflag=2/liqflag=1/iceflag=1,2,3 are available;
! (Hu & Stamnes, Ebert and Curry, Key, and Fu) are implemented.

! ------- Input -------

      integer(kind=im), intent(in) :: nlayers         ! total number of layers
      integer(kind=im), intent(in) :: inflag          ! see definitions
      integer(kind=im), intent(in) :: iceflag         ! see definitions
      integer(kind=im), intent(in) :: liqflag         ! see definitions

      real(kind=rb), intent(in) :: cldfmc(:,:)        ! cloud fraction [mcica]
                                                      !    Dimensions: (ngptsw,nlayers)
      real(kind=rb), intent(in) :: ciwpmc(:,:)        ! cloud ice water path [mcica]
                                                      !    Dimensions: (ngptsw,nlayers)
      real(kind=rb), intent(in) :: clwpmc(:,:)        ! cloud liquid water path [mcica]
                                                      !    Dimensions: (ngptsw,nlayers)
      real(kind=rb), intent(in) :: relqmc(:)          ! cloud liquid particle effective radius (microns)
                                                      !    Dimensions: (nlayers)
      real(kind=rb), intent(in) :: reicmc(:)          ! cloud ice particle effective radius (microns)
                                                      !    Dimensions: (nlayers)
                                                      ! specific definition of reicmc depends on setting of iceflag:
                                                      ! iceflag = 0: (inactive)
                                                      !
                                                      ! iceflag = 1: ice effective radius, r_ec, (Ebert and Curry, 1992),
                                                      !              r_ec range is limited to 13.0 to 130.0 microns
                                                      ! iceflag = 2: ice effective radius, r_k, (Key, Streamer Ref. Manual, 1996)
                                                      !              r_k range is limited to 5.0 to 131.0 microns
                                                      ! iceflag = 3: generalized effective size, dge, (Fu, 1996),
                                                      !              dge range is limited to 5.0 to 140.0 microns
                                                      !              [dge = 1.0315 * r_ec]
      real(kind=rb), intent(in) :: fsfcmc(:,:)        ! cloud forward scattering fraction
                                                      !    Dimensions: (ngptsw,nlayers)

! ------- Output -------

      real(kind=rb), intent(inout) :: taucmc(:,:)     ! cloud optical depth (delta scaled)
                                                      !    Dimensions: (ngptsw,nlayers)
      real(kind=rb), intent(inout) :: ssacmc(:,:)     ! single scattering albedo (delta scaled)
                                                      !    Dimensions: (ngptsw,nlayers)
      real(kind=rb), intent(inout) :: asmcmc(:,:)     ! asymmetry parameter (delta scaled)
                                                      !    Dimensions: (ngptsw,nlayers)
      real(kind=rb), intent(out) :: taormc(:,:)       ! cloud optical depth (non-delta scaled)
                                                      !    Dimensions: (ngptsw,nlayers)

! ------- Local -------

!      integer(kind=im) :: ncbands
      integer(kind=im) :: ib, lay, istr, indx, ig, icx

      real(kind=rb), parameter :: eps = 1.e-06_rb     ! epsilon
      real(kind=rb), parameter :: cldmin = 1.e-20_rb  ! minimum value for cloud quantities
      real(kind=rb) :: cwp                            ! total cloud water path
      real(kind=rb) :: radliq                         ! cloud liquid droplet radius (microns)
      real(kind=rb) :: radice                         ! cloud ice effective size (microns)
      real(kind=rb) :: factor
      real(kind=rb) :: fint

      real(kind=rb) :: taucldorig_a, taucloud_a, ssacloud_a, ffp, ffp1, ffpssa
      real(kind=rb) :: tauiceorig, scatice, ssaice, tauice, tauliqorig, scatliq, ssaliq, tauliq

      real(kind=rb) :: fdelta(ngptsw)
      real(kind=rb) :: extcoice(ngptsw), gice(ngptsw)
      real(kind=rb) :: ssacoice(ngptsw), forwice(ngptsw)
      real(kind=rb) :: extcoliq(ngptsw), gliq(ngptsw)
      real(kind=rb) :: ssacoliq(ngptsw), forwliq(ngptsw)

! Initialize

      hvrclc = '$Revision$'
      icx = 0

! Some of these initializations are done elsewhere
      do lay = 1, nlayers
         do ig = 1, ngptsw
            taormc(ig,lay) = taucmc(ig,lay)
!            taucmc(ig,lay) = 0.0_rb
!            ssacmc(ig,lay) = 1.0_rb
!            asmcmc(ig,lay) = 0.0_rb
         enddo
      enddo

! Main layer loop
      do lay = 1, nlayers

! Main g-point interval loop
         do ig = 1, ngptsw
            cwp = ciwpmc(ig,lay) + clwpmc(ig,lay)
            if (cldfmc(ig,lay) >= cldmin .and. &
               (cwp >= cldmin .or. taucmc(ig,lay) >= cldmin)) then

! (inflag=0): Cloud optical properties input directly
               if (inflag == 0) then
! Cloud optical properties already defined in taucmc, ssacmc, asmcmc are unscaled;
! Apply delta-M scaling here (using Henyey-Greenstein approximation)
                  taucldorig_a = taucmc(ig,lay)
                  ffp = fsfcmc(ig,lay)
                  ffp1 = 1.0_rb - ffp
                  ffpssa = 1.0_rb - ffp * ssacmc(ig,lay)
                  ssacloud_a = ffp1 * ssacmc(ig,lay) / ffpssa
                  taucloud_a = ffpssa * taucldorig_a

                  taormc(ig,lay) = taucldorig_a
                  ssacmc(ig,lay) = ssacloud_a
                  taucmc(ig,lay) = taucloud_a
                  asmcmc(ig,lay) = (asmcmc(ig,lay) - ffp) / (ffp1)

               elseif (inflag == 1) then
                  stop 'INFLAG = 1 OPTION NOT AVAILABLE WITH MCICA'

! (inflag=2): Separate treatement of ice clouds and water clouds.
               elseif (inflag == 2) then
                  radice = reicmc(lay)

! Calculation of absorption coefficients due to ice clouds.
                  if (ciwpmc(ig,lay) == 0.0_rb) then
                     extcoice(ig) = 0.0_rb
                     ssacoice(ig) = 0.0_rb
                     gice(ig)     = 0.0_rb
                     forwice(ig)  = 0.0_rb

! (iceflag = 1):
! Note: This option uses Ebert and Curry approach for all particle sizes similar to
! CAM3 implementation, though this is somewhat unjustified for large ice particles
                  elseif (iceflag == 1) then
#ifdef DEBUG
                     if (radice < 13.0_rb .or. radice > 130._rb) &
                         stop 'ICE RADIUS OUT OF BOUNDS'
#endif
                     ib = ngb(ig)
                     if (wavenum2(ib) > 1.43e04_rb) then
                        icx = 1
                     elseif (wavenum2(ib) > 7.7e03_rb) then
                        icx = 2
                     elseif (wavenum2(ib) > 5.3e03_rb) then
                        icx = 3
                     elseif (wavenum2(ib) > 4.0e03_rb) then
                        icx = 4
                     elseif (wavenum2(ib) >= 2.5e03_rb) then
                        icx = 5
                     endif
                     extcoice(ig) = (abari(icx) + bbari(icx)/radice)
                     ssacoice(ig) = 1._rb - &
                              cbari(icx) - dbari(icx) * radice
                     gice(ig) = ebari(icx) + fbari(icx) * radice
! Check to ensure upper limit of gice is within physical limits for large particles
                     if (gice(ig)>=1._rb) gice(ig) = 1._rb - eps
                     forwice(ig) = gice(ig)*gice(ig)
! Check to ensure all calculated quantities are within physical limits.
#ifdef DEBUG
                     if (extcoice(ig) < 0.0_rb) &
                             stop 'ICE EXTINCTION LESS THAN 0.0'
                     if (ssacoice(ig) > 1.0_rb) &
                             stop 'ICE SSA GRTR THAN 1.0'
                     if (ssacoice(ig) < 0.0_rb) &
                             stop 'ICE SSA LESS THAN 0.0'
                     if (gice(ig) > 1.0_rb) &
                             stop 'ICE ASYM GRTR THAN 1.0'
                     if (gice(ig) < 0.0_rb) &
                             stop 'ICE ASYM LESS THAN 0.0'
#endif

! For iceflag=2 option, ice particle effective radius is limited to 5.0 to 131.0 microns

                  elseif (iceflag == 2) then
#ifdef DEBUG
                     if (radice < 5.0_rb .or. radice > 131.0_rb) &
                          stop 'ICE RADIUS OUT OF BOUNDS'
#endif
                     factor = (radice - 2._rb)/3._rb
                     indx = int(factor)
                     if (indx == 43) indx = 42
                     fint = factor - real(indx,kind=rb)
                     ib = ngb(ig)
                     extcoice(ig) = extice2(indx,ib) + fint * &
                                (extice2(indx+1,ib) -  extice2(indx,ib))
                     ssacoice(ig) = ssaice2(indx,ib) + fint * &
                                (ssaice2(indx+1,ib) -  ssaice2(indx,ib))
                     gice(ig) = asyice2(indx,ib) + fint * &
                                (asyice2(indx+1,ib) -  asyice2(indx,ib))
                     forwice(ig) = gice(ig)*gice(ig)
! Check to ensure all calculated quantities are within physical limits.
#ifdef DEBUG
                     if (extcoice(ig) < 0.0_rb) &
                         stop 'ICE EXTINCTION LESS THAN 0.0'
                     if (ssacoice(ig) > 1.0_rb) &
                         stop 'ICE SSA GRTR THAN 1.0'
                     if (ssacoice(ig) < 0.0_rb) &
                         stop 'ICE SSA LESS THAN 0.0'
                     if (gice(ig) > 1.0_rb) &
                         stop 'ICE ASYM GRTR THAN 1.0'
                     if (gice(ig) < 0.0_rb) &
                         stop 'ICE ASYM LESS THAN 0.0'
#endif

! For iceflag=3 option, ice particle generalized effective size is limited to 5.0 to 140.0 microns

                  elseif (iceflag == 3) then
#ifdef DEBUG
                     if (radice < 5.0_rb .or. radice > 140.0_rb) &
                     stop 'ICE GENERALIZED EFFECTIVE SIZE OUT OF BOUNDS'
#endif
                     factor = (radice - 2._rb)/3._rb
                     indx = int(factor)
                     if (indx == 46) indx = 45
                     fint = factor - real(indx,kind=rb)
                     ib = ngb(ig)
                     extcoice(ig) = extice3(indx,ib) + fint * &
                                 (extice3(indx+1,ib) - extice3(indx,ib))
                     ssacoice(ig) = ssaice3(indx,ib) + fint * &
                                 (ssaice3(indx+1,ib) - ssaice3(indx,ib))
                     gice(ig) = asyice3(indx,ib) + fint * &
                               (asyice3(indx+1,ib) - asyice3(indx,ib))
                     fdelta(ig) = fdlice3(indx,ib) + fint * &
                                 (fdlice3(indx+1,ib) - fdlice3(indx,ib))
#ifdef DEBUG
                     if (fdelta(ig) < 0.0_rb) &
                           stop 'FDELTA LESS THAN 0.0'
                     if (fdelta(ig) > 1.0_rb) &
                           stop 'FDELTA GT THAN 1.0'
#endif
                     forwice(ig) = fdelta(ig) + 0.5_rb / ssacoice(ig)
! See Fu 1996 p. 2067
                     if (forwice(ig) > gice(ig)) &
                          forwice(ig) = gice(ig)
! Check to ensure all calculated quantities are within physical limits.
#ifdef DEBUG
                     if (extcoice(ig) < 0.0_rb) &
                        stop 'ICE EXTINCTION LESS THAN 0.0'
                     if (ssacoice(ig) > 1.0_rb) &
                        stop 'ICE SSA GRTR THAN 1.0'
                     if (ssacoice(ig) < 0.0_rb) &
                        stop 'ICE SSA LESS THAN 0.0'
                     if (gice(ig) > 1.0_rb) &
                        stop 'ICE ASYM GRTR THAN 1.0'
                     if (gice(ig) < 0.0_rb) &
                        stop 'ICE ASYM LESS THAN 0.0'
#endif

                  endif

! Calculation of absorption coefficients due to water clouds.
                  if (clwpmc(ig,lay) == 0.0_rb) then
                     extcoliq(ig) = 0.0_rb
                     ssacoliq(ig) = 0.0_rb
                     gliq(ig) = 0.0_rb
                     forwliq(ig) = 0.0_rb

                  elseif (liqflag == 1) then
                     radliq = relqmc(lay)
#ifdef DEBUG
                     if (radliq < 2.5_rb .or. radliq > 60._rb) stop &
                        'LIQUID EFFECTIVE RADIUS OUT OF BOUNDS'
#endif
                     indx = int(radliq - 1.5_rb)
                     if (indx == 0) indx = 1
                     if (indx == 58) indx = 57
                     fint = radliq - 1.5_rb - real(indx,kind=rb)
                     ib = ngb(ig)
                     extcoliq(ig) = extliq1(indx,ib) + fint * &
                                (extliq1(indx+1,ib) - extliq1(indx,ib))
                     ssacoliq(ig) = ssaliq1(indx,ib) + fint * &
                                (ssaliq1(indx+1,ib) - ssaliq1(indx,ib))
                     if (fint < 0._rb .and. ssacoliq(ig) > 1._rb) &
                                    ssacoliq(ig) = ssaliq1(indx,ib)
                     gliq(ig) = asyliq1(indx,ib) + fint * &
                               (asyliq1(indx+1,ib) - asyliq1(indx,ib))
                     forwliq(ig) = gliq(ig)*gliq(ig)
! Check to ensure all calculated quantities are within physical limits.
#ifdef DEBUG
                     if (extcoliq(ig) < 0.0_rb) &
                        stop 'LIQUID EXTINCTION LESS THAN 0.0'
                     if (ssacoliq(ig) > 1.0_rb) &
                        stop 'LIQUID SSA GRTR THAN 1.0'
                     if (ssacoliq(ig) < 0.0_rb) &
                        stop 'LIQUID SSA LESS THAN 0.0'
                     if (gliq(ig) > 1.0_rb) &
                        stop 'LIQUID ASYM GRTR THAN 1.0'
                     if (gliq(ig) < 0.0_rb) &
                        stop 'LIQUID ASYM LESS THAN 0.0'
#endif
                  endif

                  tauliqorig = clwpmc(ig,lay) * extcoliq(ig)
                  tauiceorig = ciwpmc(ig,lay) * extcoice(ig)
                  taormc(ig,lay) = tauliqorig + tauiceorig

                  ssaliq = ssacoliq(ig) * (1._rb - forwliq(ig)) / &
                          (1._rb - forwliq(ig) * ssacoliq(ig))
                  tauliq = (1._rb - forwliq(ig) * &
                            ssacoliq(ig)) * tauliqorig
                  ssaice = ssacoice(ig) * (1._rb - forwice(ig)) / &
                          (1._rb - forwice(ig) * ssacoice(ig))
                  tauice = (1._rb - forwice(ig) * ssacoice(ig)) * &
                           tauiceorig

                  scatliq = ssaliq * tauliq
                  scatice = ssaice * tauice
                  taucmc(ig,lay) = tauliq + tauice

! Ensure non-zero taucmc and scatice
                  if(taucmc(ig,lay)==0.) taucmc(ig,lay) = cldmin
                  if(scatice==0.) scatice = cldmin

                  ssacmc(ig,lay) = (scatliq + scatice) / taucmc(ig,lay)

                  if (iceflag == 3) then
! In accordance with the 1996 Fu paper, equation A.3,
! the moments for ice were calculated depending on whether using spheres
! or hexagonal ice crystals.
! Set asymetry parameter to first moment (istr=1)
                     istr = 1
                     asmcmc(ig,lay) = (1.0_rb/(scatliq+scatice))* &
                        (scatliq*(gliq(ig)**istr - forwliq(ig)) / &
                        (1.0_rb - forwliq(ig)) + &
                         scatice * ((gice(ig)-forwice(ig))/ &
                        (1.0_rb - forwice(ig)))**istr)

                  else
! This code is the standard method for delta-m scaling.
! Set asymetry parameter to first moment (istr=1)
                     istr = 1
                     asmcmc(ig,lay) = (scatliq *  &
                        (gliq(ig)**istr - forwliq(ig)) / &
                        (1.0_rb - forwliq(ig)) + &
                         scatice * (gice(ig)**istr - forwice(ig)) / &
                        (1.0_rb - forwice(ig)))/(scatliq + scatice)
                  endif

               endif

            endif

! End g-point interval loop
         enddo

! End layer loop
      enddo

      end subroutine cldprmc_sw

      end module rrtmg_sw_cldprmc

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
