!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    ICTP RegCM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with ICTP RegCM.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 
      module mod_vecbats

      use mod_runparams
      use mod_message
      use mod_bats
      use mod_lake, only : initlake
      use mod_slice
      use mod_pbldim
      use mod_date
      use mod_bndry
      use mod_drag
      use mod_main
      use mod_service
#ifdef MPP1
      use mod_mppio
#endif

      private

      public :: interf , initb , vecbats , albedov , slice1D

      real(8) , parameter :: dlogtwo = dlog(d_two)
      real(8) , parameter :: dlogten = dlog(d_10)
      real(8) , parameter :: dloglnd = dlog(zlnd)
      real(8) , parameter :: dlogocn = dlog(zoce)
      real(8) , parameter :: dlogsno = dlog(zsno)

      contains

      subroutine vecbats

!=======================================================================
!l  based on: bats version 1e          copyright 18 august 1989
!=======================================================================
!      written by robert e. dickinson and patrick j. kennedy with
! contributions from klaus bluemel, filippo giorgi, and ann henderson-
! sellers.  version 1b (1 dec 1988) greatly modifies versions 1 and 1a,
! for example:  general cleanup and stabilization of iterations;
! three-layer soil moisture; improvement and bug fixes by bluemel
! including an improved canopy model; improvement in soil temperature
! force-restore method for overlying snow; stomatal resistance
! dependence on solar radiation made dependent on lai; inclusion of a
! non-zero displacement level; elimination of soil water movement
! in frozen ground; fixing of the carbon routine;
! and derivation of leaf iteration scheme.
!
!      version 1c (1 march 1989) modifies version 1b by:
! adding a dependence of stomatal resistance on vapor pressure deficit,
! allowing for solar zenith dependence of incident light in
! calculating dependence of stomatal resistance on light;
! further streamlining the iterative calculation of leaf temperature;
! and fixing a bug in the soil moisture calculation of bats1b.
!
!      version 1d (1 june 1989) fixes two bugs introduced into 1c,
! in subroutines lftemp and tgrund, and includes minor changes
! suggested by k. bluemel.  it redefines output fluxes of sensible and
! latent heat to agree exactly with those used for soil energy balance
! to ensure conservation of energy.  some other cleanup is done also.
!
!      verson 1e (18 august 1989) makes changes consistent with 1d to
! insure water conservation.  latent heat of sublimation now applies
! to all soil terms where snow-covered or below freezing.  two minor
! bugs from misplaced lines of code are fixed.  leaf temperature
! calculation reformulated to use average longwave rather than over
! bare soil as done previously, to better match ccm1 tie-in.
!
!      vector bats (1989) was developed from bats1e by gary bates.
! it is designed such that one call to bats perfoms calculations at
! a specified arbitrary number of gridpoints instead of just one point
! as in the standard version of bats1e.  this allows a more efficient
! vectorization of the code.  vector bats was coupled to the mm4
! in februrary 1992 by gary bates.
!
!
!                  main drive program
!
! the core model begins with subroutine bndry, which is driven by
! this routine.
!
! *********************************************************************
!
!
!                      * ts (lowest model layer "midpt" at 75m)
!
!
!                         i-----------i
!                         i         ======= taf (temp air in leaves)
!                         i     tlef  i
!                         i-----------i
!                              i  i
!          ta (anemom)         i  i
!          is actually a ss    i  i
!             i   i temp (4')  i  i
!  -----------i---i------------i--i------------------------------------
!                                                  tg
! -------------------------------------------------------------------
!                                                    tgb
! ********************************************************************
!  **  type1-crop
!  **  type2-short grass                  ****************************
!  **  type3-evergreen needle leaf tree   * soils parameters are a fn of
!  **  type4-deciduous ""  """ "" "" "    * soil colour -- from 1(light)
!  **  type5- "" """"  broadleaf tree     *                  to 8(dark)
!  **  type6- evergreen   """  """"       * soil texture -- from 1(sand)
!  **  type7-tall grass                   *     thru 6(loam) to 12(clay)
!  **  type8-desert                       *
!  **  type9-tundra                       * soil drainage - 7(free),
!  **  type10-irrig crop                  *      8(poor), 9(impeded)
!  **  type11-semi-desert                 *
!  **  type12-ice                         *  (drainage to be used yet)
!  **  type13-bog or marsh                *
!  **  type14-(inland water)              *
!  **  type15-(sea)                       *****************************
!  **  type16-evgr shrub
!  **  type17-decid shrub
!  **  type18-mixed tree
!
!
!
!  flow from this driver and order of subroutines is:
!
!   initb
!   vecbats ==> soilbc
!               albedov
!                bndry   ==>   drag  ==> dragdn  ==> depth
!                             satur
!                            vcover
!                              drip
!                            lftemp   =====================> stomat
!                              drip                          frawat
!                               co2 ===> carbon             frawat
!                            tseice                           root
!                            tgrund                          satur
!                              snow                         lfdrag
!                             water                         condch
!                                                           condcq
!                                                            deriv
!
!=======================================================================
!  si version  - water fluxes are generally calculated in kg/m**2/s.
!  1000 kg/m**2/s = 1 m/s  - converted to energy units for display:
!                            1 kg/m**2/s = 2.5 x 10.e6  watts/m**2.
!  note also  1 kg/m**2/s = 1 mm/m**2/s so fluxes can be thought of
!                            in mm/s.
!=======================================================================
 
      implicit none
!
      character (len=50) :: subroutine_name='vecbats'
      integer :: idindx=0
!
      call time_begin(subroutine_name,idindx)

!---------------------------------------------------------------------

! 
! ****** calculate surface fluxes and hydrology budgets
!
      call soilbc
      call bndry
! 
      call time_end(subroutine_name,idindx)
      end subroutine vecbats
!
      subroutine initb

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!     ***  provides constant initial fields to boundary subroutine
!     ***    soil textures are set in soilbc
!     ***    soil   colors are set in albedo
!
!     ***   units are si
!
      implicit none
! 
      integer :: ill , is , itex , jll , k , nlveg
      real(8) , dimension(20) :: slmo
!
!----------------------------------------------------------------------
!     slmo   : surface moisture availability in fraction of one.
!     data slmo/0.75,0.60,0.75,0.75,0.75,
!     &          0.80,0.60,0.10,0.90,0.90,
!     &          0.10,0.95,0.90,1.00,1.00,
!     &          0.60,0.75,0.85,0.85,0.85 /
      data slmo/0.65D0 , 0.45D0 , 0.60D0 , 0.60D0 , 0.65D0 , &
                0.65D0 , 0.55D0 , 0.10D0 , 0.90D0 , 0.80D0 , &
                0.20D0 , 0.90D0 , 0.90D0 , 1.00D0 , 1.00D0 , &
                0.50D0 , 0.50D0 , 0.65D0 , 0.60D0 , 0.60D0/ !BATS land types
 
!     ****** typically resp=1.0 kg/m**2/s and changes by<10% in 10 days
!
!!
      character (len=50) :: subroutine_name='initb'
      integer :: idindx=0
!
      call time_begin(subroutine_name,idindx)
!
#ifdef MPP1
      do jll = 1 , jendx
#else
#ifdef BAND
      do jll = 1 , jx
#else
      do jll = 1 , jxm1
#endif
#endif
        do ill = 1 , iym1
          pptnc(ill,jll) = d_zero
          pptc(ill,jll) = d_zero
          veg2d(ill,jll) = idnint(mddom%satbrt(ill,jll))
        end do
        do ill = 1 , iym1
          do k = 1 , nnsg
            veg2d1(k,ill,jll) = idnint(satbrt1(k,ill,jll))
          end do
        end do
      end do

!     ******  initialize hostetler lake model
      if (lakemod == 1) call initlake
 
#ifdef MPP1
      do jll = 1 , jendx
#else
#ifdef BAND
      do jll = 1 , jx
#else
      do jll = 1 , jxm1
#endif
#endif
        do ill = 1 , iym1
          do k = 1 , nnsg
            tg2d(k,ill,jll) = sts2%tg(ill,jll)
            tgb2d(k,ill,jll) = sts2%tg(ill,jll)
            taf2d(k,ill,jll) = sts2%tg(ill,jll)
            tlef2d(k,ill,jll) = sts2%tg(ill,jll)

            if ( ocld2d(k,ill,jll) == 2 ) then
              if ( iseaice == 1 .or. lakemod == 1 ) then
                sice2d(k,ill,jll) = d_1000
                scv2d(k,ill,jll) = d_zero
              end if
              nlveg = 12
            else if ( ocld2d(k,ill,jll) == 1 ) then
              sice2d(k,ill,jll) = d_zero
              scv2d(k,ill,jll) = dmax1(scv2d(k,ill,jll),d_zero)
              nlveg = veg2d1(k,ill,jll)
            else
              sice2d(k,ill,jll) = d_zero
              scv2d(k,ill,jll) = d_zero
              nlveg = veg2d1(k,ill,jll)
            end if
!           ******  initialize soil moisture in the 3 layers
            is = idint(satbrt1(k,ill,jll))
            itex = iexsol(nlveg)
            swt2d(k,ill,jll) = deptv(nlveg)*xmopor(itex)*slmo(is)
            srw2d(k,ill,jll) = deprv(nlveg)*xmopor(itex)*slmo(is)
            ssw2d(k,ill,jll) = depuv(nlveg)*xmopor(itex)*slmo(is)

            dew2d(k,ill,jll) = d_zero
            sag2d(k,ill,jll) = d_zero
            gwet2d(k,ill,jll) = d_half
            sena2d(k,ill,jll) = d_zero
            evpa2d(k,ill,jll) = d_zero
            rnos2d(k,ill,jll) = d_zero
            rno2d(k,ill,jll) = d_zero
            ircp2d(k,ill,jll) = d_zero
          end do
        end do
      end do
 
#ifdef MPP1
      do jll = 1 , jendx
#else
#ifdef BAND
      do jll = 1 , jx
#else
      do jll = 1 , jxm1
#endif
#endif
        do ill = 1 , iym1
          fsw2d(ill,jll) = d_zero
          flw2d(ill,jll) = d_zero
          sabv2d(ill,jll) = d_zero
          sol2d(ill,jll) = d_zero
          fswa2d(ill,jll) = d_zero
          flwa2d(ill,jll) = d_zero
          prca2d(ill,jll) = d_zero
          prnca2d(ill,jll) = d_zero
          svga2d(ill,jll) = d_zero
          sina2d(ill,jll) = d_zero
        end do
      end do
 
      call time_end(subroutine_name,idindx)
      end subroutine initb
!
      subroutine interf(ivers,j,k,istart,iend,ng)

! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!  this subroutine interfaces mm42d and bats variables
!
!  ivers = 1 ,   regcm2d --> bats
!  ivers = 2 ,   bats --> regcm2d
!
      implicit none
!
      integer , intent (in) :: ivers , j , k , istart , iend , ng
!
      real(8) :: amxtem , facb , facs , fact , factuv , facv , fracb ,  &
               & fracs , fracv , hl , mmpd , rh0 , satvp , sfac ,       &
               & solvt , wpm2
      integer :: i , n , nnn
      real(4) :: real_4
!
!     ******    fbat contains nummap+1 fields to be written out to
!     iutbat ******    fields are written out in following order:
!     1.  anemon west  wind (41)
!     2.  anemom south wind (42)
!     3.  drag- surface stress, in si (48)
!     4.  ground temp (4)
!     5.  temp of foliage (7)
!     6.  anemom temp (5)
!     7.  anemom spec. humidity (15)
!     8.  upper layer soil water (26)
!     9.  root zone soil water (25)
!     10.  accum precip (19)
!     11.  accum evap (39)
!     12.  accum surf runoff (12)
!     13.  snow depth in mm h2o (30)
!     14.  accum sensible heat (40)
!     15.  accum net ir (38)
!     16.  accum net solar abs (37)
!     17.  accum downward ir
!     18.  accum solar incident at surface
!     19.  convective precipitation
!     20.  surface pressure
!     21.  pbl height
!!
      character (len=50) :: subroutine_name='interf'
      integer :: idindx=0
!
      call time_begin(subroutine_name,idindx)
 
      if ( ivers == 1 ) then ! regcm2d --> bats

        do i = istart, iend
          do n = 1 , ng
            p1d0(n,i) = (sps2%ps(i,j)+r8pt)*d_1000
            ts1d0(n,i) = thx3d(i,k,j)
            qs1d0(n,i) = qvb3d(i,k,j)/(d_one+qvb3d(i,k,j))
            qs1d(n,i) = qs1d0(n,i)
 
            hl = lh0 - lh1*(ts1d0(n,i)-tzero)
            satvp = lsvp1*dexp(lsvp2*hl*(d_one/tzero-d_one/ts1d0(n,i)))
            rh0 = dmax1(qs1d0(n,i)/(ep2*satvp/(p1d0(n,i)* &
                        0.01D0-satvp)),d_zero)
 
            ts1d(n,i) = ts1d0(n,i) - lrate*regrav*(ht1(n,i,j)- &
                                                 mddom%ht(i,j))
            p1d(n,i) = p1d0(n,i)*(ts1d(n,i)/ts1d0(n,i))
 
            hl = lh0 - lh1*(ts1d(n,i)-tzero)
            satvp = lsvp1*dexp(lsvp2*hl*(d_one/tzero-d_one/ts1d(n,i)))
            qs1d(n,i) = dmax1(rh0*ep2*satvp/(p1d(n,i)*&
                              0.01D0-satvp),d_zero)
 
            tg1d(n,i) = tg2d(n,i,j)
            rhs1d(n,i) = p1d(n,i)/(rgas*ts1d(n,i))
            prcp1d(n,i) = pptnc(i,j) + pptc(i,j)
!
!           quantities stored on 2d surface array for bats use only
!
            tgb1d(n,i) = tgb2d(n,i,j)
            taf1d(n,i) = taf2d(n,i,j)
            tlef1d(n,i) = tlef2d(n,i,j)
            tsw1d(n,i) = swt2d(n,i,j)
            rsw1d(n,i) = srw2d(n,i,j)
            ssw1d(n,i) = ssw2d(n,i,j)
            ldew1d(n,i) = dew2d(n,i,j)
            sag1d(n,i) = sag2d(n,i,j)
            scv1d(n,i) = scv2d(n,i,j)
            sice1d(n,i) = sice2d(n,i,j)
            gwet1d(n,i) = gwet2d(n,i,j)
            sent1d(n,i) = sfsta%hfx(i,j)
            evpr1d(n,i) = sfsta%qfx(i,j)
            ldoc1d(n,i) = ocld2d(n,i,j)
            ircp1d(n,i) = ircp2d(n,i,j)
            lveg(n,i) = veg2d1(n,i,j)
            oveg(n,i) = lveg(n,i)
            if ( ldoc1d(n,i) == 2 ) lveg(n,i) = 12
            amxtem = dmax1(298.0D0-tgb1d(n,i),d_zero)
            sfac = d_one - dmax1(d_zero,d_one-0.0016D0*amxtem**d_two)
            veg1d(n,i) = vegc(lveg(n,i)) - seasf(lveg(n,i))*sfac
            emiss_1d(n,i) = emiss2d(n,i,j)
            z1d(n,i) = za(i,k,j)
            z1log(n,i)  = dlog(z1d(n,i))
            z2fra(n,i)  = z1log(n,i) - dlogtwo
            z10fra(n,i) = z1log(n,i) - dlogten
            zlgocn(n,i) = z1log(n,i) - dlogocn
            zlglnd(n,i) = z1log(n,i) - dloglnd
            zlgsno(n,i) = z1log(n,i) - dlogsno
            zlgveg(n,i) = z1log(n,i) - dlog(rough(lveg(n,i)))
            zlgdis(n,i) = dlog(z1d(n,i)-displa(lveg(n,i)) / &
                               rough(lveg(n,i)))
          end do

 
          rh0 = d_zero
          do n = 1 , ng
            rh0 = rh0 + (qs1d(n,i)-qs1d0(n,i))
          end do
          rh0 = rh0/ng
          do n = 1 , ng
            qs1d(n,i) = dmax1(qs1d(n,i)-rh0,d_zero)
          end do
 
          us1d(i) = ubx3d(i,k,j)
          vs1d(i) = vbx3d(i,k,j)
          fsw1d(i) = fsw2d(i,j)
          flw1d(i) = flw2d(i,j)
          solis(i) = sol2d(i,j)
          sabveg(i) = sabv2d(i,j)
          solvt = solvd2d(i,j) + solvs2d(i,j)
          if ( solvt > d_zero ) then
            fracd(i) = solvd2d(i,j)/solvt
          else
            fracd(i) = 0.2D0
          end if
          czen(i) = dmax1(coszrs(i),d_zero)
        end do
 
      else if ( ivers == 2 ) then ! bats --> regcm2d
 
        do i = istart, iend
          sfsta%uvdrag(i,j) = d_zero
          sfsta%hfx(i,j) = d_zero
          sfsta%qfx(i,j) = d_zero
          sts2%tg(i,j) = d_zero
          sts1%tg(i,j) = d_zero
          sfsta%tgbb(i,j) = d_zero
          if ( ichem == 1 ) then
            ssw2da(i,j) = d_zero
            sdeltk2d(i,j) = d_zero
            sdelqk2d(i,j) = d_zero
            sfracv2d(i,j) = d_zero
            sfracb2d(i,j) = d_zero
            sfracs2d(i,j) = d_zero
            svegfrac2d(i,j) = d_zero
          end if

          do n = 1 , ng
            sfsta%uvdrag(i,j) = sfsta%uvdrag(i,j) + drag1d(n,i)
            sfsta%hfx(i,j) = sfsta%hfx(i,j) + sent1d(n,i)
            sfsta%qfx(i,j) = sfsta%qfx(i,j) + evpr1d(n,i)
            sts2%tg(i,j) = sts2%tg(i,j) + tg1d(n,i)
            sts1%tg(i,j) = sts1%tg(i,j) + tg1d(n,i)
            if ( ichem == 1 ) then
              ssw2da(i,j) = ssw2da(i,j) + ssw1d(n,i)
              sdeltk2d(i,j) = sdeltk2d(i,j) + delt1d(n,i)
              sdelqk2d(i,j) = sdelqk2d(i,j) + delq1d(n,i)
              sfracv2d(i,j) = sfracv2d(i,j) + sigf(n,i)
              sfracb2d(i,j) = sfracb2d(i,j) + (d_one-veg1d(n,i))    &
                            & *(d_one-scvk(n,i))
              sfracs2d(i,j) = sfracs2d(i,j) + veg1d(n,i)*wt(n,i) &
                            & + (d_one-veg1d(n,i))*scvk(n,i)
              svegfrac2d(i,j) = svegfrac2d(i,j) + veg1d(n,i)
            end if
            if ( iocnflx == 1 .or. &
                (iocnflx == 2 .and. ldoc1d(n,i) /= 0 ) ) then
              sfsta%tgbb(i,j) = sfsta%tgbb(i,j)                   &
                         + ((d_one-veg1d(n,i))*tg1d(n,i)**d_four+ &
                         veg1d(n,i)*tlef1d(n,i)**d_four)**d_rfour
            else
              sfsta%tgbb(i,j) = sfsta%tgbb(i,j) + tg1d(n,i)
            end if
            if ( ldoc1d(n,i) == 0 ) then
              ssw1d(n,i)  = dmissval
              rsw1d(n,i)  = dmissval
              tsw1d(n,i)  = dmissval
              rno1d(n,i)  = dmissval
              rnos1d(n,i) = dmissval
              scv1d(n,i)  = dmissval
            end if
          end do
          sfsta%uvdrag(i,j) = sfsta%uvdrag(i,j)/dble(ng)
          sfsta%hfx(i,j) = sfsta%hfx(i,j)/dble(ng)
          sfsta%qfx(i,j) = sfsta%qfx(i,j)/dble(ng)
          sts2%tg(i,j) = sts2%tg(i,j)/dble(ng)
          sts1%tg(i,j) = sts1%tg(i,j)/dble(ng)
          sfsta%tgbb(i,j) = sfsta%tgbb(i,j)/dble(ng)

          if ( ichem == 1 ) then
            ssw2da(i,j) = ssw2da(i,j)/dble(ng)
            sdeltk2d(i,j) = sdeltk2d(i,j)/dble(ng)
            sdelqk2d(i,j) = sdelqk2d(i,j)/dble(ng)
            sfracv2d(i,j) = sfracv2d(i,j)/dble(ng)
            sfracb2d(i,j) = sfracb2d(i,j)/dble(ng)
            sfracs2d(i,j) = sfracs2d(i,j)/dble(ng)
            svegfrac2d(i,j) = svegfrac2d(i,j)/dble(ng)
          end if
          do n = 1 , ng
            scv2d(n,i,j) = scv1d(n,i)
            tg2d(n,i,j) = tg1d(n,i)
            tgb2d(n,i,j) = tgb1d(n,i)
            taf2d(n,i,j) = taf1d(n,i)
            tlef2d(n,i,j) = tlef1d(n,i)
            swt2d(n,i,j) = tsw1d(n,i)
            srw2d(n,i,j) = rsw1d(n,i)
            ssw2d(n,i,j) = ssw1d(n,i)
            dew2d(n,i,j) = ldew1d(n,i)
            sag2d(n,i,j) = sag1d(n,i)
            sice2d(n,i,j) = sice1d(n,i)
            gwet2d(n,i,j) = gwet1d(n,i)
            ocld2d(n,i,j) = ldoc1d(n,i)
            ircp2d(n,i,j) = ircp1d(n,i)
            evpa2d(n,i,j) = evpa2d(n,i,j) + dtbat*evpr1d(n,i)
            sena2d(n,i,j) = sena2d(n,i,j) + dtbat*sent1d(n,i)
            if ( dabs(rnos1d(n,i)) > 1.0D-10 ) then
              rnos2d(n,i,j) = rnos2d(n,i,j) + rnos1d(n,i)*dtbat
            end if
            if ( dabs(rnos1d(n,i))  > 1.0D-10 .and. &
                 dabs(rno1d(n,i))   > 1.0D-10 ) then
              rno2d(n,i,j) = rno2d(n,i,j) + &
                     (rno1d(n,i)-rnos1d(n,i))*dtbat
            end if
          end do
!
!         quantities stored on 2d surface array for bats use only
!
          prca2d(i,j) = prca2d(i,j) + dtbat*pptc(i,j)
          prnca2d(i,j) = prnca2d(i,j) + dtbat*pptnc(i,j)
          if ( prnca2d(i,j) < dlowval ) prnca2d(i,j) = d_zero
          if ( prca2d(i,j) < dlowval ) prca2d(i,j) = d_zero
          flwa2d(i,j) = flwa2d(i,j) + dtbat*flw1d(i)
          flwda2d(i,j) = flwda2d(i,j) + dtbat*flwd2d(i,j)
          fswa2d(i,j) = fswa2d(i,j) + dtbat*fsw1d(i)
          svga2d(i,j) = svga2d(i,j) + dtbat*sabveg(i)
          sina2d(i,j) = sina2d(i,j) + dtbat*sinc2d(i,j)
          pptnc(i,j) = d_zero
          pptc(i,j) = d_zero
        end do

        do i = istart, iend
#ifdef MPP1
          u10m_o(j,i-1) = 0.0
          v10m_o(j,i-1) = 0.0
          tg_o(j,i-1) = 0.0
          t2m_o(j,i-1) = 0.0
          aldirs_o(j,i-1) = 0.0
          aldifs_o(j,i-1) = 0.0
          do n = 1 , ng
            if ( ldoc1d(n,i) /= 0 ) then
              fracv = sigf(n,i)
              fracb = (d_one-veg1d(n,i))*(d_one-scvk(n,i))
              fracs = veg1d(n,i)*wt(n,i) + (d_one-veg1d(n,i))*scvk(n,i)
              facv = z2fra(n,i)/zlgveg(n,i)
              facb = z2fra(n,i)/zlglnd(n,i)
              facs = z2fra(n,i)/zlgsno(n,i)
              fact = fracv*facv + fracb*facb + fracs*facs
              facv = z10fra(n,i)/zlgveg(n,i)
              facb = z10fra(n,i)/zlglnd(n,i)
              facs = z10fra(n,i)/zlgsno(n,i)
              factuv = fracv*facv + fracb*facb + fracs*facs
              u10m1d(n,i) = us1d(i)*(d_one-factuv)
              v10m1d(n,i) = vs1d(i)*(d_one-factuv)
              t2m_1d(n,i) = ts1d(n,i) - delt1d(n,i)*fact
            else 
              if ( iocnflx == 1 ) then
                fact = z2fra(n,i)/zlgocn(n,i)
                factuv = z10fra(n,i)/zlgocn(n,i)
                u10m1d(n,i) = us1d(i)*(d_one-factuv)
                v10m1d(n,i) = vs1d(i)*(d_one-factuv)
                t2m_1d(n,i) = ts1d(n,i) - delt1d(n,i)*fact
              end if
            end if
            tg_s(n,j,i-1) = real(tg1d(n,i))
            u10m_s(n,j,i-1) = real(u10m1d(n,i))
            v10m_s(n,j,i-1) = real(v10m1d(n,i))
            t2m_s(n,j,i-1) = real(t2m_1d(n,i))
 
            u10m_o(j,i-1) = u10m_o(j,i-1) + real(u10m1d(n,i))
            v10m_o(j,i-1) = v10m_o(j,i-1) + real(v10m1d(n,i))
            t2m_o(j,i-1) = t2m_o(j,i-1) + real(t2m_1d(n,i))
            tg_o(j,i-1) = tg_o(j,i-1) + real(tg1d(n,i))
            aldirs_o(j,i-1) = aldirs_o(j,i-1) + real(aldirs1d(n,i))
            aldifs_o(j,i-1) = aldifs_o(j,i-1) + real(aldifs1d(n,i))
          end do
          u10m_o(j,i-1) = u10m_o(j,i-1)/real(ng)
          v10m_o(j,i-1) = v10m_o(j,i-1)/real(ng)
          t2m_o(j,i-1) = t2m_o(j,i-1)/real(ng)
          tg_o(j,i-1) = tg_o(j,i-1)/real(ng)
          aldirs_o(j,i-1) = aldirs_o(j,i-1)/real(ng)
          aldifs_o(j,i-1) = aldifs_o(j,i-1)/real(ng)
 
          tgmx_o(j,i-1) = amax1(tgmx_o(j,i-1),tg_o(j,i-1))
          tgmn_o(j,i-1) = amin1(tgmn_o(j,i-1),tg_o(j,i-1))
          t2mx_o(j,i-1) = amax1(t2mx_o(j,i-1),t2m_o(j,i-1))
          t2mn_o(j,i-1) = amin1(t2mn_o(j,i-1),t2m_o(j,i-1))
          w10x_o(j,i-1) = amax1(w10x_o(j,i-1), &
                          sqrt(u10m_o(j,i-1)**2.0+v10m_o(j,i-1)**2.0))
          real_4 = real((sps2%ps(i,j)+r8pt)*d_10)
          psmn_o(j,i-1) = amin1(psmn_o(j,i-1),real_4)

#else
#ifdef BAND
          u10m_o(j,i-1) = 0.0
          v10m_o(j,i-1) = 0.0
          tg_o(j,i-1) = 0.0
          t2m_o(j,i-1) = 0.0
          aldirs_o(j,i-1) = 0.0
          aldifs_o(j,i-1) = 0.0
          do n = 1 , ng
            if ( ldoc1d(n,i) /= 0 ) then
              fracv = sigf(n,i)
              fracb = (d_one-veg1d(n,i))*(d_one-scvk(n,i))
              fracs = veg1d(n,i)*wt(n,i) + (d_one-veg1d(n,i))*scvk(n,i)
              facv = z2fra(n,i)/zlgveg(n,i)
              facb = z2fra(n,i)/zlglnd(n,i)
              facs = z2fra(n,i)/zlgsno(n,i)
              fact = fracv*facv + fracb*facb + fracs*facs
              facv = z10fra(n,i)/zlgveg(n,i)
              facb = z10fra(n,i)/zlglnd(n,i)
              facs = z10fra(n,i)/zlgsno(n,i)
              factuv = fracv*facv + fracb*facb + fracs*facs
              u10m1d(n,i) = us1d(i)*(d_one-factuv)
              v10m1d(n,i) = vs1d(i)*(d_one-factuv)
              t2m_1d(n,i) = ts1d(n,i) - delt1d(n,i)*fact
            else 
              if ( iocnflx == 1 ) then
                fact = z2fra(n,i)/zlgocn(n,i)
                factuv = z10fra(n,i)/zlgocn(n,i)
                u10m1d(n,i) = us1d(i)*(d_one-factuv)
                v10m1d(n,i) = vs1d(i)*(d_one-factuv)
                t2m_1d(n,i) = ts1d(n,i) - delt1d(n,i)*fact
              end if
            end if
            tg_s(n,j,i-1) = real(tg1d(n,i))
            u10m_s(n,j,i-1) = real(u10m1d(n,i))
            v10m_s(n,j,i-1) = real(v10m1d(n,i))
            t2m_s(n,j,i-1) = real(t2m_1d(n,i))

            u10m_o(j,i-1) = u10m_o(j,i-1) + real(u10m1d(n,i))
            v10m_o(j,i-1) = v10m_o(j,i-1) + real(v10m1d(n,i))
            t2m_o(j,i-1) = t2m_o(j,i-1) + real(t2m_1d(n,i))
            tg_o(j,i-1) = tg_o(j,i-1) + real(tg1d(n,i))
            aldirs_o(j,i-1) = aldirs_o(j,i-1) + real(aldirs1d(n,i))
            aldifs_o(j,i-1) = aldifs_o(j,i-1) + real(aldifs1d(n,i))
          end do
          u10m_o(j,i-1) = u10m_o(j,i-1)/real(ng)
          v10m_o(j,i-1) = v10m_o(j,i-1)/real(ng)
          t2m_o(j,i-1) = t2m_o(j,i-1)/real(ng)
          tg_o(j,i-1) = tg_o(j,i-1)/real(ng)
          aldirs_o(j,i-1) = aldirs_o(j,i-1)/real(ng)
          aldifs_o(j,i-1) = aldifs_o(j,i-1)/real(ng)
          tgmx_o(j,i-1) = amax1(tgmx_o(j,i-1),tg_o(j,i-1))
          tgmn_o(j,i-1) = amin1(tgmn_o(j,i-1),tg_o(j,i-1))
          t2mx_o(j,i-1) = amax1(t2mx_o(j,i-1),t2m_o(j,i-1))
          t2mn_o(j,i-1) = amin1(t2mn_o(j,i-1),t2m_o(j,i-1))
          w10x_o(j,i-1) = amax1(w10x_o(j,i-1), &
                           sqrt(u10m_o(j,i-1)**2.0+v10m_o(j,i-1)**2.0))
          real_4 = real((sps2%ps(i,j)+r8pt)*d_10)
          psmn_o(j,i-1) = amin1(psmn_o(j,i-1),real_4)
#else
          u10m_o(j-1,i-1) = 0.0
          v10m_o(j-1,i-1) = 0.0
          tg_o(j-1,i-1) = 0.0
          t2m_o(j-1,i-1) = 0.0
          aldirs_o(j-1,i-1) = 0.0
          aldifs_o(j-1,i-1) = 0.0
          do n = 1 , ng
            if ( ldoc1d(n,i) /= 0 ) then
              fracv = sigf(n,i)
              fracb = (d_one-veg1d(n,i))*(d_one-scvk(n,i))
              fracs = veg1d(n,i)*wt(n,i) + (d_one-veg1d(n,i))*scvk(n,i)
              facv = z2fra(n,i)/zlgveg(n,i)
              facb = z2fra(n,i)/zlglnd(n,i)
              facs = z2fra(n,i)/zlgsno(n,i)
              fact = fracv*facv + fracb*facb + fracs*facs
              facv = z10fra(n,i)/zlgveg(n,i)
              facb = z10fra(n,i)/zlglnd(n,i)
              facs = z10fra(n,i)/zlgsno(n,i)
              factuv = fracv*facv + fracb*facb + fracs*facs
              u10m1d(n,i) = us1d(i)*(d_one-factuv)
              v10m1d(n,i) = vs1d(i)*(d_one-factuv)
              t2m_1d(n,i) = ts1d(n,i) - delt1d(n,i)*fact
            else 
              if ( iocnflx == 1 ) then
                fact = z2fra(n,i)/zlgocn(n,i)
                factuv = z10fra(n,i)/zlgocn(n,i)
                u10m1d(n,i) = us1d(i)*(d_one-factuv)
                v10m1d(n,i) = vs1d(i)*(d_one-factuv)
                t2m_1d(n,i) = ts1d(n,i) - delt1d(n,i)*fact
              end if
            end if
            tg_s(n,j-1,i-1) = real(tg1d(n,i))
            u10m_s(n,j-1,i-1) = real(u10m1d(n,i))
            v10m_s(n,j-1,i-1) = real(v10m1d(n,i))
            t2m_s(n,j-1,i-1) = real(t2m_1d(n,i))

            u10m_o(j-1,i-1) = u10m_o(j-1,i-1) + real(u10m1d(n,i))
            v10m_o(j-1,i-1) = v10m_o(j-1,i-1) + real(v10m1d(n,i))
            t2m_o(j-1,i-1) = t2m_o(j-1,i-1) + real(t2m_1d(n,i))
            tg_o(j-1,i-1) = tg_o(j-1,i-1) + real(tg1d(n,i))
            aldirs_o(j-1,i-1) = aldirs_o(j-1,i-1) + real(aldirs1d(n,i))
            aldifs_o(j-1,i-1) = aldifs_o(j-1,i-1) + real(aldifs1d(n,i))
          end do
          u10m_o(j-1,i-1) = u10m_o(j-1,i-1)/real(ng)
          v10m_o(j-1,i-1) = v10m_o(j-1,i-1)/real(ng)
          t2m_o(j-1,i-1) = t2m_o(j-1,i-1)/real(ng)
          tg_o(j-1,i-1) = tg_o(j-1,i-1)/real(ng)
          aldirs_o(j-1,i-1) = aldirs_o(j-1,i-1)/real(ng)
          aldifs_o(j-1,i-1) = aldifs_o(j-1,i-1)/real(ng)
          tgmx_o(j-1,i-1) = amax1(tgmx_o(j-1,i-1),tg_o(j-1,i-1))
          tgmn_o(j-1,i-1) = amin1(tgmn_o(j-1,i-1),tg_o(j-1,i-1))
          t2mx_o(j-1,i-1) = amax1(t2mx_o(j-1,i-1),t2m_o(j-1,i-1))
          t2mn_o(j-1,i-1) = amin1(t2mn_o(j-1,i-1),t2m_o(j-1,i-1))
          w10x_o(j-1,i-1) = amax1(w10x_o(j-1,i-1), &
                      sqrt(u10m_o(j-1,i-1)**2.0+v10m_o(j-1,i-1)**2.0))
          real_4 = real((sps2%ps(i,j)+r8pt)*d_10)
          psmn_o(j-1,i-1) = amin1(psmn_o(j-1,i-1),real_4)
#endif
#endif
        end do

        if ( mod(ntime+idnint(dtmin*minph),kbats) == 0 .or.    &
            ( jyear == jyear0 .and. ktau == 0 ) .or.       & 
            ( ifrest .and. .not. done_restart ) ) then
          if ( jyear == jyear0 .and. ktau <= 1 ) then
            mmpd = secpd/dtbat
            wpm2 = d_one/dtbat
          else if ( jyear == jyear0 .and. dble(ktau*dtmin)              &
                  &  <= batfrq*minph+0.01D0 ) then
            mmpd = houpd/(batfrq-dtmin/minph)
            wpm2 = d_one/((batfrq-dtmin/minph)*secph)
          else
            mmpd = houpd/batfrq
            wpm2 = d_one/(batfrq*secph)
          end if
          do i = istart, iend
#ifdef MPP1
            drag_o(j,i-1) = 0.0
            q2m_o(j,i-1) = 0.0
            evpa_o(j,i-1) = 0.0
            sena_o(j,i-1) = 0.0
            do n = 1 , ng
              if ( ldoc1d(n,i) /= 0 ) then
                fracv = sigf(n,i)
                fracb = (d_one-veg1d(n,i))*(d_one-scvk(n,i))
                fracs = veg1d(n,i)*wt(n,i) + (d_one-veg1d(n,i))*scvk(n,i)
                facv = z2fra(n,i)/zlgveg(n,i)
                facb = z2fra(n,i)/zlglnd(n,i)
                facs = z2fra(n,i)/zlgsno(n,i)
                fact = fracv*facv + fracb*facb + fracs*facs
                q2m_1d(n,i) = qs1d(n,i) - delq1d(n,i)*fact
              else
                if ( iocnflx == 1 ) then
                  fact = z2fra(n,i)/zlgocn(n,i)
                  q2m_1d(n,i) = qs1d(n,i) - delq1d(n,i)*fact
                end if
              end if
              q2m_s(n,j,i-1) = real(q2m_1d(n,i))
              drag_s(n,j,i-1) = real(drag1d(n,i))
              evpa_s(n,j,i-1) = real(evpa2d(n,i,j)*mmpd)
              sena_s(n,j,i-1) = real(sena2d(n,i,j)*wpm2)
              tpr_s(n,j,i-1) = real((prnca2d(i,j)+prca2d(i,j))*mmpd)
              prcv_s(n,j,i-1) = real(prca2d(i,j)*mmpd)
              ps_s(n,j,i-1) = real(p1d(n,i)*0.01D0)
 
              q2m_o(j,i-1) = real(q2m_o(j,i-1) + q2m_1d(n,i))
              drag_o(j,i-1) = real(drag_o(j,i-1) + drag1d(n,i))
              evpa_o(j,i-1) = real(evpa_o(j,i-1) + evpa2d(n,i,j))
              sena_o(j,i-1) = real(sena_o(j,i-1) + sena2d(n,i,j))
            end do
            tpr_o(j,i-1) = real((prnca2d(i,j)+prca2d(i,j))*mmpd)
            q2m_o(j,i-1) = q2m_o(j,i-1)/real(ng)
            drag_o(j,i-1) = drag_o(j,i-1)/real(ng)
            evpa_o(j,i-1) = evpa_o(j,i-1)/real(ng)*real(mmpd)
            sena_o(j,i-1) = sena_o(j,i-1)/real(ng)*real(wpm2)
            flwa_o(j,i-1) = real(flwa2d(i,j)*wpm2)
            fswa_o(j,i-1) = real(fswa2d(i,j)*wpm2)
            flwd_o(j,i-1) = real(flwda2d(i,j)*wpm2)
            sina_o(j,i-1) = real(sina2d(i,j)*wpm2)
            prcv_o(j,i-1) = real(prca2d(i,j)*mmpd)
            ps_o(j,i-1) = real((sps2%ps(i,j)+r8pt)*d_10)
            zpbl_o(j,i-1) = real(sfsta%zpbl(i,j))
 
            tlef_o(j,i-1) = 0.0
            ssw_o(j,i-1) = 0.0
            rsw_o(j,i-1) = 0.0
            rnos_o(j,i-1) = 0.0
            scv_o(j,i-1) = 0.0
            nnn = 0
            do n = 1 , ng
              if ( ldoc1d(n,i) /= 0 ) then
                tlef_o(j,i-1) = tlef_o(j,i-1) + real(tlef1d(n,i))
                ssw_o(j,i-1) = ssw_o(j,i-1) + real(ssw1d(n,i))
                rsw_o(j,i-1) = rsw_o(j,i-1) + real(rsw1d(n,i))
                rnos_o(j,i-1) = rnos_o(j,i-1) + real(rnos2d(n,i,j))
                if (dabs(scv1d(n,i)) > dlowval) then
                  scv_o(j,i-1) = scv_o(j,i-1) + real(scv1d(n,i))
                end if
                tlef_s(n,j,i-1) = real(tlef1d(n,i))
                ssw_s(n,j,i-1) = real(ssw1d(n,i))
                rsw_s(n,j,i-1) = real(rsw1d(n,i))
                rnos_s(n,j,i-1) = real(rnos2d(n,i,j)*mmpd)
                if (dabs(scv1d(n,i)) > dlowval) then
                  scv_s(n,j,i-1) = real(scv1d(n,i))
                end if
                nnn = nnn + 1
              else
                tlef_s(n,j,i-1) = smissval
                ssw_s(n,j,i-1) = smissval
                rsw_s(n,j,i-1) = smissval
                rnos_s(n,j,i-1) = smissval
                scv_s(n,j,i-1) = smissval
              end if
            end do
            if ( nnn >= max0(ng/2,1) ) then
              tlef_o(j,i-1) = tlef_o(j,i-1)/real(nnn)
              ssw_o(j,i-1) = ssw_o(j,i-1)/real(nnn)
              rsw_o(j,i-1) = rsw_o(j,i-1)/real(nnn)
              rnos_o(j,i-1) = rnos_o(j,i-1)/real(nnn)*real(mmpd)
              scv_o(j,i-1) = scv_o(j,i-1)/real(nnn)
            else
              tlef_o(j,i-1) = smissval
              ssw_o(j,i-1) = smissval
              rsw_o(j,i-1) = smissval
              rnos_o(j,i-1) = smissval
              scv_o(j,i-1) = smissval
            end if
#else
#ifdef BAND
            drag_o(j,i-1) = 0.0
            q2m_o(j,i-1) = 0.0
            evpa_o(j,i-1) = 0.0
            sena_o(j,i-1) = 0.0
            do n = 1 , ng
              if ( ldoc1d(n,i) /= 0 ) then
                fracv = sigf(n,i)
                fracb = (d_one-veg1d(n,i))*(d_one-scvk(n,i))
                fracs = veg1d(n,i)*wt(n,i) + (d_one-veg1d(n,i))*scvk(n,i)
                facv = z2fra(n,i)/zlgveg(n,i)
                facb = z2fra(n,i)/zlglnd(n,i)
                facs = z2fra(n,i)/zlgsno(n,i)
                fact = fracv*facv + fracb*facb + fracs*facs
                q2m_1d(n,i) = qs1d(n,i) - delq1d(n,i)*fact
              else
                if ( iocnflx == 1 ) then
                  fact = z2fra(n,i)/zlgocn(n,i)
                  q2m_1d(n,i) = qs1d(n,i) - delq1d(n,i)*fact
                end if
              end if
              q2m_s(n,j,i-1) = real(q2m_1d(n,i))
              drag_s(n,j,i-1) = real(drag1d(n,i))
              evpa_s(n,j,i-1) = real(evpa2d(n,i,j)*mmpd)
              sena_s(n,j,i-1) = real(sena2d(n,i,j)*wpm2)
              tpr_s(n,j,i-1) = real((prnca2d(i,j)+prca2d(i,j))*mmpd)
              prcv_s(n,j,i-1) = real(prca2d(i,j)*mmpd)
              ps_s(n,j,i-1) = real(p1d(n,i)*0.01D0)

              q2m_o(j,i-1) = q2m_o(j,i-1) + real(q2m_1d(n,i))
              drag_o(j,i-1) = drag_o(j,i-1) + real(drag1d(n,i))
              evpa_o(j,i-1) = evpa_o(j,i-1) + real(evpa2d(n,i,j))
              sena_o(j,i-1) = sena_o(j,i-1) + real(sena2d(n,i,j))
            end do
            tpr_o(j,i-1) = real((prnca2d(i,j)+prca2d(i,j))*mmpd)
            q2m_o(j,i-1) = q2m_o(j,i-1)/real(ng)
            drag_o(j,i-1) = drag_o(j,i-1)/real(ng)
            evpa_o(j,i-1) = evpa_o(j,i-1)/real(ng)*real(mmpd)
            sena_o(j,i-1) = sena_o(j,i-1)/real(ng)*real(wpm2)
            flwa_o(j,i-1) = real(flwa2d(i,j)*wpm2)
            fswa_o(j,i-1) = real(fswa2d(i,j)*wpm2)
            flwd_o(j,i-1) = real(flwda2d(i,j)*wpm2)
            sina_o(j,i-1) = real(sina2d(i,j)*wpm2)
            prcv_o(j,i-1) = real(prca2d(i,j)*mmpd)
            ps_o(j,i-1) = real((sps2%ps(i,j)+r8pt)*d_10)
            zpbl_o(j,i-1) = real(sfsta%zpbl(i,j))

            tlef_o(j,i-1) = 0.0
            ssw_o(j,i-1) = 0.0
            rsw_o(j,i-1) = 0.0
            rnos_o(j,i-1) = 0.0
            scv_o(j,i-1) = 0.0
            nnn = 0
            do n = 1 , ng
              if ( ldoc1d(n,i,j) /= 0 ) then
                tlef_o(j,i-1) = tlef_o(j,i-1) + real(tlef1d(n,i))
                ssw_o(j,i-1) = ssw_o(j,i-1) + real(ssw1d(n,i))
                rsw_o(j,i-1) = rsw_o(j,i-1) + real(rsw1d(n,i))
                rnos_o(j,i-1) = rnos_o(j,i-1) + real(rnos2d(n,i,j))
                scv_o(j,i-1) = scv_o(j,i-1) + real(scv1d(n,i))
                tlef_s(n,j,i-1) = real(tlef1d(n,i))
                ssw_s(n,j,i-1) = real(ssw1d(n,i))
                rsw_s(n,j,i-1) = real(rsw1d(n,i))
                rnos_s(n,j,i-1) = real(rnos2d(n,i,j)*mmpd)
                scv_s(n,j,i-1) = real(scv1d(n,i))
                nnn = nnn + 1
              else
                tlef_s(n,j,i-1) = smissval
                ssw_s(n,j,i-1) = smissval
                rsw_s(n,j,i-1) = smissval
                rnos_s(n,j,i-1) = smissval
                scv_s(n,j,i-1) = smissval
              end if
            end do
            if ( nnn >= max0(ng/2,1) ) then
              tlef_o(j,i-1) = tlef_o(j,i-1)/real(nnn)
              ssw_o(j,i-1) = ssw_o(j,i-1)/real(nnn)
              rsw_o(j,i-1) = rsw_o(j,i-1)/real(nnn)
              rnos_o(j,i-1) = rnos_o(j,i-1)/real(nnn)*real(mmpd)
              scv_o(j,i-1) = scv_o(j,i-1)/real(nnn)
            else
              tlef_o(j,i-1) = smissval
              ssw_o(j,i-1) = smissval
              rsw_o(j,i-1) = smissval
              rnos_o(j,i-1) = smissval
              scv_o(j,i-1) = smissval
            end if
#else
            drag_o(j-1,i-1) = 0.0
            q2m_o(j-1,i-1) = 0.0
            evpa_o(j-1,i-1) = 0.0
            sena_o(j-1,i-1) = 0.0
            do n = 1 , ng
              if ( ldoc1d(n,i) /= 0 ) then
                fracv = sigf(n,i)
                fracb = (d_one-veg1d(n,i))*(d_one-scvk(n,i))
                fracs = veg1d(n,i)*wt(n,i) + (d_one-veg1d(n,i))*scvk(n,i)
                facv = z2fra(n,i)/zlgveg(n,i)
                facb = z2fra(n,i)/zlglnd(n,i)
                facs = z2fra(n,i)/zlgsno(n,i)
                fact = fracv*facv + fracb*facb + fracs*facs
                q2m_1d(n,i) = qs1d(n,i) - delq1d(n,i)*fact
              else 
                if ( iocnflx == 1 ) then
                  fact = z2fra(n,i)/zlgocn(n,i)
                  q2m_1d(n,i) = qs1d(n,i) - delq1d(n,i)*fact
                end if
              end if
              q2m_s(n,j-1,i-1) = real(q2m_1d(n,i))
              drag_s(n,j-1,i-1) = real(drag1d(n,i))
              evpa_s(n,j-1,i-1) = real(evpa2d(n,i,j)*mmpd)
              sena_s(n,j-1,i-1) = real(sena2d(n,i,j)*wpm2)
              tpr_s(n,j-1,i-1) = real((prnca2d(i,j)+prca2d(i,j))*mmpd)
              prcv_s(n,j-1,i-1) = real(prca2d(i,j)*mmpd)
              ps_s(n,j-1,i-1) = real(p1d(n,i)*0.01D0)

              q2m_o(j-1,i-1) = q2m_o(j-1,i-1) + real(q2m_1d(n,i))
              drag_o(j-1,i-1) = drag_o(j-1,i-1) + real(drag1d(n,i))
              evpa_o(j-1,i-1) = evpa_o(j-1,i-1) + real(evpa2d(n,i,j))
              sena_o(j-1,i-1) = sena_o(j-1,i-1) + real(sena2d(n,i,j))
            end do
            tpr_o(j-1,i-1) = real((prnca2d(i,j)+prca2d(i,j))*mmpd)
            q2m_o(j-1,i-1) = q2m_o(j-1,i-1)/real(ng)
            drag_o(j-1,i-1) = drag_o(j-1,i-1)/real(ng)
            evpa_o(j-1,i-1) = evpa_o(j-1,i-1)/real(ng)*real(mmpd)
            sena_o(j-1,i-1) = sena_o(j-1,i-1)/real(ng)*real(wpm2)
            flwa_o(j-1,i-1) = real(flwa2d(i,j)*wpm2)
            fswa_o(j-1,i-1) = real(fswa2d(i,j)*wpm2)
            flwd_o(j-1,i-1) = real(flwda2d(i,j)*wpm2)
            sina_o(j-1,i-1) = real(sina2d(i,j)*wpm2)
            prcv_o(j-1,i-1) = real(prca2d(i,j)*mmpd)
            ps_o(j-1,i-1) = real((sps2%ps(i,j)+r8pt)*d_10)
            zpbl_o(j-1,i-1) = real(sfsta%zpbl(i,j))

            tlef_o(j-1,i-1) = 0.0
            ssw_o(j-1,i-1) = 0.0
            rsw_o(j-1,i-1) = 0.0
            rnos_o(j-1,i-1) = 0.0
            scv_o(j-1,i-1) = 0.0
            nnn = 0
            do n = 1 , ng
              if ( ldoc1d(n,i) /= 0 ) then
                tlef_o(j-1,i-1) = tlef_o(j-1,i-1) + real(tlef1d(n,i))
                ssw_o(j-1,i-1) = ssw_o(j-1,i-1) + real(ssw1d(n,i))
                rsw_o(j-1,i-1) = rsw_o(j-1,i-1) + real(rsw1d(n,i))
                rnos_o(j-1,i-1) = rnos_o(j-1,i-1) + real(rnos2d(n,i,j))
                if (dabs(scv1d(n,i)) > dlowval) then
                  scv_o(j-1,i-1) = scv_o(j-1,i-1) + real(scv1d(n,i))
                end if
                tlef_s(n,j-1,i-1) = real(tlef1d(n,i))
                ssw_s(n,j-1,i-1) = real(ssw1d(n,i))
                rsw_s(n,j-1,i-1) = real(rsw1d(n,i))
                rnos_s(n,j-1,i-1) = real(rnos2d(n,i,j)*mmpd)
                if (dabs(scv1d(n,i)) > dlowval) then
                  scv_s(n,j-1,i-1) = real(scv1d(n,i))
                end if
                nnn = nnn + 1
              else
                tlef_s(n,j-1,i-1) = smissval
                ssw_s(n,j-1,i-1) = smissval
                rsw_s(n,j-1,i-1) = smissval
                rnos_s(n,j-1,i-1) = smissval
                scv_s(n,j-1,i-1) = smissval
              end if
            end do
            if ( nnn >= max0(ng/2,1) ) then
              tlef_o(j-1,i-1) = tlef_o(j-1,i-1)/real(nnn)
              ssw_o(j-1,i-1) = ssw_o(j-1,i-1)/real(nnn)
              rsw_o(j-1,i-1) = rsw_o(j-1,i-1)/real(nnn)
              rnos_o(j-1,i-1) = rnos_o(j-1,i-1)/real(nnn)*real(mmpd)
              scv_o(j-1,i-1) = scv_o(j-1,i-1)/real(nnn)
            else
              tlef_o(j-1,i-1) = smissval
              ssw_o(j-1,i-1) = smissval
              rsw_o(j-1,i-1) = smissval
              rnos_o(j-1,i-1) = smissval
              scv_o(j-1,i-1) = smissval
            end if
#endif
#endif
 
!           ******    reset accumulation arrays to zero
            do n = 1 , ng
              evpa2d(n,i,j) = d_zero
              rnos2d(n,i,j) = d_zero
              sena2d(n,i,j) = d_zero
            end do
            prnca2d(i,j) = d_zero
            prca2d(i,j) = d_zero
            flwa2d(i,j) = d_zero
            flwda2d(i,j) = d_zero
            fswa2d(i,j) = d_zero
            svga2d(i,j) = d_zero
            sina2d(i,j) = d_zero
          end do
        end if
!
      end if
!
      call time_end(subroutine_name,idindx)
      end subroutine interf
!
      subroutine albedov(j,iemiss)
 
      implicit none
!
      integer :: iemiss , j
      intent (in) iemiss , j
!
      real(8) :: age , albg , albgl , albgld , albgs , albgsd , albl ,  &
               & albld , albs , albsd , albzn , alwet , cf1 , cff ,     &
               & conn , cons , czeta , czf , dfalbl , dfalbs , dralbl , &
               & dralbs , fsol1 , fsol2 , sfac , sical0 , sical1 , sl , &
               & sl2 , sli , snal0 , snal1 , tdiff , tdiffs , wet
      real(8) , dimension(nnsg) :: albvl_s , albvs_s , aldifl_s ,       &
                                 & aldifs_s , aldirl_s , aldirs_s
      integer :: kolour , n , i
      character (len=50) :: subroutine_name='albedov'
      integer :: idindx=0
!
!     Albedo calculates fragmented albedos (direct and diffuse) in
!     wavelength regions split at 0.7um.
!
!     CM hands albedos to radiation package which computes
!     fsw1d(i) = net solar absorbed over full grid square
!     sabveg(i) = vegetation absorbed (full solar spectrum)
!     solis(i) = shortwave  solar incident
!
!     Here these are calculated at the end of albedo - they use only
!     direct albedos for now
!
!     in both versions :  lftemp uses sabveg
!     tgrund uses sabveg & fsw1d(i) to get
!     ground absorbed solar
!     photosynthesis uses solis - see subrouts
!     stomat and co2 (carbon)
!
!     For sea, sea-ice veg albedos are not set
!     these albedos are not treated as arrays here
!
!     (depuv/10.0)= the ratio of upper soil layer to total
!     root depth; used to compute "wet" for soil albedo
!
!     =================================================================
!     1. set initial parameters
!     =================================================================
!
!
      call time_begin(subroutine_name,idindx)
!     1.1 constants
!
!     Solar flux partitioned at wavelength of 0.7micr
      fsol1 = 0.5D0
      fsol2 = 0.5D0
!     Short and long wave albedo for new snow
      snal0 = 0.95D0
      snal1 = 0.65D0
!     Short and long wave albedo for sea ice
      sical0 = 0.6D0
      sical1 = 0.4D0
!
!     Desert seasonal albedo
!     Works for Sahara desert and generally northern emisphere
!     In souther emisphere only some points have this class
!
      if ( idesseas == 1 ) then
        if ( lmonth == 1 .or. lmonth == 2 .or. lmonth == 12 ) then
          solour(1) = 0.12D0
        endif        
        if ( lmonth == 3 .or. lmonth == 4 .or. lmonth == 5 ) then
          solour(1) = 0.15D0
        endif        
        if ( lmonth == 6 .or. lmonth == 7 .or. lmonth == 8) then
          solour(1) = 0.18D0
        endif        
        if ( lmonth == 9 .or. lmonth == 10 .or. lmonth == 11) then
          solour(1) = 0.15D0
        endif
      end if
!
!     In depth, wt is frac of grid square covered by snow;
!     depends on average snow depth, vegetation, etc.
!
      call depth
 
!     1.2 set pointers
!     ***************************************************
!     *    set n"x"k params here  in ccm but not needed *
!     ***************************************************
!
!     1.3  set default vegetation and albedo
!     do loop 50 in ccm not used here )
 
      do i = 2 , iym1
        czen(i) = dmax1(coszrs(i),d_zero)
        czeta = czen(i)
        do n = 1 , nnsg
          albgs = d_zero
          albgl = d_zero
          albgsd = d_zero
          albgld = d_zero
          albs = d_zero
          albl = d_zero
          albsd = d_zero
          albld = d_zero
 
          albvs_s(n) = d_zero
          albvl_s(n) = d_zero
 
!================================================================
!         2.   get albedo over land
!================================================================
!         can't use pointer "nalbk" here because not set - use nldock
!         instead tgb1d(i) used instead of tbelow
!
          if ( ldoc1d(n,i) == 2 ) then
            tdiffs = ts1d(n,i) - tzero
            tdiff = dmax1(tdiffs,d_zero)
            tdiffs = dmin1(tdiff,20.0D0)
            albgl = sical1 - 1.1D-2*tdiffs
            albgs = sical0 - 2.45D-2*tdiffs
            albg = fsol1*albgs + fsol2*albgl
            albgsd = albgs
            albgld = albgl
          else if ( ldoc1d(n,i) == 1 ) then
            sfac = d_one - fseas(tgb1d(n,i))
!           **********  ccm tests here on land mask for veg and soils
!c          data *** reduces albedo at low temps !!!!!should respond to
!c          moisture too the following card inactivated (commented out)
!c          (pat, 27 oct 86)
!           veg1d(i)=vegc(lveg(i))-seasf(lveg(i))*sfac
            albs = albvgs(lveg(n,i))
            albl = albvgl(lveg(n,i))
 
!----------------------------------------------------------------------
            if ( (lveg(n,i) < 12) .or. (lveg(n,i) > 15) ) then

!             2.1  bare soil albedos
!             (soil albedo depends on moisture)
              kolour = kolsol(lveg(n,i))
              wet = ssw1d(n,i)/depuv(lveg(n,i))
              alwet = dmax1((11.0D0-40.0D0*wet),d_zero)*0.01D0
              alwet = dmin1(alwet,solour(kolour))
              albg = solour(kolour) + alwet
!             if ((lveg(n,i) == 8)) albg=0.40      !Laura, cambiato il
!             DESERTO
              albgs = albg
              albgl = d_two*albg
!             higher nir albedos set diffuse albedo
              albgld = albgl
              albgsd = albgs
              albsd = albs
              albld = albl

!             Dec. 15   albzn=0.85D0+d_one/(d_one+d_10*czen(i))
!             Dec. 12, 2008
              albzn = d_one
!             Dec. 15, 2008

!             leafless hardwood canopy: no or inverse zen dep
              if ( lveg(n,i) == 5 .and. sfac < 0.1D0 ) albzn = d_one
!             multiply by zenith angle correction
              albs = albs*albzn
              albl = albl*albzn

!             albedo over vegetation after zenith angle corr
              albvs_s(n) = albs
              albvl_s(n) = albl

            else if ( lveg(n,i) == 12 ) then
 
!             2.2   permanent ice sheet
              albgs = 0.8D0
              albgsd = 0.8D0
              albgl = 0.55D0
              albgld = 0.55D0
            else
 
!             2.3  inland water, swamps, rice paddies etc.
              albg = 0.05D0/(czeta+0.15D0)
              albgs = albg
              albgsd = albg
              albgl = albg
              albgld = albg
            end if
 
          end if
! ===================================================================
!         4.  correct for snow cover
! ===================================================================
          if ( scv1d(n,i) > d_zero ) then
!           **********            snow albedo depends on  snow-age,
!           zenith angle, **********            and thickness of snow
 
!           **********            zenith angle set in zenitm
!           **********            snow albedoes for visible and ir
!           solar rad **********            visible albedo depends on
!           snow age **********            age gives reduction of
!           visible rad snow albedo **********              due to age
            cons = 0.2D0
            conn = 0.5D0
            age = (d_one-d_one/(d_one+sag1d(n,i)))
!           **********            sl helps control albedo zenith
!           dependence
            sl = d_two
            sli = d_one/sl
            sl2 = d_two*sl
!           **********            snal0= new snow albedo for vis rad,
!           sol zen le 6 **********            snal1= new snow albedo
!           for long-wave rad
            dfalbs = snal0*(d_one-cons*age)
!           **********            czf corrects albedo of new snow for
!           solar zenith
            cf1 = ((d_one+sli)/(d_one+sl2*czen(i))-sli)
            cff = dmax1(cf1,d_zero)
            czf = 0.4D0*cff*(d_one-dfalbs)
            dralbs = dfalbs + czf
            dfalbl = snal1*(d_one-conn*age)
            czf = 0.4D0*cff*(d_one-dfalbl)
            dralbl = dfalbl + czf
 
            if ( veg1d(n,i) > 0.001D0 ) then
!             **********            effective albedo over vegetation
!             with snow
              albl = (d_one-wt(n,i))*albl + dralbl*wt(n,i)
              albld = (d_one-wt(n,i))*albld + dfalbl*wt(n,i)
              albs = (d_one-wt(n,i))*albs + dralbs*wt(n,i)
              albsd = (d_one-wt(n,i))*albsd + dfalbs*wt(n,i)
            end if
 
!----------------------------------------------------------------------
!           4.1  compute albedo for snow on bare ground
!----------------------------------------------------------------------
            albgs = (d_one-scvk(n,i))*albgs + dralbs*scvk(n,i)
            albgl = (d_one-scvk(n,i))*albgl + dralbl*scvk(n,i)
            albgsd = (d_one-scvk(n,i))*albgsd + dfalbs*scvk(n,i)
            albgld = (d_one-scvk(n,i))*albgld + dfalbl*scvk(n,i)
          end if
 
!=====================================================================
!         5.  albedo over open ocean
!=====================================================================
          if ( ldoc1d(n,i) == 0 ) then
!           *********   ocean albedo depends on zenith angle
            if ( czeta >= d_zero ) then
!             **********   albedo independent of wavelength
              albg = 0.05D0/(czeta+0.15D0)
              albgs = albg
              albgl = albg
              albgsd = 0.08D0
              albgld = 0.08D0
            end if
          end if
 
!
!         ***************not part of albedo in the ccm ****************
!
          aldirs_s(n) = (d_one-veg1d(n,i))*albgs + veg1d(n,i)*albs
          aldirl_s(n) = (d_one-veg1d(n,i))*albgl + veg1d(n,i)*albl
          aldifs_s(n) = (d_one-veg1d(n,i))*albgsd + veg1d(n,i)*albsd
          aldifl_s(n) = (d_one-veg1d(n,i))*albgld + veg1d(n,i)*albld
        end do
        albvs(i) = albvs_s(1)
        albvl(i) = albvl_s(1)
        aldirs(i) = aldirs_s(1)
        aldirl(i) = aldirl_s(1)
        aldifs(i) = aldifs_s(1)
        aldifl(i) = aldifl_s(1)
        if ( iemiss == 1 ) emiss1d(i) = emiss2d(1,i,j)
        aldirs1d(1,i) = aldirs_s(1)
        aldifs1d(1,i) = aldifs_s(1)
        do n = 2 , nnsg
          albvs(i) = albvs(i) + albvs_s(n)
          albvl(i) = albvl(i) + albvl_s(n)
          aldirs(i) = aldirs(i) + aldirs_s(n)
          aldirl(i) = aldirl(i) + aldirl_s(n)
          aldifs(i) = aldifs(i) + aldifs_s(n)
          aldifl(i) = aldifl(i) + aldifl_s(n)
          if ( iemiss == 1 ) emiss1d(i) = emiss1d(i) + emiss2d(n,i,j)
          aldirs1d(n,i) = aldirs_s(n)
          aldifs1d(n,i) = aldifs_s(n)
        end do
        albvs(i) = albvs(i)/dble(nnsg)
        albvl(i) = albvl(i)/dble(nnsg)
        aldirs(i) = aldirs(i)/dble(nnsg)
        aldirl(i) = aldirl(i)/dble(nnsg)
        aldifs(i) = aldifs(i)/dble(nnsg)
        aldifl(i) = aldifl(i)/dble(nnsg)
        if ( iemiss == 1 ) emiss1d(i) = emiss1d(i)/dble(nnsg)
 
!       ******   fsw1d(i),sabveg(i),solis(i) computed in colrad
 
      end do
 
      call time_end(subroutine_name,idindx)

      contains

        function fseas(x)
          implicit none
          real(8) :: fseas
          real(8) , intent(in) :: x
          fseas = dmax1(d_zero,(d_one-0.0016D0*dmax1(298.0D0-x,d_zero)**d_two))
        end function fseas

      end subroutine albedov
!
      subroutine soilbc

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!   this subrout overwrites many of the soil constants
!   as a function of location(jlon,jlat)
!
      implicit none
!
      real(8) :: ck , dmax , dmin , dmnor , phi0 , tweak1
      integer :: itex , n , i
!
!     ================================================================
!     new soils data as a fn of texture make porosity, soil suction,
!     hydraul conduc, wilting frac variables rather than consts
!
!     for explanation of params set here see block data
!
!     relfc is the ratio of field capacity to saturated water content,
!     defined so the rate of gravitational drainage at field
!     capacity is assumed to be 2 mm/day (baver et al., 1972)
!     ===============================================================
!
!
      character (len=50) :: subroutine_name='soilbc'
      integer :: idindx=0
!
      call time_begin(subroutine_name,idindx)
 
      do i = 2 , iym1
        do n = 1 , nnsg
 
          if ( ldoc1d(n,i) /= 0 ) then
 
!           **********            lveg is set in subr. interf
            freza(lveg(n,i)) = 0.15D0*deprv(lveg(n,i))
            frezu(lveg(n,i)) = 0.15D0*depuv(lveg(n,i))
            itex = iexsol(lveg(n,i))
            texrat(n,i) = skrat(itex)
            porsl(n,i) = xmopor(itex)
            xkmx(n,i) = xmohyd(itex)
            bsw(n,i) = bee(itex)
            bfc(n,i) = 5.8D0 - bsw(n,i)*(0.8D0+0.12D0*(bsw(n,i)-d_four)* &
                      & dlog10(1.0D2*xkmx(n,i)))
 
            phi0 = xmosuc(itex)
            dmax = bsw(n,i)*phi0*xkmx(n,i)/porsl(n,i)
            dmin = 1.0D-3
            dmnor = 1550.0D0*dmin/dmax
            tweak1 = (bsw(n,i)*(bsw(n,i)-6.0D0)+10.3D0)  &
                   & /(bsw(n,i)*bsw(n,i)+40.0D0*bsw(n,i))
            ck = (d_one+dmnor)*tweak1*0.23D0/0.02356D0
            evmx0(n,i) = 1.02D0*dmax*ck/dsqrt(depuv(lveg(n,i))* &
                        & deprv(lveg(n,i)))
            gwmx0(n,i) = depuv(lveg(n,i))*porsl(n,i)
            gwmx1(n,i) = deprv(lveg(n,i))*porsl(n,i)
            gwmx2(n,i) = deptv(lveg(n,i))*porsl(n,i)
            wiltr(n,i) = xmowil(itex)
!           **********            force irrigated crop to be at field
!           capacity
            relfc(n,i) = xmofc(itex)
 
          end if
 
        end do
      end do
 
      call time_end(subroutine_name,idindx)
      end subroutine soilbc
!
      subroutine zenith(dec,alat,fjd,coszrs,frac,imax)
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
! ** used by radiation package in ccm to coszrs  *****not yet used
!                 here but could be called from albedo
!
!     this routine computes the zenith angle at each of 'imax'
!     equally-spaced longitudes at latitude 'alat', for a diurnally
!     varying sun.
!
!     input:  dec    - declination of sun in radians. computed in
!                      'solar'.
!             alat   - latitude in radians.
!             fjd    - fractional julian date in days.  the code conven-
!                      tion (see 'compjd') is that fjd=0 at noon  at
!                      greenwich meridian (1200 gmt), so the hour angle
!                      at longitude 'lon' is ha=lon+twopi*fjd.
!             imax   - number of equally spaced longitude values.
!                      it is assumed that i=1 and i=imax+1 are both
!                      at the greenwich meridian (lon=0).
!     output: coszrs - cos(z)/ where 'z' is the zenith angle at
!                      longitude 'i', latitude 'alat', and time 'fjd'.
!                           cos(z)=sin(alat)*sin(dec)+
!                                  cos(alat)*cos(dec)*cos(ha)
!                      it is assumed that the
!                      annual mean solar constant is used elsewhere
!                      in determining the solar flux.
!                      the 1/r**2 decrease of the solar flux appears
!                      in subroutine radiatn as eccf
!             frac   - not used in diurnal mode: set to 1.  in the
!                      average insolation mode, 'frac' is the daylight
!                      fraction at the point (see 'zenith'); to lowest
!                      order it should be independent of longitude.
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
      implicit none
!
      integer , intent (in) :: imax
      real(kind=8) , intent (in) :: alat , dec , fjd
      real(kind=8) , intent (out) , dimension(imax) :: coszrs , frac
!
      real(8) :: cc , cosz , dlon , ha , ss , tpifjd
      integer :: i
!
      character (len=50) :: subroutine_name='zenith'
      integer :: idindx=0
!
      call time_begin(subroutine_name,idindx)
!
!***********************************************************************
!
      ss = dsin(alat)*dsin(dec)
      cc = dcos(alat)*dcos(dec)
      dlon = twopi/imax
      tpifjd = twopi*fjd
      do i = 1 , imax
        frac(i) = d_one
        ha = (i-1)*dlon + tpifjd
!       if cosz is negative, the sun is below the horizon.
        cosz = dmax1(d_zero,ss+cc*dcos(ha))
        coszrs(i) = cosz
      end do
!
      call time_end(subroutine_name,idindx)
      end subroutine zenith
!
      subroutine slice1D(j)
 
      implicit none
!
      integer :: j
      intent (in) j
!
      real(8) :: amxtem , sfac
      integer :: n , i
! 
!     ********* For albedov
!     **********            This is taken from subroutine interf so that
!     **********            radiation can be called in tend (not vecbats).
      do i = 2 , iym1
        do n = 1 , nnsg
          ldoc1d(n,i) = ocld2d(n,i,j)
          sice1d(n,i) = sice2d(n,i,j)
          tgb1d(n,i) = tgb2d(n,i,j)
          ssw1d(n,i) = ssw2d(n,i,j)
          lveg(n,i) = veg2d1(n,i,j)
          oveg(n,i) = lveg(n,i)
          if ( ldoc1d(n,i) == 2 ) lveg(n,i) = 12
          amxtem = dmax1(298.0D0-tgb1d(n,i),d_zero)
          sfac = d_one - dmax1(d_zero,d_one-0.0016D0*amxtem**d_two)
          veg1d(n,i) = vegc(lveg(n,i)) - seasf(lveg(n,i))*sfac
          ts1d(n,i) = thx3d(i,kz,j)-6.5D-3*regrav*(ht1(n,i,j)- &
                      mddom%ht(i,j))
          scv1d(n,i) = scv2d(n,i,j)
          sag1d(n,i) = sag2d(n,i,j)
        end do
      end do

      end subroutine slice1D
!
      end module mod_vecbats
