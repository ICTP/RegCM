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

      use mod_constants
      use mod_message
      use mod_dynparam
      use mod_runparams
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

      contains

      subroutine vecbats(j,k,istart,iend,ng)

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
! Dummy arguments
!
      integer, intent(in) :: j , k , istart , iend , ng
!!
      character (len=50) :: subroutine_name='vecbats'
      integer :: idindx=0
!
      call time_begin(subroutine_name,idindx)

!---------------------------------------------------------------------

      call interf(1 , j , k , istart , iend , ng)
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
! Local variables
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
      data slmo/0.65 , 0.45 , 0.60 , 0.60 , 0.65 , 0.65 , 0.55 , 0.10 , &
         & 0.90 , 0.80 , 0.20 , 0.90 , 0.90 , 1.00 , 1.00 , 0.50 ,      &
         & 0.50 , 0.65 , 0.60 , 0.60/       !BATS land types
 
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
          pptnc(ill,jll) = 0.
          pptc(ill,jll) = 0.
!MM4      ist=nint(mddom%satbrt(ill,jll))
!MM4      if(ist.le.13)then
!MM4      ist=ist
!MM4      else
!MM4      ist=7
!MM4      end if
!MM4      veg2d(ill,jll)=vgtran(ist)
!eros     note:  when using bats dataset comment line above
          veg2d(ill,jll) = mddom%satbrt(ill,jll)
          if ( mddom%satbrt(ill,jll).gt.13.9 .and. &
               mddom%satbrt(ill,jll).lt.15.1 )  veg2d(ill,jll) = 0.
        end do
        do ill = 1 , iym1
          do k = 1 , nnsg
            veg2d1(k,ill,jll) = satbrt1(k,ill,jll)
            if ( satbrt1(k,ill,jll).gt.13.9 .and. &
                 satbrt1(k,ill,jll).lt.15.1 ) &
              veg2d1(k,ill,jll) = 0.
          end do
        end do
      end do

!     ******  initialize hostetler lake model
      if (lakemod.eq.1) call initlake
 
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
            if ( iseaice .ne. 1 ) then
              if ( veg2d1(k,ill,jll).lt.0.5 ) then
                ocld2d(k,ill,jll) = 0.
              else
                ocld2d(k,ill,jll) = 1.
              end if
            end if
            nlveg = nint(veg2d1(k,ill,jll))
            if (iseaice == 1) then
              if ( ocld2d(k,ill,jll) > 1.5 ) nlveg = 12
            end if
            if ( nlveg.eq.0 ) then
              nlveg = 15
            else
              nlveg = nlveg
            end if
!sol        itex=int(text2d(k,ill,jll))
            itex = iexsol(nlveg)
            tg2d(k,ill,jll) = sts2%tg(ill,jll)
            tgb2d(k,ill,jll) = sts2%tg(ill,jll)
            taf2d(k,ill,jll) = sts2%tg(ill,jll)
            tlef2d(k,ill,jll) = sts2%tg(ill,jll)

!           ******  initialize soil moisture in the 3 layers
            is = nint(satbrt1(k,ill,jll))
            swt2d(k,ill,jll) = deptv(nlveg)*xmopor(itex)*slmo(is)
            srw2d(k,ill,jll) = deprv(nlveg)*xmopor(itex)*slmo(is)
            ssw2d(k,ill,jll) = depuv(nlveg)*xmopor(itex)*slmo(is)

            dew2d(k,ill,jll) = 0.
            sag2d(k,ill,jll) = 0.
            scv2d(k,ill,jll) = max(snowc(k,ill,jll),0.D0)
            sice2d(k,ill,jll) = 0.
            gwet2d(k,ill,jll) = 0.5
            sena2d(k,ill,jll) = 0.
            evpa2d(k,ill,jll) = 0.
            rnos2d(k,ill,jll) = 0.
            rno2d(k,ill,jll) = 0.
            if ( sice2d(k,ill,jll).gt.0. ) then
              ocld2d(k,ill,jll) = 2.
            end if
            ircp2d(k,ill,jll) = 0.
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
          fsw2d(ill,jll) = 0.
          flw2d(ill,jll) = 0.
          sabv2d(ill,jll) = 0.
          sol2d(ill,jll) = 0.
          fswa2d(ill,jll) = 0.
          flwa2d(ill,jll) = 0.
          prca2d(ill,jll) = 0.
          prnca2d(ill,jll) = 0.
          svga2d(ill,jll) = 0.
          sina2d(ill,jll) = 0.
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
! Dummy arguments
!
      integer , intent (in) :: ivers , j , k , istart , iend , ng
!
! Local variables
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
 
      if ( ivers.eq.1 ) then ! regcm2d --> bats

        do i = istart, iend
          do n = 1 , ng
            p1d0(n,i) = (sps2%ps(i,j)+r8pt)*1000.
            z1d(n,i) = za(i,k,j)
            ts1d0(n,i) = thx3d(i,k,j)
            qs1d0(n,i) = qvb3d(i,k,j)/(1.+qvb3d(i,k,j))
            qs1d(n,i) = qs1d0(n,i)
 
            hl = lh0 - lh1*(ts1d0(n,i)-tzero)
            satvp = lsvp1*dexp(lsvp2*hl*(1./tzero-1./ts1d0(n,i)))
            rh0 = max(qs1d0(n,i)/(ep2*satvp/(p1d0(n,i)*0.01-satvp)), &
                & 0.D0)
 
            ts1d(n,i) = ts1d0(n,i) - lrate*rgti*(ht1(n,i,j)- &
                                                 mddom%ht(i,j))
            p1d(n,i) = p1d0(n,i)*(ts1d(n,i)/ts1d0(n,i))
 
            hl = lh0 - lh1*(ts1d(n,i)-tzero)
            satvp = lsvp1*dexp(lsvp2*hl*(1./tzero-1./ts1d(n,i)))
            qs1d(n,i) = max(rh0*ep2*satvp/(p1d(n,i)*0.01-satvp),0.D0)
 
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
            lveg(n,i) = nint(veg2d1(n,i,j))
            if (iseaice == 1) then
              if (ocld2d(n,i,j) > 1.5) lveg(n,i) = 12
            end if
            amxtem = max(298.-tgb1d(n,i),0.D0)
            sfac = 1. - max(0.D0,1.-0.0016*amxtem**2)
            if ( lveg(n,i).eq.0 ) then
              veg1d(n,i) = 0.
            else
              veg1d(n,i) = vegc(lveg(n,i)) - seasf(lveg(n,i))*sfac
            end if
            emiss_1d(n,i) = emiss2d(n,i,j)
          end do
 
          rh0 = 0.0D0
          do n = 1 , ng
            rh0 = rh0 + (qs1d(n,i)-qs1d0(n,i))
          end do
          rh0 = rh0/ng
          do n = 1 , ng
            qs1d(n,i) = max(qs1d(n,i)-rh0,0.0D0)
          end do
 
          us1d(i) = ubx3d(i,k,j)
          vs1d(i) = vbx3d(i,k,j)
          fsw1d(i) = fsw2d(i,j)
          flw1d(i) = flw2d(i,j)
          solis(i) = sol2d(i,j)
          sabveg(i) = sabv2d(i,j)
          solvt = solvd2d(i,j) + solvs2d(i,j)
          if ( solvt.gt.0.0 ) then
            fracd(i) = solvd2d(i,j)/solvt
          else
            fracd(i) = 0.2
          end if
          czen(i) = max(coszrs(i),0.D0)
        end do
 
      else if ( ivers.eq.2 ) then ! bats --> regcm2d
 
        do i = istart, iend
          sfsta%uvdrag(i,j) = 0.0
          sfsta%hfx(i,j) = 0.0
          sfsta%qfx(i,j) = 0.0
          sts2%tg(i,j) = 0.0
          sts1%tg(i,j) = 0.0
          sfsta%tgbb(i,j) = 0.0
!chem2
          ssw2da(i,j) = 0.0
          sdeltk2d(i,j) = 0.0
          sdelqk2d(i,j) = 0.0
          sfracv2d(i,j) = 0.0
          sfracb2d(i,j) = 0.0
          sfracs2d(i,j) = 0.0
          svegfrac2d(i,j) = 0.0
!chem2_
          do n = 1 , ng
            sfsta%uvdrag(i,j) = sfsta%uvdrag(i,j) + drag1d(n,i)
            sfsta%hfx(i,j) = sfsta%hfx(i,j) + sent1d(n,i)
            sfsta%qfx(i,j) = sfsta%qfx(i,j) + evpr1d(n,i)
            sts2%tg(i,j) = sts2%tg(i,j) + tg1d(n,i)
            sts1%tg(i,j) = sts1%tg(i,j) + tg1d(n,i)
!chem2
            ssw2da(i,j) = ssw2da(i,j) + ssw1d(n,i)
            sdeltk2d(i,j) = sdeltk2d(i,j) + delt1d(n,i)
            sdelqk2d(i,j) = sdelqk2d(i,j) + delq1d(n,i)
            sfracv2d(i,j) = sfracv2d(i,j) + sigf(n,i)
            sfracb2d(i,j) = sfracb2d(i,j) + (1.-veg1d(n,i))             &
                          & *(1.-scvk(n,i))
            sfracs2d(i,j) = sfracs2d(i,j) + veg1d(n,i)*wt(n,i)          &
                          & + (1.-veg1d(n,i))*scvk(n,i)
            svegfrac2d(i,j) = svegfrac2d(i,j) + veg1d(n,i)
!chem2_
            if ( iocnflx.eq.1 .or.                                      &
               & (iocnflx.eq.2 .and. ocld2d(n,i,j).ge.0.5 ) ) then
              sfsta%tgbb(i,j) = sfsta%tgbb(i,j)                         &
                        & + ((1.-veg1d(n,i))*tg1d(n,i)**4+veg1d(n,i)    &
                        & *tlef1d(n,i)**4)**0.25
            else
              sfsta%tgbb(i,j) = sfsta%tgbb(i,j) + tg1d(n,i)
            end if
            if ( ocld2d(n,i,j).lt.0.5 ) then
              ssw1d(n,i)  = -1.D34
              rsw1d(n,i)  = -1.D34
              tsw1d(n,i)  = -1.D34
              rno1d(n,i)  = -1.D34
              rnos1d(n,i) = -1.D34
              scv1d(n,i)  = -1.D34
            end if
          end do
          sfsta%uvdrag(i,j) = sfsta%uvdrag(i,j)/float(ng)
          sfsta%hfx(i,j) = sfsta%hfx(i,j)/float(ng)
          sfsta%qfx(i,j) = sfsta%qfx(i,j)/float(ng)
          sts2%tg(i,j) = sts2%tg(i,j)/float(ng)
          sts1%tg(i,j) = sts1%tg(i,j)/float(ng)
          sfsta%tgbb(i,j) = sfsta%tgbb(i,j)/float(ng)
!chem2
          ssw2da(i,j) = ssw2da(i,j)/float(ng)
          sdeltk2d(i,j) = sdeltk2d(i,j)/float(ng)
          sdelqk2d(i,j) = sdelqk2d(i,j)/float(ng)
          sfracv2d(i,j) = sfracv2d(i,j)/float(ng)
          sfracb2d(i,j) = sfracb2d(i,j)/float(ng)
          sfracs2d(i,j) = sfracs2d(i,j)/float(ng)
          svegfrac2d(i,j) = svegfrac2d(i,j)/float(ng)
!chem2_
          do n = 1 , ng
            snowc(n,i,j) = scv1d(n,i)
            tg2d(n,i,j) = tg1d(n,i)
            tgb2d(n,i,j) = tgb1d(n,i)
            taf2d(n,i,j) = taf1d(n,i)
            tlef2d(n,i,j) = tlef1d(n,i)
            swt2d(n,i,j) = tsw1d(n,i)
            srw2d(n,i,j) = rsw1d(n,i)
            ssw2d(n,i,j) = ssw1d(n,i)
            dew2d(n,i,j) = ldew1d(n,i)
            sag2d(n,i,j) = sag1d(n,i)
            scv2d(n,i,j) = scv1d(n,i)
            sice2d(n,i,j) = sice1d(n,i)
            gwet2d(n,i,j) = gwet1d(n,i)
            ocld2d(n,i,j) = ldoc1d(n,i)
            ircp2d(n,i,j) = ircp1d(n,i)
            evpa2d(n,i,j) = evpa2d(n,i,j) + dtbat*evpr1d(n,i)
            sena2d(n,i,j) = sena2d(n,i,j) + dtbat*sent1d(n,i)
            if ( rnos2d(n,i,j).gt.-1.D10 .and. rnos1d(n,i).gt.-1.D10 )  &
               & then
              rnos2d(n,i,j) = rnos2d(n,i,j) + rnos1d(n,i)/tau1*dtbat
            else
              rnos2d(n,i,j) = -1.D34
            end if
            if ( rno2d(n,i,j).gt.-1.D10 .and. rnos1d(n,i)               &
               & .gt.-1.D10 .and. rno1d(n,i).gt.-1.D10 ) then
              rno2d(n,i,j) = rno2d(n,i,j) + (rno1d(n,i)-rnos1d(n,i))    &
                           & /tau1*dtbat
            else
              rno2d(n,i,j) = -1.D34
            end if
          end do
!
!         quantities stored on 2d surface array for bats use only
!
          prca2d(i,j) = prca2d(i,j) + dtbat*pptc(i,j)
          prnca2d(i,j) = prnca2d(i,j) + dtbat*pptnc(i,j)
          if ( prnca2d(i,j) < 1D-30 ) prnca2d(i,j) = 0.0
          if ( prca2d(i,j) < 1D-30 ) prca2d(i,j) = 0.0
          flwa2d(i,j) = flwa2d(i,j) + dtbat*flw1d(i)
          flwda2d(i,j) = flwda2d(i,j) + dtbat*flwd2d(i,j)
          fswa2d(i,j) = fswa2d(i,j) + dtbat*fsw1d(i)
          svga2d(i,j) = svga2d(i,j) + dtbat*sabveg(i)
          sina2d(i,j) = sina2d(i,j) + dtbat*sinc2d(i,j)
          pptnc(i,j) = 0.
          pptc(i,j) = 0.
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
            if ( ocld2d(n,i,j).ge.0.5 ) then
              fracv = sigf(n,i)
              fracb = (1.-veg1d(n,i))*(1.-scvk(n,i))
              fracs = veg1d(n,i)*wt(n,i) + (1.-veg1d(n,i))*scvk(n,i)
              facv = dlog(z1(n,i)/2.)/dlog(z1(n,i)/rough(lveg(n,i)))
              facb = dlog(z1(n,i)/2.)/dlog(z1(n,i)/zlnd)
              facs = dlog(z1(n,i)/2.)/dlog(z1(n,i)/zsno)
              fact = fracv*facv + fracb*facb + fracs*facs
              facv = dlog(z1(n,i)/10.)/dlog(z1(n,i)/rough(lveg(n,i)))
              facb = dlog(z1(n,i)/10.)/dlog(z1(n,i)/zlnd)
              facs = dlog(z1(n,i)/10.)/dlog(z1(n,i)/zsno)
              factuv = fracv*facv + fracb*facb + fracs*facs
              u10m1d(n,i) = us1d(i)*(1.-factuv)
              v10m1d(n,i) = vs1d(i)*(1.-factuv)
              t2m_1d(n,i) = ts1d(n,i) - delt1d(n,i)*fact
            else if ( iocnflx.eq.1 ) then
              fact = dlog(z1(n,i)/2.)/dlog(z1(n,i)/zoce)
              factuv = dlog(z1(n,i)/10.)/dlog(z1(n,i)/zoce)
              u10m1d(n,i) = us1d(i)*(1.-factuv)
              v10m1d(n,i) = vs1d(i)*(1.-factuv)
              t2m_1d(n,i) = ts1d(n,i) - delt1d(n,i)*fact
            else
            end if
            tg_s(n,j,i-1) = tg1d(n,i)
            u10m_s(n,j,i-1) = u10m1d(n,i)
            v10m_s(n,j,i-1) = v10m1d(n,i)
            t2m_s(n,j,i-1) = t2m_1d(n,i)
 
            u10m_o(j,i-1) = u10m_o(j,i-1) + u10m1d(n,i)
            v10m_o(j,i-1) = v10m_o(j,i-1) + v10m1d(n,i)
            t2m_o(j,i-1) = t2m_o(j,i-1) + t2m_1d(n,i)
            tg_o(j,i-1) = tg_o(j,i-1) + tg1d(n,i)
            aldirs_o(j,i-1) = aldirs_o(j,i-1) + aldirs1d(n,i)
            aldifs_o(j,i-1) = aldifs_o(j,i-1) + aldifs1d(n,i)
          end do
          u10m_o(j,i-1) = u10m_o(j,i-1)/float(ng)
          v10m_o(j,i-1) = v10m_o(j,i-1)/float(ng)
          t2m_o(j,i-1) = t2m_o(j,i-1)/float(ng)
          tg_o(j,i-1) = tg_o(j,i-1)/float(ng)
          aldirs_o(j,i-1) = aldirs_o(j,i-1)/float(ng)
          aldifs_o(j,i-1) = aldifs_o(j,i-1)/float(ng)
 
          tgmx_o(j,i-1) = max(tgmx_o(j,i-1),tg_o(j,i-1))
          tgmn_o(j,i-1) = min(tgmn_o(j,i-1),tg_o(j,i-1))
          t2mx_o(j,i-1) = max(t2mx_o(j,i-1),t2m_o(j,i-1))
          t2mn_o(j,i-1) = min(t2mn_o(j,i-1),t2m_o(j,i-1))
          w10x_o(j,i-1) = max(w10x_o(j,i-1),sqrt(u10m_o(j,i-1)**2+    &
                        & v10m_o(j,i-1)**2))
          real_4 = (sps2%ps(i,j)+r8pt)*10.
          psmn_o(j,i-1) = min(psmn_o(j,i-1),real_4)

#else
#ifdef BAND
          u10m_o(j,i-1) = 0.0
          v10m_o(j,i-1) = 0.0
          tg_o(j,i-1) = 0.0
          t2m_o(j,i-1) = 0.0
          aldirs_o(j,i-1) = 0.0
          aldifs_o(j,i-1) = 0.0
          do n = 1 , ng
            if ( ocld2d(n,i,j).ge.0.5 ) then
              fracv = sigf(n,i)
              fracb = (1.-veg1d(n,i))*(1.-scvk(n,i))
              fracs = veg1d(n,i)*wt(n,i) + (1.-veg1d(n,i))*scvk(n,i)
              facv = dlog(z1(n,i)/2.)/dlog(z1(n,i)/rough(lveg(n,i)))
              facb = dlog(z1(n,i)/2.)/dlog(z1(n,i)/zlnd)
              facs = dlog(z1(n,i)/2.)/dlog(z1(n,i)/zsno)
              fact = fracv*facv + fracb*facb + fracs*facs
              facv = dlog(z1(n,i)/10.)/dlog(z1(n,i)/rough(lveg(n,i)))
              facb = dlog(z1(n,i)/10.)/dlog(z1(n,i)/zlnd)
              facs = dlog(z1(n,i)/10.)/dlog(z1(n,i)/zsno)
              factuv = fracv*facv + fracb*facb + fracs*facs
              u10m1d(n,i) = us1d(i)*(1.-factuv)
              v10m1d(n,i) = vs1d(i)*(1.-factuv)
              t2m_1d(n,i) = ts1d(n,i) - delt1d(n,i)*fact
            else if ( iocnflx.eq.1 ) then
              fact = dlog(z1(n,i)/2.)/dlog(z1(n,i)/zoce)
              factuv = dlog(z1(n,i)/10.)/dlog(z1(n,i)/zoce)
              u10m1d(n,i) = us1d(i)*(1.-factuv)
              v10m1d(n,i) = vs1d(i)*(1.-factuv)
              t2m_1d(n,i) = ts1d(n,i) - delt1d(n,i)*fact
            else
            end if
            tg_s(n,j,i-1) = tg1d(n,i)
            u10m_s(n,j,i-1) = u10m1d(n,i)
            v10m_s(n,j,i-1) = v10m1d(n,i)
            t2m_s(n,j,i-1) = t2m_1d(n,i)

            u10m_o(j,i-1) = u10m_o(j,i-1) + u10m1d(n,i)
            v10m_o(j,i-1) = v10m_o(j,i-1) + v10m1d(n,i)
            t2m_o(j,i-1) = t2m_o(j,i-1) + t2m_1d(n,i)
            tg_o(j,i-1) = tg_o(j,i-1) + tg1d(n,i)
            aldirs_o(j,i-1) = aldirs_o(j,i-1) + aldirs1d(n,i)
            aldifs_o(j,i-1) = aldifs_o(j,i-1) + aldifs1d(n,i)
          end do
          u10m_o(j,i-1) = u10m_o(j,i-1)/float(ng)
          v10m_o(j,i-1) = v10m_o(j,i-1)/float(ng)
          t2m_o(j,i-1) = t2m_o(j,i-1)/float(ng)
          tg_o(j,i-1) = tg_o(j,i-1)/float(ng)
          aldirs_o(j,i-1) = aldirs_o(j,i-1)/float(ng)
          aldifs_o(j,i-1) = aldifs_o(j,i-1)/float(ng)
          tgmx_o(j,i-1) = max(tgmx_o(j,i-1),tg_o(j,i-1))
          tgmn_o(j,i-1) = min(tgmn_o(j,i-1),tg_o(j,i-1))
          t2mx_o(j,i-1) = max(t2mx_o(j,i-1),t2m_o(j,i-1))
          t2mn_o(j,i-1) = min(t2mn_o(j,i-1),t2m_o(j,i-1))
          w10x_o(j,i-1) = max(w10x_o(j,i-1),sqrt(u10m_o(j,i-1)**&
                          & 2+v10m_o(j,i-1)**2))
          real_4 = (sps2%ps(i,j)+r8pt)*10.
          psmn_o(j,i-1) = min(psmn_o(j,i-1),real_4)
#else
          u10m_o(j-1,i-1) = 0.0
          v10m_o(j-1,i-1) = 0.0
          tg_o(j-1,i-1) = 0.0
          t2m_o(j-1,i-1) = 0.0
          aldirs_o(j-1,i-1) = 0.0
          aldifs_o(j-1,i-1) = 0.0
          do n = 1 , ng
            if ( ocld2d(n,i,j).ge.0.5 ) then
              fracv = sigf(n,i)
              fracb = (1.-veg1d(n,i))*(1.-scvk(n,i))
              fracs = veg1d(n,i)*wt(n,i) + (1.-veg1d(n,i))*scvk(n,i)
              facv = dlog(z1(n,i)/2.)/dlog(z1(n,i)/rough(lveg(n,i)))
              facb = dlog(z1(n,i)/2.)/dlog(z1(n,i)/zlnd)
              facs = dlog(z1(n,i)/2.)/dlog(z1(n,i)/zsno)
              fact = fracv*facv + fracb*facb + fracs*facs
              facv = dlog(z1(n,i)/10.)/dlog(z1(n,i)/rough(lveg(n,i)))
              facb = dlog(z1(n,i)/10.)/dlog(z1(n,i)/zlnd)
              facs = dlog(z1(n,i)/10.)/dlog(z1(n,i)/zsno)
              factuv = fracv*facv + fracb*facb + fracs*facs
              u10m1d(n,i) = us1d(i)*(1.-factuv)
              v10m1d(n,i) = vs1d(i)*(1.-factuv)
              t2m_1d(n,i) = ts1d(n,i) - delt1d(n,i)*fact
            else if ( iocnflx.eq.1 ) then
              fact = dlog(z1(n,i)/2.)/dlog(z1(n,i)/zoce)
              factuv = dlog(z1(n,i)/10.)/dlog(z1(n,i)/zoce)
              u10m1d(n,i) = us1d(i)*(1.-factuv)
              v10m1d(n,i) = vs1d(i)*(1.-factuv)
              t2m_1d(n,i) = ts1d(n,i) - delt1d(n,i)*fact
            else
            end if
            tg_s(n,j-1,i-1) = tg1d(n,i)
            u10m_s(n,j-1,i-1) = u10m1d(n,i)
            v10m_s(n,j-1,i-1) = v10m1d(n,i)
            t2m_s(n,j-1,i-1) = t2m_1d(n,i)

            u10m_o(j-1,i-1) = u10m_o(j-1,i-1) + u10m1d(n,i)
            v10m_o(j-1,i-1) = v10m_o(j-1,i-1) + v10m1d(n,i)
            t2m_o(j-1,i-1) = t2m_o(j-1,i-1) + t2m_1d(n,i)
            tg_o(j-1,i-1) = tg_o(j-1,i-1) + tg1d(n,i)
            aldirs_o(j-1,i-1) = aldirs_o(j-1,i-1) + aldirs1d(n,i)
            aldifs_o(j-1,i-1) = aldifs_o(j-1,i-1) + aldifs1d(n,i)
          end do
          u10m_o(j-1,i-1) = u10m_o(j-1,i-1)/float(ng)
          v10m_o(j-1,i-1) = v10m_o(j-1,i-1)/float(ng)
          t2m_o(j-1,i-1) = t2m_o(j-1,i-1)/float(ng)
          tg_o(j-1,i-1) = tg_o(j-1,i-1)/float(ng)
          aldirs_o(j-1,i-1) = aldirs_o(j-1,i-1)/float(ng)
          aldifs_o(j-1,i-1) = aldifs_o(j-1,i-1)/float(ng)
          tgmx_o(j-1,i-1) = max(tgmx_o(j-1,i-1),tg_o(j-1,i-1))
          tgmn_o(j-1,i-1) = min(tgmn_o(j-1,i-1),tg_o(j-1,i-1))
          t2mx_o(j-1,i-1) = max(t2mx_o(j-1,i-1),t2m_o(j-1,i-1))
          t2mn_o(j-1,i-1) = min(t2mn_o(j-1,i-1),t2m_o(j-1,i-1))
          w10x_o(j-1,i-1) = max(w10x_o(j-1,i-1),sqrt(u10m_o(j-1,i-1)**&
                          & 2+v10m_o(j-1,i-1)**2))
          real_4 = (sps2%ps(i,j)+r8pt)*10.
          psmn_o(j-1,i-1) = min(psmn_o(j-1,i-1),real_4)
#endif
#endif
        end do

        if ( mod(ntime+nint(dtmin*60.),kbats).eq.0 .or.                 &
           & (jyear.eq.jyearr .and. ktau.eq.ktaur) ) then
          if ( jyear.eq.jyear0 .and. ktau.le.1 ) then
            mmpd = 86400./dtbat
            wpm2 = 1./dtbat
          else if ( jyear.eq.jyear0 .and. dble(ktau*dtmin)              &
                  & .le.batfrq*60.+0.01 ) then
            mmpd = 24./(batfrq-dtmin/60.)
            wpm2 = 1./((batfrq-dtmin/60.)*3600.)
          else
            mmpd = 24./batfrq
            wpm2 = 1./(batfrq*3600.)
          end if
          do i = istart, iend
#ifdef MPP1
            drag_o(j,i-1) = 0.0
            q2m_o(j,i-1) = 0.0
            evpa_o(j,i-1) = 0.0
            sena_o(j,i-1) = 0.0
            do n = 1 , ng
              if ( ocld2d(n,i,j).ge.0.5 ) then
                fracv = sigf(n,i)
                fracb = (1.-veg1d(n,i))*(1.-scvk(n,i))
                fracs = veg1d(n,i)*wt(n,i) + (1.-veg1d(n,i))*scvk(n,i)
                facv = dlog(z1(n,i)/2.)/dlog(z1(n,i)/rough(lveg(n,i)))
                facb = dlog(z1(n,i)/2.)/dlog(z1(n,i)/zlnd)
                facs = dlog(z1(n,i)/2.)/dlog(z1(n,i)/zsno)
                fact = fracv*facv + fracb*facb + fracs*facs
                q2m_1d(n,i) = qs1d(n,i) - delq1d(n,i)*fact
              else if ( iocnflx.eq.1 ) then
                fact = dlog(z1(n,i)/2.)/dlog(z1(n,i)/zoce)
                q2m_1d(n,i) = qs1d(n,i) - delq1d(n,i)*fact
              else
              end if
              q2m_s(n,j,i-1) = q2m_1d(n,i)
              drag_s(n,j,i-1) = drag1d(n,i)
              evpa_s(n,j,i-1) = evpa2d(n,i,j)*mmpd
              sena_s(n,j,i-1) = sena2d(n,i,j)*wpm2
              tpr_s(n,j,i-1) = (prnca2d(i,j)+prca2d(i,j))*mmpd
              prcv_s(n,j,i-1) = prca2d(i,j)*mmpd
              ps_s(n,j,i-1) = p1d(n,i)*0.01
 
              q2m_o(j,i-1) = q2m_o(j,i-1) + q2m_1d(n,i)
              drag_o(j,i-1) = drag_o(j,i-1) + drag1d(n,i)
              evpa_o(j,i-1) = evpa_o(j,i-1) + evpa2d(n,i,j)
              sena_o(j,i-1) = sena_o(j,i-1) + sena2d(n,i,j)
            end do
            tpr_o(j,i-1) = (prnca2d(i,j)+prca2d(i,j))*mmpd
            q2m_o(j,i-1) = q2m_o(j,i-1)/float(ng)
            drag_o(j,i-1) = drag_o(j,i-1)/float(ng)
            evpa_o(j,i-1) = evpa_o(j,i-1)/float(ng)*mmpd
            sena_o(j,i-1) = sena_o(j,i-1)/float(ng)*wpm2
            flwa_o(j,i-1) = flwa2d(i,j)*wpm2
            fswa_o(j,i-1) = fswa2d(i,j)*wpm2
            flwd_o(j,i-1) = flwda2d(i,j)*wpm2
            sina_o(j,i-1) = sina2d(i,j)*wpm2
            prcv_o(j,i-1) = prca2d(i,j)*mmpd
            ps_o(j,i-1) = (sps2%ps(i,j)+r8pt)*10.
            zpbl_o(j,i-1) = sfsta%zpbl(i,j)
 
            tlef_o(j,i-1) = 0.0
            ssw_o(j,i-1) = 0.0
            rsw_o(j,i-1) = 0.0
            rnos_o(j,i-1) = 0.0
            scv_o(j,i-1) = 0.0
            nnn = 0
            do n = 1 , ng
              if ( ocld2d(n,i,j).ge.0.5 ) then
                tlef_o(j,i-1) = tlef_o(j,i-1) + tlef1d(n,i)
                ssw_o(j,i-1) = ssw_o(j,i-1) + ssw1d(n,i)
                rsw_o(j,i-1) = rsw_o(j,i-1) + rsw1d(n,i)
                rnos_o(j,i-1) = rnos_o(j,i-1) + rnos2d(n,i,j)
                if (abs(scv1d(n,i)) > 1D-30) then
                  scv_o(j,i-1) = scv_o(j,i-1) + scv1d(n,i)
                end if
                tlef_s(n,j,i-1) = tlef1d(n,i)
                ssw_s(n,j,i-1) = ssw1d(n,i)
                rsw_s(n,j,i-1) = rsw1d(n,i)
                rnos_s(n,j,i-1) = rnos2d(n,i,j)*mmpd
                if (abs(scv1d(n,i)) > 1D-30) then
                  scv_s(n,j,i-1) = scv1d(n,i)
                end if
                nnn = nnn + 1
              else
                tlef_s(n,j,i-1) = -1.D34
                ssw_s(n,j,i-1) = -1.D34
                rsw_s(n,j,i-1) = -1.D34
                rnos_s(n,j,i-1) = -1.D34
                scv_s(n,j,i-1) = -1.D34
              end if
            end do
            if ( nnn.ge.max0(ng/2,1) ) then
              tlef_o(j,i-1) = tlef_o(j,i-1)/float(nnn)
              ssw_o(j,i-1) = ssw_o(j,i-1)/float(nnn)
              rsw_o(j,i-1) = rsw_o(j,i-1)/float(nnn)
              rnos_o(j,i-1) = rnos_o(j,i-1)/float(nnn)*mmpd
              scv_o(j,i-1) = scv_o(j,i-1)/float(nnn)
            else
              tlef_o(j,i-1) = -1.D34
              ssw_o(j,i-1) = -1.D34
              rsw_o(j,i-1) = -1.D34
              rnos_o(j,i-1) = -1.D34
              scv_o(j,i-1) = -1.D34
            end if
#else
#ifdef BAND
            drag_o(j,i-1) = 0.0
            q2m_o(j,i-1) = 0.0
            evpa_o(j,i-1) = 0.0
            sena_o(j,i-1) = 0.0
            do n = 1 , ng
              if ( ocld2d(n,i,j).ge.0.5 ) then
                fracv = sigf(n,i)
                fracb = (1.-veg1d(n,i))*(1.-scvk(n,i))
                fracs = veg1d(n,i)*wt(n,i) + (1.-veg1d(n,i))*scvk(n,i)
                facv = dlog(z1(n,i)/2.)/dlog(z1(n,i)/rough(lveg(n,i)))
                facb = dlog(z1(n,i)/2.)/dlog(z1(n,i)/zlnd)
                facs = dlog(z1(n,i)/2.)/dlog(z1(n,i)/zsno)
                fact = fracv*facv + fracb*facb + fracs*facs
                q2m_1d(n,i) = qs1d(n,i) - delq1d(n,i)*fact
              else if ( iocnflx.eq.1 ) then
                fact = dlog(z1(n,i)/2.)/dlog(z1(n,i)/zoce)
                q2m_1d(n,i) = qs1d(n,i) - delq1d(n,i)*fact
              else
              end if
              q2m_s(n,j,i-1) = q2m_1d(n,i)
              drag_s(n,j,i-1) = drag1d(n,i)
              evpa_s(n,j,i-1) = evpa2d(n,i,j)*mmpd
              sena_s(n,j,i-1) = sena2d(n,i,j)*wpm2
              tpr_s(n,j,i-1) = (prnca2d(i,j)+prca2d(i,j))*mmpd
              prcv_s(n,j,i-1) = prca2d(i,j)*mmpd
              ps_s(n,j,i-1) = p1d(n,i)*0.01

              q2m_o(j,i-1) = q2m_o(j,i-1) + q2m_1d(n,i)
              drag_o(j,i-1) = drag_o(j,i-1) + drag1d(n,i)
              evpa_o(j,i-1) = evpa_o(j,i-1) + evpa2d(n,i,j)
              sena_o(j,i-1) = sena_o(j,i-1) + sena2d(n,i,j)
            end do
            tpr_o(j,i-1) = (prnca2d(i,j)+prca2d(i,j))*mmpd
            q2m_o(j,i-1) = q2m_o(j,i-1)/float(ng)
            drag_o(j,i-1) = drag_o(j,i-1)/float(ng)
            evpa_o(j,i-1) = evpa_o(j,i-1)/float(ng)*mmpd
            sena_o(j,i-1) = sena_o(j,i-1)/float(ng)*wpm2
            flwa_o(j,i-1) = flwa2d(i,j)*wpm2
            fswa_o(j,i-1) = fswa2d(i,j)*wpm2
            flwd_o(j,i-1) = flwda2d(i,j)*wpm2
            sina_o(j,i-1) = sina2d(i,j)*wpm2
            prcv_o(j,i-1) = prca2d(i,j)*mmpd
            ps_o(j,i-1) = (sps2%ps(i,j)+r8pt)*10.
            zpbl_o(j,i-1) = sfsta%zpbl(i,j)

            tlef_o(j,i-1) = 0.0
            ssw_o(j,i-1) = 0.0
            rsw_o(j,i-1) = 0.0
            rnos_o(j,i-1) = 0.0
            scv_o(j,i-1) = 0.0
            nnn = 0
            do n = 1 , ng
              if ( ocld2d(n,i,j).ge.0.5 ) then
                tlef_o(j,i-1) = tlef_o(j,i-1) + tlef1d(n,i)
                ssw_o(j,i-1) = ssw_o(j,i-1) + ssw1d(n,i)
                rsw_o(j,i-1) = rsw_o(j,i-1) + rsw1d(n,i)
                rnos_o(j,i-1) = rnos_o(j,i-1) + rnos2d(n,i,j)
                scv_o(j,i-1) = scv_o(j,i-1) + scv1d(n,i)
                tlef_s(n,j,i-1) = tlef1d(n,i)
                ssw_s(n,j,i-1) = ssw1d(n,i)
                rsw_s(n,j,i-1) = rsw1d(n,i)
                rnos_s(n,j,i-1) = rnos2d(n,i,j)*mmpd
                scv_s(n,j,i-1) = scv1d(n,i)
                nnn = nnn + 1
              else
                tlef_s(n,j,i-1) = -1.D34
                ssw_s(n,j,i-1) = -1.D34
                rsw_s(n,j,i-1) = -1.D34
                rnos_s(n,j,i-1) = -1.D34
                scv_s(n,j,i-1) = -1.D34
              end if
            end do
            if ( nnn.ge.max0(ng/2,1) ) then
              tlef_o(j,i-1) = tlef_o(j,i-1)/float(nnn)
              ssw_o(j,i-1) = ssw_o(j,i-1)/float(nnn)
              rsw_o(j,i-1) = rsw_o(j,i-1)/float(nnn)
              rnos_o(j,i-1) = rnos_o(j,i-1)/float(nnn)*mmpd
              scv_o(j,i-1) = scv_o(j,i-1)/float(nnn)
            else
              tlef_o(j,i-1) = -1.D34
              ssw_o(j,i-1) = -1.D34
              rsw_o(j,i-1) = -1.D34
              rnos_o(j,i-1) = -1.D34
              scv_o(j,i-1) = -1.D34
            end if
#else
            drag_o(j-1,i-1) = 0.0
            q2m_o(j-1,i-1) = 0.0
            evpa_o(j-1,i-1) = 0.0
            sena_o(j-1,i-1) = 0.0
            do n = 1 , ng
              if ( ocld2d(n,i,j).ge.0.5 ) then
                fracv = sigf(n,i)
                fracb = (1.-veg1d(n,i))*(1.-scvk(n,i))
                fracs = veg1d(n,i)*wt(n,i) + (1.-veg1d(n,i))*scvk(n,i)
                facv = dlog(z1(n,i)/2.)/dlog(z1(n,i)/rough(lveg(n,i)))
                facb = dlog(z1(n,i)/2.)/dlog(z1(n,i)/zlnd)
                facs = dlog(z1(n,i)/2.)/dlog(z1(n,i)/zsno)
                fact = fracv*facv + fracb*facb + fracs*facs
                q2m_1d(n,i) = qs1d(n,i) - delq1d(n,i)*fact
              else if ( iocnflx.eq.1 ) then
                fact = dlog(z1(n,i)/2.)/dlog(z1(n,i)/zoce)
                q2m_1d(n,i) = qs1d(n,i) - delq1d(n,i)*fact
              else
              end if
              q2m_s(n,j-1,i-1) = q2m_1d(n,i)
              drag_s(n,j-1,i-1) = drag1d(n,i)
              evpa_s(n,j-1,i-1) = evpa2d(n,i,j)*mmpd
              sena_s(n,j-1,i-1) = sena2d(n,i,j)*wpm2
              tpr_s(n,j-1,i-1) = (prnca2d(i,j)+prca2d(i,j))*mmpd
              prcv_s(n,j-1,i-1) = prca2d(i,j)*mmpd
              ps_s(n,j-1,i-1) = p1d(n,i)*0.01

              q2m_o(j-1,i-1) = q2m_o(j-1,i-1) + q2m_1d(n,i)
              drag_o(j-1,i-1) = drag_o(j-1,i-1) + drag1d(n,i)
              evpa_o(j-1,i-1) = evpa_o(j-1,i-1) + evpa2d(n,i,j)
              sena_o(j-1,i-1) = sena_o(j-1,i-1) + sena2d(n,i,j)
            end do
            tpr_o(j-1,i-1) = (prnca2d(i,j)+prca2d(i,j))*mmpd
            q2m_o(j-1,i-1) = q2m_o(j-1,i-1)/float(ng)
            drag_o(j-1,i-1) = drag_o(j-1,i-1)/float(ng)
            evpa_o(j-1,i-1) = evpa_o(j-1,i-1)/float(ng)*mmpd
            sena_o(j-1,i-1) = sena_o(j-1,i-1)/float(ng)*wpm2
            flwa_o(j-1,i-1) = flwa2d(i,j)*wpm2
            fswa_o(j-1,i-1) = fswa2d(i,j)*wpm2
            flwd_o(j-1,i-1) = flwda2d(i,j)*wpm2
            sina_o(j-1,i-1) = sina2d(i,j)*wpm2
            prcv_o(j-1,i-1) = prca2d(i,j)*mmpd
            ps_o(j-1,i-1) = (sps2%ps(i,j)+r8pt)*10.
            zpbl_o(j-1,i-1) = sfsta%zpbl(i,j)

            tlef_o(j-1,i-1) = 0.0
            ssw_o(j-1,i-1) = 0.0
            rsw_o(j-1,i-1) = 0.0
            rnos_o(j-1,i-1) = 0.0
            scv_o(j-1,i-1) = 0.0
            nnn = 0
            do n = 1 , ng
              if ( ocld2d(n,i,j).ge.0.5 ) then
                tlef_o(j-1,i-1) = tlef_o(j-1,i-1) + tlef1d(n,i)
                ssw_o(j-1,i-1) = ssw_o(j-1,i-1) + ssw1d(n,i)
                rsw_o(j-1,i-1) = rsw_o(j-1,i-1) + rsw1d(n,i)
                rnos_o(j-1,i-1) = rnos_o(j-1,i-1) + rnos2d(n,i,j)
                if (abs(scv1d(n,i)) > 1D-34) then
                  scv_o(j-1,i-1) = scv_o(j-1,i-1) + scv1d(n,i)
                end if
                tlef_s(n,j-1,i-1) = tlef1d(n,i)
                ssw_s(n,j-1,i-1) = ssw1d(n,i)
                rsw_s(n,j-1,i-1) = rsw1d(n,i)
                rnos_s(n,j-1,i-1) = rnos2d(n,i,j)*mmpd
                if (abs(scv1d(n,i)) > 1D-34) then
                  scv_s(n,j-1,i-1) = scv1d(n,i)
                end if
                nnn = nnn + 1
              else
                tlef_s(n,j-1,i-1) = -1.D34
                ssw_s(n,j-1,i-1) = -1.D34
                rsw_s(n,j-1,i-1) = -1.D34
                rnos_s(n,j-1,i-1) = -1.D34
                scv_s(n,j-1,i-1) = -1.D34
              end if
            end do
            if ( nnn.ge.max0(ng/2,1) ) then
              tlef_o(j-1,i-1) = tlef_o(j-1,i-1)/float(nnn)
              ssw_o(j-1,i-1) = ssw_o(j-1,i-1)/float(nnn)
              rsw_o(j-1,i-1) = rsw_o(j-1,i-1)/float(nnn)
              rnos_o(j-1,i-1) = rnos_o(j-1,i-1)/float(nnn)*mmpd
              scv_o(j-1,i-1) = scv_o(j-1,i-1)/float(nnn)
            else
              tlef_o(j-1,i-1) = -1.D34
              ssw_o(j-1,i-1) = -1.D34
              rsw_o(j-1,i-1) = -1.D34
              rnos_o(j-1,i-1) = -1.D34
              scv_o(j-1,i-1) = -1.D34
            end if
#endif
#endif
 
!           ******    reset accumulation arrays to zero
            do n = 1 , ng
              evpa2d(n,i,j) = 0.
              rnos2d(n,i,j) = 0.
              sena2d(n,i,j) = 0.
            end do
            prnca2d(i,j) = 0.
            prca2d(i,j) = 0.
            flwa2d(i,j) = 0.
            flwda2d(i,j) = 0.
            fswa2d(i,j) = 0.
            svga2d(i,j) = 0.
            sina2d(i,j) = 0.
          end do
        end if
!
      else ! end ivers test
      end if
!
      call time_end(subroutine_name,idindx)
      end subroutine interf
!
      subroutine albedov(j,iemiss)
 
      implicit none
!
! Dummy arguments
!
      integer :: iemiss , j
      intent (in) iemiss , j
!
! Local variables
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
      if (idesseas == 1) then
        if (lmonth==1.or.lmonth==2.or.lmonth==12) then
          solour(1)=0.12
        endif        
        if (lmonth==3.or.lmonth==4.or.lmonth==5) then
          solour(1)=0.15
        endif        
        if (lmonth==6.or.lmonth==7.or.lmonth==8) then
          solour(1)=0.18
        endif        
        if (lmonth==9.or.lmonth==10.or.lmonth==11) then
          solour(1)=0.15
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
        czen(i) = max(coszrs(i),0.D0)
        czeta = czen(i)
        do n = 1 , nnsg
          albgs = 0.0D0
          albgl = 0.0D0
          albgsd = 0.0D0
          albgld = 0.0D0
          albs = 0.0D0
          albl = 0.0D0
          albsd = 0.0D0
          albld = 0.0D0
 
          albvs_s(n) = 0.0D0
          albvl_s(n) = 0.0D0
 
!================================================================
!         2.   get albedo over land
!================================================================
!         can't use pointer "nalbk" here because not set - use nldock
!         instead tgb1d(i) used instead of tbelow
!
          if (iseaice == 1) then
            if ( ldoc1d(n,i).gt.1.5 ) then
              tdiffs = ts1d(n,i) - tzero
              tdiff = max(tdiffs,0.D0)
              tdiffs = min(tdiff,20.D0)
              albgl = sical1 - 1.1D-2*tdiffs
              albgs = sical0 - 2.45D-2*tdiffs
              albg = fsol1*albgs + fsol2*albgl
              albgsd = albgs
              albgld = albgl
            end if
          else if ( ldoc1d(n,i).gt.0.1D0 .and. sice1d(n,i).eq.0.D0 ) then
            sfac = 1.D0 - fseas(tgb1d(n,i))
!           **********  ccm tests here on land mask for veg and soils
!c          data *** reduces albedo at low temps !!!!!should respond to
!c          moisture too the following card inactivated (commented out)
!c          (pat, 27 oct 86)
!           veg1d(i)=vegc(lveg(i))-seasf(lveg(i))*sfac
            albs = albvgs(lveg(n,i))
            albl = albvgl(lveg(n,i))
 
!----------------------------------------------------------------------
            if ( (lveg(n,i).lt.12) .or. (lveg(n,i).gt.15) ) then
 
!             2.1  bare soil albedos
!             (soil albedo depends on moisture)
              kolour = kolsol(lveg(n,i))
              wet = ssw1d(n,i)/depuv(lveg(n,i))
              alwet = max((11.D0-40.D0*wet),0.D0)*0.01D0
              alwet = min(alwet,solour(kolour))
              albg = solour(kolour) + alwet
!             if((lveg(n,i).eq.8)) albg=0.40      !Laura, cambiato il
!             DESERTO
              albgs = albg
              albgl = 2.D0*albg
!             **********            higher nir albedos
!             **********              set diffuse albedo
              albgld = albgl
              albgsd = albgs
              albsd = albs
              albld = albl
 
!             Dec. 15   albzn=0.85+1./(1.+10.*czen(i))
!             Dec. 12, 2008
              albzn = 1.0D0
!             Dec. 15, 2008
 
!             **********            leafless hardwood canopy: no or
!             inverse zen dep
              if ( lveg(n,i).eq.5 .and. sfac.lt.0.1 ) albzn = 1.
!             **********            multiply by zenith angle correction
              albs = albs*albzn
              albl = albl*albzn
 
!             **********            albedo over vegetation after zenith
!             angle corr
              albvs_s(n) = albs
              albvl_s(n) = albl
 
            else if ( lveg(n,i).eq.12 ) then
 
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
 
          else if ( sice1d(n,i).gt.0.D0 ) then
!====================================================================
!           3.  get albedo over sea ice
!====================================================================
!           **********          albedo depends on wave-length and ts.
!           the ts **********          dependence accounts for melt
!           water puddles.
            tdiffs = ts1d(n,i) - tzero
            tdiff = max(tdiffs,0.D0)
            tdiffs = min(tdiff,20.D0)
            albgl = sical1 - 1.1D-2*tdiffs
            albgs = sical0 - 2.45D-2*tdiffs
            albg = fsol1*albgs + fsol2*albgl
            albgsd = albgs
            albgld = albgl
          end if
! ===================================================================
!         4.  correct for snow cover
! ===================================================================
          if ( scv1d(n,i).gt.0.0D0 ) then
!           **********            snow albedo depends on  snow-age,
!           zenith angle, **********            and thickness of snow
 
!           **********            zenith angle set in zenitm
!           **********            snow albedoes for visible and ir
!           solar rad **********            visible albedo depends on
!           snow age **********            age gives reduction of
!           visible rad snow albedo **********              due to age
            cons = 0.2D0
            conn = 0.5D0
            age = (1.D0-1.D0/(1.D0+sag1d(n,i)))
!           **********            sl helps control albedo zenith
!           dependence
            sl = 2.0D0
            sli = 1.D0/sl
            sl2 = 2.D0*sl
!           **********            snal0= new snow albedo for vis rad,
!           sol zen le 6 **********            snal1= new snow albedo
!           for long-wave rad
            dfalbs = snal0*(1.D0-cons*age)
!           **********            czf corrects albedo of new snow for
!           solar zenith
            cf1 = ((1.D0+sli)/(1.D0+sl2*czen(i))-sli)
            cff = max(cf1,0.D0)
            czf = 0.4D0*cff*(1.D0-dfalbs)
            dralbs = dfalbs + czf
            dfalbl = snal1*(1.D0-conn*age)
            czf = 0.4D0*cff*(1.D0-dfalbl)
            dralbl = dfalbl + czf
 
            if ( veg1d(n,i).gt.0.001D0 ) then
!             **********            effective albedo over vegetation
!             with snow
              albl = (1.D0-wt(n,i))*albl + dralbl*wt(n,i)
              albld = (1.D0-wt(n,i))*albld + dfalbl*wt(n,i)
              albs = (1.D0-wt(n,i))*albs + dralbs*wt(n,i)
              albsd = (1.D0-wt(n,i))*albsd + dfalbs*wt(n,i)
            end if
 
!----------------------------------------------------------------------
!           4.1  compute albedo for snow on bare ground
!----------------------------------------------------------------------
            albgs = (1.D0-scvk(n,i))*albgs + dralbs*scvk(n,i)
            albgl = (1.D0-scvk(n,i))*albgl + dralbl*scvk(n,i)
            albgsd = (1.D0-scvk(n,i))*albgsd + dfalbs*scvk(n,i)
            albgld = (1.D0-scvk(n,i))*albgld + dfalbl*scvk(n,i)
          end if
 
!=====================================================================
!         5.  albedo over open ocean
!=====================================================================
          if ( ldoc1d(n,i).eq.0.D0 ) then
!           *********   ocean albedo depends on zenith angle
            if ( czeta.ge.0.0D0 ) then
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
          aldirs_s(n) = (1.D0-veg1d(n,i))*albgs + veg1d(n,i)*albs
          aldirl_s(n) = (1.D0-veg1d(n,i))*albgl + veg1d(n,i)*albl
          aldifs_s(n) = (1.D0-veg1d(n,i))*albgsd + veg1d(n,i)*albsd
          aldifl_s(n) = (1.D0-veg1d(n,i))*albgld + veg1d(n,i)*albld
        end do
        albvs(i) = albvs_s(1)
        albvl(i) = albvl_s(1)
        aldirs(i) = aldirs_s(1)
        aldirl(i) = aldirl_s(1)
        aldifs(i) = aldifs_s(1)
        aldifl(i) = aldifl_s(1)
        if ( iemiss.eq.1 ) emiss1d(i) = emiss2d(1,i,j)
        aldirs1d(1,i) = aldirs_s(1)
        aldifs1d(1,i) = aldifs_s(1)
        do n = 2 , nnsg
          albvs(i) = albvs(i) + albvs_s(n)
          albvl(i) = albvl(i) + albvl_s(n)
          aldirs(i) = aldirs(i) + aldirs_s(n)
          aldirl(i) = aldirl(i) + aldirl_s(n)
          aldifs(i) = aldifs(i) + aldifs_s(n)
          aldifl(i) = aldifl(i) + aldifl_s(n)
          if ( iemiss.eq.1 ) emiss1d(i) = emiss1d(i) + emiss2d(n,i,j)
          aldirs1d(n,i) = aldirs_s(n)
          aldifs1d(n,i) = aldifs_s(n)
        end do
        albvs(i) = albvs(i)/dble(nnsg)
        albvl(i) = albvl(i)/dble(nnsg)
        aldirs(i) = aldirs(i)/dble(nnsg)
        aldirl(i) = aldirl(i)/dble(nnsg)
        aldifs(i) = aldifs(i)/dble(nnsg)
        aldifl(i) = aldifl(i)/dble(nnsg)
        if ( iemiss.eq.1 ) emiss1d(i) = emiss1d(i)/dble(nnsg)
 
!       ******   fsw1d(i),sabveg(i),solis(i) computed in colrad
 
      end do
 
      call time_end(subroutine_name,idindx)

      contains

        function fseas(x)
          implicit none
          real(8) :: fseas
          real(8) , intent(in) :: x
          fseas = max(0.D0,(1.D0-0.0016D0*max(298.D0-x,0.D0)**2.0D0))
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
! Local variables
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
 
          if ( lveg(n,i).ne.0 ) then
 
!           **********            lveg is set in subr. interf
            freza(lveg(n,i)) = 0.15*deprv(lveg(n,i))
            frezu(lveg(n,i)) = 0.15*depuv(lveg(n,i))
            itex = iexsol(lveg(n,i))
            texrat(n,i) = skrat(itex)
            porsl(n,i) = xmopor(itex)
            xkmx(n,i) = xmohyd(itex)
            bsw(n,i) = bee(itex)
            bfc(n,i) = 5.8 - bsw(n,i)*(0.8+0.12*(bsw(n,i)-4.)*          &
                      & dlog10(1.D2*xkmx(n,i)))
 
            phi0 = xmosuc(itex)
            dmax = bsw(n,i)*phi0*xkmx(n,i)/porsl(n,i)
            dmin = 1.D-3
            dmnor = 1550.*dmin/dmax
            tweak1 = (bsw(n,i)*(bsw(n,i)-6.)+10.3)                      &
                   & /(bsw(n,i)*bsw(n,i)+40.*bsw(n,i))
            ck = (1.+dmnor)*tweak1*0.23/0.02356
            evmx0(n,i) = 1.02*dmax*ck/dsqrt(depuv(lveg(n,i))*           &
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
! Dummy arguments
!
      integer , intent (in) :: imax
      real(kind=8) , intent (in) :: alat , dec , fjd
      real(kind=8) , intent (out) , dimension(imax) :: coszrs , frac
!
! Local variables
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
        frac(i) = 1.0
        ha = (i-1)*dlon + tpifjd
!       if cosz is negative, the sun is below the horizon.
        cosz = max(0.D0,ss+cc*dcos(ha))
        coszrs(i) = cosz
      end do
!
      call time_end(subroutine_name,idindx)
      end subroutine zenith
!
      subroutine slice1D(j)
 
      implicit none
!
! Dummy arguments
!
      integer :: j
      intent (in) j
!
! Local variables
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
          lveg(n,i) = nint(veg2d1(n,i,j))
          if (iseaice == 1) then
            if (ocld2d(n,i,j) > 1.5) lveg(n,i) = 12
          end if
          amxtem = max(298.-tgb1d(n,i),0.D0)
          sfac = 1. - max(0.D0,1.-0.0016*amxtem**2)
          if ( lveg(n,i).eq.0 ) then
            veg1d(n,i) = 0.
          else
            veg1d(n,i) = vegc(lveg(n,i)) - seasf(lveg(n,i))*sfac
          end if
          ts1d(n,i) = thx3d(i,kz,j)-6.5D-3*rgti*(ht1(n,i,j)- &
                                                 mddom%ht(i,j))
          scv1d(n,i) = scv2d(n,i,j)
          sag1d(n,i) = sag2d(n,i,j)
        end do
      end do
 
      end subroutine slice1D
!
      end module mod_vecbats
