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
 
module mod_bats_mtrxbats

  use mod_realkinds
  use mod_dynparam
  use mod_mpmessage
  use mod_constants
  use mod_service
  use mod_bats_common
  use mod_bats_lake
  use mod_bats_bndry
  use mod_bats_drag
  use mod_bats_mppio
  use mod_bats_zengocn
  use mod_bats_romsocn

  private

  public :: interf , initb , mtrxbats , albedov , slice1D

  contains

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
!      matrix bats (2011) was developed from vecbats by Graziano Giuliani
! it is designed such that one call to bats perfoms calculations at
! a specified arbitrary number of gridpoints in two dimensions.
! this allows a more efficient mpi parallelization of the code.
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
!   mtrxbats ==> soilbc
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
! 
  subroutine mtrxbats(jstart,jend,istart,iend,ktau)
!
    implicit none
!
    integer , intent(in) :: jstart , jend , istart , iend
    integer(8) , intent(in) :: ktau
    character (len=64) :: subroutine_name='mtrxbats'
    integer :: idindx=0
!
    call time_begin(subroutine_name,idindx)
!
!---------------------------------------------------------------------
!
!   Excange from model to BATS

    call interf(1,jstart,jend,istart,iend,ktau)

!   Calculate surface fluxes and hydrology budgets

    call soilbc(jstart,jend,istart,iend)

!   Albedo is calculated in mod_tendency

    call bndry(jstart,jend,istart,iend)

!   Zeng ocean flux model
    if ( iocnflx == 2 ) call zengocndrv(jstart,jend,istart,iend,ktau)

!   ROMS ocean model
    if ( iocnflx == 3 ) then
!     Call Zeng ocean flux model until ROMS runs once
      if (ktau <= ntcpl) then
        print*, ktau, ntcpl, "ocean - zeng"
        call zengocndrv(jstart,jend,istart,iend,ktau)
      else
!       Then update ground temperature in each coupling time step
        if (mod(ktau+1,ntcpl) == ntsrf2) then 
          print*, ktau, ntcpl, "ocean - roms"
          call romsocndrv(jstart,jend,istart,iend,ktau)
        end if
      end if
    end if

!   Hostetler lake model for every BATS timestep at lake points
    if ( llake ) then
      call lakedrv(jstart,jend,istart,iend)
    endif

!   Accumulate quantities for energy and moisture budgets

    call interf(2,jstart,jend,istart,iend,ktau)
! 
    call time_end(subroutine_name,idindx)
!
  end subroutine mtrxbats
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!     ***  provides constant initial fields to boundary subroutine
!     ***    soil textures are set in soilbc
!     ***    soil   colors are set in albedo
!
!     ***   units are si
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
  subroutine initb(jstart,jend,istart,iend)
    implicit none
    integer , intent(in) :: jstart , jend , istart , iend
! 
    integer :: i , is , itex , j , n , nlveg
!
    character (len=64) :: subroutine_name='initb'
    integer :: idindx=0
!
    call time_begin(subroutine_name,idindx)
!
    do i = istart , iend
      do j = jstart , jend
        pptnc(i,j) = d_zero
        pptc(i,j) = d_zero
        veg2d(i,j) = idnint(lndcat(i,j)+0.1D0)
        do n = 1 , nnsg
          veg2d1(n,i,j) = idnint(lndcat1(n,i,j)+0.1D0)
        end do
      end do
    end do

!   initialize hostetler lake model
    if ( llake ) call initlake(jstart,jend,istart,iend)
 
    do i = istart , iend
      do j = jstart , jend
        do n = 1 , nnsg
          tg2d(n,i,j) = tground2(i,j)
          tgb2d(n,i,j) = tground2(i,j)
          taf2d(n,i,j) = tground2(i,j)
          tlef2d(n,i,j) = tground2(i,j)

          if ( ocld2d(n,i,j) == 2 ) then
            if ( lseaice .or. llake ) then
              sice2d(n,i,j) = d_1000
              scv2d(n,i,j) = d_zero
            end if
            nlveg = 12
          else if ( ocld2d(n,i,j) == 1 ) then
            sice2d(n,i,j) = d_zero
            scv2d(n,i,j) = dmax1(scv2d(n,i,j),d_zero)
            nlveg = veg2d1(n,i,j)
          else
            sice2d(n,i,j) = d_zero
            scv2d(n,i,j) = d_zero
            nlveg = veg2d1(n,i,j)
          end if
!         Initialize soil moisture in the 3 layers
          is = idint(lndcat1(n,i,j))
          itex = iexsol(nlveg)
          swt2d(n,i,j) = deptv(nlveg)*xmopor(itex)*slmo(is)
          srw2d(n,i,j) = deprv(nlveg)*xmopor(itex)*slmo(is)
          ssw2d(n,i,j) = depuv(nlveg)*xmopor(itex)*slmo(is)

          dew2d(n,i,j) = d_zero
          sag2d(n,i,j) = d_zero
          gwet2d(n,i,j) = d_half
          sena2d(n,i,j) = d_zero
          evpa2d(n,i,j) = d_zero
          rnos2d(n,i,j) = d_zero
          rno2d(n,i,j) = d_zero
          ircp2d(n,i,j) = d_zero
        end do
      end do
    end do
 
    do i = istart , iend
      do j = jstart , jend
        fsw2d(i,j) = d_zero
        flw2d(i,j) = d_zero
        sabv2d(i,j) = d_zero
        sol2d(i,j) = d_zero
        fswa2d(i,j) = d_zero
        flwa2d(i,j) = d_zero
        prca2d(i,j) = d_zero
        prnca2d(i,j) = d_zero
        svga2d(i,j) = d_zero
        sina2d(i,j) = d_zero
      end do
    end do
! 
    call time_end(subroutine_name,idindx)
!
  end subroutine initb
!
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!  this subroutine interfaces mm42d and bats variables
!
!  ivers = 1 ,   regcm2d --> bats
!  ivers = 2 ,   bats --> regcm2d
!
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
  subroutine interf(ivers,jstart,jend,istart,iend,ktau)
!
    implicit none
!
    integer , intent (in) :: ivers , jstart , jend , istart , iend
    integer(8) , intent(in) :: ktau
!
    real(dp) :: amxtem , facb , facs , fact , factuv , facv , fracb ,  &
               fracs , fracv , hl , mmpd , rh0 , satvp , sfac ,       &
               solvt , wpm2 , p0 , qs0 , ts0
    integer :: i , j , n , nnn
    real(sp) :: real_4
    logical , save :: first_pass
!
    character (len=64) :: subroutine_name='interf'
    integer :: idindx=0
!
    data first_pass /.true./
!
    call time_begin(subroutine_name,idindx)
 
    if ( ivers == 1 ) then ! regcm2d --> bats

      do i = istart, iend
        do j = jstart, jend
          p0 = (sfps(i,j)+ptop)*d_1000
          qs0 = qvatm(i,kz,j)/(d_one+qvatm(i,kz,j))
          ts0 = thatm(i,kz,j)
          hl = lh0 - lh1*(ts0-tzero)
          satvp = lsvp1*dexp(lsvp2*hl*(d_one/tzero-d_one/ts0))
          rh0 = dmax1(qs0/(ep2*satvp/(p0*0.01D0-satvp)),d_zero)
          do n = 1 , nnsg
            qs(n,j,i) = qs0
            sts(n,j,i) = ts0-lrate*regrav*(ht1(n,i,j)-ht(i,j))
            sfcp(n,j,i) = p0*(sts(n,j,i)/ts0)
            hl = lh0 - lh1*(sts(n,j,i)-tzero)
            satvp = lsvp1*dexp(lsvp2*hl*(d_one/tzero-d_one/sts(n,j,i)))
            qs(n,j,i) = dmax1(rh0*ep2*satvp/(sfcp(n,j,i)*0.01D0-satvp),d_zero)
            tgrd(n,j,i) = tg2d(n,i,j)
            rhs(n,j,i) = sfcp(n,j,i)/(rgas*sts(n,j,i))
            prcp(n,j,i) = pptnc(i,j) + pptc(i,j)
            !
            ! quantities stored on 2d surface array for bats use only
            !
            tgbrd(n,j,i) = tgb2d(n,i,j)
            taf(n,j,i) = taf2d(n,i,j)
            tlef(n,j,i) = tlef2d(n,i,j)
            tsw(n,j,i) = swt2d(n,i,j)
            rsw(n,j,i) = srw2d(n,i,j)
            ssw(n,j,i) = ssw2d(n,i,j)
            ldew(n,j,i) = dew2d(n,i,j)
            snag(n,j,i) = sag2d(n,i,j)
            sncv(n,j,i) = scv2d(n,i,j)
            sfice(n,j,i) = sice2d(n,i,j)
            gwet(n,j,i) = gwet2d(n,i,j)
            sent(n,j,i) = hfx(i,j)
            evpr(n,j,i) = qfx(i,j)
            ldimsk(n,j,i) = ocld2d(n,i,j)
            ircp(n,j,i) = ircp2d(n,i,j)
            lveg(n,j,i) = veg2d1(n,i,j)
            oveg(n,j,i) = lveg(n,j,i)
            if ( ldimsk(n,j,i) == 2 ) lveg(n,j,i) = 12
            amxtem = dmax1(298.0D0-tgbrd(n,j,i),d_zero)
            sfac = d_one - dmax1(d_zero,d_one-0.0016D0*amxtem**d_two)
            lncl(n,j,i) = mfcv(lveg(n,j,i)) - seasf(lveg(n,j,i))*sfac
            emiss(n,j,i) = emiss2d(n,i,j)
            zh(n,j,i) = hgt(i,kz,j)
            z1log(n,j,i)  = dlog(zh(n,j,i))
            z2fra(n,j,i)  = dlog(zh(n,j,i)*d_half)
            z10fra(n,j,i) = dlog(zh(n,j,i)*d_r10)
            zlgocn(n,j,i) = dlog(zh(n,j,i)/zoce)
            zlglnd(n,j,i) = dlog(zh(n,j,i)/zlnd)
            zlgsno(n,j,i) = dlog(zh(n,j,i)/zsno)
            zlgveg(n,j,i) = dlog(zh(n,j,i)/rough(lveg(n,j,i)))
            zlgdis(n,j,i) = dlog(zh(n,j,i)-displa(lveg(n,j,i)) / &
                                           rough(lveg(n,j,i)))
          end do
 
          rh0 = d_zero
          do n = 1 , nnsg
            rh0 = rh0 + (qs(n,j,i)-qs0)
          end do
          rh0 = rh0/nnsg
          do n = 1 , nnsg
            qs(n,j,i) = dmax1(qs(n,j,i)-rh0,d_zero)
          end do
 
          usw(j,i) = uatm(i,kz,j)
          vsw(j,i) = vatm(i,kz,j)
          fsw(j,i) = fsw2d(i,j)
          flw(j,i) = flw2d(i,j)
          solis(j,i) = sol2d(i,j)
          sabveg(j,i) = sabv2d(i,j)
          solvt = solvd2d(i,j) + solvs2d(i,j)
          if ( solvt > d_zero ) then
            fracd(j,i) = solvd2d(i,j)/solvt
          else
            fracd(j,i) = 0.2D0
          end if
        end do
      end do
 
    else if ( ivers == 2 ) then ! bats --> regcm2d
 
      do i = istart, iend
        do j = jstart, jend
          uvdrag(i,j) = d_zero
          hfx(i,j) = d_zero
          qfx(i,j) = d_zero
          tground2(i,j) = d_zero
          tground1(i,j) = d_zero
          tgbb(i,j) = d_zero
          if ( lchem ) then
            ssw2da(i,j) = d_zero
            sdeltk2d(i,j) = d_zero
            sdelqk2d(i,j) = d_zero
            sfracv2d(i,j) = d_zero
            sfracb2d(i,j) = d_zero
            sfracs2d(i,j) = d_zero
            svegfrac2d(i,j) = d_zero
          end if

          do n = 1 , nnsg
            uvdrag(i,j) = uvdrag(i,j) + drag(n,j,i)
            hfx(i,j) = hfx(i,j) + sent(n,j,i)
            qfx(i,j) = qfx(i,j) + evpr(n,j,i)
            tground2(i,j) = tground2(i,j) + tgrd(n,j,i)
            tground1(i,j) = tground1(i,j) + tgrd(n,j,i)
            if ( lchem ) then
              ssw2da(i,j) = ssw2da(i,j) + ssw(n,j,i)
              sdeltk2d(i,j) = sdeltk2d(i,j) + delt(n,j,i)
              sdelqk2d(i,j) = sdelqk2d(i,j) + delq(n,j,i)
              sfracv2d(i,j) = sfracv2d(i,j) + sigf(n,j,i)
              sfracb2d(i,j) = sfracb2d(i,j)+ &
                              (d_one-lncl(n,j,i))*(d_one-scvk(n,j,i))
              sfracs2d(i,j) = sfracs2d(i,j) + &
                       lncl(n,j,i)*wt(n,j,i) + (d_one-lncl(n,j,i))*scvk(n,j,i)
              svegfrac2d(i,j) = svegfrac2d(i,j) + lncl(n,j,i)
            end if
            if ( iocnflx == 1 .or. &
                (iocnflx == 2 .and. ldimsk(n,j,i) /= 0 ) ) then
              tgbb(i,j) = tgbb(i,j) +                 &
                         ((d_one-lncl(n,j,i))*tgrd(n,j,i)**d_four +  &
                         lncl(n,j,i)*tlef(n,j,i)**d_four)**d_rfour
            else
              tgbb(i,j) = tgbb(i,j) + tgrd(n,j,i)
            end if
            if ( ldimsk(n,j,i) == 0 ) then
              ssw(n,j,i)  = dmissval
              rsw(n,j,i)  = dmissval
              tsw(n,j,i)  = dmissval
              trnof(n,j,i)  = dmissval
              srnof(n,j,i) = dmissval
              sncv(n,j,i)  = dmissval
            end if
          end do
          uvdrag(i,j) = uvdrag(i,j)*rdnnsg
          hfx(i,j) = hfx(i,j)*rdnnsg
          qfx(i,j) = qfx(i,j)*rdnnsg
          tground2(i,j) = tground2(i,j)*rdnnsg
          tground1(i,j) = tground1(i,j)*rdnnsg
          tgbb(i,j) = tgbb(i,j)*rdnnsg

          if ( lchem ) then
            ssw2da(i,j) = ssw2da(i,j)*rdnnsg
            sdeltk2d(i,j) = sdeltk2d(i,j)*rdnnsg
            sdelqk2d(i,j) = sdelqk2d(i,j)*rdnnsg
            sfracv2d(i,j) = sfracv2d(i,j)*rdnnsg
            sfracb2d(i,j) = sfracb2d(i,j)*rdnnsg
            sfracs2d(i,j) = sfracs2d(i,j)*rdnnsg
            svegfrac2d(i,j) = svegfrac2d(i,j)*rdnnsg
          end if
          do n = 1 , nnsg
            scv2d(n,i,j) = sncv(n,j,i)
            tg2d(n,i,j) = tgrd(n,j,i)
            tgb2d(n,i,j) = tgbrd(n,j,i)
            taf2d(n,i,j) = taf(n,j,i)
            tlef2d(n,i,j) = tlef(n,j,i)
            swt2d(n,i,j) = tsw(n,j,i)
            srw2d(n,i,j) = rsw(n,j,i)
            ssw2d(n,i,j) = ssw(n,j,i)
            dew2d(n,i,j) = ldew(n,j,i)
            sag2d(n,i,j) = snag(n,j,i)
            sice2d(n,i,j) = sfice(n,j,i)
            gwet2d(n,i,j) = gwet(n,j,i)
            ocld2d(n,i,j) = ldimsk(n,j,i)
            ircp2d(n,i,j) = ircp(n,j,i)
            evpa2d(n,i,j) = evpa2d(n,i,j) + dtbat*evpr(n,j,i)
            sena2d(n,i,j) = sena2d(n,i,j) + dtbat*sent(n,j,i)
            if ( dabs(srnof(n,j,i)) > 1.0D-10 ) then
              rnos2d(n,i,j) = rnos2d(n,i,j) + srnof(n,j,i)*dtbat
            end if
            if ( dabs(srnof(n,j,i))  > 1.0D-10 .and. &
                 dabs(trnof(n,j,i))   > 1.0D-10 ) then
              rno2d(n,i,j) = rno2d(n,i,j) + &
                     (trnof(n,j,i)-srnof(n,j,i))*dtbat
            end if
          end do
          !
          ! quantities stored on 2d surface array for bats use only
          !
          prca2d(i,j) = prca2d(i,j) + dtbat*pptc(i,j)
          prnca2d(i,j) = prnca2d(i,j) + dtbat*pptnc(i,j)
          flwa2d(i,j) = flwa2d(i,j) + dtbat*flw(j,i)
          flwda2d(i,j) = flwda2d(i,j) + dtbat*flwd2d(i,j)
          fswa2d(i,j) = fswa2d(i,j) + dtbat*fsw(j,i)
          svga2d(i,j) = svga2d(i,j) + dtbat*sabveg(j,i)
          sina2d(i,j) = sina2d(i,j) + dtbat*sinc2d(i,j)
          pptnc(i,j) = d_zero
          pptc(i,j) = d_zero
        end do
      end do

      do i = istart, iend
        do j = jstart, jend
          u10m_o(j,i-1) = 0.0
          v10m_o(j,i-1) = 0.0
          tg_o(j,i-1) = 0.0
          t2m_o(j,i-1) = 0.0
          q2m_o(j,i-1) = 0.0
          aldirs_o(j,i-1) = 0.0
          aldifs_o(j,i-1) = 0.0
          do n = 1 , nnsg
            if ( ldimsk(n,j,i) /= 0 ) then
              fracv = sigf(n,j,i)
              fracb = (d_one-lncl(n,j,i))*(d_one-scvk(n,j,i))
              fracs = lncl(n,j,i)*wt(n,j,i) + (d_one-lncl(n,j,i))*scvk(n,j,i)
              facv = z2fra(n,j,i)/zlgveg(n,j,i)
              facb = z2fra(n,j,i)/zlglnd(n,j,i)
              facs = z2fra(n,j,i)/zlgsno(n,j,i)
              fact = fracv*facv + fracb*facb + fracs*facs
              facv = z10fra(n,j,i)/zlgveg(n,j,i)
              facb = z10fra(n,j,i)/zlglnd(n,j,i)
              facs = z10fra(n,j,i)/zlgsno(n,j,i)
              factuv = fracv*facv + fracb*facb + fracs*facs
              u10m(n,j,i) = usw(j,i)*(d_one-factuv)
              v10m(n,j,i) = vsw(j,i)*(d_one-factuv)
              t2m(n,j,i) = sts(n,j,i) - delt(n,j,i)*fact
              q2m(n,j,i) = qs(n,j,i) - delq(n,j,i)*fact
            else 
              if ( iocnflx == 1 ) then
                fact = z2fra(n,j,i)/zlgocn(n,j,i)
                factuv = z10fra(n,j,i)/zlgocn(n,j,i)
                u10m(n,j,i) = usw(j,i)*(d_one-factuv)
                v10m(n,j,i) = vsw(j,i)*(d_one-factuv)
                t2m(n,j,i) = sts(n,j,i) - delt(n,j,i)*fact
                q2m(n,j,i) = qs(n,j,i) - delq(n,j,i)*fact
              end if
            end if
            tg_s(n,j,i-1) = real(tgrd(n,j,i))
            u10m_s(n,j,i-1) = real(u10m(n,j,i))
            v10m_s(n,j,i-1) = real(v10m(n,j,i))
            t2m_s(n,j,i-1) = real(t2m(n,j,i))
            q2m_s(n,j,i-1) = real(q2m(n,j,i))
 
            u10m_o(j,i-1) = u10m_o(j,i-1) + real(u10m(n,j,i))
            v10m_o(j,i-1) = v10m_o(j,i-1) + real(v10m(n,j,i))
            t2m_o(j,i-1) = t2m_o(j,i-1) + real(t2m(n,j,i))
            q2m_o(j,i-1) = real(q2m_o(j,i-1) + q2m(n,j,i))
            tg_o(j,i-1) = tg_o(j,i-1) + real(tgrd(n,j,i))
            aldirs_o(j,i-1) = aldirs_o(j,i-1) + real(albdirs(n,j,i))
            aldifs_o(j,i-1) = aldifs_o(j,i-1) + real(albdifs(n,j,i))
          end do
          u10m_o(j,i-1) = u10m_o(j,i-1)*rrnnsg
          v10m_o(j,i-1) = v10m_o(j,i-1)*rrnnsg
          t2m_o(j,i-1) = t2m_o(j,i-1)*rrnnsg
          q2m_o(j,i-1) = q2m_o(j,i-1)*rrnnsg
          tg_o(j,i-1) = tg_o(j,i-1)*rrnnsg
          aldirs_o(j,i-1) = aldirs_o(j,i-1)*rrnnsg
          aldifs_o(j,i-1) = aldifs_o(j,i-1)*rrnnsg
 
          tgmx_o(j,i-1) = amax1(tgmx_o(j,i-1),tg_o(j,i-1))
          tgmn_o(j,i-1) = amin1(tgmn_o(j,i-1),tg_o(j,i-1))
          t2mx_o(j,i-1) = amax1(t2mx_o(j,i-1),t2m_o(j,i-1))
          t2mn_o(j,i-1) = amin1(t2mn_o(j,i-1),t2m_o(j,i-1))
          w10x_o(j,i-1) = amax1(w10x_o(j,i-1), &
                          sqrt(u10m_o(j,i-1)**2.0+v10m_o(j,i-1)**2.0))
          real_4 = real((sfps(i,j)+ptop)*d_10)
          psmn_o(j,i-1) = amin1(psmn_o(j,i-1),real_4)
        end do
      end do

      if ( mod(ktau+1,kbats) == 0 .or. first_pass ) then
        if ( ktau <= 1 ) then
          mmpd = secpd/dtbat
          wpm2 = d_one/dtbat
        else if ( ktau+1 == kbats ) then
          mmpd = houpd/(srffrq-xdtsec/secph)
          wpm2 = d_one/((srffrq-xdtsec/secph)*secph)
        else
          mmpd = houpd/srffrq
          wpm2 = d_one/(srffrq*secph)
        end if
        do i = istart, iend
          do j = jstart, jend
            drag_o(j,i-1) = 0.0
            evpa_o(j,i-1) = 0.0
            sena_o(j,i-1) = 0.0
            do n = 1 , nnsg
              if ( ldimsk(n,j,i) /= 0 ) then
                fracv = sigf(n,j,i)
                fracb = (d_one-lncl(n,j,i))*(d_one-scvk(n,j,i))
                fracs = lncl(n,j,i)*wt(n,j,i) + (d_one-lncl(n,j,i))*scvk(n,j,i)
                facv = z2fra(n,j,i)/zlgveg(n,j,i)
                facb = z2fra(n,j,i)/zlglnd(n,j,i)
                facs = z2fra(n,j,i)/zlgsno(n,j,i)
                fact = fracv*facv + fracb*facb + fracs*facs
              else
                if ( iocnflx == 1 ) then
                  fact = z2fra(n,j,i)/zlgocn(n,j,i)
                end if
              end if
              drag_s(n,j,i-1) = real(drag(n,j,i))
              evpa_s(n,j,i-1) = real(evpa2d(n,i,j)*mmpd)
              sena_s(n,j,i-1) = real(sena2d(n,i,j)*wpm2)
              tpr_s(n,j,i-1) = real((prnca2d(i,j)+prca2d(i,j))*mmpd)
              prcv_s(n,j,i-1) = real(prca2d(i,j)*mmpd)
              ps_s(n,j,i-1) = real(sfcp(n,j,i)*0.01D0)
 
              drag_o(j,i-1) = real(drag_o(j,i-1) + drag(n,j,i))
              evpa_o(j,i-1) = real(evpa_o(j,i-1) + evpa2d(n,i,j))
              sena_o(j,i-1) = real(sena_o(j,i-1) + sena2d(n,i,j))
            end do
            tpr_o(j,i-1) = real((prnca2d(i,j)+prca2d(i,j))*mmpd)
            drag_o(j,i-1) = drag_o(j,i-1)*rrnnsg
            evpa_o(j,i-1) = evpa_o(j,i-1)*rrnnsg*real(mmpd)
            sena_o(j,i-1) = sena_o(j,i-1)*rrnnsg*real(wpm2)
            flwa_o(j,i-1) = real(flwa2d(i,j)*wpm2)
            fswa_o(j,i-1) = real(fswa2d(i,j)*wpm2)
            flwd_o(j,i-1) = real(flwda2d(i,j)*wpm2)
            sina_o(j,i-1) = real(sina2d(i,j)*wpm2)
            prcv_o(j,i-1) = real(prca2d(i,j)*mmpd)
            ps_o(j,i-1) = real((sfps(i,j)+ptop)*d_10)
            zpbl_o(j,i-1) = real(hpbl(j,i))
 
            tlef_o(j,i-1) = 0.0
            ssw_o(j,i-1) = 0.0
            rsw_o(j,i-1) = 0.0
            rnos_o(j,i-1) = 0.0
            scv_o(j,i-1) = 0.0
            nnn = 0
            do n = 1 , nnsg
              if ( ldimsk(n,j,i) /= 0 ) then
                tlef_o(j,i-1) = tlef_o(j,i-1) + real(tlef(n,j,i))
                ssw_o(j,i-1) = ssw_o(j,i-1) + real(ssw(n,j,i))
                rsw_o(j,i-1) = rsw_o(j,i-1) + real(rsw(n,j,i))
                rnos_o(j,i-1) = rnos_o(j,i-1) + real(rnos2d(n,i,j))
                scv_o(j,i-1) = scv_o(j,i-1) + real(sncv(n,j,i))
                tlef_s(n,j,i-1) = real(tlef(n,j,i))
                ssw_s(n,j,i-1) = real(ssw(n,j,i))
                rsw_s(n,j,i-1) = real(rsw(n,j,i))
                rnos_s(n,j,i-1) = real(rnos2d(n,i,j)*mmpd)
                scv_s(n,j,i-1) = real(sncv(n,j,i))
                nnn = nnn + 1
              else
                tlef_s(n,j,i-1) = smissval
                ssw_s(n,j,i-1) = smissval
                rsw_s(n,j,i-1) = smissval
                rnos_s(n,j,i-1) = smissval
                scv_s(n,j,i-1) = smissval
              end if
            end do
            if ( nnn >= max0(nnsg/2,1) ) then
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
            !
            ! reset accumulation arrays to zero
            !
            do n = 1 , nnsg
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
        end do
      end if
    end if
    first_pass = .false.
    call time_end(subroutine_name,idindx)
  end subroutine interf
!
! Albedo calculates fragmented albedos (direct and diffuse) in
! wavelength regions split at 0.7um.
!
! CM hands albedos to radiation package which computes
! fsw(i) = net solar absorbed over full grid square
! sabveg(j,i) = vegetation absorbed (full solar spectrum)
! solis(j,i) = shortwave  solar incident
!
! Here these are calculated at the end of albedo - they use only
! direct albedos for now
!
! in both versions :  lftemp uses sabveg
! tgrund uses sabveg & fsw(i) to get
! ground absorbed solar
! photosynthesis uses solis - see subrouts
! stomat and co2 (carbon)
!
! For sea, sea-ice veg albedos are not set these albedos are not
! treated as arrays here
!
! (depuv/10.0)= the ratio of upper soil layer to total root depth
! Used to compute "wet" for soil albedo
!
  subroutine albedov(imon,jstart,jend,istart,iend)
    implicit none
    integer , intent(in) :: jstart , jend , istart , iend
    integer , intent (in) :: imon
!
    real(dp) :: age , albg , albgl , albgld , albgs , albgsd , albl ,  &
               albld , albs , albsd , albzn , alwet , cf1 , cff ,     &
               conn , cons , czeta , czf , dfalbl , dfalbs , dralbl , &
               dralbs , fsol1 , fsol2 , sfac , sical0 , sical1 , sl , &
               sl2 , sli , snal0 , snal1 , tdiff , tdiffs , wet
    real(dp) , dimension(nnsg) :: albvl_s , albvs_s , aldifl_s ,       &
                                 aldifs_s , aldirl_s , aldirs_s
    integer :: kolour , n , i , j
    character (len=64) :: subroutine_name='albedov'
    integer :: idindx=0
!
    call time_begin(subroutine_name,idindx)
    !
    ! =================================================================
    ! 1. set initial parameters
    ! =================================================================
    !
    ! 1.1 constants
    !
    ! Solar flux partitioned at wavelength of 0.7micr
    fsol1 = 0.5D0
    fsol2 = 0.5D0
    ! Short and long wave albedo for new snow
    snal0 = 0.95D0
    snal1 = 0.65D0
    ! Short and long wave albedo for sea ice
    sical0 = 0.6D0
    sical1 = 0.4D0
    !
    ! Desert seasonal albedo
    ! Works for Sahara desert and generally northern emisphere
    ! In souther emisphere only some points have this class
    !
    if ( ldesseas ) then
      if ( imon == 1 .or. imon == 2 .or. imon == 12 ) then
        solour(1) = 0.12D0
      endif        
      if ( imon == 3 .or. imon == 4 .or. imon == 5 ) then
        solour(1) = 0.15D0
      endif        
      if ( imon == 6 .or. imon == 7 .or. imon == 8) then
        solour(1) = 0.18D0
      endif        
      if ( imon == 9 .or. imon == 10 .or. imon == 11) then
        solour(1) = 0.15D0
      endif
    end if

    do i = istart , iend
      do j = jstart , jend
        do n = 1 , nnsg
          lveg(n,j,i) = veg2d1(n,i,j)
        end do
      end do
    end do

    !
    ! In depth, wt is frac of grid square covered by snow;
    ! depends on average snow depth, vegetation, etc.
    !
    call depth(jstart,jend,istart,iend)
    !
    ! 1.2  set default vegetation and albedo
    ! 
    do i = istart , iend
      do j = jstart , jend
        czeta = coszrs(j,i)
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
          ! 
          !================================================================
          !       2.   get albedo over land
          !================================================================
          !
          if ( ldimsk(n,j,i) == 2 ) then
            tdiffs = sts(n,j,i) - tzero
            tdiff = dmax1(tdiffs,d_zero)
            tdiffs = dmin1(tdiff,20.0D0)
            albgl = sical1 - 1.1D-2*tdiffs
            albgs = sical0 - 2.45D-2*tdiffs
            albg = fsol1*albgs + fsol2*albgl
            albgsd = albgs
            albgld = albgl
          else if ( ldimsk(n,j,i) == 1 ) then
            sfac = d_one - fseas(tgbrd(n,j,i))
            !
            ! ccm tests here on land mask for veg and soils data
            ! reduces albedo at low temps !!!!!
            ! should respond to moisture too (commented out) (pat, 27 oct 86)
            ! lncl(i) = lncl(lveg(i)) - seasf(lveg(i)) * sfac
            !
            albs = albvgs(lveg(n,j,i))
            albl = albvgl(lveg(n,j,i))
            !---------------------------------------------------------------
            if ( (lveg(n,j,i) < 12) .or. (lveg(n,j,i) > 15) ) then
              ! 2.1  bare soil albedos
              !      (soil albedo depends on moisture)
              kolour = kolsol(lveg(n,j,i))
              wet = ssw(n,j,i)/depuv(lveg(n,j,i))
              alwet = dmax1((11.0D0-40.0D0*wet),d_zero)*0.01D0
              alwet = dmin1(alwet,solour(kolour))
              albg = solour(kolour) + alwet
              albgs = albg
              albgl = d_two*albg
              ! higher nir albedos set diffuse albedo
              albgld = albgl
              albgsd = albgs
              albsd = albs
              albld = albl
              ! Dec. 15   albzn=0.85D0+d_one/(d_one+d_10*coszrs(j,i))
              ! Dec. 12, 2008
              albzn = d_one
              ! Dec. 15, 2008
              !
              ! leafless hardwood canopy: no or inverse zen dep
              if ( lveg(n,j,i) == 5 .and. sfac < 0.1D0 ) albzn = d_one
              ! multiply by zenith angle correction
              albs = albs*albzn
              albl = albl*albzn
              ! albedo over vegetation after zenith angle corr
              albvs_s(n) = albs
              albvl_s(n) = albl
            else if ( lveg(n,j,i) == 12 ) then
              ! 2.2   permanent ice sheet
              albgs = 0.8D0
              albgsd = 0.8D0
              albgl = 0.55D0
              albgld = 0.55D0
            else
              ! 2.3  inland water, swamps, rice paddies etc.
              albg = 0.05D0/(czeta+0.15D0)
              albgs = albg
              albgsd = albg
              albgl = albg
              albgld = albg
            end if
          end if
          ! ===================================================================
          ! 4.  correct for snow cover
          ! ===================================================================
          if ( sncv(n,j,i) > d_zero ) then
            ! Snow albedo depends on snow-age, zenith angle, and thickness
            ! of snow. snow albedoes for visible and ir solar rad visible
            ! albedo depends on snow age
            ! age gives reduction of visible rad snow albedo due to age
            cons = 0.2D0
            conn = 0.5D0
            age = (d_one-d_one/(d_one+snag(n,j,i)))
            ! sl helps control albedo zenith dependence
            sl = d_two
            sli = d_one/sl
            sl2 = d_two*sl
            ! snal0= new snow albedo for vis rad, sol zen le 6
            ! snal1= new snow albedo for long-wave rad
            dfalbs = snal0*(d_one-cons*age)
            ! czf corrects albedo of new snow for solar zenith
            cf1 = ((d_one+sli)/(d_one+sl2*coszrs(j,i))-sli)
            cff = dmax1(cf1,d_zero)
            czf = 0.4D0*cff*(d_one-dfalbs)
            dralbs = dfalbs + czf
            dfalbl = snal1*(d_one-conn*age)
            czf = 0.4D0*cff*(d_one-dfalbl)
            dralbl = dfalbl + czf
            if ( lncl(n,j,i) > 0.001D0 ) then
              ! effective albedo over vegetation with snow
              albl = (d_one-wt(n,j,i))*albl + dralbl*wt(n,j,i)
              albld = (d_one-wt(n,j,i))*albld + dfalbl*wt(n,j,i)
              albs = (d_one-wt(n,j,i))*albs + dralbs*wt(n,j,i)
              albsd = (d_one-wt(n,j,i))*albsd + dfalbs*wt(n,j,i)
            end if
            !----------------------------------------------------------------
            !         4.1  compute albedo for snow on bare ground
            !----------------------------------------------------------------
            albgs = (d_one-scvk(n,j,i))*albgs + dralbs*scvk(n,j,i)
            albgl = (d_one-scvk(n,j,i))*albgl + dralbl*scvk(n,j,i)
            albgsd = (d_one-scvk(n,j,i))*albgsd + dfalbs*scvk(n,j,i)
            albgld = (d_one-scvk(n,j,i))*albgld + dfalbl*scvk(n,j,i)
          end if
          !=====================================================================
          !       5.  albedo over open ocean
          !=====================================================================

          if ( ldimsk(n,j,i) == 0 ) then
            ! ocean albedo depends on zenith angle
            if ( czeta >= d_zero ) then
              ! albedo independent of wavelength
              albg = 0.05D0/(czeta+0.15D0)
              albgs = albg
              albgl = albg
              albgsd = 0.08D0
              albgld = 0.08D0
            end if
          end if
          !
          ! not part of albedo in the ccm
          !
          aldirs_s(n) = (d_one-lncl(n,j,i))*albgs + lncl(n,j,i)*albs
          aldirl_s(n) = (d_one-lncl(n,j,i))*albgl + lncl(n,j,i)*albl
          aldifs_s(n) = (d_one-lncl(n,j,i))*albgsd + lncl(n,j,i)*albsd
          aldifl_s(n) = (d_one-lncl(n,j,i))*albgld + lncl(n,j,i)*albld
        end do
        albvs(j,i) = albvs_s(1)
        albvl(j,i) = albvl_s(1)
        aldirs(j,i) = aldirs_s(1)
        aldirl(j,i) = aldirl_s(1)
        aldifs(j,i) = aldifs_s(1)
        aldifl(j,i) = aldifl_s(1)
        albdirs(1,j,i) = aldirs_s(1)
        albdifs(1,j,i) = aldifs_s(1)
        aemiss(j,i) = emiss2d(1,i,j)
        do n = 2 , nnsg
          albvs(j,i) = albvs(j,i) + albvs_s(n)
          albvl(j,i) = albvl(j,i) + albvl_s(n)
          aldirs(j,i) = aldirs(j,i) + aldirs_s(n)
          aldirl(j,i) = aldirl(j,i) + aldirl_s(n)
          aldifs(j,i) = aldifs(j,i) + aldifs_s(n)
          aldifl(j,i) = aldifl(j,i) + aldifl_s(n)
          albdirs(n,j,i) = aldirs_s(n)
          albdifs(n,j,i) = aldifs_s(n)
          aemiss(j,i) = aemiss(j,i) + emiss2d(n,i,j)
        end do
        albvs(j,i) = albvs(j,i)*rdnnsg
        albvl(j,i) = albvl(j,i)*rdnnsg
        aldirs(j,i) = aldirs(j,i)*rdnnsg
        aldirl(j,i) = aldirl(j,i)*rdnnsg
        aldifs(j,i) = aldifs(j,i)*rdnnsg
        aldifl(j,i) = aldifl(j,i)*rdnnsg
        aemiss(j,i) = aemiss(j,i)*rdnnsg
      end do
    end do
 
    call time_end(subroutine_name,idindx)

    contains

      function fseas(x)
        implicit none
        real(dp) :: fseas
        real(dp) , intent(in) :: x
        fseas = dmax1(d_zero,(d_one-0.0016D0*dmax1(298.0D0-x,d_zero)**d_two))
      end function fseas

  end subroutine albedov
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!   this subrout overwrites many of the soil constants
!   as a function of location(jlon,jlat)
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
  subroutine soilbc(jstart,jend,istart,iend)
    implicit none
    integer , intent(in) :: jstart , jend , istart , iend
    real(dp) :: ck , dmax , dmin , dmnor , phi0 , tweak1
    integer :: itex , n , i , j
!
!   ================================================================
!   new soils data as a fn of texture make porosity, soil suction,
!   hydraul conduc, wilting frac variables rather than consts
!   relfc is the ratio of field capacity to saturated water content,
!   defined so the rate of gravitational drainage at field
!   capacity is assumed to be 2 mm/day (baver et al., 1972)
!   ===============================================================
!
    character (len=64) :: subroutine_name='soilbc'
    integer :: idindx=0
!
    call time_begin(subroutine_name,idindx)
 
    do i = istart , iend
      do j = jstart , jend
        do n = 1 , nnsg
          if ( ldimsk(n,j,i) /= 0 ) then
            ! lveg is set in subr. interf
            freza(lveg(n,j,i)) = 0.15D0*deprv(lveg(n,j,i))
            frezu(lveg(n,j,i)) = 0.15D0*depuv(lveg(n,j,i))
            itex = iexsol(lveg(n,j,i))
            texrat(n,j,i) = skrat(itex)
            porsl(n,j,i) = xmopor(itex)
            xkmx(n,j,i) = xmohyd(itex)
            bsw(n,j,i) = bee(itex)
            bfc(n,j,i) = 5.8D0 - bsw(n,j,i)*(0.8D0+0.12D0*(bsw(n,j,i)-d_four)* &
                       dlog10(1.0D2*xkmx(n,j,i)))
            phi0 = xmosuc(itex)
            dmax = bsw(n,j,i)*phi0*xkmx(n,j,i)/porsl(n,j,i)
            dmin = 1.0D-3
            dmnor = 1550.0D0*dmin/dmax
            tweak1 = (bsw(n,j,i)*(bsw(n,j,i)-6.0D0)+10.3D0) / &
                     (bsw(n,j,i)*bsw(n,j,i)+40.0D0*bsw(n,j,i))
            ck = (d_one+dmnor)*tweak1*0.23D0/0.02356D0
            evmx0(n,j,i) = 1.02D0*dmax*ck/dsqrt(depuv(lveg(n,j,i))*deprv(lveg(n,j,i)))
            gwmx0(n,j,i) = depuv(lveg(n,j,i))*porsl(n,j,i)
            gwmx1(n,j,i) = deprv(lveg(n,j,i))*porsl(n,j,i)
            gwmx2(n,j,i) = deptv(lveg(n,j,i))*porsl(n,j,i)
            wiltr(n,j,i) = xmowil(itex)
            ! force irrigated crop to be at field capacity
            relfc(n,j,i) = xmofc(itex)
          end if
        end do
      end do
    end do
 
    call time_end(subroutine_name,idindx)
  end subroutine soilbc
!
! For albedov
! This is taken from subroutine interf so that radiation can be
! called in tend (not mtrxbats).
!
  subroutine slice1D(jstart,jend,istart,iend)
    implicit none
    integer , intent(in) :: jstart , jend , istart , iend
!
    real(dp) :: amxtem , sfac
    integer :: n , i , j
    do i = istart , iend
      do j = jstart , jend
        do n = 1 , nnsg
          ldimsk(n,j,i) = ocld2d(n,i,j)
          sfice(n,j,i) = sice2d(n,i,j)
          tgbrd(n,j,i) = tgb2d(n,i,j)
          ssw(n,j,i) = ssw2d(n,i,j)
          lveg(n,j,i) = veg2d1(n,i,j)
          oveg(n,j,i) = lveg(n,j,i)
          if ( ldimsk(n,j,i) == 2 ) lveg(n,j,i) = 12
          amxtem = dmax1(298.0D0-tgbrd(n,j,i),d_zero)
          sfac = d_one - dmax1(d_zero,d_one-0.0016D0*amxtem**d_two)
          lncl(n,j,i) = mfcv(lveg(n,j,i)) - seasf(lveg(n,j,i))*sfac
          sts(n,j,i) = thatm(i,kz,j)-6.5D-3*regrav*(ht1(n,i,j)-ht(i,j))
          sncv(n,j,i) = scv2d(n,i,j)
          snag(n,j,i) = sag2d(n,i,j)
        end do
      end do
    end do
!
  end subroutine slice1D
!
end module mod_bats_mtrxbats
