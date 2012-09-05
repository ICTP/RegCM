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

  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_runparams , only : iqv , iocnflx , iocncpl
  use mod_mppparam
  use mod_mpmessage
  use mod_constants
  use mod_service
  use mod_bats_common
  use mod_bats_internal
  use mod_bats_leaftemp , only : fseas
  use mod_bats_lake
  use mod_bats_bndry
  use mod_bats_drag
  use mod_bats_mppio
  use mod_bats_zengocn

  private

  !
  ! Solar flux partitioned at wavelength of 0.7micr
  !
  real(rk8) , parameter :: fsol1 = 0.5D0
  real(rk8) , parameter :: fsol2 = 0.5D0
  !
  ! Short and long wave albedo for new snow
  !
  real(rk8) , parameter :: snal0 = 0.95D0
  real(rk8) , parameter :: snal1 = 0.65D0
  !
  ! Short and long wave albedo for sea ice
  !
  real(rk8) , parameter :: sical0 = 0.6D0
  real(rk8) , parameter :: sical1 = 0.4D0
  !

  public :: interf , initb , mtrxbats , albedobats

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
 !               albedobats
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
  subroutine mtrxbats(ktau)
!
    implicit none
!
    integer(8) , intent(in) :: ktau
    character (len=64) :: subroutine_name='mtrxbats'
    integer(ik4) :: idindx=0
!
    call time_begin(subroutine_name,idindx)
!
!---------------------------------------------------------------------
!
!   Excange from model to BATS

    call interf(1,ktau)

!   Calculate surface fluxes and hydrology budgets

    call soilbc

!   Albedo is calculated in mod_tendency

    call bndry

!   Zeng ocean flux model
    if ( iocnflx == 2 ) call zengocndrv(ktau)

!   Hostetler lake model for every BATS timestep at lake points
    if ( llake ) then
      call lakedrv
    endif

!   Accumulate quantities for energy and moisture budgets

    call interf(2,ktau)
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
  subroutine initb
    implicit none
! 
    integer(ik4) :: i , is , itex , j , n , nlveg
!
    character (len=64) :: subroutine_name='initb'
    integer(ik4) :: idindx=0
!
    call time_begin(subroutine_name,idindx)
!
!   initialize hostetler lake model
!
    if ( llake ) call initlake
 
    do i = ici1 , ici2
      do j = jci1 , jci2
        do n = 1 , nnsg
          tgrd(n,j,i) = tground2(j,i)
          tgbrd(n,j,i) = tground2(j,i)
          taf(n,j,i) = tground2(j,i)
          tlef(n,j,i) = tground2(j,i)
          if ( ldmsk1(n,j,i) == 2 ) then
            if ( lseaice .or. llake ) then
              sfice(n,j,i) = d_1000
            end if
            nlveg = 12
          else if ( ldmsk1(n,j,i) == 1 ) then
            sncv(n,j,i) = dmax1(sncv(n,j,i),d_zero)
            nlveg = iveg1(n,j,i)
          else
            nlveg = iveg1(n,j,i)
          end if
!         Initialize soil moisture in the 3 layers
          is = iveg1(n,j,i)
          itex = iexsol(nlveg)
          tsw(n,j,i) = deptv(nlveg)*xmopor(itex)*slmo(is)
          rsw(n,j,i) = deprv(nlveg)*xmopor(itex)*slmo(is)
          ssw(n,j,i) = depuv(nlveg)*xmopor(itex)*slmo(is)
          gwet(n,j,i) = d_half
        end do
      end do
    end do
! 
    call time_end(subroutine_name,idindx)
!
  end subroutine initb
!
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!  this subroutine interfaces regcm and bats variables
!
!  ivers = 1 ,   regcm --> bats
!  ivers = 2 ,   bats --> regcm
!
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
  subroutine interf(ivers,ktau)
!
    implicit none
!
    integer(ik4) , intent (in) :: ivers
    integer(8) , intent(in) :: ktau
!
    real(rk8) :: facb , facs , fact , factuv , facv , fracb ,  &
                fracs , fracv , hl , mmpd , rh0 , satvp ,     &
                solvt , wpm2 , p0 , qs0 , ts0
    integer(ik4) :: i , j , n , nnn
    real(rk4) :: real_4
!
    character (len=64) :: subroutine_name='interf'
    integer(ik4) :: idindx=0
!
    call time_begin(subroutine_name,idindx)
 
    if ( ivers == 1 ) then ! regcm --> bats

      do i = ici1 , ici2
        do j = jci1 , jci2
          p0 = (sfps(j,i)+ptop)*d_1000
          qs0 = qxatm(j,i,kz,iqv)/(d_one+qxatm(j,i,kz,iqv))
          ts0 = thatm(j,i,kz)
          hl = lh0 - lh1*(ts0-tzero)
          satvp = lsvp1*dexp(lsvp2*hl*(d_one/tzero-d_one/ts0))
          rh0 = dmax1(qs0/(ep2*satvp/(p0*0.01D0-satvp)),d_zero)
          do n = 1 , nnsg
            qs(n,j,i) = qs0
            sts(n,j,i) = ts0-lrate*regrav*(ht1(n,j,i)-ht(j,i))
            sfcp(n,j,i) = p0*(sts(n,j,i)/ts0)
            hl = lh0 - lh1*(sts(n,j,i)-tzero)
            satvp = lsvp1*dexp(lsvp2*hl*(d_one/tzero-d_one/sts(n,j,i)))
            qs(n,j,i) = dmax1(rh0*ep2*satvp/(sfcp(n,j,i)*0.01D0-satvp),d_zero)
            rhs(n,j,i) = sfcp(n,j,i)/(rgas*sts(n,j,i))
            prcp(n,j,i) = pptnc(j,i) + pptc(j,i)
            totpr(j,i) = pptnc(j,i) + pptc(j,i)
            !
            ! quantities stored on 2d surface array for bats use only
            !
            sent(n,j,i) = hfx(j,i)
            evpr(n,j,i) = qfx(j,i)
            lveg(n,j,i) = iveg1(n,j,i)
            if ( ldmsk1(n,j,i) == 2 ) lveg(n,j,i) = 12
            lncl(n,j,i) = mfcv(lveg(n,j,i)) - &
                          seasf(lveg(n,j,i))*fseas(tgbrd(n,j,i),lveg(n,j,i))
            zh(n,j,i) = hgt(j,i,kz)
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
 
          usw(j,i) = uatm(j,i,kz)
          vsw(j,i) = vatm(j,i,kz)
          solvt = solvd(j,i) + solvs(j,i)
          if ( solvt > d_zero ) then
            fracd(j,i) = solvd(j,i)/solvt
          else
            fracd(j,i) = 0.2D0
          end if
        end do
      end do

    else if ( ivers == 2 ) then ! bats --> regcm2d
 
      do i = ici1 , ici2
        do j = jci1 , jci2
          uvdrag(j,i) = d_zero
          hfx(j,i) = d_zero
          qfx(j,i) = d_zero
          tground2(j,i) = d_zero
          tground1(j,i) = d_zero
          tgbb(j,i) = d_zero
          if ( lchem ) then
            ssw2da(j,i) = d_zero
            sdeltk2d(j,i) = d_zero
            sdelqk2d(j,i) = d_zero
            sfracv2d(j,i) = d_zero
            sfracb2d(j,i) = d_zero
            sfracs2d(j,i) = d_zero
            svegfrac2d(j,i) = d_zero
          end if

          do n = 1 , nnsg
            uvdrag(j,i) = uvdrag(j,i) + drag(n,j,i)
            hfx(j,i) = hfx(j,i) + sent(n,j,i)
            qfx(j,i) = qfx(j,i) + evpr(n,j,i)
            tground2(j,i) = tground2(j,i) + tgrd(n,j,i)
            tground1(j,i) = tground1(j,i) + tgrd(n,j,i)
            if ( lchem ) then
              ssw2da(j,i) = ssw2da(j,i) + ssw(n,j,i)
              sdeltk2d(j,i) = sdeltk2d(j,i) + delt(n,j,i)
              sdelqk2d(j,i) = sdelqk2d(j,i) + delq(n,j,i)
              sfracv2d(j,i) = sfracv2d(j,i) + sigf(n,j,i)
              sfracb2d(j,i) = sfracb2d(j,i)+ &
                              (d_one-lncl(n,j,i))*(d_one-scvk(n,j,i))
              sfracs2d(j,i) = sfracs2d(j,i) + &
                       lncl(n,j,i)*wt(n,j,i) + (d_one-lncl(n,j,i))*scvk(n,j,i)
              svegfrac2d(j,i) = svegfrac2d(j,i) + lncl(n,j,i)
            end if
            if ( iocnflx == 1 .or. &
                (iocnflx == 2 .and. ldmsk1(n,j,i) /= 0 ) ) then
              tgbb(j,i) = tgbb(j,i) +                 &
                         ((d_one-lncl(n,j,i))*tgrd(n,j,i)**d_four +  &
                         lncl(n,j,i)*tlef(n,j,i)**d_four)**d_rfour
            else
              tgbb(j,i) = tgbb(j,i) + tgrd(n,j,i)
            end if
            if ( ldmsk1(n,j,i) == 0 ) then
              ssw(n,j,i)   = dmissval
              rsw(n,j,i)   = dmissval
              tsw(n,j,i)   = dmissval
              trnof(n,j,i) = dmissval
              srnof(n,j,i) = dmissval
              sncv(n,j,i)  = dmissval
            end if
          end do
          uvdrag(j,i) = uvdrag(j,i)*rdnnsg
          hfx(j,i) = hfx(j,i)*rdnnsg
          qfx(j,i) = qfx(j,i)*rdnnsg
          tground2(j,i) = tground2(j,i)*rdnnsg
          tground1(j,i) = tground1(j,i)*rdnnsg
          tgbb(j,i) = tgbb(j,i)*rdnnsg

          if ( lchem ) then
            ssw2da(j,i) = ssw2da(j,i)*rdnnsg
            sdeltk2d(j,i) = sdeltk2d(j,i)*rdnnsg
            sdelqk2d(j,i) = sdelqk2d(j,i)*rdnnsg
            sfracv2d(j,i) = sfracv2d(j,i)*rdnnsg
            sfracb2d(j,i) = sfracb2d(j,i)*rdnnsg
            sfracs2d(j,i) = sfracs2d(j,i)*rdnnsg
            svegfrac2d(j,i) = svegfrac2d(j,i)*rdnnsg
          end if
          do n = 1 , nnsg
            evpa(n,j,i) = evpa(n,j,i) + dtbat*evpr(n,j,i)
            sena(n,j,i) = sena(n,j,i) + dtbat*sent(n,j,i)
            if ( dabs(srnof(n,j,i)) > 1.0D-10 ) then
              srfrna(n,j,i) = srfrna(n,j,i) + srnof(n,j,i)*dtbat
            end if
            if ( dabs(srnof(n,j,i))  > 1.0D-10 .and. &
                 dabs(trnof(n,j,i))   > 1.0D-10 ) then
              runoff(n,j,i) = runoff(n,j,i) + &
                     (trnof(n,j,i)-srnof(n,j,i))*dtbat
            end if
          end do
          !
          ! quantities stored on 2d surface array for bats use only
          !
          prca(j,i) = prca(j,i) + dtbat*pptc(j,i)
          prnca(j,i) = prnca(j,i) + dtbat*pptnc(j,i)
          flwa(j,i) = flwa(j,i) + dtbat*flw(j,i)
          flwda(j,i) = flwda(j,i) + dtbat*flwd(j,i)
          fswa(j,i) = fswa(j,i) + dtbat*fsw(j,i)
          svga(j,i) = svga(j,i) + dtbat*sabveg(j,i)
          sina(j,i) = sina(j,i) + dtbat*sinc(j,i)
        end do
      end do

      do i = ici1 , ici2
        do j = jci1 , jci2
          fbat(j,i,u10m_o)   = 0.0
          fbat(j,i,v10m_o)   = 0.0
          fbat(j,i,tg_o)     = 0.0
          fbat(j,i,t2m_o)    = 0.0
          fbat(j,i,q2m_o)    = 0.0
          fbat(j,i,aldirs_o) = 0.0
          fbat(j,i,aldifs_o) = 0.0
          do n = 1 , nnsg
            if ( ldmsk1(n,j,i) /= 0 ) then
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
            fsub(n,j,i,tg_s)   = real(tgrd(n,j,i))
            fsub(n,j,i,u10m_s) = real(u10m(n,j,i))
            fsub(n,j,i,v10m_s) = real(v10m(n,j,i))
            fsub(n,j,i,t2m_s)  = real(t2m(n,j,i))
            fsub(n,j,i,q2m_s)  = real(q2m(n,j,i))
 
            fbat(j,i,u10m_o) = fbat(j,i,u10m_o) + real(u10m(n,j,i))
            fbat(j,i,v10m_o) = fbat(j,i,v10m_o) + real(v10m(n,j,i))
            fbat(j,i,t2m_o)  = fbat(j,i,t2m_o) + real(t2m(n,j,i))
            fbat(j,i,q2m_o)  = real(fbat(j,i,q2m_o) + q2m(n,j,i))
            fbat(j,i,tg_o)   = fbat(j,i,tg_o) + real(tgrd(n,j,i))
          end do
          fbat(j,i,u10m_o)   = fbat(j,i,u10m_o)*rrnnsg
          fbat(j,i,v10m_o)   = fbat(j,i,v10m_o)*rrnnsg
          fbat(j,i,t2m_o)    = fbat(j,i,t2m_o)*rrnnsg
          fbat(j,i,q2m_o)    = fbat(j,i,q2m_o)*rrnnsg
          fbat(j,i,tg_o)     = fbat(j,i,tg_o)*rrnnsg
          fbat(j,i,aldirs_o) = real(aldirs(j,i))
          fbat(j,i,aldifs_o) = real(aldifs(j,i))
 
          fbat(j,i,tgmx_o) = amax1(fbat(j,i,tgmx_o),fbat(j,i,tg_o))
          fbat(j,i,tgmn_o) = amin1(fbat(j,i,tgmn_o),fbat(j,i,tg_o))
          fbat(j,i,t2mx_o) = amax1(fbat(j,i,t2mx_o),fbat(j,i,t2m_o))
          fbat(j,i,t2mn_o) = amin1(fbat(j,i,t2mn_o),fbat(j,i,t2m_o))
          fbat(j,i,w10x_o) = amax1(fbat(j,i,w10x_o), &
                          sqrt(fbat(j,i,u10m_o)**2.0+fbat(j,i,v10m_o)**2.0))
          real_4 = real((pptnc(j,i)+pptc(j,i)))
          fbat(j,i,pcpx_o) = amax1(fbat(j,i,pcpx_o),real_4)
          fbat(j,i,pcpa_o) = fbat(j,i,pcpa_o) + real_4/fdaysrf
          fbat(j,i,tavg_o) = fbat(j,i,tavg_o)+fbat(j,i,t2m_o)/fdaysrf
          real_4 = real((sfps(j,i)+ptop)*d_10)
          fbat(j,i,psmn_o) = amin1(fbat(j,i,psmn_o),real_4)
          if ( fsw(j,i) > 120.0D0 ) then
            fbat(j,i,sund_o) = fbat(j,i,sund_o) + real(dtbat)
            fbat(j,i,sunt_o) = fbat(j,i,sunt_o) + real(dtbat)
          end if
          pptnc(j,i) = d_zero
          pptc(j,i)  = d_zero
        end do
      end do

      if ( mod(ktau+1,kbats) == 0 .or. ktau == 0 ) then
        if ( ktau == 0 .and. debug_level > 2 ) then
          mmpd = secpd/dtbat
          wpm2 = d_one/dtbat
        else if ( ktau+1 == kbats .and. debug_level > 2 ) then
          mmpd = houpd/(srffrq-xdtsec/secph)
          wpm2 = d_one/((srffrq-xdtsec/secph)*secph)
        else
          mmpd = houpd/srffrq
          wpm2 = d_one/(srffrq*secph)
        end if
        do i = ici1 , ici2
          do j = jci1 , jci2
            fbat(j,i,drag_o) = 0.0
            fbat(j,i,evpa_o) = 0.0
            fbat(j,i,sena_o) = 0.0
            do n = 1 , nnsg
              if ( ldmsk1(n,j,i) /= 0 ) then
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
              fsub(n,j,i,drag_s) = real(drag(n,j,i))
              fsub(n,j,i,evpa_s) = real(evpa(n,j,i)*mmpd)
              fsub(n,j,i,sena_s) = real(sena(n,j,i)*wpm2)
              fsub(n,j,i,tpr_s)  = real((prnca(j,i)+prca(j,i))*mmpd)
              fsub(n,j,i,prcv_s) = real(prca(j,i)*mmpd)
              fsub(n,j,i,ps_s)   = real(sfcp(n,j,i)*0.01D0)
 
              fbat(j,i,drag_o) = real(fbat(j,i,drag_o) + drag(n,j,i))
              fbat(j,i,evpa_o) = real(fbat(j,i,evpa_o) + evpa(n,j,i))
              fbat(j,i,sena_o) = real(fbat(j,i,sena_o) + sena(n,j,i))
            end do
            fbat(j,i,tpr_o)  = real((prnca(j,i)+prca(j,i))*mmpd)
            fbat(j,i,drag_o) = fbat(j,i,drag_o)*rrnnsg
            fbat(j,i,evpa_o) = fbat(j,i,evpa_o)*rrnnsg*real(mmpd)
            fbat(j,i,sena_o) = fbat(j,i,sena_o)*rrnnsg*real(wpm2)
            fbat(j,i,flwa_o) = real(flwa(j,i)*wpm2)
            fbat(j,i,fswa_o) = real(fswa(j,i)*wpm2)
            fbat(j,i,flwd_o) = real(flwda(j,i)*wpm2)
            fbat(j,i,sina_o) = real(sina(j,i)*wpm2)
            fbat(j,i,prcv_o) = real(prca(j,i)*mmpd)
            fbat(j,i,ps_o) = real((sfps(j,i)+ptop)*d_10)
            fbat(j,i,zpbl_o) = real(hpbl(j,i))
 
            fbat(j,i,tlef_o) = 0.0
            fbat(j,i,ssw_o)  = 0.0
            fbat(j,i,rsw_o)  = 0.0
            fbat(j,i,rnos_o) = 0.0
            fbat(j,i,scv_o) = 0.0
            nnn = 0
            do n = 1 , nnsg
              if ( ldmsk1(n,j,i) /= 0 ) then
                fbat(j,i,tlef_o) = fbat(j,i,tlef_o) + real(tlef(n,j,i))
                fbat(j,i,ssw_o)  = fbat(j,i,ssw_o) + real(ssw(n,j,i))
                fbat(j,i,rsw_o)  = fbat(j,i,rsw_o) + real(rsw(n,j,i))
                fbat(j,i,rnos_o) = fbat(j,i,rnos_o) + real(srfrna(n,j,i))
                fbat(j,i,scv_o)  = fbat(j,i,scv_o) + real(sncv(n,j,i))
                fsub(n,j,i,tlef_s) = real(tlef(n,j,i))
                fsub(n,j,i,ssw_s)  = real(ssw(n,j,i))
                fsub(n,j,i,rsw_s)  = real(rsw(n,j,i))
                fsub(n,j,i,rnos_s) = real(srfrna(n,j,i)*mmpd)
                fsub(n,j,i,scv_s)  = real(sncv(n,j,i))
                nnn = nnn + 1
              else
                fsub(n,j,i,tlef_s) = smissval
                fsub(n,j,i,ssw_s)  = smissval
                fsub(n,j,i,rsw_s)  = smissval
                fsub(n,j,i,rnos_s) = smissval
                fsub(n,j,i,scv_s)  = smissval
              end if
            end do
            if ( nnn >= max0(nnsg/2,1) ) then
              fbat(j,i,tlef_o) = fbat(j,i,tlef_o)/real(nnn)
              fbat(j,i,ssw_o)  = fbat(j,i,ssw_o)/real(nnn)
              fbat(j,i,rsw_o)  = fbat(j,i,rsw_o)/real(nnn)
              fbat(j,i,rnos_o) = fbat(j,i,rnos_o)/real(nnn)*real(mmpd)
              fbat(j,i,scv_o)  = fbat(j,i,scv_o)/real(nnn)
            else
              fbat(j,i,tlef_o) = smissval
              fbat(j,i,ssw_o)  = smissval
              fbat(j,i,rsw_o)  = smissval
              fbat(j,i,rnos_o) = smissval
              fbat(j,i,scv_o)  = smissval
            end if
            !
            ! reset accumulation arrays to zero
            !
            do n = 1 , nnsg
              evpa(n,j,i) = d_zero
              srfrna(n,j,i) = d_zero
              sena(n,j,i) = d_zero
            end do
            prnca(j,i) = d_zero
            prca(j,i) = d_zero
            flwa(j,i) = d_zero
            flwda(j,i) = d_zero
            fswa(j,i) = d_zero
            svga(j,i) = d_zero
            sina(j,i) = d_zero
          end do
        end do
      end if
    end if
    call time_end(subroutine_name,idindx)
!
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
  subroutine albedobats(imon)
    implicit none
    integer(ik4) , intent (in) :: imon
!
    real(rk8) :: age , albg , albgl , albgld , albgs , albgsd , albl ,  &
               albld , albs , albsd , albzn , alwet , cf1 , cff ,     &
               conn , cons , czeta , czf , dfalbl , dfalbs , dralbl , &
               dralbs , sfac , sl , sl2 , sli , tdiff , tdiffs , wet
    real(rk8) , dimension(nnsg) :: albvl_s , albvs_s , aldifl_s ,       &
                                 aldifs_s , aldirl_s , aldirs_s
    integer(ik4) :: kolour , n , i , j
    character (len=64) :: subroutine_name='albedobats'
    integer(ik4) :: idindx=0
!
    call time_begin(subroutine_name,idindx)
    !
    ! =================================================================
    ! 1. set initial parameters
    ! =================================================================
    !
    do i = ici1 , ici2
      do j = jci1 , jci2
        do n = 1 , nnsg
          lveg(n,j,i) = iveg1(n,j,i)
          if ( ldmsk1(n,j,i) == 2 ) lveg(n,j,i) = 12
        end do
      end do
    end do
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
    !
    ! In depth, wt is frac of grid square covered by snow;
    ! depends on average snow depth, vegetation, etc.
    !
    call depth
    !
    ! 1.2  set default vegetation and albedo
    ! 
    do i = ici1 , ici2
      do j = jci1 , jci2
        czeta = coszrs(j,i)
        do n = 1 , nnsg
          lncl(n,j,i) = mfcv(lveg(n,j,i)) - &
                        seasf(lveg(n,j,i))*fseas(tgbrd(n,j,i),lveg(n,j,i))
          sts(n,j,i) = thatm(j,i,kz)-lrate*regrav*(ht1(n,j,i)-ht(j,i))
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
          !       2.   get albedo
          !================================================================
          !
          if ( ldmsk1(n,j,i) == 0 ) then
            ! ocean albedo depends on zenith angle
            if ( czeta >= d_zero ) then
              ! albedo independent of wavelength
              albg = 0.05D0/(czeta+0.15D0)
              albgs = albg
              albgl = albg
              albgsd = 0.08D0
              albgld = 0.08D0
            else
              albg = 0.05D0
              albgs = albg
              albgl = albg
              albgsd = 0.08D0
              albgld = 0.08D0
            end if
          else if ( ldmsk1(n,j,i) == 2 ) then
            ! Ice over ocean or lake
            tdiffs = sts(n,j,i) - tzero
            tdiff = dmax1(tdiffs,d_zero)
            tdiffs = dmin1(tdiff,20.0D0)
            albgl = sical1 - 1.1D-2*tdiffs
            albgs = sical0 - 2.45D-2*tdiffs
            albg = fsol1*albgs + fsol2*albgl
            albgsd = albgs
            albgld = albgl
          else if ( ldmsk1(n,j,i) == 1 ) then
            ! Land
            sfac = d_one - fseas(tgbrd(n,j,i),lveg(n,j,i))
            !
            ! ccm tests here on land mask for veg and soils data
            ! reduces albedo at low temps !!!!!
            ! should respond to moisture too (commented out) (pat, 27 oct 86)
            ! lncl(i) = lncl(iveg1(i)) - seasf(iveg1(i)) * sfac
            !
            albs = albvgs(iveg1(n,j,i))
            albl = albvgl(iveg1(n,j,i))
            !---------------------------------------------------------------
            if ( (iveg1(n,j,i) < 12) .or. (iveg1(n,j,i) > 15) ) then
              ! 2.1  bare soil albedos
              !      (soil albedo depends on moisture)
              kolour = kolsol(iveg1(n,j,i))
              wet = ssw(n,j,i)/depuv(iveg1(n,j,i))
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
              if ( iveg1(n,j,i) == 5 .and. sfac < 0.1D0 ) albzn = d_one
              ! multiply by zenith angle correction
              albs = albs*albzn
              albl = albl*albzn
              ! albedo over vegetation after zenith angle corr
              albvs_s(n) = albs
              albvl_s(n) = albl
            else if ( iveg1(n,j,i) == 12 ) then
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
          ! ==============================================================
          ! 4.  correct for snow cover
          ! ==============================================================
          if ( ldmsk1(n,j,i) > 0 ) then
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
          end if
          !
          ! not part of albedo in the ccm
          !
          aldirs_s(n) = (d_one-lncl(n,j,i))*albgs + lncl(n,j,i)*albs
          aldirl_s(n) = (d_one-lncl(n,j,i))*albgl + lncl(n,j,i)*albl
          aldifs_s(n) = (d_one-lncl(n,j,i))*albgsd + lncl(n,j,i)*albsd
          aldifl_s(n) = (d_one-lncl(n,j,i))*albgld + lncl(n,j,i)*albld
        end do
        albvs(j,i)  = sum(albvs_s)*rdnnsg
        albvl(j,i)  = sum(albvl_s)*rdnnsg
        aldirs(j,i) = sum(aldirs_s)*rdnnsg
        aldirl(j,i) = sum(aldirl_s)*rdnnsg
        aldifs(j,i) = sum(aldifs_s)*rdnnsg
        aldifl(j,i) = sum(aldifl_s)*rdnnsg
      end do
    end do
    aemiss = sum(emiss,1)*rdnnsg
    call time_end(subroutine_name,idindx)

  end subroutine albedobats
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!   this subrout overwrites many of the soil constants
!   as a function of location(jlon,jlat)
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
  subroutine soilbc
    implicit none
    real(rk8) :: ck , dmax , dmin , dmnor , phi0 , tweak1
    integer(ik4) :: itex , n , i , j
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
    integer(ik4) :: idindx=0
!
    call time_begin(subroutine_name,idindx)
 
    do i = ici1 , ici2
      do j = jci1 , jci2
        do n = 1 , nnsg
          if ( ldmsk1(n,j,i) /= 0 ) then
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
            evmx0(n,j,i) = 1.02D0*dmax*ck / &
                 dsqrt(depuv(lveg(n,j,i))*deprv(lveg(n,j,i)))
            gwmx0(n,j,i) = depuv(lveg(n,j,i))*porsl(n,j,i)
            gwmx1(n,j,i) = deprv(lveg(n,j,i))*porsl(n,j,i)
            gwmx2(n,j,i) = deptv(lveg(n,j,i))*porsl(n,j,i)
            wiltr(n,j,i) = xmowil(itex)
            ! force irrigated crop to be at field capacity
            relfc(n,j,i) = xmofc(itex)
            ! Imported Lara Kuepper's Irrigated Crop modification from RegCM3
            ! see Kueppers et al. (2008)
            ! relaw is between field capacity and wilting point
            relaw(n,j,i) = 0.75*(xmofc(itex)-xmowil(itex))+xmowil(itex)
          end if
        end do
      end do
    end do
 
    call time_end(subroutine_name,idindx)
  end subroutine soilbc
!
end module mod_bats_mtrxbats
