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
  use mod_runparams , only : iqv , iocnflx , iocncpl , ichem , &
                             iemiss , ksrf , xmonth
  use mod_mppparam
  use mod_mpmessage
  use mod_constants
  use mod_service
  use mod_bats_common
  use mod_bats_internal
  use mod_bats_lake
  use mod_bats_bndry
  use mod_bats_drag
  use mod_bats_leaftemp
  use mod_bats_mppio
  use mod_bats_zengocn
  use mod_bats_coare
  use mod_outvars

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
    implicit none
    integer(ik8) , intent(in) :: ktau
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'mtrxbats'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
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
    if ( iocnflx == 2 ) call zengocndrv

!   Hostetler lake model for every BATS timestep at lake points
    if ( llake ) then
      call lakedrv
    endif

!   Coare bulk flux algorithm
    if ( iocncpl == 1 ) call coare3_drv

!   Accumulate quantities for energy and moisture budgets

    call interf(2,ktau)
!
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
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
    integer(ik4) :: i , itex , j , n , nlveg
    logical , parameter :: snowhack = .false.
!
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'initb'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    !
    ! initialize hostetler lake model
    !
    if ( llake ) call initlake
 
    do i = ici1 , ici2
      do j = jci1 , jci2
        do n = 1 , nnsg
          tgrd(n,j,i) = tground2(j,i)
          tgbrd(n,j,i) = tground2(j,i)
          taf(n,j,i) = tground2(j,i)
          tlef(n,j,i) = tground2(j,i)
          nlveg = iveg1(n,j,i)
          itex  = iexsol(nlveg)
          if ( ldmsk1(n,j,i) == 2 ) then
            if ( iemiss == 1 ) emiss(n,j,i) = 0.97D0
            nlveg = 12
          end if
          if ( ldmsk1(n,j,i) > 0 ) then
            if ( snowam(j,i) > d_zero .and. snowam(j,i) < dmissval ) then
              sncv(n,j,i) = snowam(j,i)
            else
              if ( snowhack ) then
                if ( ht1(n,j,i) > 2000.0D0 ) then
                  sncv(n,j,i) = 0.01D0
                else
                  sncv(n,j,i) = d_zero
                end if
              else
                sncv(n,j,i) = d_zero
              end if
            end if
          end if
!         Initialize soil moisture in the 3 layers
          tsw(n,j,i) = deptv(nlveg)*xmopor(itex)*slmo(nlveg)
          rsw(n,j,i) = deprv(nlveg)*xmopor(itex)*slmo(nlveg)
          ssw(n,j,i) = depuv(nlveg)*xmopor(itex)*slmo(nlveg)
          gwet(n,j,i) = d_half
        end do
      end do
    end do
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
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
    implicit none
    integer(ik4) , intent (in) :: ivers
    integer(ik8) , intent(in) :: ktau
!
    real(rk8) :: facb , facs , fact , factuv , facv , fracb ,  &
                fracs , fracv , hl , rh0 , satvp ,     &
                solvt , p0 , qs0 , ts0
    integer(ik4) :: i , j , n , icemsk
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'interf'
    integer(ik4) , save :: idindx = 0
    integer(ik4) :: ierr
    call time_begin(subroutine_name,idindx)
#endif
 
    if ( ivers == 1 ) then ! regcm --> bats

      call fseas(tgbrd)

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
            if ( ldmsk1(n,j,i) == 2 ) then
              lveg(n,j,i) = 12
              if ( iemiss == 1 ) emiss(n,j,i) = 0.97D0
            else if ( ldmsk1(n,j,i) == 0 ) then
              if ( iemiss == 1 ) emiss(n,j,i) = 0.995D0
            end if
            lncl(n,j,i) = mfcv(lveg(n,j,i)) - seasf(lveg(n,j,i))*aseas(n,j,i)
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

      if ( iocnflx == 1 ) then
        do i = ici1 , ici2
          do j = jci1 , jci2
            do n = 1 , nnsg
              if ( ldmsk1(n,j,i) == 0 .and. cplmsk(j,i) == 0 ) tgrd(n,j,i) = tground2(j,i)
            end do
          end do
        end do
      end if

    else if ( ivers == 2 ) then ! bats --> regcm2d
 
      ! Re-create the land sea mask to account for ice sheet melting

      if ( lseaice .or. llake ) then
        do i = ici1 , ici2
          do j = jci1 , jci2
            if ( ldmsk(j,i) == 2 ) then
              icemsk = 0
              do n = 1 , nnsg
                if ( ldmsk1(n,j,i) > 1 ) icemsk = icemsk + 1
              end do
              if ( icemsk <= nnsg/2 ) then
                ldmsk(j,i) = 0
              end if
            end if
          end do
        end do
      end if

      ! Those are needed elsewhere in the model (pbl,cumulus,chem,etc)

      do i = ici1 , ici2
        do j = jci1 , jci2
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
              if ( iocnflx == 1 .and. cplmsk(j,i) == 0 ) then
                fact = z2fra(n,j,i)/zlgocn(n,j,i)
                factuv = z10fra(n,j,i)/zlgocn(n,j,i)
                u10m(n,j,i) = usw(j,i)*(d_one-factuv)
                v10m(n,j,i) = vsw(j,i)*(d_one-factuv)
                t2m(n,j,i) = sts(n,j,i) - delt(n,j,i)*fact
                q2m(n,j,i) = qs(n,j,i) - delq(n,j,i)*fact
              end if
            end if
          end do
        end do
      end do
 
#ifdef DEBUG
      ierr = 0
      do i = ici1 , ici2
        do j = jci1 , jci2
          do n = 1 , nnsg
            if ( tgrd(n,j,i) < 150.0D0 ) then
              write(stderr,*) 'Likely error: Surface temperature too low'
              write(stderr,*) 'J   = ',global_dot_jstart+j
              write(stderr,*) 'I   = ',global_dot_istart+i
              write(stderr,*) 'VAL = ',tgrd(n,j,i)
              ierr = ierr + 1
            end if
          end do
        end do
      end do
      if ( ierr /= 0 ) then
        call fatal(__FILE__,__LINE__,'TEMP CHECK ERROR')
      end if
#endif

      uvdrag = sum(drag,1)*rdnnsg
      hfx = sum(sent,1)*rdnnsg
      qfx = sum(evpr,1)*rdnnsg
      tground2(jci1:jci2,ici1:ici2) = sum(tgrd,1)*rdnnsg
      tground1(jci1:jci2,ici1:ici2) = sum(tgrd,1)*rdnnsg
      if ( ichem == 1 ) then
        ssw2da = sum(ssw,1)*rdnnsg
        sdeltk2d = sum(delt,1)*rdnnsg
        sdelqk2d = sum(delq,1)*rdnnsg
        sfracv2d = sum(sigf,1)*rdnnsg
        sfracb2d = sum((d_one-lncl)*(d_one-scvk),1)*rdnnsg
        sfracs2d = sum(lncl*wt+(d_one-lncl)*scvk,1)*rdnnsg
        svegfrac2d = sum(lncl,1)*rdnnsg
      end if

      where ( ldmsk /= 0 )
        tgbb = sum(((d_one-lncl)*tgrd**4+lncl*tlef**4)**d_rfour,1)*rdnnsg
      else where
        tgbb = sum(tgrd,1)*rdnnsg
      end where

      ! Those are needed for output purposes

      ! Fill accumulators

      if ( ktau > 0 ) then
        if ( ifatm ) then
          if ( associated(atm_tgb_out) ) &
            atm_tgb_out = atm_tgb_out + sum(tgbrd,1)*rdnnsg
          if ( associated(atm_tsw_out) ) &
            atm_tsw_out = atm_tsw_out + sum(tsw,1)*rdnnsg
        end if
        if ( ifsrf ) then
          if ( associated(srf_evp_out) ) &
            srf_evp_out = srf_evp_out + sum(evpr,1)*rdnnsg
          if ( associated(srf_tpr_out) ) &
            srf_tpr_out = srf_tpr_out + totpr
          if ( associated(srf_prcv_out) ) &
            srf_prcv_out = srf_prcv_out + pptc
          if ( associated(srf_zpbl_out) ) &
            srf_zpbl_out = srf_zpbl_out + hpbl
          if ( associated(srf_scv_out) ) &
            srf_scv_out = srf_scv_out + sum(sncv,1)*rdnnsg
          if ( associated(srf_sund_out) ) then
            where( fsw > 120.0D0 )
              srf_sund_out = srf_sund_out + dtbat
            end where
          end if
          if ( associated(srf_runoff_out) ) then
            srf_runoff_out(:,:,1) = srf_runoff_out(:,:,1) + sum(srnof,1)*rdnnsg
            srf_runoff_out(:,:,2) = srf_runoff_out(:,:,2) + sum(trnof,1)*rdnnsg
          end if
          if ( associated(srf_sena_out) ) then
            srf_sena_out = srf_sena_out + sum(sent,1)*rdnnsg
          end if
          if ( associated(srf_flw_out) ) &
            srf_flw_out = srf_flw_out + flw
          if ( associated(srf_fsw_out) ) &
            srf_fsw_out = srf_fsw_out + fsw
          if ( associated(srf_fld_out) ) &
            srf_fld_out = srf_fld_out + flwd
          if ( associated(srf_sina_out) ) &
            srf_sina_out = srf_sina_out + sinc
        end if
        if ( ifsub ) then
          call reorder_add_subgrid(sfcp,sub_ps_out)
          if ( associated(sub_evp_out) ) &
            call reorder_add_subgrid(evpr,sub_evp_out)
          if ( associated(sub_scv_out) ) &
            call reorder_add_subgrid(sncv,sub_scv_out,mask=ldmsk1)
          if ( associated(sub_sena_out) ) &
            call reorder_add_subgrid(sent,sub_sena_out)
          if ( associated(sub_runoff_out) ) then
            call reorder_add_subgrid(srnof,sub_runoff_out,1,ldmsk1)
            call reorder_add_subgrid(trnof,sub_runoff_out,2,ldmsk1)
          end if
        end if
        if ( ifsts ) then
          if ( associated(sts_tgmax_out) ) &
            sts_tgmax_out = max(sts_tgmax_out,sum(tgrd,1)*rdnnsg)
          if ( associated(sts_tgmin_out) ) &
            sts_tgmin_out = min(sts_tgmin_out,sum(tgrd,1)*rdnnsg)
          if ( associated(sts_t2max_out) ) &
            sts_t2max_out(:,:,1) = max(sts_t2max_out(:,:,1),sum(t2m,1)*rdnnsg)
          if ( associated(sts_t2min_out) ) &
            sts_t2min_out(:,:,1) = min(sts_t2min_out(:,:,1),sum(t2m,1)*rdnnsg)
          if ( associated(sts_t2min_out) ) &
            sts_t2avg_out(:,:,1) = sts_t2avg_out(:,:,1) + sum(t2m,1)*rdnnsg
          if ( associated(sts_w10max_out) ) &
            sts_w10max_out(:,:,1) = max(sts_w10max_out(:,:,1), &
              sqrt(sum((u10m**2+v10m**2),1)*rdnnsg))
          if ( associated(sts_pcpmax_out) ) &
            sts_pcpmax_out = max(sts_pcpmax_out,totpr)
          if ( associated(sts_pcpavg_out) ) &
            sts_pcpavg_out = sts_pcpavg_out + totpr
          if ( associated(sts_psmin_out) ) &
            sts_psmin_out = min(sts_psmin_out, &
              (sfps(jci1:jci2,ici1:ici2)+ptop)*d_10)
          if ( associated(sts_sund_out) ) then
            where( fsw > 120.0D0 )
              sts_sund_out = sts_sund_out + dtbat
            end where
          end if
          if ( associated(sts_runoff_out) ) then
            sts_runoff_out(:,:,1) = sts_runoff_out(:,:,1) + sum(srnof,1)*rdnnsg
            sts_runoff_out(:,:,2) = sts_runoff_out(:,:,2) + sum(trnof,1)*rdnnsg
          end if
        end if
        if ( iflak ) then
          if ( associated(lak_tpr_out) ) &
            lak_tpr_out = lak_tpr_out + totpr
          if ( associated(lak_scv_out) ) &
            lak_scv_out = lak_scv_out + sum(sncv,1)*rdnnsg
          if ( associated(lak_sena_out) ) &
            lak_sena_out = lak_sena_out + sum(sent,1)*rdnnsg
          if ( associated(lak_flw_out) ) &
            lak_flw_out = lak_flw_out + flw
          if ( associated(lak_fsw_out) ) &
            lak_fsw_out = lak_fsw_out + fsw
          if ( associated(lak_fld_out) ) &
            lak_fld_out = lak_fld_out + flwd
          if ( associated(lak_sina_out) ) &
            lak_sina_out = lak_sina_out + sinc
          if ( associated(lak_evp_out) ) &
            lak_evp_out = lak_evp_out + sum(evpr,1)*rdnnsg
          if ( associated(lak_aveice_out) ) then
            do i = ici1 , ici2
              do j = jci1 , jci2
                do n = 1 , nnsg
                  if ( aveice(n,j,i) < dmissval ) then
                    lak_aveice_out(j,i) = lak_aveice_out(j,i) + &
                      aveice(n,j,i)*rdnnsg*d_r1000
                  else
                    lak_aveice_out(j,i) = dmissval
                  end if
                end do
              end do
            end do
          end if
        end if
      end if

      ! Those are for the output, but collected only at POINT in time

      if ( mod(ktau+1,ksrf) == 0 ) then

        if ( ifsrf ) then
          if ( associated(srf_uvdrag_out) ) &
            srf_uvdrag_out = uvdrag
          if ( associated(srf_tg_out) ) &
            srf_tg_out = tground1(jci1:jci2,ici1:ici2)
          if ( associated(srf_tlef_out) ) then
            where ( ldmsk > 0 )
              srf_tlef_out = sum(tlef,1)*rdnnsg
            elsewhere
              srf_tlef_out = dmissval
            end where
          end if
          if ( associated(srf_aldirs_out) ) &
            srf_aldirs_out = aldirs
          if ( associated(srf_aldifs_out) ) &
            srf_aldifs_out = aldifs
          if ( associated(srf_seaice_out) ) &
            srf_seaice_out = sum(sfice,1)*rdnnsg*d_r1000
          if ( associated(srf_t2m_out) ) &
            srf_t2m_out(:,:,1) = sum(t2m,1)*rdnnsg
          if ( associated(srf_q2m_out) ) &
            srf_q2m_out(:,:,1) = sum(q2m,1)*rdnnsg
          if ( associated(srf_u10m_out) ) &
            srf_u10m_out(:,:,1) = sum(u10m,1)*rdnnsg
          if ( associated(srf_v10m_out) ) &
            srf_v10m_out(:,:,1) = sum(v10m,1)*rdnnsg
          if ( associated(srf_smw_out) ) then
            srf_smw_out(:,:,1) = sum(ssw,1)*rdnnsg
            srf_smw_out(:,:,2) = sum(rsw,1)*rdnnsg
          end if
        end if

        if ( ifsub ) then
          if ( associated(sub_uvdrag_out) ) &
            call reorder_subgrid(drag,sub_uvdrag_out)
          if ( associated(sub_tg_out) ) &
            call reorder_subgrid(tgrd,sub_tg_out)
          if ( associated(sub_tlef_out) ) &
            call reorder_subgrid(tlef,sub_tlef_out,mask=ldmsk1)
          if ( llake ) then
            if ( associated(sub_tlake_out) ) &
            call reorder_subgrid(tlak,sub_tlake_out,1)
            sub_tlake_out = sub_tlake_out + tzero
          end if
          if ( associated(sub_u10m_out) ) &
            call reorder_subgrid(u10m,sub_u10m_out)
          if ( associated(sub_v10m_out) ) &
            call reorder_subgrid(v10m,sub_v10m_out)
          if ( associated(sub_t2m_out) ) &
            call reorder_subgrid(t2m,sub_t2m_out)
          if ( associated(sub_q2m_out) ) &
            call reorder_subgrid(q2m,sub_q2m_out)
          if ( associated(sub_smw_out) ) then
            call reorder_subgrid(ssw,sub_smw_out,1,ldmsk1)
            call reorder_subgrid(rsw,sub_smw_out,2,ldmsk1)
          end if
        end if

        if ( iflak ) then
          if ( associated(lak_tg_out) ) &
            lak_tg_out = tground1(jci1:jci2,ici1:ici2)
          if ( associated(lak_aldirs_out) ) &
            lak_aldirs_out = aldirs
          if ( associated(lak_aldifs_out) ) &
            lak_aldifs_out = aldifs
          if ( associated(lak_hsnow_out) ) &
            lak_hsnow_out = sum(hsnow,1)*rdnnsg
          if ( associated(lak_tlake_out) ) &
            lak_tlake_out = sum(tlak,1)*rdnnsg+tzero
        end if

      end if ! IF output time

      ! Reset accumulation from precip and cumulus
      pptnc = d_zero
      pptc  = d_zero

    end if ! Versus of the interface (1,2)

#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
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
  subroutine albedobats
    implicit none
!
    real(rk8) :: age , albg , albgl , albgld , albgs , albgsd , albl ,  &
               albld , albs , albsd , albzn , alwet , cf1 , cff ,     &
               conn , cons , czeta , czf , dfalbl , dfalbs , dralbl , &
               dralbs , sfac , sl , sl2 , sli , tdiff , tdiffs , wet
    real(rk8) , dimension(nnsg) :: albvl_s , albvs_s , aldifl_s ,       &
                                 aldifs_s , aldirl_s , aldirs_s
    integer(ik4) :: kolour , n , i , j
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'albedobats'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    !
    ! =================================================================
    ! 1. set initial parameters
    ! =================================================================
    !
    do i = ici1 , ici2
      do j = jci1 , jci2
        do n = 1 , nnsg
          lveg(n,j,i) = iveg1(n,j,i)
          if ( ldmsk1(n,j,i) == 2 ) then
            lveg(n,j,i) = 12
            if ( iemiss == 1 ) emiss(n,j,i) = 0.97D0
          else if ( ldmsk1(n,j,i) == 0 ) then
            if ( iemiss == 1 ) emiss(n,j,i) = 0.995D0
          end if
        end do
      end do
    end do
    !
    ! Desert seasonal albedo
    ! Works for Sahara desert and generally northern emisphere
    ! In souther emisphere only some points have this class
    !
    if ( ldesseas ) then
      if ( xmonth == 1 .or. xmonth == 2 .or. xmonth == 12 ) then
        solour(1) = 0.12D0
      endif        
      if ( xmonth == 3 .or. xmonth == 4 .or. xmonth == 5 ) then
        solour(1) = 0.15D0
      endif        
      if ( xmonth == 6 .or. xmonth == 7 .or. xmonth == 8) then
        solour(1) = 0.18D0
      endif        
      if ( xmonth == 9 .or. xmonth == 10 .or. xmonth == 11) then
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
    call fseas(tgbrd)
    do i = ici1 , ici2
      do j = jci1 , jci2
        czeta = coszrs(j,i)
        do n = 1 , nnsg
          lncl(n,j,i) = mfcv(lveg(n,j,i)) - seasf(lveg(n,j,i))*aseas(n,j,i)
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
            sfac = d_one - aseas(n,j,i)
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
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
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
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'soilbc'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif 
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
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine soilbc
!
end module mod_bats_mtrxbats
