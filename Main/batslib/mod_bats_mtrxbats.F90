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
              iemiss , ksrf , xmonth , rtsrf , ktau , kday
  use mod_mppparam
  use mod_mpmessage
  use mod_constants
  use mod_service
  use mod_bats_param
  use mod_bats_common
  use mod_bats_internal
  use mod_bats_lake
  use mod_bats_bndry
  use mod_bats_drag
  use mod_bats_leaftemp
  use mod_bats_zengocn
  use mod_bats_albedo
  use mod_bats_coare
  use mod_outvars
  use mod_regcm_types

  private

  public :: interf , initbats , mtrxbats , albedobats

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
!   initbats
!   mtrxbats ==> soilbc
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
  subroutine mtrxbats
    implicit none
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'mtrxbats'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
!
!---------------------------------------------------------------------
!
!   Excange from model to BATS

    call interf(1)

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
!    if ( iocncpl == 1 .and. iocnflx == 3 ) then
!      call zengocndrv
!      call coare3_drv
!    end if

!   Accumulate quantities for energy and moisture budgets

    call interf(2)
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
  subroutine initbats
    implicit none
    integer(ik4) :: i , itex , j , n , ib
    logical , parameter :: snowhack = .false.
!
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'initbats'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    if ( ktau == 0 ) then
      ib = 1
      do i = ici1 , ici2
        do j = jci1 , jci2
          do n = 1 , nnsg
            lveg(ib) = iveg1(n,j,i)
            mask(ib) = ldmsk1(n,j,i)
            lat(ib) = xlat1(n,j,i)
            tgrd(ib) = tground2(j,i)
            tgbrd(ib) = tground2(j,i)
            taf(ib) = tground2(j,i)
            tlef(ib) = tground2(j,i)
            itex  = iexsol(lveg(ib))
            if ( mask(ib) == 2 ) then
              lveg(ib) = 12
            end if
            if ( mask(ib) > 0 ) then
              if ( snowam(j,i) > d_zero .and. snowam(j,i) < dmissval ) then
                sncv(ib) = snowam(j,i)
              else
                if ( snowhack ) then
                  if ( ht1(n,j,i) > 2000.0D0 ) then
                    sncv(ib) = 0.01D0
                  else
                    sncv(ib) = d_zero
                  end if
                else
                  sncv(ib) = d_zero
                end if
              end if
            end if
            ! Initialize soil moisture in the 3 layers
            tsw(ib) = deptv(lveg(ib))*xmopor(itex)*slmo(lveg(ib))
            rsw(ib) = deprv(lveg(ib))*xmopor(itex)*slmo(lveg(ib))
            ssw(ib) = depuv(lveg(ib))*xmopor(itex)*slmo(lveg(ib))
            gwet(ib) = d_half
            ib = ib + 1
          end do
        end do
      end do
      !
      ! Calculate emission coefficients
      !
      if ( iemiss == 1 ) then
        do i = ilndbeg , ilndend
          if ( lveg(i) == 14 .or. lveg(i) == 15 ) then
            emiss(i) = 0.955D0
          else if ( lveg(i) == 8 ) then
            emiss(i) = 0.76D0
          else if ( lveg(i) == 11 ) then
            emiss(i) = 0.85D0
          else if ( lveg(i) == 12 ) then
            emiss(i) = 0.97D0
          else
            emiss(i) = 0.99D0-(albvgs(lveg(i)) + albvgl(lveg(i)))*0.1D0
          end if
        end do
      else
        emiss(:) = 1.0D0
      end if
      !
      ! Export emissivity
      !
      ib = 1
      do i = ici1, ici2
        do j = jci1, jci2
          do n = 1 , nnsg
            emiss1(n,j,i) = emiss(ib)
            ib = ib + 1
          end do
        end do
      end do
      emissivity = sum(emiss1,1)*rdnnsg
    else
      ib = 1
      do i = ici1 , ici2
        do j = jci1 , jci2
          do n = 1 , nnsg
            lat(ib) = xlat1(n,j,i)
            tgrd(ib) = tgrd1(n,j,i)
            tgbrd(ib) = tgbrd1(n,j,i)
            sncv(ib) = sncv1(n,j,i)
            gwet(ib) = gwet1(n,j,i)
            snag(ib) = snag1(n,j,i)
            sfice(ib) = sfice1(n,j,i)
            ldew(ib) = ldew1(n,j,i)
            taf(ib) = taf1(n,j,i)
            emiss(ib) = emiss1(n,j,i)
            mask(ib) = ldmsk1(n,j,i)
            tlef(ib) = tlef1(n,j,i)
            ssw(ib) = ssw1(n,j,i)
            rsw(ib) = rsw1(n,j,i)
            tsw(ib) = tsw1(n,j,i)
            ib = ib + 1
          end do
        end do
      end do
    end if
    !
    ! initialize hostetler lake model
    !
    if ( llake ) then
      ib = 1
      do i = ici1 , ici2
        do j = jci1 , jci2
          do n = 1 , nnsg
            dhlake(ib) = dhlake1(n,j,i)
            ib = ib + 1
          end do
        end do
      end do
      call initlake
    end if
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine initbats
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
  subroutine interf(ivers)
    implicit none
    integer(ik4) , intent (in) :: ivers
!
    real(rk8) :: facb , facs , fact , factuv , facv , fracb ,  &
                 fracs , fracv , hl , rh0 , satvp , solvt , p0 , qs0 , ts0
    integer(ik4) :: i , j , n , icemsk , ib , ii
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'interf'
    integer(ik4) , save :: idindx = 0
    integer(ik4) :: ierr
    call time_begin(subroutine_name,idindx)
#endif
 
    if ( ivers == 1 ) then ! regcm --> bats

      call fseas(tgbrd)

      ib = 1
      do i = ici1 , ici2
        do j = jci1 , jci2
          p0 = (sfps(j,i)+ptop)*d_1000
          qs0 = qvatm(j,i)/(d_one+qvatm(j,i))
          ts0 = thatm(j,i)
          hl = lh0 - lh1*(ts0-tzero)
          satvp = lsvp1*dexp(lsvp2*hl*(d_one/tzero-d_one/ts0))
          rh0 = dmax1(qs0/(ep2*satvp/(p0*0.01D0-satvp)),d_zero)
          solvt = swdif(j,i) + swdir(j,i)
          do n = 1 , nnsg
            mask(ib) = ldmsk1(n,j,i)
            lveg(ib) = iveg1(n,j,i)
            if ( mask(ib) == 2 ) then
              lveg(ib) = 12
            end if
            qs(ib) = qs0
            sts(ib) = ts0-lrate*regrav*(ht1(n,j,i)-ht(j,i))
            sfcp(ib) = p0*(sts(ib)/ts0)
            hl = lh0 - lh1*(sts(ib)-tzero)
            satvp = lsvp1*dexp(lsvp2*hl*(d_one/tzero-d_one/sts(ib)))
            qs(ib) = dmax1(rh0*ep2*satvp/(sfcp(ib)*0.01D0-satvp),d_zero)
            rhs(ib) = sfcp(ib)/(rgas*sts(ib))
            czenith(ib) = zencos(j,i)
            ! Average over the priod
            prcp(ib) = (ncprate(j,i) + cprate(j,i))*rtsrf
            !
            ! quantities stored on 2d surface array for bats use only
            !
            sent(ib) = hfx(j,i)
            evpr(ib) = qfx(j,i)
            lncl(ib) = mfcv(lveg(ib)) - seasf(lveg(ib))*aseas(ib)
            zh(ib) = hgt(j,i)
            z1log(ib)  = dlog(zh(ib))
            z2fra(ib)  = dlog(zh(ib)*d_half)
            z10fra(ib) = dlog(zh(ib)*d_r10)
            zlgocn(ib) = dlog(zh(ib)/zoce)
            zlglnd(ib) = dlog(zh(ib)/zlnd)
            zlgsno(ib) = dlog(zh(ib)/zsno)
            zlgveg(ib) = dlog(zh(ib)/rough(lveg(ib)))
            zlgdis(ib) = dlog(zh(ib)-displa(lveg(ib))/rough(lveg(ib)))
            if ( solvt > d_zero ) then
              fracd(ib) = swdif(j,i)/solvt
            else
              fracd(ib) = 0.2D0
            end if
            usw(ib) = uatm(j,i)
            vsw(ib) = vatm(j,i)
            swsi(ib) = solar(j,i)
            swflx(ib) = rswf(j,i)
            lwflx(ib) = rlwf(j,i)
            abswveg(ib) = vegswab(j,i)
            ib = ib + 1
          end do
        end do
      end do
 
      ib = 1
      do i = ici1 , ici2
        do j = jci1 , jci2
          rh0 = d_zero
          qs0 = qvatm(j,i)/(d_one+qvatm(j,i))
          ii = 0
          do n = 1 , nnsg
            rh0 = rh0 + (qs(ib+ii)-qs0)
            ii = ii + 1
          end do
          rh0 = rh0/nnsg
          ii = 0
          do n = 1 , nnsg
            qs(ib+ii) = dmax1(qs(ib+ii)-rh0,d_zero)
            ii = ii + 1
          end do
          ib = ib+nnsg 
        end do
      end do

      if ( iocnflx == 1 ) then
        ib = 1
        do i = ici1 , ici2
          do j = jci1 , jci2
            do n = 1 , nnsg
              if ( mask(ib) == 0 ) then
                tgrd(ib) = tground2(j,i)
              end if
              ib = ib+1
            end do
          end do
        end do
      end if

    else if ( ivers == 2 ) then ! bats --> regcm2d
 
      call bats1dto3d

      ! Re-create the land sea mask to account for ice sheet melting

      if ( lseaice .or. llake ) then
        ib = 1
        do i = ici1 , ici2
          do j = jci1 , jci2
            do n = 1 , nnsg
              ldmsk1(n,j,i) = mask(ib)
              ib = ib + 1
            end do
          end do
        end do
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

      ib = 1
      do i = ici1 , ici2
        do j = jci1 , jci2
          do n = 1 , nnsg
            if ( mask(ib) /= 0 ) then
              fracv = sigf(ib)
              fracb = (d_one-lncl(ib))*(d_one-scvk(ib))
              fracs = lncl(ib)*wt(ib) + (d_one-lncl(ib))*scvk(ib)
              facv = z2fra(ib)/zlgveg(ib)
              facb = z2fra(ib)/zlglnd(ib)
              facs = z2fra(ib)/zlgsno(ib)
              fact = fracv*facv + fracb*facb + fracs*facs
              facv = z10fra(ib)/zlgveg(ib)
              facb = z10fra(ib)/zlglnd(ib)
              facs = z10fra(ib)/zlgsno(ib)
              factuv = fracv*facv + fracb*facb + fracs*facs
              u10m(n,j,i) = usw(ib)*(d_one-factuv)
              v10m(n,j,i) = vsw(ib)*(d_one-factuv)
              t2m(n,j,i) = sts(ib) - delt(ib)*fact
              q2m(n,j,i) = qs(ib) - delq(ib)*fact
            else
              if ( iocnflx == 1 ) then
                fact = z2fra(ib)/zlgocn(ib)
                factuv = z10fra(ib)/zlgocn(ib)
                u10m(n,j,i) = usw(ib)*(d_one-factuv)
                v10m(n,j,i) = vsw(ib)*(d_one-factuv)
                t2m(n,j,i) = sts(ib) - delt(ib)*fact
                q2m(n,j,i) = qs(ib) - delq(ib)*fact
              end if
            end if
            ib = ib + 1
          end do
        end do
      end do
 
#ifdef DEBUG
      ierr = 0
      do i = ici1 , ici2
        do j = jci1 , jci2
          do n = 1 , nnsg
            if ( tgrd1(n,j,i) < 150.0D0 ) then
              write(stderr,*) 'Likely error: Surface temperature too low'
              write(stderr,*) 'J   = ',global_dot_jstart+j
              write(stderr,*) 'I   = ',global_dot_istart+i
              write(stderr,*) 'VAL = ',tgrd1(n,j,i)
              ierr = ierr + 1
            end if
          end do
        end do
      end do
      if ( ierr /= 0 ) then
        call fatal(__FILE__,__LINE__,'TEMP CHECK ERROR')
      end if
#endif
      uvdrag = sum(drag1,1)*rdnnsg
      hfx = sum(sent1,1)*rdnnsg
      qfx = sum(evpr1,1)*rdnnsg
      tground2 = sum(tgrd1,1)*rdnnsg
      tground1 = sum(tgrd1,1)*rdnnsg
      if ( ichem == 1 ) then
        ssw2da = sum(ssw1,1)*rdnnsg
        deltat = sum(delt,1)*rdnnsg
        deltaq = sum(delq,1)*rdnnsg
        sfracv2d = sum(sigf,1)*rdnnsg
        sfracb2d = sum((d_one-lncl)*(d_one-scvk),1)*rdnnsg
        sfracs2d = sum(lncl*wt+(d_one-lncl)*scvk,1)*rdnnsg
        svegfrac2d = sum(lncl,1)*rdnnsg
      end if

      ib = 1
      do i = ici1 , ici2
        do j = jci1 , jci2
          if ( ldmsk(j,i) == 0 ) then
            tgbb(j,i) = sum(tgrd(ib:ib+nnsg-1))*rdnnsg
          else
            tgbb(j,i) = sum(                                  &
             ((d_one-lncl(ib:ib+nnsg-1))*tgrd(ib:ib+nnsg-1)**4+ &
                    (lncl(ib:ib+nnsg-1))*tlef(ib:ib+nnsg-1)**4)**d_rfour)*rdnnsg
          end if
          ib = ib + nnsg
        end do
      end do

      ! Those are needed for output purposes

      ! Fill accumulators

      if ( ktau > 0 ) then
        if ( ifatm ) then
          if ( associated(atm_tgb_out) ) &
            atm_tgb_out = atm_tgb_out + sum(tgbrd1,1)*rdnnsg
          if ( associated(atm_tsw_out) ) &
            atm_tsw_out = atm_tsw_out + sum(tsw1,1)*rdnnsg
        end if
        if ( ifsrf ) then
          if ( associated(srf_evp_out) ) &
            srf_evp_out = srf_evp_out + sum(evpr1,1)*rdnnsg
          if ( associated(srf_tpr_out) ) &
            srf_tpr_out = srf_tpr_out + prcp1(1,:,:)
          if ( associated(srf_prcv_out) ) &
            srf_prcv_out = srf_prcv_out + cprate
          if ( associated(srf_zpbl_out) ) &
            srf_zpbl_out = srf_zpbl_out + hpbl
          if ( associated(srf_scv_out) ) &
            srf_scv_out = srf_scv_out + sum(sncv1,1)*rdnnsg
          if ( associated(srf_sund_out) ) then
            where( rswf > 120.0D0 )
              srf_sund_out = srf_sund_out + dtbat
            end where
          end if
          if ( associated(srf_runoff_out) ) then
            srf_runoff_out(:,:,1) = srf_runoff_out(:,:,1) + sum(srnof1,1)*rdnnsg
            srf_runoff_out(:,:,2) = srf_runoff_out(:,:,2) + sum(trnof1,1)*rdnnsg
          end if
          if ( associated(srf_sena_out) ) then
            srf_sena_out = srf_sena_out + sum(sent1,1)*rdnnsg
          end if
          if ( associated(srf_flw_out) ) &
            srf_flw_out = srf_flw_out + rlwf
          if ( associated(srf_fsw_out) ) &
            srf_fsw_out = srf_fsw_out + rswf
          if ( associated(srf_fld_out) ) &
            srf_fld_out = srf_fld_out + dwrlwf
          if ( associated(srf_sina_out) ) &
            srf_sina_out = srf_sina_out + solinc
          if ( associated(srf_snowmelt_out) ) &
            srf_snowmelt_out = srf_snowmelt_out + sum(sm,1)*rdnnsg
        end if
        if ( ifsub ) then
          call reorder_add_subgrid(ps1,sub_ps_out)
          if ( associated(sub_evp_out) ) &
            call reorder_add_subgrid(evpr1,sub_evp_out)
          if ( associated(sub_scv_out) ) &
            call reorder_add_subgrid(sncv1,sub_scv_out,mask=ldmsk1)
          if ( associated(sub_sena_out) ) &
            call reorder_add_subgrid(sent1,sub_sena_out)
          if ( associated(sub_runoff_out) ) then
            call reorder_add_subgrid(srnof1,sub_runoff_out,1,ldmsk1)
            call reorder_add_subgrid(trnof1,sub_runoff_out,2,ldmsk1)
          end if
        end if
        if ( ifsts ) then
          if ( associated(sts_tgmax_out) ) &
            sts_tgmax_out = max(sts_tgmax_out,sum(tgrd1,1)*rdnnsg)
          if ( associated(sts_tgmin_out) ) &
            sts_tgmin_out = min(sts_tgmin_out,sum(tgrd1,1)*rdnnsg)
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
            sts_pcpmax_out = max(sts_pcpmax_out,prcp1(1,:,:))
          if ( associated(sts_pcpavg_out) ) &
            sts_pcpavg_out = sts_pcpavg_out + prcp1(1,:,:)
          if ( associated(sts_psmin_out) ) &
            sts_psmin_out = min(sts_psmin_out, &
              (sfps(jci1:jci2,ici1:ici2)+ptop)*d_10)
          if ( associated(sts_sund_out) ) then
            where( rswf > 120.0D0 )
              sts_sund_out = sts_sund_out + dtbat
            end where
          end if
          if ( associated(sts_runoff_out) ) then
            sts_runoff_out(:,:,1) = sts_runoff_out(:,:,1) + sum(srnof1,1)*rdnnsg
            sts_runoff_out(:,:,2) = sts_runoff_out(:,:,2) + sum(trnof1,1)*rdnnsg
          end if
        end if
        if ( iflak ) then
          if ( associated(lak_tpr_out) ) &
            lak_tpr_out = lak_tpr_out + prcp1(1,:,:)
          if ( associated(lak_scv_out) ) &
            lak_scv_out = lak_scv_out + sum(sncv1,1)*rdnnsg
          if ( associated(lak_sena_out) ) &
            lak_sena_out = lak_sena_out + sum(sent1,1)*rdnnsg
          if ( associated(lak_flw_out) ) &
            lak_flw_out = lak_flw_out + rlwf
          if ( associated(lak_fsw_out) ) &
            lak_fsw_out = lak_fsw_out + rswf
          if ( associated(lak_fld_out) ) &
            lak_fld_out = lak_fld_out + dwrlwf
          if ( associated(lak_sina_out) ) &
            lak_sina_out = lak_sina_out + solinc
          if ( associated(lak_evp_out) ) &
            lak_evp_out = lak_evp_out + sum(evpr1,1)*rdnnsg
          if ( associated(lak_aveice_out) ) then
            lak_aveice_out = lak_aveice_out + &
              sum(sfice1*lakmsk1,1)*rdnnsg*d_r1000
          end if
        end if
      end if

      ! Those are for the output, but collected only at POINT in time

      if ( mod(ktau+1,ksrf) == 0 ) then

        if ( ifsrf ) then
          if ( associated(srf_uvdrag_out) ) &
            srf_uvdrag_out = sum(drag1,1)*rdnnsg
          if ( associated(srf_tg_out) ) &
            srf_tg_out = tground1
          if ( associated(srf_tlef_out) ) then
            where ( ldmsk > 0 )
              srf_tlef_out = sum(tlef1,1)*rdnnsg
            elsewhere
              srf_tlef_out = dmissval
            end where
          end if
          if ( associated(srf_aldirs_out) ) &
            srf_aldirs_out = swdiralb
          if ( associated(srf_aldifs_out) ) &
            srf_aldifs_out = swdifalb
          if ( associated(srf_seaice_out) ) &
            srf_seaice_out = sum(sfice1,1)*rdnnsg*d_r1000
          if ( associated(srf_t2m_out) ) &
            srf_t2m_out(:,:,1) = sum(t2m,1)*rdnnsg
          if ( associated(srf_q2m_out) ) &
            srf_q2m_out(:,:,1) = sum(q2m,1)*rdnnsg
          if ( associated(srf_u10m_out) ) &
            srf_u10m_out(:,:,1) = sum(u10m,1)*rdnnsg
          if ( associated(srf_v10m_out) ) &
            srf_v10m_out(:,:,1) = sum(v10m,1)*rdnnsg
          if ( associated(srf_smw_out) ) then
            srf_smw_out(:,:,1) = sum(ssw1,1)*rdnnsg
            srf_smw_out(:,:,2) = sum(rsw1,1)*rdnnsg
          end if
        end if

        if ( ifsub ) then
          if ( associated(sub_uvdrag_out) ) &
            call reorder_subgrid(drag1,sub_uvdrag_out)
          if ( associated(sub_tg_out) ) &
            call reorder_subgrid(tgrd1,sub_tg_out)
          if ( associated(sub_tlef_out) ) &
            call reorder_subgrid(tlef1,sub_tlef_out,mask=ldmsk1)
          if ( llake ) then
            if ( associated(sub_tlake_out) ) then
              call lake_fillvar(var_tlak,tlake,0,llakmsk1)
              call reorder_subgrid(tlake,sub_tlake_out,1)
              sub_tlake_out = sub_tlake_out + tzero
            end if
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
            call reorder_subgrid(ssw1,sub_smw_out,1,ldmsk1)
            call reorder_subgrid(rsw1,sub_smw_out,2,ldmsk1)
          end if
        end if

        if ( iflak ) then
          if ( associated(lak_tg_out) ) &
            lak_tg_out = tground1
          if ( associated(lak_aldirs_out) ) &
            lak_aldirs_out = swdiralb
          if ( associated(lak_aldifs_out) ) &
            lak_aldifs_out = swdifalb
          if ( associated(lak_hsnow_out) ) &
            lak_hsnow_out = sum(sncv1*lakmsk1,1)*rdnnsg
          if ( associated(lak_tlake_out) ) then
            call lake_fillvar(var_tlak,tlake,0,llakmsk1)
            lak_tlake_out = sum(tlake,1)*rdnnsg+tzero
          end if
        end if

      end if ! IF output time

      if ( iocncpl == 1 ) then
        ! Fill for the RTM component
        dailyrnf(:,:,1) = dailyrnf(:,:,1) + sum(srnof1,1)*rdnnsg
        dailyrnf(:,:,2) = dailyrnf(:,:,2) + (sum(trnof1,1)-sum(srnof1,1))*rdnnsg
        runoffcount = runoffcount + d_one
      end if

      ! Reset accumulation from precip and cumulus
      ncprate = d_zero
      cprate  = d_zero

    end if ! Versus of the interface (1,2)

#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine interf

  subroutine bats1dto3d
    implicit none
    integer(ik4) :: i , j , n , ib
    ib = 1
    do i = ici1 , ici2
      do j = jci1 , jci2
        do n = 1 , nnsg

          tlef1(n,j,i) = tlef(ib)
          tgrd1(n,j,i) = tgrd(ib)
          tgbrd1(n,j,i) = tgbrd(ib)
          gwet1(n,j,i) = gwet(ib)
          ldew1(n,j,i) = ldew(ib)
          snag1(n,j,i) = snag(ib)
          sncv1(n,j,i) = sncv(ib)
          sfice1(n,j,i) = sfice(ib)
          ssw1(n,j,i) = ssw(ib)
          rsw1(n,j,i) = rsw(ib)
          tsw1(n,j,i) = tsw(ib)
          emiss1(n,j,i) = emiss(ib)
          taf1(n,j,i) = taf(ib)

          ps1(n,j,i) = sfcp(ib)
          sent1(n,j,i) = sent(ib)
          evpr1(n,j,i) = evpr(ib)
          prcp1(n,j,i) = prcp(ib)
          trnof1(n,j,i) = trnof(ib)
          srnof1(n,j,i) = srnof(ib)
          drag1(n,j,i) = drag(ib)

          ib = ib + 1
        end do
      end do
    end do
  end subroutine bats1dto3d

  subroutine albedobats
    implicit none
    integer :: i , j , ib
    call interf(1)
    call albedo
    ib = 1
    do i = ici1 , ici2
      do j = jci1 , jci2
        swalb(j,i)  = sum(swal(ib:ib+nnsg-1))*rdnnsg
        lwalb(j,i)  = sum(lwal(ib:ib+nnsg-1))*rdnnsg
        swdiralb(j,i)  = sum(swdiral(ib:ib+nnsg-1))*rdnnsg
        lwdiralb(j,i)  = sum(lwdiral(ib:ib+nnsg-1))*rdnnsg
        swdifalb(j,i)  = sum(swdifal(ib:ib+nnsg-1))*rdnnsg
        lwdifalb(j,i)  = sum(lwdifal(ib:ib+nnsg-1))*rdnnsg
        ib = ib + nnsg
      end do
    end do
  end subroutine albedobats

end module mod_bats_mtrxbats
