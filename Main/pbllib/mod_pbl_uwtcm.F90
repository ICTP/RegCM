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

!**********************************************************************
!     UW PBL SCHEME                                                   *
!     Questions? Contact breth@washington.edu                         *
!                 with questions about PBL physics,                   *
!                 or tobrien@ucsc.edu for technical questions about   *
!                 the coupling of this model to RegCM4                *
!**********************************************************************

! This moist turbulence parameterization is described in detail in 
!   A New Parameterization for Shallow Cumulus Convection and Its Application to
!   Marine Subtropical Cloud-Topped Boundary Layers. Part I:
!   Description and 1D Results
!   CHRISTOPHER S. BRETHERTON, JAMES R. MCCAA, AND HERVE´ GRENIER
!   (2004), Monthly Weather Review
!
!     Order of events:
!     0.  Preliminaries:  Some things we only need to do once
!     1.  First, save initial values and calculate some derivatives.
!     2.  Calculate the buoyancy jump at flux levels
!     3.  Calculate the pbl height and length scale
!     4.  Compute diffusivity profiles
!     5.  Implicit diffusion of thetal, total water
!     6.  Implicit diffusion of momentum
!     7. Prediction of TKE
!     8. Calculation of tendencies for output to model
!
! Major Changes:
!   01/2010, TA O'Brien:
!     * Ported from MM5 code to RegCM3 code base
!   06/2010, TA O'Brien:
!     * Ported to RegCM4.0 code base
!     * Fixed momentum flux coupling
!   09/2010, TA O'Brien:
!     * Changed method for calculating N2 at surface; use
!       surface fluxes
!   01/2011-02/1011, TA O'Brien:
!     * Recoupled to RegCM4.1 code base from scratch
!     * Wrote new tridiagonal solver from scratch (GNU licence purposes)
!     * Split temp/humid solver to mod_thetal.
!     * Fixed non-convergence problem for temp/humid solver
!     * Implemented GB01/Bretherton 2010 length scale (optional)
!   07/2011 Graziano Giuliani
!     * Moved in RegCM core devel

module mod_pbl_uwtcm

  use mod_realkinds
  use mod_memutil
  use mod_dynparam
  use mod_constants
  use mod_mpmessage
  use mod_pbl_common
  use mod_pbl_thetal

  private

  real(dp) , parameter , public :: nuk = 5.0D0 ! multiplier for kethl
  integer , public :: ktmin = 3

  ! Model constants
  ! fraction of turb layer to be considered in bbls
  real(dp) , parameter :: xfr = 0.1D0
  ! see gb01 regarding the next three lines, atwo and rstbl can be tweaked
  real(dp) , parameter :: aone = 1.9D0*xfr
  ! b1/2**(3/2) from Mellor and Yamada (1982)
  real(dp) , parameter :: czero = 5.869D0
  real(dp) , parameter :: rcrit = 0.3D0
  real(dp) , parameter :: etal =  0.085D0

  ! These two parameters can be set in the regcm.in file, so don't set them
  ! here as it will overwrite the values read in from the regcm.in file.
  ! Settable model constants that need default values set
  real(dp) :: rstbl = 1.5D0
  real(dp) :: atwo = 15.0D0
        
  ! Variables that hold frequently-done calculations
  real(dp) :: rcp , rczero

  ! local variables on full levels
  real(dp) , pointer , dimension(:) :: zqx , kth , kzm , rhoxfl , &
                 tke , tkes , rrhoxfl , bbls , nsquar , richnum , &
                 bouyan , rdza , dza , svs , presfl , exnerfl ,   &
                 shear , rexnerfl , rcldb , epop , sm , sh

  ! local variables on half levels
  real(dp) , pointer , dimension(:) :: ux , vx , thx , qx , uthvx ,    &
                 zax , kethl , thlx , thlxs, thxs , tx , tvx ,         &
                 rttenx , preshl , qcx , qwx , qwxs , udzq , rrhoxhl , &
                 uxs , qxs , rhoxhl , exnerhl , rexnerhl , rdzq ,      &
                 vxs , qcxs , aimp , bimp , cimp , uimp1 , rimp1 ,     &
                 uimp2 , rimp2

  integer , pointer , dimension(:) :: isice , ktop , kbot

  ! local scalars
  real(dp) :: uflxp , vflxp , rhoxsf , tgbx , tvcon , fracz , dudz , &
              dvdz , rvls , thv0 , dthv , psbx , templ , temps ,     &
              cell , thgb , pblx , ustxsq , qfxx , hfxx , dth ,      &
              uvdragx , thvflx , kh0 , q0s

  integer :: kpbl2dx  ! Top of PBL
  integer :: kmix2dx  ! Top of mixed layer (decoupled layer)

  integer :: imethod , itbound , ilenparam , iuwvadv

  public :: init_mod_pbl_uwtcm , uwtcm , rstbl , atwo , ilenparam , iuwvadv

  contains

  subroutine init_mod_pbl_uwtcm
    implicit none

    ! If we are the head processor, output the various parameters

    if ( myid == 0 ) then
      write(*,*)'***************************************'
      write(*,*)'******** UW TCM Parameters*************'
      write(*,*)'***************************************'
      write(*,*)'rstbl = ',rstbl
      write(*,*)'atwo = ',atwo
      write(*,*)'iuwvadv = ',iuwvadv
      write(*,*)'ilenparam = ',ilenparam
      write(*,*)'***************************************'
    end if

    ! Variables that hold frequently-done calculations

    rcp = d_one/cpd
    rczero = d_one/czero

    call getmem1d(zqx,1,kzp2,'mod_uwtcm:zqx')

    ! Allocate memory local variables on full levels

    call getmem1d(kth,1,kzp1,'mod_uwtcm:kth')
    call getmem1d(kzm,1,kzp1,'mod_uwtcm:kzm')
    call getmem1d(rhoxfl,1,kzp1,'mod_uwtcm:rhoxfl')
    call getmem1d(rrhoxfl,1,kzp1,'mod_uwtcm:rrhoxfl')
    call getmem1d(tke,1,kzp1,'mod_uwtcm:tke')
    call getmem1d(tkes,1,kzp1,'mod_uwtcm:tkes')
    call getmem1d(bbls,1,kzp1,'mod_uwtcm:bbls')
    call getmem1d(nsquar,1,kzp1,'mod_uwtcm:nsquar')
    call getmem1d(richnum,1,kzp1,'mod_uwtcm:richnum')
    call getmem1d(bouyan,1,kzp1,'mod_uwtcm:bouyan')
    call getmem1d(dza,1,kzp1,'mod_uwtcm:dza')
    call getmem1d(rdza,1,kzp1,'mod_uwtcm:rdza')
    call getmem1d(svs,1,kzp1,'mod_uwtcm:svs')
    call getmem1d(presfl,1,kzp1,'mod_uwtcm:presfl')
    call getmem1d(epop,1,kzp1,'mod_uwtcm:epop')
    call getmem1d(exnerfl,1,kzp1,'mod_uwtcm:exnerfl')
    call getmem1d(rexnerfl,1,kzp1,'mod_uwtcm:rexnerfl')
    call getmem1d(shear,1,kzp1,'mod_uwtcm:shear')
    call getmem1d(rcldb,1,kzp1,'mod_uwtcm:rcldb')
    call getmem1d(sm,1,kzp1,'mod_uwtcm:sm')
    call getmem1d(sh,1,kzp1,'mod_uwtcm:sh')

    ! local variables on half levels
    call getmem1d(ux,1,kz,'mod_uwtcm:ux')
    call getmem1d(vx,1,kz,'mod_uwtcm:vx')
    call getmem1d(qx,1,kz,'mod_uwtcm:qx')
    call getmem1d(thx,1,kz,'mod_uwtcm:thx')
    call getmem1d(uthvx,1,kz,'mod_uwtcm:uthvx')
    call getmem1d(zax,1,kz,'mod_uwtcm:zax')
    call getmem1d(kethl,1,kz,'mod_uwtcm:kethl')
    call getmem1d(thlx,1,kz,'mod_uwtcm:thlx')
    call getmem1d(thlxs,1,kz,'mod_uwtcm:thlxs')
    call getmem1d(thxs,1,kz,'mod_uwtcm:thxs')
    call getmem1d(tx,1,kz,'mod_uwtcm:tx')
    call getmem1d(tvx,1,kz,'mod_uwtcm:tvx')
    call getmem1d(rttenx,1,kz,'mod_uwtcm:rttenx')
    call getmem1d(preshl,1,kz,'mod_uwtcm:preshl')
    call getmem1d(qcx,1,kz,'mod_uwtcm:qcx')
    call getmem1d(qwx,1,kz,'mod_uwtcm:qwx')
    call getmem1d(qwxs,1,kz,'mod_uwtcm:qwxs')
    call getmem1d(udzq,1,kz,'mod_uwtcm:udzq')
    call getmem1d(rrhoxhl,1,kz,'mod_uwtcm:rrhoxhl')
    call getmem1d(uxs,1,kz,'mod_uwtcm:uxs')
    call getmem1d(qxs,1,kz,'mod_uwtcm:qxs')
    call getmem1d(rhoxhl,1,kz,'mod_uwtcm:rhoxhl')
    call getmem1d(exnerhl,1,kz,'mod_uwtcm:exnerhl')
    call getmem1d(rexnerhl,1,kz,'mod_uwtcm:rexnerhl')
    call getmem1d(rdzq,1,kz,'mod_uwtcm:rdzq')
    call getmem1d(vxs,1,kz,'mod_uwtcm:vxs')
    call getmem1d(qcxs,1,kz,'mod_uwtcm:qcxs')
    call getmem1d(aimp,1,kz,'mod_uwtcm:aimp')
    call getmem1d(bimp,1,kz,'mod_uwtcm:bimp')
    call getmem1d(cimp,1,kz,'mod_uwtcm:cimp')
    call getmem1d(uimp1,1,kz,'mod_uwtcm:uimp1')
    call getmem1d(rimp1,1,kz,'mod_uwtcm:rimp1')
    call getmem1d(uimp2,1,kz,'mod_uwtcm:uimp2')
    call getmem1d(rimp2,1,kz,'mod_uwtcm:rimp2')

    call getmem1d(isice,1,kz,'mod_uwtcm:isice')
    call getmem1d(ktop,1,kz,'mod_uwtcm:ktop')
    call getmem1d(kbot,1,kz,'mod_uwtcm:kbot')

  end subroutine init_mod_pbl_uwtcm

  subroutine uwtcm
    implicit none

    integer ::  i , j , k
    integer :: ilay ! layer index
    integer :: iconv , iteration

    !Main do loop
    iloop: &
    do i = itcmstart , itcmend
      jloop: &
      do j = jtcmstart , jtcmend

!*******************************************************************************
!*******************************************************************************
!*********** Initialization of UW TCM for Current Grid Column ******************
!*******************************************************************************
!*******************************************************************************

        ! Copy in local versions of necessary variables
        psbx = sfcps(j,i)
        tgbx = tg(i,j)
        qfxx = qfx(i,j)
        hfxx = hfx(i,j)
        uvdragx = uvdrag(i,j)

        ! Integrate the hydrostatic equation to calculate the level height
        zqx(kzp1) = d_zero
        zqx(kzp1+1) = d_zero
        tke(kzp1) = tkests(i,kzp1,j)

        kinitloop: &
        do k = kz , 1 , -1
          rttenx(k) = radheatrt(j,i,k)
          cell = ptop/psbx
          zqx(k) = zqx(k+1) + rgas/egrav*tatm(j,i,k)*   &
                   log((flev(k+1)+cell)/(flev(k)+cell))
          zax(k) = d_half*(zqx(k)+zqx(k+1))
          tke(k) = tkests(i,k,j)
          tx(k)  = tatm(j,i,k)
          qx(k)  = qvatm(j,i,k)
          qcx(k) = qcatm(j,i,k)
          ux(k)  = uatm(j,i,k)
          vx(k)  = vatm(j,i,k)
          ! if ( tx(k) > tzero ) then
!         if ( tx(k) > tzero ) then
!           isice(k) = 0
!         else
!           isice(k) = 1
!         end if
          isice(k) = 0
        end do kinitloop

        ! Set all the save variables (used for determining the tendencies)
        tkes(kzp1) = tke(kzp1)
        khalfloop: &
        do k = 1 , kz
          ! pressure at half levels
          preshl(k) = hlev(k)*psbx + ptp
          ! Level spacing
          udzq(k) = zqx(k)-zqx(k+1)
          rdzq(k) = d_one/udzq(k)
          ! Exner function
          exnerhl(k)=(preshl(k)/d_100)**rovcp
          rexnerhl(k) = d_one/exnerhl(k)
          ! Potential temperature
          thx(k) = tx(k)*rexnerhl(k)
          ! Total water specific humidity
          qwx(k) = qx(k) + qcx(k)
          qwxs(k) = qwx(k)
          ! Virtual temperature and potential temperature
          tvcon = (d_one + ep1*qx(k)-qcx(k))
          tvx(k) = tx(k)*tvcon ! virtual temperature
          uthvx(k) = thx(k)*tvcon
          ! Liquid water potential temperature (accounting for ice)
          if ( isice(k) == 0 ) then
            thlx(k) = thx(k) - wlhvocp * rexnerhl(k) * qcx(k)
          else
            thlx(k) = thx(k) - wlhsocp * rexnerhl(k) * qcx(k)
          end if
          thlxs(k) = thlx(k)
          ! save initial values of winds, thx, and qx
          tkes(k) = tke(k)
          uxs(k)  = ux(k)
          vxs(k)  = vx(k)
          thxs(k) = thx(k)
          qxs(k) = qx(k)
          qcxs(k) = qcx(k)
          ! density at half levels
          rhoxhl(k)=preshl(k)*d_1000/(rgas*tvx(k))
          rrhoxhl(k) = d_one/rhoxhl(k)
        end do khalfloop

        ! Set variables that are on full levels
        presfl(1) = ptop
        kfullloop: &
        do k = 2 , kz
          ! pressure at full levels
          presfl(k) = flev(k)*psbx + ptp
          epop(k) = ep2/presfl(k)
          ! Level spacing
          dza(k) = zax(k-1)-zax(k)
          rdza(k) = d_one/dza(k)
          ! Exner function
          exnerfl(k)=(presfl(k)/d_100)**rovcp
          rexnerfl(k) = d_one/exnerfl(k)
          ! Density
          fracz = (zqx(k)-zax(k))*rdza(k)
          rhoxfl(k) = rhoxhl(k)+(rhoxhl(k-1)-rhoxhl(k))*fracz
          rrhoxfl(k) = d_one/rhoxfl(k)
          ! Wind gradient and shear
          dudz =  (ux(k-1)  - ux(k))  *rdza(k)
          dvdz =  (vx(k-1)  - vx(k))  *rdza(k)
          svs(k) = dudz*dudz + dvdz*dvdz
        end do kfullloop

        ! pressure at the surface (in centibars)
        presfl(kzp1) = psbx + ptp
        ! Surface exner function
        rexnerfl(kzp1) = (d_100/presfl(kzp1))**rovcp
        exnerfl(kzp1) = d_one/rexnerfl(kzp1)
        ! density at the surface
        rhoxsf = (presfl(kzp1)*d_1000)/(rgas*tvx(kz))
           
!*******************************************************************************
!*******************************************************************************
!*********** Calculation (and Conversion) of Surface Fluxes ********************
!*******************************************************************************
!*******************************************************************************

        ! more surface variables
        thgb = tgbx * rexnerfl(kzp1)
        ! Calculate the saturation specific humidity just above the surface
        ! Assume for now that the moisture will be available.
        ! TODO: This is a bad assumption and needs to be dealt with.
        ! This assumption will tend to cause an over-estimation of the surface
        ! virtual heat flux over land, which should cause an overestimation
        ! of boundary layer height over land.
        q0s = ep2/(presfl(kzp1)/esatw(presfl(kzp1),tgbx) - d_one)
        ! Calculate the virtual temperature right above the surface
        thv0 = thgb*(d_one+ep1*q0s)
        ! Calculate the change in virtual potential temperature from
        ! the surface to the first interface
        dthv = uthvx(kz)-thv0
        ! Calculate the change in potential temperature from the surface
        ! to the first interface
        dth = thx(kz) - thgb 
            
        ! Calculate surface momentum fluxes
        uflxp = -uvdragx*ux(kz)/rhoxsf
        vflxp = -uvdragx*vx(kz)/rhoxsf
        ustxsq = dsqrt(uflxp*uflxp + vflxp*vflxp)

        ! Estimate of the surface virtual heat flux
        thvflx = hfxx/rhoxsf*rcp*(d_one+ep1*q0s) + ep1/thgb*qfxx*rhoxsf
        ! Estimate of surface eddy diffusivity, for estimating the
        ! surface N^2 from the surface virtual heat flux
        kh0 = vonkar*d_one*dsqrt(dmax1(tkemin,3.25D0 * ustxsq))


!*******************************************************************************
!*******************************************************************************
!********************* Calculation of boundary layer height ********************
!*******************************************************************************
!*******************************************************************************
          
        ! Calculate nsquared Set N^2 based on the current potential
        ! temperature profile
        call n2(thlx,qwx,kz)
        ! nsquar(kzp1) = egrav/uthvx(kz) * dthv / zax(kz)
        ! Estimate the surface N^2 from the surface virtual heat flux
        nsquar(kzp1) = -egrav/thgb*thvflx/kh0

        ! Calculate the bulk richardson number
        richnum = nsquar/dmax1(svs,1.0D-8)

        ! Calculate the boundary layer height
        call pblhgt(kzp1,iconv)
        ! call pblhgt_tao(kzp1,iconv)

!*******************************************************************************
!*******************************************************************************
!************************* Semi-implicit Integration ***************************
!*******************************************************************************
!*******************************************************************************

        ! Calculate the diffusivity profiles, update total water and
        ! liquid water potential temperature, then re-calculate the
        ! diffusivity profiles using the updated values, and re-integrate.
        ! Also update N^2 along the way.
        myiteration: &
        do iteration = 0 , 1
          !*************************************************************
          !***** Semi-implicit calculation of diffusivity profiles *****
          !*************************************************************
          call my(kzp1,iconv)

          !*************************************************************
          !****** Implicit Diffusion of Thetal and Qtot ****************
          !*************************************************************
          ! first find the coefficients that apply to 
          ! all scalars at half-levels
          diffqandthlx: &
          do k = 1 , kz
            if ( k == 1 ) then
              aimp(k) = d_zero
            else
              aimp(k) = -(rhoxfl(k)*rrhoxhl(k))* &
                          kth(k) * dtpbl *rdzq(k)*rdza(k)
            end if
            if ( k == kz ) then
              cimp(k) = d_zero
            else
              cimp(k) = -(rhoxfl(k+1)*rrhoxhl(k))* &
                          kth(k+1) * dtpbl *rdzq(k)*rdza(k+1)
            end if
            bimp(k) = d_one - aimp(k) - cimp(k)
            ! now find right side for various scalars:
            ! no flux out top, so no (k == 1)
            if ( k == kz ) then ! at surface
              ! include surface latent heat flux
              rimp2(k) = qwx(k) + dtpbl * qfxx*rrhoxhl(k)*rdzq(kz)
              ! include surface sensible heat flux
              rimp1(k) = thlx(k) + &
                     dtpbl * hfxx*rrhoxhl(k)*rcp*rdzq(kz)*rexnerhl(kz)
            else
              rimp2(k) = qwx(k)
              rimp1(k) = thlx(k)
            end if
          end do diffqandthlx
          ! Solve total water
          call solve_tridiag(aimp,bimp,cimp,rimp2,uimp2,kz)
          ! Solve liquid water potential temperature
          call solve_tridiag(aimp,bimp,cimp,rimp1,uimp1,kz)

          ! Calculate nsquared Set N^2 based on the updated potential
          ! temperature profile (this is for the semi-implicit integration)
          call n2(uimp1,uimp2,kz)
        end do myiteration

           
        !*************************************************************
        !************ Re-calculate thx, qx, and qcx ******************
        !*************************************************************

        recalcscalar: &
        do k = 1 , kz
          ! Set thlx and qwx to their updated values
          thlx(k) = uimp1(k)
          qwx(k) = uimp2(k)
             
          ! Determine the temperature and qc and qv
          ! imethod = 1     ! Use the Brent 1973 method to solve for T
          imethod = 2       ! Use Bretherton's iterative method to solve for T
          ! imethod = 3     ! Use a finite difference method to solve for T
          itbound = 100 ! Do a maximum of 100 iterations (usually are no more
                        ! than about 30)
          temps = solve_for_t(thlx(k),qwx(k),preshl(k),thxs(k)*exnerhl(k), &
                              qwxs(k),qcxs(k),thlxs(k),itbound,imethod,    &
                              isice(k),qcx(k),qx(k))

          if ( itbound == -999 ) then
            call fatal(__FILE__,__LINE__,'UW PBL SCHEME ALGO ERROR')
          end if 

          thx(k) = temps*rexnerhl(k)
          uthvx(k)=thx(k)*(d_one + ep1*qx(k)-qcx(k))
        end do recalcscalar

        !*************************************************************
        !****** Implicit Diffusion of U and V ************************
        !*************************************************************
        diffuv: &
        do k = 1 , kz
          if ( k == 1 ) then
            aimp(k) = d_zero
          else
            aimp(k) = -(rhoxfl(k)*rrhoxhl(k))*     &
                        kzm(k) * dtpbl *rdzq(k)*rdza(k)
          end if
          if ( k == kz ) then
            cimp(k) = d_zero
          else
            cimp(k) = -(rhoxfl(k+1)*rrhoxhl(k))*   &
                        kzm(k+1) * dtpbl *rdzq(k)*rdza(k+1)
          end if
          bimp(k) = d_one - aimp(k) - cimp(k)
          ! now find right side 
          ! no flux out top, so no (k == 1)
          if ( k == kz ) then
            ! at surface
            ! include surface momentum fluxes
            rimp1(k) = ux(k) + dtpbl *              &
                       uflxp * (rhoxsf*rrhoxhl(k)) *rdzq(kz)
            rimp2(k) = vx(k) + dtpbl *              &
                       vflxp * (rhoxsf*rrhoxhl(k)) *rdzq(kz)
          else
            rimp1(k) = ux(k)
            rimp2(k) = vx(k)
          end if
        end do diffuv
        call solve_tridiag(aimp,bimp,cimp,rimp1,uimp1,kz)
        call solve_tridiag(aimp,bimp,cimp,rimp2,uimp2,kz)

        ! update the winds
        updatewind: &
        do  k = 1 , kz
          ux(k) = uimp1(k)
          vx(k) = uimp2(k)
          ! Recalculate the shear and the squared
          ! magnitude of the shear
          if ( k >= 2 ) then
            dudz =  (ux(k-1)  - ux(k))  *rdza(k)
            dvdz =  (vx(k-1)  - vx(k))  *rdza(k)
            svs(k) = dudz*dudz + dvdz*dvdz
          end if
        end do updatewind
!
!       !Re-update N^2 at the surface; this requires recalculation
!       !of the surface momentum fluxes
!
!       !Calculate surface momentum fluxes
!       uflxp = -uvdragx*ux(kz)/rhoxsf
!       vflxp = -uvdragx*vx(kz)/rhoxsf
!       ustxsq = dsqrt(uflxp*uflxp + vflxp*vflxp)
!
!       !Estimate of surface eddy diffusivity, for estimating the
!       !surface N^2 from the surface virtual heat flux
!       kh0 = vonkar*2*dsqrt(dmax1(tkemin,3.25 * ustxsq))
!
!       !Estimate the surface N^2 from the surface virtual heat flux
!       nsquar(kzp1) = -egrav/thgb*thvflx/kh0


!*******************************************************************************
!*******************************************************************************
!************************* Integration of TKE Budget Equation ******************
!*******************************************************************************
!*******************************************************************************

        ! features:
        !   a. explicit calculation of buoyancy and shear terms using
        !      time=t+1 values of thetal, qw, and winds
        !   b. surface tke is diagnosed
        !   c. semi-implicit calculation of dissipation term and
        !      implicit calculation of turbulent transfer term 
        ! first, buoyancy and shear terms
        sandb: &
        do k = 2 , kz
          ! compute buoyancy term with new values of thetal
          bouyan(k) = -kth(k) * nsquar(k)
          ! compute shear term with new values of svs
          shear(k)  = kzm(k) * svs(k)
        end do sandb
        ! Add radiative divergence contribution to the bouyancy term
        ! (only if there is cloud water at the current level)
        radib:&
        do ilay = 1 , iconv
          k = ktop(ilay)
          if ( qcx(k) > 1.0D-8 .and. k > 1 ) then
            bouyan(k) = bouyan(k) - &
                        rttenx(k)*(presfl(k+1)-presfl(k))*d_1000 * &
                        rrhoxfl(k) * rexnerfl(k) / uthvx(k) 
          end if
        end do radib
        ! tke at top is fixed
        tke(1) = d_zero
        ! diagnose tke at surface, following my 82, b1 ** (2/3) / 2 = 3.25
        !tke(kzp1) = 3.25 * ustx * ustx ! normal
        tke(kzp1) = 3.25D0 * ustxsq ! normal

        ! now the implicit calculations
        ! first find the coefficients that apply for full levels
        imptkeflux: &
        do k = 2 , kz
          if ( k == 2 ) then
            aimp(k-1) = d_zero
          else
            aimp(k-1) = -(rhoxhl(k-1)*rrhoxfl(k))*    &
                          kethl(k-1)*dtpbl*rdzq(k-1)*rdza(k)
          end if
          if ( k == kz ) then
            cimp(k-1) = d_zero
            ! account for implicit part of flux between level kz and surface
            bimp(k-1) = d_one - aimp(k-1) - cimp(k-1) + dtpbl *    &
                        ( dsqrt(tke(k))*rczero/bbls(k) +        &
                        (rhoxhl(k)*rrhoxfl(k))*kethl(k)*rdzq(k)*rdza(k) )
          else
            cimp(k-1) = -(rhoxhl(k)*rrhoxfl(k))*    &
                          kethl(k)*dtpbl*rdzq(k)*rdza(k)
            bimp(k-1) = d_one - aimp(k-1) - cimp(k-1) + dtpbl *   &
                        dsqrt(tke(k))*rczero/dmax1(bbls(k),bbls(k+1))
          end if
          ! now find right side 
          if ( k == kz ) then
            ! account for fixed part of flux between level kz and surface
            rimp1(k-1) = tke(k) + dtpbl * ( shear(k)+bouyan(k) +   &
                           tke(kzp1)*(rhoxhl(k)*rrhoxfl(k))*    &
                           kethl(k)*rdzq(k)*rdza(k) )
          else
            rimp1(k-1) = tke(k) + dtpbl * (shear(k)+bouyan(k))
          end if
        end do imptkeflux
        call solve_tridiag(aimp,bimp,cimp,rimp1,uimp1,kzm1)
        ! update the tke
        updatetke: &
        do  k = 2 , kz
          tke(k) = dmax1(uimp1(k-1),tkemin) ! background tke .001
        end do updatetke

!*******************************************************************************
!*******************************************************************************
!**************** Calculation of Tendencies for Model Output *******************
!*******************************************************************************
!*******************************************************************************

        ! Calculate the TCM tendencies for the model's prognostic variables
        ! For everything but TKE, couple the tendency (multiply by p0-ptop)
        tcmtend: &
        do k = 1 , kz
          ! Zonal wind tendency
          uuwten(i,k,j)= psbx*(ux(k)-uxs(k))*rdtpbl
          ! Meridional wind tendency
          vuwten(i,k,j)= psbx*(vx(k)-vxs(k))*rdtpbl
          ! TKE tendency
          tkeuwten(i,k,j) = (tke(k)-tkes(k))*rdtpbl
          ! Temperature tendency
          tuwten(i,k,j)= psbx*(thx(k)-thxs(k))*exnerhl(k)*rdtpbl
          ! Water vapor tendency
          qvuwten(i,k,j) = psbx*(qx(k)-qxs(k))*rdtpbl
          ! Cloud water tendency
          qcuwten(i,k,j) = psbx*(qcx(k)-qcxs(k))*rdtpbl

          ! Momentum diffusivity
          uwstateb%kzm(j,i,k) = kzm(k)
          ! Scalar diffusivity
          uwstateb%kth(j,i,k) = kth(k)
        end do tcmtend

        ! Output the diagnosed TKE
        uwstatea%srftke(j,i) = tke(kzp1)
        uwstateb%srftke(j,i) = tke(kzp1)
        ! Output the PBL top index and height
        uwstateb%zpbl(j,i) = pblx

        kpbl(j,i) = kpbl2dx

      end do jloop
    end do iloop

  end subroutine uwtcm

  subroutine solve_tridiag(a,b,c,v,x,n)
    ! n - number of equations
    ! a - sub-diagonal (means it is the diagonal below the main diagonal)
    ! b - the main diagonal
    ! c - sup-diagonal (means it is the diagonal above the main diagonal)
    ! v - right part
    ! x - the answer
    ! solves tridiagonal matrix
    ! see http://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm
    ! Written and validated by Travis A. O'Brien 01/04/11
    implicit none
    integer , intent(in) :: n
    real(dp) , dimension(n) , intent(in) :: a , b , c , v
    real(dp) , dimension(n) , intent(out) :: x
    real(dp) , dimension(n) :: bp , vp
    real(dp) :: m
    integer :: i

    bp(1) = b(1)
    vp(1) = v(1)
    ! The first pass (setting coefficients):
    firstpass: &
    do i = 2 , n
      m = a(i)/bp(i-1)
      bp(i) = b(i) - m*c(i-1)
      vp(i) = v(i) - m*vp(i-1)
    end do firstpass
    x(n) = vp(n)/bp(n)
    ! The second pass (back-substition)
    backsub: &
    do i = n-1, 1, -1
      x(i) = (vp(i) - c(i)*x(i+1))/bp(i)
    end do backsub
  end subroutine solve_tridiag 

  subroutine n2(thlxin,qwxin,kmax)
    implicit none
    integer , intent(in) :: kmax
    real(dp) , intent(in) , dimension(kmax) :: thlxin , qwxin
    ! local variables
    real(dp) :: tempv , tvbl , rcld , tvab , thvxfl , dtvdz
    integer :: k

    kloop: &
    do k = 2 , kmax
      ! buoyancy is jump in thetav across flux level/dza
      ! first, layer below, go up and see if anything condenses.
      templ = thlxin(k)*exnerfl(k)
      rvls = esatw(presfl(k),templ)*epop(k)
      temps = templ + (qwxin(k)-rvls)/(cpowlhv +    &
                       ep2*wlhv*rvls/(rgas*templ**d_two))
      rvls = esatw(presfl(k),temps)*epop(k)
      rcldb(k) = dmax1(qwxin(k)-rvls,d_zero)
      tempv = (templ + wlhvocp*rcldb(k)) *    &
              (d_one + ep1*(qwxin(k)-rcldb(k)) - rcldb(k))
      tvbl = tempv*rexnerfl(k) 
      ! now do layer above; go down to see how much evaporates
      templ = thlxin(k-1)*exnerfl(k)
      rvls = esatw(presfl(k),templ)*epop(k)
      temps = templ+(qwxin(k-1)-rvls) / &
                     (cpowlhv+ep2*wlhv*rvls/(rgas*templ**d_two))
      rvls = esatw(presfl(k),temps)*epop(k)
      rcld = dmax1(qwxin(k-1)-rvls,d_zero)
      tempv = (templ + wlhvocp*rcld) *    &
              (d_one + ep1*(qwxin(k-1)-rcld) - rcld)
      tvab = tempv*rexnerfl(k) 
      
      thvxfl= d_half * (tvab+tvbl)
      dtvdz = (tvab - tvbl) *rdza(k)
      nsquar(k) = egrav/thvxfl * dtvdz
    end do kloop
    nsquar(1) = nsquar(2)
  end subroutine n2

  subroutine my(kmax,iconv)
    ! see gb01 and mbg02
    implicit none
    integer , intent(in) :: kmax , iconv
    ! local variables
    real(dp) :: gh , a1ob1 , delthvl , elambda , bige , biga
    real(dp) , parameter :: a1 = 0.92D0 , b1 = 16.6D0 , c1 = 0.08D0 , &
                            a2 = 0.74D0 , b2 = 10.1D0
    integer :: k , ilay
    real(dp) :: kthmax

    a1ob1 = a1/b1

    ! calculate the diffusivities for momentum, thetal, qw and tke
    ! kth and kzm are at full levels.
    ! kethl is at half levels.

    kzm(kmax) = d_zero
    kth(kmax) = d_zero
    kzm(1) = d_zero
    kth(1) = d_zero
    sm(kmax) = d_one
    sh(kmax) = d_one

    kloop: &
    do k = kmax - 1, 2, -1
      gh = -bbls(k)*bbls(k)*nsquar(k)/(d_two*tke(k)+1.0D-9)
      ! gh = dmin1(gh,0.0233D0)
      ! TAO: Added the -0.28 minimum for the G function, as stated
      ! in Galperin (1988), eq 30
      gh = dmax1(dmin1(gh,0.0233D0),-0.28D0)

      sm(k) = (d_one - d_three*c1 - d_six*a1ob1 - d_three*a2*gh*   &
              ((b2-d_three*a2)*(d_one - d_six*a1ob1) -             &
                d_three*c1 * (b2 + d_six*a1))) /                   &
              ((d_one - d_three*a2*gh * (d_six*a1 + b2)) *         &
               (d_one - d_nine*a1*a2*gh))
      sh(k) = a2 * (d_one-d_six*a1ob1) / (d_one-d_three*a2*gh*(d_six*a1+b2))

      ! kzm(k) = dmin1(bbls(k)*dsqrt(2*tke(k))*sm(k),10000.0D0)
      ! kth(k) = dmin1(bbls(k)*dsqrt(2*tke(k))*sh(k),10000.0D0)

      ! Limit the diffusivity to be the vertical grid spacing squared
      ! over the time step; this implies that the entrainment rate
      ! can only be so large that the BL height would change by 
      ! one grid level over one time step -- TAO
      ! kthmax = dmin1((zax(k-1)-zax(k))**2/dt,1.d4)
          
      kthmax = 10000.0D0

      ! Calculate the diffusion coefficients
      kzm(k) = bbls(k)*dsqrt(d_two*tke(k))*sm(k)
      kth(k) = bbls(k)*dsqrt(d_two*tke(k))*sh(k)
      ! Smoothly limit kth to a maximum value
      kth(k) = d_two/mathpi*kthmax*atan(kth(k)/kthmax)
      kzm(k) = d_two/mathpi*kthmax*atan(kzm(k)/kthmax)

      kethl(k)=nuk*dsqrt(kzm(k)*kzm(k+1))
    end do kloop

    ! special case for tops of convective layers
    conv: &
    do ilay = 1 , iconv
      k = ktop(ilay)
      kethl(k) = nuk*kzm(k+1)
      if ( k >= 3 .and. nsquar(k) > 1.D-8 ) then
        kethl(k-1)=0.0D0
        delthvl = (thlx(k-2)+thx(k-2)*ep1*qwx(k-2)) -   &
                  (thlx(k) + thx(k) * ep1 * qwx(k))
        elambda = wlhvocp*rcldb(k)*rexnerhl(k)/dmax1(delthvl,0.1D0)
        bige = 0.8D0 * elambda
        biga = aone * (d_one + atwo * bige)

        ! kth(k) = dmin1(10000.0D0, biga * dsqrt(tke(k)**3)/nsquar(k)/    &
        !          dmax1(bbls(k),bbls(k+1)))
        ! Limit the diffusivity to be the vertical grid spacing squared
        ! over the time step; this implies that the entrainment rate
        ! can only be so large that the BL height would change by 
        ! one grid level over one time step -- TAO
        kthmax = dmin1((zax(k-1)-zax(k))**d_two/dtpbl,1.D4)
        ! kthmax = 10000.0D0
        kth(k) = biga * dsqrt(TKE(k)**d_three)/nsquar(k) /    &
                 dmax1(bbls(k),bbls(k+1))
        ! Smoothly limit kth to a maximum value
        kth(k) = 2*kthmax/mathpi*atan(kth(k)/kthmax)
        kzm(k) = kth(k) / sh(k+1) * sm(k+1) ! prandtl number from layer below
      end if
    end do conv

    ! need kethl at top
    kethl(1) = kethl(2)
    ! replace kethl at surface with something non-zero
    kethl(kmax-1) = nuk*d_half*kzm(kmax-1)
  end subroutine my

  subroutine pblhgt(kmax,iconv)
    ! see mbg02
    implicit none
    ! input variables
    integer , intent(in) :: kmax
    integer , intent(out) :: iconv
    ! local variables
    integer :: istabl , ibeg , ilay , nlev , k , itemp
    real(dp) :: blinf , rnnll , tkeavg , trnnll , radnnll , delthvl , &
                elambda , bige , biga , entnnll , tbbls , lambdal

    ! find noncontiguous convectively unstable layers
    iconv = 0
    istabl = 1

    findconv: &
    do k = 2 , kmax
      if ( nsquar(k) <= d_zero ) then
        if ( istabl == 1 ) then
          iconv = iconv + 1
          ktop(iconv)=k
        end if
        istabl = 0
        kbot(iconv)=k
      else
        istabl = 1
        bbls(k) = dmin1(rstbl*dsqrt(tke(k)/nsquar(k)),vonkar*zqx(k))
        ! bbls(k) = dmax1(dmin1(rstbl*dsqrt(tke(k)/nsquar(k)), &
        !           vonkar*zqx(k)),1.0D-8)
      end if
    end do findconv

    ! now see if they have sufficient buoyant convection to connect
    ibeg = 1

    bigloop: &
    do
      convlayerloop: &
      do ilay = ibeg , iconv
        blinf = xfr*(zqx(ktop(ilay)-1) - zqx(kbot(ilay)+1))
        lambdal = etal*(zqx(ktop(ilay)-1) - zqx(kbot(ilay)+1))
        ! find average n*n*l*l for layer
        rnnll = d_zero
        tkeavg = d_zero
        nlev = kbot(ilay)-ktop(ilay)+1

        avglay1: &
        do k = ktop(ilay) , kbot(ilay)
          bbls(k) = dmin1(blinf,vonkar*zqx(k))
          if ( ilenparam == 1 ) then
            bbls(k) = bbls(k)/(d_one + bbls(k)/lambdal)
          end if
          rnnll = rnnll + nsquar(k)*bbls(k)*bbls(k)
          tkeavg = tkeavg + tke(k) / nlev
        end do avglay1

        ! first extend up
        searchup1: &
        do k = ktop(ilay)-1 , 2 , -1
          ! we always go up at least one, for the entrainment interface
          ktop(ilay) = k 
          bbls(k) = dmin1(vonkar * zqx(k),blinf)
          if ( ilenparam == 1 ) then
            bbls(k) = bbls(k)/(d_one + bbls(k)/lambdal)
          end if
          trnnll = nsquar(k)*bbls(k)*bbls(k)
          ! If this is the top of the layer, stop searching upward
          if ( trnnll*nlev >= -d_half*rnnll ) exit searchup1
          if ( ilay > 1 ) then
            ! did we merge with layer above?
            if ( ktop(ilay) == kbot(ilay-1) ) then
              ibeg = ilay - 1
              ktop(ibeg) = ktop(ibeg)+1
              kbot(ibeg) = kbot(ibeg+1)
              iconv = iconv - 1
              shiftarray1: &
              do itemp = ibeg+1 , iconv
                ktop(itemp)=ktop(itemp+1)
                kbot(itemp)=kbot(itemp+1)
              end do shiftarray1
              cycle bigloop ! recompute for the new, deeper layer
                            ! (restart the do loop)
            end if
          end if
          rnnll = rnnll + trnnll
          nlev = nlev + 1
        end do searchup1
        ! add radiative/entrainment contribution to total
        k = ktop(ilay)
        radnnll = d_zero
        if ( qcx(k) > 1.0D-8 ) then
          radnnll = rttenx(k)*(presfl(k+1)-presfl(k))*d_1000/    &
                    (rhoxfl(k)*uthvx(k)*exnerfl(k))
        end if
        entnnll = d_zero
        if ( k >= 3 ) then
          delthvl = (thlx(k-2)+thx(k-2)*ep1*qwx(k-2))   &
                  - (thlx(k) + thx(k)*ep1*qwx(k))
          elambda = wlhvocp*rcldb(k)*rexnerhl(k)/dmax1(delthvl,0.1D0)
          bige = 0.8D0 * elambda
          biga = aone * (d_one + atwo * bige)
          entnnll = biga * dsqrt(tkeavg**d_three) / bbls(k)
        end if
        rnnll = rnnll + dmin1(d_zero,bbls(k)/dsqrt(dmax1(tkeavg,1.0D-8)) * &
                                            (radnnll + entnnll) )
        nlev = nlev + 1
        ! now extend down
        searchdown1: &
        do k = kbot(ilay)+1 , kmax
          tbbls = dmin1(vonkar * zqx(k),blinf)
          if ( ilenparam == 1 ) then
            tbbls = tbbls/(1 + tbbls/lambdal)
          end if
          trnnll = nsquar(k)*tbbls*tbbls
          ! is it the bottom?
          if ( trnnll*nlev >= -d_half*rnnll ) cycle convlayerloop
          ! (skip the rest of this iteration of the do loop)
          kbot(ilay)=k
          if ( ilay < iconv .and. kbot(ilay) == ktop(ilay+1) ) then
            ! did we merge with layer below?
            ktop(ilay)=ktop(ilay)+1
            kbot(ilay)=kbot(ilay+1)
            iconv = iconv - 1
            shiftarray2: &
            do itemp = ilay+1 , iconv
              ktop(itemp)=ktop(itemp+1)
              kbot(itemp)=kbot(itemp+1)
            end do shiftarray2
            cycle bigloop ! recompute for the new, deeper layer
                          ! (restart the do loop)
          end if
          rnnll = rnnll + trnnll
          bbls(k)=tbbls
          nlev = nlev + 1
        end do searchdown1
      end do convlayerloop
      exit
    end do bigloop

    setbbls: &
    do ilay = 1 , iconv
      blinf = xfr*(zqx(ktop(ilay)-1) - zqx(kbot(ilay)+1))
      lambdal = etal*(zqx(ktop(ilay)-1) - zqx(kbot(ilay)+1))
      convkloop: &
      do k = ktop(ilay) , kbot(ilay)
        bbls(k) = dmin1(vonkar * zqx(k),blinf)
        if ( ilenparam == 1 ) then
          bbls(k) = bbls(k)/(d_one + bbls(k)/lambdal)
        end if
      end do convkloop
    end do setbbls

    ! we should now have tops and bottoms for iconv layers

    if (iconv > 0 ) then
      if ( kbot(iconv) == kmax ) then
        kmix2dx = ktop(iconv)
        if ( kpbl2dx >= 0 ) then
          if ( iconv > 1 ) then
            kpbl2dx = ktop(iconv-1)
          else
            kpbl2dx = kmix2dx
          end if
        else
          kpbl2dx=-kpbl2dx
        end if
      else
        kmix2dx = kmax - 1 
        if ( kpbl2dx >= 0 ) then
          kpbl2dx = ktop(iconv)
        else
          kpbl2dx = -kpbl2dx
        end if
      end if
    else
      kmix2dx = kmax - 1
      if ( kpbl2dx >= 0 ) then
        kpbl2dx = kmix2dx
      else
        kpbl2dx = -kpbl2dx
      end if
    end if


    ! Limit the PBL height to be either the highest allowable (kt--set in
    ! mod_param), or the minimum (kmax-1; the first interface level above the
    ! surface).
    kpbl2dx = max(ktmin,kpbl2dx)
    kpbl2dx = min(kmax - 1,kpbl2dx)
    kmix2dx = max(ktmin,kmix2dx)
    kmix2dx = min(kmax - 1,kmix2dx)

    pblx = zqx(kmix2dx)
  end subroutine pblhgt
    
  ! Returns the saturation vapor pressure over water in units (cb)
  ! given the input pressure in cb and Temperature in K
  ! Modified from Buck (1981), J. App. Met. v 20
  function esatw(p,t)
    implicit none
    real(dp) , intent(in) :: p , t
    real(dp) :: esatw
    real(dp) :: dum , arg , tdum
    ! Limit T to reasonable values.  I believe that this is necessary because
    ! in the iterative calculation of T and QV from the liquid water
    ! temperature, the temperature can take on some crazy values before it
    ! converges. --TAO 01/05/2011
    tdum = dmax1(100.0D0,dmin1(399.99D0,t))
    dum = 1.0007D0 + 3.46D-5*p
    ! arg = 17.502D0*(t-tzero)/(t-32.18D0)
    arg = 17.502D0*(tdum-tzero)/(tdum-32.18D0)
    esatw = dum*0.61121D0*dexp(arg)
  end function esatw

  ! Returns the saturation vapor pressure over ice in units (cb)
  ! given the input pressure in cb and Temperature in K
  ! Modified from Buck (1981), J. App. Met. v 20
  function esati(p,t)
    implicit none
    real(dp) , intent(in) :: p , t
    real(dp) :: esati
    real(dp) :: dum , arg , tdum
    ! Limit T to reasonable values.  I believe that this is necessary because
    ! in the iterative calculation of T and QV from the liquid water
    ! temperature, the temperature can take on some crazy values before it
    ! converges. --TAO 01/05/2011
    tdum = dmax1(100.0D0,dmin1(399.99D0,t))
    dum = 1.0003D0 + 4.18D-5*p
    ! arg = 22.452D0*(t-tzero)/(t-0.6D0)
    arg = 22.452D0*(tdum-tzero)/(tdum-0.6D0)
    esati = dum*0.61115D0*dexp(arg)
  end function esati

  subroutine pblhgt_tao(kmax,iconv)
    implicit none
    ! input variables
    integer , intent(in) :: kmax
    integer , intent(out) :: iconv

    real(dp) , dimension(kmax) :: ktimesz
    real(dp) :: lambda
    logical , dimension(kmax) :: isSaturated1D , isStable1D , isBelow7001D
    logical :: foundlayer
    integer :: k

    isSaturated1D = .false.
    isStable1D = .false.
    isBelow7001D = .false.
    foundlayer = .false.

    where ( rcldb > d_zero )
     isSaturated1D = .true.
    end where

    where ( richnum > rcrit )
      isStable1D = .true.
    end where

    where ( presfl >= 70.0D0 )
      isBelow7001D = .true.
    end where

    ! First see if there is a cloud-topped boundary layer: its top will be
    ! stable, saturated, and below 700 mb
    ctsearch: &
    do k = kmax-1 , 1 , -1
      if ( isSaturated1D(k) .and. isStable1D(k) .and. isBelow7001D(k) ) then
        kmix2dx = k
        foundlayer = .true.
        exit ctsearch
      end if
    end do ctsearch

    ! If we didn't find a cloud-topped boundary layer, then find the first
    ! layer where the richardson number exceeds its threshold
    if ( .not. foundlayer ) then
      drysearch: &
      do k = kmax-1 , 1 , -1
        if ( isStable1D(k) ) then
          kmix2dx = k
          foundlayer = .true.
          exit drysearch
        end if
      end do drysearch
    end if

    ! If we still didn't find a cloud-topped boundary layer, then 
    ! set the top to be the first interface layer
    if ( .not. foundlayer ) then
      kmix2dx = kmax - 1
    end if

    ! Set the boundary layer top and the top of the convective layer
    kmix2dx = max(kmix2dx,3)
    kmix2dx = min(kmix2dx,kmax-1)
    kpbl2dx = kmix2dx
    ! Set that there is only one convective layer
    iconv = 1
    ktop(iconv) = kmix2dx
    ! Set the boundary layer height
    pblx = zqx(ktop(iconv))

    ! Smoothly interpolate the boundary layer height
    if ( (richnum(kmix2dx) >= rcrit) .and.   &
         (richnum(kmix2dx+1) < rcrit) ) then

      pblx = zqx(kmix2dx+1) + (zqx(kmix2dx) - zqx(kmix2dx+1)) * &
             ( (rcrit - richnum(kmix2dx + 1)) /                 &
               (richnum(kmix2dx) - richnum(kmix2dx + 1)) )
    end if

    ! Set the master length scale
    lambda = etal*pblx
    ktimesz = vonkar*zqx(1:kmax)
    ! Within the boundary layer, the length scale is proportional to the height
    bbls = ktimesz/(d_one+ktimesz/lambda)
    ! Otherwise use a Stability-related length scale
    do k = kmix2dx-1 , 1 , -1
      if ( nsquar(k) > d_zero ) then
        bbls(k) = dmax1(dmin1(rstbl*dsqrt(tke(k)/nsquar(k)),ktimesz(k)),1.0D-8)
      else
        bbls(k) = ktimesz(k) - ktimesz(k+1)
      end if
    end do
  end subroutine pblhgt_tao

end module mod_pbl_uwtcm
