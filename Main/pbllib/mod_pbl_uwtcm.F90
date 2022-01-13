!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
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
!   02/2012 TA O'Brien:
!     * Added vertical aerosol mixing to UW TCM

module mod_pbl_uwtcm

  use mod_intkinds
  use mod_realkinds
  use mod_memutil
  use mod_dynparam
  use mod_constants
  use mod_mpmessage
  use mod_mppparam
  use mod_pbl_common
  use mod_pbl_thetal
  use mod_runparams , only : iqv , iqc , iqi , atwo , rstbl , &
          czero , nuk , dt , rdt , ichem , ipptls
  use mod_regcm_types
  use mod_service

  implicit none

  private

  real(rkx) , public , parameter :: uwtkemin = 1.0e-3_rkx

  ! Model constants
  ! fraction of turb layer to be considered in bbls
  real(rkx) , parameter :: xfr = 0.1_rkx
  ! see gb01 regarding the next three lines, atwo and rstbl can be tweaked
  real(rkx) , parameter :: aone = 1.9_rkx*xfr
  ! real(rkx) , parameter :: rcrit = 0.3_rkx
  ! real(rkx) , parameter :: etal =  0.085_rkx
  real(rkx) , parameter :: minn2 =  1.0e-7_rkx

  !real(rkx) , parameter :: svp1 =  0.6112_rkx
  !real(rkx) , parameter :: svp1pa = d_10*svp1
  !real(rkx) , parameter :: svp2 =  17.67_rkx
  !real(rkx) , parameter :: svp3 =  29.65_rkx

  ! Variables that hold frequently-done calculations
  real(rkx) :: rczero , tkefac , b1

  ! To enable original Travis code.
  logical , parameter :: travis_code = .false.
  ! imethod = 1     ! Use the Brent 1973 method to solve for T
  ! imethod = 2     ! Use Bretherton's iterative method to solve for T
  ! imethod = 3     ! Use a finite difference method to solve for T
  integer(ik4) , parameter :: imethod = 2
  ! Do a maximum of 100 iterations (usually are no more than about 30)
  integer(ik4) , parameter :: itbound = 100

  ! Flag for implicit diffusion of ice
  logical , parameter :: implicit_ice = .true.

  public :: allocate_tcm_state
  public :: init_mod_pbl_uwtcm
  public :: uwtcm

  contains

  subroutine allocate_tcm_state(tcmstate)
    implicit none
    type(tcm_state) , intent(out) :: tcmstate
    call getmem3d(tcmstate%kzm,jci1,jci2,ici1,ici2,1,kzp1,'pbl_common:kzm')
    call getmem3d(tcmstate%kth,jci1,jci2,ici1,ici2,1,kzp1,'pbl_common:kth')
  end subroutine allocate_tcm_state

  subroutine init_mod_pbl_uwtcm
    implicit none
    rczero = d_one/czero
    tkefac = czero**(2.0_rkx/3.0_rkx)
    b1 = czero*d_two**(3.0_rkx/2.0_rkx)
  end subroutine init_mod_pbl_uwtcm

  subroutine uwtcm(m2p,p2m)
    implicit none
    type(mod_2_pbl) , intent(in) :: m2p
    type(pbl_2_mod) , intent(inout) :: p2m
    integer(ik4) ::  i , j , k , itr , ibnd
    integer(ik4) :: ilay , kpbconv , iteration
    real(rkx) :: temps , templ , deltat , rvls , pfac , rpfac , tbbls
    real(rkx) :: uflxp , vflxp , rhoxsf , tskx , tvcon , fracz , dudz , &
                 dvdz , thgb , pblx , ustxsq , qfxx , hfxx , uvdragx , &
                 thvflx , q0s , tvfac , svs , cpoxlv
    ! real(rkx) :: kh0
    real(rkx) :: thv0 , thx_t , thvx_t , dthv , dthv_t
    integer(ik4) :: kpbl2dx  ! Top of PBL
    real(rkx) , dimension(kzp2) :: zqx
    real(rkx) , dimension(kzp1) :: kth , kzm , rhoxfl , rcldb , tke , &
               tkes , bbls , nsquar , presfl , exnerfl , rexnerfl
               ! epo , richnum
    real(rkx) , dimension(kz) :: shear , buoyan , rdza , rrhoxfl ! , svs
    real(rkx) , dimension(kz) :: ux , vx , qx , thx , uthvx , zax , kethl , &
               thlx , thlxs , thxs , tx , tvx , rttenx , preshl , qcx ,     &
               qwx , qwxs , rrhoxhl , uxs , qxs , rhoxhl , exnerhl , &
               rexnerhl , rdzq , vxs , qcxs , aimp , bimp , cimp , uimp1 ,  &
               rimp1 , uimp2 , rimp2 , rlv , orlv , cp , ocp
    real(rkx) , dimension(kz) :: qix , qixs
    real(rkx) , dimension(kz,ntr) :: chix , chixs
    real(rkx) , dimension(ntr) :: chifxx
    integer(ik4) , dimension(kz) :: ktop , kbot

#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'uwtcm'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    ! Main do loop
    do i = ici1 , ici2
      do j = jci1 , jci2

!*******************************************************************************
!*******************************************************************************
!*********** Initialization of UW TCM for Current Grid Column ******************
!*******************************************************************************
!*******************************************************************************

        ! Copy in local versions of necessary variables
        if ( idynamic == 3 ) then
          pfac = dt
          rpfac = rdt
        else
          pfac = dt / m2p%psb(j,i)
          rpfac = rdt * m2p%psb(j,i)
        end if
        tskx = m2p%tg(j,i)
        qfxx = m2p%qfx(j,i)
        hfxx = m2p%hfx(j,i)
        uvdragx = m2p%uvdrag(j,i)

        ! Integrate the hydrostatic equation to calculate the level height
        ! Set variables that are on full levels
        do k = 1 , kzp1
          presfl(k) = m2p%patmf(j,i,k)
          !epop(k) = ep2/m2p%patmf(j,i,k)
          zqx(k) = m2p%zq(j,i,k)
          tke(k) = max(m2p%tkests(j,i,k),uwtkemin)
        end do
        zqx(kzp2) = d_zero

        do k = 1 , kz
          tx(k)  = m2p%tatm(j,i,k) + p2m%tten(j,i,k)*pfac
          qx(k)  = max(m2p%qxatm(j,i,k,iqv) + &
                       p2m%qxten(j,i,k,iqv)*pfac,1.0e-8_rkx)
          qcx(k) = max(m2p%qxatm(j,i,k,iqc) + &
                       p2m%qxten(j,i,k,iqc)*pfac,d_zero)
        end do

        do k = 1 , kz
          ux(k)  = m2p%uxatm(j,i,k)
          vx(k)  = m2p%vxatm(j,i,k)
          zax(k) = m2p%za(j,i,k)
          rttenx(k) = m2p%heatrt(j,i,k)
          preshl(k) = m2p%patm(j,i,k)
          rlv(k) = wlh(tx(k))
          cp(k) = cpd*(d_one-qx(k)) + cpv*qx(k)
          orlv(k) = d_one/rlv(k)
          ocp(k) = d_one/cp(k)
        end do

        if ( implicit_ice .and. ipptls > 1 ) then
          do k = 1 , kz
            qix(k) = max(m2p%qxatm(j,i,k,iqi) + &
                         p2m%qxten(j,i,k,iqi)*pfac,d_zero)
          end do
        end if

        if ( ichem == 1 ) then
          do itr = 1 , ntr
            chifxx(itr) = m2p%chifxuw(j,i,itr)
            do k = 1 , kz
              chix(k,itr) = m2p%chib(j,i,k,itr)
            end do
          end do
        end if

        ! Set all the save variables (used for determining the tendencies)
        do k = 1 , kz
          ! Level spacing
          rdzq(k) = d_one/m2p%dzq(j,i,k)
          ! Exner function
          exnerhl(k) = (preshl(k)/p00)**rovcp
          rexnerhl(k) = d_one/exnerhl(k)
          ! Potential temperature
          thx(k) = tx(k)*rexnerhl(k)
          ! Total water mixing ratio
          qwx(k) = qx(k) + qcx(k)
          ! Virtual temperature and potential temperature
          tvcon = d_one + ep1*qx(k)-qcx(k)
          tvx(k) = tx(k)*tvcon ! virtual temperature
          uthvx(k) = thx(k)*tvcon
          ! Liquid water potential temperature (accounting for ice)
          thlx(k) = thx(k) - rlv(k)*qcx(k)*ocp(k)*rexnerhl(k)
        end do

        ! save initial values
        do k = 1 , kz
          uxs(k)  = ux(k)
          vxs(k)  = vx(k)
          thxs(k) = thx(k)
          qxs(k) = qx(k)
          qwxs(k) = qwx(k)
          qcxs(k) = qcx(k)
          thlxs(k) = thlx(k)
        end do
        do k = 1 , kzp1
          tkes(k) = tke(k)
        end do

        ! density at half levels
        do k = 1 , kz
          rhoxhl(k) = preshl(k)/(rgas*tvx(k))
          rrhoxhl(k) = d_one/rhoxhl(k)
        end do

        if ( implicit_ice .and. ipptls > 1 ) then
          do k = 1 , kz
            qixs(k) = qix(k)
          end do
        end if

        if ( ichem == 1 ) then
          do itr = 1 , ntr
            do k = 1 , kz
              chixs(k,itr) = chix(k,itr)
            end do
          end do
        end if

        rdza(1) = d_zero
        rhoxfl(1) = rhoxhl(1)
        do k = 2 , kz
          ! Level spacing
          rdza(k) = d_one/(zax(k-1)-zax(k))
          ! Density
          fracz = (zqx(k)-zax(k))*rdza(k)
          rhoxfl(k) = rhoxhl(k) + (rhoxhl(k-1)-rhoxhl(k))*fracz
        end do
        rhoxfl(kzp1) = presfl(kzp1)/(rgas*tvx(kz))
        do k = 1 , kz
          rrhoxfl(k) = d_one/rhoxfl(k)
        end do

        ! Exner function
        do k = 1 , kzp1
          exnerfl(k) = (presfl(k)/p00)**rovcp
          rexnerfl(k) = d_one/exnerfl(k)
        end do
        ! Wind gradient and shear
        !do k = 1 , kz
        !  dudz = (ux(k-1) - ux(k)) * rdza(k)
        !  dvdz = (vx(k-1) - vx(k)) * rdza(k)
        !  svs(k) = dudz*dudz + dvdz*dvdz
        !end do

!*******************************************************************************
!*******************************************************************************
!*********** Calculation (and Conversion) of Surface Fluxes ********************
!*******************************************************************************
!*******************************************************************************

        ! more surface variables
        thgb = tskx * rexnerfl(kzp1)
        ! Calculate the saturation mixing ratio just above the surface
        !if ( m2p%ldmsk(j,i) == 1 ) then
        !  tvfac = d_one + ep1*m2p%q2m(j,i)
        !else
          !q0s = ep2/(presfl(kzp1)/(d_100*svp1pa * &
          !         exp(svp2*(tskx-tzero)/(tskx-svp3)))-d_one)
          q0s = pfwsat(tskx,presfl(kzp1))
          tvfac = d_one + ep1*q0s
        !end if
        ! density at the surface
        rhoxsf = presfl(kzp1)/(rgas*tvx(kz))
        ! Calculate the virtual temperature right above the surface
        thv0 = thgb * (d_one + tvfac)
        ! Calculate the change in virtual potential temperature from
        ! the surface to the first interface
        dthv = uthvx(kz) - thv0
        ! Calculate the change in potential temperature from the surface
        ! to the first interface
        ! dth = thx(kz) - thgb

        ! Calculate surface momentum fluxes
        uflxp = -uvdragx*ux(kz)/rhoxsf
        vflxp = -uvdragx*vx(kz)/rhoxsf
        ustxsq = sqrt(max(uflxp*uflxp+vflxp*vflxp,1.0e-2_rkx))
        ! Estimate of the surface virtual heat flux
        thvflx = hfxx/rhoxsf*ocp(kz)*tvfac + ep1/thgb*qfxx*rhoxsf
        ! Estimate of surface eddy diffusivity, for estimating the
        ! surface N^2 from the surface virtual heat flux
        ! kh0 = vonkar*d_one*sqrt(max(uwtkemin,tkefac*ustxsq))

!*******************************************************************************
!*******************************************************************************
!********************* Calculation of boundary layer height ********************
!*******************************************************************************
!*******************************************************************************

        ! Calculate nsquared Set N^2 based on the current potential
        ! temperature profile
        ! call n2(thlx,qwx)
        call n2(thlx,qwx)
        ! Estimate the surface N^2 from the surface virtual heat flux
        ! nsquar(kzp1) = -egrav/thgb*thvflx/kh0
        nsquar(kzp1) = egrav/uthvx(kz) * dthv / zax(kz)

        ! Calculate the bulk richardson number
        !richnum = nsquar/max(svs,1.0e-8_rkx)

        ! Calculate the boundary layer height
        call pblhgt(thlx,qwx,kpbconv)
        ! call pblhgt_tao(kpbconv)

!*******************************************************************************
!*******************************************************************************
!************************* Semi-implicit Integration ***************************
!*******************************************************************************
!*******************************************************************************

        ! Calculate the diffusivity profiles, update total water and
        ! liquid water potential temperature, then re-calculate the
        ! diffusivity profiles using the updated values, and re-integrate.
        ! Also update N^2 along the way.
        melloryamadaiteration: &
        do iteration = 1 , 2
          !*************************************************************
          !***** Semi-implicit calculation of diffusivity profiles *****
          !*************************************************************
          call melloryamada(thlx,qwx,kpbconv)
          !*************************************************************
          !****** Implicit Diffusion of Thetal and Qtot ****************
          !*************************************************************
          ! first find the coefficients that apply to
          ! all scalars at half-levels
          aimp(1) = d_zero
          do k = 2 , kz
            aimp(k) = -(rhoxfl(k)*rrhoxhl(k))*kth(k)*dt*rdzq(k)*rdza(k)
          end do
          do k = 1 , kzm1
            cimp(k) = -(rhoxfl(k+1)*rrhoxhl(k))*kth(k+1)*dt*rdzq(k)*rdza(k+1)
          end do
          cimp(kz) = d_zero
          do k = 1 , kz
            bimp(k) = d_one - aimp(k) - cimp(k)
          end do
          do k = 1 , kz
            rimp1(k) = thlx(k)
            rimp2(k) = qwx(k)
          end do
          ! include surface sensible heat flux
          rimp1(kz) = rimp1(kz) + &
                   dt*hfxx*rrhoxhl(kz)*ocp(kz)*rdzq(kz)*rexnerhl(kz)
          ! include surface latent heat flux
          rimp2(kz) = rimp2(kz) + dt*qfxx*rrhoxhl(kz)*rdzq(kz)
          ! Solve total water
          call solve_tridiag(aimp,bimp,cimp,rimp2,uimp2,kz)
          ! Solve liquid water potential temperature
          call solve_tridiag(aimp,bimp,cimp,rimp1,uimp1,kz)
          ! Calculate nsquared Set N^2 based on the updated potential
          ! temperature profile (this is for the semi-implicit integration)
          uimp2 = max(uimp2,minqq)
          call n2(uimp1,uimp2)
          thx_t = uimp1(kz) + ocp(kz)*rlv(kz)*qcx(kz)*rexnerhl(kz)
          tvcon = d_one + ep1*qx(kz)-qcx(kz)
          thvx_t = thx_t*tvcon
          dthv_t = (thvx_t-thv0)
          nsquar(kzp1) = egrav/thvx_t * dthv_t/zax(kz)
          !call pblhgt(uimp1,uimp2,kpbconv)
        end do melloryamadaiteration

        !*************************************************************
        !************ Re-calculate thx, qx, and qcx ******************
        !*************************************************************
        if ( travis_code ) then
          do k = 1 , kz
            thlx(k) = uimp1(k)
            qwx(k) = uimp2(k)
            ! Determine the temperature and qc and qv
            ibnd = itbound
            temps = solve_for_t(thlx(k),qwx(k),preshl(k),thxs(k)*exnerhl(k), &
                                qwxs(k),qcxs(k),thlxs(k),ibnd,imethod,       &
                                qx(k),qcx(k))
            if ( ibnd == -999 ) then
              call fatal(__FILE__,__LINE__, &
                           'UW PBL SCHEME ALGO ERROR')
            end if
            thx(k) = temps*rexnerhl(k)
            uthvx(k) = thx(k)*(d_one + ep1*qx(k)-qcx(k))
          end do
        else
          do k = 1 , kz
            ! Set thlx and qwx to their updated values
            thlx(k) = uimp1(k)
            qwx(k) = uimp2(k)
            templ = thlx(k)*exnerhl(k)
            temps = templ
            !rvls = ep2/(preshl(k)/(d_100*svp1pa * &
            !         exp(svp2*(temps-tzero)/(temps-svp3)))-d_one)
            rvls = pfwsat(temps,preshl(k))
            cpoxlv = cp(k)*orlv(k)
            do iteration = 1 , 3
              deltat = ((templ-temps)*cpoxlv + qwx(k)-rvls) / &
                (cpoxlv + ep2*rlv(k)*rvls/rgas/templ/templ)
              if ( abs(deltat) < 0.01_rkx ) exit
              temps = temps + deltat
              !rvls = ep2/(preshl(k)/(d_100*svp1pa * &
              !       exp(svp2*(temps-tzero)/(temps-svp3)))-d_one)
              rvls = pfwsat(temps,preshl(k))
            end do
            qcx(k) = max(qwx(k)-rvls, d_zero)
            qx(k) = qwx(k)-qcx(k)
            thx(k) = (templ + ocp(k)*rlv(k)*qcx(k))*rexnerhl(k)
            uthvx(k) = thx(k)*(d_one + ep1*qx(k)-qcx(k))
          end do
        end if

        !*************************************************************
        !****** Implicit Diffusion of U and V ************************
        !*************************************************************
        aimp(1) = d_zero
        do k = 2 , kz
          aimp(k) = -(rhoxfl(k)*rrhoxhl(k)) * kzm(k)*dt*rdzq(k)*rdza(k)
        end do
        do k = 1 , kzm1
          cimp(k) = -(rhoxfl(k+1)*rrhoxhl(k)) * kzm(k+1)*dt*rdzq(k)*rdza(k+1)
        end do
        cimp(kz) = d_zero
        do k = 1 , kz
          bimp(k) = d_one - aimp(k) - cimp(k)
        end do
        do k = 1 , kz
          rimp1(k) = ux(k)
          rimp2(k) = vx(k)
        end do
        ! at surface include surface momentum fluxes
        rimp1(kz) = rimp1(kz) + dt * uflxp * (rhoxsf*rrhoxhl(kz))*rdzq(kz)
        rimp2(kz) = rimp2(kz) + dt * vflxp * (rhoxsf*rrhoxhl(kz))*rdzq(kz)
        call solve_tridiag(aimp,bimp,cimp,rimp1,uimp1,kz)
        call solve_tridiag(aimp,bimp,cimp,rimp2,uimp2,kz)

        ! update the winds
        updatewind: &
        do  k = 1 , kz
          ux(k) = uimp1(k)
          vx(k) = uimp2(k)
        end do updatewind

        if ( implicit_ice .and. ipptls > 1 ) then
          !
          ! Implicit diffusion of cloud ice
          !
          aimp(1) = d_zero
          do  k = 2 , kz
            aimp(k) = -(rhoxfl(k)*rrhoxhl(k))*kth(k)*dt*rdzq(k)*rdza(k)
          end do
          do  k = 1 , kzm1
            cimp(k) = -(rhoxfl(k+1)*rrhoxhl(k))*kth(k+1)*dt*rdzq(k)*rdza(k+1)
          end do
          cimp(kz) = d_zero
          do  k = 1 , kz
            bimp(k) = d_one - aimp(k) - cimp(k)
          end do
          do  k = 1 , kz
            rimp1(k) = qix(k)
          end do
          call solve_tridiag(aimp,bimp,cimp,rimp1,uimp1,kz)
          do  k = 1 , kz
            qix(k) = max(uimp1(k),d_zero)
          end do
        end if

        !*************************************************************
        !****** Implicit Diffusion of tracers ************************
        !*************************************************************

        if ( ichem == 1 ) then
          ! Set the tridiagonal coefficients that apply to all of the tracers
          aimp(1) = d_zero
          do k = 2 , kz
            aimp(k) = -(rhoxfl(k)*rrhoxhl(k))*kth(k)*dt*rdzq(k)*rdza(k)
          end do
          do k = 1 , kzm1
            cimp(k) = -(rhoxfl(k+1)*rrhoxhl(k))*   &
                       kth(k+1) * dt*rdzq(k)*rdza(k+1)
          end do
          cimp(kz) = d_zero
          do k = 1 , kz
            bimp(k) = d_one - aimp(k) - cimp(k)
          end do
          !Loop through all of the tracers
          ! set the tridiagonal coefficients that are tracer-specific
          ! and solve the tridiagonal matrix for each tracer to get
          ! the tracer value implied for the next timestep
          do itr = 1 , ntr
            ! set the right side
            do k = 1 , kz
              rimp1(k) = chix(k,itr)
            end do
            ! at surface include surface momentum fluxes
            rimp1(kz) = rimp1(kz) + dt*chifxx(itr)*rrhoxhl(kz)*rdzq(kz)
            !Run the tridiagonal solver
            call solve_tridiag(aimp,bimp,cimp,rimp1,uimp1,kz)

            !Get the chi value implied for the next timestep
            do k = 1 , kz
              chix(k,itr) = max(uimp1(k),d_zero)
            end do
          end do
        end if !End tracer diffusion

        !
        ! Re-update N^2 at the surface; this requires recalculation
        ! of the surface momentum fluxes
        !
        ! Calculate surface momentum fluxes
        uflxp = -uvdragx*ux(kz)/rhoxsf
        vflxp = -uvdragx*vx(kz)/rhoxsf
        ustxsq = sqrt(max(uflxp*uflxp+vflxp*vflxp,1.0e-2_rkx))

        ! Estimate of surface eddy diffusivity, for estimating the
        ! surface N^2 from the surface virtual heat flux
        ! kh0 = vonkar*2*sqrt(max(uwtkemin,tkefac*ustxsq))

        ! Estimate the surface N^2 from the surface virtual heat flux
        dthv = uthvx(kz) - thv0
        nsquar(kzp1) = egrav/uthvx(kz) * dthv / zax(kz)
        ! nsquar(kzp1) = -egrav/thgb*thvflx/kh0
!*******************************************************************************
!*******************************************************************************
!************************* Integration of TKE Budget Equation ******************
!*******************************************************************************
!*******************************************************************************

        ! features:
        !   a. explicit calculation of buoyancy and shear terms using
        !      time = t+1 values of thetal, qw, and winds
        !   b. surface tke is diagnosed
        !   c. semi-implicit calculation of dissipation term and
        !      implicit calculation of turbulent transfer term
        ! first, buoyancy and shear terms
        buoyan(:) = d_zero
        shear(:) = d_zero
        sandb: &
        do k = 2 , kz
          ! Recalculate the shear and the squared
          ! magnitude of the shear
          dudz = (ux(k-1) - ux(k)) * rdza(k)
          dvdz = (vx(k-1) - vx(k)) * rdza(k)
          svs = dudz*dudz + dvdz*dvdz
          ! compute buoyancy term with new values of thetal
          buoyan(k) = -kth(k) * nsquar(k)
          ! compute shear term with new values of svs
          shear(k)  = kzm(k) * svs
        end do sandb
        ! Add radiative divergence contribution to the buoyancy term
        ! (only if there is cloud water at the current level)
        radib:&
        do ilay = 1 , kpbconv
          k = ktop(ilay)
          if ( k > 1 .and. k <= kz ) then
            if ( qcx(k) > 1.0e-4_rkx ) then
              buoyan(k) = buoyan(k) - rttenx(k)*(presfl(k+1)-presfl(k)) * &
                          rrhoxfl(k) * rexnerfl(k) / uthvx(k)
            end if
          end if
        end do radib
        ! tke at top is fixed
        tke(1) = d_zero
        bbls(1)= d_zero
        ! diagnose tke at surface, following my 82, b1 ** (2/3) / 2 = 3.25
        tke(kzp1) = max(tkefac*ustxsq,uwtkemin) ! normal

        ! now the implicit calculations
        ! first find the coefficients that apply for full levels
        imptkeflux: &
        do k = 2 , kz
          if ( k == 2 ) then
            aimp(k-1) = d_zero
          else
            aimp(k-1) = -(rhoxhl(k-1)*rrhoxfl(k))*    &
                          kethl(k-1)*dt*rdzq(k-1)*rdza(k)
          end if
          if ( k == kz ) then
            cimp(k-1) = d_zero
            ! account for implicit part of flux between level kz and surface
            if ( bbls(k) > d_zero ) then
              bimp(k-1) = d_one - aimp(k-1) - cimp(k-1) + dt *    &
                          ( sqrt(tke(k))*rczero/bbls(k) +         &
                          (rhoxhl(k)*rrhoxfl(k))*kethl(k)*rdzq(k)*rdza(k) )
            else
              bimp(k-1) = d_one - aimp(k-1) - cimp(k-1) + dt *    &
                          (rhoxhl(k)*rrhoxfl(k))*kethl(k)*rdzq(k)*rdza(k)
            end if
          else
            cimp(k-1) = -(rhoxhl(k)*rrhoxfl(k))*kethl(k)*dt*rdzq(k)*rdza(k)
            tbbls = max(bbls(k),bbls(k+1))
            if ( tbbls > d_zero ) then
              bimp(k-1) = d_one - aimp(k-1) - cimp(k-1) + &
                           dt * sqrt(tke(k))*rczero/tbbls
            else
              bimp(k-1) = d_one - aimp(k-1) - cimp(k-1)
            end if
          end if
          ! now find right side
          if ( k == kz ) then
            ! account for fixed part of flux between level kz and surface
            rimp1(k-1) = tke(k) + dt * ( shear(k)+buoyan(k) +   &
                           tke(kzp1)*(rhoxhl(k)*rrhoxfl(k))*    &
                           kethl(k)*rdzq(k)*rdza(k) )
          else
            rimp1(k-1) = tke(k) + dt * (shear(k)+buoyan(k))
          end if
        end do imptkeflux
        call solve_tridiag(aimp,bimp,cimp,rimp1,uimp1,kzm1)
        ! update the tke
        do  k = 2 , kz
          tke(k) = max(uimp1(k-1),uwtkemin) ! background tke .001
        end do

!*******************************************************************************
!*******************************************************************************
!**************** Calculation of Tendencies for Model Output *******************
!*******************************************************************************
!*******************************************************************************

        ! Calculate the TCM tendencies for the model's prognostic variables
        ! For everything but TKE, couple the tendency (multiply by pstar)
        do k = 1 , kz
          p2m%tten(j,i,k) = p2m%tten(j,i,k)+rpfac*((thx(k)-thxs(k))*exnerhl(k))
          p2m%qxten(j,i,k,iqv) = p2m%qxten(j,i,k,iqv)+rpfac*(qx(k)-qxs(k))
          p2m%qxten(j,i,k,iqc) = p2m%qxten(j,i,k,iqc)+rpfac*(qcx(k)-qcxs(k))
          p2m%uxten(j,i,k) = (ux(k)-uxs(k))*rdt
          p2m%vxten(j,i,k) = (vx(k)-vxs(k))*rdt
          ! Momentum diffusivity
          uwstate%kzm(j,i,k) = kzm(k)
          ! Scalar diffusivity
          uwstate%kth(j,i,k) = kth(k)
        end do

        do k = 1 , kzp1
          p2m%tketen(j,i,k) = p2m%tketen(j,i,k) + (tke(k)-tkes(k))*rdt
        end do

        if ( implicit_ice .and. ipptls > 1 ) then
          do k = 1 , kz
            p2m%qxten(j,i,k,iqi) = p2m%qxten(j,i,k,iqi)+rpfac*(qix(k)-qixs(k))
          end do
        end if

        if ( ichem == 1 ) then
          do itr = 1 , ntr
            do k = 1 , kz
              p2m%chiten(j,i,k,itr) = p2m%chiten(j,i,k,itr) + &
                         rpfac*(chix(k,itr)-chixs(k,itr))
            end do
          end do
        end if

        p2m%kpbl(j,i) = kpbl2dx
        p2m%zpbl(j,i) = pblx
      end do
    end do

#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif

    contains

#include <wlh.inc>
#include <pfesat.inc>
#include <pfwsat.inc>

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
      integer(ik4) , intent(in) :: n
      real(rkx) , dimension(:) , intent(in) :: a , b , c , v
      real(rkx) , dimension(:) , intent(out) :: x
      real(rkx) , dimension(n) :: bp , vp
      real(rkx) :: m
      integer(ik4) :: i
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

    subroutine n2(thlxin,qwxin)
      implicit none
      real(rkx) , intent(in) , dimension(kz) :: thlxin , qwxin
      ! local variables
      real(rkx) :: tvbl , rcld , tvab , thvxfl , dtvdz
      real(rkx) :: temps , templ , tempv , rvls , cpoxlv
      integer(ik4) :: k

      do k = 2 , kz
        ! buoyancy is jump in thetav across flux level/dza
        ! first, layer below, go up and see if anything condenses.
        templ = thlxin(k)*exnerfl(k)
        !rvls = d_100*svp1pa*exp(svp2*(templ-tzero)/(templ-svp3))*epop(k)
        rvls = pfwsat(templ,presfl(k))
        cpoxlv = cp(k)*orlv(k)
        temps = templ + (qwxin(k)-rvls)/(cpoxlv + &
                      ep2*rlv(k)*rvls/(rgas*templ*templ))
        !rvls = d_100*svp1pa*exp(svp2*(temps-tzero)/(temps-svp3))*epop(k)
        rvls = pfwsat(temps,presfl(k))
        rcldb(k) = max(qwxin(k)-rvls,d_zero)
        tempv = (templ + ocp(k)*rlv(k)*rcldb(k)) * &
                (d_one + ep1*(qwx(k)-rcldb(k))-rcldb(k))
        ! tempv = (templ + wlhvocp*rcldb(k)) * &
        !   (d_one + ep1*(qwx(k)-rcldb(k))-rcldb(k))
        tvbl = tempv*rexnerfl(k)
        ! now do layer above; go down to see how much evaporates
        templ = thlxin(k-1)*exnerfl(k)
        !rvls = d_100*svp1pa*exp(svp2*(templ-tzero)/(templ-svp3))*epop(k)
        rvls = pfwsat(templ,presfl(k))
        temps = templ + (qwxin(k-1)-rvls)/(cpoxlv + &
                ep2*rlv(k)*rvls/(rgas*templ*templ))
        !rvls = d_100*svp1pa*exp(svp2*(temps-tzero)/(temps-svp3))*epop(k)
        rvls = pfwsat(temps,presfl(k))
        rcld = max(qwxin(k-1)-rvls,d_zero)
        !tempv = (templ + wlhvocp*rcld) *    &
        !        (d_one + ep1*(qwxin(k-1)-rcld) - rcld)
        tempv = (templ + ocp(k)*rlv(k)*rcld) * &
                (d_one + ep1*(qwx(k-1)-rcld) - rcld)
        tvab = tempv*rexnerfl(k)

        thvxfl = d_half * (tvab+tvbl)
        dtvdz = (tvab - tvbl) * rdza(k)
        nsquar(k) = egrav/thvxfl * dtvdz
      end do
      nsquar(1) = nsquar(2)
    end subroutine n2

    subroutine melloryamada(thlxin,qwxin,kpbconv)
      implicit none
      integer(ik4) , intent(in) :: kpbconv
      real(rkx) :: gh , a1ob1 , delthvl , elambda , bige , biga
      real(rkx) , intent(in) , dimension(kz) :: thlxin , qwxin
      real(rkx) , parameter :: a1 = 0.92_rkx , c1 = 0.08_rkx , &
                               a2 = 0.74_rkx , b2 = 10.1_rkx
      integer(ik4) :: k , ilay
      real(rkx) , parameter :: kthmax = 1.0e3_rkx
      real(rkx) , parameter :: kzmmax = 1.0e3_rkx
      real(rkx) , dimension(kzp1) :: sm , sh

      a1ob1 = a1/b1

      ! calculate the diffusivities for momentum, thetal, qw and tke
      ! kth and kzm are at full levels.
      ! kethl is at half levels.

      kethl(:) = d_zero
      kzm(kzp1) = d_zero
      kth(kzp1) = d_zero
      kzm(1) = d_zero
      kth(1) = d_zero
      sm(kzp1) = d_one
      sh(kzp1) = d_one
      sm(1) = d_one
      sh(1) = d_one

      kloop: &
      do k = kz , 2, -1
        gh = -bbls(k)*bbls(k)*nsquar(k)/(d_two*tke(k)+1.0e-9_rkx)
        ! TAO: Added the -0.28 minimum for the G function, as stated
        ! in Galperin (1988), eq 30
        gh = max(min(gh,0.0233_rkx),-0.28_rkx)
        sm(k) = a1 * (d_one - 3.0_rkx*c1 - 6.0_rkx*a1ob1 - 3.0_rkx*a2*gh*   &
                     ((b2-3.0_rkx*a2)*(d_one - 6.0_rkx*a1ob1) -             &
                       3.0_rkx*c1 * (b2 + 6.0_rkx*a1))) /                   &
                     ((d_one - 3.0_rkx*a2*gh * (6.0_rkx*a1 + b2)) *         &
                      (d_one - 9.0_rkx*a1*a2*gh))
        sh(k) = a2 * (d_one-6.0_rkx*a1ob1) / &
                     (d_one-3.0_rkx*a2*gh*(6.0_rkx*a1+b2))

        ! kzm(k) = min(bbls(k)*sqrt(2*tke(k))*sm(k),10000.0_rkx)
        ! kth(k) = min(bbls(k)*sqrt(2*tke(k))*sh(k),10000.0_rkx)

        ! Limit the diffusivity to be the vertical grid spacing squared
        ! over the time step; this implies that the entrainment rate
        ! can only be so large that the BL height would change by
        ! one grid level over one time step -- TAO
        ! kthmax = max(min((zax(k-1)-zax(k))**2/dt,1.0e4_rkx),1.0e3_rkx)
        ! kthmax = 1.0e4_rkx

        ! Calculate the diffusion coefficients
        kth(k) = min(bbls(k)*sqrt(d_two*tke(k))*sh(k),kthmax)
        kzm(k) = min(bbls(k)*sqrt(d_two*tke(k))*sm(k),kzmmax)
        ! Smoothly limit kth to a maximum value
        !kth(k) = d_two/mathpi*kthmax*atan(kth(k)/kthmax)
        !kzm(k) = d_two/mathpi*kthmax*atan(kzm(k)/kzmmax)
        kethl(k) = nuk*sqrt(kzm(k)*kzm(k+1))
        kethl(k) = min(kethl(k),kzmmax)
      end do kloop

      ! special case for tops of convective layers
      conv: &
      do ilay = 1 , kpbconv
        k = ktop(ilay)
        if ( nsquar(k) >= minn2 ) then
          kethl(k) = nuk*kzm(k+1)
          if ( k >= 3 ) then
            kethl(k-1) = 0.0_rkx
            delthvl = (thlxin(k-2)+thx(k-2)*ep1*qwxin(k-2)) -   &
                      (thlxin(k) + thx(k)  *ep1*qwxin(k))
            elambda = ocp(k)*rlv(k)*rcldb(k)*rexnerhl(k)/max(delthvl,0.1_rkx)
            bige = 0.8_rkx * elambda
            biga = aone * (d_one + atwo * bige)

            ! kth(k) = min(10000.0_rkx, biga * sqrt(tke(k)**3)/nsquar(k)/    &
            !          max(bbls(k),bbls(k+1)))
            ! Limit the diffusivity to be the vertical grid spacing squared
            ! over the time step; this implies that the entrainment rate
            ! can only be so large that the BL height would change by
            ! one grid level over one time step -- TAO
            ! kthmax = max(min((zax(k-1)-zax(k))**2/dt,1.0e4_rkx),1.0e3_rkx)
            ! kthmax = 1.0e4_rkx
            kth(k) = min(kth(k), biga*sqrt(tke(k)**3)/nsquar(k)/ &
                     max(bbls(k),bbls(k+1)))
            ! kth(k) = biga * sqrt(tke(k)**3)/nsquar(k)/max(bbls(k),bbls(k+1))
            ! Smoothly limit kth to a maximum value
            ! kth(k) = d_two/mathpi*kthmax*atan(kth(k)/kthmax)
            ! prandtl number from layer below
            ! kzm(k) = kth(k) / sh(k+1) * sm(k+1)
            kzm(k) = min(kzm(k),kth(k)/sh(k+1)*sm(k+1))
          end if
        end if
      end do conv

      ! need kethl at top
      kethl(1) = kethl(2)
      ! replace kethl at surface with something non-zero
      kethl(kz) = nuk*d_half*kzm(kz)
    end subroutine melloryamada

    subroutine pblhgt(thlxin,qwxin,kpbconv)
      implicit none
      real(rkx) , intent(in) , dimension(kz) :: thlxin , qwxin
      integer(ik4) , intent(out) :: kpbconv
      integer(ik4) , dimension(kz) :: ktop_save
      integer(ik4) :: istabl , ibeg , ilay , nlev , k , itemp , kstart
      real(rkx) :: blinf , rnnll , tkeavg , trnnll , radnnll , delthvl , &
                  elambda , bige , biga , entnnll , tbbls
      integer(ik4) :: kmix2dx  ! Top of mixed layer (decoupled layer)

      ! find noncontiguous convectively unstable layers

      ktop(:) = 0
      kbot(:) = 0
      kpbl2dx = 0

      kpbconv = 0
      istabl = 1
      do k = 1 , kzp1
        bbls(k) = d_zero
      end do
      do k = 2 , kzp1
        if ( nsquar(k) <= d_zero ) then
          if ( istabl == 1 ) then
            kpbconv = kpbconv + 1
            ktop(kpbconv) = k
          end if
          istabl = 0
          kbot(kpbconv) = k
        else
          istabl = 1
          bbls(k) = min(rstbl*sqrt(tke(k)/nsquar(k)),vonkar*zqx(k))
        end if
      end do

      ! now see if they have sufficient buoyant convection to connect
      ktop_save(:) = ktop(:)
      if ( kpbconv >= 1 ) then
        ibeg = 1
        convlayerloop: &
        do ilay = ibeg , kpbconv
          blinf = xfr*(zqx(ktop(ilay)-1) - zqx(kbot(ilay)+1))
          ! find average n*n*l*l for layer
          rnnll = d_zero
          tkeavg = d_zero
          nlev = kbot(ilay)-ktop(ilay)+1
          ! Average layer
          do k = ktop(ilay) , kbot(ilay)
            bbls(k) = min(blinf,vonkar*zqx(k))
            rnnll = rnnll + nsquar(k)*bbls(k)*bbls(k)
            tkeavg = tkeavg + tke(k) / real(nlev,rkx)
          end do
          ! first extend up
          kstart = ktop(ilay) - 1
          searchup1: &
          do k = kstart , 2 , -1
            ! we always go up at least one, for the entrainment interface
            ktop(ilay) = k
            bbls(k) = min(blinf,vonkar*zqx(k))
            trnnll = nsquar(k)*bbls(k)*bbls(k)
            ! If this is the top of the layer, stop searching upward
            if ( trnnll*real(nlev,rkx) >= -d_half*rnnll ) exit searchup1
            if ( ilay > 1 ) then
              ! did we merge with layer above?
              if ( ktop(ilay) == kbot(ilay-1) ) then
                ibeg = ilay - 1
                ktop(ibeg) = ktop_save(ibeg)
                kbot(ibeg) = kbot(ibeg+1)
                kpbconv = kpbconv - 1
                do itemp = ibeg+1 , kpbconv
                  ktop(itemp) = ktop(itemp+1)
                  kbot(itemp) = kbot(itemp+1)
                  ktop_save(itemp) = ktop_save(itemp+1)
                end do
                cycle convlayerloop ! recompute for the new, deeper layer
                                    ! (restart the do loop)
              end if
            end if
            rnnll = rnnll + trnnll
            nlev = nlev + 1
          end do searchup1
          ! add radiative/entrainment contribution to total
          k = ktop(ilay)
          if ( qcx(k) > 1.0e-8_rkx ) then
            radnnll = rttenx(k)*(presfl(k+1)-presfl(k)) /  &
                      (rhoxfl(k)*uthvx(k)*exnerfl(k))
          else
            radnnll = d_zero
          end if
          entnnll = d_zero
          if ( k >= 3 ) then
            delthvl = (thlxin(k-2)+thx(k-2)*ep1*qwxin(k-2)) - &
                      (thlxin(k)  +thx(k)  *ep1*qwxin(k))
            elambda = ocp(k)*rlv(k)*rcldb(k)*rexnerhl(k)/max(delthvl,0.1_rkx)
            bige = 0.8_rkx * elambda
            biga = aone * (d_one + atwo * bige)
            entnnll = biga * sqrt(tkeavg**3) / bbls(k)
          end if
          if ( tkeavg > d_zero ) then
            rnnll = rnnll + min(d_zero,bbls(k)/sqrt(tkeavg)*(radnnll+entnnll))
          end if
          ! now extend down
          searchdown1: &
          do k = kbot(ilay)+1 , kzp1
            tbbls = min(blinf,vonkar*zqx(k))
            trnnll = nsquar(k)*tbbls*tbbls
            ! is it the bottom?
            if ( trnnll*real(nlev,rkx) >= -d_half*rnnll ) exit searchdown1
            ! (skip the rest of this iteration of the do loop)
            kbot(ilay) = k
            if ( ilay < kpbconv .and. kbot(ilay) == ktop(ilay+1) ) then
              ! did we merge with layer below?
              ktop(ilay) = ktop_save(ilay)
              kbot(ilay) = kbot(ilay+1)
              kpbconv = kpbconv - 1
              shiftarray2: &
              do itemp = ilay+1 , kpbconv
                ktop(itemp) = ktop(itemp+1)
                kbot(itemp) = kbot(itemp+1)
                ktop_save(itemp) = ktop_save(itemp+1)
              end do shiftarray2
              cycle convlayerloop ! recompute for the new, deeper layer
                                  ! (restart the do loop)
            end if
            rnnll = rnnll + trnnll
            bbls(k) = tbbls
            nlev = nlev + 1
          end do searchdown1
        end do convlayerloop
        setbbls: &
        do ilay = 1 , kpbconv
          blinf = xfr*(zqx(ktop(ilay)-1) - zqx(kbot(ilay)+1))
          convkloop: &
          do k = ktop(ilay) , kbot(ilay)
            bbls(k) = min(blinf,vonkar*zqx(k))
          end do convkloop
        end do setbbls
      end if
      ! we should now have tops and bottoms for kpbconv layers
      if ( kpbconv > 0 ) then
        if ( kbot(kpbconv) == kzp1 ) then
          kmix2dx = ktop(kpbconv)
          if ( kpbl2dx >= 0 ) then
            if ( kpbconv > 1 ) then
              kpbl2dx = ktop(kpbconv-1)
            else
              kpbl2dx = kmix2dx
            end if
          else
            kpbl2dx=-kpbl2dx
          end if
        else
          kmix2dx = kz
          if ( kpbl2dx >= 0 ) then
            kpbl2dx = ktop(kpbconv)
          else
            kpbl2dx = -kpbl2dx
          end if
        end if
      else
        ! Lowermost layer
        kmix2dx = kz
        if ( kpbl2dx >= 0 ) then
          kpbl2dx = kpbl2dx
        else
          kpbl2dx = -kpbl2dx
        end if
      end if
      pblx = zqx(kmix2dx)
    end subroutine pblhgt

    ! Returns the saturation vapor pressure over water in units (cb)
    ! given the input pressure in cb and Temperature in K
    ! Modified from Buck (1981), J. App. Met. v 20
    !function esatw(p,t)
    !  implicit none
    !  real(rkx) , intent(in) :: p , t
    !  real(rkx) :: esatw
    !  real(rkx) :: dum , arg , tdum
    !  ! Limit T to reasonable values.  I believe that this is necessary because
    !  ! in the iterative calculation of T and QV from the liquid water
    !  ! temperature, the temperature can take on some crazy values before it
    !  ! converges. --TAO 01/05/2011
    !  tdum = max(100.0_rkx,min(399.99_rkx,t))
    !  dum = 1.0007_rkx + 3.46e-5_rkx*p
    !  ! arg = 17.502_rkx*(t-tzero)/(t-32.18_rkx)
    !  arg = 17.502_rkx*(tdum-tzero)/(tdum-32.18_rkx)
    !  esatw = dum*0.61121_rkx*exp(arg)
    !end function esatw

    ! Returns the saturation vapor pressure over ice in units (cb)
    ! given the input pressure in cb and Temperature in K
    ! Modified from Buck (1981), J. App. Met. v 20
    !function esati(p,t)
    !  implicit none
    !  real(rkx) , intent(in) :: p , t
    !  real(rkx) :: esati
    !  real(rkx) :: dum , arg , tdum
    !  ! Limit T to reasonable values.  I believe that this is necessary because
    !  ! in the iterative calculation of T and QV from the liquid water
    !  ! temperature, the temperature can take on some crazy values before it
    !  ! converges. --TAO 01/05/2011
    !  tdum = max(100.0_rkx,min(399.99_rkx,t))
    !  dum = 1.0003_rkx + 4.18e-5_rkx*p
    !  ! arg = 22.452_rkx*(t-tzero)/(t-0.6_rkx)
    !  arg = 22.452_rkx*(tdum-tzero)/(tdum-0.6_rkx)
    !  esati = dum*0.61115_rkx*exp(arg)
    !end function esati

!    subroutine pblhgt_tao(kpbconv)
!      implicit none
!      integer(ik4) , intent(out) :: kpbconv
!      real(rkx) , dimension(kz) :: ktimesz
!      real(rkx) :: lambda
!      logical , dimension(kz) :: issaturated1d , isstable1d , isbelow7001d
!      logical :: foundlayer
!      integer(ik4) :: k
!      issaturated1d = .false.
!      isstable1d = .false.
!      isbelow7001d = .false.
!      foundlayer = .false.
!      where ( rcldb > d_zero )
!       issaturated1d = .true.
!      end where
!      where ( richnum > rcrit )
!        isstable1d = .true.
!      end where
!      where ( presfl >= 70000.0_rkx )
!        isbelow7001d = .true.
!      end where
!      ! First see if there is a cloud-topped boundary layer: its top will be
!      ! stable, saturated, and below 700 mb
!      do k = kzm1 , 1 , -1
!        if ( issaturated1d(k) .and. isstable1d(k) .and. isbelow7001d(k) ) then
!          kmix2dx = k
!          foundlayer = .true.
!          exit
!        end if
!      end do
!      ! If we didn't find a cloud-topped boundary layer, then find the first
!      ! layer where the richardson number exceeds its threshold
!      if ( .not. foundlayer ) then
!        do k = kzm1 , 1 , -1
!          if ( isstable1d(k) ) then
!            kmix2dx = k
!            foundlayer = .true.
!            exit
!          end if
!        end do
!      end if
!      ! If we still didn't find a cloud-topped boundary layer, then
!      ! set the top to be the first interface layer
!      if ( .not. foundlayer ) then
!        kmix2dx = kz - 1
!      end if
!      ! Set the boundary layer top and the top of the convective layer
!      kmix2dx = max(kmix2dx,3)
!      kmix2dx = min(kmix2dx,kzm1)
!      kpbl2dx = kmix2dx
!      ! Set that there is only one convective layer
!      kpbconv = 1
!      ktop(kpbconv) = kmix2dx
!      ! Set the boundary layer height
!      pblx = zqx(ktop(kpbconv))
!      ! Smoothly interpolate the boundary layer height
!      if ( (richnum(kmix2dx) >= rcrit) .and. &
!           (richnum(kmix2dx+1) < rcrit) ) then
!        pblx = zqx(kmix2dx+1) + (zqx(kmix2dx) - zqx(kmix2dx+1)) * &
!                 ( (rcrit - richnum(kmix2dx+1)) /                 &
!                   (richnum(kmix2dx) - richnum(kmix2dx+1)) )
!      end if
!      ! Set the master length scale
!      lambda = etal*pblx
!      ktimesz = vonkar*zqx(1:kz)
!      ! Within the boundary layer, the length scale is propor. to the height
!      bbls = ktimesz/(d_one+ktimesz/lambda)
!      ! Otherwise use a Stability-related length scale
!      do k = kmix2dx-1 , 1 , -1
!        if ( nsquar(k) > d_zero ) then
!          bbls(k) = max(min(rstbl*sqrt(tke(k) / &
!                        nsquar(k)),ktimesz(k)),1.0e-8_rkx)
!        else
!          bbls(k) = ktimesz(k) - ktimesz(k+1)
!        end if
!      end do
!    end subroutine pblhgt_tao

  end subroutine uwtcm

end module mod_pbl_uwtcm
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
