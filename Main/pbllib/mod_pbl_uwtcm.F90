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
          czero , nuk , dt , rdt , ichem , sigma , hsigma , dsigma , ipptls
  use mod_regcm_types
  use mod_service

  implicit none

  private

  real(rk8) , public , parameter :: tkemin = 1.0D-3

  ! Model constants
  ! fraction of turb layer to be considered in bbls
  real(rk8) , parameter :: xfr = 0.1D0
  ! see gb01 regarding the next three lines, atwo and rstbl can be tweaked
  real(rk8) , parameter :: aone = 1.9D0*xfr
  ! b1/2**(3/2) from Mellor and Yamada (1982)
  real(rk8) , parameter :: rcrit = 0.3D0
  real(rk8) , parameter :: etal =  0.085D0
  real(rk8) , parameter :: minn2 =  1.0D-7

  ! Variables that hold frequently-done calculations
  real(rk8) :: rczero

  ! local variables on full levels
  real(rk8) , pointer , dimension(:) :: zqx , kth , kzm , rhoxfl , &
                 tke , tkes , rrhoxfl , bbls , nsquar , richnum , &
                 bouyan , rdza , dza , svs , presfl , exnerfl ,   &
                 shear , rexnerfl , rcldb , epop , sm , sh, xkchi

  ! local variables on half levels
  real(rk8) , pointer , dimension(:) :: ux , vx , thx , qx , uthvx ,   &
                 zax , kethl , thlx , thlxs, thxs , tx , tvx ,         &
                 rttenx , preshl , qcx , qwx , qwxs , udzq , rrhoxhl , &
                 uxs , qxs , rhoxhl , exnerhl , rexnerhl , rdzq ,      &
                 vxs , qcxs , aimp , bimp , cimp , uimp1 , rimp1 ,     &
                 uimp2 , rimp2 , qix , qixs
  real(rk8) , pointer , dimension(:,:) :: chix , chixs

  integer(ik4) , pointer , dimension(:) :: isice , ktop , kbot

  ! local scalars
  real(rk8) :: uflxp , vflxp , rhoxsf , tgbx , tvcon , fracz , dudz , &
              dvdz , rvls , thv0 , dthv , psbx , templ , temps ,     &
              cell , thgb , pblx , ustxsq , qfxx , hfxx , dth ,      &
              uvdragx , thvflx , kh0 , q0s
  real(rk8) , pointer , dimension(:) :: chifxx

  integer(ik4) :: kpbl2dx  ! Top of PBL
  integer(ik4) :: kmix2dx  ! Top of mixed layer (decoupled layer)

  ! imethod = 1     ! Use the Brent 1973 method to solve for T
  ! imethod = 2     ! Use Bretherton's iterative method to solve for T
  ! imethod = 3     ! Use a finite difference method to solve for T
  integer(ik4) , parameter :: imethod = 2
  ! Do a maximum of 100 iterations (usually are no more than about 30)
  integer(ik4) , parameter :: itbound = 100
  ! Modding : Update tendencies only in the PBL
  logical :: update_only_in_pbl = .false.

  public :: allocate_tcm_state
  public :: init_mod_pbl_uwtcm
  public :: uwtcm
  public :: get_data_from_tcm

  contains

#include <pfesat.inc>
#include <pfqsat.inc>

  subroutine allocate_tcm_state(tcmstate)
    implicit none
    type(tcm_state) , intent(out) :: tcmstate
    call getmem3d(tcmstate%tkeps,jce1,jce2,ice1,ice2,1,kzp1,'pbl_common:tkeps')
    call getmem3d(tcmstate%advtke,jce1,jce2, &
                  ice1,ice2,1,kzp1,'pbl_common:advtke')
    call getmem3d(tcmstate%kzm,jci1,jci2,ici1,ici2,1,kzp1,'pbl_common:kzm')
    call getmem3d(tcmstate%kth,jci1,jci2,ici1,ici2,1,kzp1,'pbl_common:kth')
    call getmem2d(tcmstate%srftke,jci1,jci2,ici1,ici2,'pbl_common:srftke')
  end subroutine allocate_tcm_state

  subroutine get_data_from_tcm(p2m,atm1,atm2)
    implicit none
    type(pbl_2_mod) , intent(inout) :: p2m
    type(atmstate_a) , intent(inout) :: atm1
    type(atmstate_b) , intent(inout) :: atm2
    integer :: i , j , k , n
    ! Don't update the model variables if we are the diagnostic mode
    ! (Holtslag running, and UW updating tke)
    !
    ! Put the t and qv tendencies in to difft and diffq for
    ! application of the sponge boundary conditions (see mod_tendency)
    !
    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          p2m%diffqx(j,i,k,iqv) = p2m%diffqx(j,i,k,iqv) +  &
             p2m%qxuwten(j,i,k,iqv)
        end do
      end do
    end do

    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          p2m%difft(j,i,k) = p2m%difft(j,i,k) + p2m%tuwten(j,i,k)
        end do
      end do
    end do

    ! Put the cloud water tendency
    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          p2m%qxten(j,i,k,iqc) = p2m%qxten(j,i,k,iqc) + &
             p2m%qxuwten(j,i,k,iqc)
        end do
      end do
    end do
    if ( ipptls == 2 ) then
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            p2m%qxten(j,i,k,iqi) = p2m%qxten(j,i,k,iqi) + &
               p2m%qxuwten(j,i,k,iqi)
          end do
        end do
      end do
    end if


    ! Put the tracer tendencies in chiuwten
    ! TODO: may want to calcuate rmdr here following holtbl
    if ( ichem == 1 ) then
      do n = 1 , ntr
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
              p2m%chiten(j,i,k,n) = p2m%chiten(j,i,k,n) + &
                 chiuwten(j,i,k,n)
            end do
          end do
        end do
      end do
    end if

    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          p2m%tketen(j,i,k) = p2m%tkeuwten(j,i,k)
        end do
      end do
    end do
    !
    ! Set the surface tke (diagnosed)
    !
    do i = ici1 , ici2
      do j = jci1 , jci2
        atm1%tke(j,i,kzp1) = uwstateb%srftke(j,i)
        atm2%tke(j,i,kzp1) = uwstateb%srftke(j,i)
      end do
    end do
  end subroutine get_data_from_tcm

  subroutine init_mod_pbl_uwtcm
    implicit none

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
    call getmem1d(qix,1,kz,'mod_uwtcm:qix')
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
    call getmem1d(qixs,1,kz,'mod_uwtcm:qixs')
    call getmem1d(aimp,1,kz,'mod_uwtcm:aimp')
    call getmem1d(bimp,1,kz,'mod_uwtcm:bimp')
    call getmem1d(cimp,1,kz,'mod_uwtcm:cimp')
    call getmem1d(uimp1,1,kz,'mod_uwtcm:uimp1')
    call getmem1d(rimp1,1,kz,'mod_uwtcm:rimp1')
    call getmem1d(uimp2,1,kz,'mod_uwtcm:uimp2')
    call getmem1d(rimp2,1,kz,'mod_uwtcm:rimp2')

    if ( ichem == 1 ) then
      call getmem1d(xkchi,1,kzp1,'mod_uwtcm:xkchi')
      call getmem2d(chix,1,kz,1,ntr,'mod_uwtcm:chix')
      call getmem2d(chixs,1,kz,1,ntr,'mod_uwtcm:chixs')
      call getmem1d(chifxx,1,ntr,'mod_uwtcm:chifxx')
    end if

    call getmem1d(isice,1,kz,'mod_uwtcm:isice')
    call getmem1d(ktop,1,kz,'mod_uwtcm:ktop')
    call getmem1d(kbot,1,kz,'mod_uwtcm:kbot')

    rczero = d_one/czero
  end subroutine init_mod_pbl_uwtcm

  subroutine uwtcm(m2p,p2m)
    implicit none
    type(mod_2_pbl) , intent(in) :: m2p
    type(pbl_2_mod) , intent(inout) :: p2m
    integer(ik4) ::  i , j , k , itr , ibnd
    integer(ik4) :: ilay , kpbconv , iteration
    real(rk8) :: temps , templ , deltat , rvls
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'uwtcm'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    ! Zero out tendencies

    p2m%uuwten(:,:,:) = d_zero
    p2m%vuwten(:,:,:) = d_zero
    p2m%tkeuwten(:,:,:) = d_zero
    p2m%tuwten(:,:,:) = d_zero
    p2m%qxuwten(:,:,:,:) = d_zero
    if ( ichem == 1 ) chiuwten(:,:,:,:) = d_zero

    !Main do loop
    do i = ici1 , ici2
      do j = jci1 , jci2

!*******************************************************************************
!*******************************************************************************
!*********** Initialization of UW TCM for Current Grid Column ******************
!*******************************************************************************
!*******************************************************************************

        ! Copy in local versions of necessary variables
        psbx = m2p%psb(j,i)
        tgbx = m2p%tgb(j,i)
        qfxx = m2p%qfx(j,i)
        hfxx = m2p%hfx(j,i)
        uvdragx = m2p%uvdrag(j,i)

        ! Integrate the hydrostatic equation to calculate the level height
        ! Set variables that are on full levels
        do k = 1 , kzp1
          presfl(k) = m2p%patmf(j,i,k)
          zqx(k) = m2p%zq(j,i,k)
        end do

        tke(kzp1) = m2p%tkests(j,i,kzp1)
        do k = 1 , kz
          tx(k)  = m2p%tatm(j,i,k)
          qx(k)  = m2p%qxatm(j,i,k,iqv)
          qcx(k) = m2p%qxatm(j,i,k,iqc)
          zax(k) = m2p%za(j,i,k)
          tke(k) = m2p%tkests(j,i,k)
          ux(k)  = m2p%uxatm(j,i,k)
          vx(k)  = m2p%vxatm(j,i,k)
          rttenx(k) = m2p%heatrt(j,i,k)
        end do
        if ( ipptls == 2 ) then
          do k = 1 , kz
            qix(k) = m2p%qxatm(j,i,k,iqi)
          end do
        else
          do k = 1 , kz
            qix(k) = d_zero
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
        tkes(kzp1) = tke(kzp1)
        khalfloop: &
        do k = 1 , kz
          ! pressure at half levels
          preshl(k) = m2p%patm(j,i,k)
          ! Level spacing
          udzq(k) = m2p%dzq(j,i,k)
          rdzq(k) = d_one/udzq(k)
          ! Exner function
          exnerhl(k) = (preshl(k)/1.0D5)**rovcp
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
          thlx(k) = thx(k) - wlhvocp * rexnerhl(k) * qcx(k)
          thlxs(k) = thlx(k)
          ! save initial values of winds, thx, and qx
          tkes(k) = tke(k)
          uxs(k)  = ux(k)
          vxs(k)  = vx(k)
          thxs(k) = thx(k)
          qxs(k) = qx(k)
          qcxs(k) = qcx(k)
          qixs(k) = qix(k)
          ! density at half levels
          rhoxhl(k) = preshl(k)/(rgas*tvx(k))
          rrhoxhl(k) = d_one/rhoxhl(k)
        end do khalfloop

        if ( ichem == 1 ) then
          do itr = 1 , ntr
            do k = 1 , kz
              chixs(k,itr) = chix(k,itr)
            end do
          end do
        end if

        do k = 2 , kz
          ! pressure at full levels
          epop(k) = ep2/presfl(k)
          ! Level spacing
          dza(k) = zax(k-1)-zax(k)
          rdza(k) = d_one/dza(k)
          ! Exner function
          exnerfl(k) = (presfl(k)/1.0D5)**rovcp
          rexnerfl(k) = d_one/exnerfl(k)
          ! Density
          fracz = (zqx(k)-zax(k))*rdza(k)
          rhoxfl(k) = rhoxhl(k) + (rhoxhl(k-1)-rhoxhl(k))*fracz
          rrhoxfl(k) = d_one/rhoxfl(k)
          ! Wind gradient and shear
          dudz =  (ux(k-1) - ux(k)) * rdza(k)
          dvdz =  (vx(k-1) - vx(k)) * rdza(k)
          svs(k) = dudz*dudz + dvdz*dvdz
        end do

        ! Surface exner function
        exnerfl(kzp1) = (presfl(kzp1)/1.0D5)**rovcp
        rexnerfl(kzp1) = d_one/exnerfl(kzp1)
        ! density at the surface
        rhoxsf = presfl(kzp1)/(rgas*tvx(kz))

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
        q0s = ep2/(psbx/pfesat(tgbx) - d_one)
        ! Calculate the virtual temperature right above the surface
        thv0 = thgb * (d_one + ep1*q0s)
        ! Calculate the change in virtual potential temperature from
        ! the surface to the first interface
        dthv = uthvx(kz) - thv0
        ! Calculate the change in potential temperature from the surface
        ! to the first interface
        dth = thx(kz) - thgb

        ! Calculate surface momentum fluxes
        uflxp = -uvdragx*ux(kz)/rhoxsf
        vflxp = -uvdragx*vx(kz)/rhoxsf
        ustxsq = sqrt(uflxp*uflxp + vflxp*vflxp)

        ! Estimate of the surface virtual heat flux
        thvflx = hfxx/rhoxsf*rcpd*(d_one+ep1*q0s) + ep1/thgb*qfxx*rhoxsf
        ! Estimate of surface eddy diffusivity, for estimating the
        ! surface N^2 from the surface virtual heat flux
        kh0 = vonkar*d_one*sqrt(max(tkemin,3.25D0 * ustxsq))

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
        richnum = nsquar/max(svs,1.0D-8)

        ! Calculate the boundary layer height
        call pblhgt(kzp1,m2p%ktrop(j,i),kpbconv)
        ! call pblhgt_tao(kzp1,kpbconv)

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
        do iteration = 0 , 1
          !*************************************************************
          !***** Semi-implicit calculation of diffusivity profiles *****
          !*************************************************************
          call melloryamada(kzp1,kpbconv)

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
              aimp(k) = -(rhoxfl(k)*rrhoxhl(k))*kth(k)*dt*rdzq(k)*rdza(k)
            end if
            if ( k == kz ) then
              cimp(k) = d_zero
            else
              cimp(k) = -(rhoxfl(k+1)*rrhoxhl(k)) * &
                kth(k+1)*dt*rdzq(k)*rdza(k+1)
            end if
            bimp(k) = d_one - aimp(k) - cimp(k)
            ! now find right side for various scalars:
            ! no flux out top, so no (k == 1)
            if ( k == kz ) then ! at surface
              ! include surface latent heat flux
              rimp2(k) = qwx(k) + dt * qfxx*rrhoxhl(k)*rdzq(kz)
              ! include surface sensible heat flux
              rimp1(k) = thlx(k) + &
                     dt * hfxx*rrhoxhl(k)*rcpd*rdzq(kz)*rexnerhl(kz)
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
        end do melloryamadaiteration

        !*************************************************************
        !************ Re-calculate thx, qx, and qcx ******************
        !*************************************************************

        do k = 1 , kz
          ! Set thlx and qwx to their updated values
          thlx(k) = uimp1(k)
          qwx(k) = uimp2(k)
          templ = thlx(k)*exnerhl(k)
          temps = templ
          rvls = ep2/(preshl(k)/pfesat(temps)-d_one)
          do iteration = 1 , 3
            deltat = ((templ-temps)*cpowlhv + qwx(k)-rvls)/   &
                      (cpowlhv+ep2*wlhv*rvls/rgas/temps/temps)
            temps = temps + deltat
            rvls = ep2/(preshl(k)/pfesat(temps)-d_one)
          end do
          if ( ipptls == 2 ) then
            qcx(k) = min(1.0D-3,max(qwx(k)-rvls, d_zero))
          else
            qcx(k) = max(qwx(k)-rvls, d_zero)
          end if
          qx(k) = qwx(k) - qcx(k)
          thx(k) = (templ + qcx(k)*wlhvocp) *rexnerhl(k)
          uthvx(k) = thx(k)*(d_one + ep1*qx(k)-qcx(k))
          ! Determine the temperature and qc and qv
          !temps = solve_for_t(thlx(k),qwx(k),preshl(k),thxs(k)*exnerhl(k), &
          !                    qwxs(k),qcxs(k),thlxs(k),itbound,imethod,    &
          !                    qx(k),qcx(k))
          ! thx(k) = temps*rexnerhl(k)
          ! uthvx(k) = thx(k)*(d_one + ep1*qx(k)-qcx(k))
        end do

        !*************************************************************
        !****** Implicit Diffusion of U and V ************************
        !*************************************************************
        diffuv: &
        do k = 1 , kz
          if ( k == 1 ) then
            aimp(k) = d_zero
          else
            aimp(k) = -(rhoxfl(k)*rrhoxhl(k)) * kzm(k)*dt*rdzq(k)*rdza(k)
          end if
          if ( k == kz ) then
            cimp(k) = d_zero
          else
            cimp(k) = -(rhoxfl(k+1)*rrhoxhl(k)) * kzm(k+1)*dt*rdzq(k)*rdza(k+1)
          end if
          bimp(k) = d_one - aimp(k) - cimp(k)
          ! now find right side
          ! no flux out top, so no (k == 1)
          if ( k == kz ) then
            ! at surface include surface momentum fluxes
            rimp1(k) = ux(k) + dt * uflxp * (rhoxsf*rrhoxhl(k)) *rdzq(kz)
            rimp2(k) = vx(k) + dt * vflxp * (rhoxsf*rrhoxhl(k)) *rdzq(kz)
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
            dudz =  (ux(k-1) - ux(k)) * rdza(k)
            dvdz =  (vx(k-1) - vx(k)) * rdza(k)
            svs(k) = dudz*dudz + dvdz*dvdz
          end if
        end do updatewind

        if ( ipptls == 2 ) then
          !
          ! Implicit diffusion of cloud ice
          !
          do  k = 1 , kz
            if ( k == 1 ) then
              aimp(k) = d_zero
            else
              aimp(k) = -(rhoxfl(k)*rrhoxhl(k))*kth(k)*dt*rdzq(k)*rdza(k)
            end if
            if ( k == kz ) then
              cimp(k) = d_zero
            else
              cimp(k) = -(rhoxfl(k+1)*rrhoxhl(k))* &
                         kth(k+1) * dt*rdzq(k)*rdza(k+1)
            end if
            bimp(k) = d_one - aimp(k) - cimp(k)
            rimp1(k) = qix(k)
          end do
          call solve_tridiag(aimp,bimp,cimp,rimp1,qix,kz)
        end if
!
!       !Re-update N^2 at the surface; this requires recalculation
!       !of the surface momentum fluxes
!
!       !Calculate surface momentum fluxes
!       uflxp = -uvdragx*ux(kz)/rhoxsf
!       vflxp = -uvdragx*vx(kz)/rhoxsf
!       ustxsq = sqrt(uflxp*uflxp + vflxp*vflxp)
!
!       !Estimate of surface eddy diffusivity, for estimating the
!       !surface N^2 from the surface virtual heat flux
!       kh0 = vonkar*2*sqrt(max(tkemin,3.25 * ustxsq))
!
!       !Estimate the surface N^2 from the surface virtual heat flux
!       nsquar(kzp1) = -egrav/thgb*thvflx/kh0


        !*************************************************************
        !****** Implicit Diffusion of tracers ************************
        !*************************************************************

        if ( ichem == 1 ) then
          xkchi = kth

          !Set the tridiagonal coefficients that apply to all of the tracers
          do k = 1 , kz
            if ( k == 1 ) then
              aimp(k) = d_zero
            else
              aimp(k) = -(rhoxfl(k)*rrhoxhl(k))*     &
                          xkchi(k) * dt *rdzq(k)*rdza(k)
            end if
            if ( k == kz ) then
              cimp(k) = d_zero
            else
              cimp(k) = -(rhoxfl(k+1)*rrhoxhl(k))*   &
                          xkchi(k+1) * dt *rdzq(k)*rdza(k+1)
            end if
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
            rimp1(kz) = chix(kz,itr) + dt*chifxx(itr)*rrhoxhl(kz)*rdzq(kz)
            ! no flux out top, so no (k == 1)
            rimp1(1) = d_zero

            !Run the tridiagonal solver
            call solve_tridiag(aimp,bimp,cimp,rimp1,uimp1,kz)

            !Get the chi value implied for the next timestep
            do k = 1 , kz
              chix(k,itr) = uimp1(k)
            end do
          end do
        end if !End tracer diffusion

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
        do ilay = 1 , kpbconv
          k = ktop(ilay)
          if ( qcx(k) > minqq .and. k > 1 ) then
            bouyan(k) = bouyan(k) - rttenx(k)*(presfl(k+1)-presfl(k)) * &
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
                          kethl(k-1)*dt*rdzq(k-1)*rdza(k)
          end if
          if ( k == kz ) then
            cimp(k-1) = d_zero
            ! account for implicit part of flux between level kz and surface
            bimp(k-1) = d_one - aimp(k-1) - cimp(k-1) + dt *    &
                        ( sqrt(tke(k))*rczero/bbls(k) +        &
                        (rhoxhl(k)*rrhoxfl(k))*kethl(k)*rdzq(k)*rdza(k) )
          else
            cimp(k-1) = -(rhoxhl(k)*rrhoxfl(k))*    &
                          kethl(k)*dt*rdzq(k)*rdza(k)
            bimp(k-1) = d_one - aimp(k-1) - cimp(k-1) + dt *   &
                        sqrt(tke(k))*rczero/max(bbls(k),bbls(k+1))
          end if
          ! now find right side
          if ( k == kz ) then
            ! account for fixed part of flux between level kz and surface
            rimp1(k-1) = tke(k) + dt * ( shear(k)+bouyan(k) +   &
                           tke(kzp1)*(rhoxhl(k)*rrhoxfl(k))*    &
                           kethl(k)*rdzq(k)*rdza(k) )
          else
            rimp1(k-1) = tke(k) + dt * (shear(k)+bouyan(k))
          end if
        end do imptkeflux
        call solve_tridiag(aimp,bimp,cimp,rimp1,uimp1,kzm1)
        ! update the tke
        do  k = 2 , kz
          tke(k) = max(uimp1(k-1),tkemin) ! background tke .001
        end do

!*******************************************************************************
!*******************************************************************************
!**************** Calculation of Tendencies for Model Output *******************
!*******************************************************************************
!*******************************************************************************

        ! Calculate the TCM tendencies for the model's prognostic variables
        ! For everything but TKE, couple the tendency (multiply by pstar)
        ibnd = 1
        if ( update_only_in_pbl ) ibnd = kpbl2dx
        do k = ibnd , kz
          ! Zonal wind tendency
          p2m%uuwten(j,i,k) = psbx*(ux(k)-uxs(k))*rdt
          ! Meridional wind tendency
          p2m%vuwten(j,i,k) = psbx*(vx(k)-vxs(k))*rdt
          ! TKE tendency
          p2m%tkeuwten(j,i,k) = (tke(k)-tkes(k))*rdt
          ! Temperature tendency
          p2m%tuwten(j,i,k) = psbx*(thx(k)-thxs(k))*exnerhl(k)*rdt
          ! Water vapor tendency
          p2m%qxuwten(j,i,k,iqv) = psbx*(qx(k)-qxs(k))*rdt
          ! Cloud water tendency
          p2m%qxuwten(j,i,k,iqc) = psbx*(qcx(k)-qcxs(k))*rdt
          ! Momentum diffusivity
          uwstateb%kzm(j,i,k) = kzm(k)
          ! Scalar diffusivity
          uwstateb%kth(j,i,k) = kth(k)
        end do

        if ( ipptls == 2 ) then
          do k = ibnd , kz
            p2m%qxuwten(j,i,k,iqi) = psbx*(qix(k)-qixs(k))*rdt
          end do
        end if

        if ( ichem == 1 ) then
          do itr = 1 , ntr
            do k = ibnd , kz
              chiuwten(j,i,k,itr) = psbx*(chix(k,itr)-chixs(k,itr))*rdt
            end do
          end do
        end if

        ! Output the diagnosed TKE
        uwstatea%srftke(j,i) = tke(kzp1)
        uwstateb%srftke(j,i) = tke(kzp1)

        p2m%kpbl(j,i) = kpbl2dx
        p2m%zpbl(j,i) = pblx
      end do
    end do
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif

    contains

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
      real(rk8) , dimension(n) , intent(in) :: a , b , c , v
      real(rk8) , dimension(n) , intent(out) :: x
      real(rk8) , dimension(n) :: bp , vp
      real(rk8) :: m
      integer(ik4) :: i

      bp(1) = b(1)
      vp(1) = v(1)
      ! The first pass (setting coefficients):
      do i = 2 , n
        m = a(i)/bp(i-1)
        bp(i) = b(i) - m*c(i-1)
        vp(i) = v(i) - m*vp(i-1)
      end do
      x(n) = vp(n)/bp(n)
      ! The second pass (back-substition)
      do i = n-1, 1, -1
        x(i) = (vp(i) - c(i)*x(i+1))/bp(i)
      end do
    end subroutine solve_tridiag

    subroutine n2(thlxin,qwxin,ktmax)
      implicit none
      integer(ik4) , intent(in) :: ktmax
      real(rk8) , intent(in) , dimension(ktmax) :: thlxin , qwxin
      ! local variables
      real(rk8) :: tempv , tvbl , rcld , tvab , thvxfl , dtvdz
      integer(ik4) :: k

      do k = 2 , ktmax
        ! buoyancy is jump in thetav across flux level/dza
        ! first, layer below, go up and see if anything condenses.
        templ = thlxin(k)*exnerfl(k)
        rvls = pfesat(templ)*epop(k)
        temps = templ + (qwxin(k)-rvls)/(cpowlhv +    &
                         ep2*wlhv*rvls/(rgas*templ*templ))
        rvls = pfesat(temps)*epop(k)
        rcldb(k) = max(qwxin(k)-rvls,d_zero)
        tempv = (templ + wlhvocp*rcldb(k)) *    &
                (d_one + ep1*(qwxin(k)-rcldb(k)) - rcldb(k))
        tvbl = tempv*rexnerfl(k)
        ! now do layer above; go down to see how much evaporates
        templ = thlxin(k-1)*exnerfl(k)
        rvls = pfesat(templ)*epop(k)
        temps = templ+(qwxin(k-1)-rvls) / &
                       (cpowlhv+ep2*wlhv*rvls/(rgas*templ*templ))
        rvls = pfesat(temps)*epop(k)
        rcld = max(qwxin(k-1)-rvls,d_zero)
        tempv = (templ + wlhvocp*rcld) *    &
                (d_one + ep1*(qwxin(k-1)-rcld) - rcld)
        tvab = tempv*rexnerfl(k)

        thvxfl = d_half * (tvab+tvbl)
        dtvdz = (tvab - tvbl) *rdza(k)
        nsquar(k) = egrav/thvxfl * dtvdz
      end do
      nsquar(1) = nsquar(2)
    end subroutine n2

    subroutine melloryamada(ktmax,kpbconv)
      implicit none
      integer(ik4) , intent(in) :: ktmax , kpbconv
      real(rk8) :: gh , a1ob1 , delthvl , elambda , bige , biga
      real(rk8) , parameter :: a1 = 0.92D0 , b1 = 16.6D0 , c1 = 0.08D0 , &
                               a2 = 0.74D0 , b2 = 10.1D0
      integer(ik4) :: k , ilay
      real(rk8) , parameter :: kthmax = 1.0D3

      a1ob1 = a1/b1

      ! calculate the diffusivities for momentum, thetal, qw and tke
      ! kth and kzm are at full levels.
      ! kethl is at half levels.

      kzm(ktmax) = d_zero
      kth(ktmax) = d_zero
      kzm(1) = d_zero
      kth(1) = d_zero
      sm(ktmax) = d_one
      sh(ktmax) = d_one
      sm(1) = d_one
      sh(1) = d_one

      kloop: &
      do k = ktmax - 1, 2, -1
        gh = -bbls(k)*bbls(k)*nsquar(k)/(d_two*tke(k)+1.0D-9)
        ! gh = min(gh,0.0233D0)
        ! TAO: Added the -0.28 minimum for the G function, as stated
        ! in Galperin (1988), eq 30
        gh = max(min(gh,0.0233D0),-0.28D0)

        sm(k) = a1 * (d_one - d_three*c1 - d_six*a1ob1 - d_three*a2*gh*   &
                     ((b2-d_three*a2)*(d_one - d_six*a1ob1) -             &
                       d_three*c1 * (b2 + d_six*a1))) /                   &
                     ((d_one - d_three*a2*gh * (d_six*a1 + b2)) *         &
                      (d_one - d_nine*a1*a2*gh))
        sh(k) = a2 * (d_one-d_six*a1ob1) / (d_one-d_three*a2*gh*(d_six*a1+b2))

        ! kzm(k) = min(bbls(k)*sqrt(2*tke(k))*sm(k),10000.0D0)
        ! kth(k) = min(bbls(k)*sqrt(2*tke(k))*sh(k),10000.0D0)

        ! Limit the diffusivity to be the vertical grid spacing squared
        ! over the time step; this implies that the entrainment rate
        ! can only be so large that the BL height would change by
        ! one grid level over one time step -- TAO
        ! kthmax = max(min((zax(k-1)-zax(k))**2/dt,1.0D4),1.0D3)
        ! kthmax = 1.0D4

        ! Calculate the diffusion coefficients
        kzm(k) = min(bbls(k)*sqrt(d_two*tke(k))*sm(k),kthmax)
        kth(k) = min(bbls(k)*sqrt(d_two*tke(k))*sh(k),kthmax)
        ! Smoothly limit kth to a maximum value
        !kth(k) = d_two/mathpi*kthmax*atan(kth(k)/kthmax)
        !kzm(k) = d_two/mathpi*kthmax*atan(kzm(k)/kthmax)
        kethl(k) = nuk*sqrt(kzm(k)*kzm(k+1))
        kethl(k) = min(kethl(k),kthmax)
      end do kloop

      ! special case for tops of convective layers
      conv: &
      do ilay = 1 , kpbconv
        k = ktop(ilay)
        if ( nsquar(k) > minn2 ) then
          kethl(k) = nuk*kzm(k+1)
          if ( k >= 3 ) then
            kethl(k-1) = 0.0D0
            delthvl = (thlx(k-2)+thx(k-2)*ep1*qwx(k-2)) -   &
                      (thlx(k) + thx(k) * ep1 * qwx(k))
            elambda = wlhvocp*rcldb(k)*rexnerhl(k)/max(delthvl,0.1D0)
            bige = 0.8D0 * elambda
            biga = aone * (d_one + atwo * bige)

            ! kth(k) = min(10000.0D0, biga * sqrt(tke(k)**3)/nsquar(k)/    &
            !          max(bbls(k),bbls(k+1)))
            ! Limit the diffusivity to be the vertical grid spacing squared
            ! over the time step; this implies that the entrainment rate
            ! can only be so large that the BL height would change by
            ! one grid level over one time step -- TAO
            ! kthmax = max(min((zax(k-1)-zax(k))**2/dt,1.0D4),1.0D3)
            ! kthmax = 1.0D4
            kth(k) = min(kth(k), biga*sqrt(tke(k)**3)/nsquar(k)/ &
                     max(bbls(k),bbls(k+1)))
            ! kth(k) = biga * sqrt(tke(k)**3)/nsquar(k)/max(bbls(k),bbls(k+1))
            ! Smoothly limit kth to a maximum value
            !kth(k) = d_two/mathpi*kthmax*atan(kth(k)/kthmax)
            ! prandtl number from layer below
            ! kzm(k) = kth(k) / sh(k+1) * sm(k+1)
            kzm(k) = min(kzm(k),kth(k)/sh(k+1)*sm(k+1))
          end if
        end if
      end do conv

      ! need kethl at top
      kethl(1) = kethl(2)
      ! replace kethl at surface with something non-zero
      kethl(ktmax-1) = nuk*d_half*kzm(ktmax-1)
    end subroutine melloryamada

    subroutine pblhgt(ktmax,ktmin,kpbconv)
      implicit none
      integer(ik4) , intent(in) :: ktmax , ktmin
      integer(ik4) , intent(out) :: kpbconv
      integer(ik4) , dimension(kz) :: ktop_save
      integer(ik4) :: istabl , ibeg , ilay , nlev , k , itemp , kstart
      real(rk8) :: blinf , rnnll , tkeavg , trnnll , radnnll , delthvl , &
                  elambda , bige , biga , entnnll , tbbls

      ! find noncontiguous convectively unstable layers
      kpbconv = 0
      istabl = 1
      do k = 2 , ktmax
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
          ! bbls(k) = max(min(rstbl*sqrt(tke(k)/nsquar(k)), &
          !           vonkar*zqx(k)),1.0D-8)
        end if
      end do

      ! now see if they have sufficient buoyant convection to connect
      ktop_save = ktop
      ibeg = 1
      if ( kpbconv >= 1 ) then
        bigloop: &
        do
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
              tkeavg = tkeavg + tke(k) / nlev
            end do
            ! first extend up
            kstart = ktop(ilay) - 1
            searchup1: &
            do k = kstart , 2 , -1
              ! we always go up at least one, for the entrainment interface
              ktop(ilay) = k
              bbls(k) = min(vonkar * zqx(k),blinf)
              trnnll = nsquar(k)*bbls(k)*bbls(k)
              ! If this is the top of the layer, stop searching upward
              if ( trnnll*nlev >= -d_half*rnnll ) exit searchup1
              if ( ilay > 1 ) then
                ! did we merge with layer above?
                if ( ktop(ilay) == kbot(ilay-1) ) then
                  ibeg = ilay - 1
                  ktop(ibeg) = ktop_save(ibeg)+1
                  kbot(ibeg) = kbot(ibeg+1)
                  kpbconv = kpbconv - 1
                  do itemp = ibeg+1 , kpbconv
                    ktop(itemp) = ktop(itemp+1)
                    kbot(itemp) = kbot(itemp+1)
                    ktop_save(itemp) = ktop_save(itemp+1)
                  end do
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
            if ( qcx(k) > minqq ) then
              radnnll = rttenx(k)*(presfl(k+1)-presfl(k)) /  &
                        (rhoxfl(k)*uthvx(k)*exnerfl(k))
            end if
            entnnll = d_zero
            if ( k >= 3 ) then
              delthvl = (thlx(k-2)+thx(k-2)*ep1*qwx(k-2))   &
                      - (thlx(k) + thx(k)*ep1*qwx(k))
              elambda = wlhvocp*rcldb(k)*rexnerhl(k)/max(delthvl,0.1D0)
              bige = 0.8D0 * elambda
              biga = aone * (d_one + atwo * bige)
              entnnll = biga * sqrt(tkeavg**3) / bbls(k)
            end if
            if ( tkeavg > d_zero ) then
              rnnll = rnnll + min(d_zero,bbls(k)/sqrt(tkeavg)*(radnnll+entnnll))
            end if
            ! now extend down
            searchdown1: &
            do k = kbot(ilay)+1 , ktmax
              tbbls = min(vonkar * zqx(k),blinf)
              trnnll = nsquar(k)*tbbls*tbbls
              ! is it the bottom?
              if ( trnnll*nlev >= -d_half*rnnll ) cycle convlayerloop
              ! (skip the rest of this iteration of the do loop)
              kbot(ilay) = k
              if ( ilay < kpbconv .and. &
                   kbot(ilay) == ktop(ilay+1) ) then
                ! did we merge with layer below?
                ktop(ilay) = ktop(ilay)+1
                kbot(ilay) = kbot(ilay+1)
                kpbconv = kpbconv - 1
                shiftarray2: &
                do itemp = ilay+1 , kpbconv
                  ktop(itemp) = ktop(itemp+1)
                  kbot(itemp) = kbot(itemp+1)
                  ktop_save(itemp) = ktop_save(itemp+1)
                end do shiftarray2
                cycle bigloop ! recompute for the new, deeper layer
                              ! (restart the do loop)
              end if
              rnnll = rnnll + trnnll
              bbls(k) = tbbls
              nlev = nlev + 1
            end do searchdown1
          end do convlayerloop
          exit bigloop
        end do bigloop
      end if
      setbbls: &
      do ilay = 1 , kpbconv
        blinf = xfr*(zqx(ktop(ilay)-1) - zqx(kbot(ilay)+1))
        convkloop: &
        do k = ktop(ilay) , kbot(ilay)
          bbls(k) = min(vonkar * zqx(k),blinf)
        end do convkloop
      end do setbbls
      ! we should now have tops and bottoms for kpbconv layers
      if (kpbconv > 0 ) then
        if ( kbot(kpbconv) == ktmax ) then
          kmix2dx = ktop(kpbconv)
          if ( kpbl2dx >= 0 ) then
            if ( kpbconv > 1 ) then
              kpbl2dx = ktop(kpbconv-1)
            else
              kpbl2dx = kmix2dx
            end if
          else
            kpbl2dx = -kpbl2dx
          end if
        else
          kmix2dx = ktmax - 1
          if ( kpbl2dx >= 0 ) then
            kpbl2dx = ktop(kpbconv)
          else
            kpbl2dx = -kpbl2dx
          end if
        end if
      else
        kmix2dx = ktmax - 1
        if ( kpbl2dx >= 0 ) then
          kpbl2dx = kmix2dx
        else
          kpbl2dx = -kpbl2dx
        end if
      end if
      ! Limit the PBL height to be either the highest allowable (tropopause),
      ! or the minimum (ktmax-1; the first interface level above the surface).
      kpbl2dx = max(ktmin,kpbl2dx)
      kpbl2dx = min(ktmax-1,kpbl2dx)
      kmix2dx = max(ktmin,kmix2dx)
      kmix2dx = min(ktmax-1,kmix2dx)
      pblx = zqx(kmix2dx)
    end subroutine pblhgt

    ! Returns the saturation vapor pressure over water in units (cb)
    ! given the input pressure in cb and Temperature in K
    ! Modified from Buck (1981), J. App. Met. v 20
    function esatw(p,t)
      implicit none
      real(rk8) , intent(in) :: p , t
      real(rk8) :: esatw
      real(rk8) :: dum , arg , tdum
      ! Limit T to reasonable values.  I believe that this is necessary because
      ! in the iterative calculation of T and QV from the liquid water
      ! temperature, the temperature can take on some crazy values before it
      ! converges. --TAO 01/05/2011
      tdum = max(100.0D0,min(399.99D0,t))
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
      real(rk8) , intent(in) :: p , t
      real(rk8) :: esati
      real(rk8) :: dum , arg , tdum
      ! Limit T to reasonable values.  I believe that this is necessary because
      ! in the iterative calculation of T and QV from the liquid water
      ! temperature, the temperature can take on some crazy values before it
      ! converges. --TAO 01/05/2011
      tdum = max(100.0D0,min(399.99D0,t))
      dum = 1.0003D0 + 4.18D-5*p
      ! arg = 22.452D0*(t-tzero)/(t-0.6D0)
      arg = 22.452D0*(tdum-tzero)/(tdum-0.6D0)
      esati = dum*0.61115D0*dexp(arg)
    end function esati

    subroutine pblhgt_tao(ktmax,kpbconv)
      implicit none
      integer(ik4) , intent(in) :: ktmax
      integer(ik4) , intent(out) :: kpbconv
      real(rk8) , dimension(ktmax) :: ktimesz
      real(rk8) :: lambda
      logical , dimension(ktmax) :: isSaturated1D , isStable1D , isBelow7001D
      logical :: foundlayer
      integer(ik4) :: k
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
      where ( presfl >= 70000.0D0 )
        isBelow7001D = .true.
      end where
      ! First see if there is a cloud-topped boundary layer: its top will be
      ! stable, saturated, and below 700 mb
      do k = ktmax-1 , 1 , -1
        if ( isSaturated1D(k) .and. isStable1D(k) .and. isBelow7001D(k) ) then
          kmix2dx = k
          foundlayer = .true.
          exit
        end if
      end do
      ! If we didn't find a cloud-topped boundary layer, then find the first
      ! layer where the richardson number exceeds its threshold
      if ( .not. foundlayer ) then
        do k = ktmax-1 , 1 , -1
          if ( isStable1D(k) ) then
            kmix2dx = k
            foundlayer = .true.
            exit
          end if
        end do
      end if
      ! If we still didn't find a cloud-topped boundary layer, then
      ! set the top to be the first interface layer
      if ( .not. foundlayer ) then
        kmix2dx = ktmax - 1
      end if
      ! Set the boundary layer top and the top of the convective layer
      kmix2dx = max(kmix2dx,3)
      kmix2dx = min(kmix2dx,ktmax-1)
      kpbl2dx = kmix2dx
      ! Set that there is only one convective layer
      kpbconv = 1
      ktop(kpbconv) = kmix2dx
      ! Set the boundary layer height
      pblx = zqx(ktop(kpbconv))
      ! Smoothly interpolate the boundary layer height
      if ( (richnum(kmix2dx) >= rcrit) .and. &
           (richnum(kmix2dx+1) < rcrit) ) then
        pblx = zqx(kmix2dx+1) + (zqx(kmix2dx) - zqx(kmix2dx+1)) * &
               ( (rcrit - richnum(kmix2dx + 1)) /                 &
                 (richnum(kmix2dx) - richnum(kmix2dx + 1)) )
      end if
      ! Set the master length scale
      lambda = etal*pblx
      ktimesz = vonkar*zqx(1:ktmax)
      ! Within the boundary layer, the length scale is propor. to the height
      bbls = ktimesz/(d_one+ktimesz/lambda)
      ! Otherwise use a Stability-related length scale
      do k = kmix2dx-1 , 1 , -1
        if ( nsquar(k) > d_zero ) then
          bbls(k) = max(min(rstbl*sqrt(tke(k)/nsquar(k)),ktimesz(k)),1.0D-8)
        else
          bbls(k) = ktimesz(k) - ktimesz(k+1)
        end if
      end do
    end subroutine pblhgt_tao

  end subroutine uwtcm

end module mod_pbl_uwtcm
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
