module mod_clm_regcm
  use mod_date
  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_runparams
  use mod_ipcc_scenario , only : cgas , igh_co2
  use mod_memutil
  use mod_mppparam
  use mod_mpmessage
  use mod_constants
  use mod_sunorbit
  use mod_regcm_types
  use mod_service
  use mod_clm_initialize
  use mod_clm_driver
  use mod_clm_varctl , only : use_c13 , co2_ppmv , tcrit
  use mod_clm_varpar , only : nlevsoi
  use mod_clm_varcon , only : o2_molar_const , c13ratio , tfrz , &
                              denh2o , sb
  use mod_clm_atmlnd , only : clm_a2l , clm_l2a , adomain
  use mod_clm_decomp , only : procinfo , get_proc_bounds
  use mod_clm_megan
  private

  save

  public :: initclm45 , runclm45 , albedoclm45

  real(rk8) , dimension(:,:) , pointer :: temps

  contains

#include <pfesat.inc>
#include <pfqsat.inc>

  subroutine initclm45(lm,lms)
    implicit none
    type(lm_exchange) , intent(inout) :: lm
    type(lm_state) , intent(inout) :: lms
    integer(ik4) :: i , j , n , begg , endg , ilev
    character(len=64) :: rdate
    real(rkx) , pointer , dimension(:,:) :: p2
    real(rk8) , pointer , dimension(:) :: p1

    call getmem2d(temps,jci1,jci2,ici1,ici2,'initclm45:temps')

    allocate(adomain%xlon(lndcomm%linear_npoint_sg(myid+1)))
    allocate(adomain%xlat(lndcomm%linear_npoint_sg(myid+1)))
    allocate(adomain%topo(lndcomm%linear_npoint_sg(myid+1)))

    call glb_c2l_ss(lndcomm,lm%xlat1,adomain%xlat)
    call glb_c2l_ss(lndcomm,lm%xlon1,adomain%xlon)
    call glb_c2l_ss(lndcomm,lm%ht1,adomain%topo)
    adomain%topo = adomain%topo*regrav

    call initialize1( )

    call get_proc_bounds(begg,endg)
    deallocate(adomain%topo)
    deallocate(adomain%xlon)
    deallocate(adomain%xlat)
    allocate(adomain%snow(begg:endg))
    allocate(adomain%smoist(begg:endg))
    allocate(adomain%tgrd(begg:endg))
    allocate(adomain%luse(begg:endg))
    allocate(adomain%topo(begg:endg))
    allocate(adomain%xlon(begg:endg))
    allocate(adomain%xlat(begg:endg))
    allocate(adomain%rmoist(begg:endg,num_soil_layers))
    call glb_c2l_gs(lndcomm,lm%snowam,adomain%snow)
    call glb_c2l_gs(lndcomm,lm%smoist,adomain%smoist)
    call glb_c2l_gs(lndcomm,lm%tground2,adomain%tgrd)
    call glb_c2l_ss(lndcomm,lm%ht1,adomain%topo)
    call glb_c2l_ss(lndcomm,lm%xlat1,adomain%xlat)
    call glb_c2l_ss(lndcomm,lm%xlon1,adomain%xlon)
    adomain%topo = adomain%topo*regrav
    call glb_c2l_ss(lndcomm,lm%iveg1,adomain%luse)
    if ( replacemoist ) then
      do ilev = 1 , num_soil_layers
        call assignpnt(lm%rmoist,p2,ilev)
        call assignpnt(adomain%rmoist,p1,ilev)
        call glb_c2l_gs(lndcomm,p2,p1)
      end do
    end if

    write(rdate,'(i10)') toint10(idatex)
    call initialize2(rdate)

    ! If CLM45, the surface emissivity is not used.
    ! The CLM outputs directly to RegCM the radiant Temperature.
    ! We fill it here the output not to leave it empy, but it is not
    ! used in computing the surface Long Wave Radiation
    do n = 1 , nnsg
      do i = ici1 , ici2
        do j = jci1 , jci2
          if ( lm%ldmsk1(n,j,i) == 1 ) then
            lms%emisv(n,j,i) = lnd_sfcemiss
          end if
        end do
      end do
    end do
    if ( ktau == 0 ) then
      lms%swdiralb = 0.16_rkx
      lms%swdifalb = 0.16_rkx
      lms%lwdiralb = 0.32_rkx
      lms%lwdifalb = 0.32_rkx
      lms%swalb = (lms%swdiralb + lms%swdifalb)
      lms%lwalb = (lms%lwdiralb + lms%lwdifalb)
      do n = 1 , nnsg
        lms%tgbrd(n,:,:) = lm%tground2(:,:)
        lms%tgbb(n,:,:) = lm%tground2(:,:)
      end do
    end if
  end subroutine initclm45

  subroutine runclm45(lm,lms)
    implicit none
    type(lm_exchange) , intent(inout) :: lm
    type(lm_state) , intent(inout) :: lms
    real(rk8) :: caldayp1 , declinp1 , eccfp1 , declinp
    logical :: doalb , rstwr , nlend , nlomon
    type(rcm_time_interval) :: tdiff , triff
    type(rcm_time_and_date) :: nextt , nextr
    character(len=64) :: rdate

    call atmosphere_to_land(lm)

    ! Compute NEXT calday
    tdiff = dtsrf
    triff = dtsec
    nextt = idatex+tdiff
    nextr = idatex+triff
    declinp = declin
    caldayp1 = yeardayfrac(nextt)
    call orb_decl(real(yearpoint(nextt),rk8),eccen,mvelpp, &
                  lambm0,obliqr,declinp1,eccfp1)
    if ( ktau == 0 .or. mod(ktau+1,ntrad) == 0 ) then
      doalb = .true.
    else
      doalb = .false.
    end if
    rstwr = .false.
    nlend = .false.
    nlomon = .false.
    if ( ktau > 0 ) then
      ! Final timestep
      if ( ktau+1 == mtau ) then
        rstwr = .true.
        nlend = .true.
        if ( (lfdomonth(nextr) .and. lmidnight(nextr)) ) then
          nlomon = .true.
        end if
      else
        if ( ksav > 0 ) then
          if ( mod(ktau+1,ksav) == 0 ) then
            rstwr = .true.
            if ( (lfdomonth(nextr) .and. lmidnight(nextr)) ) then
              nlomon = .true.
            end if
          end if
        else
          if ( ksav == 0 ) then
            if ( (lfdomonth(nextr) .and. lmidnight(nextr)) ) then
              ! End of the month
              nlomon = .true.
              rstwr = .true.
            end if
          else
            if ( mod(ktau+1,-ksav) == 0 .or. &
                 (lfdomonth(nextr) .and. lmidnight(nextr)) ) then
              rstwr = .true.
              if ( (lfdomonth(nextr) .and. lmidnight(nextr)) ) then
                nlomon = .true.
              end if
            end if
          end if
        end if
      end if
      if ( rstwr .and. myid == italk ) then
        write(rdate,'(i10)') toint10(nextr)
        write (stdout,*) 'Write restart file for CLM at ', tochar(nextr)
      end if
    end if

    ! Run CLM
    call clm_drv(doalb,caldayp1,declinp1,declinp,rstwr,nlend,nlomon,rdate)
    call land_to_atmosphere(lms)

  end subroutine runclm45

  subroutine albedoclm45(lm,lms)
    implicit none
    type(lm_exchange) , intent(inout) :: lm
    type(lm_state) , intent(inout) :: lms
    real(rk4) , dimension(1:nnsg,jci1:jci2,ici1:ici2) :: lastgood
    ! Just get albedoes from clm_l2a
    clm_l2a%notused = clm_l2a%albd(:,1)
    lastgood = lms%swdiralb
    call glb_l2c_ss(lndcomm,clm_l2a%notused,lms%swdiralb)
    where ( lms%swdiralb > 0.9999_rkx )
      lms%swdiralb = lastgood
    end where
    clm_l2a%notused = clm_l2a%albd(:,2)
    lastgood = lms%lwdiralb
    call glb_l2c_ss(lndcomm,clm_l2a%notused,lms%lwdiralb)
    where ( lms%lwdiralb > 0.9999_rkx )
      lms%lwdiralb = lastgood
    end where
    clm_l2a%notused = clm_l2a%albi(:,1)
    lastgood = lms%swdifalb
    call glb_l2c_ss(lndcomm,clm_l2a%notused,lms%swdifalb)
    where ( lms%swdifalb > 0.9999_rkx )
      lms%swdifalb = lastgood
    end where
    clm_l2a%notused = clm_l2a%albi(:,2)
    lastgood = lms%lwdifalb
    call glb_l2c_ss(lndcomm,clm_l2a%notused,lms%lwdifalb)
    where ( lms%lwdifalb > 0.9999_rkx )
      lms%lwdifalb = lastgood
    end where
    ! This should be the vegetation albedo!
    lms%swalb = lms%swdiralb+lms%swdifalb
    lms%lwalb = lms%lwdiralb+lms%lwdifalb
  end subroutine albedoclm45

  subroutine atmosphere_to_land(lm)
    implicit none
    type(lm_exchange) , intent(inout) :: lm
    integer(ik4) :: begg , endg , i
    !real(rk8) :: fsnts
    !real(rk8) , parameter :: ax = -48.23_rk8
    !real(rk8) , parameter :: bx = 0.75_rk8
    !real(rk8) , parameter :: cx = 1.16_rk8
    !real(rk8) , parameter :: dx = 1.02_rk8

    real(rkx) :: satq , satp

    call get_proc_bounds(begg,endg)

    ! Fill clm_a2l
    call glb_c2l_gs(lndcomm,lm%tatm,clm_a2l%forc_t)
    call glb_c2l_gs(lndcomm,lm%uatm,clm_a2l%forc_u)
    call glb_c2l_gs(lndcomm,lm%vatm,clm_a2l%forc_v)
    call glb_c2l_gs(lndcomm,lm%qvatm,clm_a2l%forc_q)
    call glb_c2l_gs(lndcomm,lm%hgt,clm_a2l%forc_hgt)
    call glb_c2l_gs(lndcomm,lm%patm,clm_a2l%forc_pbot)
    call glb_c2l_gs(lndcomm,lm%thatm,clm_a2l%forc_th)
    call glb_c2l_gs(lndcomm,lm%rhox,clm_a2l%forc_rho)
    call glb_c2l_gs(lndcomm,lm%sfps,clm_a2l%forc_psrf)
    call glb_c2l_gs(lndcomm,lm%dwrlwf,clm_a2l%forc_lwrad)
    call glb_c2l_gs(lndcomm,lm%solar,clm_a2l%forc_solar)
    temps = (lm%cprate+lm%ncprate) * rtsrf
    call glb_c2l_gs(lndcomm,temps,clm_a2l%forc_rain)

    call glb_c2l_gs(lndcomm,lm%swdir,clm_a2l%notused)
    clm_a2l%forc_solad(:,1) = clm_a2l%notused
    call glb_c2l_gs(lndcomm,lm%lwdir,clm_a2l%notused)
    clm_a2l%forc_solad(:,2) = clm_a2l%notused
    call glb_c2l_gs(lndcomm,lm%swdif,clm_a2l%notused)
    clm_a2l%forc_solai(:,1) = clm_a2l%notused
    call glb_c2l_gs(lndcomm,lm%lwdif,clm_a2l%notused)
    clm_a2l%forc_solai(:,2) = clm_a2l%notused

    ! Compute or alias
    clm_a2l%forc_wind = sqrt(clm_a2l%forc_u**2 + clm_a2l%forc_v**2)
    !clm_a2l%forc_th = clm_a2l%forc_t*(1.0e5_rk8/clm_a2l%forc_pbot)**rovcp
    clm_a2l%forc_hgt_u = clm_a2l%forc_hgt
    clm_a2l%forc_hgt_t = clm_a2l%forc_hgt
    clm_a2l%forc_hgt_q = clm_a2l%forc_hgt

    do i = begg , endg
      satp = pfesat(real(clm_a2l%forc_t(i),rkx))
      satq = pfqsat(real(clm_a2l%forc_t(i),rkx), &
                    real(clm_a2l%forc_pbot(i),rkx),satp)
      clm_a2l%forc_rh(i) = min(max(clm_a2l%forc_q(i)/satq,real(rhmin,rk8)), &
                               real(rhmax,rk8))
      clm_a2l%forc_vp(i) = real(satp,rk8) * clm_a2l%forc_rh(i)
      clm_a2l%forc_rh(i) = clm_a2l%forc_rh(i) * 100.0_rk8
      if ( clm_a2l%forc_t(i)-tfrz < tcrit ) then
        clm_a2l%forc_snow(i) = clm_a2l%forc_rain(i)
        clm_a2l%forc_rain(i) = d_zero
        ! CLM4.5 does not allow for snow and rain together...
        !fsnts = ax * (tanh(bx*(clm_a2l%forc_t(i)-tzero-cx))-dx)
        !clm_a2l%forc_snow(i) = clm_a2l%forc_rain(i)*fsnts
        !clm_a2l%forc_rain(i) = clm_a2l%forc_rain(i)-clm_a2l%forc_snow(i)
      else
        clm_a2l%forc_snow(i) = d_zero
      end if
    end do
    ! Specific humidity
    clm_a2l%forc_q = clm_a2l%forc_q/(1.0_rk8+clm_a2l%forc_q)
    clm_a2l%rainf = clm_a2l%forc_rain+clm_a2l%forc_snow

    ! interface chemistry / surface

    if ( ichem /= 1 ) then
      clm_a2l%forc_pco2 = cgas(igh_co2,xyear)*1.e-6_rk8*clm_a2l%forc_psrf
      clm_a2l%forc_ndep = 6.34e-5_rk8
      if ( use_c13 ) then
        clm_a2l%forc_pc13o2 = c13ratio*clm_a2l%forc_pco2
      end if
      clm_a2l%forc_po2 = o2_molar_const*clm_a2l%forc_psrf
      ! deposition of aerosol == zero
      clm_a2l%forc_aer(:,:) = d_zero !
    else
      !
      ! interface with atmospheric chemistry
      ! CO2 partial pressure (Pa)
      clm_a2l%forc_pco2 = cgas(igh_co2,xyear)*1.e-6_rk8*clm_a2l%forc_psrf
      if ( use_c13 ) then
       ! C13O2 partial pressure (Pa)
       clm_a2l%forc_pc13o2 = c13ratio*clm_a2l%forc_pco2
      end if
      ! O2 partial pressure (Pa)
      clm_a2l%forc_po2 = o2_molar_const*clm_a2l%forc_psrf
      ! FAB: this is the interactive part
      ! a) snow ageing
      !    flux of species have to be consistent with snowhydrology use
      !    of indices unit is kg/m2/s
      ! b) drydeposition BC HL
      !    flux arriving through lm interface are accumulated between
      !    two surface call : needs to average with rtsrf
      ! c) dry deposition BC HL
      if ( ibchl > 0 ) then
        temps(:,:) = lm%drydepflx (jci1:jci2,ici1:ici2,ibchl) * rtsrf
        call glb_c2l_gs(lndcomm,temps,clm_a2l%notused)
        clm_a2l%forc_aer(:,1) = clm_a2l%notused
      end if
      ! drydeposition BCHB
      if ( ibchb > 0 ) then
        temps(:,:) = lm%drydepflx(jci1:jci2,ici1:ici2,ibchb) * rtsrf
        call glb_c2l_gs(lndcomm,temps,clm_a2l%notused)
        clm_a2l%forc_aer(:,2) = clm_a2l%notused
      end if
      ! wet dep BC (sum rainout and washout fluxes, sum hb amd hl)
      if ( ibchb > 0 .and. ibchl > 0 ) then
        temps(:,:) =  (lm%wetdepflx(jci1:jci2,ici1:ici2,ibchb)  &
                         +  lm%wetdepflx(jci1:jci2,ici1:ici2,ibchl)) * rtsrf
        call glb_c2l_gs(lndcomm,temps,clm_a2l%notused)
        clm_a2l%forc_aer(:,3) = clm_a2l%notused
      end if
      ! drydeposition OC HL
      if ( iochl > 0 ) then
        temps(:,:) = lm%drydepflx(jci1:jci2,ici1:ici2,iochl)*  rtsrf
        call glb_c2l_gs(lndcomm,temps,clm_a2l%notused)
        clm_a2l%forc_aer(:,4) = clm_a2l%notused
      end if
      ! drydeposition OC HB
      if ( iochb > 0 ) then
        temps(:,:) =  lm%drydepflx(jci1:jci2,ici1:ici2,iochb) *  rtsrf
        call glb_c2l_gs(lndcomm,temps,clm_a2l%notused)
        clm_a2l%forc_aer(:,5) = clm_a2l%notused
      end if
      ! wet dep OC (sum rainout and washout fluxes, sum hb and hl)
      if(iochb >0 .and. iochl >0 ) then
        temps(:,:) = (lm%wetdepflx(jci1:jci2,ici1:ici2,iochb)   &
                         +  lm%wetdepflx(jci1:jci2,ici1:ici2,iochl)) * rtsrf
        call glb_c2l_gs(lndcomm,temps,clm_a2l%notused)
        clm_a2l%forc_aer(:,6) = clm_a2l%notused
      end if
      if ( size(lm%idust) == 4 ) then
       ! wet dep dust 1
       temps(:,:) = (lm%wetdepflx(jci1:jci2,ici1:ici2,lm%idust(1)) &
                      +  lm%wetdepflx(jci1:jci2,ici1:ici2,lm%idust(1))) * rtsrf
       call glb_c2l_gs(lndcomm,temps,clm_a2l%notused)
       clm_a2l%forc_aer(:,7) = clm_a2l%notused
       ! dry dep dust 1
       temps(:,:) =  lm%drydepflx(jci1:jci2,ici1:ici2,lm%idust(1))  * rtsrf
       call glb_c2l_gs(lndcomm,temps,clm_a2l%notused)
       clm_a2l%forc_aer(:,8) = clm_a2l%notused

       ! wet dep dust 2
       temps(:,:) =(lm%wetdepflx(jci1:jci2,ici1:ici2,lm%idust(2)) &
                      +  lm%wetdepflx(jci1:jci2,ici1:ici2,lm%idust(2))) * rtsrf
       call glb_c2l_gs(lndcomm,temps,clm_a2l%notused)
       clm_a2l%forc_aer(:,9) = clm_a2l%notused
       ! dry dep dust 2
       temps(:,:) = lm%drydepflx(jci1:jci2,ici1:ici2,lm%idust(2))  * rtsrf
       call glb_c2l_gs(lndcomm,temps,clm_a2l%notused)
       clm_a2l%forc_aer(:,10) = clm_a2l%notused

       ! wet dep dust 3
       temps(:,:) = (lm%wetdepflx(jci1:jci2,ici1:ici2,lm%idust(3)) &
                      +  lm%wetdepflx(jci1:jci2,ici1:ici2,lm%idust(3))) * rtsrf
       call glb_c2l_gs(lndcomm,temps,clm_a2l%notused)
       clm_a2l%forc_aer(:,11) = clm_a2l%notused
       ! dry dep dust 3
       temps(:,:) =  lm%drydepflx(jci1:jci2,ici1:ici2,lm%idust(3))  * rtsrf
       call glb_c2l_gs(lndcomm,temps,clm_a2l%notused)
       clm_a2l%forc_aer(:,12) = clm_a2l%notused

       ! wet dep dust 4
       temps(:,:) = (lm%wetdepflx(jci1:jci2,ici1:ici2,lm%idust(4)) &
                      +  lm%wetdepflx(jci1:jci2,ici1:ici2,lm%idust(4))) * rtsrf
       call glb_c2l_gs(lndcomm,temps,clm_a2l%notused)
       clm_a2l%forc_aer(:,13) = clm_a2l%notused
       ! dry dep dust 4
       temps(:,:) = lm%drydepflx (jci1:jci2,ici1:ici2,lm%idust(4))  * rtsrf
       call glb_c2l_gs(lndcomm,temps,clm_a2l%notused)
       clm_a2l%forc_aer(:,14) = clm_a2l%notused
      end if
      ! FAB to do : treat the 12 bins case ..
      ! b : pass the nitrogen deposition flux
    end if ! end test on ichem

    if ( .true. ) then
      clm_a2l%forc_flood = 0.000_rk8
      clm_a2l%volr = 0.000_rk8
    else
      ! Runoff in input ? Chym ?
      ! clm_a2l%forc_flood  ! flood (mm/s)
      ! clm_a2l%volr        ! rof volr (m3)
    end if
  end subroutine atmosphere_to_land

  subroutine land_to_atmosphere(lms)
    implicit none
    type(lm_state) , intent(inout) :: lms
    integer(ik4) :: k , g , begg , endg
    real(rk8) , pointer , dimension(:,:,:) :: emis2d

    call get_proc_bounds(begg,endg)

    ! Get back data from clm_l2a
    call glb_l2c_ss(lndcomm,clm_l2a%t_rad,lms%tgbb)

    call glb_l2c_ss(lndcomm,clm_l2a%t_ref2m,lms%t2m)
    call glb_l2c_ss(lndcomm,clm_l2a%q_ref2m,lms%q2m)

    ! CLM gives just wind speed, assume directions are same as input.
    clm_a2l%notused = clm_l2a%u_ref10m/clm_a2l%forc_wind
    clm_l2a%notused = clm_a2l%forc_u*clm_a2l%notused
    call glb_l2c_ss(lndcomm,clm_l2a%notused,lms%u10m)
    clm_l2a%notused = clm_a2l%forc_v*clm_a2l%notused
    call glb_l2c_ss(lndcomm,clm_l2a%notused,lms%v10m)

    call glb_l2c_ss(lndcomm,clm_l2a%eflx_sh_tot,lms%sent)
    call glb_l2c_ss(lndcomm,clm_l2a%qflx_evap_tot,lms%evpr)
    clm_l2a%notused = sqrt(clm_l2a%taux**2+clm_l2a%tauy**2) / &
                      sqrt(clm_a2l%forc_u**2+clm_a2l%forc_v**2)
    call glb_l2c_ss(lndcomm,clm_l2a%notused,lms%drag)

    call glb_l2c_ss(lndcomm,clm_l2a%h2osno,lms%sncv)
    call glb_l2c_ss(lndcomm,clm_l2a%taux,lms%taux)
    call glb_l2c_ss(lndcomm,clm_l2a%tauy,lms%tauy)
    call glb_l2c_ss(lndcomm,clm_l2a%t_veg,lms%tlef)

    lms%tgrd = lms%tgbb
    lms%tgbrd = lms%tgbb

    clm_l2a%notused = 0.0_rk8
    do k = 1 , nlevsoi
      do g = begg , endg
        clm_l2a%notused(g) = max(clm_l2a%h2osoi_vol(g,k),0.0_rk8) * &
             max(clm_l2a%dzsoi(g,k),0.0_rk8)*denh2o
      end do
      call glb_l2c_ss(lndcomm,clm_l2a%notused,lms%tsw)
      lms%sw(:,:,:,k) = lms%tsw(:,:,:)
    end do
    lms%tsw(:,:,:) = sum(lms%sw,4)

    call glb_l2c_ss(lndcomm,clm_l2a%qflx_surf,lms%srnof)
    call glb_l2c_ss(lndcomm,clm_l2a%qflx_tot,lms%trnof)
    call glb_l2c_ss(lndcomm,clm_l2a%qflx_snomelt,lms%snwm)
    lms%snwm = lms%snwm * dtsrf

    ! From the input
    call glb_l2c_ss(lndcomm,clm_a2l%rainf,lms%prcp)
    call glb_l2c_ss(lndcomm,clm_a2l%forc_psrf,lms%sfcp)

    !--------------------------------------------------
    ! From land to chemistry
    ! only for Isoprene , and in kg/m^2/sec
    ! FAB add compatibility for other biogenic species /
    ! TO BE UPDATED IF THE CHEM MECHANISM CHANGES (e.g add limonene, pinene etc)
    ! passed to the chemistry scheme for the right  mechanism tracer index
    ! use temporary table emis2d for calling glb_l2c_ss

    if ( ichem == 1 .and. enable_megan_emission ) then
      allocate(emis2d(1:nnsg,jci1:jci2,ici1:ici2))
      emis2d = 0.0_rk8
      do k = 1 , shr_megan_mechcomps_n
        if (shr_megan_mechcomps(k)%name == 'ISOP' .and. iisop > 0) then
          clm_l2a%notused(:) = clm_l2a%flxvoc(:,k)
          call glb_l2c_ss(lndcomm, clm_l2a%notused, emis2d)
          lms%vocemiss(:,:,:,iisop) = real(emis2d,rk4)
        end if
        ! add compatibility for other biogenic species !! /
        !
      end do
      deallocate(emis2d)
    end if

    ! pass the CLM dust flux to regcm
    ! nb: CLM considers 4 emission bins roughly equivalent to the regcm
    ! 4 bin version.
    ! if use the regcm 12 bin, the total mass is redistributed
    ! (chemlib/mod_che_dust)
    if ( ichem == 1 .and. ichdustemd == 3 ) then
      allocate(emis2d(1:nnsg,jci1:jci2,ici1:ici2))
      emis2d = 0.0_rk8
      do k = 1 , 4
        clm_l2a%notused(:) = clm_l2a%flxdst(:,k)
        call glb_l2c_ss(lndcomm, clm_l2a%notused, emis2d)
        lms%dustemiss(:,:,:,k) = real(emis2d,rk4)
      end do
      deallocate(emis2d)
    end if

    !--------------------------------------------------

    ! Will fix
    !clm_l2a%eflx_lwrad_out
    !clm_l2a%fsa
    !clm_l2a%nee
    !clm_l2a%ram1
    !clm_l2a%rofliq
    !clm_l2a%rofice
    !clm_l2a%flxdst
    !clm_l2a%ddvel
#ifdef LCH4
    !clm_l2a%flux_ch4
#endif

  end subroutine land_to_atmosphere

end module mod_clm_regcm
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
