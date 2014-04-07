module mod_clm_regcm
  use mod_date
  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_runparams
  use mod_mppparam
  use mod_mpmessage
  use mod_constants
  use mod_sunorbit
  use mod_regcm_types
  use mod_service
  use mod_clm_initialize
  use mod_clm_driver
  use mod_clm_varctl , only : use_c13 , co2_ppmv
  use mod_clm_varpar , only : nlevsoi
  use mod_clm_varcon , only : o2_molar_const , c13ratio , tfrz , tcrit
  use mod_clm_atmlnd , only : clm_a2l , clm_l2a , adomain
  use mod_clm_decomp , only : procinfo , get_proc_bounds

  private

  save

  public :: initclm45 , runclm45 , albedoclm45

  real(rk8) , dimension(:,:) , pointer :: rprec , rsnow

  contains

  subroutine initclm45(lm,lms)
    implicit none
    type(lm_exchange) , intent(inout) :: lm
    type(lm_state) , intent(inout) :: lms
    integer(ik4) :: i , j , n , begg , endg
    character(len=64) :: rdate

    call getmem2d(rprec,jci1,jci2,ici1,ici2,'initclm45:rprec')

    allocate(adomain%xlon(lndcomm%linear_npoint_sg(myid+1)))
    allocate(adomain%xlat(lndcomm%linear_npoint_sg(myid+1)))
    allocate(adomain%topo(lndcomm%linear_npoint_sg(myid+1)))

    call glb_c2l_ss(lndcomm,lm%xlat1,adomain%xlat)
    call glb_c2l_ss(lndcomm,lm%xlon1,adomain%xlon)
    call glb_c2l_ss(lndcomm,lm%ht1,adomain%topo)

    adomain%topo = adomain%topo*regrav

    call initialize1( )

    call get_proc_bounds(begg,endg)
    allocate(adomain%snow(begg:endg))
    allocate(adomain%tgrd(begg:endg))
    call glb_c2l_gs(lndcomm,lm%snowam,adomain%snow)
    call glb_c2l_gs(lndcomm,lm%tground2,adomain%tgrd)

    write(rdate,'(i10)') toint10(idatex)
    call initialize2(rdate)

    ! Compute simple LAND emissivity for RegCM radiation
    if ( iemiss == 1 ) then
      do n = 1 , nnsg
        do i = ici1 , ici2
          do j = jci1 , jci2
            if ( lm%ldmsk1(n,j,i) == 1 ) then
              if ( lm%iveg1(n,j,i) == 8 ) then
                lms%emisv(n,j,i) = 0.76D0
              else if ( lm%iveg1(n,j,i) == 11 ) then
                lms%emisv(n,j,i) = 0.85D0
              else if ( lm%iveg1(n,j,i) == 12 ) then
                lms%emisv(n,j,i) = 0.97D0
              else
                lms%emisv(n,j,i) = 0.9995D0
              end if
            end if
          end do
        end do
      end do
    else
      do n = 1 , nnsg
        do i = ici1 , ici2
          do j = jci1 , jci2
            if ( lm%ldmsk1(n,j,i) == 1 ) then
              lms%emisv(n,j,i) = 0.9995D0
            end if
          end do
        end do
      end do
    end if
    if ( ktau == 0 ) then
      lms%swdiralb = 0.16D0
      lms%swdifalb = 0.16D0
      lms%lwdiralb = 0.32D0
      lms%lwdifalb = 0.32D0
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
    real(rk8) :: caldayp1 , declinp1 , eccfp1
    logical :: doalb , rstwr , nlend
    type(rcm_time_interval) :: tdiff
    type(rcm_time_and_date) :: nextt
    character(len=64) :: rdate

    call atmosphere_to_land(lm)

    ! Compute NEXT calday
    tdiff = dtsrf
    nextt = idatex+tdiff
    caldayp1 = yeardayfrac(nextt)
    call orb_decl(caldayp1,eccen,mvelpp,lambm0,obliqr,declinp1,eccfp1)
    if ( ktau == 0 .or. mod(ktau+1,ntrad) == 0 ) then
      doalb = .true.
    else
      doalb = .false.
    end if
    rstwr = .false.
    nlend = .false.
    if ( ktau > 0 ) then
      if ( ktau+1 == mtau ) then
        nlend = .true.
        rstwr = .true.
        write(rdate,'(i10)') toint10(nextt)
      end if
      if ( (lfdomonth(idatex) .and. lmidnight(idatex)) ) then
        rstwr = .true.
        write(rdate,'(i10)') toint10(nextt)
      end if
    end if

    ! Run CLM
    call orb_decl(caldayp1,eccen,mvelpp,lambm0,obliqr,declinp1,eccfp1)

    call clm_drv(doalb, caldayp1, declinp1, declin, rstwr, nlend, rdate)

    call land_to_atmosphere(lms)

  end subroutine runclm45

  subroutine albedoclm45(lm,lms)
    implicit none
    type(lm_exchange) , intent(inout) :: lm
    type(lm_state) , intent(inout) :: lms
    real(rk8) , dimension(1:nnsg,jci1:jci2,ici1:ici2) :: lastgood
    ! Just get albedoes from clm_l2a
    clm_l2a%notused = clm_l2a%albd(:,1)
    lastgood = lms%swdiralb
    call glb_l2c_ss(lndcomm,clm_l2a%notused,lms%swdiralb)
    where ( lms%swdiralb > 0.9999D0 )
      lms%swdiralb = lastgood
    end where
    clm_l2a%notused = clm_l2a%albd(:,2)
    lastgood = lms%lwdiralb
    call glb_l2c_ss(lndcomm,clm_l2a%notused,lms%lwdiralb)
    where ( lms%lwdiralb > 0.9999D0 )
      lms%lwdiralb = lastgood
    end where
    clm_l2a%notused = clm_l2a%albi(:,1)
    lastgood = lms%swdifalb
    call glb_l2c_ss(lndcomm,clm_l2a%notused,lms%swdifalb)
    where ( lms%swdifalb > 0.9999D0 )
      lms%swdifalb = lastgood
    end where
    clm_l2a%notused = clm_l2a%albi(:,2)
    lastgood = lms%lwdifalb
    call glb_l2c_ss(lndcomm,clm_l2a%notused,lms%lwdifalb)
    where ( lms%lwdifalb > 0.9999D0 )
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
    real(rk8) :: hl , satvp

    rprec = (lm%cprate+lm%ncprate) * rtsrf

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
    call glb_c2l_gs(lndcomm,rprec,clm_a2l%forc_rain)

    call glb_c2l_gs(lndcomm,lm%swdir,clm_a2l%notused)
    clm_a2l%forc_solad(:,1) = clm_a2l%notused
    call glb_c2l_gs(lndcomm,lm%lwdir,clm_a2l%notused)
    clm_a2l%forc_solad(:,2) = clm_a2l%notused
    call glb_c2l_gs(lndcomm,lm%swdif,clm_a2l%notused)
    clm_a2l%forc_solai(:,1) = clm_a2l%notused
    call glb_c2l_gs(lndcomm,lm%lwdif,clm_a2l%notused)
    clm_a2l%forc_solai(:,2) = clm_a2l%notused

    ! Compute or alias
    clm_a2l%forc_pbot = clm_a2l%forc_pbot * d_1000 ! In Pa
    clm_a2l%forc_wind = sqrt(clm_a2l%forc_u**2 + clm_a2l%forc_v**2)
    clm_a2l%forc_q = clm_a2l%forc_q/(1.0D0+clm_a2l%forc_q)
    clm_a2l%forc_hgt_u = clm_a2l%forc_hgt
    clm_a2l%forc_hgt_t = clm_a2l%forc_hgt
    clm_a2l%forc_hgt_q = clm_a2l%forc_hgt
    clm_a2l%forc_psrf = (clm_a2l%forc_psrf+ptop)*d_1000
    clm_a2l%rainf = clm_a2l%forc_rain+clm_a2l%forc_snow
    do i = begg , endg
      hl = lh0 - lh1*(clm_a2l%forc_t(i)-tzero)
      satvp = lsvp1*dexp(lsvp2*hl*(d_one/tzero-d_one/clm_a2l%forc_t(i)))
      clm_a2l%forc_rh(i) = max(clm_a2l%forc_q(i) / &
              (ep2*satvp/(clm_a2l%forc_pbot(i)*0.01D0-satvp)),d_zero)
      clm_a2l%forc_vp(i) = satvp * clm_a2l%forc_rh(i)
      clm_a2l%forc_rh(i) = clm_a2l%forc_rh(i) * 100.0D0
      ! Set upper limit of air temperature for snowfall at 275.65K.
      ! This cut-off was selected based on Fig. 1, Plate 3-1, of Snow
      ! Hydrology (1956).
      if ( clm_a2l%forc_t(i) < tfrz + tcrit ) then
        clm_a2l%forc_snow(i) = clm_a2l%forc_rain(i)
        clm_a2l%forc_rain(i) = 0.0D0
      else
        clm_a2l%forc_snow(i) = 0.0D0
      end if
    end do

    if ( ichem /= 1 ) then
      clm_a2l%forc_pco2 = co2_ppmv*1.D-6*clm_a2l%forc_psrf
      clm_a2l%forc_ndep = 6.34D-5 ! ?
      if ( use_c13 ) then
        clm_a2l%forc_pc13o2 = c13ratio*clm_a2l%forc_pco2
      end if
      clm_a2l%forc_po2 = o2_molar_const*clm_a2l%forc_psrf
      do i = 1 , 14
        clm_a2l%forc_aer(:,i) = clm_a2l%forc_ndep ! ????
      end do
    else
      ! Species partial pressures ! Fabien?
      ! clm_a2l%forc_pco2   ! CO2 partial pressure (Pa)
      ! clm_a2l%forc_ndep   ! nitrogen deposition rate (gN/m2/s)
      ! clm_a2l%forc_pc13o2 ! C13O2 partial pressure (Pa)
      ! clm_a2l%forc_po2    ! O2 partial pressure (Pa)
      ! clm_a2l%forc_aer    ! aerosol deposition array
    end if

    if ( .true. ) then
      clm_a2l%forc_flood = 0.000D0
      clm_a2l%volr = 0.000D0
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

    call get_proc_bounds(begg,endg)

    ! Get back data from clm_l2a
    call glb_l2c_ss(lndcomm,clm_l2a%t_rad,lms%tgbb)

    call glb_l2c_ss(lndcomm,clm_l2a%t_ref2m,lms%t2m)
    call glb_l2c_ss(lndcomm,clm_l2a%q_ref2m,lms%q2m)

    call glb_l2c_ss(lndcomm,clm_l2a%emg,lms%emisv)

    ! CLM gives just wind speed, assume directions are same as input.
    clm_a2l%notused = atan(clm_a2l%forc_v/clm_a2l%forc_u)
    clm_l2a%notused = clm_l2a%u_ref10m*cos(clm_a2l%notused)
    call glb_l2c_ss(lndcomm,clm_l2a%notused,lms%u10m)
    clm_l2a%notused = clm_l2a%u_ref10m*sin(clm_a2l%notused)
    call glb_l2c_ss(lndcomm,clm_l2a%notused,lms%v10m)

    call glb_l2c_ss(lndcomm,clm_l2a%eflx_sh_tot,lms%sent)
    call glb_l2c_ss(lndcomm,clm_l2a%qflx_evap_tot,lms%evpr)
    clm_l2a%notused = sqrt(clm_l2a%taux**2+clm_l2a%tauy**2)
    call glb_l2c_ss(lndcomm,clm_l2a%notused,lms%drag)

    call glb_l2c_ss(lndcomm,clm_l2a%h2osno,lms%sncv)
    call glb_l2c_ss(lndcomm,clm_l2a%taux,lms%taux)
    call glb_l2c_ss(lndcomm,clm_l2a%tauy,lms%tauy)

    lms%tgrd = lms%tgbb
    lms%tgbrd = lms%tgbb
    lms%tlef = lms%t2m

    clm_a2l%notused = 0.0D0
    clm_l2a%notused = 0.0D0
    do k = 1 , nlevsoi
      do g = begg , endg
        if ( clm_l2a%soidpth(g,k) < 0.10D0 ) then
          clm_l2a%notused(g) = clm_l2a%notused(g) + clm_l2a%h2osoi_liq(g,k)
        else
          clm_a2l%notused(g) = clm_a2l%notused(g) + clm_l2a%h2osoi_liq(g,k)
        end if
      end do
    end do
    call glb_l2c_ss(lndcomm,clm_a2l%notused,lms%ssw)
    call glb_l2c_ss(lndcomm,clm_l2a%notused,lms%rsw)

    ! From the input
    call glb_l2c_ss(lndcomm,clm_a2l%forc_rain,lms%prcp)
    call glb_l2c_ss(lndcomm,clm_a2l%forc_psrf,lms%sfcp)

    ! Will fix
    !clm_l2a%eflx_lwrad_out
    !clm_l2a%fsa
    !clm_l2a%nee
    !clm_l2a%ram1
    !clm_l2a%h2osoi_vol
    !clm_l2a%rofliq
    !clm_l2a%rofice
    !clm_l2a%flxdst
    !clm_l2a%ddvel
    !clm_l2a%flxvoc
#ifdef LCH4
    !clm_l2a%flux_ch4
#endif

  end subroutine land_to_atmosphere

end module mod_clm_regcm
