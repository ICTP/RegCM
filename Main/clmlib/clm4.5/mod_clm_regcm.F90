module mod_clm_regcm
  use mod_date
  use mod_stdio
  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_kdinterp
  use mod_runparams
  use mod_ipcc_scenario, only : ghgval, igh_co2, igh_ch4
  use mod_memutil
  use mod_mppparam
  use mod_ensemble
  use mod_mpmessage
  use mod_constants
  use mod_sunorbit
  use mod_regcm_types
  use mod_service
  use mod_clm_initialize
  use mod_clm_driver
  use mod_clm_varctl, only : use_c13, co2_ppmv, tcrit, nextdate
  use mod_clm_varctl, only : ndep_nochem, luse_cru, crufile, lcru_rand
  use mod_clm_varpar, only : nlevsoi, nlevgrnd
  use mod_clm_varcon, only : o2_molar_const, c13ratio, tfrz, sb
  use mod_clm_atmlnd, only : clm_a2l, clm_l2a, adomain
  use mod_clm_decomp, only : procinfo, get_proc_bounds
  use mod_clm_megan
  use mod_clm_drydep, only : n_drydep
  use netcdf
  implicit none

  private

  save

  public :: initclm45, runclm45, albedoclm45
  public :: initsaclm45, runsaclm45

  real(rkx), dimension(:), pointer, contiguous :: glon, glat
  real(rkx), dimension(:), pointer, contiguous :: ihfac
  real(rkx), dimension(:,:), pointer, contiguous :: temps
  real(rkx), dimension(:,:), pointer, contiguous :: alon, alat
  real(rkx), dimension(:,:), pointer, contiguous :: rcp
  real(rkx), dimension(:,:), pointer, contiguous :: crupre, acp0, acp1, acp2
  real(rk8), dimension(:,:,:), pointer, contiguous :: emis2d

  type(h_interpolator) :: hint

  contains

  subroutine initsaclm45(lm)
    implicit none
    type(lm_exchange), intent(inout) :: lm
    integer(ik4) :: begg, endg, ilev
    character(len=64) :: rdate
    real(rkx), pointer, contiguous, dimension(:,:) :: p2
    real(rk8), pointer, contiguous, dimension(:) :: p1

    call getmem(temps,jci1,jci2,ici1,ici2,'initclm45:temps')
    if ( ichem == 1 ) then
      call getmem(emis2d,1,nnsg,jci1,jci2,ici1,ici2,'initclm45:emis2d')
    end if

    allocate(adomain%xlon(lndcomm%linear_npoint_sg(myid+1)))
    allocate(adomain%xlat(lndcomm%linear_npoint_sg(myid+1)))
    allocate(adomain%topo(lndcomm%linear_npoint_sg(myid+1)))
    allocate(adomain%area(lndcomm%linear_npoint_sg(myid+1)))

    call glb_c2l_ss(lndcomm,lm%xlat1,adomain%xlat)
    call glb_c2l_ss(lndcomm,lm%xlon1,adomain%xlon)
    call glb_c2l_ss(lndcomm,lm%ht1,adomain%topo)
    call glb_c2l_ss(lndcomm,lm%area1,adomain%area)
    adomain%topo = adomain%topo*regrav

    nextdate = rcmtimer%idate

    call initialize1( )

    call get_proc_bounds(begg,endg)

    deallocate(adomain%topo)
    deallocate(adomain%area)
    deallocate(adomain%xlon)
    deallocate(adomain%xlat)
    allocate(adomain%snow(begg:endg))
    allocate(adomain%smoist(begg:endg))
    allocate(adomain%tgrd(begg:endg))
    allocate(adomain%iveg(begg:endg))
    allocate(adomain%itex(begg:endg))
    allocate(adomain%topo(begg:endg))
    allocate(adomain%area(begg:endg))
    allocate(adomain%xlon(begg:endg))
    allocate(adomain%xlat(begg:endg))
    allocate(adomain%rmoist(begg:endg,num_soil_layers))
    allocate(adomain%rts(begg:endg,nlevgrnd))
    call glb_c2l_gs(lndcomm,lm%snowam,adomain%snow)
    call glb_c2l_gs(lndcomm,lm%smoist,adomain%smoist)
    call glb_c2l_gs(lndcomm,lm%tg,adomain%tgrd)
    call glb_c2l_ss(lndcomm,lm%ht1,adomain%topo)
    call glb_c2l_ss(lndcomm,lm%area1,adomain%area)
    call glb_c2l_ss(lndcomm,lm%xlat1,adomain%xlat)
    call glb_c2l_ss(lndcomm,lm%xlon1,adomain%xlon)
    adomain%topo = adomain%topo*regrav
    call glb_c2l_ss(lndcomm,lm%iveg1,adomain%iveg)
    call glb_c2l_ss(lndcomm,lm%itex1,adomain%itex)
    where ( adomain%itex == 14 )
      adomain%itex = 17
    end where
    if ( replacemoist ) then
      do ilev = 1, num_soil_layers
        call assignpnt(lm%rmoist,p2,ilev)
        call assignpnt(adomain%rmoist,p1,ilev)
        call glb_c2l_gs(lndcomm,p2,p1)
      end do
    end if
    if ( replacetemp ) then
      do ilev = 1, num_soil_layers
        call assignpnt(lm%rts,p2,ilev)
        call assignpnt(adomain%rts,p1,ilev)
        call glb_c2l_gs(lndcomm,p2,p1)
      end do
      if ( num_soil_layers < nlevgrnd ) then
        do ilev = num_soil_layers+1, nlevgrnd
          adomain%rts(:,ilev) = adomain%rts(:,num_soil_layers)
        end do
      end if
    end if

    write(rdate,'(a)') tochar10(rcmtimer%idate)
    call initialize2(rdate)

  end subroutine initsaclm45

  subroutine initclm45(lm,lms)
    implicit none
    type(lm_exchange), intent(inout) :: lm
    type(lm_state), intent(inout) :: lms
    integer(ik4) :: i, j, n, begg, endg, ilev
    character(len=64) :: rdate
    real(rkx), pointer, contiguous, dimension(:,:) :: p2
    real(rk8), pointer, contiguous, dimension(:) :: p1

    call getmem(temps,jci1,jci2,ici1,ici2,'initclm45:temps')
    if ( ichem == 1 ) then
      call getmem(emis2d,1,nnsg,jci1,jci2,ici1,ici2,'initclm45:emis2d')
    end if

    allocate(adomain%xlon(lndcomm%linear_npoint_sg(myid+1)))
    allocate(adomain%xlat(lndcomm%linear_npoint_sg(myid+1)))
    allocate(adomain%topo(lndcomm%linear_npoint_sg(myid+1)))
    allocate(adomain%area(lndcomm%linear_npoint_sg(myid+1)))

    call glb_c2l_ss(lndcomm,lm%xlat1,adomain%xlat)
    call glb_c2l_ss(lndcomm,lm%xlon1,adomain%xlon)
    call glb_c2l_ss(lndcomm,lm%ht1,adomain%topo)
    call glb_c2l_ss(lndcomm,lm%area1,adomain%area)
    adomain%topo = adomain%topo*regrav

    nextdate = rcmtimer%idate

    call initialize1( )

    call get_proc_bounds(begg,endg)

    deallocate(adomain%topo)
    deallocate(adomain%area)
    deallocate(adomain%xlon)
    deallocate(adomain%xlat)
    allocate(adomain%snow(begg:endg))
    allocate(adomain%smoist(begg:endg))
    allocate(adomain%tgrd(begg:endg))
    allocate(adomain%iveg(begg:endg))
    allocate(adomain%itex(begg:endg))
    allocate(adomain%topo(begg:endg))
    allocate(adomain%area(begg:endg))
    allocate(adomain%xlon(begg:endg))
    allocate(adomain%xlat(begg:endg))
    allocate(adomain%rmoist(begg:endg,num_soil_layers))
    allocate(adomain%rts(begg:endg,nlevgrnd))
    call glb_c2l_gs(lndcomm,lm%snowam,adomain%snow)
    call glb_c2l_gs(lndcomm,lm%smoist,adomain%smoist)
    call glb_c2l_gs(lndcomm,lm%tg,adomain%tgrd)
    call glb_c2l_ss(lndcomm,lm%ht1,adomain%topo)
    call glb_c2l_ss(lndcomm,lm%area1,adomain%area)
    call glb_c2l_ss(lndcomm,lm%xlat1,adomain%xlat)
    call glb_c2l_ss(lndcomm,lm%xlon1,adomain%xlon)
    adomain%topo = adomain%topo*regrav
    call glb_c2l_ss(lndcomm,lm%iveg1,adomain%iveg)
    call glb_c2l_ss(lndcomm,lm%itex1,adomain%itex)
    where ( adomain%itex == 14 )
      adomain%itex = 17
    end where
    if ( replacemoist ) then
      do ilev = 1, num_soil_layers
        call assignpnt(lm%rmoist,p2,ilev)
        call assignpnt(adomain%rmoist,p1,ilev)
        call glb_c2l_gs(lndcomm,p2,p1)
      end do
    end if
    if ( replacetemp ) then
      do ilev = 1, num_soil_layers
        call assignpnt(lm%rts,p2,ilev)
        call assignpnt(adomain%rts,p1,ilev)
        call glb_c2l_gs(lndcomm,p2,p1)
      end do
      if ( num_soil_layers < nlevgrnd ) then
        do ilev = num_soil_layers+1, nlevgrnd
          adomain%rts(:,ilev) = adomain%rts(:,num_soil_layers)
        end do
      end if
    end if

    write(rdate,'(a)') tochar10(rcmtimer%idate)
    call initialize2(rdate)

    if ( rcmtimer%start( ) ) then
      do i = ici1, ici2
        do j = jci1, jci2
          do n = 1, nnsg
            if ( lm%ldmsk1(n,j,i) == 1 ) then
              lms%swdiralb(n,j,i) = 0.16_rkx
              lms%swdifalb(n,j,i) = 0.16_rkx
              lms%lwdiralb(n,j,i) = 0.32_rkx
              lms%lwdifalb(n,j,i) = 0.32_rkx
              lms%swalb(n,j,i) = (lms%swdiralb(n,j,i) + lms%swdifalb(n,j,i))
              lms%lwalb(n,j,i) = (lms%lwdiralb(n,j,i) + lms%lwdifalb(n,j,i))
              lms%tgbrd(n,j,i) = lm%tg(j,i)
              lms%tgbb(n,j,i) = lm%tg(j,i)
            end if
          end do
        end do
      end do
    end if
  end subroutine initclm45

  subroutine runsaclm45(lm)
    implicit none
    type(lm_exchange), intent(inout) :: lm
    real(rk8) :: caldayp1, declinp1, eccfp1, declinp
    logical :: doalb, rstwr, nlend, nlomon
    type(rcm_time_interval) :: tdiff, triff
    type(rcm_time_and_date) :: nextt, nextr
    character(len=64) :: rdate

    call atmosphere_to_land(lm)

    ! Compute NEXT calday
    tdiff = dtsrf
    triff = dtsec
    nextt = rcmtimer%idate+tdiff
    nextr = rcmtimer%idate+triff
    declinp = declin
    caldayp1 = yeardayfrac(nextt)
    call orb_decl(real(yearpoint(nextt),rk8),eccen,mvelpp, &
                  lambm0,obliqr,declinp1,eccfp1)
    if ( rcmtimer%start( ) .or. syncro_rad%will_act() ) then
      if ( debug_level > 3 .and. myid == italk ) then
        write(stdout,*) 'Updating albedo at ',trim(rcmtimer%str())
      end if
      doalb = .true.
    else
      doalb = .false.
    end if
    rstwr = .false.
    nlend = .false.
    nlomon = .false.
    if ( rcmtimer%integrating( ) ) then
      ! Final timestep
      if ( rcmtimer%next_is_endtime ) then
        rstwr = .true.
        nlend = .true.
        if ( (lfdomonth(nextr) .and. lmidnight(nextr)) ) then
          nlomon = .true.
        end if
      else
        if ( associated(alarm_out_sav) ) then
          if ( savfrq > d_zero ) then
            if ( alarm_out_sav%will_act(dtsrf) ) then
              rstwr = .true.
              if ( (lfdomonth(nextr) .and. lmidnight(nextr)) ) then
                nlomon = .true.
              end if
            end if
          else
            if ( alarm_out_sav%will_act(dtsrf) .or. &
                 (lfdomonth(nextr) .and. lmidnight(nextr)) ) then
              rstwr = .true.
              if ( (lfdomonth(nextr) .and. lmidnight(nextr)) ) then
                nlomon = .true.
              end if
            end if
          end if
        else
          if ( (lfdomonth(nextr) .and. lmidnight(nextr)) ) then
            ! End of the month
            nlomon = .true.
            rstwr = .true.
          end if
        end if
      end if
      if ( rstwr .and. myid == italk ) then
        write(rdate,'(a)') tochar10(nextr)
        write (stdout,*) 'Write restart file for CLM at ', tochar(nextr)
      end if
    end if

    ! Run CLM
    call clm_drv(doalb,caldayp1,declinp1,declinp,rstwr,nlend,nlomon,rdate)

  end subroutine runsaclm45

  subroutine runclm45(lm,lms)
    !@acc use nvtx
    implicit none
    type(lm_exchange), intent(inout) :: lm
    type(lm_state), intent(inout) :: lms
    real(rk8) :: caldayp1, declinp1, eccfp1, declinp
    logical :: doalb, rstwr, nlend, nlomon
    type(rcm_time_interval) :: tdiff, triff
    type(rcm_time_and_date) :: nextt, nextr
    character(len=64) :: rdate
    !@acc call nvtxStartRange("runclm45")
    !@acc call nvtxStartRange("atmosphere_to_land")
    call atmosphere_to_land(lm)
    !@acc call nvtxEndRange
    ! Compute NEXT calday
    tdiff = dtsrf
    triff = dtsec
    nextt = rcmtimer%idate+tdiff
    nextr = rcmtimer%idate+triff
    declinp = declin
    caldayp1 = yeardayfrac(nextt)
    call orb_decl(real(yearpoint(nextt),rk8),eccen,mvelpp, &
                  lambm0,obliqr,declinp1,eccfp1)
    if ( rcmtimer%start( ) .or. syncro_rad%will_act() ) then
      if ( debug_level > 3 .and. myid == italk ) then
        write(stdout,*) 'Updating albedo at ',trim(rcmtimer%str())
      end if
      doalb = .true.
    else
      doalb = .false.
    end if
    rstwr = .false.
    nlend = .false.
    nlomon = .false.
    if ( rcmtimer%integrating( ) ) then
      ! Final timestep
      if ( rcmtimer%next_is_endtime ) then
        rstwr = .true.
        nlend = .true.
        if ( (lfdomonth(nextr) .and. lmidnight(nextr)) ) then
          nlomon = .true.
        end if
      else
        if ( associated(alarm_out_sav) ) then
          if ( savfrq > d_zero ) then
            if ( alarm_out_sav%will_act(dtsrf) ) then
              rstwr = .true.
              if ( (lfdomonth(nextr) .and. lmidnight(nextr)) ) then
                nlomon = .true.
              end if
            end if
          else
            if ( alarm_out_sav%will_act(dtsrf) .or. &
                 (lfdomonth(nextr) .and. lmidnight(nextr)) ) then
              rstwr = .true.
              if ( (lfdomonth(nextr) .and. lmidnight(nextr)) ) then
                nlomon = .true.
              end if
            end if
          end if
        else
          if ( (lfdomonth(nextr) .and. lmidnight(nextr)) ) then
            ! End of the month
            nlomon = .true.
            rstwr = .true.
          end if
        end if
      end if
      if ( rstwr .and. myid == italk ) then
        write(rdate,'(a)') tochar10(nextr)
        write (stdout,*) 'Write restart file for CLM at ', tochar(nextr)
      end if
    end if

    ! Run CLM
    call clm_drv(doalb,caldayp1,declinp1,declinp,rstwr,nlend,nlomon,rdate)
    !@acc call nvtxStartRange("land_to_atmosphere")
    call land_to_atmosphere(lm,lms)
    !@acc call nvtxEndRange
    !@acc call nvtxEndRange
  end subroutine runclm45

  subroutine albedoclm45(lm,lms)
    implicit none
    type(lm_exchange), intent(inout) :: lm
    type(lm_state), intent(inout) :: lms
    real(rkx), dimension(1:nnsg,jci1:jci2,ici1:ici2) :: lastgood
    integer(ik4) :: i, j, n
    ! Just get albedoes from clm_l2a
    !$acc kernels
    lastgood = lms%swdiralb
    clm_l2a%notused = clm_l2a%albd(:,1)
    !$acc end kernels
    call glb_l2c_ss(lndcomm,clm_l2a%notused,lms%swdiralb)
    do concurrent ( n = 1:nnsg, j = jci1:jci2, i = ici1:ici2 )
      if ( lm%ldmsk1(n,j,i) == 1 ) then
        if ( lms%swdiralb(n,j,i) > 0.9999_rkx ) then
          lms%swdiralb(n,j,i) = lastgood(n,j,i)
        end if
      end if
    end do
    !$acc kernels
    lastgood = lms%lwdiralb
    clm_l2a%notused = clm_l2a%albd(:,2)
    !$acc end kernels
    call glb_l2c_ss(lndcomm,clm_l2a%notused,lms%lwdiralb)
    do concurrent ( n = 1:nnsg, j = jci1:jci2, i = ici1:ici2 )
      if ( lm%ldmsk1(n,j,i) == 1 ) then
        if ( lms%lwdiralb(n,j,i) > 0.9999_rkx ) then
          lms%lwdiralb(n,j,i) = lastgood(n,j,i)
        end if
      end if
    end do
    !$acc kernels
    lastgood = lms%swdifalb
    clm_l2a%notused = clm_l2a%albi(:,1)
    !$acc end kernels
    call glb_l2c_ss(lndcomm,clm_l2a%notused,lms%swdifalb)
    do concurrent ( n = 1:nnsg, j = jci1:jci2, i = ici1:ici2 )
      if ( lm%ldmsk1(n,j,i) == 1 ) then
        if ( lms%swdifalb(n,j,i) > 0.9999_rkx ) then
          lms%swdifalb(n,j,i) = lastgood(n,j,i)
        end if
      end if
    end do
    !$acc kernels
    lastgood = lms%lwdifalb
    clm_l2a%notused = clm_l2a%albi(:,2)
    !$acc end kernels
    call glb_l2c_ss(lndcomm,clm_l2a%notused,lms%lwdifalb)
    do concurrent ( n = 1:nnsg, j = jci1:jci2, i = ici1:ici2 )
      if ( lm%ldmsk1(n,j,i) == 1 ) then
        if ( lms%lwdifalb(n,j,i) > 0.9999_rkx ) then
          lms%lwdifalb(n,j,i) = lastgood(n,j,i)
        end if
      end if
    end do
    ! This should be the vegetation albedo!
    do concurrent ( n = 1:nnsg, j = jci1:jci2, i = ici1:ici2 )
      if ( lm%ldmsk1(n,j,i) == 1 ) then
        lms%swalb(n,j,i) = max(lms%swdiralb(n,j,i),lms%swdifalb(n,j,i))
        lms%lwalb(n,j,i) = max(lms%lwdiralb(n,j,i),lms%lwdifalb(n,j,i))
      end if
    end do
  end subroutine albedoclm45

  subroutine atmosphere_to_land(lm)
    implicit none
    type(lm_exchange), intent(inout) :: lm
    integer(ik4) :: begg, endg, i, j, n
    real(rkx) :: satq, lat

    call get_proc_bounds(begg,endg)

    ! Fill clm_a2l
    call glb_c2l_gs(lndcomm,lm%tatm,clm_a2l%forc_t)
    call glb_c2l_gs(lndcomm,lm%uatm,clm_a2l%forc_u)
    call glb_c2l_gs(lndcomm,lm%vatm,clm_a2l%forc_v)
    call glb_c2l_gs(lndcomm,lm%qvatm,clm_a2l%forc_q)
    call glb_c2l_gs(lndcomm,lm%hgt,clm_a2l%forc_hgt)
    call glb_c2l_gs(lndcomm,lm%patm,clm_a2l%forc_pbot)
    call glb_c2l_gs(lndcomm,lm%thatm,clm_a2l%forc_th)
    call glb_c2l_gs(lndcomm,lm%sfps,clm_a2l%forc_psrf)
    call glb_c2l_gs(lndcomm,lm%dwrlwf,clm_a2l%forc_lwrad)
    call glb_c2l_gs(lndcomm,lm%solar,clm_a2l%forc_solar)
    !$acc kernels
    temps = (lm%cprate+lm%ncprate) * syncro_srf%rw
    !$acc end kernels
    if ( luse_cru ) then
      call read_cru_pre(lm,temps)
    end if
    call glb_c2l_gs(lndcomm,temps,clm_a2l%rainf)

    call glb_c2l_gs(lndcomm,lm%swdir,clm_a2l%notused)
    !$acc kernels
    clm_a2l%forc_solad(:,1) = clm_a2l%notused
    !$acc end kernels
    call glb_c2l_gs(lndcomm,lm%lwdir,clm_a2l%notused)
    !$acc kernels
    clm_a2l%forc_solad(:,2) = clm_a2l%notused
    !$acc end kernels
    call glb_c2l_gs(lndcomm,lm%swdif,clm_a2l%notused)
    !$acc kernels
    clm_a2l%forc_solai(:,1) = clm_a2l%notused
    !$acc end kernels
    call glb_c2l_gs(lndcomm,lm%lwdif,clm_a2l%notused)
    !$acc kernels
    clm_a2l%forc_solai(:,2) = clm_a2l%notused
    !$acc end kernels

    ! Compute or alias
    !these kernels cause illegal memory access
    !!!$acc kernels
    clm_a2l%forc_wind = max(sqrt(clm_a2l%forc_u**2+clm_a2l%forc_v**2),0.1_rk8)
    clm_a2l%forc_hgt_u = clm_a2l%forc_hgt
    clm_a2l%forc_hgt_t = clm_a2l%forc_hgt
    clm_a2l%forc_hgt_q = clm_a2l%forc_hgt
    !!!$acc end kernels

    do concurrent ( i = begg:endg )
      clm_a2l%forc_rho(i) = clm_a2l%forc_pbot(i)/(rgas*clm_a2l%forc_t(i))
      satq = pfwsat(real(clm_a2l%forc_t(i),rkx), &
                    real(clm_a2l%forc_pbot(i),rkx))
      clm_a2l%forc_rh(i) = clm_a2l%forc_q(i)/satq
      clm_a2l%forc_rh(i) = min(real(rhmax,rk8), clm_a2l%forc_rh(i))
      clm_a2l%forc_rh(i) = max(real(rhmin,rk8), clm_a2l%forc_rh(i))
      clm_a2l%forc_rh(i) = clm_a2l%forc_rh(i) * 100.0_rk8
      if ( clm_a2l%forc_t(i)-tfrz < tcrit ) then
        clm_a2l%forc_snow(i) = clm_a2l%rainf(i)
        clm_a2l%forc_rain(i) = d_zero
      else
        clm_a2l%forc_snow(i) = d_zero
        clm_a2l%forc_rain(i) = clm_a2l%rainf(i)
      end if
    end do
    ! Specific humidity
    !$acc kernels
    clm_a2l%forc_q = clm_a2l%forc_q/(1.0_rk8+clm_a2l%forc_q)
    !$acc end kernels

    ! interface chemistry / surface

    if ( ichem /= 1 ) then
      do i = begg, endg
        lat = real(adomain%xlat(i),rkx)
        clm_a2l%forc_pco2(i) = &
          ghgval(igh_co2,rcmtimer%year,rcmtimer%month,lat) * &
          clm_a2l%forc_psrf(i)
#ifdef LCH4
        clm_a2l%forc_pch4(i) = &
          ghgval(igh_ch4,rcmtimer%year,rcmtimer%month,lat) * &
          clm_a2l%forc_psrf(i)
#endif
      end do
      !clm_a2l%forc_ndep = 6.34e-5_rk8
      if ( use_c13 ) then
        !$acc kernels
        clm_a2l%forc_pc13o2 = c13ratio*clm_a2l%forc_pco2
        !$acc end kernels
      end if
      !$acc kernels
      clm_a2l%forc_po2 = o2_molar_const*clm_a2l%forc_psrf
      ! deposition of aerosol == zero
      clm_a2l%forc_aer(:,:) = d_zero !
      !$acc end kernels
    else
      !
      ! interface with atmospheric chemistry
      ! CO2 partial pressure (Pa)
      do i = begg, endg
        lat = real(adomain%xlat(i),rkx)
        clm_a2l%forc_pco2(i) = &
          ghgval(igh_co2,rcmtimer%year,rcmtimer%month,lat) * &
          clm_a2l%forc_psrf(i)
#ifdef LCH4
        clm_a2l%forc_pch4(i) = &
          ghgval(igh_ch4,rcmtimer%year,rcmtimer%month,lat) * &
          clm_a2l%forc_psrf(i)
#endif
      end do
      if ( use_c13 ) then
        ! C13O2 partial pressure (Pa)
        !$acc kernels
        clm_a2l%forc_pc13o2 = c13ratio*clm_a2l%forc_pco2
        !$acc end kernels
      end if
      ! O2 partial pressure (Pa)
      !$acc kernels
      clm_a2l%forc_po2 = o2_molar_const*clm_a2l%forc_psrf
      !$acc end kernels
      ! FAB: this is the interactive part
      ! a) snow ageing
      !    flux of species have to be consistent with snowhydrology use
      !    of indices unit is kg/m2/s
      ! b) drydeposition BC HL
      !    flux arriving through lm interface are accumulated between
      !    two surface call : needs to average with syncro_srf%rw
      if ( isnowdark == 1 ) then
        ! dry deposition BC HL
        if ( nbchl > 0 ) then
          !$acc kernels
          temps = d_zero
          !$acc end kernels
          do concurrent ( j = jci1:jci2, i = ici1:ici2 )
            !$acc loop seq
            do n = 1, nbchl
              temps(j,i) = temps(j,i) + &
                lm%drydepflx(j,i,ibchl(n)) * syncro_srf%rw
            end do
          end do
          call glb_c2l_gs(lndcomm,temps,clm_a2l%notused)
          !$acc kernels
          clm_a2l%forc_aer(:,1) = clm_a2l%notused
          !$acc end kernels
        end if
        ! drydeposition BC HB
        if ( ibchb > 0 ) then
          !$acc kernels
          temps(:,:) = lm%drydepflx(jci1:jci2,ici1:ici2,ibchb) * syncro_srf%rw
          !$acc end kernels
          call glb_c2l_gs(lndcomm,temps,clm_a2l%notused)
          !$acc kernels
          clm_a2l%forc_aer(:,2) = clm_a2l%notused
          !$acc end kernels
        end if
        ! wet dep BC (sum rainout and washout fluxes, sum hb amd hl)
        if ( ibchb > 0 .and. nbchl > 0 ) then
          do concurrent ( j = jci1:jci2, i = ici1:ici2 )
            temps(j,i) = lm%wetdepflx(j,i,ibchb) * syncro_srf%rw
            !$acc loop seq
            do n = 1, nbchl
              temps(j,i) = temps(j,i) + &
                   lm%wetdepflx(j,i,ibchl(n)) * syncro_srf%rw
            end do
          end do
          call glb_c2l_gs(lndcomm,temps,clm_a2l%notused)
          !$acc kernels
          clm_a2l%forc_aer(:,3) = clm_a2l%notused
          !$acc end kernels
        end if
        ! drydeposition OC HL
        if ( nochl > 0 ) then
          !$acc kernels
          temps = d_zero
          !$acc end kernels
          do concurrent ( j = jci1:jci2, i = ici1:ici2 )
            !$acc loop seq
            do n = 1, nochl
              temps(j,i) = temps(j,i) + &
                lm%drydepflx(j,i,iochl(n)) * syncro_srf%rw
            end do
          end do
          call glb_c2l_gs(lndcomm,temps,clm_a2l%notused)
          !$acc kernels
          clm_a2l%forc_aer(:,4) = clm_a2l%notused
          !$acc end kernels
        end if
        ! drydeposition OC HB
        if ( iochb > 0 ) then
          !$acc kernels
          temps(:,:) =  lm%drydepflx(jci1:jci2,ici1:ici2,iochb) * syncro_srf%rw
          !$acc end kernels
          call glb_c2l_gs(lndcomm,temps,clm_a2l%notused)
          !$acc kernels
          clm_a2l%forc_aer(:,5) = clm_a2l%notused
          !$acc end kernels
        end if
        ! wet dep OC (sum rainout and washout fluxes, sum hb and hl)
        if ( iochb > 0 .and. nochl > 0 ) then
          do concurrent ( j = jci1:jci2, i = ici1:ici2 )
            temps(j,i) = lm%wetdepflx(j,i,iochb) * syncro_srf%rw
            !$acc loop seq
            do n = 1, nochl
              temps(j,i) = temps(j,i) + &
                     lm%wetdepflx(j,i,iochl(n)) * syncro_srf%rw
            end do
          end do
          call glb_c2l_gs(lndcomm,temps,clm_a2l%notused)
          !$acc kernels
          clm_a2l%forc_aer(:,6) = clm_a2l%notused
          !$acc end kernels
        end if
        if ( size(lm%idust) == 4 ) then
          ! wet dep dust 1
          !$acc kernels
          temps(:,:) = (lm%wetdepflx(jci1:jci2,ici1:ici2,lm%idust(1)) + &
                lm%wetdepflx(jci1:jci2,ici1:ici2,lm%idust(1))) * syncro_srf%rw
          !$acc end kernels
          call glb_c2l_gs(lndcomm,temps,clm_a2l%notused)
          !$acc kernels
          clm_a2l%forc_aer(:,7) = clm_a2l%notused
          ! dry dep dust 1
          temps(:,:) = lm%drydepflx(jci1:jci2,ici1:ici2,lm%idust(1)) * &
                       syncro_srf%rw
          !$acc end kernels
          call glb_c2l_gs(lndcomm,temps,clm_a2l%notused)
          !$acc kernels
          clm_a2l%forc_aer(:,8) = clm_a2l%notused

          ! wet dep dust 2
          temps(:,:) =(lm%wetdepflx(jci1:jci2,ici1:ici2,lm%idust(2)) + &
                 lm%wetdepflx(jci1:jci2,ici1:ici2,lm%idust(2))) * syncro_srf%rw
          !$acc end kernels
          call glb_c2l_gs(lndcomm,temps,clm_a2l%notused)
          !$acc kernels
          clm_a2l%forc_aer(:,9) = clm_a2l%notused
          ! dry dep dust 2
          temps(:,:) = lm%drydepflx(jci1:jci2,ici1:ici2,lm%idust(2)) * &
                       syncro_srf%rw
          !$acc end kernels
          call glb_c2l_gs(lndcomm,temps,clm_a2l%notused)
          !$acc kernels
          clm_a2l%forc_aer(:,10) = clm_a2l%notused

          ! wet dep dust 3
          temps(:,:) = (lm%wetdepflx(jci1:jci2,ici1:ici2,lm%idust(3)) + &
                 lm%wetdepflx(jci1:jci2,ici1:ici2,lm%idust(3))) * syncro_srf%rw
          !$acc end kernels
          call glb_c2l_gs(lndcomm,temps,clm_a2l%notused)
          !$acc kernels
          clm_a2l%forc_aer(:,11) = clm_a2l%notused
          ! dry dep dust 3
          temps(:,:) = lm%drydepflx(jci1:jci2,ici1:ici2,lm%idust(3)) * &
                       syncro_srf%rw
          !$acc end kernels
          call glb_c2l_gs(lndcomm,temps,clm_a2l%notused)
          !$acc kernels
          clm_a2l%forc_aer(:,12) = clm_a2l%notused

          ! wet dep dust 4
          temps(:,:) = (lm%wetdepflx(jci1:jci2,ici1:ici2,lm%idust(4)) + &
                  lm%wetdepflx(jci1:jci2,ici1:ici2,lm%idust(4))) * syncro_srf%rw
          !$acc end kernels
          call glb_c2l_gs(lndcomm,temps,clm_a2l%notused)
          !$acc kernels
          clm_a2l%forc_aer(:,13) = clm_a2l%notused
          ! dry dep dust 4
          temps(:,:) = lm%drydepflx(jci1:jci2,ici1:ici2,lm%idust(4)) * &
                       syncro_srf%rw
          !$acc end kernels
          call glb_c2l_gs(lndcomm,temps,clm_a2l%notused)
          !$acc kernels
          clm_a2l%forc_aer(:,14) = clm_a2l%notused
          !$acc end kernels
        end if
      end if
      if ( ichsursrc == 1 .and. ino > 0 .and. ichbion == 1 ) then
        ndep_nochem = .false.
        !$acc kernels
        temps(:,:) = (lm%drydepflx(jci1:jci2,ici1:ici2,ino) + &
                      lm%wetdepflx(jci1:jci2,ici1:ici2,ino)) * syncro_srf%rw
        !$acc end kernels
        call glb_c2l_gs(lndcomm,temps,clm_a2l%notused)
        ! to convert NO flux from mg/m2/day to g/m2/sec,
        ! we multiply by 1.15740741 * 10-8
        !$acc kernels
        clm_a2l%forc_ndep(:) = clm_a2l%notused * 1.15740741e-8_rk8
        !$acc end kernels
      end if
      ! FAB to do : treat the 12 bins case ..
      ! b : pass the nitrogen deposition flux
    end if ! end test on ichem

    if ( .true. ) then
      !$acc kernels
      clm_a2l%forc_flood = 0.0_rk8
      clm_a2l%volr = 0.0_rk8
      !$acc end kernels
    else
      ! Runoff in input ? Chym ?
      ! clm_a2l%forc_flood  ! flood (mm/s)
      ! clm_a2l%volr        ! rof volr (m3)
    end if

    contains

#include <pfwsat.inc>

  end subroutine atmosphere_to_land

  subroutine land_to_atmosphere(lm,lms)
    implicit none
    type(lm_exchange), intent(inout) :: lm
    type(lm_state), intent(inout) :: lms
    integer(ik4) :: i, j, n, k, begg, endg

    call get_proc_bounds(begg,endg)

    ! Get back data from clm_l2a
    call glb_l2c_ss(lndcomm,clm_l2a%t_rad,lms%tgbb)

    call glb_l2c_ss(lndcomm,clm_l2a%t_ref2m,lms%t2m)
    call glb_l2c_ss(lndcomm,clm_l2a%q_ref2m,lms%q2m)

    ! CLM gives just wind speed, assume directions are same as input.
    !$acc kernels
    clm_a2l%notused = clm_l2a%u_ref10m/clm_a2l%forc_wind
    clm_l2a%notused = clm_a2l%forc_u*clm_a2l%notused
    !$acc end kernels
    call glb_l2c_ss(lndcomm,clm_l2a%notused,lms%u10m)
    !$acc kernels
    clm_l2a%notused = clm_a2l%forc_v*clm_a2l%notused
    !$acc end kernels
    call glb_l2c_ss(lndcomm,clm_l2a%notused,lms%v10m)

    call glb_l2c_ss(lndcomm,clm_l2a%eflx_sh_tot,lms%sent)
    call glb_l2c_ss(lndcomm,clm_l2a%qflx_evap_tot,lms%evpr)
    call glb_l2c_ss(lndcomm,clm_l2a%ram1,lms%ram1)
    call glb_l2c_ss(lndcomm,clm_l2a%rah1,lms%rah1)
    call glb_l2c_ss(lndcomm,clm_l2a%br1,lms%br)
    !$acc kernels
    clm_l2a%notused = sqrt(clm_l2a%taux**2+clm_l2a%tauy**2)/clm_a2l%forc_wind
    !$acc end kernels
    call glb_l2c_ss(lndcomm,clm_l2a%notused,lms%drag)

    call glb_l2c_ss(lndcomm,clm_l2a%h2osno,lms%sncv)
    ! For CLM45, we collect total fraction of ground emitting dust.
    call glb_l2c_ss(lndcomm,clm_l2a%vdustfrac,lms%wt)
    call glb_l2c_ss(lndcomm,clm_l2a%taux,lms%taux)
    call glb_l2c_ss(lndcomm,clm_l2a%tauy,lms%tauy)
    !$acc kernels
    clm_l2a%zom = max(clm_l2a%zom,1.0e-4_rk8)
    !$acc end kernels
    call glb_l2c_ss(lndcomm,clm_l2a%zom,lms%zo)
    call glb_l2c_ss(lndcomm,clm_l2a%t_veg,lms%tlef)

    ! soil temperature profile
    ! note tsw is used as a temporary table
    !$acc kernels
    clm_l2a%notused = 0.0_rk8
    !$acc end kernels
    do k = 1, nlevsoi
      !$acc kernels
      clm_l2a%notused(:) = clm_l2a%tsoi(:,k)
      !$acc end kernels
      call glb_l2c_ss(lndcomm,clm_l2a%notused,lms%tsw)
      !$acc kernels
      lms%tsoi(:,:,:,k) = lms%tsw(:,:,:)
      !$acc end kernels
    end do

    do concurrent ( n = 1:nnsg, j = jci1:jci2, i = ici1:ici2 )
      if ( lm%ldmsk1(n,j,i) == 1 ) then
        lms%tgrd(n,j,i) = lms%tsoi(n,j,i,1)
        lms%tgbrd(n,j,i) = lms%tsoi(n,j,i,1)
      end if
    end do
    call glb_l2c_ss(lndcomm,clm_l2a%eflx_gnet,lms%hfso)

    !$acc kernels
    clm_l2a%notused = 0.0_rk8
    !$acc end kernels

    do k = 1, nlevsoi
      !$acc kernels
      clm_l2a%notused(:) = clm_l2a%h2osoi(:,k)
      !$acc end kernels
      call glb_l2c_ss(lndcomm,clm_l2a%notused,lms%tsw)
      !$acc kernels
      lms%sw(:,:,:,k) = lms%tsw(:,:,:)
      !$acc end kernels
    end do

    ! volumetric soil water profile (m3/m3) is saved for chem
    !$acc kernels
    clm_l2a%notused = 0.0_rk8
    !$acc end kernels
    do k = 1, nlevsoi
      !$acc kernels
      clm_l2a%notused(:) = clm_l2a%h2osoi_vol(:,k)
      !$acc end kernels
      call glb_l2c_ss(lndcomm,clm_l2a%notused,lms%tsw)
      !$acc kernels
      lms%sw_vol(:,:,:,k) = lms%tsw(:,:,:)
      !$acc end kernels
    end do
    ! tsw is finally passed as the soil water in Kg/m2 in the first
    ! 10 cm of soil
    ! a bit obsolete since volumetric water profile is passed ,
    ! but still could be potentuially usefull
    call glb_l2c_ss(lndcomm,clm_l2a%h2o10cm,lms%tsw)

    call glb_l2c_ss(lndcomm,clm_l2a%qflx_surf,lms%srnof)
    call glb_l2c_ss(lndcomm,clm_l2a%qflx_tot,lms%trnof)
    call glb_l2c_ss(lndcomm,clm_l2a%qflx_snow_melt,lms%snwm)
    !$acc kernels
    lms%snwm = lms%snwm * dtsrf
    !$acc end kernels

    ! From the input
    call glb_l2c_ss(lndcomm,clm_a2l%rainf,lms%prcp)
    call glb_l2c_ss(lndcomm,clm_a2l%forc_psrf,lms%sfcp)

    ! If CLM45, the surface emissivity is not used.
    ! The CLM outputs directly to RegCM the radiant Temperature.
    ! We fill it here the output not to leave it empy, but it is not
    ! used in computing the surface Long Wave Radiation
    !clm_l2a%notused = 1.0_rkx
    !call glb_l2c_ss(lndcomm,clm_l2a%notused,lms%emisv)
    do concurrent ( n = 1:nnsg, j = jci1:jci2, i = ici1:ici2 )
      if ( lm%ldmsk1(n,j,i) == 1 ) lms%emisv(n,j,i) = 1.0_rkx
    end do

    !--------------------------------------------------
    ! From land to chemistry
    ! only for Isoprene, and in kg/m^2/sec
    ! FAB add compatibility for other biogenic species /
    ! TO BE UPDATED IF THE CHEM MECHANISM CHANGES (e.g add limonene, pinene etc)
    ! passed to the chemistry scheme for the right  mechanism tracer index
    ! use temporary table emis2d for calling glb_l2c_ss

    if ( ichem == 1 ) then

      call glb_l2c_ss(lndcomm,clm_l2a%tlai,lms%xlai)

      if ( enable_megan_emission ) then
        emis2d = 0.0_rk8
        do k = 1, shr_megan_mechcomps_n
          if (shr_megan_mechcomps(k)%name == 'ISOP' .and. iisop > 0 ) then
            clm_l2a%notused(:) = clm_l2a%flxvoc(:,k)
            call glb_l2c_ss(lndcomm, clm_l2a%notused, emis2d)
            lms%vocemiss(:,:,:,iisop) = real(emis2d,rkx)
          end if
          ! add compatibility for other biogenic species !! /
          !
        end do
      end if

      ! pass the CLM dust flux to regcm
      ! nb: CLM considers 4 emission bins roughly equivalent to the regcm
      ! 4 bin version.
      ! if use the regcm 12 bin, the total mass is redistributed
      ! (chemlib/mod_che_dust)
      if ( ichdustemd == 3 ) then
        do k = 1, 4
          emis2d = 0.0_rk8
          clm_l2a%notused(:) = clm_l2a%flxdst(:,k)
          call glb_l2c_ss(lndcomm, clm_l2a%notused, emis2d)
          lms%dustemiss(:,:,:,k) = real(emis2d,rkx)
        end do
      end if
#ifdef LCH4
      ! factor 1.33 is to convert from kgC to kgCH4
      if ( ich4 > 0 ) then
        emis2d = 0.0_rk8
        clm_l2a%notused(:) = clm_l2a%flux_ch4(:) * 1.33_rk8
        call glb_l2c_ss(lndcomm, clm_l2a%notused, emis2d)
        lms%vocemiss(:,:,:,ich4) = real(emis2d,rkx)
      end if
#endif
     !FAB : test  drydep velocities for gas, based on w98 / improve by passing the real indices
      if (n_drydep > 0) then
        do k = 1, n_drydep
          emis2d = 0.0_rk8
          if (4 > 0 ) then
            clm_l2a%notused(:) = clm_l2a%ddvel(:,k)
            call glb_l2c_ss(lndcomm, clm_l2a%notused, emis2d)
            lms%ddepv(:,:,:,4) = real(emis2d,rkx)
          end if
        end do
      end if
    end if
    !--------------------------------------------------
    ! Will fix
    !clm_l2a%eflx_lwrad_out
    !clm_l2a%emv
    !clm_l2a%emg
    !clm_l2a%fsa
    !clm_l2a%nee
    !clm_l2a%rofliq
    !clm_l2a%rofice
    !clm_l2a%ddvel
  end subroutine land_to_atmosphere

  subroutine read_cru_pre(lm,temps)
    implicit none
    type(lm_exchange), intent(inout) :: lm
    real(rkx), dimension(:,:), pointer, contiguous, intent(inout) :: temps
    integer(ik4), save :: ncid = -1
    integer(ik4), save :: ivar = -1
    integer(ik4), save :: imon = -1
    integer(ik4), save :: iday = -1
    integer(ik4) :: istatus, idimid
    integer(ik4) :: nlat, nlon
    character(len=256) :: fname
    integer(ik4) :: iy, im, id, ih, imm, iss, ndm
    integer(ik4) :: iym1, iyp1, imm1, imp1
    integer(ik4) :: i, j
    real(rk8) :: pm, f1, f2, m1, m2
    integer(ik4), dimension(3), save :: istart, icount

    call split_idate(nextdate,iy,im,id,ih,imm,iss)

    if ( myid == iocpu ) then
      if ( ncid == -1 ) then
        call getmem(alon,jcross1,jcross2,icross1,icross2,'clm45:alon')
        call getmem(alat,jcross1,jcross2,icross1,icross2,'clm45:alat')
        call getmem(rcp,jcross1,jcross2,icross1,icross2,'clm45:rcp')
        call grid_collect(lm%xlon,alon,jci1,jci2,ici1,ici2)
        call grid_collect(lm%xlat,alat,jci1,jci2,ici1,ici2)
        fname = trim(inpglob)//pthsep//'CLM45'//pthsep// &
                'crudata'//pthsep//trim(crufile)
        istatus = nf90_open(fname,nf90_nowrite,ncid)
        if ( istatus /= nf90_noerr ) then
          write (stderr, *) nf90_strerror(istatus)
          call fatal(__FILE__,__LINE__,'ERROR OPEN CRU DATASET')
        end if
        istatus = nf90_inq_dimid(ncid,'lat',idimid)
        if ( istatus /= nf90_noerr ) then
          write (stderr, *) nf90_strerror(istatus)
          call fatal(__FILE__,__LINE__,'ERROR SEARCH DIM LAT IN CRU DATASET')
        end if
        istatus = nf90_inquire_dimension(ncid,idimid,len=nlat)
        if ( istatus /= nf90_noerr ) then
          write (stderr, *) nf90_strerror(istatus)
          call fatal(__FILE__,__LINE__,'ERROR READ DIM LAT IN CRU DATASET')
        end if
        istatus = nf90_inq_dimid(ncid,'lon',idimid)
        if ( istatus /= nf90_noerr ) then
          write (stderr, *) nf90_strerror(istatus)
          call fatal(__FILE__,__LINE__,'ERROR SEARCH DIM LON IN CRU DATASET')
        end if
        istatus = nf90_inquire_dimension(ncid,idimid,len=nlon)
        if ( istatus /= nf90_noerr ) then
          write (stderr, *) nf90_strerror(istatus)
          call fatal(__FILE__,__LINE__,'ERROR READ DIM LON IN CRU DATASET')
        end if
        call getmem(crupre,1,nlon,1,nlat,'clm45:crupre')
        call getmem(glon,1,nlon,'clm45:glon')
        call getmem(glat,1,nlat,'clm45:glat')
        call getmem(acp0,jci1,jci2,ici1,ici2,'clm45:acp0')
        call getmem(acp1,jci1,jci2,ici1,ici2,'clm45:acp1')
        call getmem(acp2,jci1,jci2,ici1,ici2,'clm45:acp2')
        if ( lcru_rand ) then
          call getmem(ihfac,1,24,'clm45:ihfac')
        end if
        istatus = nf90_inq_varid(ncid,'lat',ivar)
        if ( istatus /= nf90_noerr ) then
          write (stderr, *) nf90_strerror(istatus)
          call fatal(__FILE__,__LINE__,'ERROR SEARCH VAR LAT IN CRU DATASET')
        end if
        istatus = nf90_get_var(ncid,ivar,glat)
        if ( istatus /= nf90_noerr ) then
          write (stderr, *) nf90_strerror(istatus)
          call fatal(__FILE__,__LINE__,'ERROR READ VAR LAT FROM CRU DATASET')
        end if
        istatus = nf90_inq_varid(ncid,'lon',ivar)
        if ( istatus /= nf90_noerr ) then
          write (stderr, *) nf90_strerror(istatus)
          call fatal(__FILE__,__LINE__,'ERROR SEARCH VAR LON IN CRU DATASET')
        end if
        istatus = nf90_get_var(ncid,ivar,glon)
        if ( istatus /= nf90_noerr ) then
          write (stderr, *) nf90_strerror(istatus)
          call fatal(__FILE__,__LINE__,'ERROR READ VAR LON FROM CRU DATASET')
        end if
        istatus = nf90_inq_varid(ncid,'pre',ivar)
        if ( istatus /= nf90_noerr ) then
          write (stderr, *) nf90_strerror(istatus)
          call fatal(__FILE__,__LINE__,'ERROR SEARCH VAR PRE IN CRU DATASET')
        end if
        istart(1) = 1
        istart(2) = 1
        icount(1) = nlon
        icount(2) = nlat
        icount(3) = 1
        call h_interpolator_create(hint,glat,glon,alat,alon)
        write (stdout,*) 'Reading CRU Precipitation for month ', im
        iym1 = iy
        imm1 = im - 1
        if ( imm1 == 0 ) then
          imm1 = 12
          iym1 = iy - 1
        end if
        istart(3) = imm1
        istatus = nf90_get_var(ncid,ivar,crupre,istart,icount)
        if ( istatus /= nf90_noerr ) then
          write (stderr, *) nf90_strerror(istatus)
          call fatal(__FILE__,__LINE__,'ERROR READ VAR PRE FROM CRU DATASET')
        end if
        call h_interpolate_nn(hint,crupre,rcp)
        ndm = ndaypm(iym1,imm1,nextdate%calendar)
        rcp = rcp / real(ndm,rk8)
        call grid_distribute(rcp,acp0,jci1,jci2,ici1,ici2)
        istart(3) = im
        istatus = nf90_get_var(ncid,ivar,crupre,istart,icount)
        if ( istatus /= nf90_noerr ) then
          write (stderr, *) nf90_strerror(istatus)
          call fatal(__FILE__,__LINE__,'ERROR READ VAR PRE FROM CRU DATASET')
        end if
        call h_interpolate_nn(hint,crupre,rcp)
        ndm = ndaypm(iy,im,nextdate%calendar)
        rcp = rcp / real(ndm,rk8)
        call grid_distribute(rcp,acp1,jci1,jci2,ici1,ici2)
        iyp1 = iy
        imp1 = im + 1
        if ( imp1 == 13 ) then
          iyp1 = iy + 1
          imp1 = 1
        end if
        istart(3) = imp1
        istatus = nf90_get_var(ncid,ivar,crupre,istart,icount)
        if ( istatus /= nf90_noerr ) then
          write (stderr, *) nf90_strerror(istatus)
          call fatal(__FILE__,__LINE__,'ERROR READ VAR PRE FROM CRU DATASET')
        end if
        call h_interpolate_nn(hint,crupre,rcp)
        ndm = ndaypm(iyp1,imp1,nextdate%calendar)
        rcp = rcp / real(ndm,rk8)
        call grid_distribute(rcp,acp2,jci1,jci2,ici1,ici2)
        imon = im
      end if
      if ( imon /= im ) then
        acp0 = acp1
        acp1 = acp2
        write (stdout,*) 'Reading CRU Precipitation for month ', im
        iyp1 = iy
        imp1 = im + 1
        if ( imp1 == 13 ) then
          imp1 = 1
          iyp1 = iy + 1
        end if
        istart(3) = imp1
        istatus = nf90_get_var(ncid,ivar,crupre,istart,icount)
        if ( istatus /= nf90_noerr ) then
          write (stderr, *) nf90_strerror(istatus)
          call fatal(__FILE__,__LINE__,'ERROR READ VAR PRE FROM CRU DATASET')
        end if
        call h_interpolate_nn(hint,crupre,rcp)
        ndm = ndaypm(iyp1,imp1,nextdate%calendar)
        rcp = rcp / real(ndm,rk8)
        call grid_distribute(rcp,acp2,jci1,jci2,ici1,ici2)
        imon = im
      end if
      if ( lcru_rand ) then
        if ( iday /= id ) then
          call random_pick(24.0_rkx,ihfac,24)
          call bcast(ihfac)
          iday = id
        end if
      end if
    else
      if ( ncid == -1 ) then
        call grid_collect(lm%xlon,alon,jci1,jci2,ici1,ici2)
        call grid_collect(lm%xlat,alat,jci1,jci2,ici1,ici2)
        call getmem(acp0,jci1,jci2,ici1,ici2,'clm45:acp0')
        call getmem(acp1,jci1,jci2,ici1,ici2,'clm45:acp1')
        call getmem(acp2,jci1,jci2,ici1,ici2,'clm45:acp2')
        if ( lcru_rand ) then
          call getmem(ihfac,1,24,'clm45:ihfac')
        end if
        call grid_distribute(rcp,acp0,jci1,jci2,ici1,ici2)
        call grid_distribute(rcp,acp1,jci1,jci2,ici1,ici2)
        call grid_distribute(rcp,acp2,jci1,jci2,ici1,ici2)
        ncid = 0
        imon = im
      end if
      if ( imon /= im ) then
        acp0 = acp1
        acp1 = acp2
        call grid_distribute(rcp,acp2,jci1,jci2,ici1,ici2)
        imon = im
      end if
      if ( lcru_rand ) then
        if ( iday /= id ) then
          call bcast(ihfac)
          iday = id
        end if
      end if
    end if
    ndm = ndaypm(iy,im,nextdate%calendar)
    pm = real(id,rk8)/real(ndm,rk8)
    f1 = max(0.5_rk8 - pm,0.0_rkx)
    f2 = max(pm - 0.5_rk8,0.0_rkx)
    if ( lcru_rand ) then
      do i = ici1, ici2
        do j = jci1, jci2
          if ( acp0(j,i) > 0.0_rk8 .and. acp0(j,i) < 10000.0_rk8 ) then
            m1 = acp0(j,i) * f1 + acp1(j,i) * (1.0_rk8-f1)
            m2 = acp1(j,i) * (1.0_rk8-f2) + acp2(j,i) * f2
            temps(j,i) = (m1+m2) / (2.0_rk8*secpd) * ihfac(ih+1)
          end if
        end do
      end do
    else
      do i = ici1, ici2
        do j = jci1, jci2
          if ( acp0(j,i) > 0.0_rk8 .and. acp0(j,i) < 10000.0_rk8 ) then
            m1 = acp0(j,i) * f1 + acp1(j,i) * (1.0_rk8-f1)
            m2 = acp1(j,i) * (1.0_rk8-f2) + acp2(j,i) * f2
            temps(j,i) = (m1+m2) / (2.0_rk8*secpd)
          end if
        end do
      end do
    end if
  end subroutine read_cru_pre

end module mod_clm_regcm
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
