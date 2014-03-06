module mod_clm_regcm
  use mod_date
  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_runparams
  use mod_mppparam
  use mod_mpmessage
  use mod_constants
  use mod_regcm_types
  use mod_clm_initialize
  use mod_clm_driver
  use mod_clm_varctl , only : use_c13
  use mod_clm_atmlnd , only : clm_a2l , clm_l2a
  use mod_clm_decomp , only : procinfo

  private

  public :: initclm45 , runclm45 , albedoclm45

  real(rk8) , dimension(:,:) , pointer :: rprec , rsnow

  contains

  subroutine initclm45(cl,lm,lms)
    implicit none
    type (masked_comm) , intent(in) :: cl
    type(lm_exchange) , intent(inout) :: lm
    type(lm_state) , intent(inout) :: lms
    integer(ik4) :: i , j , n
    call initialize1(cl)
    call initialize2( )
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
                ! We should take this from CLM !
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
    call getmem2d(rprec,jci1,jci2,ici1,ici2,'initclm45:rprec')
    call getmem2d(rsnow,jci1,jci2,ici1,ici2,'initclm45:rsnow')
  end subroutine initclm45

  subroutine runclm45(lm,lms)
    implicit none
    type(lm_exchange) , intent(inout) :: lm
    type(lm_state) , intent(inout) :: lms
    integer(ik4) :: begg , endg , i
    real(rk8) :: hl , satvp
    logical :: doalb , rstwr , nlend
    character(len=64) :: rdate

    rprec = (lm%cprate+lm%ncprate) * rtsrf
    rsnow = 0.0D0

    where ( lm%tatm < 0.0D0 )
      rsnow = rprec
      rprec = 0.0D0
    end where

    begg = lbound(clm_a2l%forc_rh,1)
    endg = ubound(clm_a2l%forc_rh,1)

    ! Fill clm_a2l
    call c2l_gs(procinfo%cl,lm%tatm,clm_a2l%forc_t)
    call c2l_gs(procinfo%cl,lm%uatm,clm_a2l%forc_u)
    call c2l_gs(procinfo%cl,lm%vatm,clm_a2l%forc_v)
    call c2l_gs(procinfo%cl,lm%qvatm,clm_a2l%forc_q)
    call c2l_gs(procinfo%cl,lm%hgt,clm_a2l%forc_hgt)
    call c2l_gs(procinfo%cl,lm%patm,clm_a2l%forc_pbot)
    call c2l_gs(procinfo%cl,lm%thatm,clm_a2l%forc_th)
    ! forc_vp IS UNUSED
    call c2l_gs(procinfo%cl,lm%rhox,clm_a2l%forc_rho)
    call c2l_gs(procinfo%cl,lm%sfps,clm_a2l%forc_psrf)
    call c2l_gs(procinfo%cl,lm%dwrlwf,clm_a2l%forc_lwrad)
    call c2l_gs(procinfo%cl,lm%solar,clm_a2l%forc_solar)
    call c2l_gs(procinfo%cl,rprec,clm_a2l%forc_rain)
    call c2l_gs(procinfo%cl,rsnow,clm_a2l%forc_snow)

    call c2l_gs(procinfo%cl,lm%swdir,clm_a2l%notused)
    clm_a2l%forc_solad(:,1) = clm_a2l%notused
    call c2l_gs(procinfo%cl,lm%lwdir,clm_a2l%notused)
    clm_a2l%forc_solad(:,2) = clm_a2l%notused
    call c2l_gs(procinfo%cl,lm%swdif,clm_a2l%notused)
    clm_a2l%forc_solai(:,1) = clm_a2l%notused
    call c2l_gs(procinfo%cl,lm%lwdif,clm_a2l%notused)
    clm_a2l%forc_solai(:,2) = clm_a2l%notused

    ! Compute or alias
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
    end do

    if ( ichem /= 1 ) then
      clm_a2l%forc_pco2 = 0.039*clm_a2l%forc_psrf
      clm_a2l%forc_ndep = 6.34D-5
      if ( use_c13 ) then
        clm_a2l%forc_pc13o2 = 0.0001*clm_a2l%forc_psrf
      end if
      clm_a2l%forc_po2 = 0.2095*clm_a2l%forc_psrf
      clm_a2l%forc_aer = 0.0D0
    else
      ! Species partial pressures ! Fabien?
      ! clm_a2l%forc_pco2   ! CO2 partial pressure (Pa)
      ! clm_a2l%forc_ndep   ! nitrogen deposition rate (gN/m2/s)
      ! clm_a2l%forc_pc13o2 ! C13O2 partial pressure (Pa)
      ! clm_a2l%forc_po2    ! O2 partial pressure (Pa)
      ! clm_a2l%forc_aer    ! aerosol deposition array
    end if

    ! Runoff in input ? Chym ?
    ! clm_a2l%forc_flood  ! flood (mm/s)
    ! clm_a2l%volr        ! rof volr (m3)

    ! Run CLM
    if ( ktau == 0 .or. mod(ktau+1,ntrad) == 0 ) then
      doalb = .true.
    else
      doalb = .false.
    end if
    if ( idatex == idate2 .or. &
       (lfdomonth(idatex) .and. lmidnight(idatex)) ) then
      rstwr = .true.
      write(rdate,'(i10)') toint10(idatex)
    else
      rstwr = .false.
    end if
    ! First declin should be for NEXT time step. Stay quiet for now...
    call clm_drv(doalb, calday, declin, declin, rstwr, nlend, rdate)

    call fatal(__FILE__,__LINE__,'DONE CLM STEP!!!!')
    ! Get back data from clm_l2a

  end subroutine runclm45

  subroutine albedoclm45(lm,lms)
    implicit none
    type(lm_exchange) , intent(inout) :: lm
    type(lm_state) , intent(inout) :: lms
    ! Just get albedoes from clm_l2a
  end subroutine albedoclm45

end module mod_clm_regcm
