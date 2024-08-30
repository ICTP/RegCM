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

module mod_micro_wsm7

  use mod_realkinds
  use mod_intkinds
  use mod_constants
  use mod_dynparam
  use mod_memutil
  use mod_runparams , only : ichem , dt
  use mod_runparams , only : iqi , iqc , iqr , iqg , iqh , iqs , iqv
  use mod_runparams , only : rdt
  use mod_regcm_types

  private

  real(rkx) , parameter :: psat = 610.78_rkx

  ! maximum time step for minor loops
  real(rkx) , parameter :: dtcldcr = 120.0_rkx
  ! intercept parameter rain
  real(rkx) , parameter :: n0r = 8.e6_rkx
  ! intercept parameter graupel
  real(rkx) , parameter :: n0g = 4.e6_rkx
  ! intercept parameter hail
  real(rkx) , parameter :: n0h = 4.e4_rkx
  ! a constant for terminal velocity of rain
  real(rkx) , parameter :: avtr = 841.9_rkx
  ! a constant for terminal velocity of rain
  real(rkx) , parameter :: bvtr = 0.8_rkx
  ! a constant for terminal velocity of graupel
  real(rkx) , parameter :: avtg = 330.0_rkx
  ! a constant for terminal velocity of graupel
  real(rkx) , parameter :: bvtg = 0.8_rkx
  ! a constant for terminal velocity of hail
  real(rkx) , parameter :: avth = 285.0_rkx
  ! a constant for terminal velocity of hail
  real(rkx) , parameter :: bvth = 0.8_rkx
  ! density of graupel
  real(rkx) , parameter :: deng = 500.0_rkx
  ! density of hail
  real(rkx) , parameter :: denh = 912.0_rkx
  ! 8 microm  in contrast to 10 micro m
  real(rkx) , parameter :: r0 = 0.8e-5_rkx
  ! collection efficiency
  real(rkx) , parameter :: peaut = 0.55_rkx
  ! maritime cloud in contrast to 3.e8 in tc80
  real(rkx) , parameter :: xncr = 3.e8_rkx
  ! the dynamic viscosity kgm-1s-1
  real(rkx) , parameter :: xmyu = 1.718e-5_rkx
  ! a constant for terminal velocity of snow
  real(rkx) , parameter :: avts = 11.72_rkx
  ! a constant for terminal velocity of snow
  real(rkx) , parameter :: bvts = 0.41_rkx
  ! maximum n0s (t=-90c unlimited)
  real(rkx) , parameter :: n0smax = 1.e11_rkx
  ! limited maximum value for slope parameter of rain
  real(rkx) , parameter :: lamdarmax = 8.e4_rkx
  ! limited maximum value for slope parameter of snow
  real(rkx) , parameter :: lamdasmax = 1.e5_rkx
  ! limited maximum value for slope parameter of graupel
  real(rkx) , parameter :: lamdagmax = 6.e4_rkx
  ! limited maximum value for slope parameter of hail
  real(rkx) , parameter :: lamdahmax = 2.e4_rkx
  ! constant for the cloud-ice diamter
  real(rkx) , parameter :: dicon = 11.9_rkx
  ! limited maximum value for the cloud-ice diamter
  real(rkx) , parameter :: dimax = 500.e-6_rkx
  ! temperature dependent intercept parameter snow
  real(rkx) , parameter :: n0s = 2.e6_rkx
  ! .122 exponen factor for n0s
  real(rkx) , parameter :: alpha = 0.12_rkx
  ! constant in biggs freezing
  real(rkx) , parameter :: pfrz1 = 100.0_rkx
  ! constant in biggs freezing
  real(rkx) , parameter :: pfrz2 = 0.66_rkx
  ! minimun values for qr, qs, and qg
  real(rkx) , parameter :: qrsmin = 1.0e-9_rkx
  real(rkx) , parameter :: qcimin = 1.0e-8_rkx
  ! snow/cloud-water collection efficiency
  real(rkx) , parameter :: eacrc = 1.0_rkx

  real(rkx) , parameter :: minni = 1.0e3_rkx
  real(rkx) , parameter :: maxni = 1.0e6_rkx

  ! Hail/rain collection efficiency
  real(rkx) , parameter :: eachr = 1.0_rkx
  ! Hail/snow collection efficiency
  real(rkx) , parameter :: eachs = 1.0_rkx
  ! Hail/graupel collection efficiency
  real(rkx) , parameter :: eachg = 0.5_rkx
  ! Density of snow
  real(rkx) , parameter :: dens  =  100.0_rkx
  ! threshold amount for aggretion to occur
  real(rkx) , parameter :: qs0   =  6.e-4_rkx
  real(rkx) , parameter :: t00  = 238.16_rkx
  real(rkx) , parameter :: t01  = 273.16_rkx
  ! drag coefficient for hailsone
  real(rkx) , parameter :: cd = 0.6_rkx
  real(rkx) , parameter :: dldt  = cpv-cpw
  real(rkx) , parameter :: dldti = cpv-cpi
  real(rkx) , parameter :: xa = -dldt/rwat
  real(rkx) , parameter :: xb = xa + wlhv/(rwat*wattp)
  real(rkx) , parameter :: xai = -dldti/rwat
  real(rkx) , parameter :: xbi = xai + wlhs/(rwat*wattp)

  real(rkx) , save :: qc0 , qck1 , pidnc ,                            &
             bvtr1 , bvtr2 , bvtr3 , bvtr4 , g1pbr ,                  &
             g3pbr , g4pbr , g5pbro2 , pvtr , eacrr , pacrr ,         &
             bvtr6 , g6pbr ,                                          &
             precr1 , precr2 , roqimax , bvts1 ,                      &
             bvts2 , bvts3 , bvts4 , g1pbs , g3pbs , g4pbs ,          &
             g5pbso2 , pvts , pacrs , precs1 , precs2 , pidn0r ,      &
             pidn0s , pacrc , bvtg1 , bvtg2 , bvtg3 , bvtg4 , g1pbg , &
             g3pbg , g4pbg , g5pbgo2 , g6pbgh , pvtg , pacrg ,        &
             precg1 , precg2 , precg3 , pidn0g ,                      &
             bvth2 , bvth3 , bvth4 ,                                  &
             g3pbh , g4pbh , g5pbho2 , pvth , pacrh ,                 &
             prech1 , prech2 , prech3 , pidn0h ,                      &
             rslopermax , rslopesmax , rslopegmax , rslopehmax ,      &
             rsloperbmax , rslopesbmax , rslopegbmax , rslopehbmax ,  &
             rsloper2max , rslopes2max , rslopeg2max , rslopeh2max ,  &
             rsloper3max , rslopes3max , rslopeg3max , rslopeh3max

  public :: allocate_mod_wsm7 , init_wsm7 , wsm7

  integer(ik4) :: is , ie

  real(rkx) , dimension(:,:) , pointer :: t
  real(rkx) , dimension(:,:) , pointer :: qv
  real(rkx) , dimension(:,:,:) , pointer :: qs
  real(rkx) , dimension(:,:,:) , pointer :: rh
  real(rkx) , dimension(:,:,:) , pointer :: qci
  real(rkx) , dimension(:,:,:) , pointer :: qrs
  real(rkx) , dimension(:,:,:) , pointer :: fall
  real(rkx) , dimension(:,:) , pointer :: den
  real(rkx) , dimension(:,:) , pointer :: delz
  real(rkx) , dimension(:,:) , pointer :: p
!  real(rkx) , dimension(:,:) , pointer :: cloud_er
!  real(rkx) , dimension(:,:) , pointer :: ice_er
!  real(rkx) , dimension(:,:) , pointer :: snow_er
  real(rkx) , dimension(:) , pointer :: ptfac
  real(rkx) , dimension(:) , pointer :: rain
  real(rkx) , dimension(:) , pointer :: snow
  real(rkx) , dimension(:) , pointer :: grpl
  real(rkx) , dimension(:) , pointer :: hail

  contains

  subroutine allocate_mod_wsm7
    implicit none
    is = 1
    ie = ((ici2-ici1) + 1) * ((jci2-jci1) + 1)
    call getmem2d(t,is,ie,1,kz,'wsm7::t')
    call getmem2d(qv,is,ie,1,kz,'wsm7::qv')
    call getmem3d(qci,is,ie,1,kz,1,2,'wsm7::qci')
    call getmem3d(qrs,is,ie,1,kz,1,4,'wsm7::qrs')
    call getmem3d(qs,is,ie,1,kz,1,2,'wsm7::qs')
    call getmem3d(rh,is,ie,1,kz,1,2,'wsm7::rh')
    call getmem3d(fall,is,ie,1,kz,1,2,'wsm7::fall')
    call getmem2d(den,is,ie,1,kz,'wsm7::den')
    call getmem2d(delz,is,ie,1,kz,'wsm7::delz')
    call getmem2d(p,is,ie,1,kz,'wsm7::p')
    !call getmem2d(cloud_er,is,ie,1,kz,'wsm7::cloud_er')
    !call getmem2d(ice_er,is,ie,1,kz,'wsm7::ice_er')
    !call getmem2d(snow_er,is,ie,1,kz,'wsm7::snow_er')
    call getmem1d(ptfac,is,ie,'wsm7::ptfac')
    call getmem1d(rain,is,ie,'wsm7::rain')
    call getmem1d(snow,is,ie,'wsm7::snow')
    call getmem1d(grpl,is,ie,'wsm7::grpl')
    call getmem1d(hail,is,ie,'wsm7::hail')
  end subroutine allocate_mod_wsm7

  subroutine init_wsm7
    implicit none
    !denr=rhoh2o. den0=stdrho
    qc0  = fourt*mathpi*rhoh2o*r0**3*xncr/stdrho  ! 0.419e-3 -- .61e-3
    qck1 = 0.104_rkx*egrav*peaut / &
                  (xncr*rhoh2o)**(onet)/xmyu*stdrho**(fourt) ! 7.03
    pidnc = mathpi*rhoh2o/d_six        ! syb
    !cpv = 1885.0       ! specific heat of water vapor
    !n
    bvtr1 = d_one+bvtr
    bvtr2 = 2.5_rkx+d_half*bvtr
    bvtr3 = 3.0_rkx+bvtr
    bvtr4 = 4.0_rkx+bvtr
    bvtr6 = 6.0_rkx+bvtr
    g1pbr = rgmma(bvtr1)
    g3pbr = rgmma(bvtr3)
    g4pbr = rgmma(bvtr4)            ! 17.837825
    g6pbr = rgmma(bvtr6)
    g5pbro2 = rgmma(bvtr2)          ! 1.8273
    pvtr = avtr*g4pbr/d_six
    eacrr = 1.0_rkx
    pacrr = mathpi*n0r*avtr*g3pbr*0.25_rkx*eacrr
    precr1 = twopi*n0r*0.78_rkx
    precr2 = twopi*n0r*0.31_rkx*avtr**d_half*g5pbro2
    roqimax = 2.08e22_rkx*dimax**8
    bvts1 = d_one+bvts
    bvts2 = 2.5_rkx+d_half*bvts
    bvts3 = 3.0_rkx+bvts
    bvts4 = 4.0_rkx+bvts
    g1pbs = rgmma(bvts1)    !.8875
    g3pbs = rgmma(bvts3)
    g4pbs = rgmma(bvts4)    ! 12.0786
    g5pbso2 = rgmma(bvts2)
    pvts = avts*g4pbs/d_six
    pacrs = mathpi*n0s*avts*g3pbs*0.25_rkx
    precs1 = 4.0_rkx*n0s*0.65_rkx
    precs2 = 4.0_rkx*n0s*0.44_rkx*avts**d_half*g5pbso2
    pidn0r =  mathpi*rhoh2o*n0r
    pidn0s =  mathpi*rhosnow*n0s
    pacrc = mathpi*n0s*avts*g3pbs*0.25_rkx*eacrc

    bvtg1 = d_one+bvtg
    bvtg2 = 2.5_rkx+d_half*bvtg
    bvtg3 = 3.0_rkx+bvtg
    bvtg4 = 4.0_rkx+bvtg
    g1pbg = rgmma(bvtg1)
    g3pbg = rgmma(bvtg3)
    g4pbg = rgmma(bvtg4)
    pacrg = mathpi*n0g*avtg*g3pbg*0.25_rkx
    g5pbgo2 = rgmma(bvtg2)
    g6pbgh = rgmma(2.75_rkx)
    pvtg = avtg*g4pbg/d_six
    precg1 = twopi*n0g*0.78_rkx
    precg2 = twopi*n0g*0.31_rkx*avtg**d_half*g5pbgo2
    precg3 = twopi*n0g*0.31_rkx*g6pbgh*sqrt(sqrt(4.0_rkx*deng/3.0_rkx/cd))
    pidn0g =  mathpi*deng*n0g

    bvth2 = 2.5_rkx+d_half*bvth
    bvth3 = 3.0_rkx+bvth
    bvth4 = 4.0_rkx+bvth
    g3pbh = rgmma(bvth3)
    g4pbh = rgmma(bvth4)
    g5pbho2 = rgmma(bvth2)
    pacrh = mathpi*n0h*avth*g3pbh*0.25_rkx
    pvth = avth*g4pbh/d_six
    prech1 = twopi*n0h*0.78_rkx
    prech2 = twopi*n0h*0.31_rkx*avth**d_half*g5pbho2
    prech3 = twopi*n0h*0.31_rkx*g6pbgh*sqrt(sqrt(4.0_rkx*denh/3.0_rkx/cd))
    pidn0h = mathpi*denh*n0h

    rslopermax = d_one/lamdarmax
    rslopesmax = d_one/lamdasmax
    rslopegmax = d_one/lamdagmax
    rslopehmax = d_one/lamdahmax
    rsloperbmax = rslopermax ** bvtr
    rslopesbmax = rslopesmax ** bvts
    rslopegbmax = rslopegmax ** bvtg
    rslopehbmax = rslopehmax ** bvth
    rsloper2max = rslopermax * rslopermax
    rslopes2max = rslopesmax * rslopesmax
    rslopeg2max = rslopegmax * rslopegmax
    rslopeh2max = rslopehmax * rslopehmax
    rsloper3max = rsloper2max * rslopermax
    rslopes3max = rslopes2max * rslopesmax
    rslopeg3max = rslopeg2max * rslopegmax
    rslopeh3max = rslopeh2max * rslopehmax

    contains
    !
    ! rgmma function:  use infinite product form
    !
    pure real(rkx) function rgmma(x)
      implicit none
      real(rkx) , intent(in) :: x
      integer(ik4) , parameter :: imax = 10000
      real(rkx) :: y
      integer(ik4) :: i
      if ( abs(x-d_one) < epsilon(d_one) ) then
        rgmma = 0.0_rkx
      else
        rgmma = x * exp(m_euler*x)
        do i = 1 , imax
          y = real(i,rkx)
          rgmma = rgmma*(d_one+x/y)*exp(-x/y)
        end do
        rgmma = d_one/rgmma
      end if
    end function rgmma

  end subroutine init_wsm7

  subroutine wsm7(mo2mc,mc2mo)
    implicit none
    type(mod_2_micro) , intent(in) :: mo2mc
    type(micro_2_mod) , intent(out) :: mc2mo

    integer(ik4) :: i , j , k , kk , n
    real(rkx) :: pf1 , pf2 , pf3 , pf4 , qcw , totp

    ! to calculate effective radius for radiation
    !real(rkx) , dimension(kz) :: qv1d , t1d , p1d , qr1d , qs1d
    !real(rkx) , dimension(kz) :: den1d
    !real(rkx) , dimension(kz) :: qc1d
    !real(rkx) , dimension(kz) :: qi1d
    !real(rkx) , dimension(kz) :: re_qc , re_qi , re_qs


    if ( idynamic /= 3 ) then
      n = 1
      do i = ici1 , ici2
        do j = jci1 , jci2
          ptfac(n) = mo2mc%psb(j,i)*rdt
          n = n + 1
        end do
      end do
    else
      ptfac(:) = rdt
    end if

    do k = 1 , kz
      n = 1
      kk = kzp1-k
      do i = ici1 , ici2
        do j = jci1 , jci2
          t(n,kk) = mo2mc%t(j,i,k)
          p(n,kk) = mo2mc%phs(j,i,k)
          qv(n,kk) = mo2mc%qxx(j,i,k,iqv)
          qci(n,kk,1) = mo2mc%qxx(j,i,k,iqc)
          qci(n,kk,2) = mo2mc%qxx(j,i,k,iqi)
          qrs(n,kk,1) = mo2mc%qxx(j,i,k,iqr)
          qrs(n,kk,2) = mo2mc%qxx(j,i,k,iqs)
          qrs(n,kk,3) = mo2mc%qxx(j,i,k,iqg)
          qrs(n,kk,4) = mo2mc%qxx(j,i,k,iqh)
          delz(n,kk) = mo2mc%delz(j,i,k)
          den(n,kk) = mo2mc%rho(j,i,k)
          n = n + 1
        end do
      end do
    end do

    if ( ichem == 1 ) then
      mc2mo%remrat(:,:,:) = d_zero
    end if

    call wsm72d(dt,is,ie)

    do k = 1 , kz
      n = 1
      kk = kzp1 - k
      do i = ici1 , ici2
        do j = jci1 , jci2
          mc2mo%tten(j,i,k) = mc2mo%tten(j,i,k) + &
                  (t(n,kk)-mo2mc%t(j,i,k))*ptfac(n)
          mc2mo%qxten(j,i,k,iqv) = mc2mo%qxten(j,i,k,iqv) + &
                  (qv(n,kk)-mo2mc%qxx(j,i,k,iqv))*ptfac(n)
          mc2mo%qxten(j,i,k,iqc) = mc2mo%qxten(j,i,k,iqc) + &
                  (qci(n,kk,1)-mo2mc%qxx(j,i,k,iqc))*ptfac(n)
          mc2mo%qxten(j,i,k,iqi) = mc2mo%qxten(j,i,k,iqi) + &
                  (qci(n,kk,2)-mo2mc%qxx(j,i,k,iqi))*ptfac(n)
          mc2mo%qxten(j,i,k,iqr) = mc2mo%qxten(j,i,k,iqr) + &
                  (qrs(n,kk,1)-mo2mc%qxx(j,i,k,iqr))*ptfac(n)
          mc2mo%qxten(j,i,k,iqs) = mc2mo%qxten(j,i,k,iqs) + &
                  (qrs(n,kk,2)-mo2mc%qxx(j,i,k,iqs))*ptfac(n)
          mc2mo%qxten(j,i,k,iqg) = mc2mo%qxten(j,i,k,iqg) + &
                  (qrs(n,kk,3)-mo2mc%qxx(j,i,k,iqg))*ptfac(n)
          mc2mo%qxten(j,i,k,iqh) = mc2mo%qxten(j,i,k,iqh) + &
                  (qrs(n,kk,4)-mo2mc%qxx(j,i,k,iqh))*ptfac(n)
          n = n + 1
        end do
      end do
    end do

    if ( ichem == 1 ) then
      do k = 1 , kz
        n = 1
        kk = kzp1 - k
        do i = ici1 , ici2
          do j = jci1 , jci2
            if ( qrs(n,kk,1) > dlowval ) then
              pf1 = fall(n,kk,1)*delz(n,kk)/rhoh2o/qrs(n,kk,1)
            else
              pf1 = d_zero
            end if
            if ( qrs(n,kk,2) > dlowval ) then
              pf2 = fall(n,kk,2)*delz(n,kk)/rhoh2o/qrs(n,kk,2)
            else
              pf2 = d_zero
            end if
            if ( qrs(n,kk,3) > dlowval ) then
              pf3 = fall(n,kk,3)*delz(n,kk)/rhoh2o/qrs(n,kk,3)
            else
              pf3 = d_zero
            end if
            if ( qrs(n,kk,4) > dlowval ) then
              pf4 = fall(n,kk,4)*delz(n,kk)/rhoh2o/qrs(n,kk,4)
            else
              pf4 = d_zero
            end if
            mc2mo%remrat(j,i,k) = pf1 + pf2 + pf3 + pf4
            n = n + 1
          end do
        end do
      end do
      do i = ici1 , ici2
        do j = jci1 , jci2
          mc2mo%rembc(j,i,1) = d_zero
        end do
      end do
      do k = 2 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            mc2mo%rembc(j,i,k) = d_zero
            if ( mc2mo%remrat(j,i,k) > d_zero ) then
              do kk = 1 , k - 1
                qcw = mo2mc%qcn(j,i,k)
                mc2mo%rembc(j,i,k) = mc2mo%rembc(j,i,k) + & ![mm/hr]
                  mc2mo%remrat(j,i,kk) * qcw * &
                  (mo2mc%phs(j,i,k)-mo2mc%phs(j,i,k))*regrav
              end do
            end if
          end do
        end do
      end do
    end if

    n = 1
    do i = ici1 , ici2
      do j = jci1 , jci2
        totp = rain(n) + snow(n) + grpl(n) + hail(n)
        mc2mo%trrate(j,i) = totp * rdt
        mc2mo%rainnc(j,i) = mc2mo%rainnc(j,i) + totp
        mc2mo%snownc(j,i) = mc2mo%snownc(j,i) + snow(n)*rdt
        mc2mo%grplnc(j,i) = mc2mo%grplnc(j,i) + grpl(n)*rdt
        mc2mo%hailnc(j,i) = mc2mo%hailnc(j,i) + hail(n)*rdt
        mc2mo%lsmrnc(j,i) = mc2mo%lsmrnc(j,i) + mc2mo%trrate(j,i)
        n = n + 1
      end do
    end do

!    do n = is , ie
!      do k = 1 , kz
!        re_qc(k) = 2.51e-6_rkx
!        re_qi(k) = 10.01e-6_rkx
!        re_qs(k) = 25.e-6_rkx
!        t1d(k)  = t(n,k)
!        den1d(k)= den(n,k)
!        qc1d(k) = qci(n,k,1)
!        qi1d(k) = qci(n,k,2)
!        qs1d(k) = qrs(n,k,2)
!      end do
!      call effectrad_wsm7(t1d, qc1d, qi1d, qs1d, den1d,   &
!                          re_qc, re_qi, re_qs)
!      do k = 1 , kz
!        cloud_er(n,k) = max(2.51e-6_rkx,  min(re_qc(k),  50.e-6_rkx))
!        ice_er(n,k)   = max(10.01e-6_rkx, min(re_qi(k), 125.e-6_rkx))
!        snow_er(n,k)  = max(25.e-6_rkx,   min(re_qs(k), 999.e-6_rkx))
!      end do
!    end do
  end subroutine wsm7
  !
  ! This code is a 7-class mixed ice microphyiscs scheme (WSM7) of the
  ! single-moment microphyiscs (WSMMP,  Bae et al. 2018).
  ! The WSMMP assumes that ice nuclei
  ! number concentration is a function of temperature, and seperate assumption
  ! is developed, in which ice crystal number concentration is a function
  ! of ice amount. A theoretical background of the ice-microphysics and related
  ! processes in the WSMMPS are described in Hong et al. (2004).
  ! Production terms in the WSM6 scheme are described in Hong and Lim (2006).
  ! all units are in m.k.s. and source/sink terms in kgkg-1s-1.
  !
  !  WSM7 cloud scheme
  !
  !  Coded by Song-You Hong (Yonsei Univ.)
  !             Jimy Dudhia (NCAR) and Shu-Hua Chen (UC Davis)
  !             Summer 2002
  !
  !  Implemented by Song-You Hong (Yonsei Univ.) and Jimy Dudhia (NCAR)
  !             Summer 2003
  !
  !  Further modifications :
  !        semi-lagrangian sedimentation (JH,2010),Hong, Aug 2009
  !        ==> higher accuracy and efficient at lower resolutions
  !        reflectivity computation from Greg Thompson, Lim, Jun 2011
  !        ==> only diagnostic, but with removal of too large drops
  !        effective radius of hydrometeors, bae from Kiaps, Jan 2015
  !        ==> consistency in solar insolation of rrtmg radiation
  !
  !  Reference) Hong, Dudhia, Chen (hdc, 2004) Mon. Wea. Rev.
  !             Rutledge, Hobbs (rh83, 1983) J. Atmos. Sci.
  !             Hong and Lim (hl, 2006) J. Korean Meteor. Soc.
  !
  subroutine wsm72d(delt,ims,ime)
    implicit none
    real(rkx) , intent(in) :: delt
    integer(ik4) , intent(in) :: ims , ime

    real(rkx) , dimension(ims:ime,kz,4) :: falk , fall , work1 , &
      rslope2 , rslopeb , rslope , rslope3
    real(rkx) , dimension(ims:ime,kz) :: fallc , xl , cpm
    real(rkx) , dimension(ims:ime,kz) :: denfac , xni , denqrs1 , denqrs2
    real(rkx) , dimension(ims:ime,kz) :: denqci , n0sfac
    real(rkx) , dimension(ims:ime,kz) :: falkc , work1c , workr , workh , worka
    real(rkx) , dimension(ims:ime,kz) :: work2 , works
    real(rkx) , dimension(ims:ime) :: delqrs1 , delqrs2 , delqrs3 , &
      delqrs4 , delqi
    real(rkx) , dimension(ims:ime,kz) :: pigen , pidep , psdep , praut
    real(rkx) , dimension(ims:ime,kz) :: psaut , prevp , psevp , pracw
    real(rkx) , dimension(ims:ime,kz) :: psacw , psaci , pcond
    real(rkx) , dimension(ims:ime,kz) :: pgevp , phevp , pgdep , pvapg , &
      pvaph , phdep , pgaut , phaut , piacr , praci , &
      pracs , pracg , psacr , pgacw , pgaci , &
      pgacr , pgacs , paacw , phacw , phaci , &
      phacr , phacs , phacg , pgwet , phwet , &
      primh , psmlt , pgmlt , phmlt , pseml , &
      pgeml , pheml
    real(rkx) , dimension(ims:ime,kz) :: pgaci_w , phaci_w , qsum ,   &
                                         denqrs3 , denqrs4
    real(rkx) :: rdtcld
    real(rkx) :: supcol , supcolt , coeres , &
      supsat , dtcld , xmi , eacrs , satdt , vt2i , vt2r , vt2s , vt2g , &
      vt2h , acrfac , egi , ehi , qimax , diameter , xni0 , roqi0 , fallsum , &
      fallsum_qsi , fallsum_qg , fallsum_qh , xlwork2 , factor , source , &
      qval , xlf , pfrzdtc , pfrzdtr , supice , alpha2 , delta2 , delta3
    real(rkx) :: vt2ave
    real(rkx) :: rs0 , ghw1 , ghw2 , ghw3 , ghw4 , tr

    !real(rkx) :: tr , logtr
    ! variables for optimization
    real(rkx) :: temp
    integer(ik4) :: i , k , loop , loops , ifsat , nval

    nval = ime-ims+1
    !
    ! latent heat for phase changes and heat capacity. neglect the
    ! changes during microphysical process calculation
    ! emanuel(1994)
    !
    do k = 1 , kz
      do i = ims , ime
        cpm(i,k) = cpmcal(qv(i,k))
        xl(i,k) = wlh(t(i,k))
      end do
    end do
    do k = 1 , kz
      do i = ims , ime
        denfac(i,k) = sqrt(stdrho/den(i,k))
      end do
    end do
    !
    ! initialize the surface rain, snow , graupel , hail
    !
    do i = ims , ime
      rain(i) = d_zero
      snow(i) = d_zero
      grpl(i) = d_zero
      hail(i) = d_zero
    end do
    !
    ! compute the minor time steps.
    !
    loops = max(nint(delt/dtcldcr),1)
    dtcld = delt/real(loops,rkx)
    if ( delt <= dtcldcr ) then
      dtcld = delt
      loops = 1
    end if
    rdtcld = d_one/dtcld

    bigloop: &
    do loop = 1 , loops
      do k = 1 , kz
        do i = ims , ime
          tr = wattp/t(i,k)
          qs(i,k,1) = psat*exp(log(tr)*(xa))*exp(xb*(1.0_rkx-tr))
          qs(i,k,1) = min(qs(i,k,1),0.99_rkx*p(i,k))
          qs(i,k,1) = ep2 * qs(i,k,1) / (p(i,k) - qs(i,k,1))
          qs(i,k,1) = max(qs(i,k,1),minqq)
          rh(i,k,1) = max(qv(i,k) / qs(i,k,1),minqq)
          if ( t(i,k) < wattp ) then
            qs(i,k,2) = psat*exp(log(tr)*(xai))*exp(xbi*(1.0_rkx-tr))
          else
            qs(i,k,2) = psat*exp(log(tr)*(xa))*exp(xb*(1.0_rkx-tr))
          endif
          qs(i,k,2) = min(qs(i,k,2),0.99_rkx*p(i,k))
          qs(i,k,2) = ep2 * qs(i,k,2) / (p(i,k) - qs(i,k,2))
          qs(i,k,2) = max(qs(i,k,2),minqq)
          rh(i,k,2) = max(qv(i,k) / qs(i,k,2),minqq)
        end do
      end do
      !
      ! initialize the variables for microphysical physics
      !
      do k = 1 , kz
        do i = ims , ime
          prevp(i,k) = d_zero
          psdep(i,k) = d_zero
          pgdep(i,k) = d_zero
          phdep(i,k) = d_zero
          pvapg(i,k) = d_zero
          pvaph(i,k) = d_zero
          praut(i,k) = d_zero
          psaut(i,k) = d_zero
          pgaut(i,k) = d_zero
          phaut(i,k) = d_zero
          pracw(i,k) = d_zero
          praci(i,k) = d_zero
          piacr(i,k) = d_zero
          psaci(i,k) = d_zero
          psacw(i,k) = d_zero
          pracs(i,k) = d_zero
          pracg(i,k) = d_zero
          psacr(i,k) = d_zero
          pgacw(i,k) = d_zero
          paacw(i,k) = d_zero
          pgaci(i,k) = d_zero
          pgacr(i,k) = d_zero
          pgacs(i,k) = d_zero
          phacw(i,k) = d_zero
          phaci(i,k) = d_zero
          phacr(i,k) = d_zero
          phacs(i,k) = d_zero
          phacg(i,k) = d_zero
          pgwet(i,k) = d_zero
          phwet(i,k) = d_zero
          primh(i,k) = d_zero
          pigen(i,k) = d_zero
          pidep(i,k) = d_zero
          pcond(i,k) = d_zero
          psmlt(i,k) = d_zero
          pgmlt(i,k) = d_zero
          phmlt(i,k) = d_zero
          pseml(i,k) = d_zero
          pgeml(i,k) = d_zero
          pheml(i,k) = d_zero
          psevp(i,k) = d_zero
          pgevp(i,k) = d_zero
          phevp(i,k) = d_zero
          pgaci_w(i,k) = d_zero
          phaci_w(i,k) = d_zero
          falk(i,k,1) = d_zero
          falk(i,k,2) = d_zero
          falk(i,k,3) = d_zero
          falk(i,k,4) = d_zero
          fall(i,k,1) = d_zero
          fall(i,k,2) = d_zero
          fall(i,k,3) = d_zero
          fall(i,k,4) = d_zero
          fallc(i,k) = d_zero
          falkc(i,k) = d_zero
          xni(i,k) = 1.e3_rkx
        end do
      end do
      !
      ! ni: ice crystal number concentraiton   [hdc 5c]
      !
      do k = 1 , kz
        do i = ims , ime
          temp = den(i,k)*max(qci(i,k,2),minqq)
          temp = sqrt(sqrt(temp*temp*temp))
          xni(i,k) = min(max(5.38e7_rkx*temp,minni),maxni)
        end do
      end do
      !
      ! compute the fallout term:
      ! first, vertical terminal velocity for minor loops
      !
      call slope_wsm7(qrs,den,denfac,t, &
                      rslope,rslopeb,rslope2,rslope3,work1,ims,ime)
      do k = kz , 1, -1
        do i = ims , ime
          workr(i,k) = work1(i,k,1)
          workh(i,k) = work1(i,k,4)
          qsum(i,k) = max( (qrs(i,k,2)+qrs(i,k,3)), 1.e-15_rkx)
          if ( qsum(i,k) > 1.1e-15_rkx ) then
            worka(i,k) = (work1(i,k,2)*qrs(i,k,2) + &
                          work1(i,k,3)*qrs(i,k,3))/qsum(i,k)
          else
            worka(i,k) = 0.0_rkx
          end if
          denqrs1(i,k) = den(i,k)*qrs(i,k,1)
          denqrs2(i,k) = den(i,k)*qrs(i,k,2)
          denqrs3(i,k) = den(i,k)*qrs(i,k,3)
          denqrs4(i,k) = den(i,k)*qrs(i,k,4)
          if( qrs(i,k,1) <= d_zero ) workr(i,k) = d_zero
          if( qrs(i,k,4) <= d_zero ) workh(i,k) = d_zero
        enddo
      enddo

      call nislfv_rain_plm(nval,den,denfac,t,delz, &
                           workr,denqrs1,delqrs1,dtcld,1,1)
      call nislfv_rain_plm6(nval,den,denfac,t,delz, &
                            worka,denqrs2,denqrs3,delqrs2,delqrs3,dtcld,1,1)
      call nislfv_rain_plm(nval,den,denfac,t,delz, &
                           workh,denqrs4,delqrs4,dtcld,2,1)
      do k = 1 , kz
        do i = ims , ime
          qrs(i,k,1) = max(denqrs1(i,k)/den(i,k),d_zero)
          qrs(i,k,2) = max(denqrs2(i,k)/den(i,k),d_zero)
          qrs(i,k,3) = max(denqrs3(i,k)/den(i,k),d_zero)
          qrs(i,k,4) = max(denqrs4(i,k)/den(i,k),d_zero)
          fall(i,k,1) = denqrs1(i,k)*workr(i,k)/delz(i,k)
          fall(i,k,2) = denqrs2(i,k)*works(i,k)/delz(i,k)
          fall(i,k,3) = denqrs3(i,k)*workr(i,k)/delz(i,k)
          fall(i,k,4) = denqrs4(i,k)*works(i,k)/delz(i,k)
        end do
      end do
      do i = ims , ime
        fall(i,1,1) = delqrs1(i)/delz(i,1)*rdtcld
        fall(i,1,2) = delqrs2(i)/delz(i,1)*rdtcld
        fall(i,1,3) = delqrs3(i)/delz(i,1)*rdtcld
        fall(i,1,4) = delqrs4(i)/delz(i,1)*rdtcld

      end do
      call slope_wsm7(qrs,den,denfac,t, &
                      rslope,rslopeb,rslope2,rslope3,work1,ims,ime)
      do k = kz , 1 , -1
        do i = ims , ime
          supcol = tzero - t(i,k)
          n0sfac(i,k) = max(min(exp(alpha*supcol),n0smax/n0s),d_one)
          if ( t(i,k) > tzero ) then
            xlf = wlhf
            !
            ! psmlt: melting of snow [hl a33] [rh83 a25]
            !       (t>t0: s->r)
            !
            work2(i,k) = venfac(p(i,k),t(i,k),den(i,k))
            if ( qrs(i,k,2) > d_zero ) then
              coeres = rslope2(i,k,2)*sqrt(rslope(i,k,2)*rslopeb(i,k,2))
              psmlt(i,k) = xka(t(i,k),den(i,k))/xlf * &
                (tzero-t(i,k))*halfpi*n0sfac(i,k)*(precs1*rslope2(i,k,2) + &
                precs2*work2(i,k)*coeres)/den(i,k)
              psmlt(i,k) = min(max(psmlt(i,k)*dtcld,-qrs(i,k,2)),d_zero)
              qrs(i,k,2) = qrs(i,k,2) + psmlt(i,k)
              qrs(i,k,1) = qrs(i,k,1) - psmlt(i,k)
              t(i,k) = t(i,k) + xlf/cpm(i,k)*psmlt(i,k)
            end if
          end if
          if ( t(i,k) > tzero .and. qrs(i,k,3) > d_zero ) then
            !
            ! pgmlt: melting of graupel [HL A23]  [LFO 47]
            !       (T>T0: G->R)
            !
            coeres = rslope2(i,k,3)*sqrt(rslope(i,k,3)*rslopeb(i,k,3))
            pgmlt(i,k) = xka(t(i,k),den(i,k))/xlf*(tzero-t(i,k)) * &
              (precg1*rslope2(i,k,3)+precg2*work2(i,k)*coeres)/den(i,k)
            pgmlt(i,k) = min(max(pgmlt(i,k)*dtcld,-qrs(i,k,3)),d_zero)
            qrs(i,k,3) = qrs(i,k,3) + pgmlt(i,k)
            qrs(i,k,1) = qrs(i,k,1) - pgmlt(i,k)
            t(i,k) = t(i,k) + xlf/cpm(i,k)*pgmlt(i,k)
          end if
          if ( t(i,k) > tzero .and. qrs(i,k,4) > d_zero ) then
            !
            ! phmlt: melting of hail [BHT A22]
            !       (T>T0: H->R)
            !
            coeres = rslope2(i,k,4)*sqrt(rslope(i,k,4)*rslopeb(i,k,4))
            phmlt(i,k) = xka(t(i,k),den(i,k))/xlf*(tzero-t(i,k)) * &
              (prech1*rslope2(i,k,4)+prech2*work2(i,k)*coeres)/den(i,k)
            phmlt(i,k) = min(max(phmlt(i,k)*dtcld,-qrs(i,k,4)),d_zero)
            qrs(i,k,4) = qrs(i,k,4) + phmlt(i,k)
            qrs(i,k,1) = qrs(i,k,1) - phmlt(i,k)
            t(i,k) = t(i,k) + xlf/cpm(i,k)*phmlt(i,k)
          end if
        end do
      end do
      !
      ! vice [ms-1] : fallout of ice crystal [hdc 5a]
      !
      do k = kz , 1 , -1
        do i = ims , ime
          if ( qci(i,k,2) <= d_zero ) then
            work1c(i,k) = d_zero
          else
            xmi = den(i,k)*qci(i,k,2)/xni(i,k)
            diameter  = max(min(dicon*sqrt(xmi),dimax), 1.e-25_rkx)
            work1c(i,k) = 1.49e4_rkx*exp(log(diameter)*(1.31_rkx))
          end if
        end do
      end do
      !
      ! forward semi-laglangian scheme (jh), pcm (piecewise constant), (linear)
      !
      do k = kz , 1 , -1
        do i = ims , ime
          denqci(i,k) = den(i,k)*qci(i,k,2)
        end do
      end do
      call nislfv_rain_plm(nval,den,denfac,t,delz, &
                           work1c,denqci,delqi,dtcld,1,0)
      do k = 1 , kz
        do i = ims , ime
          qci(i,k,2) = max(denqci(i,k)/den(i,k),d_zero)
        end do
      end do
      do i = ims , ime
        fallc(i,1) = delqi(i)/delz(i,1)*rdtcld
      end do
      !
      ! rain (unit is mm/sec;kgm-2s-1: /1000*delt ===> m)==> mm for wrf
      !
      do i = ims , ime
        fallsum = fall(i,1,1)+fall(i,1,2)+fall(i,1,3)+fall(i,1,4)+fallc(i,1)
        fallsum_qsi = fall(i,1,2)+fallc(i,1)
        fallsum_qg = fall(i,1,3)
        fallsum_qh = fall(i,1,4)
        if ( fallsum > d_zero ) then
          rain(i) = fallsum*delz(i,1)/rhoh2o*dtcld*1000.0_rkx + rain(i)
        end if
        if ( fallsum_qsi > d_zero ) then
          snow(i) = fallsum_qsi*delz(i,1)/rhoh2o*dtcld*1000.0_rkx + snow(i)
        end if
        if ( fallsum_qg > d_zero ) then
          grpl(i) = fallsum_qg*delz(i,1)/rhoh2o*dtcld*1000.0_rkx + grpl(i)
        end if
        if ( fallsum_qh > d_zero ) then
          hail(i) = fallsum_qh*delz(i,1)/rhoh2o*dtcld*1000.0_rkx + hail(i)
        end if
      end do
      !
      ! pimlt: instantaneous melting of cloud ice [hl a47] [rh83 a28]
      !       (t>t0: i->c)
      !
      do k = 1 , kz
        do i = ims , ime
          supcol = tzero-t(i,k)
          xlf = max(wlhs-xl(i,k),d_zero)
          if ( supcol < d_zero ) xlf = wlhf
          if ( supcol < d_zero .and. qci(i,k,2) > d_zero ) then
            qci(i,k,1) = qci(i,k,1) + qci(i,k,2)
            t(i,k) = t(i,k) - xlf/cpm(i,k)*qci(i,k,2)
            qci(i,k,2) = d_zero
          end if
          !
          ! pihmf: homogeneous freezing of cloud water below -40c [hl a45]
          !        (t<-40c: c->i)
          !
          if ( supcol > 40.0_rkx .and. qci(i,k,1) > d_zero ) then
            qci(i,k,2) = qci(i,k,2) + qci(i,k,1)
            t(i,k) = t(i,k) + xlf/cpm(i,k)*qci(i,k,1)
            qci(i,k,1) = d_zero
          end if
          !
          ! pihtf: heterogeneous freezing of cloud water [hl a44]
          !        (t0>t>-40c: c->i)
          !
          if ( supcol > d_zero .and. qci(i,k,1) > d_zero ) then
            if ( qci(i,k,1) > 1.0e-9_rkx ) then
              supcolt = min(supcol,50.0_rkx)
              pfrzdtc = min(pfrz1*(exp(pfrz2*supcolt)-d_one) * &
                 den(i,k)/rhoh2o/xncr*qci(i,k,1)*qci(i,k,1)*dtcld,qci(i,k,1))
            else
              pfrzdtc = 0.9_rkx*qci(i,k,1)
            end if
            qci(i,k,2) = qci(i,k,2) + pfrzdtc
            t(i,k) = t(i,k) + xlf/cpm(i,k)*pfrzdtc
            qci(i,k,1) = qci(i,k,1)-pfrzdtc
          end if
          !
          ! pgfrz: freezing of rain water [hl a20] [lfo 45]
          !        (t<t0, r->g)
          !
          if ( supcol > d_zero .and. qrs(i,k,1) > d_zero ) then
            supcolt = min(supcol,50.0_rkx)
            temp = rslope3(i,k,1)
            temp = temp*temp*rslope(i,k,1)
            pfrzdtr = min(20.0_rkx*pisqr*pfrz1*n0r*rhoh2o/den(i,k) * &
                (exp(pfrz2*supcolt)-d_one)*temp*dtcld, qrs(i,k,1))
            qrs(i,k,3) = qrs(i,k,3) + pfrzdtr
            t(i,k) = t(i,k) + xlf/cpm(i,k)*pfrzdtr
            qrs(i,k,1) = qrs(i,k,1)-pfrzdtr
          end if
        end do
      end do
      !
      ! update the slope parameters for microphysics computation
      !
      call slope_wsm7(qrs,den,denfac,t, &
                      rslope,rslopeb,rslope2,rslope3,work1,ims,ime)
      ! work1:  the thermodynamic term in the denominator associated with
      !         heat conduction and vapor diffusion
      !         (ry88, y93, h85)
      ! work2: parameter associated with the ventilation effects(y93)
      !
      do k = 1 , kz
        do i = ims , ime
          work1(i,k,1) = diffac(xl(i,k),p(i,k),t(i,k),den(i,k),qs(i,k,1))
          work1(i,k,2) = diffac(wlhs,p(i,k),t(i,k),den(i,k),qs(i,k,2))
          work2(i,k) = venfac(p(i,k),t(i,k),den(i,k))
        end do
      end do
      !
      ! warm rain processes
      !
      ! - follows the processes in rh83 and lfo except for autoconcersion
      !
      do k = 1 , kz
        do i = ims , ime
          supsat = max(qv(i,k),minqq)-qs(i,k,1)
          satdt = supsat*rdtcld
          !
          ! praut: auto conversion rate from cloud to rain [hdc 16]
          !        (c->r)
          !
          if ( qci(i,k,1) > qc0 ) then
            praut(i,k) = qck1*exp(log(qci(i,k,1))*((7.0_rkx/3.0_rkx)))
            praut(i,k) = min(praut(i,k),qci(i,k,1)*rdtcld)
          end if
          !
          ! pracw: accretion of cloud water by rain [hl a40] [lfo 51]
          !        (c->r)
          !
          if ( qrs(i,k,1) > qrsmin .and. qci(i,k,1) > qcimin ) then
            pracw(i,k) = min(pacrr*rslope3(i,k,1)*rslopeb(i,k,1) * &
                             qci(i,k,1)*denfac(i,k),qci(i,k,1)*rdtcld)
          end if
          !
          ! prevp: evaporation/condensation rate of rain [hdc 14]
          !        (v->r or r->v)
          !
          if ( qrs(i,k,1) > d_zero ) then
            coeres = rslope2(i,k,1)*sqrt(rslope(i,k,1)*rslopeb(i,k,1))
            prevp(i,k) = (rh(i,k,1)-d_one)*(precr1*rslope2(i,k,1) + &
                         precr2*work2(i,k)*coeres)/work1(i,k,1)
            if ( prevp(i,k) < d_zero ) then
              prevp(i,k) = max(prevp(i,k),-qrs(i,k,1)*rdtcld)
              prevp(i,k) = max(prevp(i,k),satdt/2.0_rkx)
            else
              prevp(i,k) = min(prevp(i,k),satdt/2.0_rkx)
            end if
          end if
        end do
      end do
      !
      ! cold rain processes
      !
      ! - follows the revised ice microphysics processes in hdc
      ! - the processes same as in rh83 and rh84  and lfo behave
      !   following ice crystal hapits defined in hdc, inclduing
      !   intercept parameter for snow (n0s), ice crystal number
      !   concentration (ni), ice nuclei number concentration
      !   (n0i), ice diameter (d)
      !
      do k = 1 , kz
        do i = ims , ime
          supcol = tzero-t(i,k)
          n0sfac(i,k) = max(min(exp(alpha*supcol),n0smax/n0s),d_one)
          supsat = max(qv(i,k),minqq)-qs(i,k,2)
          satdt = supsat*rdtcld
          ifsat = 0
          !
          ! ni: ice crystal number concentraiton   [hdc 5c]
          !
          temp = den(i,k)*max(qci(i,k,2),qcimin)
          temp = sqrt(sqrt(temp*temp*temp))
          xni(i,k) = min(max(5.38e7_rkx*temp,minni),maxni)
          eacrs = exp(0.07_rkx*(-supcol))
          xmi = den(i,k)*qci(i,k,2)/xni(i,k)
          diameter  = min(dicon * sqrt(xmi),dimax)
          vt2i = 1.49e4_rkx*diameter**1.31_rkx
          vt2r=pvtr*rslopeb(i,k,1)*denfac(i,k)
          vt2s=pvts*rslopeb(i,k,2)*denfac(i,k)
          vt2g=pvtg*rslopeb(i,k,3)*denfac(i,k)
          vt2h=pvth*rslopeb(i,k,4)*denfac(i,k)
          qsum(i,k) = max( (qrs(i,k,2) + qrs(i,k,3)), 1.e-15_rkx)
          if ( qsum(i,k) > 1.1e-15_rkx) then
          vt2ave=(vt2s*qrs(i,k,2)+vt2g*qrs(i,k,3))/(qsum(i,k))
          else
          vt2ave=d_zero
          end if
          if ( supcol > d_zero .and. qci(i,k,2) > qcimin ) then
            if ( qrs(i,k,1) > qrsmin ) then
              !
              ! praci: accretion of cloud ice by rain [lf0 25]
              !        (t<t0: i->r)
              !
              acrfac = d_two*rslope3(i,k,1) + &
                       d_two*diameter*rslope2(i,k,1) + &
                       diameter**2*rslope(i,k,1)
              praci(i,k) = mathpi*qci(i,k,2)*n0r*abs(vt2r-vt2i)*acrfac*d_rfour
              ! reduce collection efficiency (suggested by B. Wilt)
              praci(i,k) = praci(i,k) * &
                min(max(d_zero,qrs(i,k,1)/qci(i,k,2)),d_one)**2
              praci(i,k) = min(praci(i,k),qci(i,k,2)/dtcld)
              !
              ! piacr: Accretion of rain by cloud ice [HL A19] [LFO 26]
              !        (t<t0: r->s or r->g)
              !
              piacr(i,k) = mathpi**2*avtr*n0r*rhoh2o*xni(i,k)*denfac(i,k) * &
                           g6pbr*rslope3(i,k,1)*rslope3(i,k,1)            * &
                           rslopeb(i,k,1)/24.0_rkx/den(i,k)
             ! reduce collection efficiency (suggested by B. Wilt)
              piacr(i,k) = piacr(i,k) * &
                min(max(d_zero,qci(i,k,2)/qrs(i,k,1)),d_one)**2
              piacr(i,k) = min(piacr(i,k),qrs(i,k,1)/dtcld)
            end if
            !
            ! psaci: Accretion of cloud ice by snow [HDC 10]
            !        (t<t0: i->s)
            !
            if( qrs(i,k,2) > qrsmin ) then
              acrfac = d_two*rslope3(i,k,2)+d_two*diameter*rslope2(i,k,2)  + &
                       diameter**2*rslope(i,k,2)
              psaci(i,k) = mathpi*qci(i,k,2)*eacrs*n0s*n0sfac(i,k) * &
                           abs(vt2ave-vt2i)*acrfac*d_rfour
              psaci(i,k) = min(psaci(i,k),qci(i,k,2)/dtcld)
            end if
            !
            ! pgaci: Accretion of cloud ice by graupel [HL A17] [LFO 41]
            !        (t<t0: i->g)
            !
            if( qrs(i,k,3) > qrsmin ) then
              egi = exp(0.07_rkx*(-supcol))
              acrfac = d_two*rslope3(i,k,3)+d_two*diameter*rslope2(i,k,3) + &
                       diameter**2*rslope(i,k,3)
              pgaci(i,k) = mathpi*egi*qci(i,k,2)*n0g * &
                abs(vt2ave-vt2i)*acrfac*d_rfour
              pgaci(i,k) = min(pgaci(i,k),qci(i,k,2)/dtcld)
            end if
            !
            ! phaci: Accretion of cloud ice by hail [BHT ]
            !        (T<T0: I->H)
            !
            if( qrs(i,k,4) > qrsmin ) then
              ehi = exp(0.07_rkx*(-supcol))
              acrfac = d_two*rslope3(i,k,4)+d_two*diameter*rslope2(i,k,4) + &
                       diameter**2*rslope(i,k,4)
              phaci(i,k) = mathpi*ehi*qci(i,k,2)*n0h*abs(vt2h-vt2i) * &
                acrfac*d_rfour
              phaci(i,k) = min(phaci(i,k),qci(i,k,2)/dtcld)
            end if
          end if
          !
          ! psacw: accretion of cloud water by snow  [hl a7] [lfo 24]
          !        (t<t0: c->s, and t>=t0: c->r)
          !
          if ( qrs(i,k,2) > qrsmin .and. qci(i,k,1) > qcimin ) then
            psacw(i,k) = min(pacrc*n0sfac(i,k)*rslope3(i,k,2)*rslopeb(i,k,2)* &
                             min(max(d_zero,qrs(i,k,2)/qci(i,k,1)),d_one)**2* &
                             qci(i,k,1)*denfac(i,k),qci(i,k,1)/dtcld)
          end if
          !
          ! pgacw: Accretion of cloud water by graupel [HL A6] [LFO 40]
          !        (T<T0: C->G, and T>=T0: C->R)
          !
          if( qrs(i,k,3) > qrsmin .and. qci(i,k,1) > minqq) then
            ! reduce collection efficiency (suggested by B. Wilt)
            pgacw(i,k) = min(pacrg*rslope3(i,k,3)*rslopeb(i,k,3)          * &
                         min(max(d_zero,qrs(i,k,3)/qci(i,k,1)),d_one)**2  * &
                         qci(i,k,1)*denfac(i,k),qci(i,k,1)/dtcld)
          end if
          !
          ! paacw: Accretion of cloud water by averaged snow/graupel
          !        (T<T0: C->G or S, and T>=T0: C->R)
          !
          if( qsum(i,k) > 1.1e-15_rkx ) then
            paacw(i,k) = (qrs(i,k,2)*psacw(i,k)+qrs(i,k,3) * &
              pgacw(i,k))/(qsum(i,k))
          end if
          !
          ! phacw: Accretion of cloud water by hail [BHT A08]
          !        (T<T0: C->H, and T>=T0: C->R)
          !
          if ( qrs(i,k,4) > qrsmin .and. qci(i,k,1) > minqq) then
            ! reduce collection efficiency (suggested by B. Wilt)
            phacw(i,k) = min(pacrh*rslope3(i,k,4)*rslopeb(i,k,4)          * &
                         min(max(d_zero,qrs(i,k,4)/qci(i,k,1)),d_one)**2  * &
                         qci(i,k,1)*denfac(i,k),qci(i,k,1)/dtcld)
          end if
          !
          ! pracs: Accretion of snow by rain [HL A11] [LFO 27]
          !         (T<T0: S->G)
          !
          if ( qrs(i,k,2) > qrsmin .and. qrs(i,k,1) > qrsmin) then
            if ( supcol > d_zero ) then
              acrfac = d_five*rslope3(i,k,2)*rslope3(i,k,2)*rslope(i,k,1)  + &
                       d_two*rslope3(i,k,2)*rslope2(i,k,2)*rslope2(i,k,1)  + &
                       d_half*rslope2(i,k,2)*rslope2(i,k,2)*rslope3(i,k,1)
              pracs(i,k) = mathpi**2*n0r*n0s*n0sfac(i,k) * &
                abs(vt2r-vt2ave)*(dens/den(i,k))*acrfac
              ! reduce collection efficiency (suggested by B. Wilt)
              pracs(i,k) = pracs(i,k)*min(max(d_zero,qrs(i,k,1) / &
                qrs(i,k,2)),d_one)**2
              pracs(i,k) = min(pracs(i,k),qrs(i,k,2)/dtcld)
            end if
            !
            ! psacr: Accretion of rain by snow [HL A10] [LFO 28]
            !         (T<T0:R->S or R->G) (T>=T0: enhance melting of snow)
            !
            acrfac = d_five*rslope3(i,k,1)*rslope3(i,k,1)*rslope(i,k,2)  + &
                     d_two*rslope3(i,k,1)*rslope2(i,k,1)*rslope2(i,k,2)  + &
                     d_half*rslope2(i,k,1)*rslope2(i,k,1)*rslope3(i,k,2)
            psacr(i,k) = mathpi**2*n0r*n0s*n0sfac(i,k) * &
              abs(vt2ave-vt2r)*(rhoh2o/den(i,k))*acrfac
            ! reduce collection efficiency (suggested by B. Wilt)
            psacr(i,k) = psacr(i,k) * &
              min(max(d_zero,qrs(i,k,2)/qrs(i,k,1)),d_one)**2
            psacr(i,k) = min(psacr(i,k),qrs(i,k,1)/dtcld)
          end if
          !
          ! pracg: Accretion of graupel by rain [BHT A17]
          !         (T<T0: G->H)
          !
          if ( qrs(i,k,3) > qrsmin .and. qrs(i,k,1) > qrsmin) then
            if ( supcol > 0 ) then
              acrfac = d_five*rslope3(i,k,3)*rslope3(i,k,3)*rslope(i,k,1)  + &
                       d_two*rslope3(i,k,3)*rslope2(i,k,3)*rslope2(i,k,1)  + &
                       d_half*rslope2(i,k,3)*rslope2(i,k,3)*rslope3(i,k,1)
              pracg(i,k) = mathpi**2*n0r*n0g * &
                abs(vt2r-vt2ave)*(deng/den(i,k))*acrfac
              ! reduce collection efficiency (suggested by B. Wilt)
              pracg(i,k) = pracg(i,k) * &
                min(max(d_zero,qrs(i,k,1)/qrs(i,k,3)),d_one)**2
              pracg(i,k) = min(pracg(i,k),qrs(i,k,3)/dtcld)
            end if
            !
            ! pgacr: Accretion of rain by graupel [HL A12] [LFO 42]
            !         (T<T0: R->G) (T>=T0: enhance melting of graupel)
            !
            acrfac = d_five*rslope3(i,k,1)*rslope3(i,k,1)*rslope(i,k,3) + &
                     d_two*rslope3(i,k,1)*rslope2(i,k,1)*rslope2(i,k,3) + &
                     d_half*rslope2(i,k,1)*rslope2(i,k,1)*rslope3(i,k,3)
            pgacr(i,k) = mathpi**2*n0r*n0g*abs(vt2ave-vt2r) * &
              (rhoh2o/den(i,k))*acrfac
            ! reduce collection efficiency (suggested by B. Wilt)
            pgacr(i,k) = pgacr(i,k)*min(max(d_zero,qrs(i,k,3) / &
              qrs(i,k,1)),d_one)**2
            pgacr(i,k) = min(pgacr(i,k),qrs(i,k,1)/dtcld)
          end if
          !
          ! pgacs: Accretion of snow by graupel [HL A13] [LFO 29]
          !        (S->G): This process is eliminated in V3.0 with the
          !        new combined snow/graupel fall speeds
          !
          if ( qrs(i,k,3) > qrsmin .and. qrs(i,k,2) > qrsmin) then
            pgacs(i,k) = 0.0_rkx
          end if
          !
          ! phacr: Accretion of rain by hail [BHT A13]
          !         (T<T0: R->H) (T>=T0: enhance melting of hail)
          !
          if ( qrs(i,k,4) > qrsmin .and. qrs(i,k,1) > qrsmin) then
            acrfac = d_five*rslope3(i,k,1)*rslope3(i,k,1)*rslope(i,k,4)  + &
                     d_two*rslope3(i,k,1)*rslope2(i,k,1)*rslope2(i,k,4)  + &
                     d_half*rslope2(i,k,1)*rslope2(i,k,1)*rslope3(i,k,4)
            phacr(i,k) = mathpi**2*n0r*n0h * &
              abs(vt2h-vt2r)*(rhoh2o/den(i,k))*acrfac
            ! reduce collection efficiency (suggested by B. Wilt)
            phacr(i,k) = phacr(i,k)*min(max(d_zero,qrs(i,k,4) / &
              qrs(i,k,1)),d_one)**2
            phacr(i,k) = min(phacr(i,k),qrs(i,k,1)/dtcld)
          end if
          !
          ! phacs: Accretion of snow by hail [BHT A14]
          !         (T<T0: S->H)
          !
          if ( qrs(i,k,4) > qrsmin .and. qrs(i,k,2) > qrsmin) then
            acrfac = d_five*rslope3(i,k,2)*rslope3(i,k,2)*rslope(i,k,4) + &
                     d_two*rslope3(i,k,2)*rslope2(i,k,2)*rslope2(i,k,4) + &
                     d_half*rslope2(i,k,2)*rslope2(i,k,2)*rslope3(i,k,4)
            phacs(i,k) = mathpi**2*eachs*n0s*n0sfac(i,k)*n0h * &
              abs(vt2h-vt2ave)*(dens/den(i,k))*acrfac
            phacs(i,k) = min(phacs(i,k),qrs(i,k,2)/dtcld)
          end if
          !
          ! phacg: Accretion of snow by hail [BHT A15]
          !         (T<T0: G->H)
          !
          if(qrs(i,k,4) > qrsmin.and.qrs(i,k,3) > qrsmin) then
            acrfac = d_five*rslope3(i,k,3)*rslope3(i,k,3)*rslope(i,k,4) + &
                     d_two*rslope3(i,k,3)*rslope2(i,k,3)*rslope2(i,k,4) + &
                     d_half*rslope2(i,k,3)*rslope2(i,k,3)*rslope3(i,k,4)
            phacg(i,k) = mathpi**2*eachg*n0g*n0h * &
              abs(vt2h-vt2ave)*(deng/den(i,k))*acrfac
            phacg(i,k) = min(phacg(i,k),qrs(i,k,3)/dtcld)
          end if
          !
          ! pgwet: wet growth of graupel [LFO 43]
          !
          !
          rs0 = psat*exp(log(wattp/tzero)*xa)*exp(xb*(d_one-wattp/tzero))
          rs0 = min(rs0,0.99_rkx*p(i,k))
          rs0 = ep2*rs0/(p(i,k)-rs0)
          rs0 = max(rs0,minqq)
          ghw1 = den(i,k)*wlhv*diffus(t(i,k),p(i,k))*(rs0-qv(i,k)) - &
                 xka(t(i,k),den(i,k))*(-supcol)
          ghw2 = den(i,k)*(wlhf+cpw*(-supcol))
          ghw3 = venfac(p(i,k),t(i,k),den(i,k)) * &
            sqrt(sqrt(egrav*den(i,k)/stdrho))
          ghw4 = den(i,k)*(wlhf-cpw*supcol+cpi*supcol)
          if ( qrs(i,k,3) > qrsmin ) then
            if( pgaci(i,k) > d_zero ) then
              egi = exp(0.07_rkx*(-supcol))
              pgaci_w(i,k) = pgaci(i,k)/egi
            else
              pgaci_w(i,k) = d_zero
            end if
            pgwet(i,k) = ghw1/ghw2*(precg1*rslope2(i,k,3)   + &
                         precg3*ghw3*rslope(i,k,4)**(2.75_rkx)  + &
                         ghw4*(pgaci_w(i,k)+pgacs(i,k)))

            pgwet(i,k) = max(pgwet(i,k), d_zero)
          end if
          !
          ! phwet: wet growth of hail [LFO 43]
          !
          !
          if ( qrs(i,k,4) > qrsmin ) then
            if ( phaci(i,k) > d_zero ) then
              ehi = exp(0.07_rkx*(-supcol))
              phaci_w(i,k) = phaci(i,k)/ehi
            else
              phaci_w(i,k) = d_zero
            end if
          end if
          phwet(i,k) = ghw1/ghw2*(prech1*rslope2(i,k,4)   + &
                       prech3*ghw3*rslope(i,k,4)**(2.75_rkx)  + &
                       ghw4*(phaci_w(i,k)+phacs(i,k)))
          phwet(i,k) = max(phwet(i,k), d_zero)
          !
          if ( phacw(i,k)+phacr(i,k) < 0.95_rkx*phwet(i,k) ) then
            phaci(i,k) = d_zero
            phacs(i,k) = d_zero
            phacg(i,k) = d_zero
          end if
          !
          if ( supcol <= d_zero ) then
            xlf = wlhf
          !
          ! pseml: Enhanced melting of snow by accretion of water [HL A34]
          !        (T>=T0: S->R)
          !
            if ( qrs(i,k,2) > d_zero ) then
              pseml(i,k) = min(max(cpw*supcol * &
                (paacw(i,k)+psacr(i,k))/xlf,-qrs(i,k,2)/dtcld),d_zero)
            end if
          !
          ! pgeml: Enhanced melting of graupel by accretion of water
          ! [HL A24] [RH84 A21-A22]
          !        (T>=T0: G->R)
          !
            if ( qrs(i,k,3) > d_zero ) then
              pgeml(i,k) = min(max(cpw*supcol * &
                (paacw(i,k)+pgacr(i,k))/xlf,-qrs(i,k,3)/dtcld),d_zero)
            end if
          !
          ! pheml: Enhanced melting of hail by accretion of water [BHT A23]
          !        (T>=T0: H->R)
          !
            if ( qrs(i,k,4) > d_zero ) then
              pheml(i,k) = min(max(cpw*supcol * &
                (phacw(i,k)+phacr(i,k))/xlf,-qrs(i,k,4)/dtcld),d_zero)
            end if
          end if
          !
         if ( supcol > d_zero ) then
            !
            ! pidep: deposition/sublimation rate of ice [hdc 9]
            !       (t<t0: v->i or i->v)
            !
           if ( qci(i,k,2) > d_zero .and. ifsat /= 1 ) then
             pidep(i,k) = d_four*diameter*xni(i,k) * &
                          (rh(i,k,2)-d_one)/work1(i,k,2)
             supice = satdt-prevp(i,k)
             if ( pidep(i,k) < d_zero ) then
               pidep(i,k) = max(max(pidep(i,k),satdt*d_half),supice)
               pidep(i,k) = max(pidep(i,k),-qci(i,k,2)*rdtcld)
             else
                pidep(i,k) = min(min(pidep(i,k),satdt*d_half),supice)
              end if
              if ( abs(prevp(i,k)+pidep(i,k)) >= abs(satdt) ) ifsat = 1
            end if
            !
            ! psdep: deposition/sublimation rate of snow [hdc 14]
            !        (v->s or s->v)
            !
            if ( qrs(i,k,2) > d_zero .and. ifsat /= 1 ) then
              coeres = rslope2(i,k,2)*sqrt(rslope(i,k,2)*rslopeb(i,k,2))
              psdep(i,k) = (rh(i,k,2)-d_one)*n0sfac(i,k) * &
                           (precs1*rslope2(i,k,2) + &
                            precs2*work2(i,k)*coeres)/work1(i,k,2)
              supice = satdt-prevp(i,k)-pidep(i,k)
              if ( psdep(i,k) < d_zero ) then
                psdep(i,k) = max(psdep(i,k),-qrs(i,k,2)*rdtcld)
                psdep(i,k) = max(max(psdep(i,k),satdt*d_half),supice)
              else
                psdep(i,k) = min(min(psdep(i,k),satdt*d_half),supice)
              end if
              if ( abs(prevp(i,k)+pidep(i,k)+psdep(i,k)) >= abs(satdt) ) then
                ifsat = 1
              end if
            end if
              !
              ! pgdep: deposition/sublimation rate of graupel [HL A21] [LFO 46]
              !        (T<T0: V->G or G->V)
              !
            if( qrs(i,k,3) > d_zero .and. ifsat /= 1 ) then
              coeres = rslope2(i,k,3)*sqrt(rslope(i,k,3)*rslopeb(i,k,3))
              pgdep(i,k) = (rh(i,k,2)-d_one)*(precg1*rslope2(i,k,3) + &
                            precg2*work2(i,k)*coeres)/work1(i,k,2)
              supice = satdt-prevp(i,k)-pidep(i,k)-psdep(i,k)
              if( pgdep(i,k) < d_zero ) then
                pgdep(i,k) = max(pgdep(i,k),-qrs(i,k,3)/dtcld)
                pgdep(i,k) = max(max(pgdep(i,k),satdt*d_half),supice)
              else
                pgdep(i,k) = min(min(pgdep(i,k),satdt*d_half),supice)
              end if
              if ( abs(prevp(i,k)+pidep(i,k) + &
                       psdep(i,k)+pgdep(i,k)) >= abs(satdt)) ifsat = 1
            end if
            !
            ! phdep: deposition/sublimation rate of hail [BHT A19]
            !        (T<T0: V->H or H->V)
            !
            if ( qrs(i,k,4) > d_zero .and. ifsat /= 1 ) then
              coeres = rslope2(i,k,4)*sqrt(rslope(i,k,4)*rslopeb(i,k,4))
              phdep(i,k) = (rh(i,k,2)-d_one)*(prech1*rslope2(i,k,4) + &
                           prech2*work2(i,k)*coeres)/work1(i,k,2)
              supice = satdt-prevp(i,k)-pidep(i,k)-psdep(i,k)-pgdep(i,k)
              if ( phdep(i,k) < d_zero ) then
                phdep(i,k) = max(phdep(i,k),-qrs(i,k,4)/dtcld)
                phdep(i,k) = max(max(phdep(i,k),satdt*d_half),supice)
              else
                phdep(i,k) = min(min(phdep(i,k),satdt*d_half),supice)
              end if
              if ( abs(prevp(i,k)+pidep(i,k)+psdep(i,k) + &
                       pgdep(i,k)+phdep(i,k)) >= abs(satdt)) ifsat = 1
            end if
           !
           ! pigen: generation(nucleation) of ice from vapor [hl a50] [hdc 7-8]
           !       (t<t0: v->i)
           !
           if ( supsat > d_zero .and. ifsat /= 1 ) then
             supice = satdt-prevp(i,k)-pidep(i,k)-psdep(i,k)
             xni0 = minni*exp(0.1_rkx*supcol)
             roqi0 = 4.92e-11_rkx*exp(log(xni0)*(1.33_rkx))
             pigen(i,k) = max(d_zero, &
                           (roqi0/den(i,k)-max(qci(i,k,2),d_zero))*rdtcld)
             pigen(i,k) = min(min(pigen(i,k),satdt),supice)
           end if
           !
           ! psaut: conversion(aggregation) of ice to snow [hdc 12]
           !       (t<t0: i->s)
           !
           if ( qci(i,k,2) > d_zero ) then
             qimax = roqimax/den(i,k)
             psaut(i,k) = max(d_zero,(qci(i,k,2)-qimax)*rdtcld)
           end if
           !
           ! pgaut: conversion(aggregation) of snow to graupel [HL A4] [LFO 37]
           !        (T<T0: QS->QG)
           !
            if ( qrs(i,k,2) > d_zero ) then
              alpha2 = 1.e-3_rkx*exp(0.09_rkx*(-supcol))
              pgaut(i,k) = min(max(d_zero, &
                alpha2*(qrs(i,k,2)-qs0)),qrs(i,k,2)/dtcld)
            end if
          end if

          !
          ! phaut: conversion(aggregation) of grauple to hail [BHT A18]
          !        (T<T0: QG->QH)
          !
          if ( qrs(i,k,3) > d_zero ) then
            alpha2 = 1.e-3_rkx*exp(0.09_rkx*(-supcol))
            phaut(i,k) = min(max(d_zero, &
              alpha2*(qrs(i,k,3)-qs0)),qrs(i,k,3)/dtcld)
          end if
         !
         ! psevp: evaporation of melting snow [hl a35] [rh83 a27]
         !       (t>t0: s->v)
         !
         if ( supcol <= d_zero ) then
           if ( qrs(i,k,2) > d_zero .and. rh(i,k,1) < d_one ) then
             coeres = rslope2(i,k,2)*sqrt(rslope(i,k,2)*rslopeb(i,k,2))
             psevp(i,k) = (rh(i,k,1)-d_one)*n0sfac(i,k) * &
               (precs1*rslope2(i,k,2)+precs2*work2(i,k)*coeres)/work1(i,k,1)
             psevp(i,k) = min(max(psevp(i,k),-qrs(i,k,2)*rdtcld),d_zero)
           end if
            !
            ! pgevp: Evaporation of melting graupel [HL A25] [RH84 A19]
            !       (T>=T0: QG->QV)
            !
            if ( qrs(i,k,3) > d_zero .and. rh(i,k,1) < d_one ) then
              coeres = rslope2(i,k,3)*sqrt(rslope(i,k,3)*rslopeb(i,k,3))
              pgevp(i,k) = (rh(i,k,1)-d_one) * &
                (precg1*rslope2(i,k,3)+precg2*work2(i,k)*coeres)/work1(i,k,1)
              pgevp(i,k) = min(max(pgevp(i,k),-qrs(i,k,3)/dtcld),d_zero)
            end if
            !
            ! phevp: Evaporation of melting hail [BHT A20]
            !       (T>=T0: QH->QV)
            !
            if ( qrs(i,k,4) > d_zero .and. rh(i,k,1) < d_one ) then
              coeres = rslope2(i,k,4)*sqrt(rslope(i,k,4)*rslopeb(i,k,4))
              phevp(i,k) = (rh(i,k,1)-d_one)*(prech1*rslope2(i,k,4) + &
                           prech2*work2(i,k)*coeres)/work1(i,k,1)
              phevp(i,k) = min(max(phevp(i,k),-qrs(i,k,4)/dtcld),d_zero)
            end if
          end if
        end do
      end do
      !
      ! check mass conservation of generation terms and feedback to the
      ! large scale
      !
      do k = 1 , kz
        do i = ims , ime
          delta2 = d_zero
          delta3 = d_zero
          if ( qrs(i,k,1) < 1.e-4_rkx .and. &
               qrs(i,k,2) < 1.e-4_rkx) delta2 = d_one
          if ( qrs(i,k,1) < 1.e-4_rkx) delta3 = d_one
          if ( t(i,k) <= tzero ) then
            !
            ! cloud water
            !
            qval = max(qci(i,k,1),qcimin)
            source = (praut(i,k)+pracw(i,k)+paacw(i,k) + &
                      paacw(i,k)+phacw(i,k))*dtcld
            if ( source > qval ) then
              factor = qval/source
              praut(i,k) = praut(i,k)*factor
              pracw(i,k) = pracw(i,k)*factor
              paacw(i,k) = paacw(i,k)*factor
              phacw(i,k) = phacw(i,k)*factor
            end if
            !
            ! cloud ice
            !
            qval = max(qci(i,k,2),qcimin)
            source = (psaut(i,k)-pigen(i,k)-pidep(i,k)+praci(i,k) + &
                      psaci(i,k)+pgaci(i,k)+phaci(i,k))*dtcld
            if ( source > qval ) then
              factor = qval/source
              psaut(i,k) = psaut(i,k)*factor
              pigen(i,k) = pigen(i,k)*factor
              pidep(i,k) = pidep(i,k)*factor
              praci(i,k) = praci(i,k)*factor
              psaci(i,k) = psaci(i,k)*factor
              pgaci(i,k) = pgaci(i,k)*factor
              phaci(i,k) = phaci(i,k)*factor
            end if
            !
            ! rain
            !
            qval = max(qrs(i,k,1),qrsmin)
            source = (-praut(i,k)-prevp(i,k)-pracw(i,k) + &
                      piacr(i,k)+psacr(i,k)+pgacr(i,k)+phacr(i,k))*dtcld
            if (source > qval) then
              factor = qval/source
              praut(i,k) = praut(i,k)*factor
              prevp(i,k) = prevp(i,k)*factor
              pracw(i,k) = pracw(i,k)*factor
              piacr(i,k) = piacr(i,k)*factor
              psacr(i,k) = psacr(i,k)*factor
              pgacr(i,k) = pgacr(i,k)*factor
              phacr(i,k) = phacr(i,k)*factor
            end if
            !
            ! snow
            !
            qval = max(qrs(i,k,2),qrsmin)
            source = -(psdep(i,k)+psaut(i,k)+paacw(i,k)+pvapg(i,k) + &
                       pvaph(i,k)+psaci(i,k)-pgaut(i,k)-pracs(i,k) * &
                       (d_one-delta2)+piacr(i,k)*delta3+praci(i,k) * &
                       delta3+psacr(i,k)*delta2-pgacs(i,k)-phacs(i,k))*dtcld
            if ( source > qval ) then
              factor = qval/source
              psdep(i,k) = psdep(i,k)*factor
              psaut(i,k) = psaut(i,k)*factor
              pgaut(i,k) = pgaut(i,k)*factor
              paacw(i,k) = paacw(i,k)*factor
              pvapg(i,k) = pvapg(i,k)*factor
              pvaph(i,k) = pvaph(i,k)*factor
              psaci(i,k) = psaci(i,k)*factor
              piacr(i,k) = piacr(i,k)*factor
              praci(i,k) = praci(i,k)*factor
              psacr(i,k) = psacr(i,k)*factor
              pracs(i,k) = pracs(i,k)*factor
              pgacs(i,k) = pgacs(i,k)*factor
              phacs(i,k) = phacs(i,k)*factor
            end if
            !
            ! graupel
            !
            qval = max(qrs(i,k,3),qrsmin)
            source = -(pgdep(i,k)+pgaut(i,k)+pgaci(i,k)+paacw(i,k) + &
                       pgacs(i,k)+piacr(i,k)*(d_one-delta3)+praci(i,k) * &
                       (d_one-delta3)+psacr(i,k)*(d_one-delta2) + &
                       pgacr(i,k)*delta2+pracs(i,k)*(d_one-delta2) - &
                       pracg(i,k)*(d_one-delta2)-phaut(i,k)-pvapg(i,k) - &
                       phacg(i,k)+primh(i,k))*dtcld
            if ( source > qval ) then
              factor = qval/source
              pgdep(i,k) = pgdep(i,k)*factor
              pgaut(i,k) = pgaut(i,k)*factor
              phaut(i,k) = phaut(i,k)*factor
              piacr(i,k) = piacr(i,k)*factor
              praci(i,k) = praci(i,k)*factor
              pracs(i,k) = pracs(i,k)*factor
              pracg(i,k) = pracg(i,k)*factor
              psacr(i,k) = psacr(i,k)*factor
              paacw(i,k) = paacw(i,k)*factor
              pgaci(i,k) = pgaci(i,k)*factor
              pgacr(i,k) = pgacr(i,k)*factor
              pgacs(i,k) = pgacs(i,k)*factor
              pvapg(i,k) = pvapg(i,k)*factor
              phacg(i,k) = phacg(i,k)*factor
              primh(i,k) = primh(i,k)*factor
            end if
!            !
!            ! hail
!            !
            qval = max(qrs(i,k,4),qrsmin)
            source = -(phdep(i,k)+phaut(i,k)                               + &
                       pgacr(i,k)*(d_one-delta2)+pracg(i,k)*(d_one-delta2) + &
                       phacw(i,k)+phacr(i,k)+phaci(i,k)+phacs(i,k)         + &
                       phacg(i,k)-pvaph(i,k)-primh(i,k))*dtcld
            if (source > qval) then
              factor = qval/source
              phdep(i,k) = phdep(i,k)*factor
              phaut(i,k) = phaut(i,k)*factor
              pracg(i,k) = pracg(i,k)*factor
              pgacr(i,k) = pgacr(i,k)*factor
              phacw(i,k) = phacw(i,k)*factor
              phaci(i,k) = phaci(i,k)*factor
              phacr(i,k) = phacr(i,k)*factor
              phacs(i,k) = phacs(i,k)*factor
              phacg(i,k) = phacg(i,k)*factor
              pvaph(i,k) = pvaph(i,k)*factor
              primh(i,k) = primh(i,k)*factor
            end if

            work2(i,k)=-(prevp(i,k)+psdep(i,k)+pgdep(i,k)+phdep(i,k) + &
                         pigen(i,k)+pidep(i,k))
           ! update
            qv(i,k) = qv(i,k)+work2(i,k)*dtcld
            qci(i,k,1) = max(qci(i,k,1)-(praut(i,k)+pracw(i,k) + &
                         paacw(i,k)+paacw(i,k)+phacw(i,k))*dtcld,d_zero)
            qrs(i,k,1) = max(qrs(i,k,1)+(praut(i,k)+pracw(i,k) + &
                            prevp(i,k)-piacr(i,k)-pgacr(i,k)   - &
                            psacr(i,k)-phacr(i,k))*dtcld,d_zero)
            qci(i,k,2) = max(qci(i,k,2)-(psaut(i,k)+praci(i,k)+psaci(i,k) + &
                            pgaci(i,k)+phaci(i,k)-pigen(i,k)-pidep(i,k))  * &
                            dtcld,d_zero)
            qrs(i,k,2) = max(qrs(i,k,2)+(psdep(i,k)+psaut(i,k)+paacw(i,k) + &
                            pvapg(i,k)+pvaph(i,k)-pgaut(i,k)              + &
                            psaci(i,k)-pgacs(i,k)-phacs(i,k)              + &
                            piacr(i,k)*delta3+praci(i,k)*delta3           + &
                            psacr(i,k)*delta2                             - &
                            pracs(i,k)*(d_one-delta2))                    * &
                            dtcld,d_zero)
            qrs(i,k,3) = max(qrs(i,k,3)+(pgdep(i,k)+pgaut(i,k)   + &
                            piacr(i,k)*(d_one-delta3)            + &
                            praci(i,k)*(d_one-delta3)            + &
                            psacr(i,k)*(d_one-delta2)            + &
                            pgacr(i,k)*delta2                    + &
                            pgaci(i,k)+paacw(i,k)                + &
                            pgacs(i,k)+primh(i,k)                + &
                            pracs(i,k)*(d_one-delta2)            - &
                            pracg(i,k)*(d_one-delta2)            - &
                            phaut(i,k)-pvapg(i,k)-phacg(i,k))    * &
                            dtcld,d_zero)
            qrs(i,k,4) = max(qrs(i,k,4)+(phdep(i,k)+phaut(i,k)          + &
                            pgacr(i,k)*(d_one-delta2)                   + &
                            pracg(i,k)*(d_one-delta2)                   + &
                            phacw(i,k)+phacr(i,k)+phaci(i,k)+phacs(i,k) + &
                            phacg(i,k)-pvaph(i,k)-primh(i,k))           * &
                            dtcld,d_zero)
           qrs(i,k,1) = max(qrs(i,k,1) + &
              (praut(i,k)+pracw(i,k)+prevp(i,k))*dtcld, d_zero)
            qci(i,k,2) = max(qci(i,k,2) - &
              (psaut(i,k)+psaci(i,k)-pigen(i,k)-pidep(i,k))*dtcld, d_zero)
            qrs(i,k,2) = max(qrs(i,k,2) + &
              (psdep(i,k)+psaut(i,k)+psaci(i,k)+psacw(i,k))*dtcld, d_zero)
            xlf = max(wlhs-xl(i,k),d_zero)
            xlwork2 = -wlhs*(psdep(i,k)+pgdep(i,k)+phdep(i,k)+pidep(i,k) + &
                       pigen(i,k))-xl(i,k)*prevp(i,k)                    - &
                       xlf*(piacr(i,k)+paacw(i,k)+paacw(i,k)+phacw(i,k)  + &
                       phacr(i,k)+pgacr(i,k)+psacr(i,k))
            t(i,k) = t(i,k)-xlwork2/cpm(i,k)*dtcld
          else   ! T > tzero
            !
            ! cloud water
            !
            qval = max(qci(i,k,1),qcimin)
            source = (praut(i,k)+pracw(i,k)+paacw(i,k) + &
                      paacw(i,k)+phacw(i,k))*dtcld
            if ( source > qval ) then
              factor = qval/source
              praut(i,k) = praut(i,k)*factor
              pracw(i,k) = pracw(i,k)*factor
              paacw(i,k) = paacw(i,k)*factor
              phacw(i,k) = phacw(i,k)*factor
            end if
            !
            ! rain
            !
            qval = max(qrs(i,k,1),qrsmin)
            source = (pseml(i,k)+pgeml(i,k)+pheml(i,k)             - &
                      pracw(i,k)-paacw(i,k)-paacw(i,k)-phacw(i,k)  - &
                      prevp(i,k)-praut(i,k))*dtcld
            if ( source > qval ) then
              factor = qval/source
              praut(i,k) = praut(i,k)*factor
              prevp(i,k) = prevp(i,k)*factor
              pracw(i,k) = pracw(i,k)*factor
              paacw(i,k) = paacw(i,k)*factor
              phacw(i,k) = phacw(i,k)*factor
              pseml(i,k) = pseml(i,k)*factor
              pgeml(i,k) = pgeml(i,k)*factor
              pheml(i,k) = pheml(i,k)*factor
            end if
            !
            ! snow
            !
            qval = max(qrs(i,k,2),qrsmin)
            source=(pgacs(i,k)+phacs(i,k)-pseml(i,k)-psevp(i,k))*dtcld
            if ( source > qval ) then
              factor = qval/source
              pgacs(i,k) = pgacs(i,k)*factor
              phacs(i,k) = phacs(i,k)*factor
              psevp(i,k) = psevp(i,k)*factor
              pseml(i,k) = pseml(i,k)*factor
            end if
            !
            ! graupel
            !
            qval = max(qrsmin,qrs(i,k,3))
            source=-(pgacs(i,k)+pgevp(i,k)+pgeml(i,k)-phacg(i,k))*dtcld
            if (source > qval) then
              factor = qval/source
              pgacs(i,k) = pgacs(i,k)*factor
              pgevp(i,k) = pgevp(i,k)*factor
              pgeml(i,k) = pgeml(i,k)*factor
              phacg(i,k) = phacg(i,k)*factor
            end if
!            !
!            ! hail
!            !
            qval = max(qrsmin,qrs(i,k,4))
            source=-(phacs(i,k)+phacg(i,k)+phevp(i,k)+pheml(i,k))*dtcld
            if (source > qval) then
              factor = qval/source
              phacs(i,k) = phacs(i,k)*factor
              phacg(i,k) = phacg(i,k)*factor
              phevp(i,k) = phevp(i,k)*factor
              pheml(i,k) = pheml(i,k)*factor
            end if
            work2(i,k)=-(prevp(i,k)+psevp(i,k)+pgevp(i,k)+phevp(i,k))
            ! update
            qv(i,k) = qv(i,k)+work2(i,k)*dtcld
            qci(i,k,1) = max(qci(i,k,1)-(praut(i,k)+pracw(i,k) + &
                     paacw(i,k)+paacw(i,k)+phacw(i,k))*dtcld,d_zero)
            qrs(i,k,1) = max(qrs(i,k,1)+(praut(i,k)+pracw(i,k)      + &
                     prevp(i,k)+paacw(i,k)+paacw(i,k)+phacw(i,k)    - &
                     pseml(i,k)-pgeml(i,k)-pheml(i,k))*dtcld,d_zero)
            qrs(i,k,2) = max(qrs(i,k,2)+(psevp(i,k)+pseml(i,k)      - &
                     pgacs(i,k)-phacs(i,k))*dtcld,d_zero)
            qrs(i,k,3) = max(qrs(i,k,3)+(pgacs(i,k)+pgevp(i,k)+pgeml(i,k) - &
                     phacg(i,k))*dtcld,d_zero)
            qrs(i,k,4) = max(qrs(i,k,4)+(phacs(i,k)+phacg(i,k)+phevp(i,k) + &
                     pheml(i,k))*dtcld,d_zero)
            xlf = max(wlhs-xl(i,k),d_zero)
            xlwork2 = -xl(i,k)*(prevp(i,k)+psevp(i,k)+pgevp(i,k) + &
                       phevp(i,k))-xlf*(pseml(i,k)+pgeml(i,k)+pheml(i,k))
            t(i,k) = t(i,k)-xlwork2/cpm(i,k)*dtcld
          end if
        end do
      end do
      do k = 1 , kz
        do i = ims , ime
          tr = wattp/t(i,k)
          qs(i,k,1) = psat*exp(log(tr)*(xa))*exp(xb*(1.0_rkx-tr))
          qs(i,k,1) = min(qs(i,k,1),0.99_rkx*p(i,k))
          qs(i,k,1) = ep2 * qs(i,k,1) / (p(i,k) - qs(i,k,1))
          qs(i,k,1) = max(qs(i,k,1),minqq)
          rh(i,k,1) = max(qv(i,k)/qs(i,k,1),minqq)
          if ( t(i,k) < wattp ) then
            qs(i,k,2) = psat*exp(log(tr)*(xai))*exp(xbi*(1.0_rkx-tr))
          else
            qs(i,k,2) = psat*exp(log(tr)*(xa))*exp(xb*(1.0_rkx-tr))
          endif
          qs(i,k,2) = min(qs(i,k,2),0.99_rkx*p(i,k))
          qs(i,k,2) = ep2 * qs(i,k,2) / (p(i,k) - qs(i,k,2))
          qs(i,k,2) = max(qs(i,k,2),minqq)
          rh(i,k,2) = max(qv(i,k)/qs(i,k,2),minqq)
        end do
      end do
      ! pcond: condensational/evaporational rate of cloud water
      !        [hl a46] [rh83 a6]
      ! if there exists additional water vapor condensated/if
      ! evaporation of cloud water is not enough to remove subsaturation
      !
      do k = 1 , kz
        do i = ims , ime
          work1(i,k,1) = ((max(qv(i,k),minqq)-qs(i,k,1))) /  &
               (d_one+(xl(i,k))*(xl(i,k))/(rwat*(cpm(i,k)))*(qs(i,k,1)) / &
               ((t(i,k))*(t(i,k))))
          work2(i,k) = qci(i,k,1)+work1(i,k,1)
          pcond(i,k) = min(max(work1(i,k,1)*rdtcld,d_zero), &
                           max(qv(i,k),minqq)*rdtcld)
          if ( qci(i,k,1) > d_zero .and. work1(i,k,1) < d_zero ) then
            pcond(i,k) = max(work1(i,k,1),-qci(i,k,1))*rdtcld
          end if
          qv(i,k) = qv(i,k) - pcond(i,k)*dtcld
          qci(i,k,1) = max(qci(i,k,1)+pcond(i,k)*dtcld,d_zero)
          t(i,k) = t(i,k) + pcond(i,k)*xl(i,k)/cpm(i,k)*dtcld
        end do
      end do
    end do bigloop

    do k = 1 , kz
      do i = ims , ime
        if ( qrs(i,k,1) < d_zero ) qrs(i,k,1) = d_zero
        if ( qrs(i,k,2) < d_zero ) qrs(i,k,2) = d_zero
          if ( qci(i,k,1) < minqq ) qci(i,k,1) = d_zero
          if ( qci(i,k,2) < minqq ) qci(i,k,2) = d_zero
      end do
    end do

  contains

#include <wlh.inc>

    pure real(rkx) function cpmcal(q)
      implicit none
      real(rkx) , intent(in) :: q
      cpmcal = cpd*(d_one-max(q,minqq)) + cpv*max(q,minqq)
    end function cpmcal

    ! diffus: diffusion coefficient of the water vapor
    pure real(rkx) function diffus(x,y)
      implicit none
      real(rkx) , intent(in) :: x , y
      diffus = 8.794e-5_rkx * exp(log(x)*(1.81_rkx)) / y
    end function diffus

    ! viscos: kinematic viscosity(m2s-1)
    pure real(rkx) function viscos(x,y)
      implicit none
      real(rkx) , intent(in) :: x , y
      viscos = 1.496e-6_rkx * (x*sqrt(x)) /(x+120.0_rkx)/y
      ! viscos = 1.496e-6_rkx *x**1.5_rkx / (x+120.0_rkx)/y
    end function viscos

    pure real(rkx) function xka(x,y)
      implicit none
      real(rkx) , intent(in) :: x , y
      xka = 1.414e3_rkx * viscos(x,y) * y
    end function xka

    pure real(rkx) function diffac(a,b,c,d,e)
      implicit none
      real(rkx) , intent(in) :: a , b , c , d , e
      diffac = d*a*a/(xka(c,d)*rwat*c*c)+d_one/(e*diffus(c,b))
    end function diffac

    pure real(rkx) function venfac(a,b,c)
      implicit none
      real(rkx) , intent(in) :: a , b , c
      venfac = exp(log((viscos(b,c)/diffus(b,a)))*((onet))) / &
                   sqrt(viscos(b,c))*sqrt(sqrt(stdrho/c))
    end function venfac

    pure real(rkx) function conden(a,b,c,d,e)
      implicit none
      real(rkx) , intent(in) :: a , b , c , d , e
      conden = (max(b,minqq)-c)/(d_one+d*d/(rwat*e)*c/(a*a))
    end function conden

  end subroutine wsm72d

  subroutine slope_wsm7(qrs,den,denfac,t,rslope,rslopeb,rslope2,rslope3,vt, &
                        ims,ime)
    implicit none
    integer(ik4) , intent(in) :: ims , ime
    real(rkx) , dimension(ims:ime,kz,4) , intent(in) :: qrs
    real(rkx) , dimension(ims:ime,kz) , intent(in) :: den , denfac , t
    real(rkx) , dimension(ims:ime,kz,4) , intent(out) :: rslope , rslopeb
    real(rkx) , dimension(ims:ime,kz,4) , intent(out) :: rslope2 , rslope3 , vt
    real(rkx) :: supcol , n0sfac
    integer(ik4) :: i , k

    do k = 1 , kz
      do i = ims , ime
        supcol = tzero-t(i,k)
        n0sfac = max(min(exp(alpha*supcol),n0smax/n0s),d_one)
        !
        ! n0s: intercept parameter for snow [m-4] [hdc 6]
        !
        if ( qrs(i,k,1) <= qrsmin ) then
          rslope(i,k,1) = rslopermax
          rslopeb(i,k,1) = rsloperbmax
          rslope2(i,k,1) = rsloper2max
          rslope3(i,k,1) = rsloper3max
        else
          rslope(i,k,1) = d_one/lamdar(qrs(i,k,1),den(i,k))
          rslopeb(i,k,1) = exp(log(rslope(i,k,1))*(bvtr))
          rslope2(i,k,1) = rslope(i,k,1)*rslope(i,k,1)
          rslope3(i,k,1) = rslope2(i,k,1)*rslope(i,k,1)
        end if
        if ( qrs(i,k,2) <= qrsmin ) then
          rslope(i,k,2) = rslopesmax
          rslopeb(i,k,2) = rslopesbmax
          rslope2(i,k,2) = rslopes2max
          rslope3(i,k,2) = rslopes3max
        else
          rslope(i,k,2) = d_one/lamdas(qrs(i,k,2),den(i,k),n0sfac)
          rslopeb(i,k,2) = exp(log(rslope(i,k,2))*(bvts))
          rslope2(i,k,2) = rslope(i,k,2)*rslope(i,k,2)
          rslope3(i,k,2) = rslope2(i,k,2)*rslope(i,k,2)
        end if
        if ( qrs(i,k,3) <= qrsmin ) then
          rslope(i,k,3) = rslopegmax
          rslopeb(i,k,3) = rslopegbmax
          rslope2(i,k,3) = rslopeg2max
          rslope3(i,k,3) = rslopeg3max
        else
          rslope(i,k,3) = d_one/lamdag(qrs(i,k,3),den(i,k))
          rslopeb(i,k,3) = rslope(i,k,3)**bvtg
          rslope2(i,k,3) = rslope(i,k,3)*rslope(i,k,3)
          rslope3(i,k,3) = rslope2(i,k,3)*rslope(i,k,3)
        end if
        if ( qrs(i,k,4) <= qrsmin) then
          rslope(i,k,4) = rslopehmax
          rslopeb(i,k,4) = rslopehbmax
          rslope2(i,k,4) = rslopeh2max
          rslope3(i,k,4) = rslopeh3max
        else
          rslope(i,k,4) = d_one/lamdah(qrs(i,k,4),den(i,k))
          rslopeb(i,k,4) = rslope(i,k,4)**bvth
          rslope2(i,k,4) = rslope(i,k,4)*rslope(i,k,4)
          rslope3(i,k,4) = rslope2(i,k,4)*rslope(i,k,4)
        end if
        vt(i,k,1) = pvtr*rslopeb(i,k,1)*denfac(i,k)
        vt(i,k,2) = pvts*rslopeb(i,k,2)*denfac(i,k)
        vt(i,k,3) = pvtg*rslopeb(i,k,3)*denfac(i,k)
        vt(i,k,4) = pvth*rslopeb(i,k,4)*denfac(i,k)
        if ( qrs(i,k,1) <= d_zero ) vt(i,k,1) = d_zero
        if ( qrs(i,k,2) <= d_zero ) vt(i,k,2) = d_zero
        if ( qrs(i,k,2) <= d_zero ) vt(i,k,3) = d_zero
        if ( qrs(i,k,2) <= d_zero ) vt(i,k,4) = d_zero
      end do
    end do

    contains

    pure real(rkx) function lamdar(x,y)
      implicit none
      real(rkx) , intent(in) :: x , y
      lamdar = sqrt(sqrt(pidn0r/(x*y)))
    end function lamdar

    pure real(rkx) function lamdas(x,y,z)
      implicit none
      real(rkx) , intent(in) :: x , y , z
      lamdas = sqrt(sqrt(pidn0s*z/(x*y)))
    end function lamdas

    pure real(rkx) function lamdag(x,y)
      implicit none
      real(rkx) , intent(in) :: x , y
      lamdag = sqrt(sqrt(pidn0g/(x*y)))
    end function lamdag

   pure real(rkx) function lamdah(x,y)
      implicit none
      real(rkx) , intent(in) :: x , y
      lamdah = sqrt(sqrt(pidn0h/(x*y)))
   end function lamdah

  end subroutine slope_wsm7

  subroutine slope_rain(qrs,den,denfac,rslope,rslopeb,rslope2,rslope3,vt)
    implicit none
    real(rkx) , dimension(kz) , intent(in) :: qrs , den , denfac
    real(rkx) , dimension(kz) , intent(out) :: rslope , rslopeb
    real(rkx) , dimension(kz) , intent(out) :: rslope2 , rslope3 , vt
    integer(ik4) :: k

    do k = 1 , kz
      if ( qrs(k) <= qrsmin ) then
        rslope(k) = rslopermax
        rslopeb(k) = rsloperbmax
        rslope2(k) = rsloper2max
        rslope3(k) = rsloper3max
      else
        rslope(k) = d_one/lamdar(qrs(k),den(k))
        rslopeb(k) = rslope(k)**bvtr
        rslope2(k) = rslope(k)*rslope(k)
        rslope3(k) = rslope2(k)*rslope(k)
      end if
      vt(k) = pvtr*rslopeb(k)*denfac(k)
      if ( qrs(k) <= d_zero ) vt(k) = d_zero
    end do

    contains
    !
    ! size distributions: (x=mixing ratio, y=air density):
    ! valid for mixing ratio > 1.e-9 kg/kg.
    !
    pure real(rkx) function lamdar(x,y)
      implicit none
      real(rkx) , intent(in) :: x , y
      lamdar = sqrt(sqrt(pidn0r/(x*y)))
    end function lamdar

  end subroutine slope_rain

  subroutine slope_snow(qrs,den,denfac,t,rslope,rslopeb,rslope2,rslope3,vt)
    implicit none
    real(rkx) , dimension(kz) , intent(in) :: t , qrs , den , denfac
    real(rkx) , dimension(kz) , intent(out) :: rslope , rslopeb
    real(rkx) , dimension(kz) , intent(out) :: rslope2 , rslope3 , vt
    real(rkx) :: n0sfac , supcol
    integer(ik4) :: k

    do k = 1 , kz
      supcol = tzero-t(k)
      !
      ! n0s: intercept parameter for snow [m-4] [hdc 6]
      !
      n0sfac = max(min(exp(alpha*supcol),n0smax/n0s),d_one)
      if ( qrs(k) <= qrsmin ) then
        rslope(k) = rslopesmax
        rslopeb(k) = rslopesbmax
        rslope2(k) = rslopes2max
        rslope3(k) = rslopes3max
      else
        rslope(k) = d_one/lamdas(qrs(k),den(k),n0sfac)
        rslopeb(k) = rslope(k)**bvts
        rslope2(k) = rslope(k)*rslope(k)
        rslope3(k) = rslope2(k)*rslope(k)
      end if
      vt(k) = pvts*rslopeb(k)*denfac(k)
      if ( qrs(k) <= d_zero ) vt(k) = d_zero
    end do

    contains
    !
    ! size distributions: (x=mixing ratio, y=air density):
    ! valid for mixing ratio > 1.e-9 kg/kg.
    !
    pure real(rkx) function lamdas(x,y,z)
      implicit none
      real(rkx) , intent(in) :: x , y , z
      lamdas = sqrt(sqrt(pidn0s*z/(x*y)))
    end function lamdas

  end subroutine slope_snow

  subroutine slope_graup(qrs,den,denfac,rslope,rslopeb,rslope2,rslope3,vt)
    implicit none
    real(rkx) , dimension(kz) , intent(in) :: qrs , den , denfac
    real(rkx) , dimension(kz) , intent(out) :: rslope , rslopeb
    real(rkx) , dimension(kz) , intent(out) :: rslope2 , rslope3 , vt
    integer(ik4) :: k
!------------------------------------------------------------------------------
!     size distributions: (x=mixing ratio, y=air density):
!     valid for mixing ratio > 1.e-9 kg/kg.
    do k = 1 , kz
      if ( qrs(k) <= qrsmin ) then
        rslope(k) = rslopegmax
        rslopeb(k) = rslopegbmax
        rslope2(k) = rslopeg2max
        rslope3(k) = rslopeg3max
      else
        rslope(k) = d_one/lamdag(qrs(k),den(k))
        rslopeb(k) = rslope(k)**bvtg
        rslope2(k) = rslope(k)*rslope(k)
        rslope3(k) = rslope2(k)*rslope(k)
      end if
      vt(k) = pvtg*rslopeb(k)*denfac(k)
      if ( qrs(k) <= d_zero ) vt(k) = d_zero
    end do

    contains
    !
    pure real(rkx) function lamdag(x,y)
      implicit none
      real(rkx) , intent(in) :: x , y
      lamdag = sqrt(sqrt(pidn0g/(x*y)))
    end function lamdag

  end subroutine slope_graup

  subroutine slope_hail(qrs,den,denfac,rslope,rslopeb,rslope2,rslope3,vt)
    implicit none
    real(rkx) , dimension(kz) , intent(in) :: qrs , den , denfac
    real(rkx) , dimension(kz) , intent(out) :: rslope , rslopeb
    real(rkx) , dimension(kz) , intent(out) :: rslope2 , rslope3 , vt
    integer(ik4) :: k
!------------------------------------------------------------------------------
!     size distributions: (x=mixing ratio, y=air density):
!     valid for mixing ratio > 1.e-9 kg/kg.
    do k = 1 , kz
      if ( qrs(k) <= qrsmin ) then
        rslope(k) = rslopehmax
        rslopeb(k) = rslopehbmax
        rslope2(k) = rslopeh2max
        rslope3(k) = rslopeh3max
      else
        rslope(k) = d_one/lamdah(qrs(k),den(k))
        rslopeb(k) = rslope(k)**bvth
        rslope2(k) = rslope(k)*rslope(k)
        rslope3(k) = rslope2(k)*rslope(k)
      end if
      vt(k) = pvth*rslopeb(k)*denfac(k)
      if ( qrs(k) <= d_zero ) vt(k) = d_zero
    end do

    contains
    !
    pure real(rkx) function lamdah(x,y)
      implicit none
      real(rkx) , intent(in) :: x , y
      lamdah = sqrt(sqrt(pidn0h/(x*y)))
    end function lamdah

  end subroutine slope_hail

  !
  ! For non-iteration semi-lagrangain forward advection for cloud
  ! with mass conservation and positive definite advection
  ! 2nd order interpolation with monotonic piecewise linear method
  ! This routine is under assumption of decfl < 1 for semi_lagrangian
  !
  ! dzl    depth of model layer in meter
  ! wwl    terminal velocity at model layer m/s
  ! rql    cloud density*mixing ration
  ! precip precipitation
  ! dt     time step
  ! id     kind of precip: 0 test case; 1 raindrop  2: snow
  ! iter   how many time to guess mean terminal velocity: 0 pure forward.
  !        0 : use departure wind for advection
  !        1 : use mean wind for advection
  !        > 1 : use mean wind after iter-1 iterations
  !
  ! Author: Hann-Ming Henry Juang <henry.juang@noaa.gov>
  !         implemented by Song-You Hong
  !
  subroutine nislfv_rain_plm(im,denl,denfacl,tkl,dzl, &
                             wwl,rql,precip,dt,id,maxiter)
    implicit none
    integer(ik4) , intent(in) :: im , id , maxiter
    real(rkx) , dimension(im,kz) , intent(in) :: denl
    real(rkx) , dimension(im,kz) , intent(in) :: denfacl
    real(rkx) , dimension(im,kz) , intent(in) :: tkl
    real(rkx) , dimension(im,kz) , intent(in) :: dzl
    real(rkx) , dimension(im,kz) , intent(in) :: wwl
    real(rkx) , dimension(im,kz) , intent(inout) :: rql
    real(rkx) , dimension(im) , intent(out) :: precip
    real(rkx) , intent(in) :: dt

    integer(ik4) :: i , k , n , m , kk , kb , kt
    real(rkx) :: tl , tl2 , qql , dql , qqd
    real(rkx) :: th , th2 , qqh , dqh
    real(rkx) :: zsum , qsum , xdim , dip , d1 , d2 , con1
    real(rkx) :: allold , decfl
    real(rkx) , dimension(kz) :: dz , ww , qq , wd , wa , was
    real(rkx) , dimension(kz) :: den , denfac , tk
    real(rkx) , dimension(kzp1) :: wi , zi , za
    real(rkx) , dimension(kz) :: qn , qr ,tmp , tmp1 , tmp2 , tmp3
    real(rkx) , dimension(kzp1) :: dza , qa , qmi , qpi
    real(rkx) , parameter :: fa1 = 9.0_rkx/16.0_rkx
    real(rkx) , parameter :: fa2 = 1.0_rkx/16.0_rkx

    precip(:) = d_zero

    i_loop : &
    do i = 1 , im
      dz(:) = dzl(i,:)
      qq(:) = rql(i,:)
      ww(:) = wwl(i,:)
      den(:) = denl(i,:)
      denfac(:) = denfacl(i,:)
      tk(:) = tkl(i,:)
      ! skip for no precipitation for all layers
      allold = d_zero
      do k = 1 , kz
        allold = allold + qq(k)
      end do
      if ( allold <= d_zero ) then
        cycle i_loop
      end if
      !
      ! compute interface values
      !
      zi(1) = d_zero
      do k = 1 , kz
        zi(k+1) = zi(k)+dz(k)
      end do
      !
      ! save departure wind
      !
      wd(:) = ww(:)
      n = 1
      do
        ! plm is 2nd order, we can use 2nd order wi or 3rd order wi
        ! 2nd order interpolation to get wi
        wi(1) = ww(1)
        wi(kzp1) = ww(kz)
        do k = 2 , kz
          wi(k) = (ww(k)*dz(k-1)+ww(k-1)*dz(k))/(dz(k-1)+dz(k))
        end do
        ! 3rd order interpolation to get wi
        wi(1) = ww(1)
        wi(2) = d_half*(ww(2)+ww(1))
        do k = 3 , kzm1
          wi(k) = fa1*(ww(k)+ww(k-1))-fa2*(ww(k+1)+ww(k-2))
        end do
        wi(kz) = d_half*(ww(kz)+ww(kz-1))
        wi(kzp1) = ww(kz)
        !
        ! terminate of top of raingroup
        !
        do k = 2 , kz
          if ( abs(ww(k)) < epsilon(d_one) ) wi(k) = ww(k-1)
        end do
        !
        ! diffusivity of wi
        !
        con1 = 0.05_rkx
        do k = kz , 1 , -1
          decfl = (wi(k+1)-wi(k))*dt/dz(k)
          if ( decfl > con1 ) then
            wi(k) = wi(k+1) - con1*dz(k)/dt
          end if
        end do
        ! compute arrival point
        do k = 1 , kzp1
          za(k) = zi(k) - wi(k)*dt
        end do
        do k = 1 , kz
          dza(k) = za(k+1)-za(k)
        end do
        dza(kzp1) = zi(kzp1) - za(kzp1)
        !
        ! compute deformation at arrival point
        !
        do k = 1 , kz
          qa(k) = qq(k)*dz(k)/dza(k)
          qr(k) = qa(k)/den(k)
        end do
        qa(kzp1) = d_zero
        if ( n <= maxiter ) then
          !
          ! compute arrival terminal velocity, and estimate mean
          ! terminal velocity then back to use mean terminal velocity
          !
          if ( id == 1 ) then
            call slope_rain(qr,den,denfac,tmp,tmp1,tmp2,tmp3,wa)
          else
            call slope_snow(qr,den,denfac,tk,tmp,tmp1,tmp2,tmp3,wa)
          end if
          if ( n >= 2 ) wa(1:kz) = d_half*(wa(1:kz)+was(1:kz))
          do k = 1 , kz
            ! mean wind is average of departure and new arrival winds
            ww(k) = d_half * ( wd(k)+wa(k) )
          end do
          was(:) = wa(:)
          n = n + 1
        else
          exit
        end if
      end do
      !
      ! estimate values at arrival cell interface with monotone
      !
      do k = 2 , kz
        d1 = qa(k+1)-qa(k)
        d2 = qa(k)-qa(k-1)
        if ( d1 < 1.0e-20_rkx ) d1 = d_zero
        if ( d2 < 1.0e-20_rkx ) d2 = d_zero
        dip = d1 / (dza(k+1)+dza(k))
        xdim = d2 / (dza(k-1)+dza(k))
        if ( dip*xdim <= d_zero ) then
          qmi(k) = qa(k)
          qpi(k) = qa(k)
        else
          qpi(k) = qa(k) + d_half*(dip+xdim)*dza(k)
          qmi(k) = d_two*qa(k) - qpi(k)
          if( qpi(k) < d_zero .or. qmi(k) < d_zero ) then
            qpi(k) = qa(k)
            qmi(k) = qa(k)
          end if
        end if
      end do
      qpi(1) = qa(1)
      qmi(1) = qa(1)
      qmi(kzp1) = qa(kzp1)
      qpi(kzp1) = qa(kzp1)
      !
      ! interpolation to regular point
      !
      qn = d_zero
      kb = 1
      kt = 1
      intp : &
      do k = 1 , kz
        kb = max(kb-1,1)
        kt = max(kt-1,1)
        ! find kb and kt
        if ( zi(k) >= za(kzp1) ) then
          exit intp
        else
          find_kb : &
          do kk = kb , kz
            if ( zi(k) <= za(kk+1) ) then
              kb = kk
              exit find_kb
            else
              cycle find_kb
            end if
          end do find_kb
          find_kt : &
          do kk = kt , kz
            if ( zi(k+1) <= za(kk) ) then
              kt = kk
              exit find_kt
            else
              cycle find_kt
            end if
          end do find_kt
          kt = kt - 1
          ! compute q with piecewise constant method
          if ( kt == kb ) then
            tl = (zi(k)-za(kb)) / dza(kb)
            th = (zi(k+1)-za(kb)) / dza(kb)
            tl2 = tl*tl
            th2 = th*th
            qqd = d_half*(qpi(kb)-qmi(kb))
            qqh = qqd*th2+qmi(kb)*th
            qql = qqd*tl2+qmi(kb)*tl
            qn(k) = (qqh-qql)/(th-tl)
          else if ( kt > kb ) then
            tl = (zi(k)-za(kb))/dza(kb)
            tl2 = tl*tl
            qqd = d_half*(qpi(kb)-qmi(kb))
            qql = qqd*tl2+qmi(kb)*tl
            dql = qa(kb)-qql
            zsum = (d_one-tl)*dza(kb)
            qsum = dql*dza(kb)
            if ( kt-kb > 1 ) then
              do m = kb+1 , kt-1
                zsum = zsum + dza(m)
                qsum = qsum + qa(m) * dza(m)
              end do
            end if
            th = (zi(k+1)-za(kt))/dza(kt)
            th2 = th*th
            qqd = d_half*(qpi(kt)-qmi(kt))
            dqh = qqd*th2+qmi(kt)*th
            zsum = zsum + th*dza(kt)
            qsum = qsum + dqh*dza(kt)
            qn(k) = qsum/zsum
          end if
          cycle intp
        end if
      end do intp
      !
      ! rain out
      !
      sum_precip: &
      do k = 1 , kz
        if ( za(k) < d_zero .and. za(k+1) < d_zero ) then
          precip(i) = precip(i) + qa(k)*dza(k)
          cycle sum_precip
        else if ( za(k) < d_zero .and. za(k+1) >= d_zero ) then
          precip(i) = precip(i) + qa(k)*(d_zero-za(k))
          exit sum_precip
        end if
        exit sum_precip
      end do sum_precip
      !
      ! replace the new values
      !
      rql(i,:) = qn(:)
    end do i_loop
  end subroutine nislfv_rain_plm

  subroutine nislfv_rain_plm6(im,denl,denfacl,tkl,dzl,wwl,rql, &
                              rql2,precip1,precip2,dt,id,maxiter)
    implicit none
    integer(ik4) , intent(in) :: im , id , maxiter
    real(rkx) , dimension(im,kz) , intent(in) :: denl
    real(rkx) , dimension(im,kz) , intent(in) :: denfacl
    real(rkx) , dimension(im,kz) , intent(in) :: tkl
    real(rkx) , dimension(im,kz) , intent(in) :: dzl
    real(rkx) , dimension(im,kz) , intent(in) :: wwl
    real(rkx) , dimension(im,kz) , intent(inout) :: rql , rql2
    real(rkx) , dimension(im) , intent(out) :: precip1 , precip2
    real(rkx) , intent(in) :: dt

    integer(ik4) :: i , k , n , m , kk , kb , kt , ist
    real(rkx) :: tl , tl2 , qql , dql , qqd
    real(rkx) :: th , th2 , qqh , dqh
    real(rkx) :: zsum , qsum , xdim , dip , d1 , d2 , con1
    real(rkx) :: allold , decfl
    real(rkx) , dimension(im) :: precip
    real(rkx) , dimension(kz) :: dz , ww , qq , wd , wa , wa2 , was
    real(rkx) , dimension(kz) :: den , denfac , tk
    real(rkx) , dimension(kzp1) :: wi , zi , za
    real(rkx) , dimension(kz) :: qn , qr , qr2 , qq2 , tmp , tmp1 , tmp2 , tmp3
    real(rkx) , dimension(kzp1) :: dza , qa , qa2 , qmi , qpi
    real(rkx) , parameter :: fa1 = 9.0_rkx/16.0_rkx
    real(rkx) , parameter :: fa2 = 1.0_rkx/16.0_rkx

    precip(:) = d_zero
    precip1(:) = d_zero
    precip2(:) = d_zero

    j_loop : &
    do i = 1 , im
      dz(:) = dzl(i,:)
      qq(:) = rql(i,:)
      qq2(:) = rql2(i,:)
      ww(:) = wwl(i,:)
      den(:) = denl(i,:)
      denfac(:) = denfacl(i,:)
      tk(:) = tkl(i,:)
      ! skip for no precipitation for all layers
      allold = d_zero
      do k=1 , kz
        allold = allold + qq(k) + qq2(k)
      enddo
      if ( allold <= d_zero ) then
        cycle j_loop
      end if
      ! compute interface values
      zi(1) = d_zero
      do k = 1 , kz
        zi(k+1) = zi(k)+dz(k)
      end do
      ! save departure wind
      wd(:) = ww(:)
      n=1
      do !inf
        ! plm is 2nd order, we can use 2nd order wi or 3rd order wi
        ! 2nd order interpolation to get wi
        wi(1) = ww(1)
        wi(kzp1) = ww(kz)
        do k = 2 , kz
          wi(k) = (ww(k)*dz(k-1)+ww(k-1)*dz(k))/(dz(k-1)+dz(k))
        end do
        ! 3rd order interpolation to get wi
        wi(1) = ww(1)
        wi(2) = d_half*(ww(2)+ww(1))
        do k = 3 , kzm1
          wi(k) = fa1*(ww(k)+ww(k-1))-fa2*(ww(k+1)+ww(k-2))
        end do
        wi(kz) = d_half*(ww(kz)+ww(kz-1))
        wi(kzp1) = ww(kz)
        !
        ! terminate of top of raingroup
        !
        do k = 2 , kz
          if ( abs(ww(k)) < epsilon(d_one) ) wi(k) = ww(k-1)
        end do
        !diffusivity of wi
        con1 = 0.05_rkx
        do k = kz , 1 , -1
          decfl = (wi(k+1)-wi(k))*dt/dz(k)
          if ( decfl > con1 ) then
            wi(k) = wi(k+1) - con1*dz(k)/dt
          end if
        end do
        ! compute arrival point
        do k = 1 , kzp1
          za(k) = zi(k) - wi(k)*dt
        end do
        do k = 1 , kz
          dza(k) = za(k+1)-za(k)
        end do
        dza(kzp1) = zi(kzp1) - za(kzp1)
        ! compute deformation at arrival point
        do k = 1 , kz
          qa(k) = qq(k)*dz(k)/dza(k)
          qa2(k) = qq2(k)*dz(k)/dza(k)
          qr(k) = qa(k)/den(k)
          qr2(k) = qa2(k)/den(k)
        end do
        qa(kzp1) = d_zero
        qa2(kzp1) = d_zero
        if ( n <= maxiter ) then
          !
          ! compute arrival terminal velocity, and estimate mean
          ! terminal velocity then back to use mean terminal velocity
          !
          call slope_snow(qr,den,denfac,tk,tmp,tmp1,tmp2,tmp3,wa)
          call slope_graup(qr2,den,denfac,tmp,tmp1,tmp2,tmp3,wa2)
          do k = 1 , kz
            tmp(k) = max((qr(k)+qr2(k)), 1.e-15_rkx)
            if ( tmp(k) > 1.1e-15_rkx ) then
              wa(k) = (wa(k)*qr(k) + wa2(k)*qr2(k))/tmp(k)
            else
              wa(k) = d_zero
            end if
          end do
          if ( n >= 2 ) wa(1:kz) = d_half*(wa(1:kz)+was(1:kz))
          do k = 1 , kz
           ! mean wind is average of departure and new arrival winds
            ww(k) = d_half * ( wd(k)+wa(k) )
          end do
          was(:) = wa(:)
          n = n + 1
        else
          exit
        end if
        ist_loop : &
        do ist = 1 , 2
          if ( ist == 2 ) then
            qa(:) = qa2(:)
          else
            exit
          end if
        !end do ist_loop
        precip(i) = d_zero
      ! estimate values at arrival cell interface with monotone
      do k = 2 , kz
        d1 = qa(k+1)-qa(k)
        d2 = qa(k)-qa(k-1)
        if ( d1 < 1.0e-20_rkx ) d1 = d_zero
        if ( d2 < 1.0e-20_rkx ) d2 = d_zero
        dip = d1 / (dza(k+1)+dza(k))
        xdim = d2 / (dza(k-1)+dza(k))
        if ( dip*xdim <= d_zero ) then
          qmi(k) = qa(k)
          qpi(k) = qa(k)
        else
          qpi(k) = qa(k) + d_half*(dip+xdim)*dza(k)
          qmi(k) = d_two*qa(k) - qpi(k)
          if( qpi(k) < d_zero .or. qmi(k) < d_zero ) then
            qpi(k) = qa(k)
            qmi(k) = qa(k)
          end if
        end if
      end do
      qpi(1) = qa(1)
      qmi(1) = qa(1)
      qmi(kzp1) = qa(kzp1)
      qpi(kzp1) = qa(kzp1)
          ! interpolation to regular point
      !
      qn = d_zero
      kb = 1
      kt = 1
     intp : &
      do k = 1 , kz
        kb = max(kb-1,1)
        kt = max(kt-1,1)
        ! find kb and kt
        if ( zi(k) >= za(kzp1) ) then
          exit intp
        else
          find_kb : &
          do kk = kb , kz
            if ( zi(k) <= za(kk+1) ) then
              kb = kk
              exit find_kb
            else
              cycle find_kb
            end if
          end do find_kb
          find_kt : &
          do kk = kt , kz
            if ( zi(k+1) <= za(kk) ) then
              kt = kk
              exit find_kt
            else
              cycle find_kt
            end if
          end do find_kt
          kt = kt - 1
          ! compute q with piecewise constant method
          if ( kt == kb ) then
            tl = (zi(k)-za(kb)) / dza(kb)
            th = (zi(k+1)-za(kb)) / dza(kb)
            tl2 = tl*tl
            th2 = th*th
            qqd = d_half*(qpi(kb)-qmi(kb))
            qqh = qqd*th2+qmi(kb)*th
            qql = qqd*tl2+qmi(kb)*tl
            qn(k) = (qqh-qql)/(th-tl)
          else if ( kt > kb ) then
            tl = (zi(k)-za(kb))/dza(kb)
            tl2 = tl*tl
            qqd = d_half*(qpi(kb)-qmi(kb))
            qql = qqd*tl2+qmi(kb)*tl
           dql = qa(kb)-qql
            zsum = (d_one-tl)*dza(kb)
            qsum = dql*dza(kb)
            if ( kt-kb > 1 ) then
              do m = kb+1 , kt-1
                zsum = zsum + dza(m)
                qsum = qsum + qa(m) * dza(m)
              end do
            end if
            th = (zi(k+1)-za(kt))/dza(kt)
            th2 = th*th
            qqd = d_half*(qpi(kt)-qmi(kt))
            dqh = qqd*th2+qmi(kt)*th
            zsum = zsum + th*dza(kt)
            qsum = qsum + dqh*dza(kt)
            qn(k) = qsum/zsum
          end if
          cycle intp
        end if
      end do intp
      ! rain out
      !
      sum_precip1: &
      do k = 1 , kz
        if ( za(k) < d_zero .and. za(k+1) < d_zero ) then
          precip(i) = precip(i) + qa(k)*dza(k)
          cycle sum_precip1
        else if ( za(k) < d_zero .and. za(k+1) >= d_zero ) then
          precip(i) = precip(i) + qa(k)*(d_zero-za(k))
          exit sum_precip1
        end if
      end do sum_precip1
      if(ist.eq.1) then
        rql(i,:) = qn(:)
        precip1(i) = precip(i)
      else
        rql2(i,:) = qn(:)
        precip2(i) = precip(i)
      end if
    end do ist_loop
   end do
  end do j_loop

  end subroutine nislfv_rain_plm6
  !
  !  Compute radiation effective radii of cloud water, ice, and snow for
  !  single-moment microphysics.
  !  These are entirely consistent with microphysics assumptions, not
  !  constant or otherwise ad hoc as is internal to most radiation
  !  schemes.
  !  Coded and implemented by Soo Ya Bae, KIAPS, January 2015.
  !
!  subroutine effectrad_wsm7(t,qc,qi,qs,rho,re_qc,re_qi,re_qs)
!    implicit none
!    real(rkx) , dimension(kz) , intent(in) :: t
!    real(rkx) , dimension(kz) , intent(in) :: qc
!    real(rkx) , dimension(kz) , intent(in) :: qi
!    real(rkx) , dimension(kz) , intent(in) :: qs
!    real(rkx) , dimension(kz) , intent(in) :: rho
!    real(rkx) , dimension(kz) , intent(inout) :: re_qc
!    real(rkx) , dimension(kz) , intent(inout) :: re_qi
!    real(rkx) , dimension(kz) , intent(inout) :: re_qs
!
!    integer(ik4) :: k
!    real(rkx) , dimension(kz) :: ni
!    real(rkx) , dimension(kz) :: rqc
!    real(rkx) , dimension(kz) :: rqi
!    real(rkx) , dimension(kz) :: rni
!    real(rkx) , dimension(kz) :: rqs
!    real(rkx) :: temp , lamdac , lamdas , supcol , n0sfac
!    ! diameter of ice in m
!    real(rkx) :: diai
!    logical :: has_qc , has_qi , has_qs
!
!    !..minimum microphys values
!    real(rkx) , parameter :: r1 = 1.e-12_rkx
!    real(rkx) , parameter :: r2 = 1.e-6_rkx
!    !..mass power law relations:  mass = am*d**bm
!    real(rkx) , parameter :: bm_r = 3.0_rkx
!    real(rkx) , parameter :: obmr = 1.0_rkx/bm_r
!    real(rkx) , parameter :: nc0  = 3.e8_rkx
!
!    has_qc = .false.
!    has_qi = .false.
!    has_qs = .false.
!
!    do k = 1 , kz
!      ! for cloud
!      rqc(k) = max(r1, qc(k)*rho(k))
!      if ( rqc(k) > r1 ) has_qc = .true.
!      ! for ice
!      rqi(k) = max(r1, qi(k)*rho(k))
!      temp = (rho(k)*qi(k))
!      temp = sqrt(sqrt(temp*temp*temp))
!      ni(k) = min(max(5.38e7_rkx*temp,minni),maxni)
!      rni(k)= max(r2, ni(k)*rho(k))
!      if ( rqi(k) > r1 .and. rni(k) > r2 ) has_qi = .true.
!      ! for snow
!      rqs(k) = max(r1, qs(k)*rho(k))
!      if (rqs(k) > r1) has_qs = .true.
!    end do
!
!    if ( has_qc ) then
!      do k = 1 , kz
!        if ( rqc(k) <= r1 ) cycle
!        lamdac   = (pidnc*nc0/rqc(k))**obmr
!        re_qc(k) =  max(2.51e-6_rkx,min(1.5_rkx*(d_one/lamdac),50.e-6_rkx))
!      end do
!    end if
!
!    if ( has_qi ) then
!      do k = 1 , kz
!        if ( rqi(k) <= r1 .or. rni(k) <= r2 ) cycle
!        diai = 11.9_rkx*sqrt(rqi(k)/ni(k))
!        re_qi(k) = max(10.01e-6_rkx,min(0.75_rkx*0.163_rkx*diai,125.e-6_rkx))
!      end do
!    end if
!
!    if ( has_qs ) then
!      do k = 1 , kz
!        if ( rqs(k) <= r1 ) cycle
!        supcol = tzero-t(k)
!        n0sfac = max(min(exp(alpha*supcol),n0smax/n0s),d_one)
!        lamdas = sqrt(sqrt(pidn0s*n0sfac/rqs(k)))
!        re_qs(k) = max(25.e-6_rkx,min(0.5_rkx*(d_one/lamdas), 999.e-6_rkx))
!      end do
!    end if
!  end subroutine effectrad_wsm7

end module mod_micro_wsm7

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
