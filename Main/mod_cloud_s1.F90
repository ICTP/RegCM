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
 
module mod_cloud_s1

  use mod_realkinds
  use mod_dynparam
  use mod_mpmessage
  use mod_memutil
  use mod_atm_interface , only : atmstate , slice , surfstate
  use mod_runparams , only : iqqv => iqv                        !vapor
  use mod_runparams , only : iqql => iqc                        !liquid   
  use mod_runparams , only : iqqr => iqr                        !rain
  use mod_runparams , only : iqqi => iqi                        !ice
  use mod_runparams , only : iqqs => iqs                        !snow 
  use mod_pbl_common

  private

  integer(ik4) , parameter :: nmax = 5                               !number of progn.variables 
  real(rk8) , parameter :: ptsphy = 60.0                         !time step
  real(rk8) , parameter :: rg = 9.81                             !gravitation's const. 
  real(rk8) :: w2                                                !used for the conservation of water 
  real(rk8) :: w1                                                !used for the conservation of water
  real(rk8) :: difftemp        
  real(rk8) :: slht                                              ! sublimation latent heat
  real(rk8) :: zfall
  real(rk8) , pointer , dimension(:,:,:) :: pres                 ! from atms
  real(rk8) , pointer , dimension(:,:,:) :: zt                   ! from atms
  real(rk8) , pointer , dimension(:,:,:) :: zeta                 ! from atms
  real(rk8) , pointer , dimension(:,:,:) :: dzeta                ! from atms
  real(rk8) , pointer , dimension(:,:,:,:) :: zqxx               ! from atms
  real(rk8) , pointer , dimension(:,:,:) :: rhob3d               ! from atms 
  real(rk8) , pointer , dimension(:,:,:) :: omega                ! from atms 
  real(rk8) , pointer , dimension(:,:,:) :: heatrt               ! radiation heat rate
  real(rk8) , pointer , dimension(:,:,:) :: ztten                ! tendency of temperature
  real(rk8) , pointer , dimension(:,:,:) :: zqdetr               ! conv. detrained water
  real(rk8) , pointer , dimension(:,:,:,:) :: zqxten             ! tendency of zqx
  real(rk8) , pointer , dimension(:,:,:,:) :: zpfplsx            !j,i,k+1,n ! generalized precipitation flux
  integer(ik4) , parameter :: nqx = 5             

  public :: allocate_mod_cloud_s1 , init_cloud_s1 , microphys

  
  integer(ik4) , pointer , dimension(:) :: kphase                     ! marker for water phase of each species
                                                                 ! 0=vapour, 1=liquid, 2=ice  
  integer(ik4) , pointer , dimension(:) :: imelt                      ! marks melting linkage for ice categories
                                                                 ! ice->liquid, snow->rain

  real(rk8) :: zphases, zmelt
 
  real(rk8) , pointer , dimension(:,:) :: zfrzmax
  real(rk8) , pointer , dimension(:,:) :: zdp                     ! dp
  real(rk8) , pointer , dimension(:,:) :: zgdp                    ! g/dp
  real(rk8) , pointer , dimension(:,:) :: zdtgdp                  ! dt * g/dp
  real(rk8) , pointer , dimension(:,:) :: zrdtgdp                 ! dp / (dt * g)
  real(rk8) , pointer , dimension(:,:) :: zicetot
  real(rk8) , pointer , dimension(:,:) :: zmeltmax
  
  integer(ik4) , pointer , dimension(:,:,:) :: jindex2
  real(rk8) , pointer , dimension(:,:,:) :: ztl
  real(rk8) , pointer , dimension(:,:,:) :: ztln
  real(rk8) , pointer , dimension(:,:,:) :: dqsatdt
  real(rk8) , pointer , dimension(:,:,:) :: satvp
  real(rk8) , pointer , dimension(:,:,:) :: zqxfg
  real(rk8) , pointer , dimension(:,:,:) :: zfallsink
  real(rk8) , pointer , dimension(:,:,:) :: zfallsrce
  real(rk8) , pointer , dimension(:,:,:) :: zfluxq
  real(rk8) , pointer , dimension(:,:,:) :: zratio
  real(rk8) , pointer , dimension(:,:,:) :: zsinksum
  real(rk8) , pointer , dimension(:,:,:) :: totalw1
  real(rk8) , pointer , dimension(:,:,:) :: totalw2
  real(rk8) , pointer , dimension(:,:,:) :: diff
  real(rk8) , pointer , dimension(:,:,:) :: zconvsrce
  real(rk8) , pointer , dimension(:,:,:) :: zconvsink
  real(rk8) , pointer , dimension(:,:,:) :: zfoeew
  real(rk8) , pointer , dimension(:,:,:) :: zqsliq                ! liquid water saturation
  real(rk8) , pointer , dimension(:,:,:) :: zqsice                ! ice water saturation
  real(rk8) , pointer , dimension(:,:,:) :: zfoeeliq

  real(rk8) , pointer , dimension(:,:,:,:) :: wvflux
  real(rk8) , pointer , dimension(:,:,:,:) :: vqx
  real(rk8) , pointer , dimension(:,:,:,:) :: zqx
  real(rk8) , public  , pointer, dimension(:,:,:,:) :: zqxn
  real(rk8) , pointer , dimension(:,:,:,:) :: zqlhs
  real(rk8) , pointer , dimension(:,:,:,:) :: zsolqa
  real(rk8) , pointer , dimension(:,:,:,:) :: zsolqb
  integer(ik4) , pointer , dimension(:,:,:,:) :: jindex1
  
  contains

  subroutine allocate_mod_cloud_s1
    implicit none
    call getmem1d(imelt,1,nqx,'clouds1:imelt')
    call getmem1d(kphase,1,nqx,'clouds1:kphase')

    call getmem2d(zicetot,jci1,jci2,ici1,ici2,'clouds1:zicetot')
    call getmem2d(zmeltmax,jci1,jci2,ici1,ici2,'clouds1:zmeltmax')
    call getmem2d(zdp,jci1,jci2,ici1,ici2,'clouds1:zdp')
    call getmem2d(zgdp,jci1,jci2,ici1,ici2,'clouds1:zgdp')
    call getmem2d(zdtgdp,jci1,jci2,ici1,ici2,'clouds1:zdtgdp')
    call getmem2d(zrdtgdp,jci1,jci2,ici1,ici2,'clouds1:zrdtgdp')
    call getmem2d(zfrzmax,jci1,jci2,ici1,ici2,'clouds1:zfrzmax')

    call getmem3d(zqxfg,jci1,jci2,ici1,ici2,1,nqx,'clouds1:zqxfg')
    call getmem3d(zconvsrce,jci1,jci2,ici1,ici2,1,kz,'clouds1:zconvsrce')
    call getmem3d(zconvsink,jci1,jci2,ici1,ici2,1,kz,'clouds1:zconvsink')
    call getmem3d(zfoeew,jci1,jci2,ici1,ici2,1,kz,'clouds1:zfoeew') 
    call getmem3d(zfoeeliq,jci1,jci2,ici1,ici2,1,kz,'clouds1:zfoeeliq')
    call getmem3d(zqsice,jci1,jci2,ici1,ici2,1,kz,'clouds1:zqsice')
    call getmem3d(zqsliq,jci1,jci2,ici1,ici2,1,kz,'clouds1:zqsliq')
    call getmem3d(jindex2,jci1,jci2,ici1,ici2,1,nqx,'clouds1:jindex2')
    call getmem3d(totalw1,jci1,jci2,ici1,ici2,1,kz,'clouds1:totalw1')
    call getmem3d(totalw2,jci1,jci2,ici1,ici2,1,kz,'clouds1:totalw2')
    call getmem3d(zfallsink,jci1,jci2,ici1,ici2,1,nqx,'clouds1:zfallsink')
    call getmem3d(zfallsrce,jci1,jci2,ici1,ici2,1,nqx,'clouds1:zfallsrce')
    call getmem3d(zfluxq,jci1,jci2,ici1,ici2,1,nqx,'clouds1:zfluxq')
    call getmem3d(zratio,jci1,jci2,ici1,ici2,1,nqx,'clouds1:zratio')
    call getmem3d(zsinksum,jci1,jci2,ici1,ici2,1,nqx,'clouds1:zsinksum')
    call getmem3d(ztl,jci1,jci2,ici1,ici2,1,kz,'clouds1:ztl')
    call getmem3d(ztln,jci1,jci2,ici1,ici2,1,kz,'clouds1:ztln')
    call getmem3d(dqsatdt,jci1,jci2,ici1,ici2,1,kz,'clouds1:dqsatdt')
    call getmem3d(satvp,jci1,jci2,ici1,ici2,1,kz,'clouds1:satvp')
    call getmem3d(diff,jci1,jci2,ici1,ici2,1,kz,'clouds1:diff')

    call getmem4d(vqx,jci1,jci2,ici1,ici2,1,kz,1,nqx,'clouds1:vqx')
    call getmem4d(zqx,jci1,jci2,ici1,ici2,1,kz,1,nqx,'clouds1:zqx')
    call getmem4d(zqxn,jce1,jce2,ice1,ice2,1,kz,1,nqx,'clouds1:zqxn')
    call getmem4d(zqlhs,jci1,jci2,ici1,ici2,1,nqx,1,nqx,'clouds1:zqlhs')
    call getmem4d(zsolqa,jci1,jci2,ici1,ici2,1,nqx,1,nqx,'clouds1:zsolqa')
    call getmem4d(zsolqb,jci1,jci2,ici1,ici2,1,nqx,1,nqx,'clouds1:zsolqb')
    call getmem4d(jindex1,jci1,jci2,ici1,ici2,1,nqx,1,nqx,'clouds1:jindex1')
    call getmem4d(wvflux,jci1,jci2,ici1,ici2,1,kz,1,nqx,'clouds1:wvflux')
    call getmem4d(zpfplsx,jce1,jce2,ice1,ice2,1,kz+1,1,nqx,'clouds1:zpfplsx')

  end subroutine allocate_mod_cloud_s1

  subroutine init_cloud_s1(atms,aten,qdetr,heatrt)
    implicit none
    type(slice) , intent(in) :: atms
    type(atmstate) , intent(in) :: aten
    real(rk8) , pointer , intent(in), dimension(:,:,:) :: qdetr, heatrt
    call assignpnt(atms%pb3d,pres)
    call assignpnt(atms%tb3d,zt)
    call assignpnt(atms%za,zeta)
    call assignpnt(atms%dzq,dzeta)
    call assignpnt(atms%qxb3d,zqxx)
    call assignpnt(atms%rhob3d,rhob3d)
    call assignpnt(heatrt,radheatrt)
    call assignpnt(aten%qx,zqxten)
    call assignpnt(aten%t,ztten)
    call assignpnt(qdetr,zqdetr)
  end subroutine init_cloud_s1

  subroutine microphys(omega,jstart,jend,istart,iend)
    implicit none
    integer(ik4) , intent(in) :: jstart , jend , istart , iend
    integer(ik4) :: i , j , k , n , m 
    integer(ik4) :: iqqi , iqql , iqqr , iqqs , iqqv , jn , jo , kautoconv
    logical :: lmicro
    real(rk8) , pointer , dimension(:,:,:) , intent(in) :: omega             
    real(rk8) :: zexplicit
 
    ! local real variables for autoconversion rate constants
    real(rk8) :: alpha1                                              ! coefficient autoconversion
    real(rk8) , parameter :: zauto_rate_khair = 0.355D0              ! microphysical terms
    real(rk8) , parameter :: zauto_expon_khair = 1.47D0
    real(rk8) , parameter :: zauto_rate_sundq = 0.5D-3
    real(rk8) , parameter :: zauto_rate_kesl = 1.D-3                 !giusto!
    real(rk8) , parameter :: zauto_rate_klepi = 0.5D-3
    real(rk8) , parameter :: zautocrit = 5.D-4                       !giusto!
    real(rk8) , parameter :: zt0 = 273.16
    real(rk8) , parameter :: zepsec = 1.D-10
    real(rk8) , parameter :: qi0 = 1.0D-3                            !g g^(-1)
    real(rk8) , parameter :: vlht = 2.26D6                           !J/kg evaporation latent heat  
    real(rk8) , parameter :: rcpd = 1005.                            !J/(kg·K) Cp_dry  
    real(rk8) , parameter :: retv = 0.60                             !rv/rd-1   rv = 461.5, rd = 287.05              
    real(rk8) , parameter :: rd=287.05      
    
    ! local real constants and variables for condensation
    real(rk8) :: zcond  
    real(rk8) :: zdtdp  
    real(rk8) :: rovcp

    ! local real constants and variables for freezing
    real(rk8) :: zfrz
    real(rk8) :: ralsdcp
    real(rk8) :: ralvdcp
    real(rk8) :: zrldcp

    ! local real constants and variables for melting
    real(rk8) :: zsubsat
    real(rk8) :: ztdiff
    real(rk8) :: zcons1

    ! local real constants for evaporation
    real(rk8) , parameter :: kevap = 0.100D-02                       ! Raindrop evap rate coef
    real(rk8) , parameter :: rlmin = 1.D-8
    
    ! Numerical fit to wet bulb temperature
    real(rk8) , parameter :: ztw1 = 1329.31
    real(rk8) , parameter :: ztw2 = 0.0074615
    real(rk8) , parameter :: ztw3 = 0.85D5
    real(rk8) , parameter :: ztw4 = 40.637
    real(rk8) , parameter :: ztw5 = 275.0
    real(rk8) , parameter :: rtaumel=0.01188 

    ! Define species tags, hard wire for now
    iqqv = 1    ! vapour
    iqql = 2    ! liquid cloud water
    iqqr = 3    ! rain water
    iqqi = 4    ! ice
    iqqs = 5    ! snow
 
    ! Define which autoconversion paramaterization to be chosen
    !KAUTOCONV=1 ! Klein & Pincus (2000)
    !kautoconv = 2 ! Khairoutdinov and Kogan (2000)
    KAUTOCONV=3 ! Kessler (1969)
    !KAUTOCONV=4 ! Sundqvist
 
    ! Define species phase, 0=vapour, 1=liquid, 2=ice
    kphase(iqqv) = 0
    kphase(iqql) = 1
    kphase(iqqr) = 1
    kphase(iqqi) = 2
    kphase(iqqs) = 2

    ! Set up melting/freezing index,
    ! if an ice category melts/freezes, where does it go?

    imelt(iqqv)=-99
    imelt(iqql)=iqqi
    imelt(iqqr)=iqqs
    imelt(iqqi)=iqql
    imelt(iqqs)=iqqr


    !set initial value of the liquid temperature, ............
    ztl(jci1:jci2,ici1:ici2,1:kz) = zt(jci1:jci2,ici1:ici2,1:kz)
 
    do n = 1 , nqx
      do k = 1 , kz
        do j = jstart , jend
          do i = istart , iend
            if ( kphase(n) == 1 ) then
              ztl(j,i,k) = ztl(j,i,k) - vlht*zqxx(j,i,k,n)/rcpd
            end if
            if ( kphase(n) == 2 ) then
              ztl(j,i,k) = ztl(j,i,k) - slht*zqxx(j,i,k,n)/rcpd
            end if
          end do
        end do
      end do
    end do
    
    !-----------------
    ! loop over levels
    do k = 1 , kz - 1
    !-----------------
    ! First guess microphysics
    !---------------------------------
      do n = 1 , nqx
        do j = jstart , jend
          do i = istart , iend
            zqxfg(j,i,n)=zqxx(j,i,k,n)
          end do
        end do
      end do

 
     ! derived variables needed                                                  
      do j = jstart, jend                                                  
        do i = istart, iend                                                
          zdp(j,i)=pres(j,i,k+1)-pres(j,i,k)                        !dp                 
          zgdp(j,i)=rg/zdp(j,i)                                     !g/dp                 
          zdtgdp(j,i)=ptsphy*zgdp(j,i)                              !(dt*g)/dp                  
          zrdtgdp(j,i)=zdp(j,i)*(1.0/(ptsphy*rg))
        end do                                                             
      end do     

    ! set the fall velocities
      do j = jstart , jend
        do i = istart , iend
          vqx(j,i,k,iqqv) = d_zero !*sqrt(ZQX(JL,JK,IQV))
          vqx(j,i,k,iqql) = 0.1D0  !*sqrt(ZQX(JL,JK,IQL))
          vqx(j,i,k,iqqr) = 4.0D0  !*sqrt(ZQX(JL,JK,IQR))
          vqx(j,i,k,iqqi) = 0.15D0 !*sqrt(ZQX(JL,JK,IQI))
          vqx(j,i,k,iqqs) = 1.0D0  !*sqrt(ZQX(JL,JK,IQS))
        end do
      end do


    ! calculate cloud total water at the beginning
      do j = jstart , jend
        do i = istart , iend
          totalw1(j,i,k)= zqxx(j,i,k,iqqv)+zqxx(j,i,k,iqql)+zqxx(j,i,k,iqqr) + &
                          zqxx(j,i,k,iqqi)+zqxx(j,i,k,iqqs)
        end do
      end do
      w1=sum(totalw1)

    ! reset matrix so missing pathways are set
      zsolqb(:,:,:,:)  = d_zero  !_JPRB
      zsolqa(:,:,:,:)  = d_zero  !_JPRB
      zfluxq(:,:,:)    = d_zero  !_JPRB
      zfallsrce(:,:,:) = d_zero  !_JPRB
      zfallsink(:,:,:) = d_zero  !_JPRB
      zpfplsx(:,:,:,:) = d_zero 
      zconvsrce(:,:,:) = d_zero
      zconvsink(:,:,:) = d_zero
    !-------------------------------------------------------
    ! SOURCE/SINK array for implicit and explicit terms
    !-------------------------------------------------------
    !
    ! a POSITIVE value entered into the arrays is a...
    !
    !             Source of this variable
    !             |
    !             |   Sink of this variable
    !             |   |
    !             V   V
    ! ZSOLQA/B:q(IQa,IQb)
    !
    ! Thus if ZSOLQA/B(IQL,IQV)=K where K>0 then this is
    ! a source of IQL and a sink of IQV
    !
    ! put 'magic' source terms such as PLUDE from
    ! detrainment into explicit source/sink array diagnognal ! ZSOLQA(IQL,IQL)=PLUDE
    !--------------------------------------------------------
    ! Calculate the saturation mixing ratio
      do j = jstart , jend
        do i = istart , iend
          !Teton's formula for the saturation mixing ratio:
          if ( zt(j,i,k) > d_zero ) then
            satvp(j,i,k) = satw(zt(j,i,k),pres(j,i,k))
            dqsatdt(j,i,k) = dqsatdtw(zt(j,i,k),satvp(j,i,k))
          else
            satvp(j,i,k) = satc(zt(j,i,k),pres(j,i,k))
            dqsatdt(j,i,k) = dqsatdtc(zt(j,i,k),satvp(j,i,k))
          end if
        end do
      end do
 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                         CONDENSATION                       ! 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      rovcp=rd/rcpd
      do j = jstart , jend
        ! these are explicit sources and sinks
        do i = istart , iend
          zdtdp = rovcp*zt(j,i,k)/pres(j,i,k)
          zcond = dqsatdt(j,i,k)*((omega(j,i,k)*zdtdp)+radheatrt(j,i,k))
          zsolqa(j,i,iqqv,iqql) = zsolqa(j,i,iqqv,iqql) - zcond
          zsolqa(j,i,iqql,iqqv) = zsolqa(j,i,iqql,iqqv) + zcond
        end do
      end do 
    ! define the microphysics
    ! the matrix will be sparse is this a problem ?
    ! (X,Y) means a sink of X and a source of Y
    ! for the implementation I will use flexible pointers
    ! such that it will be written (IQR,IQG) to indicate graupel to rain
    ! and the parametrization can have different variables switched on
    ! and off.
    ! each of these is a parametrization for a microphysical process.
      lmicro = .true.
    !LMICRO=.false.
 
      if ( lmicro ) then
      ! turn off microphysics
 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                         AUTOCONVERSION                     ! 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    

    !warm clouds
    
        do j = jstart , jend
          do i = istart , iend
            ! if liquid > 0 then B_lv=0, no source of liquid from vapor
            if ( zqxx(j,i,k,iqql) > d_zero ) then    !why if B=0 by default?
              zsolqb(j,i,iqql,iqqv) = d_zero 
              select case (kautoconv)
              case (1)          !Klein & Pincus (2000)
                zsolqb(j,i,iqqr,iqql) = zsolqb(j,i,iqqr,iqql) + &
                      ptsphy*zauto_rate_klepi * (zqxx(j,i,k,iqql)**(3.3D0))
                zsolqa(j,i,iqqr,iqql) = d_zero
              case (2)          ! Khairoutdinov and Kogan (2000)
                zsolqb(j,i,iqqr,iqql) = zsolqb(j,i,iqqr,iqql) + &
                      ptsphy*zauto_rate_khair *         &
                      (zqxx(j,i,k,iqql)**(zauto_expon_khair))
                case (3)          !Kessler(1969)
                  if ( zqxx(j,i,k,iqql) > zautocrit ) then
                    zsolqa(j,i,iqqr,iqql) = zsolqa(j,i,iqqr,iqql) - &
                                            zauto_rate_kesl*zautocrit
                    zsolqa(j,i,iqql,iqqr) = zsolqa(j,i,iqql,iqqr) + &
                                            zauto_rate_kesl*zautocrit
                    zsolqb(j,i,iqqr,iqql) = zsolqb(j,i,iqqr,iqql) + &
                                            ptsphy*zauto_rate_kesl
                  else
                    zsolqb(j,i,iqqr,iqql) = d_zero
                    zsolqa(j,i,iqql,iqqr) = d_zero
                  end if
                case (4)           !Sundqvist
                  zsolqb(j,i,iqqr,iqql) = zsolqb(j,i,iqqr,iqql) + & 
                         ptsphy*zauto_rate_sundq* &
                         (d_one-dexp(-(zqxx(j,i,k,iqql)/zautocrit)**2))
                  zsolqa(j,i,iqqr,iqql) = zsolqa(j,i,iqqr,iqql) + d_zero
              end select
            end if !(ZQX(JL,JK,IQL)>0.0)

    !cold clouds
    !  Snow Autoconversion rate follow Lin et al. 1983
            if (zt(j,i,k) <=  zt0) then
              if (zqxx(j,i,k,iqqi) > zepsec) then
                alpha1 = ptsphy*1.0E-3*exp(0.025*(zt(j,i,k)-zt0))
                zsolqa(j,i,iqqs,iqqi)=zsolqa(j,i,iqqs,iqqi)-alpha1*qi0
                zsolqa(j,i,iqqi,iqqs)=zsolqa(j,i,iqqi,iqqs)+alpha1*qi0
                zsolqb(j,i,iqqs,iqqi)=zsolqb(j,i,iqqs,iqqi)+ptsphy*alpha1
              end if
            end if
   
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                         EVAPORATION                        ! 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
       
            zsolqb(j,i,iqqv,iqqr) = zsolqb(j,i,iqqv,iqqr) + &
                                ptsphy*kevap*max((satvp(j,i,k)-zqxx(j,i,k,iqqv)),d_zero)
    ! IF (ZQX(JL,JK,IQR)>0.0.AND.ZQX(JL,JK,IQV)<SATVP(JL,JK)) THEN
    !         ZSOLQB(JL,IQV,IQR)=KEVAP*(SATVP(JL,JK)-ZQX(JL,JK,IQV))
    ! ELSE
    !        ZSOLQB(JL,IQV,IQR)=0.0
    ! END IF
          end do
        end do  
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                         FREEZING                           ! 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !Freezing of rain. 
    !All rain freezes in a timestep if the temperature is below 0 C
    !calculate sublimation latent heat
        do j = jstart , jend
          do i = istart , iend
            slht=(2834.1 -0.29*zt(j,i,k)-0.004*(zt(j,i,k)**2)) !sublimation latent heat
          end do
        end do

        ralsdcp=slht/rcpd
        ralvdcp=vlht/rcpd
        zrldcp  = 1.0/(ralsdcp-ralvdcp)

        do j = jstart, jend
          do i = istart, iend
            zfrzmax(j,i) = max((zt0-zt(j,i,k))*zrldcp,d_zero)
          end do
        end do

          
        do j = jstart, jend
          do i = istart, iend
            if (zfrzmax(j,i) > zepsec .AND. zqxfg(j,i,iqqr) > zepsec) then
              zfrz=min(zqxfg(j,i,iqqr),zfrzmax(j,i))
              zsolqa(j,i,iqqs,iqqr)=zsolqa(j,i,iqqs,iqqr) + zfrz
              zsolqa(j,i,iqqr,iqqs)=zsolqa(j,i,iqqr,iqqs) - zfrz
            end if 
          end do
        end do
          
    !Freezing of liquid
        do j = jstart, jend
          do i = istart, iend
            if (zfrzmax(j,i) > zepsec .AND. zqxfg(j,i,iqql) > zepsec) then
              zfrz=min(zqxfg(j,i,iqql),zfrzmax(j,i))
              zsolqa(j,i,iqqi,iqql)=zsolqa(j,i,iqqi,iqql) + zfrz
              zsolqa(j,i,iqql,iqqi)=zsolqa(j,i,iqql,iqqi) - zfrz
            end if
          end do
        end do


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                         MELTING                                  !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !The melting of ice and snow are treated explicitly.
    !First water and ice saturation are found
    !---------------------------------------------
    ! ice saturation T<273K
    ! liquid water saturation for T>273K
    !---------------------------------------------

        do j = jstart, jend
          do i = istart, iend 
            zphases = phases(zt(j,i,k))    
            zfoeew(j,i,k)=min(foeewm(zt(j,i,k))/pres(j,i,k),0.5)
            zfoeew(j,i,k)=min(0.5,zfoeew(j,i,k))
            zqsice(j,i,k)=zfoeew(j,i,k)/(1.0-retv*zfoeew(j,i,k))
     
           !----------------------------------
           ! liquid water saturation
           !----------------------------------
           zfoeeliq(j,i,k)=min(foeeliq(zt(j,i,k))/pres(j,i,k),0.5)
           zqsliq(j,i,k)=zfoeeliq(j,i,k)
           zqsliq(j,i,k)=zqsliq(j,i,k)/(1.0-retv*zqsliq(j,i,k))
   
          end do
        end do
        
    ! MELTING OF SNOW and ICE
        do j = jstart, jend
          do i = istart, iend
            zicetot(j,i)=zqxfg(j,i,iqqi)+zqxfg(j,i,iqqs)
            zmeltmax(j,i)=0.0 
            if (zicetot(j,i) > zepsec .and. zt(j,i,k) > zt0) then
              zsubsat = max(zqsice(j,i,k)-zqx(j,i,k,iqqv),0.0)
   
              ! Calculate difference between dry-bulb (zt)  and the temperature
              ! at which the wet-bulb=0degC  
              ! Melting only occurs if the wet-bulb temperature >0
              ! i.e. warming of ice particle due to melting > cooling
              ! due to evaporation.     
              ztdiff = zt(j,i,k)-zt0-zsubsat*(ztw1+ztw2*(pres(j,i,k)-ztw3)-ztw4*(zt(j,i,k)-ztw5))
              ! Ensure ZCONS1 is positive so that ZMELTMAX=0 if ZTDMTW0<0
              zcons1 = abs(ptsphy*(1.0+0.5*ztdiff)/rtaumel)  
              zmeltmax(j,i) = max(ztdiff*zcons1*zrldcp,0.0)
            end if
          end do
        end do    

        ! Loop over frozen hydrometeors (ice, snow)
        do n = 1, nqx
          if (kphase(n) == 2) then
            m = imelt(n) 
            do j = jstart, jend
              do i = istart, iend
                if (zmeltmax(j,i)>zepsec .and. zicetot(j,i) > zepsec) then 
                  zphases = zqxfg(j,i,n)/zicetot(j,i)
                  zmelt = min(zqxfg(j,i,n),zphases*zmeltmax(j,i))
                  !zqxfg(j,i,n)     = zqxfg(j,i,n)-zmelt
                  !zqxfg(j,i,n)     = zqxfg(j,i,n)+zmelt
                  !zsolqa(j,i,m,n) = zsolqa(j,i,m,n)+zmelt
                  !zsolqa(j,i,n,m) = zsolqa(j,i,n,m)-zmelt
                end if
              end do 
            end do
          end if
        end do
         
 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                         EVAPORATION                              !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
    ! Evaporate very small amounts of liquid and ice
        do j = jstart, jend
          do i = istart, iend
            if ( zqxx(j,i,k,iqql) < rlmin ) then
              zsolqa(j,i,iqqv,iqql) = zqxx(j,i,k,iqql)
              zsolqa(j,i,iqql,iqqv) = -zqxx(j,i,k,iqql)
            end if
 
            if ( zqxx(j,i,k,iqqi) < rlmin ) then
              zsolqa(j,i,iqqv,iqqi) = zqxx(j,i,k,iqqi)
              zsolqa(j,i,iqqi,iqqv) = -zqxx(j,i,k,iqqi)
            end if
        
          end do ! IY
        end do !JPX
      
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                 DETRAINMENT FROM CONVECTION
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
     
        do j = jstart, jend
          do i = istart, iend
            zphases=phases(zt(j,i,k))
            !zqdetr(j,i,k)=zqdetr(j,i,k)*zdtgdp(j,i)
            if(zqdetr(j,i,k)>rlmin) then
              zconvsrce(j,i,iqql) = zphases*zqdetr(j,i,k)
              !zconvsrce(j,i,iqqi) = (1.0-zphases)*zqdetr(j,i,k)
              !zsolqa(j,i,iqql,iqql) = zsolqa(j,i,iqql,iqql)+ zphases*zqdetr(j,i,k)*zdtgdp(j,i)
              zsolqa(j,i,iqql,iqql) = zsolqa(j,i,iqql,iqql)+zconvsrce(j,i,iqql)*zdtgdp(j,i)
              !zconvsrce(j,i,iqql)
              zsolqa(j,i,iqqi,iqqi) = zsolqa(j,i,iqqi,iqqi)+(1.0-zphases)*zqdetr(j,i,k)*zdtgdp(j,i)
              !zconvsrce(j,i,iqqi)
            else
              zqdetr(j,i,k)=0.0
            end if 
          end do
        end do 
   
      end if ! lmicro
 
  !!!!!!!!!!!!! sedimentation/falling of microphysical species
  ! do n = 1 , nqx
  !   do j = jstart , jend
  !     do i = istart , iend
  !       zfallsink(j,i,n) = ptsphy*vqx(j,i,k,n)/dzeta(j,i,k)
  !       if ( k > 1 ) then
  !         zfallsrce(j,i,n) = ptsphy*vqx(j,i,k,n) * &
  !                          zqxx(j,i,k-1,n)/dzeta(j,i,k)
  !       end if
  !     end do
  !   end do
  ! end do
      

  ! source from layer above
      do n = 1 , nqx
        do j = jstart , jend
          do i = istart , iend   
            if (k > 1) then
              zfallsrce(j,i,n) = zpfplsx(j,i,k,n)*zdtgdp(j,i)
              zsolqa(j,i,n,n)  = zsolqa(j,i,n,n)+zfallsrce(j,i,n)
              zqxfg(j,i,n)   = zqxfg(j,i,n)+zfallsrce(j,i,n)
            end if
            zfall=vqx(j,i,k,n)*rhob3d(j,i,k)
            zfallsink(j,i,n)=zdtgdp(j,i)*zfall 
          end do
        end do
      end do

  ! sink to next layer, constant fall speed
  !-----------------------------------------------
  ! solvers...
  !-----------------------------------------------
 
  !----------------------------------------------------------
  !truncate sum of explicit sinks to size of bin
  ! this approach is inaccurate, but conserves -
  ! prob best can do with explicit (i.e. not implicit!) terms
  !----------------------------------------------------------
      jindex1(:,:,:,:) = 0
      zsinksum(:,:,:) = d_zero
      do n = 1 , nqx
        do jn = 1 , nqx
          do j = jstart , jend
            do i = istart , iend
              if ( zsolqa(j,i,jn,n) < d_zero ) then
                jindex1(j,i,jn,n) = 1
               zsinksum(j,i,jn) = zsinksum(j,i,jn) - zsolqa(j,i,jn,n)
             end if
            end do
          end do
        end do
      end do
 
      do n = 1 , nqx
        do j = jstart , jend
          do i = istart , iend
            zratio(j,i,n) = d_one   !JPRB
            if ( zsinksum(j,i,n) > d_zero ) then
              zratio(j,i,n) = zqxx(j,i,k,n)/zsinksum(j,i,n)
              zratio(j,i,n) = max(min(zratio(j,i,n),d_one),d_zero)
            end if
          end do
        end do
      end do
      do n = 1 , nqx
        do jn = 1 , nqx
          do j = jstart , jend
            do i = istart , iend
              if ( jindex1(j,i,jn,n) == 1 ) then
                zsolqa(j,i,jn,n) = zsolqa(j,i,jn,n)*zratio(j,i,jn)
                zsolqa(j,i,n,jn) = zsolqa(j,i,n,jn)*zratio(j,i,jn)
              end if
            end do
          end do
        end do
      end do
  ! set the LHS of equation
      do n = 1 , nqx
        do jn = 1 , nqx
          do j = jstart , jend
            do i = istart , iend
            ! diagonals: microphysical sink terms+transport
              if ( jn == n ) then
                zqlhs(j,i,jn,n) = d_one + zfallsink(j,i,n)
                do jo = 1 , nqx
                  zqlhs(j,i,jn,n) = zqlhs(j,i,jn,n) + zsolqb(j,i,jo,jn)
                end do
                ! non-diagonals: microphysical source terms
              else
              ! here is the delta T - missing from doc.
              zqlhs(j,i,jn,n) = -zsolqb(j,i,jn,n)
              end if
            end do
          end do
        end do
      end do
 
  ! set the RHS of equation
      do n = 1 , nqx
        do j = jstart , jend
          do i = istart , iend
            !   sum the explicit source and sink
            zexplicit = d_zero
            do jn = 1 , nqx
              ! negative, since summed over 1st index
              zexplicit = zexplicit - zsolqa(j,i,jn,n)
            end do
          zqxn(j,i,k,n) = zqxx(j,i,k,n) + zexplicit
          ! note *** density missing*** for the moment...
          zfluxq(j,i,n) = zfallsrce(j,i,n)+ zconvsrce(j,i,n)   !(delz at k-1)
          zqxn(j,i,k,n) = zqxn(j,i,k,n) + zfluxq(j,i,n)
          !    ZQXN(JL,JM)=MAX(ZQXN(JL,JM),0.0)
          end do
        end do
      end do
 
  ! solve by LU decomposition:
  ! dummy test
  !ZQLHS(:,:,:)=0.0
  !DO JM=1,NQX
  !  ZQLHS(1,JM,JM)=1.0
  !ENDDO
  !ZQLHS(1,3,3)=ZQLHS(1,3,3)+0.47
  !ZQLHS(1,5,3)=ZQLHS(1,5,3)-0.47
  !print*,'prima subroutine'
      call ludcmp(jstart,jend,istart,iend,zqlhs,jindex2)
      call lubksb(jstart,jend,istart,iend,zqlhs,jindex2,zqxn)
 
  ! calculate fluxes in and out of box for conservation of TL
      do n = 1 , nqx
        do j = jstart , jend
          do i = istart , iend
            !   Note ***explicit PLUDE*** fluxes missing for moment
            zfluxq(j,i,n) = zfluxq(j,i,n) - ptsphy*vqx(j,i,k,n) * &
                            zqxn(j,i,k,n)/dzeta(j,i,n)
            if ( kphase(n) == 1 ) then
              zt(j,i,k) = zt(j,i,k) - vlht*zfluxq(j,i,n)/rcpd
              ! liquid terms
              zt(j,i,k) = zt(j,i,k)+vlht*(zqxn(j,i,k,n)-zqxx(j,i,k,n))/rcpd
            end if
            if ( kphase(n) == 2 ) then
              zt(j,i,k) = zt(j,i,k) + slht*zfluxq(j,i,n)/rcpd  !ice terms
            end if
          end do
        end do
      end do
      

      ! generalized precipitation flux
      do n=1,nqx
        do j = jstart , jend
          do i = istart , iend
            zpfplsx(j,i,k+1,n) = zfallsink(j,i,n)*zqxn(j,i,k,n)*zrdtgdp(j,i)
          end do
        end do
      end do

  ! copy the result into the Q array
  ! do n = 1 , nqx
  !   do j = jstart , jend
  !     zqx(j,k,n) = zqxn(j,k,n)
  !   end do
  ! end do

  ! recalcalate ZTL to check for conservation
  !ztln(jci1:jci2,ici1:ici2,1:k) = ztl(jci1:jci2,ici1:ici2,1:k)

      do j = jstart, jend
        do i = istart, iend 
          ztln(j,i,k) = ztl(j,i,k)
        end do
      end do

      do n = 1 , nqx
        do j = jstart , jend
          do i = istart , iend
            if ( kphase(n) == 1 ) then
              ztln(j,i,k) = ztln(j,i,k) - vlht*zqxx(j,i,k,n)/rcpd
            end if
            if ( kphase(n) == 2 ) then
              ztln(j,i,k) = ztln(j,i,k) - slht*zqxx(j,i,k,n)/rcpd
            end if
          end do
        end do
      end do
      do j = jstart , jend
        do i = istart , iend
        !print*,'ztln',ztln(j,i,k),'differenza',ztln(j,i,k)-ztl(j,i,k)             
        !print*,'detr',detr(2,2,:)
          ztten(j,i,k)=ztten(j,i,k)+0.001
        !print*,'ztl',ztl(2,2,18)
        !print*,'ztln',ztln(2,2,18) 
        end do
      end do
      do n = 1 , nqx
        do j = jstart , jend
          do i = istart , iend
            zqxten(j,i,k,n) = zqxn(j,i,k,n)-zqxx(j,i,k,n)
            wvflux(j,i,k,n) = rhob3d(j,i,k)*zqxn(j,i,k,n)*vqx(j,i,k,n)
          end do
        end do 
      end do 
      do j = jstart , jend
        do i = istart , iend
          difftemp = ztln(j,i,k) - ztl(j,i,k)
        end do
      end do 
  ! print*,'fine loop ztten'

  !do j = jstart , jend
  !do i = istart , iend
  !ztten(j,i,k)=ztten(j,i,k)+0.1
  !print*,'ztten2',ztten(j,i,k)
  !end do
  !end do
    
  ! calculate cloud total water at the end
      do j = jstart , jend
        do i = istart , iend
          totalw2(j,i,k)= zqxx(j,i,k,iqqv)+zqxx(j,i,k,iqql)+zqxx(j,i,k,iqqr) + &
                          zqxx(j,i,k,iqqi)+zqxx(j,i,k,iqqs)
        end do
      end do
      w2=sum(totalw2)
    end do   ! kz : end of vertical loop
!print*,' diff qt=', w2-w1  
!print*, 'ultimo step'
 
contains

     real(rk8) function phases(zt)
       implicit none
       real(rk8) , intent(in):: zt
       real(rk8) , parameter :: rtice =  250.16                     !zt0 - 23.0
       real(rk8) , parameter :: rtwat_rtice_r=1./23.0
       real(rk8) , parameter :: zt0 = 273.16
       phases = min(1.0,((max(rtice,min(zt0,zt))-rtice)*rtwat_rtice_r)**2)
     end function phases

     real(rk8) function foeewm(zt)
       implicit none
       real(rk8) , intent(in):: zt
       real(rk8) , parameter :: r2es =  611.21                     !
       real(rk8) , parameter :: r3les = 17.502
       real(rk8) , parameter :: r3ies = 22.587
       real(rk8) , parameter :: r4les = 32.19
       real(rk8) , parameter :: r4ies = -0.7
       real(rk8) , parameter :: zt0 = 273.16
       foeewm = r2es*(phases(zt)*exp(r3les*(zt-zt0)/(zt-r4les))+&
       &(1.0-phases(zt))*exp(r3ies*(zt-zt0)/(zt-r4ies)))
       end function foeewm

     real(rk8) function foeeliq(zt)
       implicit none
       real(rk8) , intent(in):: zt
       real(rk8) , parameter :: rtice =  250.16                     !zt0 - 23.0
       real(rk8) , parameter :: rtwat_rtice_r=1./23.0
       real(rk8) , parameter :: zt0 = 273.16
       foeeliq = min(1.0,((max(rtice,min(zt0,zt))-rtice)*rtwat_rtice_r)**2)
     end function foeeliq

     real(rk8) function foeeice(zt)
       implicit none
       real(rk8) , intent(in):: zt
       real(rk8) , parameter :: rtice =  250.16                     !zt0 - 23.0
       real(rk8) , parameter :: rtwat_rtice_r=1./23.0
       real(rk8) , parameter :: zt0 = 273.16
       foeeice = min(1.0,((max(rtice,min(zt0,zt))-rtice)*rtwat_rtice_r)**2)
     end function foeeice

     real(rk8) function dqsatdtc(zt,satc)
       implicit none
       real(rk8) , intent(in) :: satc , zt
       dqsatdtc = satc*(6150.0D0/(zt**2))
     end function dqsatdtc
      
     real(rk8) function dqsatdtw(zt,satw)
       implicit none
       real(rk8) , intent(in) :: satw , zt
       dqsatdtw = satw*(4097.99D0/((zt-32.15D0)**2))
     end function dqsatdtw

     real(rk8) function satc(zt,xp)
       implicit none
       real(rk8) , intent(in) :: xp , zt
       ! saturation mixing ratio for ice in hPa (A Description of the Fifth-Generation Penn State/NCAR
       !Mesoscale Model (MM5)Georg A. Grell, Jimy Dudhia, David R. Stauffer, NCAR TECHNICAL NOTE Dec 1994)
       satc = (3.79D0/xp)*dexp(22.514D0-6150.0D0/zt) 
     end function satc

     real(rk8) function satw(zt,xp)
       implicit none
       real(rk8) , intent(in) :: xp , zt
       !  ! saturation mixing ratio for water in hPa (A Description of the Fifth-Generation Penn State/NCAR
       !Mesoscale Model (MM5)Georg A. Grell, Jimy Dudhia, David R. Stauffer, NCAR TECHNICAL NOTE Dec 1994)
       satw = (3.79D0/xp)*dexp(17.67D0*(zt-tzero)/(zt-29.65D0))
     end function satw
     
  end subroutine microphys

  subroutine lubksb(jstart,jend,istart,iend,aam,indx,bbm)
    implicit none
    integer(ik4) , intent(in) :: jstart , jend , istart , iend
    real(rk8) , pointer , intent(in) , dimension(:,:,:,:) :: aam
    integer(ik4) , pointer , intent(in) , dimension(:,:,:) :: indx
    real(rk8) , pointer , intent(inout) , dimension(:,:,:,:) :: bbm
    integer(ik4) :: i , j , ii , jj , k , ll , m
    real(rk8) :: xsum

    ! SOLVES THE SET OF N LINEAR EQUATIONS A * X = B.
    ! HERE A IS INPUT, NOT AS THE MATRIX A BUT RATHER AS
    ! ITS LU DECOMPOSITION, DETERMINED BY THE ROUTINE LUDCMP.
    ! INDX IS INPUT AS THE PERMUTATION VECTOR RETURNED BY LUDCMP.
    ! B(1:N) IS INPUT AS THE RIGHT-HAND SIDE VECTOR B,
    ! AND RETURNS WITH THE SOLUTION VECTOR X. A, N, NP,
    ! AND INDX ARE NOT MODI ED BY THIS ROUTINE AND CAN BE
    ! LEFT IN PLACE FOR SUCCESSIVE CALLS WITH DI
     
    do k = 1 , kz
      do j = jstart , jend
        do i = istart , iend
          ii = 0
          ! WHEN II IS SET TO A POSITIVE VALUE, IT WILL BECOME
          ! THE INDEX OF THE  RST NONVANISHING ELEMENT OF B.
          ! WE NOW DO THE FORWARD SUBSTITUTION, EQUATION (2.3.6).
          ! THE ONLY NEW WRINKLE IS TO UNSCRAMBLE THE PERMUTATION AS WE GO.
          do m = 1 , nqx
            ll = indx(j,i,m)
            xsum = bbm(j,i,k,ll)
            bbm(j,i,k,ll) = bbm(j,i,k,m)
            if ( ii /= 0 ) then
              do jj = ii , m - 1
                xsum = xsum - aam(j,i,m,jj)*bbm(j,i,k,jj)
              end do
            else if ( dabs(xsum) > dlowval ) then
                ii = m
            end if
            bbm(j,i,k,m) = xsum
          end do
          do m = nqx , 1 , -1 ! NOW WE DO THE BACKSUBSTITUTION, EQUATION (2.3.7)
            xsum = bbm(j,i,k,m)
            do jj = m + 1 , nqx
              xsum = xsum - aam(j,i,m,jj)*bbm(j,i,k,jj)
            end do
            ! STORE A COMPONENT OF THE SOLUTION VECTOR X.
            bbm(j,i,k,m) = xsum/aam(j,i,m,m)
          end do
        end do
      end do
    end do
  end subroutine lubksb

  subroutine ludcmp(jstart,jend,istart,iend,aam,indx)
    implicit none
!
    integer(ik4) , intent(in) :: jstart , jend , istart , iend
    real(rk8) , pointer , intent(inout) , dimension(:,:,:,:) :: aam
    integer(ik4) , pointer , intent(out) , dimension(:,:,:) :: indx
!
    real(rk8) :: aamax , dum , xsum
    integer(ik4) :: i , j , k , imax , n , m
    real(rk8) , dimension(nmax) :: vv
!
    do j = jstart , jend
      do i = istart , iend
        do m = 1 , nqx !LOOP OVER ROWS TO GET THE IMPLICIT SCALING INFORMATION.
          aamax = d_zero
          do n = 1 , nqx
            if ( dabs(aam(j,i,m,n)) > aamax ) aamax = dabs(aam(j,i,m,n))
          end do
          if ( dabs(aamax) < dlowval ) then
            call fatal(__FILE__,__LINE__,'SINGULAR MATRIX')
          end if ! SINGULAR MATRIX
          vv(m) = d_one/aamax !SAVE THE SCALING.
        end do
        do n = 1 , nqx
          ! THIS IS THE LOOP OVER COLUMNS OF CROUT S METHOD.
          do m = 1 , n - 1
            !THIS IS EQUATION (2.3.12) EXCEPT FOR I = J.
            xsum = aam(j,i,m,n)
            do k = 1 , m - 1
              xsum = xsum - aam(j,i,m,k)*aam(j,i,k,n)
            end do
            aam(j,i,m,n) = xsum
          end do
          aamax = d_zero ! INITIALIZE FOR THE SEARCH FOR LARGEST PIVOT ELEMENT.
          do m = n , nqx ! THIS IS I = J OF EQUATION (2.3.12)
            ! AND I = J+1. . .N OF EQUATION (2.3.13).
            xsum = aam(j,i,m,n)
            do k = 1 , n - 1
              xsum = xsum - aam(j,i,m,k)*aam(j,i,k,n)
            end do
            aam(j,i,m,n) = xsum
            dum = vv(m)*dabs(xsum)   ! FIGURE OF MERIT FOR THE PIVOT.
            if ( dum >= aamax ) then ! IS IT BETTER THAN THE BEST SO FAR?
              imax = m
              aamax = dum
            end if
          end do
          if ( n /= imax ) then
            ! DO WE NEED TO INTERCHANGE ROWS?
            do k = 1 , nqx ! YES, DO SO...
              dum = aam(j,i,imax,k)
              aam(j,i,imax,k) = aam(j,i,n,k)
              aam(j,i,n,k) = dum
            end do
            ! D=-D !...AND CHANGE THE PARITY OF D.
            vv(imax) = vv(n) ! ALSO INTERCHANGE THE SCALE FACTOR.
          end if
          indx(j,i,n) = imax
          if ( dabs(aam(j,i,n,n)) < dlowval ) aam(j,i,n,n) = dlowval
          if ( n /= nqx ) then
            dum = d_one/aam(j,i,n,n)
            do m = n + 1 , nqx
              aam(j,i,m,n) = aam(j,i,m,n)*dum
            end do
          end if
        end do
      end do
    end do
  end subroutine ludcmp

end module mod_cloud_s1
