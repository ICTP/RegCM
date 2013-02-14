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
  use mod_runparams , only : sigma
  use mod_runparams , only : dt
  use mod_pbl_common
  use mod_constants
  use mod_precip , only : fcc
  use mod_runparams , only : rtsrf
  use mod_service
  private

  integer(ik4) , parameter :: nmax = 5                           ! number of progn.variables 
  integer(ik4) , parameter :: nqx = 5                             
  real(rk8) :: zfall                                             ! constant fall speed   
  real(rk8) :: zqtmst                                            ! 1/dt
  real(rk8) , pointer , dimension(:,:,:) :: pres                 ! from atms
  real(rk8) , pointer , dimension(:,:,:) :: zt                   ! from atms
  real(rk8) , pointer , dimension(:,:,:) :: zeta                 ! from atms
  real(rk8) , pointer , dimension(:,:,:) :: dzeta                ! from atms
  real(rk8) , pointer , dimension(:,:,:,:) :: zqxx               ! from atms
  real(rk8) , pointer , dimension(:,:,:) :: rhob3d               ! from atms 
  real(rk8) , pointer , dimension(:,:,:) :: omega                ! from atms 
  real(rk8) , pointer , dimension(:,:,:) :: heatrt               ! radiation heat rate
  real(rk8) , pointer , dimension(:,:,:) :: ztten                ! tendency of temperature
  real(rk8) , pointer , dimension(:,:,:) :: qdetr                ! conv. detrained water
  real(rk8) , pointer , dimension(:,:,:,:) :: zqxten             ! tendency of zqx
  real(rk8) , pointer , dimension(:,:) :: psf , rainnc, lsmrnc
  
  public :: allocate_mod_cloud_s1 , init_cloud_s1 , microphys

! Total water and enthalpy budget diagnostics variables
  integer(ik4) , pointer , dimension(:) :: kphase                ! marker for water phase of each species
                                                                 ! 0=vapour, 1=liquid, 2=ice  
  integer(ik4) , pointer , dimension(:) :: imelt                 ! marks melting linkage for ice categories
                                                                 ! ice->liquid, snow->rain
  real(rk8) :: zalfaw
  real(rk8) :: zphases, zmelt, zice
  real(rk8) :: ztmpl,ztmpi, ztnew,zqe,zrain
  real(rk8) , pointer , dimension(:,:,:) :: papf
  real(rk8) , pointer , dimension(:,:,:):: zsumh0,zsumq0
  real(rk8) , pointer , dimension(:,:,:) :: zsumh1,zsumq1
  real(rk8) , pointer , dimension(:,:,:) :: zerrorq,zerrorh
  real(rk8) , pointer , dimension(:,:,:):: ztentkeep
  real(rk8) , pointer , dimension(:,:,:,:) :: ztenkeep
 
! Mass variables 
  real(rk8) , pointer , dimension(:,:) :: zdp                     ! dp
  real(rk8) , pointer , dimension(:,:) :: zgdp                    ! g/dp
  real(rk8) , pointer , dimension(:,:) :: zdtgdp                  ! dt * g/dp
  real(rk8) , pointer , dimension(:,:) :: zrdtgdp                 ! dp / (dt * g)
  
! Microphysics 
  integer(ik4) , pointer , dimension(:,:,:)   :: jindex2          ! index variable
  integer(ik4) , pointer , dimension(:,:,:,:) :: jindex1          ! index variable
  real(rk8) , pointer , dimension(:,:) :: zfrzmax
  real(rk8) , pointer , dimension(:,:) :: zicetot
  real(rk8) , pointer , dimension(:,:) :: zmeltmax
  real(rk8) , pointer , dimension(:,:) :: prcflxw
  real(rk8) , pointer , dimension(:,:) :: prcflxc
  real(rk8) , pointer , dimension(:,:,:) :: dqsatdt
  real(rk8) , pointer , dimension(:,:,:) :: satvp
  real(rk8) , pointer , dimension(:,:,:) :: satice          
! for sedimentation source/sink terms
  real(rk8) , pointer , dimension(:,:,:) :: zfallsink
  real(rk8) , pointer , dimension(:,:,:) :: zfallsrce
! for convection detrainment source and subsidence source/sink terms
  real(rk8) , pointer , dimension(:,:,:) :: zconvsrce
  real(rk8) , pointer , dimension(:,:,:) :: zconvsink
  
  real(rk8) , pointer , dimension(:,:) :: zliqcld
  real(rk8) , pointer , dimension(:,:) :: zicecld
  real(rk8) , pointer , dimension(:,:,:) :: zfluxq                ! fluxes convergence of species 
  real(rk8) , pointer , dimension(:,:,:) :: zratio
  real(rk8) , pointer , dimension(:,:,:) :: zsinksum
  real(rk8) , pointer , dimension(:,:,:) :: zfoeew                 
  real(rk8) , pointer , dimension(:,:,:) :: zliq                   
  real(rk8) , pointer , dimension(:,:,:) :: zqsliq                ! liquid water saturation
  real(rk8) , pointer , dimension(:,:,:) :: zqsice                ! ice water saturation
  real(rk8) , pointer , dimension(:,:,:) :: zfoeeliq              ! melting 
  real(rk8) , pointer , dimension(:,:,:) :: zpfplsl                !liq+rain sedim flux
  real(rk8) , pointer , dimension(:,:,:) :: zpfplsn                !ice+snow sedim flux
  real(rk8) , pointer , dimension(:,:,:) :: zttendc               !decoupled temperature tendency 
  real(rk8) , pointer , dimension(:,:,:) :: zqdetr                ! detrainment from tiedtke scheme 
  real(rk8) , pointer , dimension(:,:,:,:) :: zvqx      ! fall speeds of three categories
  real(rk8) , pointer , dimension(:,:,:,:) :: zqlhs     ! n x n matrix storing the LHS of implicit solver
  real(rk8) , pointer , dimension(:,:,:,:) :: zsolqa    ! explicit sources and sinks
  real(rk8) , pointer , dimension(:,:,:,:) :: zsolqb    ! implicit sources and sinks
  real(rk8) , pointer , dimension(:,:,:,:) :: zqxtendc  ! decoupled mixing ratios tendency   
  real(rk8) , pointer , dimension(:,:,:,:) :: zpfplsx            ! j,i,n ! generalized precipitation flux
  real(rk8) , public  , pointer, dimension(:,:,:,:) :: zqxn  ! new values for zqxx at time+1
 
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
    call getmem2d(prcflxw,jci1,jci2,ici1,ici2,'clouds1:prcflxw')
    call getmem2d(prcflxc,jci1,jci2,ici1,ici2,'clouds1:prcflxc')
    call getmem3d(zttendc,jci1,jci2,ici1,ici2,1,kz,'clouds1:zttendc')
    call getmem3d(zconvsrce,jci1,jci2,ici1,ici2,1,nqx,'clouds1:zconvsrce')
    call getmem3d(zconvsink,jci1,jci2,ici1,ici2,1,nqx,'clouds1:zconvsink')
    call getmem3d(zfoeew,jci1,jci2,ici1,ici2,1,kz,'clouds1:zfoeew') 
    call getmem3d(zfoeeliq,jci1,jci2,ici1,ici2,1,kz,'clouds1:zfoeeliq')
    call getmem3d(zqsice,jci1,jci2,ici1,ici2,1,kz,'clouds1:zqsice')
    call getmem3d(zqsliq,jci1,jci2,ici1,ici2,1,kz,'clouds1:zqsliq')
    call getmem3d(jindex2,jci1,jci2,ici1,ici2,1,nqx,'clouds1:jindex2')
    call getmem3d(zfallsink,jci1,jci2,ici1,ici2,1,nqx,'clouds1:zfallsink')
    call getmem3d(zfallsrce,jci1,jci2,ici1,ici2,1,nqx,'clouds1:zfallsrce')
    call getmem3d(zfluxq,jci1,jci2,ici1,ici2,1,nqx,'clouds1:zfluxq')
    call getmem3d(zratio,jci1,jci2,ici1,ici2,1,nqx,'clouds1:zratio')
    call getmem3d(zsinksum,jci1,jci2,ici1,ici2,1,nqx,'clouds1:zsinksum')
    call getmem3d(dqsatdt,jci1,jci2,ici1,ici2,1,kz,'clouds1:dqsatdt')
    call getmem3d(satvp,jci1,jci2,ici1,ici2,1,kz,'clouds1:satvp')
    call getmem3d(satice,jci1,jci2,ici1,ici2,1,kz,'clouds1:satice')
    call getmem3d(zpfplsl,jci1,jci2,ici1,ici2,1,kz+1,'clouds1:zpfplsl')
    call getmem3d(zpfplsn,jci1,jci2,ici1,ici2,1,kz+1,'clouds1:zpfplsn')
    call getmem3d(ztentkeep,jci1,jci2,ici1,ici2,1,kz,'clouds1:ztentkeep')
    call getmem3d(zsumq0,jci1,jci2,ici1,ici2,1,kz,'clouds1:zsumq0')
    call getmem3d(zsumh0,jci1,jci2,ici1,ici2,1,kz,'clouds1:zsumh0')
    call getmem3d(zsumq1,jci1,jci2,ici1,ici2,1,kz,'clouds1:zsumq1')
    call getmem3d(zsumh1,jci1,jci2,ici1,ici2,1,kz,'clouds1:zsumh1')
    call getmem3d(zerrorq,jci1,jci2,ici1,ici2,1,kz,'clouds1:zerrorq')
    call getmem3d(zerrorh,jci1,jci2,ici1,ici2,1,kz,'clouds1:zerrorh')
    call getmem3d(papf,jci1,jci2,ici1,ici2,1,kz+1,'clouds1:papf')
    call getmem3d(zliq,jci1,jci2,ici1,ici2,1,kz+1,'clouds1:zliq')
    call getmem2d(zliqcld,jci1,jci2,ici1,ici2,'clouds1:zliqcld')
    call getmem2d(zicecld,jci1,jci2,ici1,ici2,'clouds1:zicecld')
    call getmem3d(zqdetr,jci1,jci2,ici1,ici2,1,kz,'clouds1:zqdetr')
    call getmem4d(zqxtendc,jci1,jci2,ici1,ici2,1,kz,1,nqx,'clouds1:zqxtendc')
    call getmem4d(zvqx,jci1,jci2,ici1,ici2,1,kz,1,nqx,'clouds1:zvqx')
    call getmem4d(zqxn,jce1,jce2,ice1,ice2,1,kz,1,nqx,'clouds1:zqxn')
    call getmem4d(zqlhs,jci1,jci2,ici1,ici2,1,nqx,1,nqx,'clouds1:zqlhs')
    call getmem4d(zsolqa,jci1,jci2,ici1,ici2,1,nqx,1,nqx,'clouds1:zsolqa')
    call getmem4d(zsolqb,jci1,jci2,ici1,ici2,1,nqx,1,nqx,'clouds1:zsolqb')
    call getmem4d(jindex1,jci1,jci2,ici1,ici2,1,nqx,1,nqx,'clouds1:jindex1')
    call getmem4d(zpfplsx,jce1,jce2,ice1,ice2,1,kz+1,1,nqx,'clouds1:zpfplsx')
    call getmem4d(ztenkeep,jce1,jce2,ice1,ice2,1,kz,1,nqx,'clouds1:ztenkeep')
  end subroutine allocate_mod_cloud_s1

  subroutine init_cloud_s1(atms,aten,heatrt,sfs,q_detr,pptnc)
    implicit none
    type(slice) , intent(in) :: atms
    type(atmstate) , intent(in) :: aten
    type(surfstate) , intent(in) :: sfs
    real(rk8) , pointer , intent(in), dimension(:,:,:) :: heatrt, q_detr
    real(rk8) , pointer , dimension(:,:) :: pptnc
    call assignpnt(atms%pb3d,pres)
    call assignpnt(atms%tb3d,zt)
    call assignpnt(atms%za,zeta)
    call assignpnt(atms%dzq,dzeta)
    call assignpnt(atms%qxb3d,zqxx)
    call assignpnt(atms%rhob3d,rhob3d)
    call assignpnt(heatrt,radheatrt)
    call assignpnt(aten%qx,zqxten)
    call assignpnt(aten%t,ztten)
    call assignpnt(q_detr,qdetr)
    call assignpnt(sfs%psb,psf)
    call assignpnt(sfs%psb,sfcps)
    call assignpnt(sfs%rainnc,rainnc)
    call assignpnt(pptnc,lsmrnc)

  end subroutine init_cloud_s1

  subroutine microphys(omega,jstart,jend,istart,iend)
    implicit none
    integer(ik4) , intent(in) :: jstart , jend , istart , iend
    integer(ik4) :: i , j , k , n , m 
    integer(ik4) :: iqqi , iqql , iqqr , iqqs , iqqv , jn , jo , kautoconv
    logical :: lmicro, budget
    real(rk8) , pointer , dimension(:,:,:) , intent(in) :: omega             
    real(rk8) :: zexplicit
 
    ! local real variables for autoconversion rate constants
    real(rk8) :: alpha1                                              ! coefficient autoconversion cold cloud
    real(rk8) :: ztmpa
    real(rk8) , parameter :: zauto_rate_khair = 0.355D0              ! microphysical terms
    real(rk8) , parameter :: zauto_expon_khair = 1.47D0
    real(rk8) , parameter :: zauto_rate_sundq = 0.5D-3
    real(rk8) , parameter :: zauto_rate_kesl = 1.D-3                 !giusto!
    real(rk8) , parameter :: zauto_rate_klepi = 0.5D-3
    real(rk8) , parameter :: zautocrit = 5.D-4                       !giusto!
    real(rk8) , parameter :: zepsec = 1.D-10
    real(rk8) , parameter :: qi0 = 1.0D-3                            !g g^(-1)
    real(rk8) , parameter :: retv = 0.60                             !rv/rd-1   rv = 461.5, rd = 287.05              
   
    ! local real constants and variables for condensation
    real(rk8) :: zcond  
    real(rk8) :: zdtdp  

    ! local real constants and variables for freezing
    real(rk8) :: zfrz
    real(rk8) :: zrldcp

    ! local real constants and variables for melting
    real(rk8) :: zsubsat
    real(rk8) :: ztdiff
    real(rk8) :: zcons1

    real(rk8) :: rovcp
     
    ! constant for converting the fluxes unit measures
    real(rk8) :: prainx 
     
    ! local real constants for evaporation
    real(rk8) , parameter :: kevap = 0.100D-02                       ! Raindrop evap rate coef
    real(rk8) , parameter :: rlmin = 1.D-8
    
    ! Numerical fit to wet bulb temperature
    real(rk8) , parameter :: ztw1 = 1329.31
    real(rk8) , parameter :: ztw2 = 0.0074615
    real(rk8) , parameter :: ztw3 = 0.85D5
    real(rk8) , parameter :: ztw4 = 40.637
    real(rk8) , parameter :: ztw5 = 275.0
    real(rk8) , parameter :: rtaumel=11880.0 

#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'microphys'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif


    ! Define species tags
    iqqv = 1    ! vapour
    iqql = 2    ! liquid cloud water
    iqqr = 3    ! rain water
    iqqi = 4    ! ice
    iqqs = 5    ! snow
 
    ! Choose the autoconversion paramaterization 
    !KAUTOCONV = 1 ! Klein & Pincus (2000)
    !kautoconv = 2 ! Khairoutdinov and Kogan (2000)
     KAUTOCONV = 3 ! Kessler (1969)
    ! KAUTOCONV = 4 ! Sundqvist
 
    ! Define species phase, 0=vapour, 1=liquid, 2=ice
    kphase(iqqv) = 0
    kphase(iqql) = 1
    kphase(iqqr) = 1
    kphase(iqqi) = 2
    kphase(iqqs) = 2

    ! Set up melting/freezing index, 
    ! if an ice category melts/freezes, where does it go?
    
    imelt(iqqv) = -99
    imelt(iqql) = iqqi
    imelt(iqqr) = iqqs
    imelt(iqqi) = iqql
    imelt(iqqs) = iqqr

!Define zliq the function for mixed phase
do k = 1 , kz
  do j = jstart , jend
    do i = istart , iend
      zliq(j,i,k) = phases(zt(j,i,k)) 
    end do
  end do
end do

! Total water and enthalpy budget on/off
budget = .true.
                                
!-----------------------------------
! initialization for cloud variables
! -------------------------------------

!Define pressure at full levels (half levels for ECMWF)
do k = 1 , kz+1
  do i = istart , iend
    do j = jstart , jend
      papf(j,i,k) = (sigma(k)*sfcps(j,i)+ptop)*d_1000          ! (Pa)
    end do
  end do
end do

!Convert pressure at half levels in Pa
do k = 1 , kz
  do i = istart , iend
    do j = jstart , jend
      pres(j,i,k) = pres(j,i,k)*d_1000        !   (Pa)
    end do
  end do
end do

!Define a new array for detrainment
do k = 1 , kz
  do i = istart , iend
    do j = jstart , jend
      zqdetr(j,i,k) = qdetr(j,i,k)
    end do
  end do
end do

! Budget reset errors and variables
zerrorq(:,:,:) = d_zero
zerrorh(:,:,:) = d_zero
zsumh0(:,:,:) = d_zero
zsumq0(:,:,:) = d_zero
zsumh1(:,:,:) = d_zero
zsumq1(:,:,:) = d_zero
ztentkeep(:,:,:) = d_zero
ztenkeep(:,:,:,:) = d_zero

! Decouple tendencies
do k = 1 , kz
  do n = 1 , nqx
    do j = jstart , jend
      do i = istart , iend
        zqxtendc(j,i,k,n) = zqxten(j,i,k,n)/psf(j,i)
      end do
    end do
  end do
end do

do k = 1 , kz
  do j = jstart , jend
    do i = istart , iend
      zttendc(j,i,k) = ztten(j,i,k)/psf(j,i)
    end do
  end do
end do

! Record the tendencies
do k = 1 , kz
  do j = jstart , jend
    do i = istart , iend
      do n = 1 , nqx
        ztenkeep(j,i,k,n) = zqxtendc(j,i,k,n)               
      end do
      ztentkeep(j,i,k) = zttendc(j,i,k)
    enddo
  enddo
enddo

zrldcp  = d_zero/(wlhfocp-wlhvocp)                 !

!-------------------------------------
! Initial enthalpy and total water diagnostics 
!-------------------------------------

if ( budget ) then
! starting budget 
! initialize the flux arrays
do k = 1 , kz
  do j = jstart , jend
    do i = istart , iend 
      ztnew = zt(j,i,k)+dt*(zttendc(j,i,k)-ztentkeep(j,i,k))              ![ztnew]=K      
      if (k == 1) then                              !top of the atm
        zsumq0(j,i,k) = d_zero ! total water
        zsumh0(j,i,k) = d_zero ! liquid water temperature
      else
        zsumq0(j,i,k) = zsumq0(j,i,k-1)             
        zsumh0(j,i,k) = zsumh0(j,i,k-1)
      endif
    
      do n = 1 , nqx
          if (kphase(n) == 1) ztnew=ztnew-wlhvocp*(zqxx(j,i,k,n)+ &
                                    & (zqxtendc(j,i,k,n)-ztenkeep(j,i,k,n))*dt)
          if (kphase(n) == 2) ztnew=ztnew-wlhfocp*(zqxx(j,i,k,n)+ &
                                    & (zqxtendc(j,i,k,n)-ztenkeep(j,i,k,n))*dt)
          zsumq0(j,i,k)=zsumq0(j,i,k)+ (zqxx(j,i,k,n)+(zqxtendc(j,i,k,n)-ztenkeep(j,i,k,n))*dt)* &
                        & (papf(j,i,k+1)-papf(j,i,k))*regrav
      end do

      ztnew = ztnew-wlhvocp*ztmpl-wlhfocp*ztmpi  
      zsumq0(j,i,k) = zsumq0(j,i,k)+(ztmpl+ztmpi)*(papf(j,i,k+1)-papf(j,i,k))*regrav    !(kg/m^2)
           
      ! Detrained water treated here
      zqe = zqdetr(j,i,k)*dt*egrav/(papf(j,i,k+1)-papf(j,i,k))                  !adimensionale?
      if (zqe > rlmin) then
        zsumq0(j,i,k) = zsumq0(j,i,k)+zqdetr(j,i,k)*dt                         ![zqdetr]=kg/(m^2*s)
        zalfaw = phases(zt(j,i,k))
        ztnew = ztnew-(wlhvocp*zalfaw+wlhfocp*(d_one-zalfaw))*zqe
      endif
      zsumh0(j,i,k) = zsumh0(j,i,k)+(papf(j,i,k+1)-papf(j,i,k))*ztnew
    end do 
  end do
end do 

do k = 1 , kz
  do j = jstart , jend
    do i = istart , iend
      zsumh0(j,i,k) = zsumh0(j,i,k)/(papf(j,i,k+1)-papf(j,i,1))
    end do
  end do
end do

end if !budget

!----------------------------------------------------------------------
!                       START OF VERTICAL LOOP
!----------------------------------------------------------------------
!-----------------
! Loop over levels
do k = 1 , kz 
  ! Derived variables needed                                                  
  do j = jstart, jend                                                  
    do i = istart, iend                                                
      zdp(j,i) = papf(j,i,k+1)-papf(j,i,k)                        !dp                 
      zgdp(j,i) = egrav/zdp(j,i)                                  !g/dp  =(1/m)               
      zdtgdp(j,i) = dt*zgdp(j,i)                              !(dt*g)/dp =(dt/m)                  
      zrdtgdp(j,i) = zdp(j,i)*(d_one/(dt*egrav))              !dp/(gdt)=m/dt
      zqtmst = 1/dt                                           !1/dt    
    end do                                                             
  end do     
  ! Set the fall velocities
  do j = jstart , jend
    do i = istart , iend
      zvqx(j,i,k,iqqv) = d_zero !*sqrt(ZQX(JL,JK,IQV))
      zvqx(j,i,k,iqql) = 0.1D0  !*sqrt(ZQX(JL,JK,IQL))
      zvqx(j,i,k,iqqr) = 4.0D0  !*sqrt(ZQX(JL,JK,IQR))
      zvqx(j,i,k,iqqi) = 0.15D0 !*sqrt(ZQX(JL,JK,IQI))
      zvqx(j,i,k,iqqs) = 1.0D0  !*sqrt(ZQX(JL,JK,IQS))
    end do
  end do

  ! Reset matrix so missing pathways are set
  zsolqb(:,:,:,:)  = d_zero  !_JPRB
  zsolqa(:,:,:,:)  = d_zero  !_JPRB
  zfluxq(:,:,:)    = d_zero  !_JPRB
  zfallsrce(:,:,:) = d_zero  !_JPRB
  zfallsink(:,:,:) = d_zero  !_JPRB
  zconvsrce(:,:,:) = d_zero
  zconvsink(:,:,:) = d_zero
  zratio(:,:,:) = d_zero

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
      if ( zt(j,i,k) > tzero ) then
         satvp(j,i,k) = satw(zt(j,i,k),pres(j,i,k))
         dqsatdt(j,i,k) = dqsatdtw(zt(j,i,k),satvp(j,i,k))
      else
        satvp(j,i,k) = satc(zt(j,i,k),pres(j,i,k))
        satice(j,i,k) = satvp(j,i,k) 
        dqsatdt(j,i,k) = dqsatdtc(zt(j,i,k),satvp(j,i,k))
      end if
    end do
  end do
 
  ! Define the microphysics
  ! the matrix will be sparse is this a problem ?
  ! (X,Y) means a sink of X and a source of Y
  ! for the implementation I will use flexible pointers
  ! such that it will be written (IQR,IQG) to indicate graupel to rain
  ! and the parametrization can have different variables switched on
  ! and off.
  ! each of these is a parametrization for a microphysical process.
  ! lmicro = .false.
  lmicro = .true.
 
  if ( lmicro ) then
    ! Turn off microphysics
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                         CONDENSATION                       ! 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    rovcp = rgas/cpd
    do j = jstart , jend
      do i = istart , iend
        zdtdp = rovcp*zt(j,i,k)/pres(j,i,k)
        zcond = dqsatdt(j,i,k)*((omega(j,i,k)*zdtdp)+radheatrt(j,i,k))
        if (zcond < rlmin) zcond = d_zero
        zsolqa(j,i,iqqv,iqql) = zsolqa(j,i,iqqv,iqql) - zcond
        zsolqa(j,i,iqql,iqqv) = zsolqa(j,i,iqql,iqqv) + zcond
      end do
    end do 
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                         AUTOCONVERSION                     ! 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    !Warm clouds
    do j = jstart , jend
      do i = istart , iend
        ztmpa=d_one/max(fcc(j,i,k),zepsec) 
        zliqcld(j,i) = zqxx(j,i,k,iqql)*ztmpa
        zicecld(j,i) = zqxx(j,i,k,iqqi)*ztmpa
        if ( zliqcld(j,i) > zepsec ) then    !why if B=0 by default?
          zsolqb(j,i,iqql,iqqv) = d_zero 
          select case (kautoconv)
          case (1)          !Klein & Pincus (2000)
            zsolqb(j,i,iqqr,iqql) = zsolqb(j,i,iqqr,iqql) + &
                                    dt*zauto_rate_klepi * (zqxx(j,i,k,iqql)**(3.3D0))
            zsolqa(j,i,iqqr,iqql) = d_zero
          case (2)          ! Khairoutdinov and Kogan (2000)
            zsolqb(j,i,iqqr,iqql) = zsolqb(j,i,iqqr,iqql) + &
                                    dt*zauto_rate_khair *         &
                                   (zqxx(j,i,k,iqql)**(zauto_expon_khair))
          case (3)          !Kessler(1969)
            if ( zqxx(j,i,k,iqql) > zautocrit ) then
              zsolqa(j,i,iqqr,iqql) = zsolqa(j,i,iqqr,iqql) - zauto_rate_kesl*zautocrit
              zsolqa(j,i,iqql,iqqr) = zsolqa(j,i,iqql,iqqr) + zauto_rate_kesl*zautocrit
              zsolqb(j,i,iqqr,iqql) = zsolqb(j,i,iqqr,iqql) + dt*zauto_rate_kesl
            end if
          case (4)           !Sundqvist
            zsolqb(j,i,iqqr,iqql) = zsolqb(j,i,iqqr,iqql) + & 
                                    dt*zauto_rate_sundq* &
                                   (d_one-dexp(-(zqxx(j,i,k,iqql)/zautocrit)**2))
            zsolqa(j,i,iqqr,iqql) = zsolqa(j,i,iqqr,iqql) + d_zero
          end select
        end if !(ZQX(JL,JK,IQL)>0.0)
       
        ! Cold clouds
        ! Snow Autoconversion rate follow Lin et al. 1983
        if (zt(j,i,k) <=  tzero .and. zicecld(j,i) > zepsec ) then
          alpha1 = dt*1.0D-3*exp(0.025*(zt(j,i,k)-tzero))
          zsolqa(j,i,iqqs,iqqi)=zsolqa(j,i,iqqs,iqqi)-alpha1*qi0
          zsolqa(j,i,iqqi,iqqs)=zsolqa(j,i,iqqi,iqqs)+alpha1*qi0
          zsolqb(j,i,iqqs,iqqi)=zsolqb(j,i,iqqs,iqqi)+dt*alpha1
        end if
       
      
       ! if (zliqcld(j,i) > zepsec) then
       !   if (zt(j,i,k) <= tzero) then
       !     zsolqb
       ! end if    
   
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !                         EVAPORATION                        ! 
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!        
        zsolqb(j,i,iqqv,iqqr) = zsolqb(j,i,iqqv,iqqr)+&
                                & dt*kevap*max((satvp(j,i,k)-zqxx(j,i,k,iqqv)),d_zero)
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
        zfrzmax(j,i) = max((tzero-zt(j,i,k))*zrldcp,d_zero)
      end do
    end do
          
    do j = jstart , jend
      do i = istart , iend
        if (zfrzmax(j,i) > zepsec .AND. zqxx(j,i,k,iqqr) > zepsec) then
         zfrz = min(zqxx(j,i,k,iqqr),zfrzmax(j,i))
          zsolqa(j,i,iqqs,iqqr) = zsolqa(j,i,iqqs,iqqr) + zfrz
          zsolqa(j,i,iqqr,iqqs) = zsolqa(j,i,iqqr,iqqs) - zfrz
        end if 
      end do
    end do
          
    ! Freezing of liquid
    do j = jstart , jend
      do i = istart , iend
        if (zfrzmax(j,i) > zepsec .AND. zqxx(j,i,k,iqql) > zepsec) then
          zfrz = min(zqxx(j,i,k,iqql),zfrzmax(j,i))
          zsolqa(j,i,iqqi,iqql) = zsolqa(j,i,iqqi,iqql) + zfrz
          zsolqa(j,i,iqql,iqqi) = zsolqa(j,i,iqql,iqqi) - zfrz
        end if
      end do
    end do

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                         MELTING                                  !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! The melting of ice and snow are treated explicitly.
    ! First water and ice saturation are found
    !---------------------------------------------
    ! ice saturation T<273K
    ! liquid water saturation for T>273K
    !---------------------------------------------

    do j = jstart , jend
      do i = istart , iend 
        zphases = phases(zt(j,i,k))    
        zfoeew(j,i,k) = min(foeewm(zt(j,i,k))/pres(j,i,k),d_half)
        zfoeew(j,i,k) = min(d_half,zfoeew(j,i,k))
        zqsice(j,i,k) = zfoeew(j,i,k)/(d_one-retv*zfoeew(j,i,k))
        !----------------------------------
        ! liquid water saturation
        !----------------------------------
        zfoeeliq(j,i,k) = min(foeeliq(zt(j,i,k))/pres(j,i,k),d_half)
        zqsliq(j,i,k) = zfoeeliq(j,i,k)
        zqsliq(j,i,k) = zqsliq(j,i,k)/(d_one-retv*zqsliq(j,i,k))
      end do
    end do
        
    ! MELTING OF SNOW and ICE
    do j = jstart, jend
      do i = istart, iend
        zicetot(j,i) = zqxx(j,i,k,iqqi)+zqxx(j,i,k,iqqs)
        zmeltmax(j,i) = d_zero
        if (zicetot(j,i) > zepsec .and. zt(j,i,k) > tzero) then
          !Calculate subsaturation 
          zsubsat = max(zqsice(j,i,k)-zqxx(j,i,k,iqqv),d_zero)
          ! Calculate difference between dry-bulb (zt)  and the temperature
          ! at which the wet-bulb=0degC  
          ! Melting only occurs if the wet-bulb temperature >0
          ! i.e. warming of ice particle due to melting > cooling
          ! due to evaporation.     
          ztdiff = zt(j,i,k)-tzero-zsubsat*&
                   &(ztw1+ztw2*(pres(j,i,k)-ztw3)-ztw4*(zt(j,i,k)-ztw5))
          ! Ensure ZCONS1 is positive so that ZMELTMAX=0 if ZTDMTW0<0
          zcons1 = abs(dt*(d_one+d_half*ztdiff)/rtaumel)  
          zmeltmax(j,i) = max(ztdiff*zcons1*zrldcp,d_zero)
        end if
      end do
    end do    

    ! Loop over frozen hydrometeors (kphase==2 (ice, snow))
    do n = 1, nqx
      if (kphase(n) == 2) then
        m = imelt(n) 
        do j = jstart, jend
          do i = istart, iend
            if (zmeltmax(j,i) > zepsec .and. zicetot(j,i) > d_zero) then !zepsec
              zphases = zqxx(j,i,k,n)/zicetot(j,i)
              zmelt = min(d_one,zphases*zmeltmax(j,i))
              zsolqa(j,i,m,n) = zsolqa(j,i,m,n)+zmelt
              zsolqa(j,i,n,m) = zsolqa(j,i,n,m)-zmelt
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
 

  end if!lmicro      

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !                 DETRAINMENT FROM CONVECTION
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  if (k < kz .and. k >= 1) then
    do j = jstart , jend
      do i = istart , iend
        zalfaw = zliq(j,i,k)
        zice = d_one-zalfaw
        zqdetr(j,i,k) = qdetr(j,i,k)*zdtgdp(j,i)
        if (zqdetr(j,i,k) > rlmin) then
          zconvsrce(j,i,iqql) = zalfaw*zqdetr(j,i,k)
          zconvsrce(j,i,iqqi) = zice*zqdetr(j,i,k)
          zsolqa(j,i,iqql,iqql) = zsolqa(j,i,iqql,iqql)+zconvsrce(j,i,iqql)
          zsolqa(j,i,iqqi,iqqi) = zsolqa(j,i,iqqi,iqqi)+zconvsrce(j,i,iqqi)
        else
          zqdetr(j,i,k) = d_zero
        end if 
      end do
    end do 
  end if
  !!!!!!!!!!!!! sedimentation/falling of microphysical species
  !     now that rain, snow, graupel species are prognostic
  !     the precipitation flux can be defined directly level by level

  do n = 1 , nqx
    do j = jstart , jend
      do i = istart , iend
        ! Source from layer above 
        if ( k > 1 ) then
          zfallsrce(j,i,n) = zpfplsx(j,i,k,n)*zdtgdp(j,i)
          zsolqa(j,i,n,n) = zsolqa(j,i,n,n)+zfallsrce(j,i,n)
        end if
        ! Sink to next layer, constant fall speed
        zfall = zvqx(j,i,k,n)*rhob3d(j,i,k)
        zfallsink(j,i,n) = zdtgdp(j,i)*zfall
      end do
    end do
  end do
      
  !----------------------------------------------------------
  ! Truncate sum of explicit sinks to size of bin
  ! this approach is inaccurate, but conserves -
  ! prob best can do with explicit (i.e. not implicit!) terms
  !----------------------------------------------------------
  jindex1(:,:,:,:) = 0
  zsinksum(:,:,:) = d_zero

   do n = 1 , nqx
    do jn = 1 , nqx
      do j = jstart , jend
        do i = istart , iend
          if ( zsolqa(j,i,jn,n) < d_zero ) then                           !se A<0 jindex=1 e zsinksum = -A
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
        if ( zsinksum(j,i,n) > d_zero ) zratio(j,i,n) = zqxx(j,i,k,n)/zsinksum(j,i,n)
          zratio(j,i,n) = max(min(zratio(j,i,n),d_one),d_zero)
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
  ! Set the LHS of equation
  do n = 1 , nqx
    do jn = 1 , nqx
      do j = jstart , jend
        do i = istart , iend
          ! Diagonals: microphysical sink terms+transport
          if ( jn == n ) then
            zqlhs(j,i,jn,n) = d_one + zfallsink(j,i,n)
            do jo = 1 , nqx
              zqlhs(j,i,jn,n) = zqlhs(j,i,jn,n) + zsolqb(j,i,jo,jn)
            end do
            ! Non-diagonals: microphysical source terms
          else
            ! Here is the delta T - missing from doc.
            zqlhs(j,i,jn,n) = -zsolqb(j,i,jn,n)
          end if
        end do
      end do
    end do
  end do

  ! Set the RHS of equation
  do n = 1 , nqx
    do j = jstart , jend
      do i = istart , iend
        !   Sum the explicit source and sink
        zexplicit = d_zero
        do jn = 1 , nqx
          ! Positive, since summed over 2nd index
          zexplicit = zexplicit + zsolqa(j,i,n,jn)
        end do
        zqxn(j,i,k,n) = zqxx(j,i,k,n) + zexplicit
      end do
    end do
  end do

  ! Solve by LU decomposition:
  ! dummy test
  ! ZQLHS(:,:,:)=0.0
  ! DO JM=1,NQX
  ! ZQLHS(1,JM,JM)=1.0
  ! ENDDO
  ! ZQLHS(1,3,3)=ZQLHS(1,3,3)+0.47
  ! ZQLHS(1,5,3)=ZQLHS(1,5,3)-0.47

  call ludcmp(jstart,jend,istart,iend,zqlhs,jindex2)
  call lubksb(jstart,jend,istart,iend,zqlhs,jindex2,zqxn)
  !------------------------------------------------------------------------
  !  Precipitation/sedimentation fluxes to next level
  !     diagnostic precipitation fluxes
  !     It is this scaled flux that must be used for source to next layer
  !------------------------------------------------------------------------
  ! Generalized precipitation flux
  do n = 1 , nqx
    do j = jstart , jend
      do i = istart , iend
        zpfplsx(j,i,k+1,n) = zfallsink(j,i,n)*zqxn(j,i,k,n)*zrdtgdp(j,i) !this will be the source for the k 
      end do                                                             !kg/m2/s 
    end do
  end do

  ! Calculate fluxes in and out of box for conservation of TL
  do n = 1 , nqx
    do j = jstart , jend
      do i = istart , iend
        zfluxq(j,i,n)=zconvsrce(j,i,n)+ zfallsrce(j,i,n) -&
                    &(zfallsink(j,i,n)+zconvsink(j,i,n))*zqxn(j,i,k,n)
      end do
    end do
  end do

  ! Calculate the water variables tendencies
  do n = 1 , nqx
    do j = jstart , jend
      do i = istart , iend
        zqxtendc(j,i,k,n) = zqxtendc(j,i,k,n) + (zqxn(j,i,k,n)-zqxx(j,i,k,n))*zqtmst
      end do
    end do 
  end do
     
  ! Calculate the temperature tendencies
  do n = 1 , nqx
    do j = jstart , jend
      do i = istart , iend
        if ( kphase(n) == 1 ) then
          zttendc(j,i,k) = zttendc(j,i,k) + wlhvocp*(zqxn(j,i,k,n)-zqxx(j,i,k,n)-zfluxq(j,i,n))*zqtmst
        end if
        if ( kphase(n) == 2 ) then
          zttendc(j,i,k) = zttendc(j,i,k) + wlhfocp*(zqxn(j,i,k,n)-zqxx(j,i,k,n)-zfluxq(j,i,n))*zqtmst
        end if
      end do
    end do 
  end do

  ! Couple tendencies with pressure
  do n = 1 , nqx
    do j = jstart , jend
      do i = istart , iend
        zqxten(j,i,k,n) = zqxtendc(j,i,k,n)*psf(j,i) 
      end do
    end do
  end do

  do j = jstart , jend
    do i = istart , iend
      ztten(j,i,k) = zttendc(j,i,k)*psf(j,i)
    end do
  end do
 
end do   ! kz : end of vertical loop

!-------------------------------------
! Final enthalpy and total water diagnostics 
!-------------------------------------

if (budget) then
  ! Initialize the flux arrays
  do k = 1 , kz
    do j = jstart , jend
      do i = istart , iend 
        ztnew = zt(j,i,k)+dt*(zttendc(j,i,k)-ztentkeep(j,i,k))      
        if (k == 1) then
          zsumq1(j,i,k) = d_zero ! total water
          zsumh1(j,i,k) = d_zero ! liquid water temperature
        else
          zsumq1(j,i,k) = zsumq1(j,i,k-1)
          zsumh1(j,i,k) = zsumh1(j,i,k-1)
        endif

        ! cld vars
        do n = 1 , nqx
          if (kphase(n) == 1) ztnew = ztnew-wlhvocp*(zqxx(j,i,k,n)+ &
                                      & (zqxtendc(j,i,k,n)-ztenkeep(j,i,k,n))*dt)
          if (kphase(n) == 2) ztnew = ztnew-wlhfocp*(zqxx(j,i,k,n)+ &
                                      & (zqxtendc(j,i,k,n)-ztenkeep(j,i,k,n))*dt)
          zsumq1(j,i,k) = zsumq1(j,i,k)+ (zqxx(j,i,k,n)+(zqxtendc(j,i,k,n)-ztenkeep(j,i,k,n))*dt)* &
                      & (papf(j,i,k+1)-papf(j,i,k))*regrav
        end do
      
        zsumh1(j,i,k) = zsumh1(j,i,k)+(papf(j,i,k+1)-papf(j,i,k))*ztnew
      
        zrain = d_zero         
    
        do n = 1,nqx
          zrain = zrain+dt*zpfplsx(j,i,k+1,n) 
        end do
        zerrorq(j,i,k) = zsumq1(j,i,k)+zrain-zsumq0(j,i,k)
      end do
    end do
  end do

  do k = 1 , kz
    do i = istart, iend
      do j = jstart , jend
        zdtgdp(j,i) = dt*egrav/(papf(j,i,k+1)-papf(j,i,k))
        zrain = d_zero
        do n = 1 , nqx
          if (kphase(n) == 1) zrain = zrain+wlhvocp*zdtgdp(j,i)*zpfplsx(j,i,k+1,n)* &                    !k+1?
                                      & (papf(j,i,k+1)-papf(j,i,k))
          if (kphase(n) == 2) zrain = zrain+wlhfocp*zdtgdp(j,i)*zpfplsx(j,i,k+1,n)* &
                                      & (papf(j,i,k+1)-papf(j,i,k))
        end do
        zsumh1(j,i,k) = (zsumh1(j,i,k)-zrain)/(papf(j,i,k+1)-papf(j,i,1))
        zerrorh(j,i,k) = zsumh1(j,i,k)-zsumh0(j,i,k)
      enddo
    enddo
  end do

  do k = 1 , kz
    do i = istart, iend
      do j = jstart, jend
        do n=1,nqx
          if (abs(zerrorq(j,i,kz))>1.D-12.OR.abs(zerrorh(j,i,kz))>1.D-12) then
            print*,'not conserved' 
          end if
        end do
      end do
    end do
  end do
end if!budget
 
! Sum fluxes over the levels
! Initialize fluxes
prcflxw(:,:) = d_zero
prcflxc(:,:) = d_zero
zpfplsl(:,:,:) = d_zero
zpfplsn(:,:,:) = d_zero

!--------------------------------------------------------------------
! Copy general precip arrays back into PFP arrays for GRIB archiving
! Add rain and liquid fluxes, ice and snow fluxes
!--------------------------------------------------------------------
!Rain+liquid, snow+ice
do k = 1 , kz+1
  do j = jstart , jend
    do i = istart , iend
      do n = 1 , nqx
        if ( kphase(n) == 1 ) then
          zpfplsl(j,i,k) = zpfplsl(j,i,k) + zpfplsx(j,i,k,n)
        end if
        if ( kphase(n) == 2 ) then
          zpfplsn(j,i,k) = zpfplsn(j,i,k) + zpfplsx(j,i,k,n)
        end if
      end do
    end do
  end do
end do

!--------------------------------------------------------------------
!Convert the accumlated precipitation to appropriate units for
!the surface physics and the output
!--------------------------------------------------------------------

do j = jstart , jend
  do i = istart , iend
    do k = 1 , kz
      prainx = zpfplsl(j,i,k+1)*dt
      if ( prainx > dlowval ) then
        rainnc(j,i) =  rainnc(j,i) + prainx   !mm
        lsmrnc(j,i) =  lsmrnc(j,i) + zpfplsl(j,i,k+1)*rtsrf
      end if
     ! rainnc(j,i) = prcflxc(j,i) + zpfplsn(j,i,k)   !mm/s
      
   !   print*,'cold', prcflxc(j,i)
    end do
  end do
end do

#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif

contains

     real(rk8) function phases(zt)
       implicit none
       real(rk8) , intent(in):: zt
       real(rk8) , parameter :: rtice =  250.16D0               !tzero - 23.0
       real(rk8) , parameter :: rtwat_rtice_r=d_one/23.0D0      !1.0_JPRB/(RTWAT-RTICE)
       phases = min(d_one,((max(rtice,min(tzero,zt))-rtice)*rtwat_rtice_r)**2)
     end function phases

     real(rk8) function foeewm(zt)
       implicit none
       real(rk8) , intent(in):: zt
       real(rk8) , parameter :: r2es =  611.21D0*rgow
       real(rk8) , parameter :: r3les = 17.502D0
       real(rk8) , parameter :: r3ies = 22.587D0
       real(rk8) , parameter :: r4les = 32.19D0
       real(rk8) , parameter :: r4ies = -0.7D0
       foeewm = r2es*(phases(zt)*exp(r3les*(zt-tzero)/(zt-r4les))+&
       &(d_one-phases(zt))*exp(r3ies*(zt-tzero)/(zt-r4ies)))
       end function foeewm

     real(rk8) function foeeliq(zt)
       implicit none
       real(rk8) , intent(in):: zt
      ! real(rk8) , parameter :: rtice =  250.16D0 
      ! real(rk8) , parameter :: rtwat_rtice_r=d_one/23.0D0
       real(rk8) , parameter :: r2es =  611.21D0*rgow
       real(rk8) , parameter :: r3les = 17.502D0
       real(rk8) , parameter :: r4les = 32.19D0
       foeeliq = r2es*exp(r3les*(zt-tzero)/(zt-r4les))
     end function foeeliq

     real(rk8) function foeeice(zt)
       implicit none
       real(rk8) , intent(in):: zt
      ! real(rk8) , parameter :: rtice =  250.16D0
       real(rk8) , parameter :: r3ies = 22.587D0
       real(rk8) , parameter :: r2es =  611.21D0*rgow
       real(rk8) , parameter :: r4ies = -0.7D0
       real(rk8) , parameter :: rtwat_rtice_r=d_one/23.0D0
        foeeice = r2es*exp(r3ies*(zt-tzero)/(zt-r4ies))
     end function foeeice

     real(rk8) function satc(zt,xp)
       implicit none
       real(rk8) , intent(in) :: xp , zt
       ! saturation mixing ratio for ice in Pa (A Description of the Fifth-Generation Penn State/NCAR
       !Mesoscale Model (MM5)Georg A. Grell, Jimy Dudhia, David R. Stauffer, NCAR TECHNICAL NOTE Dec 1994)
       satc = (3.79D2/xp)*dexp(22.514D0-6150.0D0/zt) 
     end function satc

     real(rk8) function satw(zt,xp)
       implicit none
       real(rk8) , intent(in) :: xp , zt
       !  ! saturation mixing ratio for water in Pa (A Description of the Fifth-Generation Penn State/NCAR
       !Mesoscale Model (MM5)Georg A. Grell, Jimy Dudhia, David R. Stauffer, NCAR TECHNICAL NOTE Dec 1994)
       satw = (3.8014D2/xp)*dexp(17.67D0*(zt-tzero)/(zt-29.65D0))
     end function satw
     
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

!#ifdef DEBUG
!  call time_end(subroutine_name,idindx)
!#endif

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
    ! AND INDX ARE NOT MODIFIED BY THIS ROUTINE AND CAN BE
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

! solves A x = b-------------> LU x = b---------> Ly=b
!                                                 Ux=y
!

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
