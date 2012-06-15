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
!
  use mod_runparams
  use mod_mpmessage
  use mod_memutil
  use mod_atm_interface , only : atmstate , slice , surfstate
  use mod_mppparam , only : iqqv => iqv
  use mod_mppparam , only : iqqc => iqc

  use mod_pbl_common

 

  private

  integer , parameter :: nmax = 5
  

  real(dp) , parameter :: ptsphy = 1000.0
  real(dp) , parameter :: lvl = 2.50E6
  real(dp) , parameter :: lvs = 2.83E6
  real(dp) , parameter :: dpotdz = -1.0E-3
  real(dp) , parameter :: vervel = 0.0015
  real(dp) , parameter :: qrad = 1.0E-2
  real(dp) :: w2 
  real(dp) :: w1
 
  real(dp) , pointer , dimension(:,:,:) :: pres   ! from atms
  real(dp) , pointer , dimension(:,:,:) :: zt     ! from atms
  real(dp) , pointer , dimension(:,:,:) :: zeta   ! from atms
  real(dp) , pointer , dimension(:,:,:) :: dzeta  ! from atms
  real(dp) , pointer , dimension(:,:,:,:) :: zqxx ! from atms
  real(dp) , pointer , dimension(:,:,:) :: rhob3d ! from atms added by R
  real(dp) , pointer , dimension(:,:,:) :: omega  ! from atms added by R
  real(dp) , pointer , dimension(:,:,:) :: heatrt ! added by R


  integer , parameter :: nqx = 5

  public :: allocate_mod_cloud_s1 , init_cloud_s1 , microphys

  real(dp) , pointer , dimension(:,:,:) :: ztl
  real(dp) , pointer , dimension(:,:,:) :: dqsatdt
  real(dp) , pointer , dimension(:,:,:) :: satvp
  real(dp) , pointer , dimension(:,:,:,:) :: wvflux
  real(dp) , pointer , dimension(:,:,:,:) :: vqx
  real(dp) , pointer , dimension(:,:,:,:) :: zqx
  real(dp) , pointer , dimension(:,:,:) :: zqxfg
  real(dp) , public  , pointer, dimension(:,:,:,:) :: zqxn
  real(dp) , pointer , dimension(:,:,:) :: zfallsink
  real(dp) , pointer , dimension(:,:,:) :: zfallsrce
  real(dp) , pointer , dimension(:,:,:) :: zfluxq
  real(dp) , pointer , dimension(:,:,:) :: zratio
  real(dp) , pointer , dimension(:,:,:) :: zsinksum
  real(dp) , pointer , dimension(:,:,:,:) :: zqlhs
  real(dp) , pointer , dimension(:,:,:,:) :: zsolqa
  real(dp) , pointer , dimension(:,:,:,:) :: zsolqb
  real(dp) , pointer , dimension(:,:,:) :: totalw1
  real(dp) , pointer , dimension(:,:,:) :: totalw2
  real(dp) , pointer , dimension(:,:,:) :: diff
  real(dp) , pointer , dimension(:,:) :: zfrzmax
  integer , pointer , dimension(:) :: kphase                     ! marker for water phase of each species
                                                                 ! 0=vapour, 1=liquid, 2=ice  
  integer , pointer , dimension(:) :: imelt                      ! marks melting linkage for ice categories
                                                                 ! ice->liquid, snow->rain
  integer , pointer , dimension(:,:,:,:) :: jindex1
  integer , pointer , dimension(:,:,:) :: jindex2
   


  contains

  subroutine allocate_mod_cloud_s1
    implicit none
    call getmem3d(ztl,jci1,jci2,ici1,ici2,1,kz,'clouds1:ztl')
    call getmem3d(dqsatdt,jci1,jci2,ici1,ici2,1,kz,'clouds1:dqsatdt')
    call getmem3d(satvp,jci1,jci2,ici1,ici2,1,kz,'clouds1:satvp')
    call getmem4d(vqx,jci1,jci2,ici1,ici2,1,kz,1,nqx,'clouds1:vqx')
    call getmem4d(zqx,jci1,jci2,ici1,ici2,1,kz,1,nqx,'clouds1:zqx')
    call getmem3d(zqxfg,jci1,jci2,ici1,ici2,1,nqx,'clouds1:zqxfg')
    call getmem4d(zqxn,jce1,jce2,ice1,ice2,1,kz,1,nqx,'clouds1:zqxn')
    call getmem2d(zfrzmax,jci1,jci2,ici1,ici2,'clouds1:zfrzmax')
    call getmem3d(zfallsink,jci1,jci2,ici1,ici2,1,nqx,'clouds1:zfallsink')
    call getmem3d(zfallsrce,jci1,jci2,ici1,ici2,1,nqx,'clouds1:zfallsrce')
    call getmem3d(zfluxq,jci1,jci2,ici1,ici2,1,nqx,'clouds1:zfluxq')
    call getmem3d(zratio,jci1,jci2,ici1,ici2,1,nqx,'clouds1:zratio')
    call getmem3d(zsinksum,jci1,jci2,ici1,ici2,1,nqx,'clouds1:zsinksum')
    call getmem4d(zqlhs,jci1,jci2,ici1,ici2,1,nqx,1,nqx,'clouds1:zqlhs')
    call getmem4d(zsolqa,jci1,jci2,ici1,ici2,1,nqx,1,nqx,'clouds1:zsolqa')
    call getmem4d(zsolqb,jci1,jci2,ici1,ici2,1,nqx,1,nqx,'clouds1:zsolqb')
    call getmem4d(jindex1,jci1,jci2,ici1,ici2,1,nqx,1,nqx,'clouds1:jindex1')
    call getmem3d(jindex2,jci1,jci2,ici1,ici2,1,nqx,'clouds1:jindex2')
    call getmem1d(kphase,1,nqx,'clouds1:kphase')
    call getmem1d(imelt,1,nqx,'clouds1:imelt')
    call getmem3d(totalw1,jci1,jci2,ici1,ici2,1,kz,'clouds1:totalw1')
    call getmem3d(totalw2,jci1,jci2,ici1,ici2,1,kz,'clouds1:totalw2')
    call getmem3d(diff,jci1,jci2,ici1,ici2,1,kz,'clouds1:diff')
    call getmem4d(wvflux,jci1,jci2,ici1,ici2,1,kz,1,nqx,'clouds1:wvflux')
    !call getmem3d(omega,jce1,jce2,ice1,ice2,1,kz,'atm_interface:omega')

    



  end subroutine allocate_mod_cloud_s1

  subroutine init_cloud_s1(atms)
    implicit none
    type(slice) , intent(in) :: atms
    call assignpnt(atms%pb3d,pres)
    call assignpnt(atms%tb3d,zt)
    call assignpnt(atms%za,zeta)
    call assignpnt(atms%dzq,dzeta)
    call assignpnt(atms%qxb3d,zqxx)
    call assignpnt(atms%rhob3d,rhob3d)
    call assignpnt(heatrt,radheatrt)
  end subroutine init_cloud_s1

  subroutine microphys(omega,jstart,jend,istart,iend)
    implicit none
    integer , intent(in) :: jstart , jend , istart , iend
    real(dp) , pointer , dimension(:,:,:) , intent(in) :: omega  !added by R

    integer :: i , j , k , n , m
    integer :: iqi , iql , iqr , iqs , iqv , jn , jo , kautoconv
    logical :: lmicro
    real(dp) :: zcond , zdtdp , zexplicit
    real(dp) :: zfrz
    real(dp) :: alpha1 !coefficient autoconversion
    real(dp) :: slht !sublimation latent heat 
    real(dp) :: ralsdcp
    real(dp) :: ralvdcp
    real(dp) :: zrldcp
!
    ! local real variables for autoconversion rate constants
    real(dp) , parameter :: zauto_rate_khair = 0.355D0  ! microphysical terms
    real(dp) , parameter :: zauto_expon_khair = 1.47D0
    real(dp) , parameter :: zauto_rate_sundq = 0.5D-3
    real(dp) , parameter :: zauto_rate_kesl = 1.D-3         !giusto!
    real(dp) , parameter :: zauto_rate_klepi = 0.5D-3
    real(dp) , parameter :: zautocrit = 5.D-4               !giusto!
    real(dp) , parameter :: zt0 = 273.16
    real(dp) , parameter :: zepsec = 1.E-10
    real(dp) , parameter :: qi0 = 1.0E-3 !g g^(-1)
    real(dp) , parameter :: vlht = 2260. !kJ/kg evaporation latent heat  
    real(dp) , parameter :: rcpd =  1005. !J/(kg·K) Cp_dry (dry air calorific capacity at constant pressure) 
 

    ! local real constants for evaporation
    real(dp) , parameter :: kevap = 0.100D-02    !! Raindrop evap rate coef
    real(dp) , parameter :: rlmin = 1.D-8
 
     
    ! Define species tags, hard wire for now
    iqv = 1    ! vapour
    iql = 2    ! liquid cloud water
    iqr = 3    ! rain water
    iqi = 4    ! ice
    iqs = 5    ! snow
 
    ! Define which autoconversion paramaterization to be chosen
    !KAUTOCONV=1 ! Klein & Pincus (2000)
    !kautoconv = 2 ! Khairoutdinov and Kogan (2000)
    KAUTOCONV=3 ! Kessler (1969)
    !KAUTOCONV=4 ! Sundqvist
 
    ! Define species phase, 0=vapour, 1=liquid, 2=ice
    kphase(iqv) = 0
    kphase(iql) = 1
    kphase(iqr) = 1
    kphase(iqi) = 2
    kphase(iqs) = 2


    ! Set up melting/freezing index,
    ! if an ice category melts/freezes, where does it go?

    imelt(iqv)=-99
    imelt(iql)=iqi
    imelt(iqr)=iqs
    imelt(iqi)=iql
    imelt(iqs)=iqr

   
    
   
    !zqxn(jce1:jce2,ice1:ice2,1:kz,1:nqx)= d_zero
    zqx(jci1:jci2,ici1:ici2,1:kz,iqv) = zqxx(jci1:jci2,ici1:ici2,1:kz,iqqv)
    zqx(jci1:jci2,ici1:ici2,1:kz,iql) = zqxx(jci1:jci2,ici1:ici2,1:kz,iqqc)
    ztl(jci1:jci2,ici1:ici2,1:kz)     = zt(jci1:jci2,ici1:ici2,1:kz)
  
    do n = 1 , nqx
      do k = 1 , kz
    !print*, omega(2,2,k), 'level', k
        do j = jstart , jend
          do i = istart , iend
            if ( kphase(n) == 1 ) then
              ztl(j,i,k) = ztl(j,i,k) - lvl*zqx(j,i,k,n)*rcpd
            end if
            if ( kphase(n) == 2 ) then
              ztl(j,i,k) = ztl(j,i,k) - lvs*zqx(j,i,k,n)*rcpd
            end if
          end do
        end do
      end do
    end do
    
    !-----------------
    ! loop over levels
    do k = 1 , kz - 1
    !-----------------

    !---------------------------------
    ! First guess microphysics
    !---------------------------------
    do n = 1 , nqx
      do j = jstart , jend
        do i = istart , iend
          zqxfg(j,i,n)=zqx(j,i,k,n)
        end do
      end do
    end do


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!set the fall velocities!!!!!!!!!!
      do j = jstart , jend
        do i = istart , iend
          vqx(j,i,k,iqv) = d_zero !*sqrt(ZQX(JL,JK,IQV))
          vqx(j,i,k,iql) = 0.1D0  !*sqrt(ZQX(JL,JK,IQL))
          vqx(j,i,k,iqr) = 2.0D0  !*sqrt(ZQX(JL,JK,IQR))
          vqx(j,i,k,iqi) = d_zero !*sqrt(ZQX(JL,JK,IQI))
          vqx(j,i,k,iqs) = d_zero !*sqrt(ZQX(JL,JK,IQS))
        end do
      end do
      




      ! calculate cloud total water at the beginning

      do j = jstart , jend
        do i = istart , iend
          totalw1(j,i,k)= zqx(j,i,k,iqv)+zqx(j,i,k,iql)+zqx(j,i,k,iqr) + &
                          zqx(j,i,k,iqi)+zqx(j,i,k,iqs)
        end do
      end do
      w1=sum(totalw1)

        ! reset matrix so missing pathways are set
        zsolqb(:,:,:,:)  = d_zero  !_JPRB
        zsolqa(:,:,:,:)  = d_zero  !_JPRB
        zfluxq(:,:,:)    = d_zero  !_JPRB
        zfallsrce(:,:,:) = d_zero  !_JPRB
        zfallsink(:,:,:) = d_zero  !_JPRB
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
        ! detrainment into explicit source/sink array diagnognal
        ! ZSOLQA(IQL,IQL)=PLUDE
        !--------------------------------------------------------
        !!!Calculate the saturation mixing ratio
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
 
      !!!!!!!!!!!!!!!!!!CONDENSATION!!!!!!!!!!!!!!!!!!!!
 
      do j = jstart , jend
        ! these are explicit sources and sinks
        ! for example from condensation processes
        do i = istart , iend
          zdtdp = rovcp*zt(j,i,k)/pres(j,i,k)
          zcond = dqsatdt(j,i,k)*((omega(j,i,k)*zdtdp)+radheatrt(j,i,k))
          zsolqa(j,i,iqv,iql) = zsolqa(j,i,iqv,iql) - zcond
          zsolqa(j,i,iql,iqv) = zsolqa(j,i,iql,iqv) + zcond
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
        !!!!!!!!!!!!!!!!!!!AUTOCONVERSION!!!!!!!!!!!!!!!!!
        do j = jstart , jend
          do i = istart , iend
            if ( zqx(j,i,k,iql) > d_zero ) then
              ! se liquid > 0 allora B_lv=0, no source
              zsolqb(j,i,iql,iqv) = d_zero ! of liquid from vapor
              select case (kautoconv)
                case (1)          !Klein & Pincus (2000)
                  zsolqb(j,i,iqr,iql) = zsolqb(j,i,iqr,iql) + &
                        ptsphy*zauto_rate_klepi * (zqx(j,i,k,iql)**(3.3D0))
                  zsolqa(j,i,iqr,iql) = d_zero
                case (2)          ! Khairoutdinov and Kogan (2000)
                  zsolqb(j,i,iqr,iql) = zsolqb(j,i,iqr,iql) + &
                        ptsphy*zauto_rate_khair *         &
                        (zqx(j,i,k,iql)**(zauto_expon_khair))
                case (3)          !Kessler(1969)
                  if ( zqx(j,i,k,iql) > zautocrit ) then
                    zsolqa(j,i,iqr,iql) = zsolqa(j,i,iqr,iql) - &
                                       zauto_rate_kesl*zautocrit
                    zsolqa(j,i,iql,iqr) = zsolqa(j,i,iql,iqr) + &
                                       zauto_rate_kesl*zautocrit
                    zsolqb(j,i,iqr,iql) = zsolqb(j,i,iqr,iql) + &
                                       ptsphy*zauto_rate_kesl
                  else
                    zsolqb(j,i,iqr,iql) = d_zero
                    zsolqa(j,i,iql,iqr) = d_zero
                  end if
                case (4)           !Sundqvist
                  zsolqb(j,i,iqr,iql) = zsolqb(j,i,iqr,iql) + & 
                           ptsphy*zauto_rate_sundq* &
                          (d_one-dexp(-(zqx(j,i,k,iql)/zautocrit)**d_two))
                  zsolqa(j,i,iqr,iql) = zsolqa(j,i,iqr,iql) + d_zero
              end select
            end if !(ZQX(JL,JK,IQL)>0.0)
 
          !!!!!!!!!!!!!!!!!!!!EVAPORATION (R-->V) !!!!!!!!!!!!!!!!!!!!!
 
            zsolqb(j,i,iqv,iqr) = zsolqb(j,i,iqv,iqr) + &
                      ptsphy*kevap*max((satvp(j,i,k)-zqx(j,i,k,iqv)),d_zero)
            ! IF (ZQX(JL,JK,IQR)>0.0.AND.ZQX(JL,JK,IQV)<SATVP(JL,JK)) THEN
            !         ZSOLQB(JL,IQV,IQR)=KEVAP*(SATVP(JL,JK)-ZQX(JL,JK,IQV))
            ! ELSE
            !        ZSOLQB(JL,IQV,IQR)=0.0
            ! END IF

          !  Snow Autoconversion rate follow Lin et al. 1983
          print*,'temperatura', zt(j,i,k)
            if (zt(j,i,k) <=  zt0) then
              if (zqx(j,i,k,iqi) > zepsec) then
                alpha1 = ptsphy*1.0E-3*exp(0.025*(zt(j,i,k)-zt0))
                zsolqa(j,i,iqs,iqi)=zsolqa(j,i,iqs,iqi)-alpha1*qi0
                zsolqa(j,i,iqi,iqs)=zsolqa(j,i,iqi,iqs)+alpha1*qi0
                zsolqb(j,i,iqs,iqi)=zsolqb(j,i,iqs,iqi)+ptsphy*alpha1
              print*,'ice', zsolqa(j,i,iqi,iqs)
              end if
            end if
         end do
        end do

        !Freezing of rain

        do j = jstart , jend
          do i = istart , iend
             slht=(2834.1 -0.29*zt(j,i,k)-0.004*(zt(j,i,k)**2)) !sublimation latent heat
             !print*, 'slht', slht
          end do
        end do

        ralsdcp=slht/rcpd
        ralvdcp=vlht/rcpd
        zrldcp  = 1.0/(ralsdcp-ralvdcp)

        do j = jstart, jend
          do i = istart, iend
            zfrzmax(j,i) = max((zt0-zt(j,i,k))*zrldcp,0.0)
            !print*, 'zfrmax', zfrzmax(j,i)
          end do
        end do

          
          do j = jstart, jend
            do i = istart, iend
              if (zfrzmax(j,i) > zepsec .AND. zqxfg(j,i,iqr) > zepsec) then
                zfrz=min(zqxfg(j,i,iqr),zfrzmax(j,i))
                zsolqa(j,i,iqs,iqr)=zsolqa(j,i,iqs,iqr) + zfrz
                zsolqa(j,i,iqr,iqs)=zsolqa(j,i,iqr,iqs) - zfrz
                !print*, 'sink of rain', zsolqa(j,i,iqs,iqr)
                !print*, 'source of rain', zsolqa(j,i,iqr,iqs)      
              end if 
            end do
          end do
          

            !!!!!!!!!!!!!!!!!!! Evaporate very small amounts of liquid and ice
        do j = jstart, jend
          do i = istart, iend
 
            if ( zqx(j,i,k,iql) < rlmin ) then
              zsolqa(j,i,iqv,iql) = zqx(j,i,k,iql)
              zsolqa(j,i,iql,iqv) = -zqx(j,i,k,iql)
            end if
 
            if ( zqx(j,i,k,iqi) < rlmin ) then
              zsolqa(j,i,iqv,iqi) = zqx(j,i,k,iqi)
              zsolqa(j,i,iqi,iqv) = -zqx(j,i,k,iqi)
            end if
        
          end do ! IY
        end do !JPX
             !!!!!!!
      end if ! lmicro
 
      !!!!!!!!!!!!! sedimentation/falling of microphysical species
      do n = 1 , nqx
        do j = jstart , jend
          do i = istart , iend
            zfallsink(j,i,n) = ptsphy*vqx(j,i,k,n)/dzeta(j,i,k)
            if ( k > 1 ) then
              zfallsrce(j,i,n) = ptsphy*vqx(j,i,k,n) * &
                               zqx(j,i,k-1,n)/dzeta(j,i,k)
            end if
          end do
        end do
      end do
 
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
              zratio(j,i,n) = zqx(j,i,k,n)/zsinksum(j,i,n)
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
            zqxn(j,i,k,n) = zqx(j,i,k,n) + zexplicit
            ! note ***density missing*** for the moment...
            zfluxq(j,i,n) = zfallsrce(j,i,n)   !(delz at k-1)
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
 
      call ludcmp(jstart,jend,istart,iend,zqlhs,jindex2)
      call lubksb(jstart,jend,istart,iend,zqlhs,jindex2,zqxn)
 
      ! calculate fluxes in and out of box for conservation of TL
      do n = 1 , nqx
        do j = jstart , jend
          do i = istart , iend
            !   Note ***explicit PLUDE*** fluxes missing for moment
            zfluxq(j,i,n) = zfluxq(j,i,n) - ptsphy*vqx(j,i,k,n) * &
                            zqxn(j,i,k,n)/dzeta(j,i,n)
            ! if ( kphase(n) == 1 ) then
            !   zt(j,i,k) = zt(j,i,k) - lvl*zfluxq(j,i,n)*rcpd
            !   ! liquid terms
            !   zt(j,i,k) = zt(j,i,k)+lvl*(zqxn(j,i,k,n)-zqx(j,i,k,n))*rcpd
            ! end if
            ! if ( kphase(n) == 2 ) then
            !   ztl(k) = ztl(k) + lvs*zfluxq(n)*rcpd  !ice terms
            ! end if
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
      ztl(jci1:jci2,ici1:ici2,k) = zt(jci1:jci2,ici1:ici2,k)
      do n = 1 , nqx
        do j = jstart , jend
          do i = istart , iend
            if ( kphase(n) == 1 ) then
              ztl(j,i,k) = ztl(j,i,k) - lvl*zqx(j,i,k,n)*rcpd
            end if
            if ( kphase(n) == 2 ) then
              ztl(j,i,k) = ztl(j,i,k) - lvs*zqx(j,i,k,n)*rcpd
            end if
          end do
        end do
      end do

      do n = 1 , nqx
        do j = jstart , jend
          do i = istart , iend
            zqx(j,i,k,n) = zqxn(j,i,k,n)
            wvflux(j,i,k,n) = rhob3d(j,i,k)*zqxn(j,i,k,n)*vqx(j,i,k,n)
          end do
       end do 
     end do 

  ! calculate cloud total water at the end
      do j = jstart , jend
        do i = istart , iend
          totalw2(j,i,k)= zqx(j,i,k,iqv)+zqx(j,i,k,iql)+zqx(j,i,k,iqr) + &
                          zqx(j,i,k,iqi)+zqx(j,i,k,iqs)
        end do
      end do
      w2=sum(totalw2)
   end do   ! kz : end of vertical loop
   !print*,' diff qt=', w2-w1  

 
    contains

     real(dp) function dqsatdtc(zt,satc)
       implicit none
       real(dp) , intent(in) :: satc , zt
       dqsatdtc = satc*(6150.0D0/(zt**d_two))
     end function dqsatdtc
      
     real(dp) function dqsatdtw(zt,satw)
       implicit none
       real(dp) , intent(in) :: satw , zt
       dqsatdtw = satw*(4097.99D0/((zt-32.15D0)**d_two))
     end function dqsatdtw

     real(dp) function satc(zt,xp)
       implicit none
       real(dp) , intent(in) :: xp , zt
       ! saturation mixing ratio in hPa
       satc = (3.79D0/xp)*dexp(22.514D0-6150.0D0/zt) 
     end function satc

     real(dp) function satw(zt,xp)
       implicit none
       real(dp) , intent(in) :: xp , zt
       ! saturation mixing ratio in hPa
       satw = (3.79D0/xp)*dexp(17.5D0*(zt-tzero)/(zt-32.15D0))
     end function satw
      
  end subroutine microphys

  subroutine lubksb(jstart,jend,istart,iend,aam,indx,bbm)
    implicit none
    integer , intent(in) :: jstart , jend , istart , iend
    real(dp) , pointer , intent(in) , dimension(:,:,:,:) :: aam
    integer , pointer , intent(in) , dimension(:,:,:) :: indx
    real(dp) , pointer , intent(inout) , dimension(:,:,:,:) :: bbm
    integer :: i , j , ii , jj , k , ll , m
    real(dp) :: xsum

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
    integer , intent(in) :: jstart , jend , istart , iend
    real(dp) , pointer , intent(inout) , dimension(:,:,:,:) :: aam
    integer , pointer , intent(out) , dimension(:,:,:) :: indx
!
    real(dp) :: aamax , dum , xsum
    integer :: i , j , k , imax , n , m
    real(dp) , dimension(nmax) :: vv
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
