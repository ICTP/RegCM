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
!

! fsolmon@ictp.it: this is the interface for calling rrtm

module mod_rrtmg_driver
!
  use mod_realkinds
  use mod_dynparam
  use mod_rad_common
  use mod_rad_scenarios
  use mod_rad_tracer
  use rrtmg_sw_rad
  use mcica_subcol_gen_sw
  use parrrsw
  use rrsw_wvn
  use parrrtm
  use rrtmg_lw_rad
  use mod_rad_outrad

private
!
  public :: rrtmg_driver

  contains
!
  subroutine rrtmg_driver(ktau,iyear,eccf,j,lout)

  implicit none

  integer(8) , intent(in) :: ktau
  integer , intent(in) :: iyear
  integer , intent(in) :: j
  logical, intent(in) :: lout
  real(dp), intent(in) :: eccf
! these are temporary passed to colmod to test RRTM SW only 
! when RRTMLW is implemented , radout will be called in this module. The driver will be called from tendency.
! 

!local variables

real(dp), dimension (iym1)  :: solin,frsa,sabtp,clrst,clrss, firtp, frla, clrlt, clrls , empty1,srfrad,sols,soll,solsd,solld,slwd
real(dp), dimension (iym1,kz)  :: qrs  ,qrl, clwp_int, cld_int,empty2


integer:: dyofyr , inflgsw, iceflgsw, liqflgsw, icld, idrv, permuteseed,irng,iplon,ns,k,kj, &
          inflglw ,iceflglw,liqflglw
 
 real(dp) :: adjes

 real(dp) , dimension(iym1,kzp1) :: plev,tlev,swuflx,swdflx,swuflxc,swdflxc,&
                                    lwuflx ,lwdflx, lwuflxc ,lwdflxc ,&
                                    swddiruviflx, swddifuviflx ,& 
                                    swddirpirflx, swddifpirflx,swdvisflx

 real(dp) , dimension(iym1,kz) :: play,tlay,h2ovmr , o3vmr   ,co2vmr  ,ch4vmr  ,n2ovmr ,o2vmr , & 
                                 cfc11vmr, cfc12vmr,  cfc22vmr, ccl4vmr,    &
                                 reicmcl, relqmcl,swhr,swhrc,ciwp,clwp,rei,rel,cldf,  &
                                 lwhr,lwhrc,duflx_dt,duflxc_dt
                               
 real(dp) , dimension(iym1):: tsfc, psfc

 
 real(dp) , dimension(ngptsw,iym1,kz) :: cldfmcl, taucmcl, ssacmcl, asmcmcl, fsfcmcl, &
                                          ciwpmcl, clwpmcl 


! spectral dependant quantities
 real(dp) , dimension(iym1,kz,nbndsw) :: tauaer, ssaaer, asmaer, ecaer
 real(dp), dimension(nbndsw,iym1,kz) :: tauc,ssac,asmc,fsfc      


real(dp), dimension(iym1,nbndlw) :: emis_surf 
real(dp), dimension(iym1) :: asdir 
real(dp), dimension(iym1) :: asdif 
real(dp), dimension(iym1) :: aldir 
real(dp), dimension(iym1) :: aldif 
real(dp), dimension(iy) :: czen 
real(dp), dimension(nbndlw,iym1,kz) :: tauc_lw 
real(dp),dimension( iym1,kz,nbndlw) :: tauaer_lw 

real(dp) , dimension(ngptlw,iym1,kz) :: cldfmcl_lw, taucmcl_lw, &
                                          ciwpmcl_lw, clwpmcl_lw 
!-----------------------------------------------------------------------


iplon=1 ! not effectively used

!adjustment of earthsun distance:
! if dyofyr = 0, we pass directly the ecc. factor in adjes, 
! otherwise it is calculated in RRTM as a function of day of the year.
dyofyr = 0
adjes = eccf

! options for optical properties of cloud calculation

!inflgsw  = 0 : use the optical properties calculated in perp_dat_rrtm (same as standard radiation)
!inflgsw  = 2 : use RRTM option to calculate cloud optical prperties from water path and cloud drop radius
!check the diferent param available in rrtmg_rad (iceflgsw / liqflgsw  should be nameliste interactive) 
inflgsw  = 2 
iceflgsw = 3
liqflgsw = 1

! IN LW : optical depth is calculated internally 
!from water path and cloud radius / tauc_LW is not requested
inflglw  = 2
iceflglw = 3
liqflglw = 1
tauc_lw = dlowval
!
!McICA parameteres
icld  = 1 ! Cloud Overlapp hypothesis ( should be interactive) 
!
irng = 1 ! mersenne twister random generator for McICA
! Cloud Overlapp hypothesis flag ( return if 0)


call prep_dat_rrtm (j,iyear,inflgsw,plev,play,tlev,tlay,tsfc,psfc,         &
                     cldf,h2ovmr,o3vmr,n2ovmr,co2vmr,ch4vmr,o2vmr, & 
                     cfc11vmr,cfc12vmr, cfc22vmr, ccl4vmr, &
                     tauc,ssac,asmc,fsfc,ciwp,clwp,rei,rel)


! Call to the shortwave radiation code as soon one element of coszen is > 0.
! 

 if ( maxval (coszen) > dlowval) then

! generates cloud properties:
permuteseed = 1
call mcica_subcol_sw( iplon, iym1, kz, icld, permuteseed, irng, play, &
                       cldf , ciwp, clwp, rei, rel, tauc, ssac, asmc, fsfc, &
                       cldfmcl, ciwpmcl, clwpmcl, reicmcl, relqmcl, &
                       taucmcl, ssacmcl, asmcmcl, fsfcmcl)

! for now initialise aerosol OP to zero:
tauaer = d_zero
ssaaer = d_one
asmaer = 0.85D0
ecaer =0.78D0

asdir = swdiralb(j,:)
asdif = swdifalb(j,:)
aldir = lwdiralb(j,:)
aldif = lwdifalb(j,:)
czen = coszen(j,:)

     call rrtmg_sw &
            (iym1    , kz   ,icld    , &
             play    ,plev   ,tlay    ,tlev    ,tsfc   , &
             h2ovmr , o3vmr   ,co2vmr  ,ch4vmr  ,n2ovmr ,o2vmr , &
             asdir  , asdif   , aldir  ,aldif   , &
             czen  ,adjes   ,dyofyr  ,solcon    , &
             inflgsw ,iceflgsw,liqflgsw,cldfmcl , &
             taucmcl ,ssacmcl ,asmcmcl ,fsfcmcl , &
             ciwpmcl ,clwpmcl ,reicmcl ,relqmcl , &
             tauaer  ,ssaaer  ,asmaer  ,ecaer   , &
             swuflx  ,swdflx  ,swhr    ,swuflxc ,swdflxc ,swhrc, &
             swddiruviflx, swddifuviflx ,& 
             swddirpirflx, swddifpirflx, swdvisflx     )

      
end if ! end shortwave call 


! LW call :

  permuteseed = 150
  call  mcica_subcol_lw(iplon, iym1, kz, icld, permuteseed, irng, play, &
                       cldf, ciwp, clwp, rei, rel, tauc_lw, cldfmcl_lw, &
                       ciwpmcl_lw, clwpmcl_lw, reicmcl, relqmcl, taucmcl_lw)

tauaer_lw = d_zero

!  provisoire
do k = 1,nbndlw
emis_surf(:,k) = 0.97D0
! = emsvt(:)
end do 

  idrv = 0
  call  rrtmg_lw &
            (iym1    ,kz    ,icld    ,idrv    , &
             play    ,plev    ,tlay    ,tlev    ,tsfc    , & 
             h2ovmr  ,o3vmr   ,co2vmr  ,ch4vmr  ,n2ovmr  ,o2vmr , &
             cfc11vmr,cfc12vmr,cfc22vmr,ccl4vmr ,emis_surf   , &
             inflglw ,iceflglw,liqflglw,cldfmcl_lw , &
             taucmcl_lw ,ciwpmcl_lw ,clwpmcl_lw ,reicmcl ,relqmcl , &
             tauaer_lw  , &
             lwuflx ,lwdflx,lwhr, lwuflxc ,lwdflxc,lwhrc, &    
             duflx_dt,duflxc_dt )


! Output and interface 
! 
! solin  - instantaneous incident solar
! sabtp  - total column absorbed solar flux
! frsa   - surface absorbed solar flux
! clrst  - clear sky total column abs solar flux
! clrss  - clear sky surface absorbed solar flux
! qrs    - solar heating rate
! firtp  - net up flux top of model (up-dwn flx)
! frla   - longwave cooling of surface (up-dwn flx)
! clrlt  - clr sky net up flx top of model (up-dwn f
! clrls  - clr sky lw cooling of srf (up-dwn flx)
! qrl    - longwave cooling rate
! slwd   - surface longwave down flux
! srfrad - surface radiative heating flux (frsa+slwd)
! h2ommr - ozone mixing ratio
! cld    - cloud fractional cover
! clwp   - cloud liquid water path
! soll   - Downward solar rad onto surface (lw direct)
! solld  - Downward solar rad onto surface (lw diffuse)
! sols   - Downward solar rad onto surface (sw direct)
! solsd  - Downward solar rad onto surface (sw diffuse)
! solar = visible band only downard SW radiation 
!
!EES  next 3 added, they are calculated in radcsw : but not used further
! fsnirt   - Near-IR flux absorbed at toa
! fsnrtc   - Clear sky near-IR flux absorbed at toa
! fsnirtsq - Near-IR flux absorbed at toa >= 0.7 microns
! fsds     - Flux Shortwave Downwelling Surface
!

solin(:) = swdflx(:,kzp1)
frsa(:) =   swdflx(:,1) - swuflx(:,1)
sabtp(:) =  swdflx(:,kzp1) -  swuflx(:,kzp1)
clrst(:) =   swdflxc(:,kzp1) -  swuflxc(:,kzp1)
clrss(:) =   swdflxc(:,1) -  swuflxc(:,1)

firtp(:) =  -d_one*(lwdflx(:,kzp1) -  lwuflx(:,kzp1))
frla(:) =   -d_one*(lwdflx(:,1) - lwuflx(:,1))
!!$
clrlt(:) =  -d_one* (lwdflxc(:,kzp1) -  lwuflxc(:,kzp1))
clrls(:) =   -d_one*(lwdflxc(:,1) -  lwuflxc(:,1))

! coupling with  BATS
!  abveg set to frsa (as in standard version : potential inconsieny if soil fraction is large)
abveg(j,:) = frsa(:) 
! solar is normally the visible band only total incident surface flux
solar(j,:) = swdvisflx(:,1)
! surface SW incident
sols(:) =  swddiruviflx(:,1)
solsd(:) =  swddifuviflx(:,1)
soll(:) =  swddirpirflx(:,1)
solld(:) =  swddifpirflx(:,1)
!LW incident 
slwd(:) = lwdflx(:,1)

!3d heating rate back on regcm grid and converted to K.S-1
do k=1,kz
kj = kzp1-k
qrs(:,kj) = swhr(:,k) / secpd
qrl(:,kj) = lwhr(:,k) / secpd

cld_int(:,kj) = cldf(:,k) !ouptut : these are in cloud diagnostics  
clwp_int(:,kj) = clwp(:,k)
end do 


! Finally call radout for coupling to BATS/CLM/ATM and outputing fields in RAD files
! these are diagnostics that are passed to radout but not used furthermore
empty1 = dmissval
empty2 = dmissval
!!$
!!$
    call radout (lout,solin,sabtp,frsa,clrst,clrss,qrs,firtp,frla,  &
                    clrlt,clrls,qrl,slwd,srfrad,sols,soll,solsd,solld, &
                    empty1,empty1,empty1,empty1,empty1,empty1,j,empty2,cld_int,clwp_int)



end subroutine  rrtmg_driver




subroutine prep_dat_rrtm (j,iyear,inflagsw,plev,play,tlev,tlay,ts,ps,         &
                          cldf,h2ovmr,o3vmr,n2ovmr,co2vmr,ch4vmr,o2vmr,&
                          cfc11vmr,cfc12vmr,cfc22vmr, ccl4vmr,                   &                  
                          tauc,ssac,asmc,fsfc,ciwp,clwp,rei,rel )
!
    implicit none
!
    integer :: j
    integer , intent(in) :: iyear, inflagsw
    real(dp) , dimension(iym1) :: alat , coslat ,  ps , ts
    real(dp) , dimension(iym1,kzp1) ::  plev,tlev
    real(dp) , dimension(iym1,kz):: play  , tlay,  &
                                   h2ovmr,o3vmr ,n2ovmr,co2vmr,ch4vmr,o2vmr, &
                                   cfc11vmr,cfc12vmr, cfc22vmr, ccl4vmr, &
                                   rei,rel,fice, &
                                   wcl,wci,gcl,gci,fcl,fci, tauxcl,tauxci
 

 real(dp)  :: tauc (nbndsw,iym1,kz)        ! in-cloud optical depth
!!$                                                      !    Dimensions: (nbndsw,ncol,nlay)
 real(dp)  :: ssac(nbndsw,iym1,kz)        ! in-cloud single scattering albedo (non-delta scaled)
!!$                                                      !    Dimensions: (nbndsw,ncol,nlay) 
 real(dp)  :: asmc (nbndsw,iym1,kz)       ! in-cloud asymmetry parameter (non-delta scaled)
!!$                                                      !    Dimensions: (nbndsw,ncol,nlay)
 real(dp)  :: fsfc(nbndsw,iym1,kz)       ! in-cloud forward scattering fraction (non-delta scaled)
!!$                                                      !    Dimensions: (nbndsw,ncol,nlay)
 real(dp)  :: ciwp(iym1,kz)           ! in-cloud ice water path
!!$                                                      !    Dimensions: (ncol,nlay)
 real(dp)  :: clwp(iym1,kz)        ! in-cloud liquid water path
!!$                                                      !    Dimensions: (ncol,nlay)
real(dp)  :: cldf(iym1,kz)


    intent (in) j
    intent (out)  o3vmr ,n2ovmr,co2vmr,ch4vmr,o2vmr ,cfc11vmr,cfc12vmr, cfc22vmr, ccl4vmr, ts , cldf , h2ovmr ,& 
       plev , play , ps , tlay, tauc, ssac, asmc, fsfc, ciwp,clwp, rei,  rel
!
    
      
    real(dp) , dimension(iym1,kz) ::  h2ommr , o3mmr, n2ommr,co2mmr,ch4mmr, cfc11mmr,cfc12mmr
    real(dp) , dimension(iym1,kz) :: deltaz
    real(dp) :: c287,ccvtem,clwtem,w1,w2,wmin,wmax
    integer :: i , k , kj , ncldm1, indxsl,ns,n,nrr,ii0,ii1,ii2
    integer , dimension(iym1) :: ioro
    real(dp) , parameter :: lowcld = 1.0D-30
    real(8) , parameter :: verynearone = 0.999999D0
    real(dp) :: tmp1l,tmp2l,tmp3l, tmp1i,tmp2i,tmp3i


!
!     Set index for cloud particle properties based on the wavelength,
!     according to A. Slingo (1989) equations 1-3:
!     Use index 1 (0.25 to 0.69 micrometers) for visible
!     Use index 2 (0.69 - 1.19 micrometers) for near-infrared
!     Use index 3 (1.19 to 2.38 micrometers) for near-infrared
!     Use index 4 (2.38 to 4.00 micrometers) for near-infrared
!
!     Note that the minimum wavelength is encoded (with 0.001, 0.002,
!     0.003) in order to specify the index appropriate for the
!     near-infrared cloud absorption properties
!
!the 14 RRTM SW bands intervall are given in wavenumber (wavenum1/2): 
!equivalence in micrometer are
!3.07-3.84 / 2.50-3.07 / 2.15-2.5 / 1.94-2.15 / 1.62-1.94 / 1.29-1.62 / 1.24-1.29 /
!0.778-1.24 / 0.625-0.778 / 0.441-0.625 / 0.344-0.441 / 0.263-0.344 / 0.200-0.263 / 
!3.846-12.195 

 integer, dimension (nbndsw) :: indsl
   
! A. Slingo's data for cloud particle radiative properties
! (from 'A GCM Parameterization for the Shortwave Properties of Water
! Clouds' JAS vol. 46 may 1989 pp 1419-1427)
!
! abarl    - A coefficient for extinction optical depth
! bbarl    - B coefficient for extinction optical depth
! cbarl    - C coefficient for single particle scat albedo
! dbarl    - D coefficient for single particle scat albedo
! ebarl    - E coefficient for asymmetry parameter
! fbarl    - F coefficient for asymmetry parameter
!
! Caution... A. Slingo recommends no less than 4.0 micro-meters nor
! greater than 20 micro-meters
!/nnr

! ice water coefficients (Ebert and Curry,1992, JGR, 97, 3831-3836)
!
! abari    - a coefficient for extinction optical depth
! bbari    - b coefficient for extinction o


real(dp) , dimension(4) ::  abari , abarl , bbari , bbarl , cbari , &
                            cbarl , dbari , dbarl , ebari , ebarl , &
                            fbari , fbarl


real(dp) :: abarii , abarli , bbarii , bbarli , cbarii , cbarli ,  &
               dbarii , dbarli , ebarii , ebarli , fbarii ,   &
               fbarli 

  data abarl / 2.817D-02 ,  2.682D-02 , 2.264D-02 , 1.281D-02/
  data bbarl / 1.305D+00 ,  1.346D+00 , 1.454D+00 , 1.641D+00/
  data cbarl /-5.620D-08 , -6.940D-06 , 4.640D-04 , 0.201D+00/
  data dbarl / 1.630D-07 ,  2.350D-05 , 1.240D-03 , 7.560D-03/
  data ebarl / 0.829D+00 ,  0.794D+00 , 0.754D+00 , 0.826D+00/
  data fbarl / 2.482D-03 ,  4.226D-03 , 6.560D-03 , 4.353D-03/
! 
  data abari / 3.4480D-03 , 3.4480D-03 , 3.4480D-03 , 3.44800D-03/
  data bbari / 2.4310D+00 , 2.4310D+00 , 2.4310D+00 , 2.43100D+00/
  data cbari / 1.0000D-05 , 1.1000D-04 , 1.8610D-02 , 0.46658D+00/
  data dbari / 0.0000D+00 , 1.4050D-05 , 8.3280D-04 , 2.05000D-05/
  data ebari / 0.7661D+00 , 0.7730D+00 , 0.7940D+00 , 0.95950D+00/
  data fbari / 5.8510D-04 , 5.6650D-04 , 7.2670D-04 , 1.07600D-04/

! define index pointing on appropriate parameter in slingo's table for eachRRTM SW band
!
data indsl /4,4,3,3,3,3,3,2,2,1,1,1,1,4 /

!
   !
    do i = 1 , iym1
      alat(i) = xlat(i,j)*degrad
      coslat(i) = dcos(alat(i))
    end do


! CONVENTION : RRTMG driver takes layering form botom to TOA. 
! regcm consider Top to bottom


!   surface pressure and scaled pressure, from which level are computed, RRTM SW takes pressure in mb,hpa
    do i = 1 , iym1
      ps(i) = (sfps(i,j)+ptp)*d_10
     do k = 1 , kz
        kj = kzp1 - k        
        play(i,kj) = (sfps(i,j)*hlev(k)+ptp)*d_10
      end do
    end do
!
!   convert pressures from mb to pascals and define interface pressures:
!

    do k = 1 , kzp1
      do i = 1 , iym1
        kj = kzp1 - k  +1   
        plev(i,kj) = (sfps(i,j)*flev(k)+ptp)*d_10
      end do
    end do
!

!   ground temperature
!
    do i = 1 , iym1
      ts(i) = tground(i,j)
    end do
!

!   air temperatures
!
    do k = 1 , kz
     kj = kzp1 - k 
     do i = 1 , iym1
        tlay(i,kj) = tatms(i,k,j)
      end do
    end do
!
! air temperature at the interface (same method as in vertical advection routine)
          c287 = 0.287D+00
          do k = 2 , kz
            kj = kzp1 - k   +1        
            do i = 1,iym1
            w1 =  (hlev(kj) - flev(kj)) / (hlev(kj) - hlev(kj-1))
            w2 =  (flev(kj) - hlev(kj-1) ) / (hlev(kj) - hlev(kj-1))
            
             if (k < kz-1) then    
              tlev(i,k) =  w1*tlay(i,k-1) * ( plev(i,k)/ play(i,k-1) )**c287 + &
                           w2*tlay(i,k)   * ( plev(i,k)/ play(i,k) )**c287     
             else 
!linear interpolation for last upper tropospheric levels
               tlev(i,k) =  w1*tlay(i,k-1)  + &
                           w2*tlay(i,k)     
             end if 
            end do
           end do
          
         do i = 1,iym1
          tlev(i,1) = ts(i)
          tlev(i,kzp1) = tlay(i,kz)
         end do

!   h2o volume mixing ratio
!
    do k = 1 , kz 
        kj = kzp1 - k       
     do i = 1 , iym1
        h2ommr(i,kj) = dmax1(1.0D-7,qvatms(i,k,j))
        h2ovmr(i,kj) = h2ommr(i,kj) * ep2
      end do
    end do
!
!   o3 volume mixing ratio (already on the right grid)
!
    do k = 1 , kz
      do i = 1 , iym1
!        kj = kzp1 - k
        o3vmr(i,k) = o3prof(i,k,j) * amo/amd
      end do
    end do
!
! other gas (n2o,ch4)
!
  call trcmix  (play,alat,coslat,n2ommr,ch4mmr,cfc11mmr,cfc12mmr)

 do k = 1 , kz
 do i = 1 , iym1

  n2ovmr (i,k) = n2ommr(i,k) * (44.D0/amd)
  ch4vmr (i,k) =  ch4mmr (i,k) * (16.D0/amd)
  co2vmr(i,k)  = cgas(2,iyear)*1.0D-6
  o2vmr(i,k) =  0.23143D0 * (32.D0/amd)
  cfc11vmr(i,k) =  cfc11mmr (i,k) / (163.1278D0 /amd)
  cfc12vmr(i,k) =  cfc12mmr (i,k) / (175.1385D0 /amd)
!
! No data FOR NOW : IMPROVE !!
   cfc22vmr(i,k) = d_zero
   ccl4vmr (i,k) =d_zero 

  end do 
 end do


!cloud fraction and cloud liquid waterpath calculation:
!=  as in STANDARD SCHEME for now (getdat) : We need to improve this 
!according to new cloud microphysics!

!   fractional cloud cover (dependent on relative humidity)
!
!   qc   = gary's mods for clouds/radiation tie-in to exmois
    do k = 1 , kz
        kj = kzp1 - k
      do i = 1 , iym1   
        ccvtem = d_zero   !cqc mod

        cldf(i,kj) = dmax1(cldfra(i,k)*0.9999999D0,ccvtem)

        cldf(i,kj) = dmin1(cldf(i,kj),0.9999999D0)
!
!       convert liquid water content into liquid water path, i.e.
!       multiply b deltaz
        clwtem = cldlwc(i,kj) ! put on the right grid !

!deltaz,clwp are on the right grid since plev and tlay are
!care pressure is on botom/toa grid
        deltaz(i,k) = rgas*tlay(i,k)*(plev(i,k) - &
                        plev(i,k+1))/(egrav*play(i,k))
        clwp(i,k) = clwtem * deltaz(i,k)

        if ( dabs(cldf(i,k)) < lowcld ) then
          cldf(i,k) = d_zero
          clwp(i,k) = d_zero
        end if
      end do
    end do
   

!!$!   set cloud fractional cover at top model level = 0
!!$    do i = 1 , iym1
!!$      cld(i,kzp1) = d_zero
!!$      clwp(i,kzp1) = d_zero
!!$      cld(i,kz) = d_zero      
!!$      clwp(i,kz) = d_zero
!!$    end do
!
!  partition Take into accountice water path 
  


!   set cloud fractional cover at bottom (ncld) model levels = 0
!
    ncldm1 = ncld - 1
    do k = 1 ,ncldm1
      do i = 1 , iym1
        cldf(i,k) = d_zero
        clwp(i,k) = d_zero
      end do
    end do
!
!maximum cloud fraction
!----------------------------------------------------------------------
    do i = 1 , iym1
      do k = 1 , kz
        if ( cldf(i,k) > 0.999D0 ) cldf(i,k) = 0.999D0
      end do
    end do
!
!
! CLOUD Properties:
!
! cloud effective radius 
!   NB: orography types are specified in the following
!
    do i = 1 , iym1
      ii0 = 0
      ii1 = 0
      ii2 = 0
      do n = 1 , nnsg
        if ( lndocnicemsk(n,i,j) == 2 ) then
          ii2 = ii2 + 1
        else if ( lndocnicemsk(n,i,j) == 1 ) then
          ii1 = ii1 + 1
        else
          ii0 = ii0 + 1
        end if
      end do
      if ( ii0 >= ii1 .and. ii0 >= ii2 ) ioro(i) = 0
      if ( ii1 > ii0 .and. ii1 >= ii2 ) ioro(i) = 1
      if ( ii2 > ii0 .and. ii2 > ii1 ) ioro(i) = 2
    end do

!
    call cldefr_rrtm(ioro,tlay,rel,rei,fice,sfps(:,j),play)
!
!
!  partition of total  water path betwwen liquide and ice.
! ( waiting for prognostic ice !) 
!
     do i = 1 , iym1
      do k = 1 , kz
        ciwp (i,k) =  clwp (i,k) *   fice (i,k)
        clwp (i,k) =  clwp (i,k) * (1- fice (i,k))
        ! now clwp is liquide only !
     end do
    end do

! 
! Cloud optical properties(tau,ssa,g,f) :
! 2 options :
! inflgsw == 0 :  treated as the standard radiation scheme and passed to McICA/ RRTM
! inflgsw == 2 : direcly calculated within RRTMsw 



!initialise and  begin spectral loop
!
        tauc  =  d_zero
        ssac  =  verynearone
        asmc  =  0.850D0
        fsfc  =  0.725D0
           
if (inflagsw == 0) then 

    do ns = 1 ,nbndsw
 
!      nrr = ns + jpb1 -1
!     Set cloud extinction optical depth, single scatter albedo,
!     asymmetry parameter, and forward scattered fraction:
!
      abarli = abarl(indsl(ns))
      bbarli = bbarl(indsl(ns))
      cbarli = cbarl(indsl(ns))
      dbarli = dbarl(indsl(ns))
      ebarli = ebarl(indsl(ns))
      fbarli = fbarl(indsl(ns))
!
      abarii = abari(indsl(ns))
      bbarii = bbari(indsl(ns))
      cbarii = cbari(indsl(ns))
      dbarii = dbari(indsl(ns))
      ebarii = ebari(indsl(ns))
      fbarii = fbari(indsl(ns))
!
      do k = 1 , kz
        do i =  1,iym1
!
         if (clwp(i,k) < dlowval .and.  ciwp(i,k)< dlowval) cycle 

!           liquid
!
            tmp1l = abarli + bbarli/rel(i,k)
            tmp2l = d_one - cbarli - dbarli*rel(i,k)
            tmp3l = fbarli*rel(i,k)
!
!           ice
!
            tmp1i = abarii + bbarii/rei(i,k)
            tmp2i = d_one - cbarii - dbarii*rei(i,k)
            tmp3i = fbarii*rei(i,k)
!
!! cloud optical properties extinction optical depth : IN_CLOUD quantities !

            tauxcl(i,k) = clwp(i,k)*tmp1l
            tauxci(i,k) = ciwp(i,k)*tmp1i
   
            wcl(i,k) = dmin1(tmp2l,verynearone)
            gcl(i,k) = ebarli + tmp3l
            fcl(i,k) = gcl(i,k)*gcl(i,k)
!
            wci(i,k) = dmin1(tmp2i,verynearone)
            gci(i,k) = ebarii + tmp3i
            fci(i,k) = gci(i,k)*gci(i,k)
!

            tauc (ns,i,k) =  tauxcl(i,k) + tauxci(i,k) 
        
            ssac (ns,i,k) =  (tauxcl(i,k) * wcl(i,k) + tauxci(i,k) * wci(i,k) ) /  tauc (ns,i,k)
            asmc (ns,i,k) =  (tauxcl(i,k) * gcl(i,k) + tauxci(i,k) * gci(i,k) ) /  tauc (ns,i,k)
            fsfc (ns,i,k) =  (tauxcl(i,k) * fcl(i,k) + tauxci(i,k) * fci(i,k) ) /  tauc (ns,i,k)
           
          end do
       
      end do
!
!     Set reflectivities

     end do ! spectral loop

end if  ! inflagsw

!
  end subroutine prep_dat_rrtm
!
!
  subroutine cldefr_rrtm(ioro,t,rel,rei,fice,ps,pmid)
!
! for now we use for RRTM the same param as in standard rad 
    implicit none
!
    real(dp) , dimension(iym1,kz) :: fice , pmid , rei , rel , t
    integer , dimension(iym1) :: ioro
    real(dp) , dimension(iym1) :: ps
    intent (in) ioro , pmid , ps , t
    intent (out) fice , rei , rel
!
    integer :: i , k
    real(dp) :: pnrml , rliq , weight
!
!   reimax - maximum ice effective radius
    real(dp) , parameter :: reimax = 30.0D0
!   rirnge - range of ice radii (reimax - 10 microns)
    real(dp) , parameter :: rirnge = 20.0D0
!   pirnge - nrmlzd pres range for ice particle changes
    real(dp) , parameter :: pirnge = 0.4D0
!   picemn - normalized pressure below which rei=reimax
    real(dp) , parameter :: picemn = 0.4D0
!   Temperatures in K (263.16 , 243.16)
    real(dp) , parameter :: minus10 = wattp-d_10
    real(dp) , parameter :: minus30 = wattp-(d_three*d_10)
!
    do k = 1 , kz
      do i = 1 , iym1
!
!       Define liquid drop size
!
        if ( ioro(i) /= 1 ) then
!
!         Effective liquid radius over ocean and sea ice
!
          rliq = d_10
        else
!
!         Effective liquid radius over land
!
          rliq = d_five+d_five* & 
                  dmin1(d_one,dmax1(d_zero,(minus10-t(i,k))*0.05D0))
        end if
!
        rel(i,k) = rliq
!fil
!       test radius = d_10
!       rel(i,k) = d_10
!fil
!+      rei(i,k) = 30.0
!
!       Determine rei as function of normalized pressure
!
        pnrml = pmid(i,k)/ps(i)
        weight = dmax1(dmin1((pnrml-picemn)/pirnge,d_one),d_zero)
        rei(i,k) = reimax - rirnge*weight
!
!       Define fractional amount of cloud that is ice
!
!       if warmer than -10 degrees C then water phase
!
        if ( t(i,k) > minus10 ) fice(i,k) = d_zero
!
!       if colder than -10 degrees C but warmer than -30 C mixed phase
!
        if ( t(i,k) <= minus10 .and. t(i,k) >= minus30 ) &
          fice(i,k) = (minus10-t(i,k))/20.0D0
!
!       if colder than -30 degrees C then ice phase
!
        if ( t(i,k) < minus30 ) fice(i,k) = d_one
!
!       Turn off ice radiative properties by setting fice = 0.0
!
!fil    no-ice test
!       fice(i,k) = d_zero
!
      end do
    end do
!
  end subroutine cldefr_rrtm
!


end module mod_rrtmg_driver
