module mod_clm_pftvarcon

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: pftvarcon
!
! !DESCRIPTION:
! Module containing vegetation constants and method to
! read and initialize vegetation (PFT) constants.
!
! !USES:
  use mod_realkinds
  use mod_mpmessage
  use mod_nchelper
  use mod_dynparam
  use mod_mppparam
  use mod_clm_varpar  , only : mxpft, numrad, ivis, inir
!
! !PUBLIC TYPES:
  implicit none
  save
!
! Vegetation type constants
!
  character(len=40) pftname(0:mxpft) !PFT description

  integer :: noveg                  !value for not vegetated 
  integer :: ndllf_evr_tmp_tree     !value for Needleleaf evergreen temperate tree
  integer :: ndllf_evr_brl_tree     !value for Needleleaf evergreen boreal tree
  integer :: ndllf_dcd_brl_tree     !value for Needleleaf deciduous boreal tree
  integer :: nbrdlf_evr_trp_tree    !value for Broadleaf evergreen tropical tree
  integer :: nbrdlf_evr_tmp_tree    !value for Broadleaf evergreen temperate tree
  integer :: nbrdlf_dcd_trp_tree    !value for Broadleaf deciduous tropical tree
  integer :: nbrdlf_dcd_tmp_tree    !value for Broadleaf deciduous temperate tree
  integer :: nbrdlf_dcd_brl_tree    !value for Broadleaf deciduous boreal tree
  integer :: ntree                  !value for last type of tree
  integer :: nbrdlf_evr_shrub       !value for Broadleaf evergreen shrub
  integer :: nbrdlf_dcd_tmp_shrub   !value for Broadleaf deciduous temperate shrub
  integer :: nbrdlf_dcd_brl_shrub   !value for Broadleaf deciduous boreal shrub
  integer :: nc3_arctic_grass       !value for C3 arctic grass
  integer :: nc3_nonarctic_grass    !value for C3 non-arctic grass
  integer :: nc4_grass              !value for C4 grass
  integer :: npcropmin              !value for first crop
  integer :: ncorn                  !value for corn, rain fed (rf)
  integer :: ncornirrig             !value for corn, irrigated (ir)
  integer :: nscereal               !value for spring temperate cereal (rf)
  integer :: nscerealirrig          !value for spring temperate cereal (ir)
  integer :: nwcereal               !value for winter temperate cereal (rf)
  integer :: nwcerealirrig          !value for winter temperate cereal (ir)
  integer :: nsoybean               !value for soybean (rf)
  integer :: nsoybeanirrig          !value for soybean (ir)
  integer :: npcropmax              !value for last prognostic crop in list
  integer :: nc3crop                !value for generic crop (rf)
  integer :: nc3irrig               !value for irrigated generic crop (ir)

  real(rk8):: dleaf(0:mxpft)       !characteristic leaf dimension (m)
  real(rk8):: c3psn(0:mxpft)       !photosynthetic pathway: 0. = c4, 1. = c3
  real(rk8):: xl(0:mxpft)          !leaf/stem orientation index
  real(rk8):: rhol(0:mxpft,numrad) !leaf reflectance: 1=vis, 2=nir
  real(rk8):: rhos(0:mxpft,numrad) !stem reflectance: 1=vis, 2=nir
  real(rk8):: taul(0:mxpft,numrad) !leaf transmittance: 1=vis, 2=nir
  real(rk8):: taus(0:mxpft,numrad) !stem transmittance: 1=vis, 2=nir
  real(rk8):: z0mr(0:mxpft)        !ratio of momentum roughness length to canopy top height (-)
  real(rk8):: displar(0:mxpft)     !ratio of displacement height to canopy top height (-)
  real(rk8):: roota_par(0:mxpft)   !CLM rooting distribution parameter [1/m]
  real(rk8):: rootb_par(0:mxpft)   !CLM rooting distribution parameter [1/m]
  real(rk8):: crop(0:mxpft)        ! crop pft: 0. = not crop, 1. = crop pft
  real(rk8):: irrigated(0:mxpft)   ! irrigated pft: 0. = not, 1. = irrigated
  real(rk8):: smpso(0:mxpft)       !soil water potential at full stomatal opening (mm)
  real(rk8):: smpsc(0:mxpft)       !soil water potential at full stomatal closure (mm)
  real(rk8):: fnitr(0:mxpft)       !foliage nitrogen limitation factor (-)
  ! begin new pft parameters for CN code
  real(rk8):: slatop(0:mxpft)      !SLA at top of canopy [m^2/gC]
  real(rk8):: dsladlai(0:mxpft)    !dSLA/dLAI [m^2/gC]
  real(rk8):: leafcn(0:mxpft)      !leaf C:N [gC/gN]
  real(rk8):: flnr(0:mxpft)        !fraction of leaf N in Rubisco [no units]
  real(rk8):: woody(0:mxpft)       !woody lifeform flag (0 or 1)
  real(rk8):: lflitcn(0:mxpft)      !leaf litter C:N (gC/gN)
  real(rk8):: frootcn(0:mxpft)      !fine root C:N (gC/gN)
  real(rk8):: livewdcn(0:mxpft)     !live wood (phloem and ray parenchyma) C:N (gC/gN)
  real(rk8):: deadwdcn(0:mxpft)     !dead wood (xylem and heartwood) C:N (gC/gN)
  real(rk8):: grperc(0:mxpft)       !growth respiration parameter
  real(rk8):: grpnow(0:mxpft)       !growth respiration parameter
  real(rk8):: rootprof_beta(0:mxpft)   !CLM rooting distribution parameter for C and N inputs [unitless]

! for crop
  real(rk8):: graincn(0:mxpft)      !grain C:N (gC/gN)
  real(rk8):: mxtmp(0:mxpft)        !parameter used in accFlds
  real(rk8):: baset(0:mxpft)        !parameter used in accFlds
  real(rk8):: declfact(0:mxpft)     !parameter used in CNAllocation
  real(rk8):: bfact(0:mxpft)        !parameter used in CNAllocation
  real(rk8):: aleaff(0:mxpft)       !parameter used in CNAllocation
  real(rk8):: arootf(0:mxpft)       !parameter used in CNAllocation
  real(rk8):: astemf(0:mxpft)       !parameter used in CNAllocation
  real(rk8):: arooti(0:mxpft)       !parameter used in CNAllocation
  real(rk8):: fleafi(0:mxpft)       !parameter used in CNAllocation
  real(rk8):: allconsl(0:mxpft)     !parameter used in CNAllocation
  real(rk8):: allconss(0:mxpft)     !parameter used in CNAllocation
  real(rk8):: ztopmx(0:mxpft)       !parameter used in CNVegStructUpdate
  real(rk8):: laimx(0:mxpft)        !parameter used in CNVegStructUpdate
  real(rk8):: gddmin(0:mxpft)       !parameter used in CNPhenology
  real(rk8):: hybgdd(0:mxpft)       !parameter used in CNPhenology
  real(rk8):: lfemerg(0:mxpft)      !parameter used in CNPhenology
  real(rk8):: grnfill(0:mxpft)      !parameter used in CNPhenology
  integer :: mxmat(0:mxpft)        !parameter used in CNPhenology
  integer :: mnNHplantdate(0:mxpft)!minimum planting date for NorthHemisphere (YYYYMMDD)
  integer :: mxNHplantdate(0:mxpft)!maximum planting date for NorthHemisphere (YYYYMMDD)
  integer :: mnSHplantdate(0:mxpft)!minimum planting date for SouthHemisphere (YYYYMMDD)
  integer :: mxSHplantdate(0:mxpft)!maximum planting date for SouthHemisphere (YYYYMMDD)
  real(rk8):: planttemp(0:mxpft)    !planting temperature used in CNPhenology (K)
  real(rk8):: minplanttemp(0:mxpft) !mininum planting temperature used in CNPhenology (K)
  real(rk8):: froot_leaf(0:mxpft)   !allocation parameter: new fine root C per new leaf C (gC/gC) 
  real(rk8):: stem_leaf(0:mxpft)    !allocation parameter: new stem c per new leaf C (gC/gC)
  real(rk8):: croot_stem(0:mxpft)   !allocation parameter: new coarse root C per new stem C (gC/gC)
  real(rk8):: flivewd(0:mxpft)      !allocation parameter: fraction of new wood that is live (phloem and ray parenchyma) (no units)
  real(rk8):: fcur(0:mxpft)         !allocation parameter: fraction of allocation that goes to currently displayed growth, remainder to storage
  real(rk8):: fcurdv(0:mxpft)       !alternate fcur for use with cndv
  real(rk8):: lf_flab(0:mxpft)      !leaf litter labile fraction
  real(rk8):: lf_fcel(0:mxpft)      !leaf litter cellulose fraction
  real(rk8):: lf_flig(0:mxpft)      !leaf litter lignin fraction
  real(rk8):: fr_flab(0:mxpft)      !fine root litter labile fraction
  real(rk8):: fr_fcel(0:mxpft)      !fine root litter cellulose fraction
  real(rk8):: fr_flig(0:mxpft)      !fine root litter lignin fraction
  real(rk8):: leaf_long(0:mxpft)    !leaf longevity (yrs)
  real(rk8):: evergreen(0:mxpft)    !binary flag for evergreen leaf habit (0 or 1)
  real(rk8):: stress_decid(0:mxpft) !binary flag for stress-deciduous leaf habit (0 or 1)
  real(rk8):: season_decid(0:mxpft) !binary flag for seasonal-deciduous leaf habit (0 or 1)
  real(rk8):: pconv(0:mxpft)        !proportion of deadstem to conversion flux
  real(rk8):: pprod10(0:mxpft)      !proportion of deadstem to 10-yr product pool
  real(rk8):: pprod100(0:mxpft)     !proportion of deadstem to 100-yr product pool
  real(rk8):: pprodharv10(0:mxpft)  !harvest mortality proportion of deadstem to 10-yr pool
  ! pft paraemeters for fire code
  real(rk8):: cc_leaf(0:mxpft)
  real(rk8):: cc_lstem(0:mxpft)
  real(rk8):: cc_dstem(0:mxpft)
  real(rk8):: cc_other(0:mxpft)
  real(rk8):: fm_leaf(0:mxpft)
  real(rk8):: fm_lstem(0:mxpft)
  real(rk8):: fm_dstem(0:mxpft)
  real(rk8):: fm_other(0:mxpft)
  real(rk8):: fm_root(0:mxpft)
  real(rk8):: fm_lroot(0:mxpft)
  real(rk8):: fm_droot(0:mxpft)
  real(rk8):: fsr_pft(0:mxpft)
  real(rk8):: fd_pft(0:mxpft)
  ! pft parameters for crop code
  real(rk8):: fertnitro(0:mxpft)    !fertilizer
  real(rk8):: fleafcn(0:mxpft)      !C:N during grain fill; leaf
  real(rk8):: ffrootcn(0:mxpft)     !C:N during grain fill; fine root
  real(rk8):: fstemcn(0:mxpft)      !C:N during grain fill; stem

  ! pft parameters for CNDV code
  ! from LPJ subroutine pftparameters
  real(rk8) pftpar20(0:mxpft)       !tree maximum crown area (m2)
  real(rk8) pftpar28(0:mxpft)       !min coldest monthly mean temperature
  real(rk8) pftpar29(0:mxpft)       !max coldest monthly mean temperature
  real(rk8) pftpar30(0:mxpft)       !min growing degree days (>= 5 deg C)
  real(rk8) pftpar31(0:mxpft)       !upper limit of temperature of the warmest month (twmax)
  real(rk8), parameter :: reinickerp = 1.6D0 !parameter in allometric equation
  real(rk8), parameter :: dwood  = 2.5D5   !cn wood density (gC/m3); lpj:2.0e5
  real(rk8), parameter :: allom1 = 100.0D0   !parameters in
  real(rk8), parameter :: allom2 =  40.0D0   !...allometric
  real(rk8), parameter :: allom3 =   0.5D0   !...equations
  real(rk8), parameter :: allom1s = 250.0D0  !modified for shrubs by
  real(rk8), parameter :: allom2s =   8.0D0  !X.D.Z
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: pftconrd ! Read and initialize vegetation (PFT) constants
!
! !REVISION HISTORY:
! Created by Sam Levis (put into module form by Mariana Vertenstein)
! 10/21/03, Peter Thornton: Added new variables for CN code
! 06/24/09, Erik Kluzek: Add indices for all pft types, and add expected_pftnames array and comparision
! 09/17/10, David Lawrence: Modified code to read in netCDF pft physiology file
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: pftconrd
!
! !INTERFACE:
  subroutine pftconrd
!
! !DESCRIPTION:
! Read and initialize vegetation (PFT) constants
!
! !USES:
    use mod_clm_varctl, only : fpftcon
    use mod_clm_varcon, only : tfrz
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
! routine initialize in module initializeMod
!
! !REVISION HISTORY:
! Created by Gordon Bonan
!F. Li and S. Levis (11/06/12)
!
! !LOCAL VARIABLES:
!EOP
    integer :: i,n              ! loop indices
    integer :: ier              ! error code
    integer :: ncid
    integer :: npft             ! number of pfts on pft-physiology file
    logical :: readv            ! read variable in or not
    character(len=32) :: subname = 'pftconrd'              ! subroutine name
    !
    ! Expected PFT names: The names expected on the fpftcon file and the order they are expected to be in.
    ! NOTE: similar types are assumed to be together, first trees (ending with broadleaf_deciduous_boreal_tree
    !       then shrubs, ending with broadleaf_deciduous_boreal_shrub, then grasses starting with c3_arctic_grass
    !       and finally crops, ending with soybean
    ! DO NOT CHANGE THE ORDER -- WITHOUT MODIFYING OTHER PARTS OF THE CODE WHERE THE ORDER MATTERS!
    !
    character(len=40), parameter :: expected_pftnames(0:mxpft) = (/ &
                 'not_vegetated                      '  &
               , 'needleleaf_evergreen_temperate_tree'  &
               , 'needleleaf_evergreen_boreal_tree   '  &
               , 'needleleaf_deciduous_boreal_tree   '  &
               , 'broadleaf_evergreen_tropical_tree  '  &
               , 'broadleaf_evergreen_temperate_tree '  &
               , 'broadleaf_deciduous_tropical_tree  '  &
               , 'broadleaf_deciduous_temperate_tree '  &
               , 'broadleaf_deciduous_boreal_tree    '  &
               , 'broadleaf_evergreen_shrub          '  &
               , 'broadleaf_deciduous_temperate_shrub'  &
               , 'broadleaf_deciduous_boreal_shrub   '  &
               , 'c3_arctic_grass                    '  &
               , 'c3_non-arctic_grass                '  &
               , 'c4_grass                           '  &
               , 'c3_crop                            '  &
               , 'c3_irrigated                       '  &
               , 'corn                               '  &
               , 'irrigated_corn                     '  &
               , 'spring_temperate_cereal            '  &
               , 'irrigated_spring_temperate_cereal  '  &
               , 'winter_temperate_cereal            '  &
               , 'irrigated_winter_temperate_cereal  '  &
               , 'soybean                            '  &
               , 'irrigated_soybean                  '  &
    /)
!-----------------------------------------------------------------------

    ! Set specific vegetation type values

    if (myid == italk) then
       write(stdout,*) 'Attempting to read PFT physiological data .....'
    end if
    call openfile_withname(fpftcon,ncid)
    call ncd_inqdim(ncid,'pft',dlen=npft)
    call read_var1d_static(ncid,'pftname',pftname)
    call read_var1d_static(ncid,'z0mr',npft,z0mr)
    call read_var1d_static(ncid,'displar',npft,displar)
    call read_var1d_static(ncid,'dleaf',npft,dleaf)
    call read_var1d_static(ncid,'c3psn',npft,c3psn)
    call read_var1d_static(ncid,'rholvis',npft,rhol(:,ivis))
    call read_var1d_static(ncid,'rholnir',npft,rhol(:,inir))
    call read_var1d_static(ncid,'rhosvis',npft,rhos(:,ivis))
    call read_var1d_static(ncid,'rhosnir',npft,rhos(:,inir))
    call read_var1d_static(ncid,'taulvis',npft,taul(:,ivis))
    call read_var1d_static(ncid,'taulnir',npft,taul(:,inir))
    call read_var1d_static(ncid,'tausvis',npft,taus(:,ivis))
    call read_var1d_static(ncid,'tausnir',npft,taus(:,inir))
    call read_var1d_static(ncid,'xl',npft,xl)
    call read_var1d_static(ncid,'roota_par',npft,roota_par)
    call read_var1d_static(ncid,'rootb_par',npft,rootb_par)
    call read_var1d_static(ncid,'slatop',npft,slatop)
    call read_var1d_static(ncid,'dsladlai',npft,dsladlai)
    call read_var1d_static(ncid,'leafcn',npft,leafcn)
    call read_var1d_static(ncid,'flnr',npft,flnr)
    call read_var1d_static(ncid,'smpso',npft,smpso)
    call read_var1d_static(ncid,'smpsc',npft,smpsc)
    call read_var1d_static(ncid,'fnitr',npft,fnitr)
    call read_var1d_static(ncid,'woody',npft,woody)
    call read_var1d_static(ncid,'lflitcn',npft,lflitcn)
    call read_var1d_static(ncid,'frootcn',npft,frootcn)
    call read_var1d_static(ncid,'livewdcn',npft,livewdcn)
    call read_var1d_static(ncid,'deadwdcn',npft,deadwdcn)
    call read_var1d_static(ncid,'grperc',npft,grperc)
    call read_var1d_static(ncid,'grpnow',npft,grpnow)
    call read_var1d_static(ncid,'froot_leaf',npft,froot_leaf)
    call read_var1d_static(ncid,'stem_leaf',npft,stem_leaf)
    call read_var1d_static(ncid,'croot_stem',npft,croot_stem)
    call read_var1d_static(ncid,'flivewd',npft,flivewd)
    call read_var1d_static(ncid,'fcur',npft,fcur)
    call read_var1d_static(ncid,'fcurdv',npft,fcurdv)
    call read_var1d_static(ncid,'lf_flab',npft,lf_flab)
    call read_var1d_static(ncid,'lf_fcel',npft,lf_fcel)
    call read_var1d_static(ncid,'lf_flig',npft,lf_flig)
    call read_var1d_static(ncid,'fr_flab',npft,fr_flab)
    call read_var1d_static(ncid,'fr_fcel',npft,fr_fcel)
    call read_var1d_static(ncid,'fr_flig',npft,fr_flig)
    call read_var1d_static(ncid,'leaf_long',npft,leaf_long)
    call read_var1d_static(ncid,'evergreen',npft,evergreen)
    call read_var1d_static(ncid,'stress_decid',npft,stress_decid)
    call read_var1d_static(ncid,'season_decid',npft,season_decid)
    call read_var1d_static(ncid,'pftpar20',npft,pftpar20)
    call read_var1d_static(ncid,'pftpar28',npft,pftpar28)
    call read_var1d_static(ncid,'pftpar29',npft,pftpar29)
    call read_var1d_static(ncid,'pftpar30',npft,pftpar30)
    call read_var1d_static(ncid,'pftpar31',npft,pftpar31)
    call read_var1d_static(ncid,'fertnitro',npft,fertnitro)
    call read_var1d_static(ncid,'fleafcn',npft,fleafcn)
    call read_var1d_static(ncid,'ffrootcn',npft,ffrootcn)
    call read_var1d_static(ncid,'fstemcn',npft,fstemcn)
#ifdef VERTSOILC
    call read_var1d_static(ncid,'rootprof_beta',npft,rootprof_beta)
#endif
    call read_var1d_static(ncid,'pconv',npft,pconv)
    call read_var1d_static(ncid,'pprod10',npft,pprod10)
    call read_var1d_static(ncid,'pprodharv10',npft,pprodharv10)
    call read_var1d_static(ncid,'pprod100',npft,pprod100)
    call read_var1d_static(ncid,'graincn',npft,graincn)
    call read_var1d_static(ncid,'mxtmp',npft,mxtmp)
    call read_var1d_static(ncid,'baset',npft,baset)
    call read_var1d_static(ncid,'declfact',npft,declfact)
    call read_var1d_static(ncid,'bfact',npft,bfact)
    call read_var1d_static(ncid,'aleaff',npft,aleaff)
    call read_var1d_static(ncid,'arootf',npft,arootf)
    call read_var1d_static(ncid,'astemf',npft,astemf)
    call read_var1d_static(ncid,'arooti',npft,arooti)
    call read_var1d_static(ncid,'fleafi',npft,fleafi)
    call read_var1d_static(ncid,'allconsl',npft,allconsl)
    call read_var1d_static(ncid,'allconss',npft,allconss)
    call read_var1d_static(ncid,'crop',npft,crop)
    call read_var1d_static(ncid,'irrigated',npft,irrigated)
    call read_var1d_static(ncid,'ztopmx',npft,ztopmx)
    call read_var1d_static(ncid,'laimx',npft,laimx)
    call read_var1d_static(ncid,'gddmin',npft,gddmin)
    call read_var1d_static(ncid,'hybgdd',npft,hybgdd)
    call read_var1d_static(ncid,'lfemerg',npft,lfemerg)
    call read_var1d_static(ncid,'grnfill',npft,grnfill)
    call read_var1d_static(ncid,'mxmat',npft,mxmat)
    call read_var1d_static(ncid,'cc_leaf',npft,cc_leaf)
    call read_var1d_static(ncid,'cc_lstem',npft,cc_lstem)
    call read_var1d_static(ncid,'cc_dstem',npft,cc_dstem)
    call read_var1d_static(ncid,'cc_other',npft,cc_other)
    call read_var1d_static(ncid,'fm_leaf',npft,fm_leaf)
    call read_var1d_static(ncid,'fm_lstem',npft,fm_lstem)
    call read_var1d_static(ncid,'fm_dstem',npft,fm_dstem)
    call read_var1d_static(ncid,'fm_other',npft,fm_other)
    call read_var1d_static(ncid,'fm_root',npft,fm_root)
    call read_var1d_static(ncid,'fm_lroot',npft,fm_lroot)
    call read_var1d_static(ncid,'fm_droot',npft,fm_droot)
    call read_var1d_static(ncid,'fsr_pft',npft,fsr_pft)
    call read_var1d_static(ncid,'fd_pft',npft,fd_pft)
    call read_var1d_static(ncid,'planting_temp',npft,planttemp)
    call read_var1d_static(ncid,'min_planting_temp',npft,minplanttemp)
    call read_var1d_static(ncid,'min_NH_planting_date',npft,mnNHplantdate)
    call read_var1d_static(ncid,'min_SH_planting_date',npft,mnSHplantdate)
    call read_var1d_static(ncid,'max_NH_planting_date',npft,mxNHplantdate)
    call read_var1d_static(ncid,'max_SH_planting_date',npft,mxSHplantdate)
    call closefile(ncid)

    do i = 0, mxpft
       if ( trim(adjustl(pftname(i))) /= trim(expected_pftnames(i)) )then
          write(stderr,*)'pftconrd: pftname is NOT what is expected, name = ', &
                        trim(pftname(i)), ', expected name = ', trim(expected_pftnames(i))
          call fatal(__FILE__,__LINE__,'pftconrd: bad name for pft on fpftcon dataset' )
       end if
       if ( trim(pftname(i)) == 'not_vegetated'                       ) noveg               = i
       if ( trim(pftname(i)) == 'needleleaf_evergreen_temperate_tree' ) ndllf_evr_tmp_tree  = i
       if ( trim(pftname(i)) == 'needleleaf_evergreen_boreal_tree'    ) ndllf_evr_brl_tree  = i
       if ( trim(pftname(i)) == 'needleleaf_deciduous_boreal_tree'    ) ndllf_dcd_brl_tree  = i
       if ( trim(pftname(i)) == 'broadleaf_evergreen_tropical_tree'   ) nbrdlf_evr_trp_tree  = i
       if ( trim(pftname(i)) == 'broadleaf_evergreen_temperate_tree'  ) nbrdlf_evr_tmp_tree  = i
       if ( trim(pftname(i)) == 'broadleaf_deciduous_tropical_tree'   ) nbrdlf_dcd_trp_tree  = i
       if ( trim(pftname(i)) == 'broadleaf_deciduous_temperate_tree'  ) nbrdlf_dcd_tmp_tree  = i
       if ( trim(pftname(i)) == 'broadleaf_deciduous_boreal_tree'     ) nbrdlf_dcd_brl_tree  = i
       if ( trim(pftname(i)) == 'broadleaf_evergreen_shrub'           ) nbrdlf_evr_shrub     = i
       if ( trim(pftname(i)) == 'broadleaf_deciduous_temperate_shrub' ) nbrdlf_dcd_tmp_shrub = i
       if ( trim(pftname(i)) == 'broadleaf_deciduous_boreal_shrub'    ) nbrdlf_dcd_brl_shrub = i
       if ( trim(pftname(i)) == 'c3_arctic_grass'                     ) nc3_arctic_grass     = i
       if ( trim(pftname(i)) == 'c3_non-arctic_grass'                 ) nc3_nonarctic_grass  = i
       if ( trim(pftname(i)) == 'c4_grass'                            ) nc4_grass            = i
       if ( trim(pftname(i)) == 'c3_crop'                             ) nc3crop              = i
       if ( trim(pftname(i)) == 'c3_irrigated'                        ) nc3irrig               = i
       if ( trim(pftname(i)) == 'corn'                                ) ncorn                = i
       if ( trim(pftname(i)) == 'irrigated_corn'                      ) ncornirrig           = i
       if ( trim(pftname(i)) == 'spring_temperate_cereal'             ) nscereal             = i
       if ( trim(pftname(i)) == 'irrigated_spring_temperate_cereal'   ) nscerealirrig        = i
       if ( trim(pftname(i)) == 'winter_temperate_cereal'             ) nwcereal             = i
       if ( trim(pftname(i)) == 'irrigated_winter_temperate_cereal'   ) nwcerealirrig        = i
       if ( trim(pftname(i)) == 'soybean'                             ) nsoybean             = i
       if ( trim(pftname(i)) == 'irrigated_soybean'                   ) nsoybeanirrig        = i
    end do

    ntree                = nbrdlf_dcd_brl_tree  ! value for last type of tree
    npcropmin            = ncorn                ! first prognostic crop
    npcropmax            = nsoybeanirrig        ! last prognostic crop in list

#if (defined CNDV)
    fcur(:) = fcurdv(:)
#endif

    !
    ! Do some error checking
    !
    if ( npcropmax /= mxpft )then
       call fatal(__FILE__,__LINE__,trim(subname)//' ERROR: npcropmax is NOT the last value' )
    end if
    do i = 0, mxpft
       if (      irrigated(i) == 1.0D0  .and. (i == nc3irrig .or. &
                                                i == ncornirrig .or. &
                                                i == nscerealirrig .or. &
                                                i == nwcerealirrig .or. &
                                                i == nsoybeanirrig) )then
          ! correct
       else if ( irrigated(i) == 0.0D0 )then
          ! correct
       else
          call fatal(__FILE__,__LINE__,trim(subname)//' ERROR: irrigated has wrong values' )
       end if
       if (      crop(i) == 1.0D0 .and. (i >= nc3crop .and. i <= npcropmax) )then
          ! correct
       else if ( crop(i) == 0.0D0 )then
          ! correct
       else
          call fatal(__FILE__,__LINE__,trim(subname)//' ERROR: crop has wrong values' )
       end if
       if ( (i /= noveg) .and. (i < npcropmin) .and. &
            abs(pconv(i)+pprod10(i)+pprod100(i) - 1.0D0) > 1.D-7 )then
          call fatal(__FILE__,__LINE__,trim(subname)//' ERROR: pconv+pprod10+pprod100 do NOT sum to one.' )
       end if
       if ( pprodharv10(i) > 1.0D0 .or. pprodharv10(i) < 0.0D0 )then
          call fatal(__FILE__,__LINE__,trim(subname)//' ERROR: pprodharv10 outside of range.' )
       end if
    end do

    if (myid == italk) then
       write(stdout,*) 'Successfully read PFT physiological data'
       write(stdout,*)
    end if

  end subroutine pftconrd

end module mod_clm_pftvarcon
