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

module mod_bats_internal
!
  use mod_realkinds
  use mod_memutil
  use mod_dynparam
!
  private

  real(dp) , pointer , dimension(:,:) :: gwatr , rnof , rsubsr , rsubss , &
        rsubst , rsur , wflux1 , wflux2 , wfluxc , xkmx1 , xkmx2 , xkmxr
  real(dp) , pointer , dimension(:,:) :: sold
  real(dp) , pointer , dimension(:,:) :: bb , bcoef , cc , depann , depdiu , &
        deprat , fct2 , hs , rscsa , rscsd , ska , skd , sks , swtrta , swtrtd
  real(dp) , pointer , dimension(:,:) :: cari
  real(dp) , pointer , dimension(:,:) :: lfta , lftb
  real(dp) , pointer , dimension(:,:) :: lftra , lftrs
  real(dp) , pointer , dimension(:,:) :: cdrd , vpdc
  real(dp) , pointer , dimension(:,:) :: rppq , efe
  real(dp) , pointer , dimension(:,:) :: dcd , etrc
  real(dp) , pointer , dimension(:,:) :: qsatld
  real(dp) , pointer , dimension(:,:) :: dels
  real(dp) , pointer , dimension(:,:) :: efpot , tbef
  real(dp) , pointer , dimension(:,:) :: fsol0 , fsold
  real(dp) , pointer , dimension(:,:) :: radf , rmini
  real(dp) , pointer , dimension(:,:) :: trup , trupd
  real(dp) , pointer , dimension(:,:) :: cdrmin , dlstaf
  real(dp) , pointer , dimension(:,:) :: rib , rib1
  real(dp) , pointer , dimension(:,:) :: p1d0 , qs1d0 , ts1d0
  real(dp) , pointer , dimension(:,:) :: emiss1d , pbp1d , prcp1d
  real(dp) , pointer , dimension(:,:) :: qg1d , qs1d , resp1d , rhs1d
  real(dp) , pointer , dimension(:,:) :: taf1d , ts1d , z1d
  real(dp) , pointer , dimension(:,:) :: bfc , bsw
  real(dp) , pointer , dimension(:,:) :: evmx0 , fdry , fwet
  real(dp) , pointer , dimension(:,:) :: gwmx0 , gwmx1 , gwmx2
  real(dp) , pointer , dimension(:,:) :: porsl , relfc , rnet
  real(dp) , pointer , dimension(:,:) :: texrat , vegt , wiltr , wt , xkmx
  real(dp) , pointer , dimension(:,:) :: cdr , cdrn , cdrx , cf ,     &
         cgrnd , cgrndl , cgrnds , clead , densi , efpr , eg ,       &
         etr , etrrun , evaps , evapw , fevpg , flnet , flneto ,     &
         fseng , htvp , ps , pw , qice , qsatl , rhosw , ribd ,      &
         rlai , rpp , scrat , scvk , sdrop , seasb , sigf , sm ,     &
         tm , uaf , vspda , wata , watr , watt , watu , wta , xlai , &
         xlsai , xrun , z1log , z2fra , z10fra , zlgocn , zlglnd ,   &
         zlgsno , zlgveg , zlgdis
  real(dp) , pointer , dimension(:,:) :: cn1 , rgr , wta0 , wtaq0 ,   &
         wtg , wtg0 , wtg2 , wtga , wtgaq , wtgl , wtglq ,  wtgq ,   &
         wtgq0 , wtl0 , wtlh , wtlq , wtlq0 , wtshi , wtsqi , df
  real(dp) , pointer , dimension(:,:) :: ribl , ribn
  integer , pointer , dimension(:,:) :: ldoc1d , lveg , oveg
!
  public :: gwatr , rnof , rsubsr , rsubss , rsubst , rsur , sold ,     &
            wflux1 , wflux2 , wfluxc , xkmx1 , xkmx2 , xkmxr  , bb ,    &
            bcoef , cc , depann , depdiu , deprat , fct2 , hs , rscsa , &
            rscsd , ska , skd , sks , swtrta , swtrtd , cari
  public :: lfta , lftb , lftra , lftrs , cdrd , vpdc , rppq , efe ,    &
            dcd , etrc , qsatld , dels , efpot , tbef , fsol0 , fsold , &
            radf , rmini , trup , trupd , cdrmin , dlstaf , rib  , rib1
  public :: p1d0 , qs1d0 , ts1d0 , emiss1d , pbp1d , prcp1d , qg1d , qs1d , &
            resp1d , rhs1d , taf1d , ts1d , z1d , bfc , bsw
  public :: evmx0 , fdry , fwet , gwmx0 , gwmx1 , gwmx2 , porsl , relfc ,  &
            rnet , texrat , vegt , wiltr , wt , xkmx , cdr , cdrn , cdrx , &
            cf , cgrnd , cgrndl , cgrnds , clead , densi , efpr , eg ,     &
            etr , etrrun , evaps , evapw , fevpg , flnet , flneto ,        &
            fseng , htvp , ps , pw , qice , qsatl , rhosw , ribd ,         &
            rlai , rpp , scrat , scvk , sdrop , seasb , sigf , sm ,        &
            tm , uaf , vspda , wata , watr , watt , watu , wta , xlai ,    &
            xlsai , xrun , z1log , z2fra , z10fra , zlgocn , zlglnd ,      &
            zlgsno , zlgveg , zlgdis
  public :: cn1 , rgr , wta0 , wtaq0 , wtg , wtg0 , wtg2 , wtga , wtgaq , &
            wtgl , wtglq ,  wtgq , wtgq0 , wtl0 , wtlh , wtlq , wtlq0 ,   &
            wtshi , wtsqi , df , ribl , ribn
  public :: ldoc1d , lveg , oveg

  public :: allocate_mod_bats_internal

  contains

  subroutine allocate_mod_bats_internal
    implicit none
    call getmem2d(lfta,1,nnsg,1,iym1,'bats_internal:lfta')
    call getmem2d(lftb,1,nnsg,1,iym1,'bats_internal:lftb')
    call getmem2d(lftra,1,nnsg,1,iym1,'bats_internal:lftra')
    call getmem2d(lftrs,1,nnsg,1,iym1,'bats_internal:lftrs')
    call getmem2d(cdrd,1,nnsg,1,iym1,'bats_internal:cdrd')
    call getmem2d(vpdc,1,nnsg,1,iym1,'bats_internal:vpdc')
    call getmem2d(rppq,1,nnsg,1,iym1,'bats_internal:rppq')
    call getmem2d(efe,1,nnsg,1,iym1,'bats_internal:efe')
    call getmem2d(qsatld,1,nnsg,1,iym1,'bats_internal:qsatld')
    call getmem2d(dcd,1,nnsg,1,iym1,'bats_internal:dcd')
    call getmem2d(etrc,1,nnsg,1,iym1,'bats_internal:etrc')
    call getmem2d(dels,1,nnsg,1,iym1,'bats_internal:dels')
    call getmem2d(efpot,1,nnsg,1,iym1,'bats_internal:efpot')
    call getmem2d(tbef,1,nnsg,1,iym1,'bats_internal:tbef')
    call getmem2d(fsol0,1,nnsg,1,iym1,'bats_internal:fsol0')
    call getmem2d(fsold,1,nnsg,1,iym1,'bats_internal:fsold')
    call getmem2d(radf,1,nnsg,1,iym1,'bats_internal:radf')
    call getmem2d(rmini,1,nnsg,1,iym1,'bats_internal:rmini')
    call getmem2d(trup,1,nnsg,1,iym1,'bats_internal:trup')
    call getmem2d(trupd,1,nnsg,1,iym1,'bats_internal:trupd')
    call getmem2d(cdrmin,1,nnsg,1,iym1,'bats_internal:cdrmin')
    call getmem2d(dlstaf,1,nnsg,1,iym1,'bats_internal:dlstaf')
    call getmem2d(rib,1,nnsg,1,iym1,'bats_internal:rib')
    call getmem2d(rib1,1,nnsg,1,iym1,'bats_internal:rib1')
    call getmem2d(gwatr,1,nnsg,1,iym1,'bats_internal:gwatr')
    call getmem2d(rnof,1,nnsg,1,iym1,'bats_internal:rnof')
    call getmem2d(rsubsr,1,nnsg,1,iym1,'bats_internal:rsubsr')
    call getmem2d(rsubss,1,nnsg,1,iym1,'bats_internal:rsubss')
    call getmem2d(rsubst,1,nnsg,1,iym1,'bats_internal:rsubst')
    call getmem2d(rsur,1,nnsg,1,iym1,'bats_internal:rsur')
    call getmem2d(wflux1,1,nnsg,1,iym1,'bats_internal:wflux1')
    call getmem2d(wflux2,1,nnsg,1,iym1,'bats_internal:wflux2')
    call getmem2d(wfluxc,1,nnsg,1,iym1,'bats_internal:wfluxc')
    call getmem2d(xkmx1,1,nnsg,1,iym1,'bats_internal:xkmx1')
    call getmem2d(xkmx2,1,nnsg,1,iym1,'bats_internal:xkmx2')
    call getmem2d(xkmxr,1,nnsg,1,iym1,'bats_internal:xkmxr')
    call getmem2d(sold,1,nnsg,1,iym1,'bats_internal:sold')
    call getmem2d(bb,1,nnsg,1,iym1,'bats_internal:bb')
    call getmem2d(bcoef,1,nnsg,1,iym1,'bats_internal:bcoef')
    call getmem2d(cc,1,nnsg,1,iym1,'bats_internal:cc')
    call getmem2d(depann,1,nnsg,1,iym1,'bats_internal:depann')
    call getmem2d(depdiu,1,nnsg,1,iym1,'bats_internal:depdiu')
    call getmem2d(deprat,1,nnsg,1,iym1,'bats_internal:deprat')
    call getmem2d(fct2,1,nnsg,1,iym1,'bats_internal:fct2')
    call getmem2d(hs,1,nnsg,1,iym1,'bats_internal:hs')
    call getmem2d(rscsa,1,nnsg,1,iym1,'bats_internal:rscsa')
    call getmem2d(rscsd,1,nnsg,1,iym1,'bats_internal:rscsd')
    call getmem2d(ska,1,nnsg,1,iym1,'bats_internal:ska')
    call getmem2d(skd,1,nnsg,1,iym1,'bats_internal:skd')
    call getmem2d(sks,1,nnsg,1,iym1,'bats_internal:sks')
    call getmem2d(swtrta,1,nnsg,1,iym1,'bats_internal:swtrta')
    call getmem2d(swtrtd,1,nnsg,1,iym1,'bats_internal:swtrtd')
    call getmem2d(cari,1,nnsg,1,iym1,'bats_internal:cari')
    call getmem2d(p1d0,1,nnsg,1,iym1,'bats_internal:p1d0')
    call getmem2d(qs1d0,1,nnsg,1,iym1,'bats_internal:qs1d0')
    call getmem2d(ts1d0,1,nnsg,1,iym1,'bats_internal:ts1d0')
    call getmem2d(emiss1d,1,nnsg,1,iym1,'bats_internal:emiss1d')
    call getmem2d(pbp1d,1,nnsg,1,iym1,'bats_internal:pbp1d')
    call getmem2d(prcp1d,1,nnsg,1,iym1,'bats_internal:prcp1d')
    call getmem2d(qg1d,1,nnsg,1,iym1,'bats_internal:qg1d')
    call getmem2d(qs1d,1,nnsg,1,iym1,'bats_internal:qs1d')
    call getmem2d(resp1d,1,nnsg,1,iym1,'bats_internal:resp1d')
    call getmem2d(rhs1d,1,nnsg,1,iym1,'bats_internal:rhs1d')
    call getmem2d(taf1d,1,nnsg,1,iym1,'bats_internal:taf1d')
    call getmem2d(ts1d,1,nnsg,1,iym1,'bats_internal:ts1d')
    call getmem2d(z1d,1,nnsg,1,iym1,'bats_internal:z1d')
    call getmem2d(bfc,1,nnsg,1,iym1,'bats_internal:bfc')
    call getmem2d(bsw,1,nnsg,1,iym1,'bats_internal:bsw')
    call getmem2d(evmx0,1,nnsg,1,iym1,'bats_internal:evmx0')
    call getmem2d(gwmx0,1,nnsg,1,iym1,'bats_internal:gwmx0')
    call getmem2d(gwmx1,1,nnsg,1,iym1,'bats_internal:gwmx1')
    call getmem2d(gwmx2,1,nnsg,1,iym1,'bats_internal:gwmx2')
    call getmem2d(fdry,1,nnsg,1,iym1,'bats_internal:fdry')
    call getmem2d(fwet,1,nnsg,1,iym1,'bats_internal:fwet')
    call getmem2d(porsl,1,nnsg,1,iym1,'bats_internal:porsl')
    call getmem2d(relfc,1,nnsg,1,iym1,'bats_internal:relfc')
    call getmem2d(rnet,1,nnsg,1,iym1,'bats_internal:rnet')
    call getmem2d(texrat,1,nnsg,1,iym1,'bats_internal:texrat')
    call getmem2d(vegt,1,nnsg,1,iym1,'bats_internal:vegt')
    call getmem2d(wiltr,1,nnsg,1,iym1,'bats_internal:wiltr')
    call getmem2d(wt,1,nnsg,1,iym1,'bats_internal:wt')
    call getmem2d(xkmx,1,nnsg,1,iym1,'bats_internal:xkmx')
    call getmem2d(cdr,1,nnsg,1,iym1,'bats_internal:cdr')
    call getmem2d(cdrn,1,nnsg,1,iym1,'bats_internal:cdrn')
    call getmem2d(cdrx,1,nnsg,1,iym1,'bats_internal:cdrx')
    call getmem2d(cf,1,nnsg,1,iym1,'bats_internal:cf')
    call getmem2d(cgrnd,1,nnsg,1,iym1,'bats_internal:cgrnd')
    call getmem2d(cgrndl,1,nnsg,1,iym1,'bats_internal:cgrndl')
    call getmem2d(cgrnds,1,nnsg,1,iym1,'bats_internal:cgrnds')
    call getmem2d(clead,1,nnsg,1,iym1,'bats_internal:clead')
    call getmem2d(densi,1,nnsg,1,iym1,'bats_internal:densi')
    call getmem2d(efpr,1,nnsg,1,iym1,'bats_internal:efpr')
    call getmem2d(eg,1,nnsg,1,iym1,'bats_internal:eg')
    call getmem2d(etr,1,nnsg,1,iym1,'bats_internal:etr')
    call getmem2d(etrrun,1,nnsg,1,iym1,'bats_internal:etrrun')
    call getmem2d(evaps,1,nnsg,1,iym1,'bats_internal:evaps')
    call getmem2d(evapw,1,nnsg,1,iym1,'bats_internal:evapw')
    call getmem2d(fevpg,1,nnsg,1,iym1,'bats_internal:fevpg')
    call getmem2d(flnet,1,nnsg,1,iym1,'bats_internal:flnet')
    call getmem2d(flneto,1,nnsg,1,iym1,'bats_internal:flneto')
    call getmem2d(fseng,1,nnsg,1,iym1,'bats_internal:fseng')
    call getmem2d(htvp,1,nnsg,1,iym1,'bats_internal:htvp')
    call getmem2d(ps,1,nnsg,1,iym1,'bats_internal:ps')
    call getmem2d(pw,1,nnsg,1,iym1,'bats_internal:pw')
    call getmem2d(qice,1,nnsg,1,iym1,'bats_internal:qice')
    call getmem2d(qsatl,1,nnsg,1,iym1,'bats_internal:qsatl')
    call getmem2d(rhosw,1,nnsg,1,iym1,'bats_internal:rhosw')
    call getmem2d(ribd,1,nnsg,1,iym1,'bats_internal:ribd')
    call getmem2d(rlai,1,nnsg,1,iym1,'bats_internal:rlai')
    call getmem2d(rpp,1,nnsg,1,iym1,'bats_internal:rpp')
    call getmem2d(scrat,1,nnsg,1,iym1,'bats_internal:scrat')
    call getmem2d(scvk,1,nnsg,1,iym1,'bats_internal:scvk')
    call getmem2d(sdrop,1,nnsg,1,iym1,'bats_internal:sdrop')
    call getmem2d(seasb,1,nnsg,1,iym1,'bats_internal:seasb')
    call getmem2d(sigf,1,nnsg,1,iym1,'bats_internal:sigf')
    call getmem2d(sm,1,nnsg,1,iym1,'bats_internal:sm')
    call getmem2d(tm,1,nnsg,1,iym1,'bats_internal:tm')
    call getmem2d(uaf,1,nnsg,1,iym1,'bats_internal:uaf')
    call getmem2d(vspda,1,nnsg,1,iym1,'bats_internal:vspda')
    call getmem2d(wata,1,nnsg,1,iym1,'bats_internal:wata')
    call getmem2d(watr,1,nnsg,1,iym1,'bats_internal:watr')
    call getmem2d(watt,1,nnsg,1,iym1,'bats_internal:watt')
    call getmem2d(watu,1,nnsg,1,iym1,'bats_internal:watu')
    call getmem2d(wta,1,nnsg,1,iym1,'bats_internal:wta')
    call getmem2d(xlai,1,nnsg,1,iym1,'bats_internal:xlai')
    call getmem2d(xlsai,1,nnsg,1,iym1,'bats_internal:xlsai')
    call getmem2d(xrun,1,nnsg,1,iym1,'bats_internal:xrun')
    call getmem2d(cn1,1,nnsg,1,iym1,'bats_internal:cn1')
    call getmem2d(df,1,nnsg,1,iym1,'bats_internal:df')
    call getmem2d(rgr,1,nnsg,1,iym1,'bats_internal:rgr')
    call getmem2d(wta0,1,nnsg,1,iym1,'bats_internal:wta0')
    call getmem2d(wtaq0,1,nnsg,1,iym1,'bats_internal:wtaq0')
    call getmem2d(wtg,1,nnsg,1,iym1,'bats_internal:wtg')
    call getmem2d(wtg0,1,nnsg,1,iym1,'bats_internal:wtg0')
    call getmem2d(wtg2,1,nnsg,1,iym1,'bats_internal:wtg2')
    call getmem2d(wtga,1,nnsg,1,iym1,'bats_internal:wtga')
    call getmem2d(wtgaq,1,nnsg,1,iym1,'bats_internal:wtgaq')
    call getmem2d(wtgl,1,nnsg,1,iym1,'bats_internal:wtgl')
    call getmem2d(wtglq,1,nnsg,1,iym1,'bats_internal:wtglq')
    call getmem2d(wtgq,1,nnsg,1,iym1,'bats_internal:wtgq')
    call getmem2d(wtgq0,1,nnsg,1,iym1,'bats_internal:wtgq0')
    call getmem2d(wtl0,1,nnsg,1,iym1,'bats_internal:wtl0')
    call getmem2d(wtlh,1,nnsg,1,iym1,'bats_internal:wtlh')
    call getmem2d(wtlq,1,nnsg,1,iym1,'bats_internal:wtlq')
    call getmem2d(wtlq0,1,nnsg,1,iym1,'bats_internal:wtlq0')
    call getmem2d(wtshi,1,nnsg,1,iym1,'bats_internal:wtshi')
    call getmem2d(wtsqi,1,nnsg,1,iym1,'bats_internal:wtsqi')
    call getmem2d(z2fra,1,nnsg,1,iym1,'bats_internal:z2fra')
    call getmem2d(z10fra,1,nnsg,1,iym1,'bats_internal:z10fra')
    call getmem2d(z1log,1,nnsg,1,iym1,'bats_internal:z1log')
    call getmem2d(zlgocn,1,nnsg,1,iym1,'bats_internal:zlgocn')
    call getmem2d(zlglnd,1,nnsg,1,iym1,'bats_internal:zlglnd')
    call getmem2d(zlgsno,1,nnsg,1,iym1,'bats_internal:zlgsno')
    call getmem2d(zlgveg,1,nnsg,1,iym1,'bats_internal:zlgveg')
    call getmem2d(zlgdis,1,nnsg,1,iym1,'bats_internal:zlgdis')
    call getmem2d(ribl,1,nnsg,1,iym1,'bats_internal:ribl')
    call getmem2d(ribn,1,nnsg,1,iym1,'bats_internal:ribn')

    call getmem2d(ldoc1d,1,nnsg,1,iym1,'bats_internal:ldoc1d')
    call getmem2d(lveg,1,nnsg,1,iym1,'bats_internal:lveg')
    call getmem2d(oveg,1,nnsg,1,iym1,'bats_internal:oveg')

!
  end subroutine allocate_mod_bats_internal
!
end module mod_bats_internal
