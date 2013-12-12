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
  use mod_intkinds
  use mod_realkinds
  use mod_memutil
  use mod_dynparam

  implicit none

  public

  real(rk8) , pointer , dimension(:,:,:) :: rnof , rsubst , rsur , &
    wflux1 , wflux2 , wfluxc , bb , cc , bcoef , deprat , fct2
  real(rk8) , pointer , dimension(:,:,:) :: lfta , lftb , lftra , lftrs
  real(rk8) , pointer , dimension(:,:,:) :: cdrd , vpdc
  real(rk8) , pointer , dimension(:,:,:) :: rppq , efe
  real(rk8) , pointer , dimension(:,:,:) :: dcd , etrc
  real(rk8) , pointer , dimension(:,:,:) :: dels , aseas
  real(rk8) , pointer , dimension(:,:,:) :: rib , swsi
  real(rk8) , pointer , dimension(:,:,:) :: qgrd , qs , resp , rhs
  real(rk8) , pointer , dimension(:,:,:) :: sts , zh
  real(rk8) , pointer , dimension(:,:,:) :: bfc , bsw
  real(rk8) , pointer , dimension(:,:,:) :: evmx0 , fdry , fwet
  real(rk8) , pointer , dimension(:,:,:) :: gwmx0 , gwmx1 , gwmx2
  real(rk8) , pointer , dimension(:,:,:) :: porsl , relfc , relaw , rnet
  real(rk8) , pointer , dimension(:,:,:) :: texrat , vegt , wiltr , wt , xkmx
  real(rk8) , pointer , dimension(:,:,:) :: cdr , cdrn , cf ,  &
    cgrnd , cgrndl , cgrnds , clead , efpr , eg , etr , etrrun , &
    evaps , evapw , fevpg , flnet , flneto , fseng , htvp , ps , &
    pw , qsatl , rhosw , ribd , rlai , rpp , scrat ,      &
    scvk , sdrop , seasb , sigf , sm , tm , uaf , vspda , wata , &
    watr , watt , watu , wta , xlai , xlsai , xrun , z1log ,     &
    z2fra , z10fra , zlgocn , zlglnd ,  zlgsno , zlgveg , zlgdis
  real(rk8) , pointer , dimension(:,:,:) :: cn1 , rgr , wta0 , wtaq0 ,   &
    wtg , wtg0 , wtg2 , wtga , wtgaq , wtgl , wtglq ,  wtgq ,   &
    wtgq0 , wtl0 , wtlh , wtlq , wtlq0 , wtshi , wtsqi , df , fracd
  real(rk8) , pointer , dimension(:,:,:) :: usw , vsw , czenith
  integer(ik4) , pointer , dimension(:,:,:) :: lveg

  contains

  subroutine allocate_mod_bats_internal
    implicit none
    call getmem3d(lfta,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:lfta')
    call getmem3d(lftb,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:lftb')
    call getmem3d(lftra,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:lftra')
    call getmem3d(lftrs,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:lftrs')
    call getmem3d(cdrd,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:cdrd')
    call getmem3d(czenith,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:czenith')
    call getmem3d(vpdc,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:vpdc')
    call getmem3d(rppq,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:rppq')
    call getmem3d(efe,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:efe')
    call getmem3d(dcd,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:dcd')
    call getmem3d(etrc,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:etrc')
    call getmem3d(dels,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:dels')
    call getmem3d(aseas,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:aseas')
    call getmem3d(rib,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:rib')
    call getmem3d(rnof,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:rnof')
    call getmem3d(rsubst,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:rsubst')
    call getmem3d(rsur,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:rsur')
    call getmem3d(wflux1,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:wflux1')
    call getmem3d(wflux2,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:wflux2')
    call getmem3d(wfluxc,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:wfluxc')
    call getmem3d(bb,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:bb')
    call getmem3d(bcoef,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:bcoef')
    call getmem3d(cc,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:cc')
    call getmem3d(deprat,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:deprat')
    call getmem3d(fct2,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:fct2')
    call getmem3d(qgrd,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:qgrd')
    call getmem3d(qs,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:qs')
    call getmem3d(resp,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:resp')
    call getmem3d(rhs,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:rhs')
    call getmem3d(sts,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:sts')
    call getmem3d(zh,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:zh')
    call getmem3d(bfc,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:bfc')
    call getmem3d(bsw,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:bsw')
    call getmem3d(evmx0,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:evmx0')
    call getmem3d(gwmx0,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:gwmx0')
    call getmem3d(gwmx1,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:gwmx1')
    call getmem3d(gwmx2,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:gwmx2')
    call getmem3d(fdry,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:fdry')
    call getmem3d(fwet,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:fwet')
    call getmem3d(porsl,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:porsl')
    call getmem3d(relfc,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:relfc')
    call getmem3d(relaw,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:relaw')
    call getmem3d(rnet,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:rnet')
    call getmem3d(texrat,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:texrat')
    call getmem3d(vegt,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:vegt')
    call getmem3d(wiltr,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:wiltr')
    call getmem3d(wt,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:wt')
    call getmem3d(xkmx,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:xkmx')
    call getmem3d(cdr,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:cdr')
    call getmem3d(cdrn,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:cdrn')
    call getmem3d(cf,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:cf')
    call getmem3d(cgrnd,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:cgrnd')
    call getmem3d(cgrndl,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:cgrndl')
    call getmem3d(cgrnds,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:cgrnds')
    call getmem3d(clead,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:clead')
    call getmem3d(efpr,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:efpr')
    call getmem3d(eg,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:eg')
    call getmem3d(etr,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:etr')
    call getmem3d(etrrun,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:etrrun')
    call getmem3d(evaps,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:evaps')
    call getmem3d(evapw,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:evapw')
    call getmem3d(fevpg,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:fevpg')
    call getmem3d(flnet,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:flnet')
    call getmem3d(flneto,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:flneto')
    call getmem3d(fracd,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:fracd')
    call getmem3d(fseng,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:fseng')
    call getmem3d(htvp,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:htvp')
    call getmem3d(ps,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:ps')
    call getmem3d(pw,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:pw')
    call getmem3d(qsatl,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:qsatl')
    call getmem3d(rhosw,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:rhosw')
    call getmem3d(ribd,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:ribd')
    call getmem3d(rlai,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:rlai')
    call getmem3d(rpp,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:rpp')
    call getmem3d(scrat,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:scrat')
    call getmem3d(scvk,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:scvk')
    call getmem3d(sdrop,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:sdrop')
    call getmem3d(seasb,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:seasb')
    call getmem3d(sigf,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:sigf')
    call getmem3d(sm,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:sm')
    call getmem3d(tm,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:tm')
    call getmem3d(uaf,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:uaf')
    call getmem3d(vspda,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:vspda')
    call getmem3d(wata,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:wata')
    call getmem3d(watr,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:watr')
    call getmem3d(watt,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:watt')
    call getmem3d(watu,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:watu')
    call getmem3d(wta,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:wta')
    call getmem3d(xlai,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:xlai')
    call getmem3d(xlsai,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:xlsai')
    call getmem3d(xrun,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:xrun')
    call getmem3d(cn1,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:cn1')
    call getmem3d(df,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:df')
    call getmem3d(rgr,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:rgr')
    call getmem3d(swsi,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:swsi')
    call getmem3d(usw,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:usw')
    call getmem3d(vsw,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:vsw')
    call getmem3d(wta0,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:wta0')
    call getmem3d(wtaq0,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:wtaq0')
    call getmem3d(wtg,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:wtg')
    call getmem3d(wtg0,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:wtg0')
    call getmem3d(wtg2,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:wtg2')
    call getmem3d(wtga,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:wtga')
    call getmem3d(wtgaq,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:wtgaq')
    call getmem3d(wtgl,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:wtgl')
    call getmem3d(wtglq,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:wtglq')
    call getmem3d(wtgq,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:wtgq')
    call getmem3d(wtgq0,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:wtgq0')
    call getmem3d(wtl0,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:wtl0')
    call getmem3d(wtlh,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:wtlh')
    call getmem3d(wtlq,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:wtlq')
    call getmem3d(wtlq0,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:wtlq0')
    call getmem3d(wtshi,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:wtshi')
    call getmem3d(wtsqi,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:wtsqi')
    call getmem3d(z2fra,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:z2fra')
    call getmem3d(z10fra,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:z10fra')
    call getmem3d(z1log,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:z1log')
    call getmem3d(zlgocn,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:zlgocn')
    call getmem3d(zlglnd,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:zlglnd')
    call getmem3d(zlgsno,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:zlgsno')
    call getmem3d(zlgveg,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:zlgveg')
    call getmem3d(zlgdis,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:zlgdis')
    call getmem3d(lveg,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:lveg')
  end subroutine allocate_mod_bats_internal
!
end module mod_bats_internal
