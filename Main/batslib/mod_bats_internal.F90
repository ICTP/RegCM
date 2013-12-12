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

  real(rk8) , pointer , dimension(:,:,:) :: delq
  real(rk8) , pointer , dimension(:,:,:) :: rnof
  real(rk8) , pointer , dimension(:,:,:) :: rsubst
  real(rk8) , pointer , dimension(:,:,:) :: rsur
  real(rk8) , pointer , dimension(:,:,:) :: wflux1
  real(rk8) , pointer , dimension(:,:,:) :: wflux2
  real(rk8) , pointer , dimension(:,:,:) :: wfluxc
  real(rk8) , pointer , dimension(:,:,:) :: bb
  real(rk8) , pointer , dimension(:,:,:) :: cc
  real(rk8) , pointer , dimension(:,:,:) :: bcoef
  real(rk8) , pointer , dimension(:,:,:) :: deprat
  real(rk8) , pointer , dimension(:,:,:) :: fct2
  real(rk8) , pointer , dimension(:,:,:) :: lfta
  real(rk8) , pointer , dimension(:,:,:) :: lftb
  real(rk8) , pointer , dimension(:,:,:) :: lftra
  real(rk8) , pointer , dimension(:,:,:) :: lftrs
  real(rk8) , pointer , dimension(:,:,:) :: cdrd
  real(rk8) , pointer , dimension(:,:,:) :: vpdc
  real(rk8) , pointer , dimension(:,:,:) :: rppq
  real(rk8) , pointer , dimension(:,:,:) :: efe
  real(rk8) , pointer , dimension(:,:,:) :: dcd
  real(rk8) , pointer , dimension(:,:,:) :: etrc
  real(rk8) , pointer , dimension(:,:,:) :: dels
  real(rk8) , pointer , dimension(:,:,:) :: aseas
  real(rk8) , pointer , dimension(:,:,:) :: rib
  real(rk8) , pointer , dimension(:,:,:) :: swsi
  real(rk8) , pointer , dimension(:,:,:) :: qgrd
  real(rk8) , pointer , dimension(:,:,:) :: qs
  real(rk8) , pointer , dimension(:,:,:) :: resp
  real(rk8) , pointer , dimension(:,:,:) :: rhs
  real(rk8) , pointer , dimension(:,:,:) :: sts
  real(rk8) , pointer , dimension(:,:,:) :: zh
  real(rk8) , pointer , dimension(:,:,:) :: bfc
  real(rk8) , pointer , dimension(:,:,:) :: bsw
  real(rk8) , pointer , dimension(:,:,:) :: evmx0
  real(rk8) , pointer , dimension(:,:,:) :: fdry
  real(rk8) , pointer , dimension(:,:,:) :: fwet
  real(rk8) , pointer , dimension(:,:,:) :: gwmx0
  real(rk8) , pointer , dimension(:,:,:) :: gwmx1
  real(rk8) , pointer , dimension(:,:,:) :: gwmx2
  real(rk8) , pointer , dimension(:,:,:) :: porsl
  real(rk8) , pointer , dimension(:,:,:) :: relfc
  real(rk8) , pointer , dimension(:,:,:) :: relaw
  real(rk8) , pointer , dimension(:,:,:) :: rnet
  real(rk8) , pointer , dimension(:,:,:) :: texrat
  real(rk8) , pointer , dimension(:,:,:) :: vegt
  real(rk8) , pointer , dimension(:,:,:) :: wiltr
  real(rk8) , pointer , dimension(:,:,:) :: wt
  real(rk8) , pointer , dimension(:,:,:) :: xkmx
  real(rk8) , pointer , dimension(:,:,:) :: cdr
  real(rk8) , pointer , dimension(:,:,:) :: cdrn
  real(rk8) , pointer , dimension(:,:,:) :: cf
  real(rk8) , pointer , dimension(:,:,:) :: cgrnd
  real(rk8) , pointer , dimension(:,:,:) :: cgrndl
  real(rk8) , pointer , dimension(:,:,:) :: cgrnds
  real(rk8) , pointer , dimension(:,:,:) :: clead
  real(rk8) , pointer , dimension(:,:,:) :: efpr
  real(rk8) , pointer , dimension(:,:,:) :: eg
  real(rk8) , pointer , dimension(:,:,:) :: etr
  real(rk8) , pointer , dimension(:,:,:) :: etrrun
  real(rk8) , pointer , dimension(:,:,:) :: evaps
  real(rk8) , pointer , dimension(:,:,:) :: evapw
  real(rk8) , pointer , dimension(:,:,:) :: fevpg
  real(rk8) , pointer , dimension(:,:,:) :: flnet
  real(rk8) , pointer , dimension(:,:,:) :: flneto
  real(rk8) , pointer , dimension(:,:,:) :: fseng
  real(rk8) , pointer , dimension(:,:,:) :: htvp
  real(rk8) , pointer , dimension(:,:,:) :: ps
  real(rk8) , pointer , dimension(:,:,:) :: pw
  real(rk8) , pointer , dimension(:,:,:) :: qsatl
  real(rk8) , pointer , dimension(:,:,:) :: rhosw
  real(rk8) , pointer , dimension(:,:,:) :: ribd
  real(rk8) , pointer , dimension(:,:,:) :: rlai
  real(rk8) , pointer , dimension(:,:,:) :: rpp
  real(rk8) , pointer , dimension(:,:,:) :: scrat
  real(rk8) , pointer , dimension(:,:,:) :: scvk
  real(rk8) , pointer , dimension(:,:,:) :: sdrop
  real(rk8) , pointer , dimension(:,:,:) :: seasb
  real(rk8) , pointer , dimension(:,:,:) :: sigf
  real(rk8) , pointer , dimension(:,:,:) :: sm
  real(rk8) , pointer , dimension(:,:,:) :: tm
  real(rk8) , pointer , dimension(:,:,:) :: uaf
  real(rk8) , pointer , dimension(:,:,:) :: vspda
  real(rk8) , pointer , dimension(:,:,:) :: wata
  real(rk8) , pointer , dimension(:,:,:) :: watr
  real(rk8) , pointer , dimension(:,:,:) :: watt
  real(rk8) , pointer , dimension(:,:,:) :: watu
  real(rk8) , pointer , dimension(:,:,:) :: wta
  real(rk8) , pointer , dimension(:,:,:) :: xlai
  real(rk8) , pointer , dimension(:,:,:) :: xlsai
  real(rk8) , pointer , dimension(:,:,:) :: xrun
  real(rk8) , pointer , dimension(:,:,:) :: z1log
  real(rk8) , pointer , dimension(:,:,:) :: z2fra
  real(rk8) , pointer , dimension(:,:,:) :: z10fra
  real(rk8) , pointer , dimension(:,:,:) :: zlgocn
  real(rk8) , pointer , dimension(:,:,:) :: zlglnd
  real(rk8) , pointer , dimension(:,:,:) :: zlgsno
  real(rk8) , pointer , dimension(:,:,:) :: zlgveg
  real(rk8) , pointer , dimension(:,:,:) :: zlgdis
  real(rk8) , pointer , dimension(:,:,:) :: cn1
  real(rk8) , pointer , dimension(:,:,:) :: rgr
  real(rk8) , pointer , dimension(:,:,:) :: wta0
  real(rk8) , pointer , dimension(:,:,:) :: wtaq0
  real(rk8) , pointer , dimension(:,:,:) :: wtg
  real(rk8) , pointer , dimension(:,:,:) :: wtg0
  real(rk8) , pointer , dimension(:,:,:) :: wtg2
  real(rk8) , pointer , dimension(:,:,:) :: wtga
  real(rk8) , pointer , dimension(:,:,:) :: wtgaq
  real(rk8) , pointer , dimension(:,:,:) :: wtgl
  real(rk8) , pointer , dimension(:,:,:) :: wtglq
  real(rk8) , pointer , dimension(:,:,:) ::  wtgq
  real(rk8) , pointer , dimension(:,:,:) :: wtgq0
  real(rk8) , pointer , dimension(:,:,:) :: wtl0
  real(rk8) , pointer , dimension(:,:,:) :: wtlh
  real(rk8) , pointer , dimension(:,:,:) :: wtlq
  real(rk8) , pointer , dimension(:,:,:) :: wtlq0
  real(rk8) , pointer , dimension(:,:,:) :: wtshi
  real(rk8) , pointer , dimension(:,:,:) :: wtsqi
  real(rk8) , pointer , dimension(:,:,:) :: df
  real(rk8) , pointer , dimension(:,:,:) :: fracd
  real(rk8) , pointer , dimension(:,:,:) :: usw
  real(rk8) , pointer , dimension(:,:,:) :: vsw
  real(rk8) , pointer , dimension(:,:,:) :: czenith
  real(rk8) , pointer , dimension(:,:,:) :: swflx
  real(rk8) , pointer , dimension(:,:,:) :: lwflx
  real(rk8) , pointer , dimension(:,:,:) :: abswveg
  integer(ik4) , pointer , dimension(:,:,:) :: lveg

  contains

  subroutine allocate_mod_bats_internal
    implicit none
    call getmem3d(delq,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:delq')
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
    call getmem3d(swflx,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:swflx')
    call getmem3d(lwflx,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:lwflx')
    call getmem3d(abswveg,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:abswveg')
    call getmem3d(lveg,1,nnsg,jci1,jci2,ici1,ici2,'bats_internal:lveg')
  end subroutine allocate_mod_bats_internal
!
end module mod_bats_internal
