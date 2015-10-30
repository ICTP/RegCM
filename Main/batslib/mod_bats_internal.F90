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
  use mod_regcm_types

  implicit none

  public

  real(rk8) :: dtbat

  integer(ik4) :: nlandp
  integer(ik4) :: ilndbeg , ilndend

  real(rk8) , pointer , dimension(:) :: abswveg
  real(rk8) , pointer , dimension(:) :: aseas
  real(rk8) , pointer , dimension(:) :: bseas
  real(rk8) , pointer , dimension(:) :: bb
  real(rk8) , pointer , dimension(:) :: bcoef
  real(rk8) , pointer , dimension(:) :: bfc
  real(rk8) , pointer , dimension(:) :: bsw
  real(rk8) , pointer , dimension(:) :: cc
  real(rk8) , pointer , dimension(:) :: cdr
  real(rk8) , pointer , dimension(:) :: cdrd
  real(rk8) , pointer , dimension(:) :: cdrn
  real(rk8) , pointer , dimension(:) :: cdrx
  real(rk8) , pointer , dimension(:) :: cf
  real(rk8) , pointer , dimension(:) :: cgrnd
  real(rk8) , pointer , dimension(:) :: cgrndl
  real(rk8) , pointer , dimension(:) :: cgrnds
  real(rk8) , pointer , dimension(:) :: cn1
  real(rk8) , pointer , dimension(:) :: czenith
  real(rk8) , pointer , dimension(:) :: dcd
  real(rk8) , pointer , dimension(:) :: delq
  real(rk8) , pointer , dimension(:) :: dels
  real(rk8) , pointer , dimension(:) :: delt
  real(rk8) , pointer , dimension(:) :: deprat
  real(rk8) , pointer , dimension(:) :: df
  real(rk8) , pointer , dimension(:) :: drag
  real(rk8) , pointer , dimension(:) :: dzh
  real(rk8) , pointer , dimension(:) :: efe
  real(rk8) , pointer , dimension(:) :: efpr
  real(rk8) , pointer , dimension(:) :: eg
  real(rk8) , pointer , dimension(:) :: emiss
  real(rk8) , pointer , dimension(:) :: etr
  real(rk8) , pointer , dimension(:) :: etrc
  real(rk8) , pointer , dimension(:) :: etrrun
  real(rk8) , pointer , dimension(:) :: evaps
  real(rk8) , pointer , dimension(:) :: evapw
  real(rk8) , pointer , dimension(:) :: evmx0
  real(rk8) , pointer , dimension(:) :: evpr
  real(rk8) , pointer , dimension(:) :: fact10
  real(rk8) , pointer , dimension(:) :: fact2
  real(rk8) , pointer , dimension(:) :: fct2
  real(rk8) , pointer , dimension(:) :: fdry
  real(rk8) , pointer , dimension(:) :: fevpg
  real(rk8) , pointer , dimension(:) :: flnet
  real(rk8) , pointer , dimension(:) :: flneto
  real(rk8) , pointer , dimension(:) :: fracd
  real(rk8) , pointer , dimension(:) :: fseng
  real(rk8) , pointer , dimension(:) :: fwet
  real(rk8) , pointer , dimension(:) :: gwet
  real(rk8) , pointer , dimension(:) :: gwmx0
  real(rk8) , pointer , dimension(:) :: gwmx1
  real(rk8) , pointer , dimension(:) :: gwmx2
  real(rk8) , pointer , dimension(:) :: ht
  real(rk8) , pointer , dimension(:) :: hts
  real(rk8) , pointer , dimension(:) :: htvp
  real(rk8) , pointer , dimension(:) :: lat
  real(rk8) , pointer , dimension(:) :: ldew
  real(rk8) , pointer , dimension(:) :: lftra
  real(rk8) , pointer , dimension(:) :: lftrs
  real(rk8) , pointer , dimension(:) :: lncl
  real(rk8) , pointer , dimension(:) :: lwal
  real(rk8) , pointer , dimension(:) :: lwdifal
  real(rk8) , pointer , dimension(:) :: lwdiral
  real(rk8) , pointer , dimension(:) :: lwflx
  real(rk8) , pointer , dimension(:) :: porsl
  real(rk8) , pointer , dimension(:) :: prcp
  real(rk8) , pointer , dimension(:) :: ps
  real(rk8) , pointer , dimension(:) :: pw
  real(rk8) , pointer , dimension(:) :: qgrd
  real(rk8) , pointer , dimension(:) :: qs
  real(rk8) , pointer , dimension(:) :: qsatl
  ! real(rk8) , pointer , dimension(:) :: relaw
  real(rk8) , pointer , dimension(:) :: relfc
  real(rk8) , pointer , dimension(:) :: resp
  real(rk8) , pointer , dimension(:) :: rgr
  real(rk8) , pointer , dimension(:) :: rhosw
  real(rk8) , pointer , dimension(:) :: rhs
  real(rk8) , pointer , dimension(:) :: rib
  real(rk8) , pointer , dimension(:) :: ribd
  real(rk8) , pointer , dimension(:) :: rlai
  real(rk8) , pointer , dimension(:) :: rnet
  real(rk8) , pointer , dimension(:) :: rnof
  real(rk8) , pointer , dimension(:) :: rpp
  real(rk8) , pointer , dimension(:) :: rppq
  real(rk8) , pointer , dimension(:) :: rsubst
  real(rk8) , pointer , dimension(:) :: rsur
  real(rk8) , pointer , dimension(:) :: rsw
  real(rk8) , pointer , dimension(:) :: scrat
  real(rk8) , pointer , dimension(:) :: scvk
  real(rk8) , pointer , dimension(:) :: sdrop
  real(rk8) , pointer , dimension(:) :: sent
  real(rk8) , pointer , dimension(:) :: sfcp
  real(rk8) , pointer , dimension(:) :: sigf
  real(rk8) , pointer , dimension(:) :: sm
  real(rk8) , pointer , dimension(:) :: snag
  real(rk8) , pointer , dimension(:) :: sncv
  real(rk8) , pointer , dimension(:) :: srnof
  real(rk8) , pointer , dimension(:) :: ssw
  real(rk8) , pointer , dimension(:) :: sts
  real(rk8) , pointer , dimension(:) :: swal
  real(rk8) , pointer , dimension(:) :: swdifal
  real(rk8) , pointer , dimension(:) :: swdiral
  real(rk8) , pointer , dimension(:) :: swflx
  real(rk8) , pointer , dimension(:) :: swsi
  real(rk8) , pointer , dimension(:) :: taf
  real(rk8) , pointer , dimension(:) :: texrat
  real(rk8) , pointer , dimension(:) :: tgbb
  real(rk8) , pointer , dimension(:) :: tgbrd
  real(rk8) , pointer , dimension(:) :: tgrd
  real(rk8) , pointer , dimension(:) :: tlef
  real(rk8) , pointer , dimension(:) :: tm
  real(rk8) , pointer , dimension(:) :: trnof
  real(rk8) , pointer , dimension(:) :: tsw
  real(rk8) , pointer , dimension(:) :: uaf
  real(rk8) , pointer , dimension(:) :: usw
  real(rk8) , pointer , dimension(:) :: vegt
  real(rk8) , pointer , dimension(:) :: vpdc
  real(rk8) , pointer , dimension(:) :: vspda
  real(rk8) , pointer , dimension(:) :: vsw
  real(rk8) , pointer , dimension(:) :: wata
  real(rk8) , pointer , dimension(:) :: watr
  real(rk8) , pointer , dimension(:) :: watt
  real(rk8) , pointer , dimension(:) :: watu
  real(rk8) , pointer , dimension(:) :: wflux1
  real(rk8) , pointer , dimension(:) :: wflux2
  real(rk8) , pointer , dimension(:) :: wfluxc
  real(rk8) , pointer , dimension(:) :: wiltr
  real(rk8) , pointer , dimension(:) :: wt
  real(rk8) , pointer , dimension(:) :: wta
  real(rk8) , pointer , dimension(:) :: wta0
  real(rk8) , pointer , dimension(:) :: wtaq0
  real(rk8) , pointer , dimension(:) :: wtg
  real(rk8) , pointer , dimension(:) :: wtg0
  real(rk8) , pointer , dimension(:) :: wtg2
  real(rk8) , pointer , dimension(:) :: wtga
  real(rk8) , pointer , dimension(:) :: wtgaq
  real(rk8) , pointer , dimension(:) :: wtgl
  real(rk8) , pointer , dimension(:) :: wtglq
  real(rk8) , pointer , dimension(:) :: wtgq
  real(rk8) , pointer , dimension(:) :: wtgq0
  real(rk8) , pointer , dimension(:) :: wtl0
  real(rk8) , pointer , dimension(:) :: wtlh
  real(rk8) , pointer , dimension(:) :: wtlq
  real(rk8) , pointer , dimension(:) :: wtlq0
  real(rk8) , pointer , dimension(:) :: wtshi
  real(rk8) , pointer , dimension(:) :: wtsqi
  real(rk8) , pointer , dimension(:) :: xkmx
  real(rk8) , pointer , dimension(:) :: xlai
  real(rk8) , pointer , dimension(:) :: xlsai
  real(rk8) , pointer , dimension(:) :: xrun
  real(rk8) , pointer , dimension(:) :: z10fra
  real(rk8) , pointer , dimension(:) :: z1log
  real(rk8) , pointer , dimension(:) :: z2fra
  real(rk8) , pointer , dimension(:) :: zh
  real(rk8) , pointer , dimension(:) :: zlgdis
  real(rk8) , pointer , dimension(:) :: zlglnd
  real(rk8) , pointer , dimension(:) :: zlgsno
  real(rk8) , pointer , dimension(:) :: zlgveg

  integer(ik4) , pointer , dimension(:) :: lveg

  real(rk8) , pointer , dimension(:) :: p0
  real(rk8) , pointer , dimension(:) :: qs0
  real(rk8) , pointer , dimension(:) :: ts0
  real(rk8) , pointer , dimension(:) :: swd0
  real(rk8) , pointer , dimension(:) :: swf0
  real(rk8) , pointer , dimension(:) :: ncpr0
  real(rk8) , pointer , dimension(:) :: cpr0

  contains

  subroutine allocate_mod_bats_internal(cl)
    implicit none
    type (masked_comm) , intent(in) :: cl
    nlandp = cl%linear_npoint_sg(myid+1)
    ilndbeg = 1
    ilndend = nlandp
    call getmem1d(abswveg,1,nlandp,'bats_internal:abswveg')
    call getmem1d(aseas,1,nlandp,'bats_internal:aseas')
    call getmem1d(bseas,1,nlandp,'bats_internal:bseas')
    call getmem1d(bb,1,nlandp,'bats_internal:bb')
    call getmem1d(bcoef,1,nlandp,'bats_internal:bcoef')
    call getmem1d(bfc,1,nlandp,'bats_internal:bfc')
    call getmem1d(bsw,1,nlandp,'bats_internal:bsw')
    call getmem1d(cc,1,nlandp,'bats_internal:cc')
    call getmem1d(cdr,1,nlandp,'bats_internal:cdr')
    call getmem1d(cdrd,1,nlandp,'bats_internal:cdrd')
    call getmem1d(cdrn,1,nlandp,'bats_internal:cdrn')
    call getmem1d(cdrx,1,nlandp,'bats_internal:cdrx')
    call getmem1d(cf,1,nlandp,'bats_internal:cf')
    call getmem1d(cgrnd,1,nlandp,'bats_internal:cgrnd')
    call getmem1d(cgrndl,1,nlandp,'bats_internal:cgrndl')
    call getmem1d(cgrnds,1,nlandp,'bats_internal:cgrnds')
    call getmem1d(cn1,1,nlandp,'bats_internal:cn1')
    call getmem1d(czenith,1,nlandp,'bats_internal:czenith')
    call getmem1d(dcd,1,nlandp,'bats_internal:dcd')
    call getmem1d(delq,1,nlandp,'bats_internal:delq')
    call getmem1d(dels,1,nlandp,'bats_internal:dels')
    call getmem1d(delt,1,nlandp,'bats_internal:delt')
    call getmem1d(deprat,1,nlandp,'bats_internal:deprat')
    call getmem1d(df,1,nlandp,'bats_internal:df')
    call getmem1d(drag,1,nlandp,'bats_internal:drag')
    call getmem1d(dzh,1,nlandp,'bats_internal:dzh')
    call getmem1d(efe,1,nlandp,'bats_internal:efe')
    call getmem1d(efpr,1,nlandp,'bats_internal:efpr')
    call getmem1d(eg,1,nlandp,'bats_internal:eg')
    call getmem1d(emiss,1,nlandp,'bats_internal:emiss')
    call getmem1d(etr,1,nlandp,'bats_internal:etr')
    call getmem1d(etrc,1,nlandp,'bats_internal:etrc')
    call getmem1d(etrrun,1,nlandp,'bats_internal:etrrun')
    call getmem1d(evaps,1,nlandp,'bats_internal:evaps')
    call getmem1d(evapw,1,nlandp,'bats_internal:evapw')
    call getmem1d(evmx0,1,nlandp,'bats_internal:evmx0')
    call getmem1d(evpr,1,nlandp,'bats_internal:evpr')
    call getmem1d(fact10,1,nlandp,'bats_internal:fact10')
    call getmem1d(fact2,1,nlandp,'bats_internal:fact2')
    call getmem1d(fct2,1,nlandp,'bats_internal:fct2')
    call getmem1d(fdry,1,nlandp,'bats_internal:fdry')
    call getmem1d(fevpg,1,nlandp,'bats_internal:fevpg')
    call getmem1d(flnet,1,nlandp,'bats_internal:flnet')
    call getmem1d(flneto,1,nlandp,'bats_internal:flneto')
    call getmem1d(fracd,1,nlandp,'bats_internal:fracd')
    call getmem1d(fseng,1,nlandp,'bats_internal:fseng')
    call getmem1d(fwet,1,nlandp,'bats_internal:fwet')
    call getmem1d(gwet,1,nlandp,'bats_internal:gwet')
    call getmem1d(gwmx0,1,nlandp,'bats_internal:gwmx0')
    call getmem1d(gwmx1,1,nlandp,'bats_internal:gwmx1')
    call getmem1d(gwmx2,1,nlandp,'bats_internal:gwmx2')
    call getmem1d(ht,1,nlandp,'bats_internal:ht')
    call getmem1d(hts,1,nlandp,'bats_internal:hts')
    call getmem1d(htvp,1,nlandp,'bats_internal:htvp')
    call getmem1d(lat,1,nlandp,'bats_internal:lat')
    call getmem1d(ldew,1,nlandp,'bats_internal:ldew')
    call getmem1d(lftra,1,nlandp,'bats_internal:lftra')
    call getmem1d(lftrs,1,nlandp,'bats_internal:lftrs')
    call getmem1d(lncl,1,nlandp,'bats_internal:lncl')
    call getmem1d(lwal,1,nlandp,'bats_internal:lwalb')
    call getmem1d(lwdifal,1,nlandp,'bats_internal:lwdifal')
    call getmem1d(lwdiral,1,nlandp,'bats_internal:lwdiral')
    call getmem1d(lwflx,1,nlandp,'bats_internal:lwflx')
    call getmem1d(porsl,1,nlandp,'bats_internal:porsl')
    call getmem1d(prcp,1,nlandp,'bats_internal:prcp')
    call getmem1d(ps,1,nlandp,'bats_internal:ps')
    call getmem1d(pw,1,nlandp,'bats_internal:pw')
    call getmem1d(qgrd,1,nlandp,'bats_internal:qgrd')
    call getmem1d(qs,1,nlandp,'bats_internal:qs')
    call getmem1d(qsatl,1,nlandp,'bats_internal:qsatl')
    ! call getmem1d(relaw,1,nlandp,'bats_internal:relaw')
    call getmem1d(relfc,1,nlandp,'bats_internal:relfc')
    call getmem1d(resp,1,nlandp,'bats_internal:resp')
    call getmem1d(rgr,1,nlandp,'bats_internal:rgr')
    call getmem1d(rhosw,1,nlandp,'bats_internal:rhosw')
    call getmem1d(rhs,1,nlandp,'bats_internal:rhs')
    call getmem1d(rib,1,nlandp,'bats_internal:rib')
    call getmem1d(ribd,1,nlandp,'bats_internal:ribd')
    call getmem1d(rlai,1,nlandp,'bats_internal:rlai')
    call getmem1d(rnet,1,nlandp,'bats_internal:rnet')
    call getmem1d(rnof,1,nlandp,'bats_internal:rnof')
    call getmem1d(rpp,1,nlandp,'bats_internal:rpp')
    call getmem1d(rppq,1,nlandp,'bats_internal:rppq')
    call getmem1d(rsubst,1,nlandp,'bats_internal:rsubst')
    call getmem1d(rsur,1,nlandp,'bats_internal:rsur')
    call getmem1d(rsw,1,nlandp,'bats_internal:rsw')
    call getmem1d(scrat,1,nlandp,'bats_internal:scrat')
    call getmem1d(scvk,1,nlandp,'bats_internal:scvk')
    call getmem1d(sdrop,1,nlandp,'bats_internal:sdrop')
    call getmem1d(sent,1,nlandp,'bats_internal:sent')
    call getmem1d(sfcp,1,nlandp,'bats_internal:sfcp')
    call getmem1d(sigf,1,nlandp,'bats_internal:sigf')
    call getmem1d(sm,1,nlandp,'bats_internal:sm')
    call getmem1d(snag,1,nlandp,'bats_internal:snag')
    call getmem1d(sncv,1,nlandp,'bats_internal:sncv')
    call getmem1d(srnof,1,nlandp,'bats_internal:srnof')
    call getmem1d(ssw,1,nlandp,'bats_internal:ssw')
    call getmem1d(sts,1,nlandp,'bats_internal:sts')
    call getmem1d(swal,1,nlandp,'bats_internal:swal')
    call getmem1d(swdifal,1,nlandp,'bats_internal:swdifal')
    call getmem1d(swdiral,1,nlandp,'bats_internal:swdiral')
    call getmem1d(swflx,1,nlandp,'bats_internal:swflx')
    call getmem1d(swsi,1,nlandp,'bats_internal:swsi')
    call getmem1d(taf,1,nlandp,'bats_internal:taf')
    call getmem1d(texrat,1,nlandp,'bats_internal:texrat')
    call getmem1d(tgbb,1,nlandp,'bats_internal:tgbb')
    call getmem1d(tgbrd,1,nlandp,'bats_internal:tgbrd')
    call getmem1d(tgrd,1,nlandp,'bats_internal:tgrd')
    call getmem1d(tlef,1,nlandp,'bats_internal:tlef')
    call getmem1d(tm,1,nlandp,'bats_internal:tm')
    call getmem1d(trnof,1,nlandp,'bats_internal:trnof')
    call getmem1d(tsw,1,nlandp,'bats_internal:tsw')
    call getmem1d(uaf,1,nlandp,'bats_internal:uaf')
    call getmem1d(usw,1,nlandp,'bats_internal:usw')
    call getmem1d(vegt,1,nlandp,'bats_internal:vegt')
    call getmem1d(vpdc,1,nlandp,'bats_internal:vpdc')
    call getmem1d(vspda,1,nlandp,'bats_internal:vspda')
    call getmem1d(vsw,1,nlandp,'bats_internal:vsw')
    call getmem1d(wata,1,nlandp,'bats_internal:wata')
    call getmem1d(watr,1,nlandp,'bats_internal:watr')
    call getmem1d(watt,1,nlandp,'bats_internal:watt')
    call getmem1d(watu,1,nlandp,'bats_internal:watu')
    call getmem1d(wflux1,1,nlandp,'bats_internal:wflux1')
    call getmem1d(wflux2,1,nlandp,'bats_internal:wflux2')
    call getmem1d(wfluxc,1,nlandp,'bats_internal:wfluxc')
    call getmem1d(wiltr,1,nlandp,'bats_internal:wiltr')
    call getmem1d(wt,1,nlandp,'bats_internal:wt')
    call getmem1d(wta0,1,nlandp,'bats_internal:wta0')
    call getmem1d(wta,1,nlandp,'bats_internal:wta')
    call getmem1d(wtaq0,1,nlandp,'bats_internal:wtaq0')
    call getmem1d(wtg0,1,nlandp,'bats_internal:wtg0')
    call getmem1d(wtg,1,nlandp,'bats_internal:wtg')
    call getmem1d(wtg2,1,nlandp,'bats_internal:wtg2')
    call getmem1d(wtga,1,nlandp,'bats_internal:wtga')
    call getmem1d(wtgaq,1,nlandp,'bats_internal:wtgaq')
    call getmem1d(wtgl,1,nlandp,'bats_internal:wtgl')
    call getmem1d(wtglq,1,nlandp,'bats_internal:wtglq')
    call getmem1d(wtgq0,1,nlandp,'bats_internal:wtgq0')
    call getmem1d(wtgq,1,nlandp,'bats_internal:wtgq')
    call getmem1d(wtl0,1,nlandp,'bats_internal:wtl0')
    call getmem1d(wtlh,1,nlandp,'bats_internal:wtlh')
    call getmem1d(wtlq0,1,nlandp,'bats_internal:wtlq0')
    call getmem1d(wtlq,1,nlandp,'bats_internal:wtlq')
    call getmem1d(wtshi,1,nlandp,'bats_internal:wtshi')
    call getmem1d(wtsqi,1,nlandp,'bats_internal:wtsqi')
    call getmem1d(xkmx,1,nlandp,'bats_internal:xkmx')
    call getmem1d(xlai,1,nlandp,'bats_internal:xlai')
    call getmem1d(xlsai,1,nlandp,'bats_internal:xlsai')
    call getmem1d(xrun,1,nlandp,'bats_internal:xrun')
    call getmem1d(z10fra,1,nlandp,'bats_internal:z10fra')
    call getmem1d(z1log,1,nlandp,'bats_internal:z1log')
    call getmem1d(z2fra,1,nlandp,'bats_internal:z2fra')
    call getmem1d(zh,1,nlandp,'bats_internal:zh')
    call getmem1d(zlgdis,1,nlandp,'bats_internal:zlgdis')
    call getmem1d(zlglnd,1,nlandp,'bats_internal:zlglnd')
    call getmem1d(zlgsno,1,nlandp,'bats_internal:zlgsno')
    call getmem1d(zlgveg,1,nlandp,'bats_internal:zlgveg')
    call getmem1d(lveg,1,nlandp,'bats_internal:lveg')
    call getmem1d(p0,1,nlandp,'bats_internal:p0')
    call getmem1d(qs0,1,nlandp,'bats_internal:qs0')
    call getmem1d(ts0,1,nlandp,'bats_internal:ts0')
    call getmem1d(swd0,1,nlandp,'bats_internal:swd0')
    call getmem1d(swf0,1,nlandp,'bats_internal:swf0')
    call getmem1d(ncpr0,1,nlandp,'bats_internal:ncpr0')
    call getmem1d(cpr0,1,nlandp,'bats_internal:cpr0')
  end subroutine allocate_mod_bats_internal
!
end module mod_bats_internal
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
