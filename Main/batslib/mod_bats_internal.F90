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

  real(rkx) :: dtbat

  integer(ik4) :: nlandp
  integer(ik4) :: ilndbeg , ilndend

  real(rkx) , pointer , dimension(:) :: abswveg
  real(rkx) , pointer , dimension(:) :: aseas
  real(rkx) , pointer , dimension(:) :: bseas
  real(rkx) , pointer , dimension(:) :: bb
  real(rkx) , pointer , dimension(:) :: bcoef
  real(rkx) , pointer , dimension(:) :: bfc
  real(rkx) , pointer , dimension(:) :: bsw
  real(rkx) , pointer , dimension(:) :: cc
  real(rkx) , pointer , dimension(:) :: cdr
  real(rkx) , pointer , dimension(:) :: cdrd
  real(rkx) , pointer , dimension(:) :: cdrn
  real(rkx) , pointer , dimension(:) :: cdrx
  real(rkx) , pointer , dimension(:) :: cf
  real(rkx) , pointer , dimension(:) :: cgrnd
  real(rkx) , pointer , dimension(:) :: cgrndl
  real(rkx) , pointer , dimension(:) :: cgrnds
  real(rkx) , pointer , dimension(:) :: cn1
  real(rkx) , pointer , dimension(:) :: czenith
  real(rkx) , pointer , dimension(:) :: dcd
  real(rkx) , pointer , dimension(:) :: delq
  real(rkx) , pointer , dimension(:) :: dels
  real(rkx) , pointer , dimension(:) :: delt
  real(rkx) , pointer , dimension(:) :: deprat
  real(rkx) , pointer , dimension(:) :: df
  real(rkx) , pointer , dimension(:) :: drag
  real(rkx) , pointer , dimension(:) :: dzh
  real(rkx) , pointer , dimension(:) :: efe
  real(rkx) , pointer , dimension(:) :: efpr
  real(rkx) , pointer , dimension(:) :: eg
  real(rkx) , pointer , dimension(:) :: emiss
  real(rkx) , pointer , dimension(:) :: etr
  real(rkx) , pointer , dimension(:) :: etrc
  real(rkx) , pointer , dimension(:) :: etrrun
  real(rkx) , pointer , dimension(:) :: evaps
  real(rkx) , pointer , dimension(:) :: evapw
  real(rkx) , pointer , dimension(:) :: evmx0
  real(rkx) , pointer , dimension(:) :: evpr
  real(rkx) , pointer , dimension(:) :: fact10
  real(rkx) , pointer , dimension(:) :: fact2
  real(rkx) , pointer , dimension(:) :: fct2
  real(rkx) , pointer , dimension(:) :: fdry
  real(rkx) , pointer , dimension(:) :: fevpg
  real(rkx) , pointer , dimension(:) :: flnet
  real(rkx) , pointer , dimension(:) :: flneto
  real(rkx) , pointer , dimension(:) :: fracd
  real(rkx) , pointer , dimension(:) :: fseng
  real(rkx) , pointer , dimension(:) :: fwet
  real(rkx) , pointer , dimension(:) :: gwet
  real(rkx) , pointer , dimension(:) :: gwmx0
  real(rkx) , pointer , dimension(:) :: gwmx1
  real(rkx) , pointer , dimension(:) :: gwmx2
  real(rkx) , pointer , dimension(:) :: ht
  real(rkx) , pointer , dimension(:) :: hts
  real(rkx) , pointer , dimension(:) :: htvp
  real(rkx) , pointer , dimension(:) :: lat
  real(rkx) , pointer , dimension(:) :: ldew
  real(rkx) , pointer , dimension(:) :: lftra
  real(rkx) , pointer , dimension(:) :: lftrs
  real(rkx) , pointer , dimension(:) :: lncl
  real(rkx) , pointer , dimension(:) :: lwal
  real(rkx) , pointer , dimension(:) :: lwdifal
  real(rkx) , pointer , dimension(:) :: lwdiral
  real(rkx) , pointer , dimension(:) :: lwflx
  real(rkx) , pointer , dimension(:) :: porsl
  real(rkx) , pointer , dimension(:) :: prcp
  real(rkx) , pointer , dimension(:) :: ps
  real(rkx) , pointer , dimension(:) :: pw
  real(rkx) , pointer , dimension(:) :: qgrd
  real(rkx) , pointer , dimension(:) :: qs
  real(rkx) , pointer , dimension(:) :: qsatl
  ! real(rkx) , pointer , dimension(:) :: relaw
  real(rkx) , pointer , dimension(:) :: relfc
  real(rkx) , pointer , dimension(:) :: resp
  real(rkx) , pointer , dimension(:) :: rgr
  real(rkx) , pointer , dimension(:) :: rhosw
  real(rkx) , pointer , dimension(:) :: rhs
  real(rkx) , pointer , dimension(:) :: rib
  real(rkx) , pointer , dimension(:) :: ribd
  real(rkx) , pointer , dimension(:) :: rlai
  real(rkx) , pointer , dimension(:) :: rnet
  real(rkx) , pointer , dimension(:) :: rnof
  real(rkx) , pointer , dimension(:) :: rpp
  real(rkx) , pointer , dimension(:) :: rppq
  real(rkx) , pointer , dimension(:) :: rsubst
  real(rkx) , pointer , dimension(:) :: rsur
  real(rkx) , pointer , dimension(:) :: rsw
  real(rkx) , pointer , dimension(:) :: scrat
  real(rkx) , pointer , dimension(:) :: scvk
  real(rkx) , pointer , dimension(:) :: sdrop
  real(rkx) , pointer , dimension(:) :: sent
  real(rkx) , pointer , dimension(:) :: sfcp
  real(rkx) , pointer , dimension(:) :: sigf
  real(rkx) , pointer , dimension(:) :: sm
  real(rkx) , pointer , dimension(:) :: snag
  real(rkx) , pointer , dimension(:) :: sncv
  real(rkx) , pointer , dimension(:) :: srnof
  real(rkx) , pointer , dimension(:) :: ssw
  real(rkx) , pointer , dimension(:) :: sts
  real(rkx) , pointer , dimension(:) :: swal
  real(rkx) , pointer , dimension(:) :: swdifal
  real(rkx) , pointer , dimension(:) :: swdiral
  real(rkx) , pointer , dimension(:) :: swflx
  real(rkx) , pointer , dimension(:) :: swsi
  real(rkx) , pointer , dimension(:) :: taf
  real(rkx) , pointer , dimension(:) :: texrat
  real(rkx) , pointer , dimension(:) :: tgbb
  real(rkx) , pointer , dimension(:) :: tgbrd
  real(rkx) , pointer , dimension(:) :: tgrd
  real(rkx) , pointer , dimension(:) :: tlef
  real(rkx) , pointer , dimension(:) :: tm
  real(rkx) , pointer , dimension(:) :: trnof
  real(rkx) , pointer , dimension(:) :: tsw
  real(rkx) , pointer , dimension(:) :: uaf
  real(rkx) , pointer , dimension(:) :: usw
  real(rkx) , pointer , dimension(:) :: vegt
  real(rkx) , pointer , dimension(:) :: vpdc
  real(rkx) , pointer , dimension(:) :: vspda
  real(rkx) , pointer , dimension(:) :: vsw
  real(rkx) , pointer , dimension(:) :: watr
  real(rkx) , pointer , dimension(:) :: watt
  real(rkx) , pointer , dimension(:) :: watu
  real(rkx) , pointer , dimension(:) :: wflux1
  real(rkx) , pointer , dimension(:) :: wflux2
  real(rkx) , pointer , dimension(:) :: wfluxc
  real(rkx) , pointer , dimension(:) :: wiltr
  real(rkx) , pointer , dimension(:) :: wt
  real(rkx) , pointer , dimension(:) :: wta
  real(rkx) , pointer , dimension(:) :: wta0
  real(rkx) , pointer , dimension(:) :: wtaq0
  real(rkx) , pointer , dimension(:) :: wtg
  real(rkx) , pointer , dimension(:) :: wtg0
  real(rkx) , pointer , dimension(:) :: wtg2
  real(rkx) , pointer , dimension(:) :: wtga
  real(rkx) , pointer , dimension(:) :: wtgaq
  real(rkx) , pointer , dimension(:) :: wtgl
  real(rkx) , pointer , dimension(:) :: wtglq
  real(rkx) , pointer , dimension(:) :: wtgq
  real(rkx) , pointer , dimension(:) :: wtgq0
  real(rkx) , pointer , dimension(:) :: wtl0
  real(rkx) , pointer , dimension(:) :: wtlh
  real(rkx) , pointer , dimension(:) :: wtlq
  real(rkx) , pointer , dimension(:) :: wtlq0
  real(rkx) , pointer , dimension(:) :: wtshi
  real(rkx) , pointer , dimension(:) :: wtsqi
  real(rkx) , pointer , dimension(:) :: xkmx
  real(rkx) , pointer , dimension(:) :: xlai
  real(rkx) , pointer , dimension(:) :: xlsai
  real(rkx) , pointer , dimension(:) :: z10fra
  real(rkx) , pointer , dimension(:) :: z1log
  real(rkx) , pointer , dimension(:) :: z2fra
  real(rkx) , pointer , dimension(:) :: zh
  real(rkx) , pointer , dimension(:) :: zo
  real(rkx) , pointer , dimension(:) :: zlgdis
  real(rkx) , pointer , dimension(:) :: zlglnd
  real(rkx) , pointer , dimension(:) :: zlgsno
  real(rkx) , pointer , dimension(:) :: zlgveg

  integer(ik4) , pointer , dimension(:) :: lveg

  real(rkx) , pointer , dimension(:) :: p0
  real(rkx) , pointer , dimension(:) :: qs0
  real(rkx) , pointer , dimension(:) :: ts0
  real(rkx) , pointer , dimension(:) :: swd0
  real(rkx) , pointer , dimension(:) :: swf0
  real(rkx) , pointer , dimension(:) :: ncpr0
  real(rkx) , pointer , dimension(:) :: cpr0

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
    call getmem1d(z10fra,1,nlandp,'bats_internal:z10fra')
    call getmem1d(z1log,1,nlandp,'bats_internal:z1log')
    call getmem1d(z2fra,1,nlandp,'bats_internal:z2fra')
    call getmem1d(zh,1,nlandp,'bats_internal:zh')
    call getmem1d(zo,1,nlandp,'bats_internal:zo')
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

  subroutine bats_psat(t,e)
    implicit none
    real(rkx) , intent(in) :: t
    real(rkx) , intent(out) :: e
    if ( t < tzero ) then
      e = c1es * exp(c3ies*(t-tzero)/(t-c4ies))
    else
      e = c1es * exp(c3les*(t-tzero)/(t-c4les))
    end if
  end subroutine bats_psat

  subroutine bats_satur(t,p,e,qs)
    implicit none
    real(rkx) , intent(in) :: t , p
    real(rkx) , intent(out) :: e , qs
    call bats_psat(t,e)
    qs = ep2 * e/(p-0.378_rkx*e)
  end subroutine bats_satur

  subroutine bats_qsdt(t,qs,qsdt)
    implicit none
    real(rkx) , intent(in) :: t , qs
    real(rkx) , intent(out) :: qsdt
    ! Eqn. 83 from BATS manual
    if ( t < tzero ) then
      qsdt = qs * (c3ies * (tzero-c4ies) * (d_one/(t-c4ies))**d_two)
    else
      qsdt = qs * (c3les * (tzero-c4les) * (d_one/(t-c4les))**d_two)
    end if
  end subroutine bats_qsdt

end module mod_bats_internal

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
