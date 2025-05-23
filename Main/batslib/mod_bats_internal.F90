!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    Use of this source code is governed by an MIT-style license that can
!    be found in the LICENSE file or at
!
!         https://opensource.org/licenses/MIT.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

module mod_bats_internal
!
  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_memutil
  use mod_dynparam
  use mod_regcm_types

  implicit none

  public

  real(rkx) :: dtbat

  integer(ik4) :: nlandp
  integer(ik4) :: ilndbeg, ilndend

  real(rkx), pointer, contiguous, dimension(:) :: abswveg
  real(rkx), pointer, contiguous, dimension(:) :: aseas
  real(rkx), pointer, contiguous, dimension(:) :: bseas
  real(rkx), pointer, contiguous, dimension(:) :: bb
  real(rkx), pointer, contiguous, dimension(:) :: bcoef
  real(rkx), pointer, contiguous, dimension(:) :: bfc
  real(rkx), pointer, contiguous, dimension(:) :: bsw
  real(rkx), pointer, contiguous, dimension(:) :: cc
  real(rkx), pointer, contiguous, dimension(:) :: cdr
  real(rkx), pointer, contiguous, dimension(:) :: cdrd
  real(rkx), pointer, contiguous, dimension(:) :: cdrn
  real(rkx), pointer, contiguous, dimension(:) :: cdrx
  real(rkx), pointer, contiguous, dimension(:) :: cf
  real(rkx), pointer, contiguous, dimension(:) :: cgrnd
  real(rkx), pointer, contiguous, dimension(:) :: cgrndl
  real(rkx), pointer, contiguous, dimension(:) :: cgrnds
  real(rkx), pointer, contiguous, dimension(:) :: cn1
  real(rkx), pointer, contiguous, dimension(:) :: czenith
  real(rkx), pointer, contiguous, dimension(:) :: dcd
  real(rkx), pointer, contiguous, dimension(:) :: delq
  real(rkx), pointer, contiguous, dimension(:) :: delt
  real(rkx), pointer, contiguous, dimension(:) :: deprat
  real(rkx), pointer, contiguous, dimension(:) :: df
  real(rkx), pointer, contiguous, dimension(:) :: drag
  real(rkx), pointer, contiguous, dimension(:) :: dzh
  real(rkx), pointer, contiguous, dimension(:) :: efe
  real(rkx), pointer, contiguous, dimension(:) :: efpr
  real(rkx), pointer, contiguous, dimension(:) :: eg
  real(rkx), pointer, contiguous, dimension(:) :: emiss
  real(rkx), pointer, contiguous, dimension(:) :: etr
  real(rkx), pointer, contiguous, dimension(:) :: etrc
  real(rkx), pointer, contiguous, dimension(:) :: etrrun
  real(rkx), pointer, contiguous, dimension(:) :: evaps
  real(rkx), pointer, contiguous, dimension(:) :: evapw
  real(rkx), pointer, contiguous, dimension(:) :: evmx0
  real(rkx), pointer, contiguous, dimension(:) :: evpr
  real(rkx), pointer, contiguous, dimension(:) :: fact10
  real(rkx), pointer, contiguous, dimension(:) :: fact2
  real(rkx), pointer, contiguous, dimension(:) :: fct2
  real(rkx), pointer, contiguous, dimension(:) :: fdry
  real(rkx), pointer, contiguous, dimension(:) :: fevpg
  real(rkx), pointer, contiguous, dimension(:) :: flnet
  real(rkx), pointer, contiguous, dimension(:) :: flneto
  real(rkx), pointer, contiguous, dimension(:) :: fracd
  real(rkx), pointer, contiguous, dimension(:) :: fseng
  real(rkx), pointer, contiguous, dimension(:) :: fwet
  real(rkx), pointer, contiguous, dimension(:) :: gwet
  real(rkx), pointer, contiguous, dimension(:) :: gwmx0
  real(rkx), pointer, contiguous, dimension(:) :: gwmx1
  real(rkx), pointer, contiguous, dimension(:) :: gwmx2
  real(rkx), pointer, contiguous, dimension(:) :: ht
  real(rkx), pointer, contiguous, dimension(:) :: hts
  real(rkx), pointer, contiguous, dimension(:) :: htvp
  real(rkx), pointer, contiguous, dimension(:) :: lat
  real(rkx), pointer, contiguous, dimension(:) :: ldew
  real(rkx), pointer, contiguous, dimension(:) :: lftra
  real(rkx), pointer, contiguous, dimension(:) :: lftrs
  real(rkx), pointer, contiguous, dimension(:) :: lncl
  real(rkx), pointer, contiguous, dimension(:) :: lwal
  real(rkx), pointer, contiguous, dimension(:) :: lwdifal
  real(rkx), pointer, contiguous, dimension(:) :: lwdiral
  real(rkx), pointer, contiguous, dimension(:) :: lwflx
  real(rkx), pointer, contiguous, dimension(:) :: porsl
  real(rkx), pointer, contiguous, dimension(:) :: prcp
  real(rkx), pointer, contiguous, dimension(:) :: ps
  real(rkx), pointer, contiguous, dimension(:) :: pw
  real(rkx), pointer, contiguous, dimension(:) :: qgrd
  real(rkx), pointer, contiguous, dimension(:) :: qs
  real(rkx), pointer, contiguous, dimension(:) :: qsatl
  ! real(rkx), pointer, contiguous, dimension(:) :: relaw
  real(rkx), pointer, contiguous, dimension(:) :: relfc
  real(rkx), pointer, contiguous, dimension(:) :: resp
  real(rkx), pointer, contiguous, dimension(:) :: rgr
  real(rkx), pointer, contiguous, dimension(:) :: rhosw
  real(rkx), pointer, contiguous, dimension(:) :: rhs
  real(rkx), pointer, contiguous, dimension(:) :: rib
  real(rkx), pointer, contiguous, dimension(:) :: ribd
  real(rkx), pointer, contiguous, dimension(:) :: rlai
  real(rkx), pointer, contiguous, dimension(:) :: rnet
  real(rkx), pointer, contiguous, dimension(:) :: rpp
  real(rkx), pointer, contiguous, dimension(:) :: rppq
  real(rkx), pointer, contiguous, dimension(:) :: rsubst
  real(rkx), pointer, contiguous, dimension(:) :: rsw
  real(rkx), pointer, contiguous, dimension(:) :: scrat
  real(rkx), pointer, contiguous, dimension(:) :: scvk
  real(rkx), pointer, contiguous, dimension(:) :: sdrop
  real(rkx), pointer, contiguous, dimension(:) :: sent
  real(rkx), pointer, contiguous, dimension(:) :: sfcp
  real(rkx), pointer, contiguous, dimension(:) :: sigf
  real(rkx), pointer, contiguous, dimension(:) :: sm
  real(rkx), pointer, contiguous, dimension(:) :: snag
  real(rkx), pointer, contiguous, dimension(:) :: sncv
  real(rkx), pointer, contiguous, dimension(:) :: srnof
  real(rkx), pointer, contiguous, dimension(:) :: ssw
  real(rkx), pointer, contiguous, dimension(:) :: sts
  real(rkx), pointer, contiguous, dimension(:) :: swal
  real(rkx), pointer, contiguous, dimension(:) :: swdifal
  real(rkx), pointer, contiguous, dimension(:) :: swdiral
  real(rkx), pointer, contiguous, dimension(:) :: swflx
  real(rkx), pointer, contiguous, dimension(:) :: swsi
  real(rkx), pointer, contiguous, dimension(:) :: taf
  real(rkx), pointer, contiguous, dimension(:) :: texrat
  real(rkx), pointer, contiguous, dimension(:) :: tgbb
  real(rkx), pointer, contiguous, dimension(:) :: tgbrd
  real(rkx), pointer, contiguous, dimension(:) :: tgrd
  real(rkx), pointer, contiguous, dimension(:) :: tlef
  real(rkx), pointer, contiguous, dimension(:) :: tm
  real(rkx), pointer, contiguous, dimension(:) :: trnof
  real(rkx), pointer, contiguous, dimension(:) :: tsw
  real(rkx), pointer, contiguous, dimension(:) :: uaf
  real(rkx), pointer, contiguous, dimension(:) :: usw
  real(rkx), pointer, contiguous, dimension(:) :: vegt
  real(rkx), pointer, contiguous, dimension(:) :: vpdc
  real(rkx), pointer, contiguous, dimension(:) :: vspda
  real(rkx), pointer, contiguous, dimension(:) :: vsw
  real(rkx), pointer, contiguous, dimension(:) :: watr
  real(rkx), pointer, contiguous, dimension(:) :: watt
  real(rkx), pointer, contiguous, dimension(:) :: watu
  real(rkx), pointer, contiguous, dimension(:) :: wflux1
  real(rkx), pointer, contiguous, dimension(:) :: wflux2
  real(rkx), pointer, contiguous, dimension(:) :: wfluxc
  real(rkx), pointer, contiguous, dimension(:) :: wiltr
  real(rkx), pointer, contiguous, dimension(:) :: wt
  real(rkx), pointer, contiguous, dimension(:) :: wta
  real(rkx), pointer, contiguous, dimension(:) :: wta0
  real(rkx), pointer, contiguous, dimension(:) :: wtaq0
  real(rkx), pointer, contiguous, dimension(:) :: wtg
  real(rkx), pointer, contiguous, dimension(:) :: wtg0
  real(rkx), pointer, contiguous, dimension(:) :: wtg2
  real(rkx), pointer, contiguous, dimension(:) :: wtga
  real(rkx), pointer, contiguous, dimension(:) :: wtgaq
  real(rkx), pointer, contiguous, dimension(:) :: wtgl
  real(rkx), pointer, contiguous, dimension(:) :: wtglq
  real(rkx), pointer, contiguous, dimension(:) :: wtgq
  real(rkx), pointer, contiguous, dimension(:) :: wtgq0
  real(rkx), pointer, contiguous, dimension(:) :: wtl0
  real(rkx), pointer, contiguous, dimension(:) :: wtlh
  real(rkx), pointer, contiguous, dimension(:) :: wtlq
  real(rkx), pointer, contiguous, dimension(:) :: wtlq0
  real(rkx), pointer, contiguous, dimension(:) :: wtshi
  real(rkx), pointer, contiguous, dimension(:) :: wtsqi
  real(rkx), pointer, contiguous, dimension(:) :: xkmx
  real(rkx), pointer, contiguous, dimension(:) :: xlai
  real(rkx), pointer, contiguous, dimension(:) :: xlsai
  real(rkx), pointer, contiguous, dimension(:) :: z10fra
  real(rkx), pointer, contiguous, dimension(:) :: z1log
  real(rkx), pointer, contiguous, dimension(:) :: z2fra
  real(rkx), pointer, contiguous, dimension(:) :: zh
  real(rkx), pointer, contiguous, dimension(:) :: zo
  real(rkx), pointer, contiguous, dimension(:) :: zlgdis
  real(rkx), pointer, contiguous, dimension(:) :: zlglnd
  real(rkx), pointer, contiguous, dimension(:) :: zlgsno
  real(rkx), pointer, contiguous, dimension(:) :: zlgveg

  integer(ik4), pointer, contiguous, dimension(:) :: lveg
  integer(ik4), pointer, contiguous, dimension(:) :: ltex

  real(rkx), pointer, contiguous, dimension(:) :: p0
  real(rkx), pointer, contiguous, dimension(:) :: qs0
  real(rkx), pointer, contiguous, dimension(:) :: ts0
  real(rkx), pointer, contiguous, dimension(:) :: swd0
  real(rkx), pointer, contiguous, dimension(:) :: swf0
  real(rkx), pointer, contiguous, dimension(:) :: ncpr0
  real(rkx), pointer, contiguous, dimension(:) :: cpr0

  contains

  subroutine allocate_mod_bats_internal(cl)
    implicit none
    type (masked_comm), intent(in) :: cl
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
    call getmem1d(rpp,1,nlandp,'bats_internal:rpp')
    call getmem1d(rppq,1,nlandp,'bats_internal:rppq')
    call getmem1d(rsubst,1,nlandp,'bats_internal:rsubst')
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
    call getmem1d(ltex,1,nlandp,'bats_internal:ltex')
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
    real(rkx), intent(in) :: t
    real(rkx), intent(out) :: e
    if ( t < tzero ) then
      e = c1es * exp(c3ies*(t-tzero)/(t-c4ies))
    else
      e = c1es * exp(c3les*(t-tzero)/(t-c4les))
    end if
  end subroutine bats_psat

  subroutine bats_satur(t,p,e,qs)
    implicit none
    real(rkx), intent(in) :: t, p
    real(rkx), intent(out) :: e, qs
    call bats_psat(t,e)
    qs = ep2 * e/(p-0.378_rkx*e)
  end subroutine bats_satur

  subroutine bats_qsdt(t,qs,qsdt)
    implicit none
    real(rkx), intent(in) :: t, qs
    real(rkx), intent(out) :: qsdt
    ! Eqn. 83 from BATS manual
    if ( t < tzero ) then
      qsdt = qs * (c3ies * (tzero-c4ies) * (d_one/(t-c4ies))**d_two)
    else
      qsdt = qs * (c3les * (tzero-c4les) * (d_one/(t-c4les))**d_two)
    end if
  end subroutine bats_qsdt

end module mod_bats_internal

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
