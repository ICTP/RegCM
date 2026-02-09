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
    call getmem(abswveg,1,nlandp,'bats_internal:abswveg')
    call getmem(aseas,1,nlandp,'bats_internal:aseas')
    call getmem(bseas,1,nlandp,'bats_internal:bseas')
    call getmem(bb,1,nlandp,'bats_internal:bb')
    call getmem(bcoef,1,nlandp,'bats_internal:bcoef')
    call getmem(bfc,1,nlandp,'bats_internal:bfc')
    call getmem(bsw,1,nlandp,'bats_internal:bsw')
    call getmem(cc,1,nlandp,'bats_internal:cc')
    call getmem(cdr,1,nlandp,'bats_internal:cdr')
    call getmem(cdrd,1,nlandp,'bats_internal:cdrd')
    call getmem(cdrn,1,nlandp,'bats_internal:cdrn')
    call getmem(cdrx,1,nlandp,'bats_internal:cdrx')
    call getmem(cf,1,nlandp,'bats_internal:cf')
    call getmem(cgrnd,1,nlandp,'bats_internal:cgrnd')
    call getmem(cgrndl,1,nlandp,'bats_internal:cgrndl')
    call getmem(cgrnds,1,nlandp,'bats_internal:cgrnds')
    call getmem(cn1,1,nlandp,'bats_internal:cn1')
    call getmem(czenith,1,nlandp,'bats_internal:czenith')
    call getmem(dcd,1,nlandp,'bats_internal:dcd')
    call getmem(delq,1,nlandp,'bats_internal:delq')
    call getmem(delt,1,nlandp,'bats_internal:delt')
    call getmem(deprat,1,nlandp,'bats_internal:deprat')
    call getmem(df,1,nlandp,'bats_internal:df')
    call getmem(drag,1,nlandp,'bats_internal:drag')
    call getmem(dzh,1,nlandp,'bats_internal:dzh')
    call getmem(efe,1,nlandp,'bats_internal:efe')
    call getmem(efpr,1,nlandp,'bats_internal:efpr')
    call getmem(eg,1,nlandp,'bats_internal:eg')
    call getmem(emiss,1,nlandp,'bats_internal:emiss')
    call getmem(etr,1,nlandp,'bats_internal:etr')
    call getmem(etrc,1,nlandp,'bats_internal:etrc')
    call getmem(etrrun,1,nlandp,'bats_internal:etrrun')
    call getmem(evaps,1,nlandp,'bats_internal:evaps')
    call getmem(evapw,1,nlandp,'bats_internal:evapw')
    call getmem(evmx0,1,nlandp,'bats_internal:evmx0')
    call getmem(evpr,1,nlandp,'bats_internal:evpr')
    call getmem(fact10,1,nlandp,'bats_internal:fact10')
    call getmem(fact2,1,nlandp,'bats_internal:fact2')
    call getmem(fct2,1,nlandp,'bats_internal:fct2')
    call getmem(fdry,1,nlandp,'bats_internal:fdry')
    call getmem(fevpg,1,nlandp,'bats_internal:fevpg')
    call getmem(flnet,1,nlandp,'bats_internal:flnet')
    call getmem(flneto,1,nlandp,'bats_internal:flneto')
    call getmem(fracd,1,nlandp,'bats_internal:fracd')
    call getmem(fseng,1,nlandp,'bats_internal:fseng')
    call getmem(fwet,1,nlandp,'bats_internal:fwet')
    call getmem(gwet,1,nlandp,'bats_internal:gwet')
    call getmem(gwmx0,1,nlandp,'bats_internal:gwmx0')
    call getmem(gwmx1,1,nlandp,'bats_internal:gwmx1')
    call getmem(gwmx2,1,nlandp,'bats_internal:gwmx2')
    call getmem(ht,1,nlandp,'bats_internal:ht')
    call getmem(hts,1,nlandp,'bats_internal:hts')
    call getmem(htvp,1,nlandp,'bats_internal:htvp')
    call getmem(lat,1,nlandp,'bats_internal:lat')
    call getmem(ldew,1,nlandp,'bats_internal:ldew')
    call getmem(lftra,1,nlandp,'bats_internal:lftra')
    call getmem(lftrs,1,nlandp,'bats_internal:lftrs')
    call getmem(lncl,1,nlandp,'bats_internal:lncl')
    call getmem(lwal,1,nlandp,'bats_internal:lwalb')
    call getmem(lwdifal,1,nlandp,'bats_internal:lwdifal')
    call getmem(lwdiral,1,nlandp,'bats_internal:lwdiral')
    call getmem(lwflx,1,nlandp,'bats_internal:lwflx')
    call getmem(porsl,1,nlandp,'bats_internal:porsl')
    call getmem(prcp,1,nlandp,'bats_internal:prcp')
    call getmem(ps,1,nlandp,'bats_internal:ps')
    call getmem(pw,1,nlandp,'bats_internal:pw')
    call getmem(qgrd,1,nlandp,'bats_internal:qgrd')
    call getmem(qs,1,nlandp,'bats_internal:qs')
    call getmem(qsatl,1,nlandp,'bats_internal:qsatl')
    ! call getmem(relaw,1,nlandp,'bats_internal:relaw')
    call getmem(relfc,1,nlandp,'bats_internal:relfc')
    call getmem(resp,1,nlandp,'bats_internal:resp')
    call getmem(rgr,1,nlandp,'bats_internal:rgr')
    call getmem(rhosw,1,nlandp,'bats_internal:rhosw')
    call getmem(rhs,1,nlandp,'bats_internal:rhs')
    call getmem(rib,1,nlandp,'bats_internal:rib')
    call getmem(ribd,1,nlandp,'bats_internal:ribd')
    call getmem(rlai,1,nlandp,'bats_internal:rlai')
    call getmem(rnet,1,nlandp,'bats_internal:rnet')
    call getmem(rpp,1,nlandp,'bats_internal:rpp')
    call getmem(rppq,1,nlandp,'bats_internal:rppq')
    call getmem(rsubst,1,nlandp,'bats_internal:rsubst')
    call getmem(rsw,1,nlandp,'bats_internal:rsw')
    call getmem(scrat,1,nlandp,'bats_internal:scrat')
    call getmem(scvk,1,nlandp,'bats_internal:scvk')
    call getmem(sdrop,1,nlandp,'bats_internal:sdrop')
    call getmem(sent,1,nlandp,'bats_internal:sent')
    call getmem(sfcp,1,nlandp,'bats_internal:sfcp')
    call getmem(sigf,1,nlandp,'bats_internal:sigf')
    call getmem(sm,1,nlandp,'bats_internal:sm')
    call getmem(snag,1,nlandp,'bats_internal:snag')
    call getmem(sncv,1,nlandp,'bats_internal:sncv')
    call getmem(srnof,1,nlandp,'bats_internal:srnof')
    call getmem(ssw,1,nlandp,'bats_internal:ssw')
    call getmem(sts,1,nlandp,'bats_internal:sts')
    call getmem(swal,1,nlandp,'bats_internal:swal')
    call getmem(swdifal,1,nlandp,'bats_internal:swdifal')
    call getmem(swdiral,1,nlandp,'bats_internal:swdiral')
    call getmem(swflx,1,nlandp,'bats_internal:swflx')
    call getmem(swsi,1,nlandp,'bats_internal:swsi')
    call getmem(taf,1,nlandp,'bats_internal:taf')
    call getmem(texrat,1,nlandp,'bats_internal:texrat')
    call getmem(tgbb,1,nlandp,'bats_internal:tgbb')
    call getmem(tgbrd,1,nlandp,'bats_internal:tgbrd')
    call getmem(tgrd,1,nlandp,'bats_internal:tgrd')
    call getmem(tlef,1,nlandp,'bats_internal:tlef')
    call getmem(tm,1,nlandp,'bats_internal:tm')
    call getmem(trnof,1,nlandp,'bats_internal:trnof')
    call getmem(tsw,1,nlandp,'bats_internal:tsw')
    call getmem(uaf,1,nlandp,'bats_internal:uaf')
    call getmem(usw,1,nlandp,'bats_internal:usw')
    call getmem(vegt,1,nlandp,'bats_internal:vegt')
    call getmem(vpdc,1,nlandp,'bats_internal:vpdc')
    call getmem(vspda,1,nlandp,'bats_internal:vspda')
    call getmem(vsw,1,nlandp,'bats_internal:vsw')
    call getmem(watr,1,nlandp,'bats_internal:watr')
    call getmem(watt,1,nlandp,'bats_internal:watt')
    call getmem(watu,1,nlandp,'bats_internal:watu')
    call getmem(wflux1,1,nlandp,'bats_internal:wflux1')
    call getmem(wflux2,1,nlandp,'bats_internal:wflux2')
    call getmem(wfluxc,1,nlandp,'bats_internal:wfluxc')
    call getmem(wiltr,1,nlandp,'bats_internal:wiltr')
    call getmem(wt,1,nlandp,'bats_internal:wt')
    call getmem(wta0,1,nlandp,'bats_internal:wta0')
    call getmem(wta,1,nlandp,'bats_internal:wta')
    call getmem(wtaq0,1,nlandp,'bats_internal:wtaq0')
    call getmem(wtg0,1,nlandp,'bats_internal:wtg0')
    call getmem(wtg,1,nlandp,'bats_internal:wtg')
    call getmem(wtg2,1,nlandp,'bats_internal:wtg2')
    call getmem(wtga,1,nlandp,'bats_internal:wtga')
    call getmem(wtgaq,1,nlandp,'bats_internal:wtgaq')
    call getmem(wtgl,1,nlandp,'bats_internal:wtgl')
    call getmem(wtglq,1,nlandp,'bats_internal:wtglq')
    call getmem(wtgq0,1,nlandp,'bats_internal:wtgq0')
    call getmem(wtgq,1,nlandp,'bats_internal:wtgq')
    call getmem(wtl0,1,nlandp,'bats_internal:wtl0')
    call getmem(wtlh,1,nlandp,'bats_internal:wtlh')
    call getmem(wtlq0,1,nlandp,'bats_internal:wtlq0')
    call getmem(wtlq,1,nlandp,'bats_internal:wtlq')
    call getmem(wtshi,1,nlandp,'bats_internal:wtshi')
    call getmem(wtsqi,1,nlandp,'bats_internal:wtsqi')
    call getmem(xkmx,1,nlandp,'bats_internal:xkmx')
    call getmem(xlai,1,nlandp,'bats_internal:xlai')
    call getmem(xlsai,1,nlandp,'bats_internal:xlsai')
    call getmem(z10fra,1,nlandp,'bats_internal:z10fra')
    call getmem(z1log,1,nlandp,'bats_internal:z1log')
    call getmem(z2fra,1,nlandp,'bats_internal:z2fra')
    call getmem(zh,1,nlandp,'bats_internal:zh')
    call getmem(zo,1,nlandp,'bats_internal:zo')
    call getmem(zlgdis,1,nlandp,'bats_internal:zlgdis')
    call getmem(zlglnd,1,nlandp,'bats_internal:zlglnd')
    call getmem(zlgsno,1,nlandp,'bats_internal:zlgsno')
    call getmem(zlgveg,1,nlandp,'bats_internal:zlgveg')
    call getmem(lveg,1,nlandp,'bats_internal:lveg')
    call getmem(ltex,1,nlandp,'bats_internal:ltex')
    call getmem(p0,1,nlandp,'bats_internal:p0')
    call getmem(qs0,1,nlandp,'bats_internal:qs0')
    call getmem(ts0,1,nlandp,'bats_internal:ts0')
    call getmem(swd0,1,nlandp,'bats_internal:swd0')
    call getmem(swf0,1,nlandp,'bats_internal:swf0')
    call getmem(ncpr0,1,nlandp,'bats_internal:ncpr0')
    call getmem(cpr0,1,nlandp,'bats_internal:cpr0')
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
