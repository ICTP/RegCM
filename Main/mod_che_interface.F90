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

module mod_che_interface
!
  use mod_realkinds
  use mod_regcm_types
  use mod_che_common
  use mod_che_cumtran
  use mod_che_dust
  use mod_che_indices
  use mod_che_mppio
  use mod_che_ncio
  use mod_che_param
  use mod_che_drydep
  use mod_che_bdyco
  use mod_che_emission
  use mod_che_carbonaer
  use mod_che_species
  use mod_che_tend
  use mod_che_start
!
  public
!
  contains 
!
#if (defined CLM)
  subroutine init_chem(atms,mddom,sfs,xpsb,ba_cr,fcc,cldfra,          &
                       rembc,remrat,coszrs,svegfrac2d,sxlai2d,        &
                       sfracv2d,sfracb2d,sfracs2d,solis,sdelt,     &
                       sdelq,ssw2da,convpr,icutop,icubot,taucldsp, &
                       voc_em,voc_em1,voc_em2,dep_vels)
#else
  subroutine init_chem(atms,mddom,sfs,xpsb,ba_cr,fcc,cldfra,      &
                       rembc,remrat,coszrs,svegfrac2d,sxlai2d,    &
                       sfracv2d,sfracb2d,sfracs2d,solis,sdelt, &
                       sdelq,ssw2da,convpr,icutop,icubot,      &
                       taucldsp)
#endif

    ! this routine define the pointer interface between the chem module and
    ! the rest of the model
    ! It also call startchem which is the chemistry initialisation routine

    implicit none

    real(rk8), pointer, dimension(:,:,:),intent(in) :: fcc
    real(rk8), pointer, dimension(:,:) :: svegfrac2d , solis , sdelt , &
             sdelq , ssw2da , sfracv2d , sfracb2d , sfracs2d , sxlai2d
    real(rk8), pointer, dimension(:,:,:) :: cldfra , rembc , remrat , convpr
    real(rk8), pointer, dimension(:,:,:,:) :: taucldsp
    integer(ik4) , pointer , dimension(:,:)  :: icutop , icubot
    type(slice) , intent(in) :: atms
    type(domain) , intent(in) :: mddom
    type(surfstate) , intent(in) :: sfs
    type(v2dbound) , intent(in) :: xpsb
    type(bound_area) , intent(in) :: ba_cr
    real(rk8) , pointer , dimension(:,:) :: coszrs

#if (defined CLM)
    real(rk8), pointer :: voc_em(:,:), voc_em1(:,:), voc_em2(:,:), dep_vels(:,:,:)
#endif

    call assignpnt(icutop,kcumtop)
    call assignpnt(icubot,kcumbot)
    call assignpnt(atms%tb3d,ctb3d)
    call assignpnt(atms%qxb3d,cqxb3d)
    call assignpnt(atms%rhob3d,crhob3d)
    ! wind at cell center
    call assignpnt(atms%ubx3d,cubx3d)
    call assignpnt(atms%vbx3d,cvbx3d)
    call assignpnt(atms%chib3d,chib3d)
    call assignpnt(atms%rhb3d,crhb3d)
!   
    call assignpnt(mddom%lndcat,clndcat)
    call assignpnt(mddom%xlat,cxlat)
    call assignpnt(mddom%ht,cht)
    call assignpnt(mddom%iveg,cveg2d) 
    call assignpnt(sfs%psb,cpsb)
    call assignpnt(xpsb%b0,psbb0)
    call assignpnt(xpsb%b1,psbb1)
    call assignpnt(sfs%tgb,ctg)
    call assignpnt(sfs%uvdrag,cuvdrag)
    call assignpnt(convpr,cconvpr)
    call assignpnt(fcc,cfcc)
    call assignpnt(cldfra,ccldfra)
    call assignpnt(rembc,crembc)
    call assignpnt(remrat,cremrat)
    call assignpnt(solis,csol2d)
    call assignpnt(svegfrac2d,cvegfrac)
    call assignpnt(sfracv2d,csfracv2d)
    call assignpnt(sfracb2d,csfracb2d)
    call assignpnt(sfracs2d,csfracs2d)
    call assignpnt(sxlai2d,cxlai2d)

    call assignpnt(sdelt,csdeltk2d) 
    call assignpnt(sdelq,csdelqk2d)

    call assignpnt(atms%za,cza)
    call assignpnt(atms%dzq,cdzq)
    call assignpnt(coszrs,czen)
    call assignpnt(ssw2da,cssw2da)   

    call assignpnt(taucldsp, ctaucld)

    cba%havebound = ba_cr%havebound
    call assignpnt(ba_cr%bsouth, cba%bsouth)   
    call assignpnt(ba_cr%bnorth, cba%bnorth)   
    call assignpnt(ba_cr%beast, cba%beast)   
    call assignpnt(ba_cr%bwest, cba%bwest)   
    call assignpnt(ba_cr%ibnd, cba%ibnd)   

    cba%ns = ba_cr%ns
    cba%nn = ba_cr%nn
    cba%ne = ba_cr%ne
    cba%nw = ba_cr%nw
    cba%nsp = ba_cr%nsp

#if (defined CLM)
#if (defined VOC)
    call assignpnt(voc_em,cvoc_em)
    call assignpnt(voc_em1,cvoc_em1)
    call assignpnt(voc_em2,cvoc_em2)
#endif
    call assignpnt(dep_vels,cdep_vels)
#endif

  end subroutine init_chem
!
end module mod_che_interface
