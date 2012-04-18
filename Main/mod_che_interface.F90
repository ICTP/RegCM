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

module mod_che_interface
!
  use mod_realkinds
  use mod_atm_interface , only : slice , domain, surfstate,bound_area
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
  subroutine init_chem(ifrest, idirect,dt,dx,chemfrq,dtrad,dsigma,atms, &
                       mddom,sfs,ba_cr, fcc,cldfra,rembc,remrat,a,anudg,   &
                       twt,ptop,coszrs,iveg,svegfrac2d,solis,sdeltk2d,     &
                       sdelqk2d,ssw2da,icutop,icubot)

! this routine define the pointer interface between the chem module and the rest of the model
! It also call startchem which is the chemistry initialisation routine

    implicit none
    logical, intent(in) :: ifrest
    integer , intent(in) :: idirect
    real(dp) , intent(in) :: dt , chemfrq , dtrad,dx 

    real(dp) , pointer , dimension(:) , intent(in) :: dsigma ! dsigma
    real(dp), pointer, dimension(:,:,:),intent(in) :: fcc
    real(dp), pointer, dimension(:,:) :: svegfrac2d , solis , sdeltk2d , &
                                         sdelqk2d , ssw2da , twt
    real(dp), pointer, dimension(:,:,:) :: cldfra , rembc , remrat
    integer , pointer , dimension(:,:) :: icutop , icubot, iveg
    type(slice) , intent(in) :: atms
    type(domain), intent(in):: mddom
    type (surfstate) , intent(in) :: sfs
    type(bound_area) , intent(in) :: ba_cr
    real(dp) , pointer , dimension(:) :: a , anudg
    real(dp) , pointer , dimension(:,:) :: coszrs
    real(dp) :: ptop

    ichdir = idirect

    chfrq = chemfrq
    rafrq = dtrad
    dtche = dt
    crdxsq = d_one/(dx*dx)    
    chptop = ptop

    call assignpnt(dsigma,cdsigma)
    call assignpnt(icutop,kcumtop)
    call assignpnt(icubot,kcumbot)
    call assignpnt(atms%tb3d,ctb3d)
    call assignpnt(atms%qvb3d,cqvb3d)
    call assignpnt(atms%qcb3d,cqcb3d)
    call assignpnt(atms%rhob3d,crhob3d)
!wind at cell center
    call assignpnt(atms%ubx3d,cubx3d)
    call assignpnt(atms%vbx3d,cvbx3d)
    call assignpnt(atms%chib3d,chib3d)
!   
    call assignpnt(mddom%lndcat,clndcat)
    call assignpnt(mddom%ht,cht)
    call assignpnt(sfs%psb,cpsb)
    call assignpnt(sfs%tgb,ctg)
    call assignpnt(sfs%uvdrag,cuvdrag)
    
    call assignpnt(fcc,cfcc)
    call assignpnt(cldfra,ccldfra)
    call assignpnt(rembc,crembc)
    call assignpnt(remrat,cremrat)
    call assignpnt(solis,csol2d)
    call assignpnt(svegfrac2d,cvegfrac)
    call assignpnt(sdeltk2d,csdeltk2d) 
    call assignpnt(sdelqk2d,csdelqk2d)
    call assignpnt(iveg,cveg2d) 
!call assignpnt(rough,crough) 
!call assignpnt(iexsol,ciexsol) 
!call assignpnt(xmopor,cxmopor) 
!call assignpnt(depuv,cdepuv) 
    call assignpnt(atms%za,cza)
    call assignpnt(atms%dzq,cdzq)
    call assignpnt(a,hlev)
    call assignpnt(anudg,canudg)
    call assignpnt(twt,ctwt)
    call assignpnt(coszrs,czen)
    call assignpnt(ssw2da,cssw2da)   

!!$!
  cba_cr%dotflag = ba_cr%dotflag
  cba_cr%havebound = ba_cr%havebound
 call assignpnt(ba_cr%bsouth, cba_cr%bsouth)   
 call assignpnt(ba_cr%bnorth, cba_cr%bnorth)   
 call assignpnt(ba_cr%beast, cba_cr%beast)   
 call assignpnt(ba_cr%bwest, cba_cr%bwest)   
 call assignpnt(ba_cr%ibnd, cba_cr%ibnd)   

 cba_cr%ns = ba_cr%ns
 cba_cr%nn = ba_cr%nn
 cba_cr%ne = ba_cr%ne
 cba_cr%nw = ba_cr%nw
 cba_cr%nsp = ba_cr%nsp
!!$


!$
! Peform chemistry initialisation

!   call start_chem(ifrest)

  end subroutine init_chem
!
end module mod_che_interface
