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
  use mod_atm_interface , only : slice , surfpstate, domain, surftstate
  use mod_che_common
  use mod_che_cumtran
  use mod_che_dust
  use mod_che_indices
  use mod_che_mppio
  use mod_che_ncio
  use mod_che_param
  use mod_che_drydep
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
  subroutine init_chem(ifrest, idirect,dt,rdxsq,chemfrq,dtrad,calday,dsigma,atms,&
                       sps2,mddom,sts2,fcc,cldfra,rembc,remrat,a,anudg,za,dzq,twt,&
                       ptop,coszrs,veg2d,svegfrac2d,solis,sdeltk2d,sdelqk2d,ssw2da,icutop,icubot)

! this routine define the pointer interface between the chem module and the rest of the model
! It also call startchem which is the chemistry initialisation routine

    implicit none
    logical, intent(in) :: ifrest
    integer , intent(in) :: idirect
    real(dp) , intent(in) :: dt , chemfrq , dtrad, calday,rdxsq

    real(dp) , pointer , dimension(:) , intent(in) :: dsigma ! dsigma
    real(dp), pointer, dimension(:,:,:),intent(in) :: fcc,za,dzq
    real(dp), pointer, dimension(:,:) :: rembc,svegfrac2d,solis,sdeltk2d,sdelqk2d,ssw2da,remrat,twt
    real(dp), pointer, dimension(:,:,:) :: cldfra
    integer , pointer , dimension(:,:) :: icutop , icubot, veg2d
    type(surfpstate) , intent(in) :: sps2
    type(slice) , intent(in) :: atms
    type(domain), intent(in):: mddom
    type(surftstate), intent(in) :: sts2
    real(dp) , pointer , dimension(:) :: a,anudg
    real(dp) , pointer , dimension(:,:) :: coszrs
    real(dp) :: ptop

    ichdir = idirect

    chfrq = chemfrq
    rafrq = dtrad
    dtche = dt
    crdxsq = rdxsq    
    ccalday=calday
    chptop = ptop



    call assignpnt(dsigma,cdsigma)
    call assignpnt(icutop,kcumtop)
    call assignpnt(icubot,kcumbot)
    call assignpnt(sps2%ps,cpsb)
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
    call assignpnt(sts2%tg,ctg)

    call assignpnt(fcc,cfcc)
    call assignpnt(cldfra,ccldfra)
    call assignpnt(rembc,crembc)
    call assignpnt(remrat,cremrat)
    call assignpnt(solis,csol2d)
    call assignpnt(svegfrac2d,cvegfrac)
    call assignpnt(sdeltk2d,csdeltk2d) 
    call  assignpnt(sdelqk2d,csdelqk2d)
call assignpnt(veg2d,cveg2d) 
!call assignpnt(rough,crough) 
!call  assignpnt(iexsol,ciexsol) 
!call assignpnt(xmopor,cxmopor) 
!call assignpnt(depuv,cdepuv) 
    call assignpnt(za,cza)
    call assignpnt(a,hlev)
    call assignpnt(anudg,canudg)
    call assignpnt(dzq,cdzq)
    call assignpnt(twt,ctwt)
    call assignpnt(coszrs,czen)
    call assignpnt(ssw2da,cssw2da)   
!!$
! Peform chemistry initialisation

    call start_chem(ifrest)



  end subroutine init_chem
!
end module mod_che_interface
