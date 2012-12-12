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

module mod_cu_interface
!
! Link atmospheric model and cumulus schemes
!
  use mod_constants
  use mod_realkinds
  use mod_cu_common
  use mod_cu_tiedtke
  use mod_cu_tables
  use mod_cu_bm
  use mod_cu_em
  use mod_cu_kuo
  use mod_cu_grell
  use mod_atm_interface , only : atmstate , slice , surfstate , domain
!
  contains

  subroutine init_cuscheme(mddom,atm1,aten,atms,chiten,sfs,qdot,pptc, &
                           ldmsk,cldfra,cldlwc,ktrop)
    implicit none
    type(domain) , intent(in) :: mddom
    type(atmstate) , intent(in) :: atm1 , aten
    type(slice) , intent(in) :: atms
    real(rk8) , pointer , intent(in) , dimension(:,:,:,:) :: chiten
    type(surfstate) , intent(in) :: sfs
    real(rk8) , pointer , intent(in) , dimension(:,:,:) :: qdot
    real(rk8) , pointer , intent(in) , dimension(:,:) :: pptc
    integer(ik4) , pointer , intent(in) , dimension(:,:) :: ldmsk , ktrop
    real(rk8) , pointer , dimension(:,:,:) :: cldlwc , cldfra
!
    call assignpnt(mddom%ht,sfhgt)
    call assignpnt(atm1%t,ptatm)
    call assignpnt(atm1%u,puatm)
    call assignpnt(atm1%v,pvatm)
    call assignpnt(atm1%qx,pvqxtm)
    call assignpnt(atms%tb3d,tas)
    call assignpnt(atms%ubx3d,uas)
    call assignpnt(atms%vbx3d,vas)
    call assignpnt(atms%pb3d,pas)
    call assignpnt(atms%qsb3d,qsas)
    call assignpnt(atms%qxb3d,qxas)
    call assignpnt(atms%chib3d,chias)
    call assignpnt(atms%za,hgt)
    call assignpnt(aten%t,tten)
    call assignpnt(aten%u,uten)
    call assignpnt(aten%v,vten)
    call assignpnt(aten%qx,qxten)
    call assignpnt(chiten,tchiten)
    call assignpnt(sfs%psa,psfcps)
    call assignpnt(sfs%psb,sfcps)
    call assignpnt(sfs%rainc,rainc)
    call assignpnt(sfs%qfx,qfx)
    call assignpnt(qdot,svv)
    call assignpnt(pptc,lmpcpc)
    call assignpnt(ldmsk,lmask)
    call assignpnt(cldfra,rcldfra)
    call assignpnt(cldlwc,rcldlwc)
    call assignpnt(ktrop,rktrop)
  end subroutine init_cuscheme
!
end module mod_cu_interface
