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

module mod_ocn_common
!
! Storage for Surface (BATS and shared by CLM) variables
!
  use mod_intkinds
  use mod_realkinds
  use mod_regcm_types
  use mod_ocn_internal
  use mod_ocn_bats
  use mod_ocn_coare
  use mod_ocn_lake
  use mod_ocn_zeng
  use mod_ocn_albedo
  use mod_runparams , only : ktau , iemiss , rtsrf
  use mod_mppparam

  implicit none

  private

  logical :: lcoup

  public :: initocn , vecocn , albedoocn
  public :: llake , ldcsst , lseaice , lcoup
  public :: allocate_mod_ocn_internal

  contains

  subroutine vecocn(lm,lms)
    implicit none
    type(lm_exchange) , intent(inout) :: lm
    type(lm_state) , intent(inout) :: lms
    call ocn_interf(lm,lms,1)
    select case ( iocnflx )
      case (1)
        call ocnbats
      case (2)
        call zengocndrv
      case (3)
        call coare3_drv
    end select
    if ( llake ) call lakedrv
    call seaice
    call ocn_interf(lm,lms,2)
    call ocn_albedo
  end subroutine vecocn

  subroutine initocn(lm,lms)
    implicit none
    type(lm_exchange) , intent(inout) :: lm
    type(lm_state) , intent(inout) :: lms
    ! Set up Masks for lake and sea ice
    if ( lcoup ) then
      call c2l_gs(ocncomm,lm%icplmsk,icpl)
    end if
    if ( llake .or. lseaice ) then
      call c2l_ss(ocncomm,lm%ldmsk1,xmask)
      where ( xmask == 2 )
        mask = 2
      else where
        mask = 1
      end where
      if ( llake ) then
        call c2l_ss(ocncomm,lm%dhlake1,dhlake)
        call c2l_ss(ocncomm,lm%iveg1,ilake)
        where ( ilake == 14 )
          ilake = 1
        else where
          ilake = 0
        end where
        where ( ilake == 1 )
          mask = mask + 2
        end where
        call allocate_mod_ocn_lake
      end if
    else
      mask = 1
    end if
    if ( ktau == 0 ) then
      call c2l_ss(ocncomm,lm%xlat1,lat)
      call c2l_gs(ocncomm,lm%tground2,tgb)
      tgrd = tgb
      tgbrd = tgb
      if ( llake .or. lseaice ) then
        call c2l_gs(ocncomm,lm%snowam,sncv)
      end if
      if ( iemiss == 1 ) then
        where ( mask == 1 .or. mask == 3 )
          emiss = 0.955D0
        else where
          emiss = 0.97D0
        end where
      else
        emiss = 0.9995D0
      end if
      call l2c_ss(ocncomm,emiss,lms%emisv)
    else
      call c2l_ss(ocncomm,lm%xlat1,lat)
      call c2l_ss(ocncomm,lms%tgrd,tgrd)
      call c2l_ss(ocncomm,lms%tgbrd,tgbrd)
      call c2l_ss(ocncomm,lms%emisv,emiss)
      call c2l_gs(ocncomm,lm%qfx,evpr)
      call c2l_gs(ocncomm,lm%hfx,sent)
      if ( ldcsst ) then
        call c2l_ss(ocncomm,lms%deltas,deltas)
        call c2l_ss(ocncomm,lms%tdeltas,tdeltas)
        call c2l_ss(ocncomm,lms%tskin,tskin)
      end if
      if ( llake .or. lseaice ) then
        call c2l_ss(ocncomm,lms%sfice,sfice)
        call c2l_ss(ocncomm,lms%snag,snag)
        call c2l_ss(ocncomm,lms%sncv,sncv)
        call c2l_ss(ocncomm,lms%scvk,scvk)
      end if
      if ( llake ) then
        call c2l_ss(ocncomm,lms%eta,laketa)
        call c2l_ss(ocncomm,lms%hi,lakhi)
        call c2l_ss(ocncomm,lms%tlake,laktlake)
        call lake_fillvar(var_eta,laketa,1)
        call lake_fillvar(var_hi,lakhi,1)
        call lake_fillvar(var_tlak,laktlake,1)
      end if
    end if
    if ( llake ) then
      call initlake
      call l2c_ss(ocncomm,lakmsk,lms%lakmsk)
    end if
  end subroutine initocn

  subroutine ocn_interf(lm,lms,ivers)
    implicit none
    type(lm_exchange) , intent(inout) :: lm
    type(lm_state) , intent(inout) :: lms
    integer , intent(in) :: ivers
    if ( ivers == 1 ) then
      ! RegCM -> OCN
      if ( llake .or. lseaice ) then
        call c2l_ss(ocncomm,lm%ldmsk1,xmask)
        where ( xmask == 2 )
          mask = 2
        else where
          mask = 1
        end where
        call c2l_ss(ocncomm,lms%sfice,sfice)
      end if
      call c2l_gs(ocncomm,lm%hgt,ht)
      call c2l_gs(ocncomm,lm%uatm,usw)
      call c2l_gs(ocncomm,lm%vatm,vsw)
      call c2l_gs(ocncomm,lm%tatm,sts)
      call c2l_gs(ocncomm,lm%tground2,tgb)
      call c2l_gs(ocncomm,lm%hpbl,hpbl)
      call c2l_gs(ocncomm,lm%qvatm,qv)
      call c2l_gs(ocncomm,lm%rhox,rhox)
      call c2l_gs(ocncomm,lm%rlwf,rlwf)
      call c2l_gs(ocncomm,lm%rswf,rswf)
      call c2l_gs(ocncomm,lm%dwrlwf,dwrlwf)
      call c2l_gs(ocncomm,lm%zencos,czenith)
      call c2l_gs(ocncomm,lm%sfps,sfps)
      call c2l_gs(ocncomm,lm%hfx,sent)
      call c2l_gs(ocncomm,lm%qfx,evpr)
      call c2l_gs(ocncomm,lm%cprate,cprate)
      call c2l_gs(ocncomm,lm%ncprate,ncprate)
      prcp = (cprate+ncprate) * rtsrf
      sfps = (sfps+ptop)*d_1000
      if ( ldcsst ) then
        call c2l_ss(ocncomm,lms%deltas,deltas)
        call c2l_ss(ocncomm,lms%tdeltas,tdeltas)
        call c2l_ss(ocncomm,lms%tskin,tskin)
        call c2l_ss(ocncomm,lms%sst,sst)
      end if
      if ( llake ) then
        where ( ilake == 1 )
          mask = mask + 2
        end where
      end if
    else
      ! OCN -> RegCM
      call l2c_ss(ocncomm,tgb,lms%tgbb)
      call l2c_ss(ocncomm,tgrd,lms%tgrd)
      call l2c_ss(ocncomm,tgbrd,lms%tgbrd)
      call l2c_ss(ocncomm,tgrd,lms%tlef)
      call l2c_ss(ocncomm,sent,lms%sent)
      call l2c_ss(ocncomm,evpr,lms%evpr)
      call l2c_ss(ocncomm,drag,lms%drag)
      call l2c_ss(ocncomm,u10m,lms%u10m)
      call l2c_ss(ocncomm,v10m,lms%v10m)
      call l2c_ss(ocncomm,taux,lms%taux)
      call l2c_ss(ocncomm,tauy,lms%tauy)
      call l2c_ss(ocncomm,t2m,lms%t2m)
      call l2c_ss(ocncomm,q2m,lms%q2m)
      call l2c_ss(ocncomm,sfps,lms%sfcp)
      call l2c_ss(ocncomm,prcp,lms%prcp)
      if ( ldcsst ) then
        call l2c_ss(ocncomm,deltas,lms%deltas)
        call l2c_ss(ocncomm,tdeltas,lms%tdeltas)
        call l2c_ss(ocncomm,tskin,lms%tskin)
      end if
      if ( llake .or. lseaice ) then
        xmask = mask
        where ( xmask == 1 .or. xmask == 3 ) xmask = 0
        if ( llake ) then
          call lake_fillvar(var_eta,laketa,0)
          call lake_fillvar(var_hi,lakhi,0)
          call lake_fillvar(var_tlak,laktlake,0)
          call l2c_ss(ocncomm,laketa,lms%eta)
          call l2c_ss(ocncomm,lakhi,lms%hi)
          call l2c_ss(ocncomm,laktlake,lms%tlake)
          where ( xmask == 4 ) xmask = 2
        end if
        call l2c_ss(ocncomm,xmask,lm%ldmsk1)
        call l2c_ss(ocncomm,sfice,lms%sfice)
        call l2c_ss(ocncomm,snag,lms%snag)
        call l2c_ss(ocncomm,sncv,lms%sncv)
        call l2c_ss(ocncomm,scvk,lms%scvk)
        call l2c_ss(ocncomm,sm,lms%snwm)
        ! Emissivity for surface ice
        if ( iemiss == 1 ) then
          where ( xmask == 0 )
            emiss = 0.955D0
          else where
            emiss = 0.97D0
          end where
          call l2c_ss(ocncomm,emiss,lms%emisv)
        end if
      end if
    end if
  end subroutine ocn_interf

  subroutine albedoocn(lm,lms)
    implicit none
    type(lm_exchange) , intent(inout) :: lm
    type(lm_state) , intent(inout) :: lms
    call ocn_interf(lm,lms,1)
    call ocn_albedo
    call l2c_ss(ocncomm,swdiral,lms%swdiralb)
    call l2c_ss(ocncomm,lwdiral,lms%lwdiralb)
    call l2c_ss(ocncomm,swdifal,lms%swdifalb)
    call l2c_ss(ocncomm,lwdifal,lms%lwdifalb)
  end subroutine albedoocn

end module mod_ocn_common
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
