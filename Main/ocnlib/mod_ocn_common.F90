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
  use mod_constants
  use mod_regcm_types
  use mod_ocn_internal
  use mod_ocn_bats
  use mod_ocn_coare
  use mod_ocn_lake
  use mod_ocn_zeng
  use mod_ocn_albedo
  use mod_runparams , only : rcmtimer , syncro_srf , iwavcpl
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
    if ( lseaice ) call seaice
    call ocn_interf(lm,lms,2)
  end subroutine vecocn

  subroutine initocn(lm,lms)
    implicit none
    type(lm_exchange) , intent(inout) :: lm
    type(lm_state) , intent(inout) :: lms
    integer(ik4) :: n
    ! Set up Masks for lake and sea ice
    if ( lcoup ) then
      call c2l_gs(ocncomm,lm%icplmsk,icpl)
    end if
    call c2l_ss(ocncomm,lm%xlat1,lat)
    call c2l_ss(ocncomm,lm%iveg1,omask)
    call c2l_gs(ocncomm,lm%zencos,czenith)
    if ( llake .or. lseaice ) then
      call c2l_ss(ocncomm,lm%ldmsk1,mask)
      if ( llake ) then
        call c2l_ss(ocncomm,lm%dhlake1,dhlake)
        ilake(:) = 0
        do n = iocnbeg , iocnend
          if ( omask(n) == 14 ) then
            ilake(n) = 1
            if ( mask(n) == 2 ) then
              mask(n) = 4
            else
              mask(n) = 3
            end if
          else
            if ( mask(n) == 2 ) then
              mask(n) = 2
            else
              mask(n) = 1
            end if
          end if
        end do
        call allocate_mod_ocn_lake
      else
        do n = iocnbeg , iocnend
          if ( mask(n) == 2 ) then
            mask(n) = 2
          else
            mask(n) = 1
          end if
        end do
      end if
    else
      mask(:) = 1
    end if
    if ( llake .or. lseaice ) then
      call c2l_ss(ocncomm,lms%sfice,sfice)
      call c2l_ss(ocncomm,lms%sncv,sncv)
      call c2l_ss(ocncomm,lms%snag,snag)
    end if
    if ( rcmtimer%start( ) ) then
      call c2l_gs(ocncomm,lm%tg,tgb)
      tgrd = tgb
      tgbrd = tgb
      t2m = tgb
      um10 = 1.0_rkx ! Assume a mean of 1m/s wind for init.
    else
      call c2l_ss(ocncomm,lms%tgrd,tgrd)
      call c2l_ss(ocncomm,lms%tgbrd,tgbrd)
      call c2l_ss(ocncomm,lms%um10,um10)
      call c2l_ss(ocncomm,lms%tlef,t2m)
      call c2l_gs(ocncomm,lm%qfx,evpr)
      call c2l_gs(ocncomm,lm%hfx,sent)
      if ( ldcsst ) then
        call c2l_ss(ocncomm,lms%sst,sst)
        call c2l_ss(ocncomm,lms%deltas,deltas)
        call c2l_ss(ocncomm,lms%tdeltas,tdeltas)
        call c2l_ss(ocncomm,lms%tskin,tskin)
      end if
      if ( llake ) then
        call c2l_ss(ocncomm,lms%eta,laketa)
        call c2l_ss(ocncomm,lms%hi,lakhi)
        call c2l_ss(ocncomm,lms%tlake,laktlake)
        call lake_fillvar(var_eta,laketa,1)
        call lake_fillvar(var_hi,lakhi,1)
        call lake_fillvar(var_tlak,laktlake,1)
      end if
      if ( llake .or. lseaice ) then
        call c2l_ss(ocncomm,lm%ldmsk1,mask)
        if ( llake ) then
          do n = iocnbeg , iocnend
            if ( ilake(n) == 1 ) then
              if ( mask(n) == 2 ) then
                mask(n) = 4
              else
                mask(n) = 3
              end if
            else
              if ( mask(n) == 2 ) then
                mask(n) = 2
              else
                mask(n) = 1
              end if
            end if
          end do
        else
          do n = iocnbeg , iocnend
            if ( mask(n) == 2 ) then
              mask(n) = 2
            else
              mask(n) = 1
            end if
          end do
        end if
      end if
    end if
    if ( llake ) then
      call initlake
      lms%lakmsk = .false.
      call l2c_ss(ocncomm,lakmsk,lms%lakmsk)
    end if
    emiss(:) = ocn_sfcemiss
    call l2c_ss(ocncomm,emiss,lms%emisv)
  end subroutine initocn

  subroutine ocn_interf(lm,lms,ivers)
    implicit none
    type(lm_exchange) , intent(inout) :: lm
    type(lm_state) , intent(inout) :: lms
    integer(ik4) , intent(in) :: ivers
    integer(ik4) :: i , j , i1 , i2 , n
    if ( ivers == 1 ) then
      ! RegCM -> OCN
      if ( llake .or. lseaice ) then
        call c2l_ss(ocncomm,lm%ldmsk1,mask)
        if ( llake ) then
          do n = iocnbeg , iocnend
            if ( ilake(n) == 1 ) then
              if ( mask(n) == 2 ) then
                mask(n) = 4
              else
                mask(n) = 3
              end if
            else
              if ( mask(n) == 2 ) then
                mask(n) = 2
              else
                mask(n) = 1
              end if
            end if
          end do
        else
          do n = iocnbeg , iocnend
            if ( mask(n) == 2 ) then
              mask(n) = 2
            else
              mask(n) = 1
            end if
          end do
        end if
      end if
      if ( lseaice .or. llake ) then
        call c2l_ss(ocncomm,lms%sfice,sfice)
      end if
      call c2l_gs(ocncomm,lm%hgt,ht)
      call c2l_gs(ocncomm,lm%uatm,usw)
      call c2l_gs(ocncomm,lm%vatm,vsw)
      usw = sign(1.0_rkx,usw)*max(abs(usw),0.001_rkx)
      vsw = sign(1.0_rkx,vsw)*max(abs(vsw),0.001_rkx)
      call c2l_gs(ocncomm,lm%tatm,tatm)
      call c2l_gs(ocncomm,lm%patm,patm)
      call c2l_gs(ocncomm,lm%sfta,sfta)
      call c2l_ss(ocncomm,lms%tlef,t2m)
      call c2l_gs(ocncomm,lm%tg,tgb)
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
      prcp = (cprate+ncprate) * syncro_srf%rw
      if ( ldcsst ) then
        call c2l_ss(ocncomm,lms%sst,sst)
        call c2l_ss(ocncomm,lms%deltas,deltas)
        call c2l_ss(ocncomm,lms%tdeltas,tdeltas)
        call c2l_ss(ocncomm,lms%tskin,tskin)
      end if
      if ( iwavcpl == 1 ) then
        call c2l_gs(ocncomm,lm%zo,zoo)
        call c2l_gs(ocncomm,lm%ustar,ustr)
      end if
    else
      ! OCN -> RegCM
      call l2c_ss(ocncomm,tgb,lms%tgbb)
      call l2c_ss(ocncomm,tgrd,lms%tgrd)
      call l2c_ss(ocncomm,tgbrd,lms%tgbrd)
      call l2c_ss(ocncomm,t2m,lms%tlef)
      call l2c_ss(ocncomm,sent,lms%sent)
      call l2c_ss(ocncomm,evpr,lms%evpr)
      call l2c_ss(ocncomm,drag,lms%drag)
      call l2c_ss(ocncomm,u10m,lms%u10m)
      call l2c_ss(ocncomm,v10m,lms%v10m)
      call l2c_ss(ocncomm,ustr,lms%ustar)
      call l2c_ss(ocncomm,zoo,lms%zo)
      call l2c_ss(ocncomm,rhoa,lms%rhoa)
      call l2c_ss(ocncomm,taux,lms%taux)
      call l2c_ss(ocncomm,tauy,lms%tauy)
      call l2c_ss(ocncomm,t2m,lms%t2m)
      call l2c_ss(ocncomm,q2m,lms%q2m)
      call l2c_ss(ocncomm,sfps,lms%sfcp)
      call l2c_ss(ocncomm,prcp,lms%prcp)
      call l2c_ss(ocncomm,um10,lms%um10)
      call l2c_ss(ocncomm,ram1,lms%ram1)
      call l2c_ss(ocncomm,rah1,lms%rah1)
      call l2c_ss(ocncomm,br,lms%br)
      if ( ldcsst ) then
        call l2c_ss(ocncomm,deltas,lms%deltas)
        call l2c_ss(ocncomm,tdeltas,lms%tdeltas)
        call l2c_ss(ocncomm,tskin,lms%tskin)
      end if
      if ( llake .or. lseaice ) then
        do n = iocnbeg , iocnend
          if ( mask(n) == 2 .or. mask(n) == 4 ) then
            mask(n) = 2
          else
            mask(n) = 0
          end if
        end do
        call l2c_ss(ocncomm,mask,lm%ldmsk1)
        if ( llake ) then
          call lake_fillvar(var_eta,laketa,0)
          call lake_fillvar(var_hi,lakhi,0)
          call lake_fillvar(var_tlak,laktlake,0)
          call l2c_ss(ocncomm,laketa,lms%eta)
          call l2c_ss(ocncomm,lakhi,lms%hi)
          call l2c_ss(ocncomm,laktlake,lms%tlake)
        end if
        call l2c_ss(ocncomm,sfice,lms%sfice)
        call l2c_ss(ocncomm,snag,lms%snag)
        call l2c_ss(ocncomm,sncv,lms%sncv)
        call l2c_ss(ocncomm,sm,lms%snwm)
        do i = ici1 , ici2
          do j = jci1 , jci2
            if ( lm%ldmsk(j,i) /= 1 ) then
              i1 = count(lm%ldmsk1(:,j,i) == 0)
              i2 = count(lm%ldmsk1(:,j,i) == 2)
              if ( i1 > i2 ) then
                lm%ldmsk(j,i) = 0
              else
                lm%ldmsk(j,i) = 2
              end if
            end if
          end do
        end do
      end if
    end if
  end subroutine ocn_interf

  subroutine albedoocn(lm,lms)
    implicit none
    type(lm_exchange) , intent(inout) :: lm
    type(lm_state) , intent(inout) :: lms
    integer(ik4) :: i , j , n
    call ocn_albedo
    call l2c_ss(ocncomm,swdiral,lms%swdiralb)
    call l2c_ss(ocncomm,lwdiral,lms%lwdiralb)
    call l2c_ss(ocncomm,swdifal,lms%swdifalb)
    call l2c_ss(ocncomm,lwdifal,lms%lwdifalb)
    do i = ici1 , ici2
      do j = jci1 , jci2
        do n = 1 , nnsg
          if ( lm%ldmsk1(n,j,i) /= 1 ) then
            lms%swalb(n,j,i) = max(lms%swdiralb(n,j,i),lms%swdifalb(n,j,i))
            lms%lwalb(n,j,i) = max(lms%lwdiralb(n,j,i),lms%lwdifalb(n,j,i))
          end if
        end do
      end do
    end do
  end subroutine albedoocn

end module mod_ocn_common
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
