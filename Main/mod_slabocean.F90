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
 
module mod_slabocean
!
!
  use mod_runparams
  use mod_memutil
  use mod_atm_interface , only :  surfstate , ts1 , domain
  use mod_mppparam
  use mod_constants
  use mod_outvars
!
  private
!
  real(rk8) :: mixed_layer_salin = d_100/d_three
  integer(ik4) , parameter :: monperyear = 12
!
  ! the actual prognotic sst pointing on tg2
  real(rk8) , pointer , dimension(:,:) :: sstemp
  real(ik8) , pointer , dimension(:,:) :: ohfx , oqfx , ofsw , oflw
  real(rk8) , pointer , dimension(:,:) :: olndcat
  real(rk8) , pointer , dimension(:,:,:) :: qflux_restore_sst
  integer(ik4) , dimension(monperyear) , target :: stepcount
  real(rk8) , pointer , dimension (:,:) :: qflux_sst , qflux_adj , net_hflx , &
    hflx , qflb0 , qflb1 , qflbt
  integer(ik4) , pointer , dimension (:,:) :: ocmask

  real(rk8) :: dtocean
!
  public :: allocate_mod_slabocean , init_slabocean , update_slabocean
  public :: fill_slaboc_outvars
  public :: qflb0 , qflb1 , qflbt , qflux_restore_sst , stepcount

  contains
!
    subroutine allocate_mod_slabocean
      implicit none
      call getmem3d(qflux_restore_sst,jci1,jci2,ici1,ici2,1,monperyear,&
                    'slab_ocean:qflux_restore_sst')
      call getmem2d(qflux_sst,jci1,jci2,ici1,ici2,'slab_ocean:qflux_sst')
      call getmem2d(qflux_adj,jci1,jci2,ici1,ici2,'slab_ocean:qflux_adj')
      call getmem2d(qflb0,jci1,jci2,ici1,ici2,'slab_ocean:qflb0')
      call getmem2d(qflb1,jci1,jci2,ici1,ici2,'slab_ocean:qflb1')
      call getmem2d(qflbt,jci1,jci2,ici1,ici2,'slab_ocean:qflbt')
      call getmem2d(net_hflx,jci1,jci2,ici1,ici2,'slab_ocean:net_hflx')
      call getmem2d(hflx,jci1,jci2,ici1,ici2,'slab_ocean:hflx')
      stepcount(:) = 0
    end subroutine allocate_mod_slabocean
!
    subroutine init_slabocean(sfs,ldmsk,fsw,flw)
      implicit none
      ! interface for regcm variable / slab ocean
      type(surfstate) , intent(in) :: sfs
      integer(ik4) , pointer , intent(in) , dimension(:,:) :: ldmsk
      real(rk8) , pointer , intent(in) , dimension(:,:) :: fsw , flw

      call assignpnt(sfs%tgb,sstemp)
      call assignpnt(ldmsk,ocmask)
      call assignpnt(sfs%hfx,ohfx)
      call assignpnt(sfs%qfx,oqfx)
      call assignpnt(fsw,ofsw)
      call assignpnt(flw,oflw) 
      dtocean =  dtsrf
    end subroutine init_slabocean

    subroutine update_slabocean
      implicit none
      ! mlcp is the heat capacity of the mixed layer [J / m3 / deg C] * m
      real(rk8) :: mlcp
      integer(ik4) :: i , j
#ifdef DEBUG
      real(rk8) , dimension(5) :: pval , pval1
#endif
      ! water heat capacity ~ 4 J/g/K         
      mlcp = mixed_layer_depth*4.0D6

      if ( do_restore_sst ) then
        stepcount(xmonth) = stepcount(xmonth)+1
        do i = ici1 , ici2
          do j = jci1 , jci2
            if ( ocmask(j,i) == 0 ) then
              qflux_sst(j,i) = (ts1(j,i) - sstemp(j,i)) * &
                mlcp / (sst_restore_timescale * 86400.0D0) ! w/m2
              qflux_restore_sst(j,i,xmonth) = &
                qflux_restore_sst(j,i,xmonth) + qflux_sst(j,i)
            end if
          end do
        end do
      else if ( do_qflux_adj ) then
        !
        ! Find the current climatological heat flux adjustment (qflux_adj).  
        ! The qflux_adj can be added to the heat flux to find the net heat
        ! flux anomaly into the mixed layer for the coupled mixed layer expt.
        !
        ! Note that any ice model contributions need to be added in, which
        ! is handled in an off-line calculation to derive the final qflux
        ! adjustment. These include ice restoring and ice lid contributions.
        !
        qflux_adj = qflb0 + dtocean*qflbt
      end if    
      !
      ! energy budget in the mixed layer including the q flux therm
      !
#ifdef DEBUG
      if ( mod(ktau+1,krep) == 0 ) then
        pval(1) = maxval(ofsw)
        pval(2) = maxval(oflw)
        pval(3) = maxval(ohfx)
        pval(4) = maxval(2.26D6*oqfx)
        pval(5) = maxval(qflux_adj)
        do i = 1 , size(pval)
          call maxall(pval(i),pval1(i))
        end do
        if ( myid == italk ) then
          write(stdout,'(a,f12.5)') ' $$$ SLAB :  max fsw = ', pval1(1)
          write(stdout,'(a,f12.5)') ' $$$         max flw = ', pval1(2)
          write(stdout,'(a,f12.5)') ' $$$         max hfx = ', pval1(3)
          write(stdout,'(a,f12.5)') ' $$$         max qfx = ', pval1(4)
          write(stdout,'(a,f12.5)') ' $$$         max adj = ', pval1(5)
        end if
      end if
#endif
      do i = ici1 , ici2
        do j = jci1 , jci2
          if ( ocmask(j,i) == 0 ) then
            ! The following are some key equations for this model:
            ! flux from or to the atmosphere ( convention  = positive downward)
            ! multiply evaporation by latent heat of evaporation     
            hflx(j,i) = ofsw(j,i)-oflw(j,i)-ohfx(j,i)-wlhv*oqfx(j,i)
            ! account for the retaured or adjustment flux term
            net_hflx(j,i) = hflx(j,i) + qflux_adj(j,i) + qflux_sst(j,i)
            ! update the prognostic seasurface temperature (pointing on tg2)
            sstemp(j,i) = sstemp(j,i) + dtocean*net_hflx(j,i)/mlcp
          end if
        end do
      end do
    end subroutine update_slabocean

    subroutine fill_slaboc_outvars
      implicit none
      integer(ik4) :: imon
      if ( associated(slab_qflx_out) ) then
        do imon = 1 , monperyear
          if ( stepcount(imon) /= 0 ) then
            where ( ocmask == 0 )
              slab_qflx_out(:,:,imon) = &
                qflux_restore_sst(:,:,imon)/stepcount(imon)
            else where
              slab_qflx_out(:,:,imon) = dmissval
            end where
          else
            where ( ocmask == 0 )
              slab_qflx_out(:,:,imon) = d_zero
            else where
              slab_qflx_out(:,:,imon) = dmissval
            end where
          end if
        end do
      end if
    end subroutine fill_slaboc_outvars
!
end module mod_slabocean
