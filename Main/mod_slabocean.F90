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
  use mod_atm_interface , only :  surfstate, ts1,domain
  use mod_mppparam
  use mod_constants
  use mod_outvars


!
  private
!
  real(rk8)    :: mixed_layer_depth 
  real (rk8)   :: mixed_layer_salin = 33.333
  integer(ik4) :: islab_ocean
  logical :: do_qflux_adj         ! TK mod
  logical :: do_restore_sst       ! TK mod
  real(rk8)    :: sst_restore_timescale    ! TK mod; restoring time scale for SST in days
!
  real(rk8) , pointer , dimension(:,:) :: sstemp ! the actual prognotic sst pointing on tg2
  real(ik8) , pointer , dimension(:,:) :: ohfx, oqfx, ofsw, oflw
  real(rk8) , pointer , dimension(:,:) :: olndcat
  real(rk8), pointer , dimension (:,:) ::  qflux_restore_sst , qflux_adj,  net_hflx , hflx 
  integer(ik4), pointer , dimension (:,:) :: ocmask


 integer(ik8) :: dtocean
!
  public :: allocate_mod_slabocean , init_slabocean, update_slabocean
public :: islab_ocean,do_qflux_adj,do_restore_sst,mixed_layer_depth,  sst_restore_timescale
  contains
!
    subroutine allocate_mod_slabocean
      implicit none

        call getmem2d(qflux_restore_sst,jci1,jci2,ici1,ici2, &
                      'slab_ocean:qflux_restore_sst')

        call getmem2d(qflux_adj,jci1,jci2,ici1,ici2, &
                      'slab_ocean:qflux_adj')

        call getmem2d(net_hflx,jci1,jci2,ici1,ici2, &
                      'slab_ocean:net_hflx')

        call getmem2d(hflx,jci1,jci2,ici1,ici2, &
                      'slab_ocean:hflx')

    end subroutine allocate_mod_slabocean
!
    subroutine init_slabocean(sfs,ldmsk,fsw,flw)
      implicit none

! interface for regcm varaibele / slab ocean   
      type(surfstate) , intent(in) :: sfs
      integer(ik4) , pointer , intent(in) , dimension(:,:) :: ldmsk
      real(rk8) , pointer , intent(in) , dimension(:,:) :: fsw,flw
      
      call assignpnt(sfs%tgb, sstemp)
      call assignpnt(ldmsk , ocmask)
      call assignpnt(sfs%hfx ,ohfx)     
      call assignpnt(sfs%qfx , oqfx)     
      call assignpnt(fsw ,ofsw)     
      call assignpnt(flw , oflw)     



       dtocean =  dtsrf
       
 
    end subroutine init_slabocean


    subroutine update_slabocean
      implicit none
      
      real(rk8) :: mlcp

!      mlcp is the heat capacity of the mixed layer [J / m3 / deg C] * m
! water heat capacity ~ 4 J/g/K         
       mlcp = mixed_layer_depth*4e6

     qflux_restore_sst(:,:) = 0.0

      if (do_restore_sst) then

      where ( ocmask == 0 )
            qflux_restore_sst = (ts1(jci1:jci2, ici1:ici2) - sstemp(jci1:jci2, ici1:ici2)) &
                                 * mlcp / (sst_restore_timescale * 86400) ! w/m2                  
      end where


      if (do_qflux_adj) then
       
  !   Find the current climatological heat flux adjustment (qflux_adj).  
  !   The qflux_adj can be added to the heat flux to find the net
  !   heat flux anomaly into the mixed layer for the coupled mixed layer expt.
  
  !   Note that any ice model contributions need to be added in, which
  !   is handled in an off-line calculation to derive the final qflux
  !   adjustment.  These include ice restoring and ice lid contributions.

  !   Obtain indices for local array:
      
!ICI LIRE LE QFLUX ADJ FROM THE EXTERNAL FILE !!!!
!        call get_heat_flux_adj (Ocean_state%ocean_time, qflux_adj)       

      else
        qflux_adj = 0.0
      end if    

     endif

!
! energy budget in the mixed layer including the q flux therm
!
      print* , 'HE', maxval(ofsw ), maxval(oflw), maxval(  ohfx) , maxval(2.26D6 * oqfx), & 
                   maxval(qflux_restore_sst)
          

   where ( ocmask == 0 )

! The following are some key equations for this model:
! flux from or to the atmosphere ( convention  = positive downward)  
! multiply evaporation by latent heat of evaporation     
       hflx = ofsw - oflw - ohfx - wlhv * oqfx 

  

! account for the retaured or adjustment flux term
       net_hflx = hflx + qflux_adj + qflux_restore_sst

! update the prognostic seasurface temperature (pointing on tg2)

       sstemp(jci1:jci2, ici1:ici2)  = sstemp(jci1:jci2, ici1:ici2) + dtocean*net_hflx/mlcp
     

      end where






 
    end subroutine update_slabocean



   subroutine fill_slaboc_outvars
      implicit none


      if ( associated(slab_qflx_out) ) then

        slab_qflx_out  =  qflux_restore_sst(jci1:jci2,ici1:ici2)

      end if


   end subroutine fill_slaboc_outvars






!
end module mod_slabocean
