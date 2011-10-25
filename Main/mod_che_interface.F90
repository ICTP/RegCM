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
  use mod_atm_interface , only : slice , surfpstate
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
!
  public
!
  contains 
!
  subroutine init_chem(idirect,dt,chemfrq,dtrad,dsigma,atms,sps1, &
                       sps2,a,ptop,coszrs,icutop,icubot)
    implicit none
    integer , intent(in) :: idirect
    real(dp) , intent(in) :: dt , chemfrq , dtrad
    real(dp) , pointer , dimension(:) , intent(in) :: dsigma ! dsigma
    integer , pointer , dimension(:,:) :: icutop , icubot
    type(surfpstate) , intent(in) :: sps1 , sps2
    type(slice) , intent(in) :: atms
    real(dp) , pointer , dimension(:) :: a
    real(dp) , pointer , dimension(:,:) :: coszrs
    real(dp) :: ptop
    integer :: ibin , jbin , kbin , itr

    ichdir = idirect

    chfrq = chemfrq
    rafrq = dtrad
    dtche = dt
    chptop = ptop

    call assignpnt(dsigma,chlevs)
    call assignpnt(icutop,kcumtop)
    call assignpnt(icubot,kcumbot)
    call assignpnt(sps1%ps,psfcp)
    call assignpnt(sps2%ps,sfcp)
    call assignpnt(atms%tb3d,tatm)
    call assignpnt(atms%qvb3d,qvatm)
    call assignpnt(a,hlev)
    call assignpnt(coszrs,czen)
        
    iso2  = 0
    iso4  = 0
    idms  = 0
    imsa  = 0
    ibchl = 0
    ibchb = 0
    iochl = 0
    iochb = 0

    io3    =  0
    ino    =  0
    ino2   =  0
    ino3   =  0
    ioh    =  0
    iho2   =  0
    ih2o2  =  0
    ihno2  =  0
    ihno3  =  0
    ihno4  =  0
    isulf  =  0
    ih2so4 =  0
    ihono  =  0
    in2o5  =  0
    ihc    =  0
    ihcr   =  0
    ic2h4  =  0
    ico    =  0
    ihcho  =  0
    iald2  =  0
    iethe  =  0
    ic2h6  =  0
    ic3h8  =  0
    iisop  =  0
    itolue =  0
    ixyl   =  0
    inh3   =  0
    ipan   =  0
    irooh  =  0
    iacet  =  0
    ibenz  =  0
    inox   =  0
    ihox   =  0
    isox   =  0
    ich4   =  0
    ieoh   =  0
    imoh   =  0
    iaco2  =  0
    ico2   =  0
    in2o   =  0
    ipar   =  0
    iolt   =  0
    ioli   =  0

    ! abt
    ! For initializing megan tracer biogenic voc mask  
    ! Mask not equal to zero when any MEGAN species is defined as a
    ! tracer within regcm.in  
    ! If not equal to zero then that compound will be used as a surface
    ! emission from MEGAN and not from inventory (see surftnd.F for use)
    ! abt

#ifdef CLM
#if (defined VOC)
    if ( ibvoc == 1 ) then
      if (.not. allocated(bvoc_trmask)) then
        allocate(bvoc_trmask(ntr))
        bvoc_trmask(:) = 0
      endif
    endif
#endif
#endif

    ibin = 0
    jbin = 0
    kbin = 0

    ! initialise tracer indices and important parameter

    do itr = 1 , ntr
      if ( chtrname(itr) == 'SO2' ) iso2 = itr
      if ( chtrname(itr) == 'DMS')  idms = itr
      if ( chtrname(itr) == 'SO4' ) then
        ! sulfate index is added to carb vector for treatment in
        ! drydep and wetdep
        ! sulfate effective diameter and bin is taken equal to ochl
        kbin = kbin + 1
        iso4 = itr
        icarb(kbin) = itr
        carbsiz(kbin,1) = chreffochl - 0.1D0
        carbsiz(kbin,2) = chreffochl + 0.1D0
      end if
      if ( chtrname(itr) == 'BC_HL' ) then 
        kbin = kbin + 1
        ibchl = itr
        icarb(kbin) = itr
        carbsiz(kbin,1) = chreffbchl - 0.06D0
        carbsiz(kbin,2) = chreffbchl + 0.06D0
      end if
      if ( chtrname(itr) == 'BC_HB' ) then 
        kbin = kbin + 1
        ibchb = itr
        icarb(kbin) = itr
        carbsiz(kbin,1) = chreffbc - 0.01D0
        carbsiz(kbin,2) = chreffbc + 0.01D0
      end if
      if ( chtrname(itr) == 'OC_HL' ) then
        kbin = kbin + 1
        iochl = itr
        icarb(kbin) = itr
        carbsiz(kbin,1) = chreffochl - 0.1D0
        carbsiz(kbin,2) = chreffochl + 0.1D0
      end if
      if ( chtrname(itr) == 'OC_HB' ) then
        kbin = kbin + 1
        iochb = itr
        icarb(kbin) = itr 
        carbsiz(kbin,1) = chreffoc - 0.07D0
        carbsiz(kbin,2) = chreffoc + 0.07D0
      end if
      if ( chtrname(itr)(1:4) == 'DUST') then
        ibin = ibin + 1
        idust(ibin) = itr
      end if
      if ( chtrname(itr)(1:4) == 'SSLT') then
        jbin = jbin + 1
        isslt(jbin) = itr
      end if
      if ( chtrname(itr) == 'O3'    ) io3       = itr
      if ( chtrname(itr) == 'NO'    ) ino       = itr
      if ( chtrname(itr) == 'NO2'   ) ino2      = itr
      if ( chtrname(itr) == 'NO3'   ) ino3      = itr
      if ( chtrname(itr) == 'OH'    ) ioh       = itr
      if ( chtrname(itr) == 'HO2'   ) iho2      = itr
      if ( chtrname(itr) == 'H2O2'  ) ih2o2     = itr
      if ( chtrname(itr) == 'HNO2'  ) ihno2     = itr
      if ( chtrname(itr) == 'HNO3'  ) ihno3     = itr
      if ( chtrname(itr) == 'HNO4'  ) ihno4     = itr
      if ( chtrname(itr) == 'SULF'  ) isulf     = itr
      if ( chtrname(itr) == 'SO4'   ) iso4      = itr
      if ( chtrname(itr) == 'H2SO4' ) ih2so4    = itr
      if ( chtrname(itr) == 'HONO'  ) ihono     = itr
      if ( chtrname(itr) == 'N2O5'  ) in2o5     = itr
      if ( chtrname(itr) == 'HC'    ) ihc       = itr
      if ( chtrname(itr) == 'HCR'   ) ihcr      = itr
      if ( chtrname(itr) == 'C2H4'  ) ic2h4     = itr
      if ( chtrname(itr) == 'OLT'   ) iolt      = itr
      if ( chtrname(itr) == 'OLI'   ) ioli      = itr
      if ( chtrname(itr) == 'ALK4'  ) ialk4     = itr
      if ( chtrname(itr) == 'ALK7'  ) ialk7     = itr
      if ( chtrname(itr) == 'CO'    ) ico       = itr
      if ( chtrname(itr) == 'HCHO'  ) ihcho     = itr
      if ( chtrname(itr) == 'ALD2'  ) iald2     = itr
      if ( chtrname(itr) == 'ETHE'  ) iethe     = itr
      if ( chtrname(itr) == 'C2H6'  ) ic2h6     = itr
      if ( chtrname(itr) == 'C3H8'  ) ic3h8     = itr
      if ( chtrname(itr) == 'ISOP'  ) iisop     = itr
      if ( chtrname(itr) == 'TOLUE' ) itolue    = itr
      if ( chtrname(itr) == 'XYL'   ) ixyl      = itr
      if ( chtrname(itr) == 'NH3'   ) inh3      = itr
      if ( chtrname(itr) == 'PAN'   ) ipan      = itr
      if ( chtrname(itr) == 'ROOH'  ) irooh     = itr
      if ( chtrname(itr) == 'ACET'  ) iacet     = itr
      if ( chtrname(itr) == 'BENZ'  ) ibenz     = itr
      if ( chtrname(itr) == 'CH4'   ) ich4      = itr
      if ( chtrname(itr) == 'MOH'   ) imoh      = itr
      if ( chtrname(itr) == 'EOH'   ) ieoh      = itr
      if ( chtrname(itr) == 'ACO2'  ) iaco2     = itr
      if ( chtrname(itr) == 'CO2'   ) ico2      = itr
      if ( chtrname(itr) == 'DMS'   ) idms      = itr
      if ( chtrname(itr) == 'NOX'   ) inox      = itr
      if ( chtrname(itr) == 'HOX'   ) ihox      = itr
      if ( chtrname(itr) == 'SOX'   ) isox      = itr
      if ( chtrname(itr) == 'PAR'   ) ipar      = itr

      !abt *** Check to make sure SO4 is not defined twice as SULF or SO4 in
      !    *** regcm.in namelist.  If both are defined then STOP
      if ( iso4 /= 0 .and. isulf /= 0 ) then
        write(*,*)"******* ERROR: Defined both SO4 and SULF"
        write(*,*)"*******        in regcm.in              "
        write(*,*)"*******        must choose one b/c they "
        write(*,*)"*******        both represent Sulfate   "
        call fatal(__FILE__,__LINE__,'Defined both SO4 and SULF')
      end if

      ! abt
      ! Added below to determine which MEGAN biogenic emission species
      ! will be passed to the gas phase mechanism
      ! commented out lines correspond to species not advected but potentially
      ! used in chemistry mechanism.  Uncomment to give potential to advect
!!$
!!$#if (defined VOC)
!!$   if ( ibvoc == 1 ) then
!!$     if ( chtrname(itr) == 'ISOP'  ) bvoc_trmask(itr) = 1
!!$     if ( chtrname(itr) == 'APIN'  ) bvoc_trmask(itr) = 7
!!$     if ( chtrname(itr) == 'LIMO'  ) bvoc_trmask(itr) = 4
!!$   end if
!!$#endif
!!$
      ! abt above added 
    end do

  end subroutine init_chem
!
end module mod_che_interface
