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

module mod_che_start 

  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_mpmessage
  use mod_che_common
  use mod_che_indices
  use mod_che_bdyco
  use mod_che_wetdep
  use mod_che_carbonaer
  use mod_che_ncio
  use mod_che_mppio
  use mod_che_bdyco
  use mod_cbmz_init1
  use mod_che_dust
  use mod_che_sox
  use mod_che_seasalt
  use mod_mppparam  

  implicit none

  private

  real(rk8) , parameter :: solso4 = 0.9D0

  public  :: start_chem

  contains

  !--------------------------------------------------------------------------

  subroutine start_chem
    implicit none
    integer(ik4) :: i , j , k , itr , ibin , jbin , kbin
    integer(ik4) :: lyear , lmonth , lday , lhour

    ! A : Intialise chemistry tracer indices         

    iso2  = 0
    iso4  = 0
    idms  = 0
    imsa  = 0
    ibchl = 0
    ibchb = 0
    iochl = 0
    iochb = 0
    idust = 0
    isslt = 0
    icarb = 0

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
    ircooh =  0
    ico2   =  0
    in2o   =  0
    ipar   =  0
    iolt   =  0
    ioli   =  0
    imgly  =  0
    icres  =  0
    iopen  =  0
    iisoprd = 0
    ionit   = 0
    ihcooh  = 0
    ich3ooh = 0
    iethooh = 0
    irooh   = 0
    ixo2 =    0

    ipollen = 0


    !abt *** For initializing megan tracer biogenic voc mask  
    !    *** Mask not equal to zero when any MEGAN species is
    !    *** defined as a tracer within regcm.in  
    !    *** If not equal to zero then that compound will be
    !    *** used as a surface emission from MEGAN and not
    !    *** from inventory (see surftnd.F for use)

#if (defined VOC && defined CLM)
    if ( igaschem == 1 ) then
      bvoc_trmask(:) = 0
    end if
#endif

    ibin = 0
    jbin = 0
    kbin = 0

    do itr = 1 , ntr
      if ( chtrname(itr) == 'SO2' ) iso2 = itr
      if ( chtrname(itr) == 'DMS' ) idms = itr
      if ( chtrname(itr) == 'SO4' ) then
        ! sulfate index is added to carb vector for treatment in drydep
        ! and wetdep sulfate effective diameter and bin is taken equal to ochl
        kbin = kbin + 1
        iso4 = itr
        icarb(kbin) = itr
        carbed(kbin) = reffochl
        chtrsol(iso4) = solso4
      end if
      if ( chtrname(itr) == 'BC_HL' ) then 
        kbin = kbin + 1
        ibchl = itr
        icarb(kbin) = itr
        carbed(kbin) = reffbchl
        chtrsol(itr) = solbchl
      end if
      if ( chtrname(itr) == 'BC_HB' ) then 
        kbin = kbin + 1
        ibchb = itr
        icarb(kbin) = itr
        carbed(kbin) = reffbc
        chtrsol(itr) = solbc   
      end if
      if ( chtrname(itr) == 'OC_HL' ) then
        kbin = kbin + 1
        iochl = itr
        icarb(kbin) = itr
        carbed(kbin) = reffochl 
        chtrsol(itr) = soloc             
      end if
      if ( chtrname(itr) == 'OC_HB' ) then
        kbin = kbin + 1
        iochb = itr
        icarb(kbin) = itr 
        carbed(kbin) = reffoc 
        chtrsol(itr) = solochl  
      end if
      if ( chtrname(itr)(1:4) ==  'DUST') then
        ibin = ibin + 1
        idust(ibin) = itr
        chtrsol(itr) = soldust(ibin)   
      end if
      if ( chtrname(itr)(1:4) ==  'SSLT') then
        jbin = jbin + 1
        isslt(jbin) = itr
        chtrsol(itr) = solsslt(jbin)   
      end if

      ! gas phas species (CBMZ)

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
      if ( chtrname(itr) == 'RCOOH' ) ircooh    = itr
      if ( chtrname(itr) == 'CO2'   ) ico2      = itr
      if ( chtrname(itr) == 'DMS'   ) idms      = itr
      if ( chtrname(itr) == 'NOX'   ) inox      = itr
      if ( chtrname(itr) == 'HOX'   ) ihox      = itr
      if ( chtrname(itr) == 'SOX'   ) isox      = itr
      if ( chtrname(itr) == 'PAR'   ) ipar      = itr
      if ( chtrname(itr) == 'MGLY'  ) imgly     = itr
      if ( chtrname(itr) == 'CRES'  ) icres     = itr
      if ( chtrname(itr) == 'OPEN'  ) iopen     = itr
      if ( chtrname(itr) == 'ISOPRD') iisoprd   = itr 
      if ( chtrname(itr) == 'ONIT'  ) ionit     = itr
      if ( chtrname(itr) == 'HCOOH' ) ihcooh    = itr
      if ( chtrname(itr) == 'RCOOH' ) ircooh    = itr
      if ( chtrname(itr) == 'CH3OOH') ich3ooh   = itr 
      if ( chtrname(itr) == 'ETHOOH') iethooh   = itr 
      if ( chtrname(itr) == 'ROOH'  ) irooh     = itr
      if ( chtrname(itr) == 'HONO'  ) ihono     = itr
      if ( chtrname(itr) == 'HNO4'  ) ihno4     = itr
      if ( chtrname(itr) == 'XO2'   ) ixo2      = itr
      if ( chtrname(itr) == 'POLLEN') ipollen   = itr

      !abt *** Check to make sure SO4 is not defined twice as SULF or SO4 in
      !    *** regcm.in namelist.  If both are defined then STOP
      if ( iso4 /= 0 .and. isulf /= 0 ) then
        if ( myid == italk ) then
          write(*,*) "******* ERROR: Defined both SO4 and SULF"
          write(*,*) "*******        in regcm.in              "
          write(*,*) "*******        must choose one b/c they "
          write(*,*) "*******        both represent Sulfate   "
        end if
        call fatal(__FILE__,__LINE__,'CHEM CANNOT START')
      end if

#if (defined VOC && defined CLM)
      !abt *** Added below to determine which MEGAN biogenic emission species
      !    *** will be passed to the gas phase mechanism
      !    *** commented out lines correspond to species not advected but
      !    *** potentially used in chemistry mechanism.
      !    *** Uncomment to give potential to advect    
      if ( igaschem == 1 ) then
        if ( chtrname(itr) == 'ISOP'  ) bvoc_trmask(itr) = 1
      end if
#endif

    end do

    ! define now correspndance between boundary species indices and
    ! tracer indices
    ! must be absoutely consistent with ch_bdy  / depends on chem mechanism

    ichbdy2trac(:) = 0 

    select case (chemsimtype)
      case ('DUST')
        do i = 1 , ibin
          ichbdy2trac(i) = idust(i)
        end do
      case ('SSLT')
        do i = 1 , jbin
          ichbdy2trac(i) = isslt(i)
        end do
      case ('CARB')        
        ichbdy2trac(1) = ibchb
        ichbdy2trac(2) = ibchl
        ichbdy2trac(3) = iochb
        ichbdy2trac(4) = iochl
      case ('SULF')
        ichbdy2trac(1) = iso2
        ichbdy2trac(2) = iso4
      case ('SUCA')
        ichbdy2trac(1) = ibchb
        ichbdy2trac(2) = ibchl
        ichbdy2trac(3) = iochb
        ichbdy2trac(4) = iochl
        ichbdy2trac(5) = iso2
        ichbdy2trac(6) = iso4
      case ('AERO')
        ichbdy2trac(1) = ibchb
        ichbdy2trac(2) = ibchl
        ichbdy2trac(3) = iochb
        ichbdy2trac(4) = iochl
        ichbdy2trac(5) = iso2
        ichbdy2trac(6) = iso4
        itr = 6
        do i = 1 , jbin
          ichbdy2trac(itr+i) = isslt(i)
        end do
        itr = itr + jbin
        do i = 1 , ibin
          ichbdy2trac(itr+i) = idust(i)
        end do
      case ('CBMZ')
        ichbdy2trac(1)  = io3
        ichbdy2trac(2)  = ino
        ichbdy2trac(3)  = ino2
        ichbdy2trac(4)  = ihno3
        ichbdy2trac(5)  = in2o5
        ichbdy2trac(6)  = ih2o2
        ichbdy2trac(7)  = ich4
        ichbdy2trac(8)  = ico
        ichbdy2trac(9)  = ihcho
        ichbdy2trac(10) = imoh
        ichbdy2trac(11) = ieoh
        ichbdy2trac(12) = iethe
        ichbdy2trac(13) = ic2h6
        ichbdy2trac(14) = iald2
        ichbdy2trac(15) = iacet
        ichbdy2trac(16) = ioli
        ! ichbdy2trac(chbc_ivar(17)) = bigalk is no used here !!
        ichbdy2trac(17) = iolt
        ichbdy2trac(18) = ic3h8
        ichbdy2trac(19) = iisop
        ichbdy2trac(20) = itolue
        ichbdy2trac(21) = ipan
        ichbdy2trac(22) = iso2
        ichbdy2trac(23) = iso4
        ichbdy2trac(24) = idms
      case ('DCCB')
        ichbdy2trac(1)  = io3
        ichbdy2trac(2)  = ino
        ichbdy2trac(3)  = ino2
        ichbdy2trac(4)  = ihno3
        ichbdy2trac(5)  = in2o5
        ichbdy2trac(6)  = ih2o2
        ichbdy2trac(7)  = ich4
        ichbdy2trac(8)  = ico
        ichbdy2trac(9)  = ihcho
        ichbdy2trac(10) = imoh
        ichbdy2trac(11) = ieoh
        ichbdy2trac(12) = iethe
        ichbdy2trac(13) = ic2h6
        ichbdy2trac(14) = iald2
        ichbdy2trac(15) = iacet
        ichbdy2trac(16) = ioli
        ichbdy2trac(17) = iolt
        ichbdy2trac(18) = ic3h8
        ichbdy2trac(19) = iisop
        ichbdy2trac(20) = itolue
        ichbdy2trac(21) = ipan
        ichbdy2trac(22) = iso2
        ichbdy2trac(23) = iso4
        ichbdy2trac(24) = idms
        do i = 1 , ibin
          ichbdy2trac(i+24) = idust(i)
        end do
        ichbdy2trac(ibin+24+1) = ibchb
        ichbdy2trac(ibin+24+2) = ibchl
        ichbdy2trac(ibin+24+3) = iochb
        ichbdy2trac(ibin+24+4) = iochl
    end select

    if ( idust(1) > 0 ) then
      ! fisrt activate dust initialization
      if ( myid == italk ) write(stdout,*) 'Calling inidust'
      call inidust
    end if

    !*** abt added for wet deposition scheme
    ! if ( .not.allocated(chevap) ) allocate(chevap(iy,kz))
    ! if ( .not.allocated(checum) ) allocate(checum(iy,kz))

    !*** Initialize accumulation factor for output diagnostics 
    if ( ifrest ) then
      ! Care to the 0.5 factor added (leap frog related)
      cdiagf =  dt / (3600.0D0 * chemfrq) * d_half
    else
      cdiagf =  dt / (3600.0D0 * chemfrq)
    end if

    if ( igaschem == 1 ) then
      open(26,file='TUVGRID2', status='old', err=900)
      open(25,file='REACTION.DAT_CBMZ', status='old', err=901)
      open(27,file='cbmz_chemmech.out', status='replace', err=902)
902   continue
      call chemread
      call hvread
      call cheminit 
    end if

    call init_mod_che_ncio(chemsimtype)
    call che_init_bdy

    ! Finally initialise chia and chib to chib0 over the whole domain

    if ( .not. ifrest ) then
      do k = 1 , kz
        do i = ice1 , ice2
          do j = jce1 , jce2
            chia(j,i,k,:) = chib0(j,i,k,:)
            chib(j,i,k,:) = chib0(j,i,k,:)
          end do
        end do
      end do
    end if
!
    return

900 write(stderr,*) 'Cannot open required file TUVGRID2.'
    call fatal(__FILE__,__LINE__,'TUVGRID2 NOT FOUND')
901 write(stderr,*) 'Cannot open required file REACTION.DAT_CBMZ.'
    call fatal(__FILE__,__LINE__,'REACTION.DAT_CBMZ NOT FOUND')

  end subroutine start_chem

end module mod_che_start
