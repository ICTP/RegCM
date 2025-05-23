!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    Use of this source code is governed by an MIT-style license that can
!    be found in the LICENSE file or at
!
!         https://opensource.org/licenses/MIT.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

module mod_che_start

  use mod_intkinds
  use mod_runparams
  use mod_realkinds
  use mod_dynparam
  use mod_mpmessage
  use mod_che_common
  use mod_che_indices
  use mod_che_bdyco
  use mod_che_wetdep
  use mod_che_chemistry
  use mod_che_carbonaer
  use mod_che_ncio
  use mod_che_mppio
  use mod_che_bdyco
  use mod_che_dust
  use mod_che_sox
  use mod_che_seasalt
  use mod_mppparam
  use mod_che_molwg
  use mod_che_bionit
  use mod_cbmz_hvread

  implicit none

  private

  public  :: start_chem

  contains

  !--------------------------------------------------------------------------

  subroutine start_chem
    implicit none
    integer(ik4) :: i, j, k, n, itr, ibin, jbin, kbin, mmin, mbin
    integer(ik4) :: ipunit
    character(len=8) :: minamesav

    ! A : Intialise chemistry tracer indices

    dtchsolv = dtche

    iso2  = 0
    iso4  = 0
    isulf = 0
    idms  = 0
    imsa  = 0
    ibchl = 0
    nbchl = 0
    ibchb = 0
    nochl = 0
    iochl = 0
    iochb = 0
    idust = 0
    isslt = 0
    icarb = 0
    ianh4 = 0
    iano3 = 0
    ism1  = 0
    ism2  = 0

    ino     = 0
    ino2    = 0
    in2o5   = 0
    ihno2   = 0
    ihno3   = 0
    ihno4   = 0
    io3     = 0
    ih2o2   = 0
    ico     = 0
    iso2    = 0
    idms    = 0
    ih2so4  = 0
    ich4    = 0
    ic2h6   = 0
    ipar    = 0
    ich3oh  = 0
    ihcho   = 0
    iald2   = 0
    iaone   = 0
    ieth    = 0
    iolet   = 0
    iolei   = 0
    itol    = 0
    ixyl    = 0
    iisop   = 0
    ionit   = 0
    ipan    = 0
    ihcooh  = 0
    ircooh  = 0
    ich3ooh = 0
    iethooh = 0
    irooh   = 0
    imgly   = 0
    iisoprd = 0
    iisopn  = 0
    iopen   = 0
    icres   = 0
    ipollen = 0

    if ( nmine > 0 ) imine(:,:) = 0

    !abt *** For initializing megan tracer biogenic voc mask
    !    *** Mask not equal to zero when any MEGAN species is
    !    *** defined as a tracer within regcm.in
    !    *** If not equal to zero then that compound will be
    !    *** used as a surface emission from MEGAN and not
    !    *** from inventory (see surftnd.F for use)

#if defined CLM45 || (defined VOC && defined CLM)
    if ( igaschem == 1 ) then
      bvoc_trmask(:) = 0
    end if
#endif

    trac%indcbmz(:) = -1
    trac%mw(:) = d_zero

    ibin = 0
    jbin = 0
    kbin = 0
    mmin = 0
    mbin = 0
    do itr = 1, ntr
      if ( chtrname(itr) == 'SO2' ) then
        iso2 = itr
        chtrsol(iso2) = solso2
      end if
      if ( chtrname(itr) == 'DMS' ) idms = itr
      !
      ! For sulfuric acid and sulfate aer, we define here the same tracer
      ! index for compatibility with gas-phase chemistry options, and simple
      ! sulfate aer option. This might change with adding explicit sulfate
      ! aq chemistry.
      ! We consider here that all h2so4 partition in aerosol phase
      !
      if ( chtrname(itr) == 'SO4' .or. chtrname(itr) == 'H2SO4'    ) then
        ! sulfate index is added to carb vector for treatment in drydep
        ! and wetdep sulfate effective diameter and bin is taken equal to ochl
        kbin = kbin + 1
        ih2so4 = itr
        iso4 = itr
        icarb(kbin) = itr
        carbed(kbin) = reffochl(1)
        chtrsol(iso4) = solso4
      end if
      if ( chtrname(itr) == 'ANO3' ) then
        kbin = kbin + 1
        iano3 = itr
        icarb(kbin) = itr
        carbed(kbin) = reffochl(1)
        chtrsol(iano3) = solso4
      end if
      if ( chtrname(itr) == 'ANH4' ) then
        kbin = kbin + 1
        ianh4 = itr
        icarb(kbin) = itr
        carbed(kbin) = reffochl(1)
        chtrsol(ianh4) = solso4
      end if
      if ( chtrname(itr)(1:5) == 'BC_HL' ) then
        kbin = kbin + 1
        icarb(kbin) = itr
        nbchl = nbchl + 1
        ibchl(nbchl) = itr
        carbed(kbin) = reffochl(nbchl)
        chtrsol(itr) = solochl(nbchl)
      end if
      if ( chtrname(itr) == 'BC_HB' ) then
        kbin = kbin + 1
        ibchb = itr
        icarb(kbin) = itr
        carbed(kbin) = reffbc
        chtrsol(itr) = solbc
      end if
      if ( chtrname(itr)(1:5) == 'OC_HL' ) then
        kbin = kbin + 1
        icarb(kbin) = itr
        nochl = nochl + 1
        iochl(nochl) = itr
        carbed(kbin) = reffochl(nochl)
        chtrsol(itr) = solochl(nochl)
      end if
      if ( chtrname(itr) == 'OC_HB' ) then
        kbin = kbin + 1
        iochb = itr
        icarb(kbin) = itr
        carbed(kbin) = reffoc
        chtrsol(itr) = soloc
      end if
      if ( chtrname(itr) == 'SM1' ) then
        kbin = kbin + 1
        ism1 = itr
        icarb(kbin) = itr
        carbed(kbin) = reffsm1
        chtrsol(itr) = solsm1
      end if
      if ( chtrname(itr) == 'SM2' ) then
        kbin = kbin + 1
        ism2 = itr
        icarb(kbin) = itr
        carbed(kbin) = reffsm2
        chtrsol(itr) = solsm2
      end if
      if ( chtrname(itr)(1:4) ==  'DUST') then
        ibin = ibin + 1
        idust(ibin) = itr
      end if
      if ( chtrname(itr)(1:4) ==  'DU12') then
        ibin = ibin + 1
        idust(ibin) = itr
      end if
      if ( chtrname(itr)(1:4) ==  'SSLT') then
        jbin = jbin + 1
        isslt(jbin) = itr
        chtrsol(itr) = solsslt(jbin)
      end if
      ! mineralogical dust tracer
      if ( any(mine_name == chtrname(itr)(1:4)) ) then
        if ( chtrname(itr)(1:4) /= minamesav(1:4) ) then
          mbin = 0
          mmin = mmin+1
        end if
        mbin = mbin + 1
        imine(mbin,mmin) = itr
        minamesav = chtrname(itr)
      end if
!      if ( chtrname(itr)(1:3) == 'DMI' ) then
!        if ( chtrname(itr)(1:4) /= minamesav(1:4) ) then
!          mbin = 0
!          mmin = mmin+1
!        end if
!        mbin = mbin + 1
!        imine(mbin,mmin) = itr
!        minamesav = chtrname(itr)
!      end if

      ! gas phas species (CBMZ),
      ! max configuration : number of tracer = number of species

!      trac%indcbmz(:) = -1
!      trac%mw(:) = d_zero
      do n = 1,totsp
        if ( chtrname(itr) == cbmzspec(n) )then
          ! index of the tracer in the CBMZ list of species
          trac%indcbmz(itr) = n
          ! correponding molecular weight
          trac%mw(itr) = mw_cbmz(n)
        end if
      end do

      if ( myid == italk ) then
        if ( itr == 1 )  write(stdout,*) 'tracer', ' cbmz index',' molw'
        write(stdout,*) chtrname(itr),  trac%indcbmz(itr),  trac%mw(itr)
      end if

      !!$ Define also some specific indices for practical purpose
      !  Not all the possible cbmz species have a tracer index:
      !  however this information is also contained in trac%indcbmz table
      !
      !CBMZ mechanims
      if ( chtrname(itr) == 'NO'     ) ino         = itr
      if ( chtrname(itr) == 'NO2'    ) ino2        = itr
      if ( chtrname(itr) == 'N2O5'   ) in2o5       = itr
      if ( chtrname(itr) == 'HNO2'   ) ihno2       = itr
      if ( chtrname(itr) == 'HNO3'   ) ihno3       = itr
      if ( chtrname(itr) == 'HNO4'   ) ihno4       = itr
      if ( chtrname(itr) == 'O3'     ) io3         = itr
      if ( chtrname(itr) == 'H2O2'   ) ih2o2       = itr
      if ( chtrname(itr) == 'CO'     ) ico         = itr
      if ( chtrname(itr) == 'SO2'    ) iso2        = itr
      if ( chtrname(itr) == 'DMS'    ) idms        = itr
      if ( chtrname(itr) == 'H2SO4'  ) ih2so4      = itr
      if ( chtrname(itr) == 'CH4'    ) ich4        = itr
      if ( chtrname(itr) == 'C2H6'   ) ic2h6       = itr
      if ( chtrname(itr) == 'PAR'    ) ipar        = itr
      if ( chtrname(itr) == 'CH3OH'  ) ich3oh      = itr
      if ( chtrname(itr) == 'HCHO'   ) ihcho       = itr
      if ( chtrname(itr) == 'ALD2'   ) iald2       = itr
      if ( chtrname(itr) == 'AONE'   ) iaone       = itr
      if ( chtrname(itr) == 'ETH'    ) ieth        = itr
      if ( chtrname(itr) == 'OLET'   ) iolet       = itr
      if ( chtrname(itr) == 'OLEI'   ) iolei       = itr
      if ( chtrname(itr) == 'TOL'    ) itol        = itr
      if ( chtrname(itr) == 'XYL'    ) ixyl        = itr
      if ( chtrname(itr) == 'ISOP'   ) iisop       = itr
      if ( chtrname(itr) == 'ONIT'   ) ionit       = itr
      if ( chtrname(itr) == 'PAN'    ) ipan        = itr
      if ( chtrname(itr) == 'HCOOH'  ) ihcooh      = itr
      if ( chtrname(itr) == 'RCOOH'  ) ircooh      = itr
      if ( chtrname(itr) == 'CH3OOH' ) ich3ooh     = itr
      if ( chtrname(itr) == 'ETHOOH' ) iethooh     = itr
      if ( chtrname(itr) == 'ROOH'   ) irooh       = itr
      if ( chtrname(itr) == 'MGLY'   ) imgly       = itr
      if ( chtrname(itr) == 'ISOPRD' ) iisoprd     = itr
      if ( chtrname(itr) == 'ISOPN'  ) iisopn      = itr
      if ( chtrname(itr) == 'OPEN'   ) iopen       = itr
      if ( chtrname(itr) == 'CRES'   ) icres       = itr
      !!$
      if ( chtrname(itr) == 'NH3'   ) inh3       = itr

      if ( chtrname(itr) == 'POLLEN') ipollen   = itr

! special case of biogenic options
#if defined CLM45 || (defined VOC && defined CLM)
      !abt *** Added below to determine which MEGAN biogenic emission species
      !    *** will be passed to the gas phase mechanism
      !    *** commented out lines correspond to species not advected but
      !    *** potentially used in chemistry mechanism.
      !    *** Uncomment to give potential to advect
      if ( igaschem == 1 ) then
        if ( chtrname(itr) == 'ISOP'  ) bvoc_trmask(itr) = 1
        if ( chtrname(itr) == 'APIN'  ) bvoc_trmask(itr) = 1
        if ( chtrname(itr) == 'LIMO'  ) bvoc_trmask(itr) = 1
      end if
#endif

    end do

    cfdout = alarm_out_che%rw

    if ( count(icarb > 0) > 0 ) then
      call carb_init( )
    end if

    !
    ! define now correspndance between boundary species indices and
    ! determine tracer indices corresponding to ch boundary conditions
    !
    ichbdy2trac(:) = 0
    itr = 1
    if( igaschem == 1 ) then
      do n = 1, n_chbcvar
        do i = 1, ntr
          if (chbcname(n) == chtrname(i)) then
            ichbdy2trac(itr) = i
            itr = itr + 1
          end if
        end do
      end do
    end if
    !
    ! look also in aerosol bc and pile them after. (Expect when using DU12)
    !
    if ( iaerosol == 1 .and. chemsimtype(1:4) .ne. 'DU12' ) then
      do n = 1, size(aeaero)
        do i = 1, ntr
          if ( aeaero(n) == chtrname(i) ) then
            ichbdy2trac(itr) = i
            itr = itr + 1
          end if
        end do
      end do
    end if
    !
    ! Only when using DU12
    !
    if ( iaerosol == 1 .and. chemsimtype(1:4) == 'DU12' ) then
      do n = 1, size(aedu12)
        do i = 1, ntr
          if ( aedu12(n) == chtrname(i) ) then
            ichbdy2trac(itr) = i
            itr = itr + 1
          end if
        end do
      end do
    end if
    if ( myid == italk ) then
      write(stdout,*) 'tracer index coreesponding to bdy species '
      call iprntv(ichbdy2trac,size(ichbdy2trac),'ichbdy2trac')
    end if
!!$  FAB : work on that later
!!$    do itr = 1,ntr
!!$       do n = 1,n_chbcvar
!!$       if ( chtrname(itr) == chbcname(n))then
!!$        trac%indchbdy(itr) = n  ! index of the tracer in the chbdy list
!!$       end if
!!$      end do
!!$        print*,'test', itr, chtrname(itr), trac%indchbdy(itr)
!!$     end do
#ifndef CLM45
    if ( ichdustemd == 3 ) then
      call fatal(__FILE__,__LINE__, &
        'ichdustemd == 3 valid only if CLM45 is active.')
    end if
#endif
    if ( nbin > 0 .or. ichbion == 1 ) then
      ! activate dust initialization
      if ( myid == italk ) write(stdout,'(a)',advance='no') ' Calling inidust'
      call inidust
      if ( myid == italk ) write(stdout,*) '. Success.'
    end if

    if ( ichbion == 1 ) call ini_bionit

    if ( igaschem == 1 ) then
      open(newunit=ipunit,file='TUVGRID2', status='old', err=900)
      call hvread(ipunit)
    end if

    call init_mod_che_ncio(chemsimtype)
    call che_init_bdy

    ! Finally initialise tracer MXR to chib0 over the whole domain

    if ( .not. ifrest .or. ichecold == 1 ) then
      if ( myid == 0 ) then
        write(stdout,*) '#################################'
        write(stdout,*) '       Initializing tracers'
        write(stdout,*) '#################################'
      end if
      if ( idynamic == 3 ) then
        do n = 1, ntr
          do k = 1, kz
            do i = ice1, ice2
              do j = jce1, jce2
                chemt(j,i,k,n) = max(chib0(j,i,k,n),d_zero)
              end do
            end do
          end do
        end do
      else
        do n = 1, ntr
          do k = 1, kz
            do i = ice1, ice2
              do j = jce1, jce2
                chia(j,i,k,n) = max(chib0(j,i,k,n),d_zero)
                chib(j,i,k,n) = max(chib0(j,i,k,n),d_zero)
              end do
            end do
          end do
        end do
      end if
    end if

    return

900 write(stderr,*) 'Cannot open required file TUVGRID2.'
    call fatal(__FILE__,__LINE__,'TUVGRID2 NOT FOUND')

  end subroutine start_chem

end module mod_che_start
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
