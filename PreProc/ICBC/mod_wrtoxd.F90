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

module mod_wrtoxd

  use mod_dynparam
  use mod_date
  use mod_grid
  use mod_realkinds
  use mod_memutil
  use mod_message
  use mod_stdio
  use mod_ncstream_types
  use mod_ncstream

  private

  public :: chv4, oxv4, aev4, mw
  public :: nchsp, noxsp, naesp, ncbmz
  public :: chspec, oxspec, aespec, cbmzspec

  public :: init_outoxd, close_outoxd
  public :: newfile_ch_icbc, newfile_ox_icbc, newfile_ae_icbc, newfile_ae_icbc1
  public :: write_ch_icbc, write_ox_icbc, write_ae_icbc, write_ae_icbc1

  character(len=256) :: ofname
  character(len=8) :: chtype

  ! species in Mozart output
  integer(ik4), parameter :: nchsp = 40
  ! cbmz species in chemistry lateral boundaries
  integer(ik4), parameter :: ncbmz = 33
  integer(ik4), parameter :: noxsp = 5
  !aero species in chemistry lateral boundaries
  integer(ik4), parameter :: naero = 12

  integer(ik4) :: naesp = -1

  character(len=8), dimension(nchsp) :: chspec   ! Names of Mozart species
  character(len=8), dimension(ncbmz) :: cbmzspec ! Name of CBMZ species
  character(len=8), dimension(noxsp) :: oxspec
  character(len=8), dimension(naero) :: aerospec

  real(rkx), dimension(nchsp) :: mw

  character(len=8), pointer, contiguous, dimension(:) :: aespec
  character(len=8), target, dimension(4) :: aedust
  character(len=8), target, dimension(12) :: aedu12
  character(len=8), target, dimension(8) :: aedccb
  character(len=8), target, dimension(4) :: aesslt
  character(len=8), target, dimension(8) :: aeduss
  character(len=8), target, dimension(5) :: aecarb
  character(len=8), target, dimension(2) :: aesulf
  character(len=8), target, dimension(7) :: aesuca
  character(len=8), target, dimension(15) :: aeaero

  real(rkx), pointer, contiguous, dimension(:,:,:,:) :: chv4
  real(rkx), pointer, contiguous, dimension(:,:,:,:) :: oxv4
  real(rkx), pointer, contiguous, dimension(:,:,:,:) :: aev4

  data oxspec / 'OH', 'HO2', 'O3', 'NO3', 'H2O2' /

  data chspec / 'NO      ','NO2     ','N2O5    ','HNO3    ','HO2NO2  ',  &
                'O3      ','H2O2    ','SO2     ','SO4     ','CH4     ',  &
                'CH2O    ','CH3OH   ','PAN     ','C2H6    ','C3H8    ',  &
                'BIGALK  ','C2H4    ','C3H6    ','BIGENE  ','TOLUENE ',  &
                'ISOP    ','CH3CHO  ','CH3COOH ','GLYALD  ','CH3OOH  ',  &
                'C2H5OOH ','CH3COCH3','HYAC    ','CH3COCHO','ONIT    ',  &
                'MEK     ','MVK     ','MACR    ','HYDRALD ','BIGALD  ',  &
                'ISOPNO3 ','ONITR   ','CRESOL  ','CO      ','DMS     ' /


 data cbmzspec / 'O3      ','NO      ','NO2     ','HNO3    ','HNO4    ', &
                 'N2O5    ','H2O2    ','CH4     ','CO      ','SO2     ', &
                 'H2SO4   ','DMS     ','PAR     ','C2H6    ','ETH     ', &
                 'OLET    ','OLEI    ','TOL     ','XYL     ','ISOP    ', &
                 'CRES    ','OPEN    ','ISOPN   ','ISOPRD  ','ONIT    ', &
                 'MGLY    ','AONE    ','PAN     ','CH3OOH  ','ETHOOH  ', &
                 'ALD2    ','HCHO    ','CH3OH   '/
! This is usedi only for CAMS preproc, in the future try to harmonize..
 data aerospec /'BC_HL ','BC_HB ','OC_HL ','OC_HB ', &
                'SO2   ','SO4   ','DUST01','DUST02', &
                'DUST03','DUST04','SSLT01','SSLT02'/


  integer, parameter :: maxaeout = 16

  data aedust / 'DST01', 'DST02', 'DST03', 'DST04' /
  data aedu12 / 'D1201', 'D1202', 'D1203', 'D1204', &
                'D1205', 'D1206', 'D1207', 'D1208', &
                'D1209', 'D1210', 'D1211', 'D1212' /
  data aesslt / 'SSLT01', 'SSLT02', 'SSLT03', 'SSLT04' /
  data aeduss / 'DST01', 'DST02', 'DST03', 'DST04', &
                'SSLT01', 'SSLT02', 'SSLT03', 'SSLT04' /
  data aecarb / 'CB1', 'CB2', 'OC1', 'SOA', 'OC2' /
  data aesulf / 'SO2', 'SO4' /
  data aesuca / 'CB1', 'CB2', 'OC1', 'SOA', 'OC2', 'SO2', 'SO4' /
  data aeaero / 'CB1', 'CB2', 'OC1', 'SOA', 'OC2', 'SO2', 'SO4', &
                'SSLT01', 'SSLT02', 'SSLT03', 'SSLT04', 'DST01',      &
                'DST02', 'DST03', 'DST04' /
  data aedccb / 'CB1', 'CB2', 'OC1', 'OC2' ,'DST01', 'DST02',  &
                'DST03', 'DST04' /

  integer(ik4) :: ioc2, isoa
  integer(ik4) :: isslt1, isslt2, isslt3, isslt4

  logical :: sum_soa_to_oc2
  logical :: sum_sslt_bins

  type(nc_output_stream), save :: ncoutch
  type(nc_output_stream), save :: ncoutox
  type(nc_output_stream), save :: ncoutae

  type(ncvariable2d_mixed), save, dimension(2) :: v2dvar_base
  type(ncvariable3d_mixed), save, dimension(ncbmz) :: v3dvar_ch
  type(ncvariable3d_mixed), save, dimension(noxsp) :: v3dvar_ox
  type(ncvariable3d_mixed), save, dimension(maxaeout) :: v3dvar_ae

  data sum_soa_to_oc2 /.false./
  data sum_sslt_bins  /.false./

  contains

  subroutine init_outoxd(chemsimtype)
    implicit none
    character(len=8), intent(in) :: chemsimtype
    logical :: doaero, dochem, dooxcl
    integer(ik4) :: i
    data doaero /.false./
    data dochem /.false./
    data dooxcl /.false./
    chtype = chemsimtype
    select case ( chemsimtype )
      case ( 'DUST' )
        naesp = 4
        aespec => aedust
        doaero = .true.
      case ( 'DU12' )
        naesp = 12
        aespec => aedu12
        doaero = .true.
      case ( 'SSLT' )
        naesp = 4
        aespec => aesslt
        doaero = .true.
        sum_sslt_bins = .true.
      case ( 'DUSS' )
        naesp = 8
        aespec => aeduss
        doaero = .true.
        sum_sslt_bins = .true.
      case ( 'CARB' )
        naesp = 5
        aespec => aecarb
        doaero = .true.
        dooxcl = .true.
        sum_soa_to_oc2 = .true.
      case ( 'SULF' )
        naesp = 2
        aespec => aesulf
        doaero = .true.
        dooxcl = .true.
      case ( 'SUCA', 'SUCE' )
        naesp = 7
        aespec => aesuca
        doaero = .true.
        dooxcl = .true.
        sum_soa_to_oc2 = .true.
      case ( 'AERO' )
        naesp = 15
        aespec => aeaero
        doaero = .true.
        dooxcl = .true.
        sum_sslt_bins = .true.
        sum_soa_to_oc2 = .true.
      case ( 'CBMZ' )
        dochem = .true.
      case ('DCCB')
        naesp = 15
        aespec => aeaero
        doaero = .true.
        dooxcl = .true.
        dochem=.true.
        sum_sslt_bins = .true.
      case default
    call die('init_outoxd','Unknown chemsimtype')
    end select
    if ( doaero ) then
      if ( naesp > maxaeout ) then
        write(stderr,*) 'Added more species without increasing maxaeout'
        call die('init_outoxd','INCRESE maxaeout AND RECOMPILE')
      end if
      call getmem4d(aev4,1,jx,1,iy,1,kz,1,naesp,'mod_wrtoxd:aev4')
    end if
    if ( dochem ) call getmem4d(chv4,1,jx,1,iy,1,kz,1,ncbmz,'mod_wrtoxd:chv4')
    if ( dooxcl ) call getmem4d(oxv4,1,jx,1,iy,1,kz,1,noxsp,'mod_wrtoxd:oxv4')
    if ( sum_soa_to_oc2 ) then
      ioc2 = -1
      isoa = -1
      do i = 1, naesp
        if ( aespec(i) == 'SOA' ) isoa = i
        if ( aespec(i) == 'OC2' ) ioc2 = i
      end do
      if ( isoa < 0 .or. ioc2 < 0 ) then
        call fatal(__FILE__,__LINE__, &
                   'Logical error: Search SOA error')
      end if
    end if
    if ( sum_sslt_bins ) then
      isslt1 = -1
      isslt2 = -1
      isslt3 = -1
      isslt4 = -1
      do i = 1, naesp
        if ( aespec(i) == 'SSLT01' ) isslt1 = i
        if ( aespec(i) == 'SSLT02' ) isslt2 = i
        if ( aespec(i) == 'SSLT03' ) isslt3 = i
        if ( aespec(i) == 'SSLT04' ) isslt4 = i
      end do
      if ( isslt1 < 0 .or. isslt2 < 0 .or. isslt3 < 0 .or. isslt4 < 0 ) then
        call fatal(__FILE__,__LINE__, &
                   'Logical error: Search SSLT error.')
      end if
    end if
    v2dvar_base(1)%vname = 'xlon'
    v2dvar_base(1)%vunit = 'degrees_east'
    v2dvar_base(1)%long_name = 'Longitude on Cross Points'
    v2dvar_base(1)%standard_name = 'longitude'
    v2dvar_base(2)%vname = 'xlat'
    v2dvar_base(2)%vunit = 'degrees_north'
    v2dvar_base(2)%long_name = 'Latitude on Cross Points'
    v2dvar_base(2)%standard_name = 'latitude'
    v2dvar_base(1)%rval => xlon
    v2dvar_base(2)%rval => xlat
  end subroutine init_outoxd

  subroutine close_outoxd
    implicit none
    call outstream_dispose(ncoutch)
    call outstream_dispose(ncoutox)
    call outstream_dispose(ncoutae)
  end subroutine close_outoxd

  subroutine newfile_ch_icbc(idate1)
    implicit none
    type(rcm_time_and_date), intent(in) :: idate1
    type(ncoutstream_params) :: opar
    integer(ik4) :: ivar

    call outstream_dispose(ncoutch)
    write (ofname,'(a,a,a,a,a,a)') trim(dirglob), pthsep, trim(domname), &
      '_CHBC.', trim(tochar10(idate1)), '.nc'
    opar%fname = ofname
    opar%pname = 'chem_icbc'
    opar%zero_date = idate1
    opar%l_bound = .true.
    call outstream_setup(ncoutch,opar)
    call outstream_addatt(ncoutch, &
      ncattribute_string('simulation_type',chtype))
    call outstream_addvar(ncoutch,v2dvar_base(1))
    call outstream_addvar(ncoutch,v2dvar_base(2))
    do ivar = 1, ncbmz
      v3dvar_ch(ivar)%vname = cbmzspec(ivar)
      v3dvar_ch(ivar)%vunit = 'kg kg-1'
      v3dvar_ch(ivar)%long_name = trim(cbmzspec(ivar))//' Volume Mixing Ratio'
      v3dvar_ch(ivar)%standard_name = &
        'mass_fraction_of_'//trim(cbmzspec(ivar))//'_in_air'
      v3dvar_ch(ivar)%lrecords = .true.
      v3dvar_ch(ivar)%is_slice = .true.
      v3dvar_ch(ivar)%rval_slice => chv4
      call outstream_addvar(ncoutch,v3dvar_ch(ivar))
    end do
    call outstream_enable(ncoutch,sigmah)
    call outstream_writevar(ncoutch,v2dvar_base(1))
    call outstream_writevar(ncoutch,v2dvar_base(2))
  end subroutine newfile_ch_icbc

  subroutine newfile_ox_icbc(idate1)
    implicit none
    type(rcm_time_and_date), intent(in) :: idate1
    type(ncoutstream_params) :: opar
    integer(ik4) :: ivar

    call outstream_dispose(ncoutox)
    write (ofname,'(a,a,a,a,a,a)') trim(dirglob), pthsep, trim(domname), &
      '_OXBC.', trim(tochar10(idate1)), '.nc'
    opar%fname = ofname
    opar%pname = 'chem_icbc'
    opar%zero_date = idate1
    opar%l_bound = .true.
    call outstream_setup(ncoutox,opar)
    call outstream_addatt(ncoutox, &
      ncattribute_string('simulation_type',chtype))
    call outstream_addvar(ncoutox,v2dvar_base(1))
    call outstream_addvar(ncoutox,v2dvar_base(2))
    do ivar = 1, noxsp
      v3dvar_ox(ivar)%vname = oxspec(ivar)
      v3dvar_ox(ivar)%vunit = 'kg kg-1'
      v3dvar_ox(ivar)%long_name = trim(oxspec(ivar))//' Volume Mixing Ratio'
      v3dvar_ox(ivar)%standard_name = &
        'mass_fraction_of_'//trim(oxspec(ivar))//'_in_air'
      v3dvar_ox(ivar)%lrecords = .true.
      v3dvar_ox(ivar)%is_slice = .true.
      v3dvar_ox(ivar)%rval_slice => oxv4
      call outstream_addvar(ncoutox,v3dvar_ox(ivar))
    end do
    call outstream_enable(ncoutox,sigmah)
    call outstream_writevar(ncoutox,v2dvar_base(1))
    call outstream_writevar(ncoutox,v2dvar_base(2))
  end subroutine newfile_ox_icbc


   subroutine newfile_ae_icbc(idate1)
    implicit none
    type(rcm_time_and_date), intent(in) :: idate1
    type(ncoutstream_params) :: opar
    integer(ik4) :: ivar

    call outstream_dispose(ncoutae)
    write (ofname,'(a,a,a,a,a,a)') trim(dirglob), pthsep, trim(domname), &
      '_AEBC.', trim(tochar10(idate1)), '.nc'
    opar%fname = ofname
    opar%pname = 'chem_icbc'
    opar%zero_date = idate1
    opar%l_bound = .true.
    call outstream_setup(ncoutae,opar)
    call outstream_addatt(ncoutae, &
      ncattribute_string('simulation_type',chtype))
    call outstream_addvar(ncoutae,v2dvar_base(1))
    call outstream_addvar(ncoutae,v2dvar_base(2))
    do ivar = 1, naero
      v3dvar_ae(ivar)%vname = aerospec(ivar)
      v3dvar_ae(ivar)%vunit = 'kg kg-1'
      v3dvar_ae(ivar)%long_name = trim(aerospec(ivar))//' Mass Mixing Ratio'
      v3dvar_ae(ivar)%standard_name = &
        'mass_fraction_of_'//trim(aerospec(ivar))//'_in_air'
      v3dvar_ae(ivar)%lrecords = .true.
      v3dvar_ae(ivar)%is_slice = .true.
      v3dvar_ae(ivar)%rval_slice => aev4
      call outstream_addvar(ncoutae,v3dvar_ae(ivar))
    end do
    call outstream_enable(ncoutae,sigmah)
    call outstream_writevar(ncoutae,v2dvar_base(1))
    call outstream_writevar(ncoutae,v2dvar_base(2))
  end subroutine newfile_ae_icbc


  subroutine newfile_ae_icbc1(idate1)
    implicit none
    type(rcm_time_and_date), intent(in) :: idate1
    type(ncoutstream_params) :: opar
    character(len=8) :: specname
    integer(ik4) :: ivar

    call outstream_dispose(ncoutae)
    write (ofname,'(a,a,a,a,a,a)') trim(dirglob), pthsep, trim(domname), &
      '_AEBC.', trim(tochar10(idate1)), '.nc'
    opar%fname = ofname
    opar%pname = 'chem_icbc'
    opar%zero_date = idate1
    opar%l_bound = .true.
    call outstream_setup(ncoutae,opar)
    call outstream_addatt(ncoutae, &
      ncattribute_string('simulation_type',chtype))
    call outstream_addvar(ncoutae,v2dvar_base(1))
    call outstream_addvar(ncoutae,v2dvar_base(2))
    do ivar = 1, naesp
      if ( aespec(ivar) == 'SOA' ) cycle
      if ( aespec(ivar) == 'SSLT03' ) cycle
      if ( aespec(ivar) == 'SSLT04' ) cycle
      if ( aespec(ivar)(1:3) == 'DST' ) then
        specname = 'DUST'//aespec(ivar)(4:5)
      else if ( aespec(ivar)(1:3) == 'D12' ) then
        specname = 'DUST'//aespec(ivar)(4:5)
      else if ( aespec(ivar)(1:3) == 'OC1' ) then
        specname = 'OC_HB'
      else if ( aespec(ivar)(1:3) == 'OC2' ) then
        specname = 'OC_HL'
      else if ( aespec(ivar)(1:3) == 'CB1' ) then
        specname = 'BC_HB'
      else if ( aespec(ivar)(1:3) == 'CB2' ) then
        specname = 'BC_HL'
      else
        specname = aespec(ivar)
      end if
      v3dvar_ae(ivar)%vname = specname
      v3dvar_ae(ivar)%vunit = 'kg kg-1'
      v3dvar_ae(ivar)%long_name = trim(specname)//' Mass Mixing Ratio'
      v3dvar_ae(ivar)%standard_name = &
        'mass_fraction_of_'//trim(specname)//'_in_air'
      v3dvar_ae(ivar)%lrecords = .true.
      v3dvar_ae(ivar)%is_slice = .true.
      v3dvar_ae(ivar)%rval_slice => aev4
      call outstream_addvar(ncoutae,v3dvar_ae(ivar))
    end do
    call outstream_enable(ncoutae,sigmah)
    call outstream_writevar(ncoutae,v2dvar_base(1))
    call outstream_writevar(ncoutae,v2dvar_base(2))
  end subroutine newfile_ae_icbc1

  subroutine write_ch_icbc(idate)
    implicit none
    type(rcm_time_and_date), intent(in) :: idate
    integer(ik4):: ivar
    call outstream_addrec(ncoutch,idate)
    do ivar = 1, ncbmz
      call outstream_writevar(ncoutch,v3dvar_ch(ivar),is=ivar)
    end do
    write (stdout ,*) 'Write ch_icbc : ', tochar(idate)
  end subroutine write_ch_icbc

  subroutine write_ae_icbc(idate)
    implicit none
    type(rcm_time_and_date), intent(in) :: idate
    integer(ik4):: ivar
    call outstream_addrec(ncoutae,idate)
    do ivar = 1, naero
      call outstream_writevar(ncoutae,v3dvar_ae(ivar),is=ivar)
    end do
    write (stdout ,*) 'Write ae_icbc : ', tochar(idate)
  end subroutine write_ae_icbc

  subroutine write_ox_icbc(idate)
    implicit none
    type(rcm_time_and_date), intent(in) :: idate
    integer(ik4) :: ivar
    call outstream_addrec(ncoutox,idate)
    do ivar = 1, noxsp
      call outstream_writevar(ncoutox,v3dvar_ox(ivar),is=ivar)
    end do
    write (stdout ,*) 'Write ox_icbc : ', tochar(idate)
  end subroutine write_ox_icbc
!
  subroutine write_ae_icbc1(idate)
    implicit none
    type(rcm_time_and_date), intent(in) :: idate
    integer(ik4) :: ivar

    if ( sum_sslt_bins ) then
      aev4(:,:,:,isslt1) = aev4(:,:,:,isslt1) + aev4(:,:,:,isslt2)
      aev4(:,:,:,isslt2) = aev4(:,:,:,isslt3) + aev4(:,:,:,isslt4)
    end if

    if ( sum_soa_to_oc2 ) then
      aev4(:,:,:,ioc2) = aev4(:,:,:,ioc2) + aev4(:,:,:,isoa)
    end if

    call outstream_addrec(ncoutae,idate)
    do ivar = 1, naesp
      if ( aespec(ivar) == 'SOA' ) cycle
      if ( aespec(ivar) == 'SSLT03' ) cycle
      if ( aespec(ivar) == 'SSLT04' ) cycle
      call outstream_writevar(ncoutae,v3dvar_ae(ivar),is=ivar)
    end do
    write (stdout ,*) 'Write ae_icbc : ', tochar(idate)
  end subroutine write_ae_icbc1

end module mod_wrtoxd
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
