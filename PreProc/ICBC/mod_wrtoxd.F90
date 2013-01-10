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

  public :: chv4 , oxv4 , aev4 , mw
  public :: nchsp , noxsp , naesp
  public :: chspec , oxspec , aespec

  public :: init_outoxd , close_outoxd
  public :: newfile_ch_icbc , newfile_ox_icbc , newfile_ae_icbc
  public :: write_ch_icbc , write_ox_icbc , write_ae_icbc

  character(len=256) :: ofname
  character(len=8) :: chtype

  integer(ik4) , parameter :: nchsp = 25
  integer(ik4) , parameter :: noxsp = 5
  integer(ik4) :: naesp = -1

  character(len=8) , dimension(nchsp) :: chspec
  character(len=8) , dimension(noxsp) :: oxspec
  real(rk8) , dimension(nchsp) :: mw

  character(len=8) , pointer , dimension(:) :: aespec
  character(len=8) , target , dimension(4) :: aedust
  character(len=8) , target , dimension(8) :: aedccb
  character(len=8) , target , dimension(4) :: aesslt
  character(len=8) , target , dimension(5) :: aecarb
  character(len=8) , target , dimension(2) :: aesulf
  character(len=8) , target , dimension(7) :: aesuca
  character(len=8) , target , dimension(15) :: aeaero

  real(rk8) , pointer , dimension(:,:,:,:) :: chv4
  real(rk8) , pointer , dimension(:,:,:,:) :: oxv4
  real(rk8) , pointer , dimension(:,:,:,:) :: aev4
  
  data oxspec / 'OH' , 'HO2' , 'O3' , 'NO3' , 'H2O2' /
  data chspec / 'O3' , 'NO' , 'NO2' , 'HNO3' , 'N2O5' , 'H2O2' , 'CH4' , &
                'CO' , 'CH2O' , 'CH3OH' , 'C2H5OH' , 'C2H4' , 'C2H6' ,   &
                'CH3CHO' , 'CH3COCH3' , 'BIGENE' , 'BIGALK' , 'C3H6' ,   &
                'C3H8' , 'ISOP' , 'TOLUENE' , 'PAN' , 'SO2' , 'SO4' , 'DMS' /
 
  data mw / 48.0D0,  30.0D0,  46.0D0,  63.0D0, 108.0D0,  34.0D0,  16.0D0, & 
            28.0D0,  30.0D0,  32.0D0,  46.0D0,  28.0D0,  30.0D0,  44.0D0, &
            58.0D0,  56.0D0,  72.0D0,  42.0D0,  44.0D0,  68.0D0,  92.0D0, &
           121.0D0,  64.0D0,  96.0D0,  62.0D0 /

  integer , parameter :: maxaeout = 16

  data aedust / 'DST01', 'DST02', 'DST03', 'DST04' /
  data aesslt / 'SSLT01' , 'SSLT02', 'SSLT03', 'SSLT04' /
  data aecarb / 'CB1' , 'CB2' , 'OC1' , 'SOA' , 'OC2' /
  data aesulf / 'SO2' , 'SO4' /
  data aesuca / 'CB1' , 'CB2' , 'OC1' , 'SOA' , 'OC2' , 'SO2' , 'SO4' /
  data aeaero / 'CB1' , 'CB2' , 'OC1' , 'SOA' , 'OC2' , 'SO2' , 'SO4' , &
                'SSLT01' , 'SSLT02', 'SSLT03', 'SSLT04' , 'DST01',      &
                'DST02', 'DST03', 'DST04' /
  data aedccb / 'CB1' , 'CB2' , 'OC1' , 'OC2' ,'DST01', 'DST02',  &
                'DST03', 'DST04' /

  integer(ik4) :: ioc2 , isoa
  integer(ik4) :: isslt1 , isslt2 , isslt3 , isslt4

  logical :: sum_soa_to_oc2
  logical :: sum_sslt_bins

  type(nc_output_stream) , save :: ncoutch
  type(nc_output_stream) , save :: ncoutox
  type(nc_output_stream) , save :: ncoutae

  type(ncvariable2d_real) , save , dimension(2) :: v2dvar_base
  type(ncvariable3d_real) , save , dimension(nchsp) :: v3dvar_ch
  type(ncvariable3d_real) , save , dimension(noxsp) :: v3dvar_ox
  type(ncvariable3d_real) , save , dimension(maxaeout) :: v3dvar_ae

  data sum_soa_to_oc2 /.false./
  data sum_sslt_bins  /.false./

  contains

  subroutine init_outoxd(chemsimtype)
    implicit none
    character(len=8) , intent(in) :: chemsimtype
    logical :: doaero , dochem , dooxcl
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
      case ( 'SSLT' )
        naesp = 4
        aespec => aesslt
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
      case ( 'SUCA' )
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
      case ( 'DCCB' )
        naesp = 8
        aespec => aedccb
        doaero = .true.
        dooxcl = .true.
        dochem = .true.
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
    if ( dochem ) call getmem4d(chv4,1,jx,1,iy,1,kz,1,nchsp,'mod_wrtoxd:chv4')
    if ( dooxcl ) call getmem4d(oxv4,1,jx,1,iy,1,kz,1,noxsp,'mod_wrtoxd:oxv4')
    if ( sum_soa_to_oc2 ) then
      ioc2 = -1
      isoa = -1
      do i = 1 , naesp
        if ( aespec(i) == 'SOA' ) isoa = i
        if ( aespec(i) == 'OC2' ) ioc2 = i
      end do
      if ( isoa < 0 .or. ioc2 < 0 ) then
        call fatal(__FILE__,__LINE__,'Logical error: Search SOA error')
      end if
    end if
    if ( sum_sslt_bins ) then
      isslt1 = -1
      isslt2 = -1
      isslt3 = -1
      isslt4 = -1
      do i = 1 , naesp
        if ( aespec(i) == 'SSLT01' ) isslt1 = i
        if ( aespec(i) == 'SSLT02' ) isslt2 = i
        if ( aespec(i) == 'SSLT03' ) isslt3 = i
        if ( aespec(i) == 'SSLT04' ) isslt4 = i
      end do
      if ( isslt1 < 0 .or. isslt2 < 0 .or. isslt3 < 0 .or. isslt4 < 0 ) then
        call fatal(__FILE__,__LINE__,'Logical error: Search SSLT error.')
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
    type(rcm_time_and_date) , intent(in) :: idate1
    type(ncoutstream_params) :: opar
    integer(ik4) :: ivar

    call outstream_dispose(ncoutch)
    write (ofname,'(a,a,a,a,i10,a)') trim(dirglob), pthsep, trim(domname), &
      '_CHBC.', toint10(idate1), '.nc'
    opar%fname = ofname
    opar%pname = 'chem_icbc'
    opar%zero_date = idate1
    opar%l_bound = .true.
    call outstream_setup(ncoutch,opar)
    call outstream_addatt(ncoutch, &
      ncattribute_string('simulation_type',chtype))
    call outstream_addvar(ncoutch,v2dvar_base(1))
    call outstream_addvar(ncoutch,v2dvar_base(2))
    do ivar = 1 , nchsp
      v3dvar_ch(ivar)%vname = chspec(ivar)
      v3dvar_ch(ivar)%vunit = 'kg kg-1'
      v3dvar_ch(ivar)%long_name = trim(chspec(ivar))//' Volume Mixing Ratio'
      v3dvar_ch(ivar)%standard_name = &
        'mass_fraction_of_'//trim(chspec(ivar))//'_in_air'
      v3dvar_ch(ivar)%lrecords = .true.
      v3dvar_ch(ivar)%is_slice = .true.
      v3dvar_ch(ivar)%rval_slice => chv4
      call outstream_addvar(ncoutch,v3dvar_ch(ivar))
    end do
    call outstream_enable(ncoutch,sigma2)
    call outstream_writevar(ncoutch,v2dvar_base(1))
    call outstream_writevar(ncoutch,v2dvar_base(2))
  end subroutine newfile_ch_icbc

  subroutine newfile_ox_icbc(idate1)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate1
    type(ncoutstream_params) :: opar
    integer(ik4) :: ivar

    call outstream_dispose(ncoutox)
    write (ofname,'(a,a,a,a,i10,a)') trim(dirglob), pthsep, trim(domname), &
      '_OXBC.', toint10(idate1), '.nc'
    opar%fname = ofname
    opar%pname = 'chem_icbc'
    opar%zero_date = idate1
    opar%l_bound = .true.
    call outstream_setup(ncoutox,opar)
    call outstream_addatt(ncoutox, &
      ncattribute_string('simulation_type',chtype))
    call outstream_addvar(ncoutox,v2dvar_base(1))
    call outstream_addvar(ncoutox,v2dvar_base(2))
    do ivar = 1 , noxsp
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
    call outstream_enable(ncoutox,sigma2)
    call outstream_writevar(ncoutox,v2dvar_base(1))
    call outstream_writevar(ncoutox,v2dvar_base(2))
  end subroutine newfile_ox_icbc

  subroutine newfile_ae_icbc(idate1)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate1
    type(ncoutstream_params) :: opar
    character(len=8) :: specname
    integer(ik4) :: ivar

    call outstream_dispose(ncoutae)
    write (ofname,'(a,a,a,a,i10,a)') trim(dirglob), pthsep, trim(domname), &
      '_AEBC.', toint10(idate1), '.nc'
    opar%fname = ofname
    opar%pname = 'chem_icbc'
    opar%zero_date = idate1
    opar%l_bound = .true.
    call outstream_setup(ncoutae,opar)
    call outstream_addatt(ncoutae, &
      ncattribute_string('simulation_type',chtype))
    call outstream_addvar(ncoutae,v2dvar_base(1))
    call outstream_addvar(ncoutae,v2dvar_base(2))
    do ivar = 1 , naesp
      if ( aespec(ivar) == 'SOA' ) cycle
      if ( aespec(ivar) == 'SSLT03' ) cycle
      if ( aespec(ivar) == 'SSLT04' ) cycle
      if ( aespec(ivar)(1:3) == 'DST' ) then
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
      v3dvar_ae(ivar)%long_name = trim(specname)//' Volume Mixing Ratio'
      v3dvar_ae(ivar)%standard_name = &
        'mass_fraction_of_'//trim(specname)//'_in_air'
      v3dvar_ae(ivar)%lrecords = .true.
      v3dvar_ae(ivar)%is_slice = .true.
      v3dvar_ae(ivar)%rval_slice => aev4
      call outstream_addvar(ncoutae,v3dvar_ae(ivar))
    end do
    call outstream_enable(ncoutae,sigma2)
    call outstream_writevar(ncoutae,v2dvar_base(1))
    call outstream_writevar(ncoutae,v2dvar_base(2))
  end subroutine newfile_ae_icbc

  subroutine write_ch_icbc(idate)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate
    integer(ik4):: ivar
    call outstream_addrec(ncoutch,idate)
    do ivar = 1 , nchsp
      call outstream_writevar(ncoutch,v3dvar_ch(ivar),is=ivar)
    end do
    write (stdout ,*) 'Write ch_icbc : ', tochar(idate)
  end subroutine write_ch_icbc

  subroutine write_ox_icbc(idate)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate
    integer(ik4) :: ivar
    call outstream_addrec(ncoutox,idate)
    do ivar = 1 , noxsp
      call outstream_writevar(ncoutox,v3dvar_ox(ivar),is=ivar)
    end do
    write (stdout ,*) 'Write ox_icbc : ', tochar(idate)
  end subroutine write_ox_icbc
!
  subroutine write_ae_icbc(idate)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate
    integer(ik4) :: ivar

    if ( sum_sslt_bins ) then
      aev4(:,:,:,isslt1) = aev4(:,:,:,isslt1) + aev4(:,:,:,isslt2)
      aev4(:,:,:,isslt2) = aev4(:,:,:,isslt3) + aev4(:,:,:,isslt4)
    end if

    if ( sum_soa_to_oc2 ) then
      aev4(:,:,:,ioc2) = aev4(:,:,:,ioc2) + aev4(:,:,:,isoa)
    end if

    call outstream_addrec(ncoutae,idate)
    do ivar = 1 , naesp
      if ( aespec(ivar) == 'SOA' ) cycle
      if ( aespec(ivar) == 'SSLT03' ) cycle
      if ( aespec(ivar) == 'SSLT04' ) cycle
      call outstream_writevar(ncoutae,v3dvar_ae(ivar),is=ivar)
    end do
    write (stdout ,*) 'Write ae_icbc : ', tochar(idate)
  end subroutine write_ae_icbc

end module mod_wrtoxd
