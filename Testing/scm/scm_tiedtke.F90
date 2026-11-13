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

! Single column test of the Tiedtke cumulus scheme (icup = 5).
!
! Usage: scm_tiedtke [namelist file]       (default: scm_tiedtke.in)
!
! The columns of the profile file are laid out as a 1 x ncol grid and the
! scheme is called once through its model entry point, tiedtkedrv, as the
! MOLOCH core does (idynamic = 3, NoTo microphysics, no chemistry). The
! namelist has a &scmparam group (files, time step, grid spacing) and the
! model groups &tiedtkeparam and &cldparam (only kfac_deep, kfac_shal,
! k2_const), with the RegCM defaults. The output file is plain text, read
! by scm_tiedtke.py: a header with the namelists used, then per column a
! COLUMN line and one line per level, top to bottom.

program scm_tiedtke
  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_dynparam
  use mod_runparams
  use mod_memutil
  use mod_regcm_types
  use mod_cu_common
  use mod_cu_tiedtke, only : allocate_mod_cu_tiedtke, tiedtkedrv, &
                             pmean, ktype, pmfu
  implicit none

  character(len=256) :: nmlfile, infile, outfile
  real(rkx) :: dxkm
  integer(ik4) :: ncol, nk, j, k, iu, ios
  character(len=32), allocatable, dimension(:) :: col_name
  character(len=10), allocatable, dimension(:) :: col_date
  real(rkx), allocatable, dimension(:) :: col_lat, col_lon, col_ps, &
    col_shf, col_lhf, col_zpbl, col_cape, col_crr, col_topo
  integer(ik4), allocatable, dimension(:) :: col_land
  type(mod_2_cum) :: m2c
  real(rkx), pointer, contiguous, dimension(:,:,:) :: omega
  real(rkx) :: rain, heat, moist, dp
  character(len=512) :: line

  namelist /scmparam/ infile, outfile, dtcum, dxkm
  namelist /tiedtkeparam/ entrdd, entrpen_lnd, entrpen_ocn, &
    detrpen_lnd, detrpen_ocn, entshalp, rcuc_lnd, rcuc_ocn, &
    rcpec_lnd, rcpec_ocn, rhebc_lnd, rhebc_ocn, rprc_ocn, rprc_lnd, &
    revap_lnd, revap_ocn, lmfpen, lmfmid, lmfdd, lepcld, lmfdudv,   &
    lmfscv, lmfuvdis, lmftrac, lmfsmooth, lmfwstar, lmfcapebl
  namelist /cldparam/ kfac_deep, kfac_shal, k2_const

  ! Defaults: as in Main/mod_params.F90
  infile = 'scm_profiles.txt'
  outfile = 'scm_tiedtke.out'
  dtcum = 300.0_rkx
  dxkm = 25.0_rkx
  entrpen_lnd = 1.75e-3_rkx
  entrpen_ocn = 1.75e-3_rkx
  detrpen_lnd = 0.75e-4_rkx
  detrpen_ocn = 0.75e-4_rkx
  entrdd = 3.0e-4_rkx
  entshalp = 2.0_rkx
  rcuc_lnd = 0.05_rkx
  rcuc_ocn = 0.05_rkx
  rcpec_lnd = 5.55e-5_rkx
  rcpec_ocn = 5.55e-5_rkx
  rhebc_lnd = 0.8_rkx
  rhebc_ocn = 0.8_rkx
  rprc_lnd = 1.4e-3_rkx
  rprc_ocn = 1.4e-3_rkx
  revap_lnd = 1.0e-5_rkx
  revap_ocn = 1.0e-5_rkx
  lmfpen = .true.
  lmfmid = .true.
  lmfdd = .true.
  lepcld = .true.
  lmfdudv = .true.
  lmfscv = .true.
  lmfuvdis = .true.
  lmftrac = .true.
  lmfsmooth = .false.
  lmfwstar = .false.
  lmfcapebl = .true.
  k2_const = 500.0_rkx
  kfac_shal = 0.07_rkx
  kfac_deep = 0.14_rkx

  nmlfile = 'scm_tiedtke.in'
  if ( command_argument_count( ) > 0 ) call get_command_argument(1,nmlfile)
  open(newunit=iu,file=nmlfile,status='old',action='read',iostat=ios)
  if ( ios == 0 ) then
    read(iu,nml=scmparam,iostat=ios)
    if ( ios > 0 ) call scm_error('reading scmparam in '//trim(nmlfile))
    rewind(iu)
    read(iu,nml=tiedtkeparam,iostat=ios)
    if ( ios > 0 ) call scm_error('reading tiedtkeparam in '//trim(nmlfile))
    rewind(iu)
    read(iu,nml=cldparam,iostat=ios)
    if ( ios > 0 ) call scm_error('reading cldparam in '//trim(nmlfile))
    close(iu)
  else if ( command_argument_count( ) > 0 ) then
    call scm_error('cannot open '//trim(nmlfile))
  end if

  ! Model configuration seen by the scheme
  idynamic = 3
  icup_lnd = 5
  icup_ocn = 5
  icup(:) = 5
  ipptls = 2
  nqx = 5
  ichem = 0
  iaerosol = 0
  iindirect = 0
  icumcloud = 1
  ntr = 0
  dx = dxkm * 1000.0_rkx

  ! Read the column profiles
  open(newunit=iu,file=infile,status='old',action='read',iostat=ios)
  if ( ios /= 0 ) call scm_error('cannot open '//trim(infile))
  do
    read(iu,'(a)') line
    if ( line(1:1) /= '#' ) exit
  end do
  read(line,*) ncol, nk
  kz = nk
  kzp1 = nk + 1
  kzm1 = nk - 1
  call set_grid(ncol)
  call memory_init
  allocate(col_name(ncol), col_date(ncol), col_lat(ncol), col_lon(ncol), &
    col_ps(ncol), col_shf(ncol), col_lhf(ncol), col_zpbl(ncol),           &
    col_cape(ncol), col_crr(ncol), col_topo(ncol), col_land(ncol))
  call getmem(cuscheme,jci1,jci2,ici1,ici2,'scm:cuscheme')
  cuscheme(:,:) = 5
  call init_mod_cumulus
  call allocate_mod_cu_tiedtke
  call allocate_m2c
  allocate(omega(ncol,1,nk))
  do j = 1, ncol
    read(iu,*) col_name(j), col_lat(j), col_lon(j), col_land(j), &
      col_date(j), col_ps(j), col_shf(j), col_lhf(j), col_zpbl(j), &
      col_cape(j), col_crr(j), col_topo(j)
    do k = 1, nk
      read(iu,*) m2c%zas(j,1,k), m2c%pas(j,1,k), m2c%tas(j,1,k), &
        m2c%qxas(j,1,k,iqv), m2c%uas(j,1,k), m2c%vas(j,1,k),       &
        omega(j,1,k), m2c%tten(j,1,k), m2c%qxten(j,1,k,iqv)
    end do
    do k = 1, nk + 1
      read(iu,*) m2c%zfs(j,1,k), m2c%pasf(j,1,k)
    end do
    ! specific humidity to mixing ratio, as in the model
    do k = 1, nk
      m2c%qxas(j,1,k,iqv) = m2c%qxas(j,1,k,iqv)/(d_one-m2c%qxas(j,1,k,iqv))
    end do
    m2c%hfx(j,1) = col_shf(j)
    m2c%qfx(j,1) = col_lhf(j)/wlhv
    m2c%zpbl(j,1) = col_zpbl(j)
    m2c%ldmsk(j,1) = col_land(j)
    m2c%psb(j,1) = col_ps(j)
  end do
  close(iu)
  do k = 1, nk
    pmean(k) = sum(m2c%pas(:,1,k))/real(ncol,rkx)
  end do
  avg_ww(:,:,:) = omega

  ! One call of the scheme
  cu_prate(:,:) = d_zero
  cu_srate(:,:) = d_zero
  cu_ktop(:,:) = 0
  cu_kbot(:,:) = 0
  cu_tten(:,:,:) = d_zero
  cu_qten(:,:,:,:) = d_zero
  cu_uten(:,:,:) = d_zero
  cu_vten(:,:,:) = d_zero
  cu_cldfrc(:,:,:) = d_zero
  cu_qdetr(:,:,:) = d_zero
  cu_raincc(:,:,:) = d_zero
  total_precip_points = 0
  call tiedtkedrv(m2c)

  ! Output
  open(newunit=iu,file=outfile,status='replace',action='write')
  write(iu,'(a)') '# RegCM Tiedtke single column test'
  write(iu,'(a,a)') '# profiles: ', trim(infile)
  write(iu,'(a,f10.2,a,f10.3)') '# dtcum (s) ', dtcum, '  dx (km) ', dxkm
  write(iu,nml=tiedtkeparam)
  write(iu,nml=cldparam)
  write(iu,'(a)') '# COLUMN name land ktype kctop kcbot rain(mm/day) '// &
    'snow(mm/day) heating(W/m2) moistening(W/m2) era5_cape '// &
    'era5_conv_rain(mm/day)'
  write(iu,'(a)') '# z(m) p(hPa) mfu(kg/m2/s) dT(K/day) dq(g/kg/day) '// &
    'dqc(g/kg/day) dqi(g/kg/day) du(m/s/day) dv(m/s/day) cldfrc '// &
    'qdetr(g/kg)'
  write(iu,'(a,2i6)') 'NCOL ', ncol, nk
  do j = 1, ncol
    heat = d_zero
    moist = d_zero
    do k = 1, nk
      dp = m2c%pasf(j,1,k+1) - m2c%pasf(j,1,k)
      heat = heat + cpd*cu_tten(j,1,k)*dp*regrav
      moist = moist + wlhv*cu_qten(j,1,k,iqv)*dp*regrav
    end do
    rain = cu_prate(j,1)*86400.0_rkx
    write(iu,'(a,1x,a,i3,3i4,4es14.5,2es14.5)') 'COLUMN', &
      trim(col_name(j)), col_land(j), ktype(j), cu_ktop(j,1), cu_kbot(j,1), &
      rain, cu_srate(j,1)*86400.0_rkx, heat, moist, col_cape(j), col_crr(j)
    do k = 1, nk
      write(iu,'(f10.2,f10.3,9es14.5)') m2c%zas(j,1,k), &
        m2c%pas(j,1,k)*0.01_rkx, pmfu(j,k),                      &
        cu_tten(j,1,k)*86400.0_rkx, cu_qten(j,1,k,iqv)*8.64e7_rkx, &
        cu_qten(j,1,k,iqc)*8.64e7_rkx, cu_qten(j,1,k,iqi)*8.64e7_rkx, &
        cu_uten(j,1,k)*86400.0_rkx, cu_vten(j,1,k)*86400.0_rkx,     &
        cu_cldfrc(j,1,k), cu_qdetr(j,1,k)*1000.0_rkx
    end do
  end do
  close(iu)
  write(*,'(a,i4,a,a)') ' scm_tiedtke: ', ncol, ' columns written to ', &
    trim(outfile)

  contains

  subroutine scm_error(msg)
    implicit none
    character(len=*), intent(in) :: msg
    write(*,'(a,a)') ' scm_tiedtke error: ', msg
    stop 1
  end subroutine scm_error

  ! A 1 x n grid: every index range of the model domain is 1:n x 1:1
  subroutine set_grid(n)
    implicit none
    integer(ik4), intent(in) :: n
    jci1 = 1 ; jci2 = n ; ici1 = 1 ; ici2 = 1
    jce1 = 1 ; jce2 = n ; ice1 = 1 ; ice2 = 1
    jdi1 = 1 ; jdi2 = n ; idi1 = 1 ; idi2 = 1
    jce1ga = 1 ; jce2ga = n ; ice1ga = 1 ; ice2ga = 1
    jce1gb = 1 ; jce2gb = n ; ice1gb = 1 ; ice2gb = 1
  end subroutine set_grid

  subroutine allocate_m2c
    implicit none
    call getmem(m2c%pas,1,ncol,1,1,1,nk,'scm:pas')
    call getmem(m2c%pasf,1,ncol,1,1,1,nk+1,'scm:pasf')
    call getmem(m2c%zas,1,ncol,1,1,1,nk,'scm:zas')
    call getmem(m2c%zfs,1,ncol,1,1,1,nk+1,'scm:zfs')
    call getmem(m2c%tas,1,ncol,1,1,1,nk,'scm:tas')
    call getmem(m2c%uas,1,ncol,1,1,1,nk,'scm:uas')
    call getmem(m2c%vas,1,ncol,1,1,1,nk,'scm:vas')
    call getmem(m2c%qxas,1,ncol,1,1,1,nk,1,nqx,'scm:qxas')
    call getmem(m2c%tten,1,ncol,1,1,1,nk,'scm:tten')
    call getmem(m2c%uten,1,ncol,1,1,1,nk,'scm:uten')
    call getmem(m2c%vten,1,ncol,1,1,1,nk,'scm:vten')
    call getmem(m2c%qxten,1,ncol,1,1,1,nk,1,nqx,'scm:qxten')
    call getmem(m2c%qfx,1,ncol,1,1,'scm:qfx')
    call getmem(m2c%hfx,1,ncol,1,1,'scm:hfx')
    call getmem(m2c%zpbl,1,ncol,1,1,'scm:zpbl')
    call getmem(m2c%psb,1,ncol,1,1,'scm:psb')
    call getmem(m2c%ldmsk,1,ncol,1,1,'scm:ldmsk')
    m2c%qxas(:,:,:,:) = d_zero
    m2c%qxten(:,:,:,:) = d_zero
    m2c%uten(:,:,:) = d_zero
    m2c%vten(:,:,:) = d_zero
  end subroutine allocate_m2c

end program scm_tiedtke

! The library calls myabort (Main/abort.F90) on fatal errors.
subroutine myabort
  implicit none
  stop 2
end subroutine myabort
