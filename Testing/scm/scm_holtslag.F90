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

! Single column test of the Holtslag boundary layer scheme (ibltyp = 1).
!
! Usage: scm_holtslag [namelist file]      (default: scm_holtslag.in)
!
! The columns of the profile file are laid out as a 1 x ncol grid and the
! scheme is called through its model entry point, holtbl, as the MOLOCH
! core does (idynamic = 3, NoTo microphysics, no chemistry). With
! nsteps > 1 the columns are integrated in time with the boundary layer
! tendencies only, the surface fluxes held fixed and a Coriolis term
! that relaxes the wind towards its initial (geostrophic) profile.
! The namelist has a &scmparam group (files, time step, number of steps)
! and the model group &holtslagparam, with the RegCM defaults. The output
! file is plain text, read by scm_holtslag.py: a header with the namelists
! used, then per column a COLUMN line and one line per level, top to
! bottom.

program scm_holtslag
  use mod_intkinds, only : ik4
  use mod_realkinds, only : rkx
  use mod_constants, only : cpd, d_one, d_two, d_zero, degrad, eomeg, ep1, &
    p00, regrav, rgas, rkx, rovcp, vonkar, wlhv
  use mod_dynparam, only : ice1, ice1ga, ice1gb, ice2, ice2ga, ice2gb, ici1, &
    ici1ga, ici2, ici2ga, idi1, idi2, idynamic, jce1, jce1ga, jce1gb, jce2, &
    jce2ga, jce2gb, jci1, jci1ga, jci2, jci2ga, jdi1, jdi2, kz, kzm1, kzp1, &
    ntr
  use mod_runparams, only : dt, holtth10iter, ibltyp, ichem, ifaholt, &
    ifaholtth10, ipptls, iqc, iqv, lmoistn2, lsctop, maxprandtl, nqx, rdt, &
    ricr_lnd, ricr_ocn, zhnew_fac
  use mod_memutil, only : getmem, memory_init
  use mod_regcm_types, only : ik4, mod_2_pbl, pbl_2_mod, rkx
  use mod_pbl_common, only : ricr, kmxpbl
  use mod_pbl_holtbl, only : allocate_mod_pbl_holtbl, holtbl
  implicit none

  character(len=256) :: nmlfile, infile, outfile
  real(rkx) :: dtpbl
  integer(ik4) :: nsteps
  integer(ik4) :: ncol, nk, j, k, n, iu, ios, istep
  character(len=32), allocatable, dimension(:) :: col_name
  character(len=10), allocatable, dimension(:) :: col_date
  real(rkx), allocatable, dimension(:) :: col_lat, col_lon, col_ps, &
    col_shf, col_lhf, col_zpbl, col_cape, col_crr, col_topo, col_z0
  integer(ik4), allocatable, dimension(:) :: col_land
  real(rkx), allocatable, dimension(:,:) :: t0, q0, u0, v0, &
    dt1, dq1, dqc1, du1, dv1
  type(mod_2_pbl) :: m2p
  type(pbl_2_mod) :: p2m
  real(rkx) :: heat, moist, dp, wspd, ust, unew, vnew
  character(len=512) :: line

  namelist /scmparam/ infile, outfile, dtpbl, nsteps
  namelist /holtslagparam/ ricr_ocn, ricr_lnd, zhnew_fac, &
    ifaholtth10, ifaholt, holtth10iter, lmoistn2, maxprandtl, lsctop

  ! Defaults: as in Main/mod_params.F90
  infile = 'scm_profiles.txt'
  outfile = 'scm_holtslag.out'
  dtpbl = 240.0_rkx
  nsteps = 1
  ricr_ocn = 0.275_rkx
  ricr_lnd = 0.250_rkx
  zhnew_fac = 0.25_rkx
  ifaholt = 0
  ifaholtth10 = 2
  holtth10iter = 1
  lmoistn2 = .false.
  maxprandtl = 1.0_rkx
  lsctop = .true.

  nmlfile = 'scm_holtslag.in'
  if ( command_argument_count( ) > 0 ) call get_command_argument(1,nmlfile)
  open(newunit=iu,file=nmlfile,status='old',action='read',iostat=ios)
  if ( ios == 0 ) then
    read(iu,nml=scmparam,iostat=ios)
    if ( ios > 0 ) call scm_error('reading scmparam in '//trim(nmlfile))
    rewind(iu)
    read(iu,nml=holtslagparam,iostat=ios)
    if ( ios > 0 ) call scm_error('reading holtslagparam in '//trim(nmlfile))
    close(iu)
  else if ( command_argument_count( ) > 0 ) then
    call scm_error('cannot open '//trim(nmlfile))
  end if
  nsteps = max(nsteps,1)

  ! Model configuration seen by the scheme
  idynamic = 3
  ibltyp = 1
  ipptls = 2
  nqx = 5
  ichem = 0
  ntr = 0
  dt = dtpbl
  rdt = d_one/dt

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
    col_cape(ncol), col_crr(ncol), col_topo(ncol), col_land(ncol),        &
    col_z0(ncol))
  allocate(t0(ncol,nk), q0(ncol,nk), u0(ncol,nk), v0(ncol,nk), &
    dt1(ncol,nk), dq1(ncol,nk), dqc1(ncol,nk), du1(ncol,nk), dv1(ncol,nk))
  call getmem(ricr,jci1,jci2,ici1,ici2,'scm:ricr')
  call getmem(kmxpbl,jci1,jci2,ici1,ici2,'scm:kmxpbl')
  call allocate_mod_pbl_holtbl
  call allocate_m2p
  do j = 1, ncol
    read(iu,*) col_name(j), col_lat(j), col_lon(j), col_land(j), &
      col_date(j), col_ps(j), col_shf(j), col_lhf(j), col_zpbl(j), &
      col_cape(j), col_crr(j), col_topo(j)
    do k = 1, nk
      read(iu,*) m2p%za(j,1,k), m2p%patm(j,1,k), m2p%tatm(j,1,k), &
        m2p%qxatm(j,1,k,iqv), m2p%uxatm(j,1,k), m2p%vxatm(j,1,k)
    end do
    do k = 1, nk + 1
      read(iu,*) dp, m2p%patmf(j,1,k)
    end do
    ! specific humidity to mixing ratio, as in the model
    do k = 1, nk
      m2p%qxatm(j,1,k,iqv) = m2p%qxatm(j,1,k,iqv)/(d_one-m2p%qxatm(j,1,k,iqv))
    end do
    m2p%hfx(j,1) = col_shf(j)
    m2p%qfx(j,1) = col_lhf(j)/wlhv
    m2p%ldmsk(j,1) = col_land(j)
    m2p%coriol(j,1) = d_two*eomeg*sin(col_lat(j)*degrad)
    ! Roughness length 0.1 m on land, 1 mm on the ocean
    col_z0(j) = merge(0.1_rkx,1.0e-3_rkx,col_land(j) == 1)
    ! Critical Richardson number and highest boundary layer level, as in
    ! mod_params and mod_slice
    ricr(j,1) = merge(ricr_lnd,ricr_ocn,col_land(j) == 1)
    kmxpbl(j,1) = kz
    do k = kzm1, 2, -1
      if ( m2p%za(j,1,k) >= 5500.0_rkx ) exit
      kmxpbl(j,1) = k
    end do
  end do
  close(iu)
  t0 = m2p%tatm(:,1,:)
  q0 = m2p%qxatm(:,1,:,iqv)
  u0 = m2p%uxatm(:,1,:)
  v0 = m2p%vxatm(:,1,:)

  do istep = 1, nsteps
    call update_derived
    p2m%tten(:,:,:) = d_zero
    p2m%uten(:,:,:) = d_zero
    p2m%vten(:,:,:) = d_zero
    p2m%qxten(:,:,:,:) = d_zero
    call holtbl(m2p,p2m)
    if ( istep == 1 ) then
      dt1 = p2m%tten(:,1,:)
      dq1 = p2m%qxten(:,1,:,iqv)
      dqc1 = p2m%qxten(:,1,:,iqc)
      du1 = p2m%uten(:,1,:)
      dv1 = p2m%vten(:,1,:)
    end if
    if ( nsteps > 1 ) then
      ! Advance the state: boundary layer tendencies, and the Coriolis
      ! force with the initial wind as the geostrophic wind
      do j = 1, ncol
        do k = 1, nk
          m2p%tatm(j,1,k) = m2p%tatm(j,1,k) + dt*p2m%tten(j,1,k)
          do n = 1, nqx
            m2p%qxatm(j,1,k,n) = max(m2p%qxatm(j,1,k,n) + &
              dt*p2m%qxten(j,1,k,n),d_zero)
          end do
          unew = m2p%uxatm(j,1,k) + dt*(p2m%uten(j,1,k) + &
            m2p%coriol(j,1)*(m2p%vxatm(j,1,k)-v0(j,k)))
          vnew = m2p%vxatm(j,1,k) + dt*(p2m%vten(j,1,k) - &
            m2p%coriol(j,1)*(m2p%uxatm(j,1,k)-u0(j,k)))
          m2p%uxatm(j,1,k) = unew
          m2p%vxatm(j,1,k) = vnew
        end do
      end do
    end if
  end do

  ! Output
  open(newunit=iu,file=outfile,status='replace',action='write')
  write(iu,'(a)') '# RegCM Holtslag single column test'
  write(iu,'(a,a)') '# profiles: ', trim(infile)
  write(iu,'(a,f10.2,a,i8)') '# dt (s) ', dt, '  nsteps ', nsteps
  write(iu,nml=holtslagparam)
  write(iu,'(a)') '# COLUMN name land kpbl zpbl(m) ustar(m/s) '// &
    'heating(W/m2) moistening(W/m2) shf(W/m2) lhf(W/m2) era5_zpbl(m)'
  write(iu,'(a)') '# z(m) p(hPa) dT(K/day) dq(g/kg/day) dqc(g/kg/day) '// &
    'du(m/s/day) dv(m/s/day) T-T0(K) q-q0(g/kg) u-u0(m/s) v-v0(m/s)'
  write(iu,'(a)') '# tendencies of the first call, differences after '// &
    'nsteps (zpbl, kpbl and ustar of the last call)'
  write(iu,'(a,2i6)') 'NCOL ', ncol, nk
  do j = 1, ncol
    heat = d_zero
    moist = d_zero
    do k = 1, nk
      dp = m2p%patmf(j,1,k+1) - m2p%patmf(j,1,k)
      heat = heat + cpd*dt1(j,k)*dp*regrav
      moist = moist + wlhv*dq1(j,k)*dp*regrav
    end do
    wspd = sqrt(m2p%uxatm(j,1,nk)**2+m2p%vxatm(j,1,nk)**2)
    ust = sqrt(m2p%uvdrag(j,1)*wspd/m2p%rhox2d(j,1))
    write(iu,'(a,1x,a,i3,i4,f10.2,7es14.5)') 'COLUMN', &
      trim(col_name(j)), col_land(j), p2m%kpbl(j,1), p2m%zpbl(j,1), &
      ust, heat, moist, col_shf(j), col_lhf(j), col_zpbl(j)
    do k = 1, nk
      write(iu,'(f10.2,f10.3,9es14.5)') m2p%za(j,1,k),      &
        m2p%patm(j,1,k)*0.01_rkx, dt1(j,k)*86400.0_rkx,       &
        dq1(j,k)*8.64e7_rkx, dqc1(j,k)*8.64e7_rkx,            &
        du1(j,k)*86400.0_rkx, dv1(j,k)*86400.0_rkx,           &
        m2p%tatm(j,1,k)-t0(j,k), (m2p%qxatm(j,1,k,iqv)-q0(j,k))*1.0e3_rkx, &
        m2p%uxatm(j,1,k)-u0(j,k), m2p%vxatm(j,1,k)-v0(j,k)
    end do
  end do
  close(iu)
  write(*,'(a,i4,a,a)') ' scm_holtslag: ', ncol, ' columns written to ', &
    trim(outfile)

  contains

  subroutine scm_error(msg)
    implicit none
    character(len=*), intent(in) :: msg
    write(*,'(a,a)') ' scm_holtslag error: ', msg
    stop 1
  end subroutine scm_error

  ! A 1 x n grid: every index range of the model domain is 1:n x 1:1
  subroutine set_grid(n)
    implicit none
    integer(ik4), intent(in) :: n
    jci1 = 1 ; jci2 = n ; ici1 = 1 ; ici2 = 1
    jce1 = 1 ; jce2 = n ; ice1 = 1 ; ice2 = 1
    jdi1 = 1 ; jdi2 = n ; idi1 = 1 ; idi2 = 1
    jci1ga = 1 ; jci2ga = n ; ici1ga = 1 ; ici2ga = 1
    jce1ga = 1 ; jce2ga = n ; ice1ga = 1 ; ice2ga = 1
    jce1gb = 1 ; jce2gb = n ; ice1gb = 1 ; ice2gb = 1
  end subroutine set_grid

  ! Fields the model computes from the prognostic state: potential
  ! temperature, surface air density, ground temperature (same potential
  ! temperature as the lowest level, used only with ifaholtth10 /= 2) and
  ! the surface drag coefficient from a neutral log law friction velocity.
  subroutine update_derived
    implicit none
    real(rkx) :: ps, tv, zl
    do j = 1, ncol
      do k = 1, nk
        m2p%thatm(j,1,k) = m2p%tatm(j,1,k) * &
          (p00/m2p%patm(j,1,k))**rovcp
      end do
      ps = m2p%patmf(j,1,nk+1)
      tv = m2p%tatm(j,1,nk)*(d_one+ep1*m2p%qxatm(j,1,nk,iqv) / &
        (d_one+m2p%qxatm(j,1,nk,iqv)))
      m2p%rhox2d(j,1) = ps/(rgas*tv)
      m2p%tg(j,1) = m2p%thatm(j,1,nk)*(ps/p00)**rovcp
      wspd = max(sqrt(m2p%uxatm(j,1,nk)**2+m2p%vxatm(j,1,nk)**2),0.1_rkx)
      zl = max(m2p%za(j,1,nk),d_two*col_z0(j))
      ust = vonkar*wspd/log(zl/col_z0(j))
      m2p%uvdrag(j,1) = m2p%rhox2d(j,1)*ust*ust/wspd
    end do
  end subroutine update_derived

  subroutine allocate_m2p
    implicit none
    call getmem(m2p%za,1,ncol,1,1,1,nk,'scm:za')
    call getmem(m2p%heatrt,1,ncol,1,1,1,nk,'scm:pheatrt')
    m2p%heatrt(:,:,:) = d_zero
    ! large scale cloud fraction of the model (cldfrac)
    call getmem(m2p%fcc,1,ncol,1,1,1,nk,'scm:pfcc')
    m2p%fcc(:,:,:) = d_zero
    call getmem(m2p%patm,1,ncol,1,1,1,nk,'scm:patm')
    call getmem(m2p%patmf,1,ncol,1,1,1,nk+1,'scm:patmf')
    call getmem(m2p%tatm,1,ncol,1,1,1,nk,'scm:tatm')
    call getmem(m2p%thatm,1,ncol,1,1,1,nk,'scm:thatm')
    call getmem(m2p%uxatm,1,ncol,1,1,1,nk,'scm:uxatm')
    call getmem(m2p%vxatm,1,ncol,1,1,1,nk,'scm:vxatm')
    call getmem(m2p%qxatm,1,ncol,1,1,1,nk,1,nqx,'scm:qxatm')
    call getmem(m2p%hfx,1,ncol,1,1,'scm:hfx')
    call getmem(m2p%qfx,1,ncol,1,1,'scm:qfx')
    call getmem(m2p%uvdrag,1,ncol,1,1,'scm:uvdrag')
    call getmem(m2p%rhox2d,1,ncol,1,1,'scm:rhox2d')
    call getmem(m2p%tg,1,ncol,1,1,'scm:tg')
    call getmem(m2p%coriol,1,ncol,1,1,'scm:coriol')
    call getmem(m2p%ldmsk,1,ncol,1,1,'scm:ldmsk')
    call getmem(p2m%tten,1,ncol,1,1,1,nk,'scm:tten')
    call getmem(p2m%uten,1,ncol,1,1,1,nk,'scm:uten')
    call getmem(p2m%vten,1,ncol,1,1,1,nk,'scm:vten')
    call getmem(p2m%qxten,1,ncol,1,1,1,nk,1,nqx,'scm:qxten')
    call getmem(p2m%zpbl,1,ncol,1,1,'scm:zpbl')
    call getmem(p2m%kpbl,1,ncol,1,1,'scm:kpbl')
    m2p%qxatm(:,:,:,:) = d_zero
  end subroutine allocate_m2p

end program scm_holtslag

! The library calls myabort (Main/abort.F90) on fatal errors.
subroutine myabort
  implicit none
  stop 2
end subroutine myabort
