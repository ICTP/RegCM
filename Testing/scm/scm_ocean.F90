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

! Stand alone test of the ocean surface flux schemes: Zeng (iocnflux = 1,
! with izengrough) and COARE 3.5 (iocnflux = 2).
!
! Usage: scm_ocean [namelist file]         (default: scm_ocean.in)
!
! The samples of the input file are the ocean points of the scheme, called
! once through its model entry point (zengocndrv, coare3_drv) on
! open water (no sea ice, no lake, no diurnal SST scheme, no wave model).
! The inputs are what the model passes from its lowest level and from the
! surface (mod_ocn_common): height of the level, wind, temperature and
! mixing ratio at the level, pressure at the level and at the surface, the
! bulk SST, the PBL height and the net shortwave and longwave radiation at
! the surface. As in mod_slice, the surface air temperature passed to the
! schemes (sfta) is the temperature of the level brought dry adiabatically
! to the surface pressure, the air density the surface pressure over
! R T of the level. The 3-day mean 10 m wind (used by izengrough 3)
! is the wind of the level.
!
! Input file: comment lines starting with #, then the number of samples,
! then one line per sample:
!   z(m) u v(m/s) T(K) qv(kg/kg, mixing ratio) p(Pa) ps(Pa) sst(K)
!   zpbl(m) swnet lwnet(W/m2, lwnet upward positive)
! Output file: a header, then one line per sample:
!   shf lhf(W/m2, upward) tau(N/m2) ustar(m/s) u10(m/s) t2m(K)
!   q2m(g/kg) drag(kg/m2/s)
! The latent heat flux is the evaporation times the constant wlhv, the
! same for all the schemes.

program scm_ocean
  use mod_intkinds, only : ik4
  use mod_realkinds, only : rkx
  use mod_constants, only : d_zero, rgas, rkx, rovcp, wlhv
  use mod_runparams, only : iocnflux, iwavcpl, izengrough
  use mod_memutil, only : getmem, memory_init
  use mod_ocn_internal, only : br, drag, dtocn, evpr, getmem, hpbl, ht, icpl, &
    ik4, iocnbeg, iocnend, iocnflux, ldcsst, llake, lseaice, mask, &
    memory_init, nocnp, patm, q2m, qv, rah1, ram1, rhoa, rhox, rkx, rlwf, &
    rswf, sent, sfice, sfps, sfta, t2m, tatm, taux, tauy, tgb, tgbrd, tgrd, &
    u10m, um10, ustr, usw, v10m, vsw, zoo
  use mod_ocn_zeng, only : zengocndrv
  use mod_ocn_coare, only : coare3_drv
  implicit none

  character(len=256) :: nmlfile, infile, outfile
  integer(ik4) :: n, i, iu, ios
  character(len=512) :: line
  real(rkx) :: wspd

  namelist /scmparam/ infile, outfile, iocnflux, izengrough

  infile = 'scm_ocean.txt'
  outfile = 'scm_ocean.out'
  ! Defaults: as in Main/mod_params.F90
  iocnflux = 1
  izengrough = 1

  nmlfile = 'scm_ocean.in'
  if ( command_argument_count( ) > 0 ) call get_command_argument(1,nmlfile)
  open(newunit=iu,file=nmlfile,status='old',action='read',iostat=ios)
  if ( ios == 0 ) then
    read(iu,nml=scmparam,iostat=ios)
    if ( ios > 0 ) call scm_error('reading scmparam in '//trim(nmlfile))
    close(iu)
  else if ( command_argument_count( ) > 0 ) then
    call scm_error('cannot open '//trim(nmlfile))
  end if

  iwavcpl = 0
  dtocn = 600.0_rkx
  ldcsst = .false.
  llake = .false.
  lseaice = .false.

  open(newunit=iu,file=infile,status='old',action='read',iostat=ios)
  if ( ios /= 0 ) call scm_error('cannot open '//trim(infile))
  do
    read(iu,'(a)') line
    if ( line(1:1) /= '#' ) exit
  end do
  read(line,*) n
  call memory_init
  nocnp = n
  iocnbeg = 1
  iocnend = n
  call allocate_fields
  do i = 1, n
    read(iu,*) ht(i), usw(i), vsw(i), tatm(i), qv(i), patm(i), &
      sfps(i), tgb(i), hpbl(i), rswf(i), rlwf(i)
  end do
  close(iu)
  ! as in mod_ocn_common and mod_slice
  usw = sign(1.0_rkx,usw)*max(abs(usw),0.001_rkx)
  vsw = sign(1.0_rkx,vsw)*max(abs(vsw),0.001_rkx)
  do i = 1, n
    sfta(i) = tatm(i)*(sfps(i)/patm(i))**rovcp
    rhox(i) = sfps(i)/(rgas*tatm(i))
    um10(i) = sqrt(usw(i)**2+vsw(i)**2)
  end do
  mask(:) = 1
  tgrd(:) = tgb(:)
  tgbrd(:) = tgb(:)

  select case ( iocnflux )
    case (1)
      call zengocndrv
    case (2)
      call coare3_drv
    case default
      call scm_error('iocnflux must be 1 or 2')
  end select

  open(newunit=iu,file=outfile,status='replace',action='write')
  write(iu,'(a,2i3)') '# RegCM ocean flux test, iocnflux izengrough ', &
    iocnflux, izengrough
  write(iu,'(a)') '# shf lhf(W/m2) tau(N/m2) ustar(m/s) u10(m/s) '// &
    't2m(K) q2m(g/kg) drag(kg/m2/s)'
  write(iu,'(i10)') n
  do i = 1, n
    wspd = sqrt(u10m(i)**2+v10m(i)**2)
    write(iu,'(8es14.6)') sent(i), evpr(i)*wlhv, &
      sqrt(taux(i)**2+tauy(i)**2), ustr(i), wspd, t2m(i), &
      q2m(i)*1000.0_rkx, drag(i)
  end do
  close(iu)
  write(*,'(a,i8,a,a)') ' scm_ocean: ', n, ' samples written to ', &
    trim(outfile)

  contains

  subroutine scm_error(msg)
    implicit none
    character(len=*), intent(in) :: msg
    write(*,'(a,a)') ' scm_ocean error: ', msg
    stop 1
  end subroutine scm_error

  subroutine allocate_fields
    implicit none
    call getmem(mask,1,n,'scm:mask')
    call getmem(icpl,1,n,'scm:icpl')
    call getmem(ht,1,n,'scm:ht')
    call getmem(usw,1,n,'scm:usw')
    call getmem(vsw,1,n,'scm:vsw')
    call getmem(tatm,1,n,'scm:tatm')
    call getmem(qv,1,n,'scm:qv')
    call getmem(patm,1,n,'scm:patm')
    call getmem(sfps,1,n,'scm:sfps')
    call getmem(sfta,1,n,'scm:sfta')
    call getmem(tgb,1,n,'scm:tgb')
    call getmem(tgrd,1,n,'scm:tgrd')
    call getmem(tgbrd,1,n,'scm:tgbrd')
    call getmem(hpbl,1,n,'scm:hpbl')
    call getmem(rswf,1,n,'scm:rswf')
    call getmem(rlwf,1,n,'scm:rlwf')
    call getmem(rhox,1,n,'scm:rhox')
    call getmem(um10,1,n,'scm:um10')
    call getmem(sent,1,n,'scm:sent')
    call getmem(evpr,1,n,'scm:evpr')
    call getmem(drag,1,n,'scm:drag')
    call getmem(ustr,1,n,'scm:ustr')
    call getmem(zoo,1,n,'scm:zoo')
    call getmem(taux,1,n,'scm:taux')
    call getmem(tauy,1,n,'scm:tauy')
    call getmem(u10m,1,n,'scm:u10m')
    call getmem(v10m,1,n,'scm:v10m')
    call getmem(t2m,1,n,'scm:t2m')
    call getmem(q2m,1,n,'scm:q2m')
    call getmem(rhoa,1,n,'scm:rhoa')
    call getmem(ram1,1,n,'scm:ram1')
    call getmem(rah1,1,n,'scm:rah1')
    call getmem(br,1,n,'scm:br')
    call getmem(sfice,1,n,'scm:sfice')
    icpl(:) = 0
    sfice(:) = d_zero
  end subroutine allocate_fields

end program scm_ocean

! The library calls myabort (Main/abort.F90) on fatal errors.
subroutine myabort
  implicit none
  stop 2
end subroutine myabort
