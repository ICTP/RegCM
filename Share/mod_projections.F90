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

module mod_projections

  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_memutil

#ifdef F2008
  use, intrinsic :: iso_fortran_env
#endif

  implicit none

  private

  public :: rounder

  type anyprojparams
    character(len=6) :: pcode
    real(rk8) :: ds
    real(rk8) :: clat, clon
    real(rk8) :: plat, plon
    real(rk8) :: trlat1, trlat2
    integer(ik4) :: nlat, nlon
    logical :: staggerx, staggery
    logical :: rotparam
  end type anyprojparams

  type regcm_projdata
    character(len=6) :: code
    real(rk8) :: stdlon, stdlat
    real(rk8) :: truelat1, truelat2, tl2r, ctl1r
    real(rk8) :: chi1, chi2, xct1, tanchi1h, schi1, tchi1
    real(rk8) :: colat1, colat2, nfac
    real(rk8) :: rsw, rebydx, hemi
    real(rk8) :: reflon, dlon, dlat, scale_top
    real(rk8) :: polei, polej
    real(rk8) :: rlon0, rlat0, phi0, lam0
    real(rk8) :: xoff, yoff
    real(rk8) :: zsinpol, zcospol, zlampol
    real(rk8) :: pollam, polcphi, polsphi
    real(rk8) :: conefac = 1.0_rk8
    real(rk8) :: rconefac
    integer(ik4) :: nlat, nlon
    logical :: lamtan
    logical :: skiprot
  end type regcm_projdata

  type regcm_projection
    private
    type(regcm_projdata), pointer :: p => null()
    procedure(transform), pass(pj), pointer, public :: llij => NULL()
    procedure(transform), pass(pj), pointer, public :: ijll => NULL()
    procedure(map_factor), pass(pj), pointer, public :: mapfac => NULL()
    procedure(rotate2), pass(pj), pointer, public :: uvrotate2 => NULL()
    procedure(rotate2), pass(pj), pointer, public :: uvbkrotate2 => NULL()
    procedure(rotate3), pass(pj), pointer, public :: uvrotate3 => NULL()
    procedure(rotate3), pass(pj), pointer, public :: uvbkrotate3 => NULL()

    real(rk8), pointer, contiguous, dimension(:,:) :: f1 => null()
    real(rk8), pointer, contiguous, dimension(:,:) :: f2 => null()
    real(rk8), pointer, contiguous, dimension(:,:) :: f3 => null()
    real(rk8), pointer, contiguous, dimension(:,:) :: f4 => null()

    contains

      procedure, public :: initialize
      procedure, public :: rotation_angle
      procedure, public :: wind_rotate
      procedure, public :: wind2_rotate
      procedure, public :: wind_antirotate
      procedure, public :: wind2_antirotate
      procedure, public :: rl00
      procedure, public :: conefac
      procedure, public :: destruct

  end type regcm_projection

  abstract interface
    pure subroutine transform(pj,a,b,c,d)
      import
      implicit none
      class(regcm_projection), intent(in) :: pj
#ifdef SINGLE_PRECISION_REAL
#ifdef F2008
      integer, parameter :: rk4  = REAL32
#else
      integer, parameter :: rk4  = kind(1.0)
#endif
      integer, parameter :: rkx = rk4
#else
#ifdef F2008
      integer, parameter :: rk8  = REAL64
#else
      integer, parameter :: rk8  = selected_real_kind(2*precision(1.0))
#endif
      integer, parameter :: rkx = rk8
#endif
      real(rkx), intent(in) :: a, b
      real(rkx), intent(out) :: c, d
    end subroutine transform
  end interface

  abstract interface
    subroutine rotate3(pj,u,v)
      import
      implicit none
      class(regcm_projection), intent(in) :: pj
#ifdef SINGLE_PRECISION_REAL
#ifdef F2008
      integer, parameter :: rk4  = REAL32
#else
      integer, parameter :: rk4  = kind(1.0)
#endif
      integer, parameter :: rkx = rk4
#else
#ifdef F2008
      integer, parameter :: rk8  = REAL64
#else
      integer, parameter :: rk8  = selected_real_kind(2*precision(1.0))
#endif
      integer, parameter :: rkx = rk8
#endif
      real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: u, v
    end subroutine rotate3
  end interface

  abstract interface
    subroutine rotate2(pj,u,v)
      import
      implicit none
      class(regcm_projection), intent(in) :: pj
#ifdef SINGLE_PRECISION_REAL
#ifdef F2008
      integer, parameter :: rk4  = REAL32
#else
      integer, parameter :: rk4  = kind(1.0)
#endif
      integer, parameter :: rkx = rk4
#else
#ifdef F2008
      integer, parameter :: rk8  = REAL64
#else
      integer, parameter :: rk8  = selected_real_kind(2*precision(1.0))
#endif
      integer, parameter :: rkx = rk8
#endif
      real(rkx), pointer, contiguous, dimension(:,:), intent(inout) :: u, v
    end subroutine rotate2
  end interface

  abstract interface
    pure subroutine map_factor(pj,a,b,c)
      import
      implicit none
      class(regcm_projection), intent(in) :: pj
#ifdef SINGLE_PRECISION_REAL
#ifdef F2008
      integer, parameter :: rk4  = REAL32
#else
      integer, parameter :: rk4  = kind(1.0)
#endif
      integer, parameter :: rkx = rk4
#else
#ifdef F2008
      integer, parameter :: rk8  = REAL64
#else
      integer, parameter :: rk8  = selected_real_kind(2*precision(1.0))
#endif
      integer, parameter :: rkx = rk8
#endif
      real(rkx), pointer, contiguous, dimension(:,:), intent(in) :: a, b
      real(rkx), pointer, contiguous, dimension(:,:), intent(inout) :: c
    end subroutine map_factor
  end interface

  public :: anyprojparams, regcm_projection

  contains

  subroutine rl00(pj,lat0,lon0)
    implicit none
    class(regcm_projection), intent(in) :: pj
    real(rkx), intent(out) :: lat0, lon0
    lat0 = real(pj%p%rlat0,rkx)
    lon0 = real(pj%p%rlon0,rkx)
  end subroutine rl00

  real(rk8) function conefac(pj)
    implicit none
    class(regcm_projection), intent(in) :: pj
    conefac = pj%p%conefac
  end function conefac

  subroutine wind_rotate(pj,u,v)
    implicit none
    class(regcm_projection), intent(in) :: pj
    real(rkx), dimension(:,:,:), pointer, contiguous, intent(inout) :: u, v
    select case(pj%p%code)
      case ('LAMCON')
        call rotate3_lc(pj,u,v)
      case ('ROTMER')
        call rotate3_rc(pj,u,v)
      case ('POLSTR')
        call rotate3_ps(pj,u,v)
      case ('ROTLLR')
        call rotate3_rl(pj,u,v)
      case default
        ! No action
    end select
  end subroutine wind_rotate

  subroutine wind_antirotate(pj,u,v)
    implicit none
    class(regcm_projection), intent(in) :: pj
    real(rkx), dimension(:,:,:), pointer, contiguous, intent(inout) :: u, v
    select case(pj%p%code)
      case ('LAMCON')
        call backrotate3_lc(pj,u,v)
      case ('ROTMER')
        call backrotate3_rc(pj,u,v)
      case ('POLSTR')
        call backrotate3_ps(pj,u,v)
      case ('ROTLLR')
        call backrotate3_rl(pj,u,v)
      case default
        ! No action
    end select
  end subroutine wind_antirotate

  subroutine wind2_rotate(pj,u,v)
    implicit none
    class(regcm_projection), intent(in) :: pj
    real(rkx), dimension(:,:), pointer, contiguous, intent(inout) :: u, v
    select case(pj%p%code)
      case ('LAMCON')
        call rotate2_lc(pj,u,v)
      case ('ROTMER')
        call rotate2_rc(pj,u,v)
      case ('POLSTR')
        call rotate2_ps(pj,u,v)
      case ('ROTLLR')
        call rotate2_rl(pj,u,v)
      case default
        ! No action
    end select
  end subroutine wind2_rotate

  subroutine wind2_antirotate(pj,u,v)
    implicit none
    class(regcm_projection), intent(in) :: pj
    real(rkx), dimension(:,:), pointer, contiguous, intent(inout) :: u, v
    select case(pj%p%code)
      case ('LAMCON')
        call backrotate2_lc(pj,u,v)
      case ('ROTMER')
        call backrotate2_rc(pj,u,v)
      case ('POLSTR')
        call backrotate2_ps(pj,u,v)
      case ('ROTLLR')
        call backrotate2_rl(pj,u,v)
      case default
        ! No action
    end select
  end subroutine wind2_antirotate

  subroutine rotation_angle(pj,lon,lat,p)
    implicit none
    class(regcm_projection), intent(in) :: pj
    real(rkx), dimension(:), intent(in) :: lon
    real(rkx), dimension(:), intent(in) :: lat
    real(rkx), dimension(:,:), intent(out) :: p
    integer(ik4) :: i, j, nlon, nlat
    real(rk8) :: f1, f2
    nlon = size(lon)
    nlat = size(lat)
    if ( size(p,1) /= nlon .or. size(p,2) /= nlat ) then
      return
    end if
    select case (pj%p%code)
      case ('LAMCON')
        do i = 1, nlat
          do j = 1, nlon
            p(j,i) = real(uvrot_lc(pj,lon(j)),rkx)
          end do
        end do
      case ('ROTMER')
        do i = 1, nlat
          do j = 1, nlon
            call uvrot_rc(pj,lat(i),lon(j),f1,f2)
            p(j,i) = real(atan2(f2,f1),rkx)
          end do
        end do
      case ('POLSTR')
        do i = 1, nlat
          do j = 1, nlon
            p(j,i) = real(uvrot_ps(pj,lon(j)),rkx)
          end do
        end do
      case default
        p(:,:) = 0.0_rkx
    end select
  end subroutine rotation_angle

  subroutine destruct(pj)
    implicit none
    class(regcm_projection), intent(inout) :: pj
    if ( associated(pj%p) ) deallocate(pj%p)
  end subroutine destruct

  subroutine initialize(pj,pjpara)
    implicit none
    class(regcm_projection), intent(out) :: pj
    type(anyprojparams), intent(in) :: pjpara
    real(rk8) :: ci, cj, ci1, cj1
    if ( .not. associated(pj%p) ) then
      allocate(pj%p)
    end if
    !
    ! Offset of the domain centre from the first grid point: the centre is
    ! in the middle of the grid, at index ci+1, cj+1 (half a grid point
    ! further for a staggered grid). The projections working with the grid
    ! index (LAMCON, POLSTR, NORMER, ROTMER) take it as ci1, cj1; the
    ! regular latitude-longitude grids take the offset ci, cj.
    !
    ci = real(pjpara%nlon,rk8) * 0.5_rk8
    cj = real(pjpara%nlat,rk8) * 0.5_rk8
    if ( .not. pjpara%staggerx ) ci = ci - 0.5_rk8
    if ( .not. pjpara%staggery ) cj = cj - 0.5_rk8
    ci1 = ci + 1.0_rk8
    cj1 = cj + 1.0_rk8
    pj%p%code = pjpara%pcode
    pj%p%nlon = pjpara%nlon
    pj%p%nlat = pjpara%nlat
    select case (pjpara%pcode)
      case ('LAMCON')
        call setup_lcc(pj,pjpara%clon,pjpara%clat,ci1,cj1,pjpara%ds, &
                       pjpara%clon,pjpara%trlat1,pjpara%trlat2,pjpara%rotparam)
        pj%llij => llij_lc
        pj%ijll => ijll_lc
        pj%uvrotate2 => rotate2_lc
        pj%uvrotate3 => rotate3_lc
        pj%uvbkrotate2 => backrotate2_lc
        pj%uvbkrotate3 => backrotate3_lc
        pj%mapfac => mapfac_lc
      case ('POLSTR')
        call setup_plr(pj,pjpara%clon,pjpara%clat,ci1,cj1, &
                       pjpara%ds,pjpara%clon,pjpara%rotparam)
        pj%llij => llij_ps
        pj%ijll => ijll_ps
        pj%uvrotate2 => rotate2_ps
        pj%uvrotate3 => rotate3_ps
        pj%uvbkrotate2 => backrotate2_ps
        pj%uvbkrotate3 => backrotate3_ps
        pj%mapfac => mapfac_ps
      case ('NORMER')
        call setup_mrc(pj,pjpara%clon,pjpara%clat,ci1,cj1, &
                       pjpara%ds)
        pj%llij => llij_mc
        pj%ijll => ijll_mc
        pj%uvrotate2 => rotate2_mc
        pj%uvrotate3 => rotate3_mc
        pj%uvbkrotate2 => rotate2_mc
        pj%uvbkrotate3 => rotate3_mc
        pj%mapfac => mapfac_mc
      case ('ROTMER')
        call setup_rmc(pj,pjpara%clon,pjpara%clat,ci1,cj1, &
                       pjpara%ds,pjpara%plon,pjpara%plat,pjpara%rotparam)
        pj%llij => llij_rc
        pj%ijll => ijll_rc
        pj%uvrotate2 => rotate2_rc
        pj%uvrotate3 => rotate3_rc
        pj%uvbkrotate2 => backrotate2_rc
        pj%uvbkrotate3 => backrotate3_rc
        pj%mapfac => mapfac_rc
      case ('ROTLLR')
        call setup_rll(pj,pjpara%clon,pjpara%clat,ci,cj, &
                       pjpara%ds,pjpara%plon,pjpara%plat,pjpara%rotparam)
        pj%llij => llij_rl
        pj%ijll => ijll_rl
        pj%uvrotate2 => rotate2_rl
        pj%uvrotate3 => rotate3_rl
        pj%uvbkrotate2 => backrotate2_rl
        pj%uvbkrotate3 => backrotate3_rl
        pj%mapfac => mapfac_rl
      case default
        call setup_ll(pj,pjpara%clon,pjpara%clat,ci,cj,pjpara%ds)
        pj%llij => llij_ll
        pj%ijll => ijll_ll
        pj%uvrotate2 => rotate2_mc
        pj%uvrotate3 => rotate3_mc
        pj%uvbkrotate2 => rotate2_mc
        pj%uvbkrotate3 => rotate3_mc
        pj%mapfac => mapfac_ll
    end select
  end subroutine initialize

  subroutine setup_ll(pj,clon,clat,ci,cj,ds)
    implicit none
    type(regcm_projection), intent(inout) :: pj
    real(rk8), intent(in) :: ci, cj, clat, clon, ds
    pj%p%dlon = raddeg_r8 * ds / earthrad_r8
    pj%p%dlat = pj%p%dlon
    pj%p%rlat0 = clat - cj*pj%p%dlat
    pj%p%rlon0 = clon - ci*pj%p%dlon
    if ( pj%p%rlat0 >  deg90 )  pj%p%rlat0 =  deg180 - pj%p%rlat0
    if ( pj%p%rlat0 < -deg90 )  pj%p%rlat0 = -deg180 - pj%p%rlat0
    if ( pj%p%rlon0 >  deg180 ) pj%p%rlon0 = pj%p%rlon0 - deg360
    if ( pj%p%rlon0 < -deg180 ) pj%p%rlon0 = pj%p%rlon0 + deg360
  end subroutine setup_ll

  pure subroutine ijll_ll(pj,i,j,lat,lon)
    implicit none
    class(regcm_projection), intent(in) :: pj
    real(rkx), intent(in) :: i, j
    real(rkx), intent(out) :: lat, lon
    lat = real(pj%p%rlat0 + (j-1.0_rk8) * pj%p%dlat,rkx)
    lon = real(pj%p%rlon0 + (i-1.0_rk8) * pj%p%dlon,rkx)
    ! Across a pole: reflect the latitude, go to the opposite meridian
    if ( lat >  90.0_rkx ) then
      lat =  180.0_rkx - lat
      lon = lon + 180.0_rkx
    end if
    if ( lat < -90.0_rkx ) then
      lat = -180.0_rkx - lat
      lon = lon + 180.0_rkx
    end if
    if ( lon >  180.0_rkx ) lon = lon - 360.0_rkx
    if ( lon < -180.0_rkx ) lon = lon + 360.0_rkx
  end subroutine ijll_ll

  pure subroutine llij_ll(pj,lat,lon,i,j)
    implicit none
    class(regcm_projection), intent(in) :: pj
    real(rkx), intent(in) :: lat, lon
    real(rkx), intent(out) :: i, j
    i = real(( lon - pj%p%rlon0 ) / pj%p%dlon + 1.0_rk8,rkx)
    j = real(( lat - pj%p%rlat0 ) / pj%p%dlat + 1.0_rk8,rkx)
  end subroutine llij_ll

  subroutine setup_rll(pj,clon,clat,ci,cj,ds,plon,plat,luvrot)
    implicit none
    type(regcm_projection), intent(inout) :: pj
    real(rk8), intent(in) :: ci, cj, clat, clon, plon, plat, ds
    logical, intent(in) :: luvrot
    real(rk8) :: phi, lam, plam, pphi, clam, cphi, dlam
    real(rk8) :: delta, arg1, arg2
    real(rkx) :: lon, lat, ri, rj
    integer(ik4) :: i, j
    pj%p%skiprot = .false.
    if ( abs(plat-90.0_rk8) < 0.001 ) pj%p%skiprot = .true.
    pj%p%dlon = raddeg_r8 * ds / earthrad_r8
    pj%p%dlat = pj%p%dlon
    pphi = degrad_r8*plat
    plam = degrad_r8*plon
    call get_equator(pphi,plam,pj%p%phi0,pj%p%lam0)
    cphi = degrad_r8 * clat
    clam = degrad_r8 * clon
    call geo_to_rotated(pj%p%phi0,pj%p%lam0,cphi,clam, &
                        pj%p%rlat0,pj%p%rlon0)
    ! Rotated longitude of the domain centre
    pj%p%xoff = raddeg_r8*pj%p%rlon0
    pj%p%rlat0 = raddeg_r8*pj%p%rlat0 - cj*pj%p%dlat
    pj%p%rlon0 = raddeg_r8*pj%p%rlon0 - ci*pj%p%dlon
    if ( pj%p%rlon0 >  deg180 ) pj%p%rlon0 = pj%p%rlon0 - deg360
    if ( pj%p%rlon0 < -deg180 ) pj%p%rlon0 = pj%p%rlon0 + deg360
    if ( luvrot ) then
      if ( .not. pj%p%skiprot ) then
        call getmem(pj%f1,1,pj%p%nlon,1,pj%p%nlat,'projections:f1')
        call getmem(pj%f2,1,pj%p%nlon,1,pj%p%nlat,'projections:f2')
        do j = 1, pj%p%nlat
          do i = 1, pj%p%nlon
            ri = i
            rj = j
            call ijll_rl(pj,ri,rj,lat,lon)
            lam = degrad_r8*lon
            phi = degrad_r8*lat
            dlam = plam - lam
            arg1 = cos(pphi)*sin(dlam)
            arg2 = cos(phi)*sin(pphi) - sin(phi)*cos(pphi)*cos(dlam)
            delta = atan2(arg1,arg2)
            pj%f1(i,j) = cos(delta)
            pj%f2(i,j) = sin(delta)
          end do
        end do
      end if
    end if
  end subroutine setup_rll

  pure subroutine ijll_rl(pj,i,j,lat,lon)
    class(regcm_projection), intent(in) :: pj
    real(rkx), intent(in) :: i, j
    real(rkx), intent(out) :: lat, lon
    real(rk8) :: phi, lam, glat, glon
    phi = pj%p%rlat0 + (j-1.0_rk8) * pj%p%dlat
    lam = pj%p%rlon0 + (i-1.0_rk8) * pj%p%dlon
    ! Across the rotated pole: reflect the latitude, opposite meridian
    if ( phi >  deg90 ) then
      phi =  deg180 - phi
      lam = lam + deg180
    end if
    if ( phi < -deg90 ) then
      phi = -deg180 - phi
      lam = lam + deg180
    end if
    phi = degrad_r8 * phi
    lam = degrad_r8 * lam
    call rotated_to_geo(pj%p%phi0,pj%p%lam0,phi,lam,glat,glon)
    glon = raddeg_r8*glon
    if ( glon >  deg180 ) glon = glon - deg360
    if ( glon < -deg180 ) glon = glon + deg360
    lat = real(raddeg_r8*glat,rkx)
    lon = real(glon,rkx)
  end subroutine ijll_rl

  pure subroutine llij_rl(pj,lat,lon,i,j)
    implicit none
    class(regcm_projection), intent(in) :: pj
    real(rkx), intent(in) :: lat, lon
    real(rkx), intent(out) :: i, j
    real(rk8) :: rlat, rlon, phi, lam
    phi = lat
    lam = lon
    if ( phi >  deg90 ) then
      phi =  deg180 - phi
      lam = lam + deg180
    end if
    if ( phi < -deg90 ) then
      phi = -deg180 - phi
      lam = lam + deg180
    end if
    phi = degrad_r8 * phi
    lam = degrad_r8 * lam
    call geo_to_rotated(pj%p%phi0,pj%p%lam0,phi,lam,rlat,rlon)
    rlon = raddeg_r8*rlon
    rlat = raddeg_r8*rlat
    ! The domain may cross the rotated longitude 180: take the rotated
    ! longitude within 180 degrees of the domain centre
    rlon = pj%p%xoff + modulo(rlon - pj%p%xoff + deg180, deg360) - deg180
    i = real(( rlon - pj%p%rlon0 ) / pj%p%dlon + 1.0_rk8,rkx)
    j = real(( rlat - pj%p%rlat0 ) / pj%p%dlat + 1.0_rk8,rkx)
  end subroutine llij_rl

  !
  ! Rotation between geographic and rotated coordinates (radians). The
  ! rotated frame is the geographic one rotated by phi0 around the axis
  ! through the longitudes lam0 +/- 90. With phi0 = 0 (pole at the
  ! geographic pole) the coordinates are the same.
  !
  pure subroutine geo_to_rotated(phi0,lam0,phi,lam,rlat,rlon)
    implicit none
    real(rk8), intent(in) :: phi0, lam0, phi, lam
    real(rk8), intent(out) :: rlat, rlon
    real(rk8) :: x, y, z
    if ( abs(phi0) <= 1.0e-7_rk8 ) then
      rlat = phi
      rlon = lam
      return
    end if
    x = cos(phi0)*cos(phi)*cos(lam-lam0) + sin(phi0)*sin(phi)
    y = cos(phi)*sin(lam-lam0)
    z = -sin(phi0)*cos(phi)*cos(lam-lam0) + cos(phi0)*sin(phi)
    rlat = asin(max(-1.0_rk8,min(1.0_rk8,z)))
    if ( abs(x) < dlowval_r8 .and. abs(y) < dlowval_r8 ) then
      rlon = 0.0_rk8
    else
      rlon = atan2(y,x)
    end if
  end subroutine geo_to_rotated

  pure subroutine rotated_to_geo(phi0,lam0,rlat,rlon,phi,lam)
    implicit none
    real(rk8), intent(in) :: phi0, lam0, rlat, rlon
    real(rk8), intent(out) :: phi, lam
    real(rk8) :: x, y, z
    if ( abs(phi0) <= 1.0e-7_rk8 ) then
      phi = rlat
      lam = rlon
      return
    end if
    x = cos(phi0)*cos(rlat)*cos(rlon) - sin(phi0)*sin(rlat)
    y = cos(rlat)*sin(rlon)
    z = sin(phi0)*cos(rlat)*cos(rlon) + cos(phi0)*sin(rlat)
    phi = asin(max(-1.0_rk8,min(1.0_rk8,z)))
    if ( abs(x) < dlowval_r8 .and. abs(y) < dlowval_r8 ) then
      lam = lam0
    else
      lam = lam0 + atan2(y,x)
    end if
  end subroutine rotated_to_geo

  subroutine setup_lcc(pj,clon,clat,ci,cj,ds,slon,trlat1,trlat2,luvrot)
    implicit none
    type(regcm_projection), intent(inout) :: pj
    real(rk8), intent(in) :: ci, cj, slon, clat, clon, ds, &
                              trlat1, trlat2
    logical, intent(in) :: luvrot
    real(rk8) :: arg, deltalon1, tl1r, tl2r, pf3
    real(rkx) :: ri, rj, lat, lon
    integer(ik4) :: i, j

    pj%p%stdlon = slon
    pj%p%stdlat = clat
    pj%p%truelat1 = trlat1
    pj%p%truelat2 = trlat2
    tl1r = pj%p%truelat1*degrad_r8
    tl2r = pj%p%truelat2*degrad_r8
    pj%p%colat1  = degrad_r8*(deg90 - pj%p%truelat1)
    pj%p%colat2  = degrad_r8*(deg90 - pj%p%truelat2)
    pj%p%nfac = (log(sin(pj%p%colat1))         - log(sin(pj%p%colat2))) / &
                (log(tan(pj%p%colat1*0.5_rk8)) - log(tan(pj%p%colat2*0.5_rk8)))
    if ( pj%p%stdlat > 0.0_rk8 ) then
      pj%p%hemi =  1.0_rk8
    else
      pj%p%hemi = -1.0_rk8
    end if
    pj%p%rebydx = earthrad_r8 / ds
    if ( abs(pj%p%truelat1-pj%p%truelat2) > 0.1_rk8 ) then
      pj%p%conefac = log10(cos(tl1r)) - log10(cos(tl2r))
      pj%p%conefac = pj%p%conefac / &
                (log10(tan((deg45-abs(pj%p%truelat1)/2.0_rk8)*degrad_r8)) - &
                 log10(tan((deg45-abs(pj%p%truelat2)/2.0_rk8)*degrad_r8)))
      pj%p%lamtan = .false.
    else
      pj%p%conefac = sin(abs(tl1r))
      pj%p%lamtan = .true.
    end if
    pj%p%rconefac = 1.0_rk8/pj%p%conefac
    deltalon1 = clon - pj%p%stdlon
    if ( deltalon1 >  deg180 ) deltalon1 = deltalon1 - deg360
    if ( deltalon1 < -deg180 ) deltalon1 = deltalon1 + deg360
    pj%p%ctl1r = cos(tl1r)
    pj%p%xct1 = tan((deg90*pj%p%hemi-pj%p%truelat1)*degrad_r8*0.5_rk8)
    pj%p%rsw = pj%p%rebydx * pj%p%ctl1r*pj%p%rconefac * &
        (tan((deg90*pj%p%hemi-clat)*degrad_r8*0.5_rk8)/pj%p%xct1)**pj%p%conefac
    arg = pj%p%conefac*(deltalon1*degrad_r8)
    pj%p%polei = pj%p%hemi*ci - pj%p%hemi * pj%p%rsw * sin(arg)
    pj%p%polej = pj%p%hemi*cj + pj%p%rsw * cos(arg)
    pj%p%chi1 = (deg90 - pj%p%hemi*pj%p%truelat1)*degrad_r8
    pj%p%chi2 = (deg90 - pj%p%hemi*pj%p%truelat2)*degrad_r8
    pj%p%tanchi1h = tan(pj%p%chi1*0.5_rk8)
    pj%p%tchi1 = tan(pj%p%chi1)
    pj%p%schi1 = sin(pj%p%chi1)
    if ( luvrot ) then
      call getmem(pj%f1,1,pj%p%nlon,1,pj%p%nlat,'projections:f1')
      call getmem(pj%f2,1,pj%p%nlon,1,pj%p%nlat,'projections:f2')
      call getmem(pj%f3,1,pj%p%nlon,1,pj%p%nlat,'projections:f3')
      do j = 1, pj%p%nlat
        do i = 1, pj%p%nlon
          ri = i
          rj = j
          call ijll_lc(pj,ri,rj,lat,lon)
          pf3 = uvrot_lc(pj,lon)
          pj%f3(i,j) = pf3
          pj%f1(i,j) = cos(pf3)
          pj%f2(i,j) = sin(pf3)
        end do
      end do
    end if
  end subroutine setup_lcc

  pure subroutine ijll_lc(pj,i,j,lat,lon)
    implicit none
    class(regcm_projection), intent(in) :: pj
    real(rkx), intent(in) :: i, j
    real(rkx), intent(out) :: lat, lon
    real(rk8) :: chi
    real(rk8) :: inew, jnew, xx, yy, r2, r

    inew = pj%p%hemi * i
    jnew = pj%p%hemi * j
    xx = inew - pj%p%polei
    yy = pj%p%polej - jnew
    r2 = (xx*xx + yy*yy)
    r = sqrt(r2)/pj%p%rebydx
    if ( abs(r2) < dlowval_r8 ) then
      lat = real(pj%p%hemi*deg90,rkx)
      lon = real(pj%p%stdlon,rkx)
    else
      lon = real(pj%p%stdlon + &
           raddeg_r8 * atan2(pj%p%hemi*xx,yy)*pj%p%rconefac,rkx)
      if ( pj%p%lamtan ) then
        chi = 2.0_rk8 * &
          atan(((r/pj%p%tchi1)**pj%p%rconefac)*pj%p%tanchi1h)
      else
        chi = 2.0_rk8 * &
          atan(((r*pj%p%conefac/pj%p%schi1)**pj%p%rconefac)*pj%p%tanchi1h)
      end if
      lat = real((deg90-chi*raddeg_r8)*pj%p%hemi,rkx)
    end if
    if ( lon >  180.0_rkx ) lon = lon - 360.0_rkx
    if ( lon < -180.0_rkx ) lon = lon + 360.0_rkx
  end subroutine ijll_lc

  pure subroutine llij_lc(pj,lat,lon,i,j)
    implicit none
    class(regcm_projection), intent(in) :: pj
    real(rkx), intent(in) :: lat, lon
    real(rkx), intent(out) :: i, j
    real(rk8) :: arg, deltalon, rm

    deltalon = lon - pj%p%stdlon
    if ( deltalon > +deg180 ) deltalon = deltalon - deg360
    if ( deltalon < -deg180 ) deltalon = deltalon + deg360
    rm = pj%p%rebydx * pj%p%ctl1r * pj%p%rconefac * &
        (tan((deg90*pj%p%hemi-lat)*degrad_r8*0.5_rk8)/pj%p%xct1)**pj%p%conefac
    arg = pj%p%conefac*(deltalon*degrad_r8)
    i = real(pj%p%hemi*(pj%p%polei + pj%p%hemi * rm * sin(arg)),rkx)
    j = real(pj%p%hemi*(pj%p%polej - rm * cos(arg)),rkx)
  end subroutine llij_lc

  subroutine setup_plr(pj,clon,clat,ci,cj,ds,slon,luvrot)
    implicit none
    type(regcm_projection), intent(inout) :: pj
    real(rk8), intent(in) :: clat, clon, cj, ci, ds, slon
    logical, intent(in) :: luvrot
    real(rk8) :: ala1, alo1, pf3
    real(rkx) :: lat, lon, ri, rj
    integer(ik4) :: i, j

    pj%p%stdlon = slon
    pj%p%stdlat = clat
    if ( pj%p%stdlat > 0.0_rk8 ) then
      pj%p%hemi = 1.0_rk8
    else
      pj%p%hemi = -1.0_rk8
    end if
    pj%p%rebydx = earthrad_r8 / ds
    pj%p%reflon = pj%p%stdlon + deg90
    ala1 = clat*degrad_r8
    alo1 = (clon-pj%p%reflon)*degrad_r8
    pj%p%scale_top = 1.0_rk8 + pj%p%hemi * sin(ala1)
    pj%p%rsw = pj%p%rebydx * &
         cos(ala1)*pj%p%scale_top/(1.0_rk8+pj%p%hemi*sin(ala1))
    pj%p%polei = ci - pj%p%rsw * cos(alo1)
    pj%p%polej = cj - pj%p%hemi * pj%p%rsw * sin(alo1)
    if ( luvrot ) then
      call getmem(pj%f1,1,pj%p%nlon,1,pj%p%nlat,'projections:f1')
      call getmem(pj%f2,1,pj%p%nlon,1,pj%p%nlat,'projections:f2')
      call getmem(pj%f3,1,pj%p%nlon,1,pj%p%nlat,'projections:f3')
      do j = 1, pj%p%nlat
        do i = 1, pj%p%nlon
          ri = i
          rj = j
          call ijll_ps(pj,ri,rj,lat,lon)
          pf3 = uvrot_ps(pj,lon)
          pj%f3(i,j) = pf3
          pj%f1(i,j) = cos(pf3)
          pj%f2(i,j) = sin(pf3)
        end do
      end do
    end if
  end subroutine setup_plr

  pure subroutine llij_ps(pj,lat,lon,i,j)
    implicit none
    class(regcm_projection), intent(in) :: pj
    real(rkx), intent(in) :: lat, lon
    real(rkx), intent(out) :: i, j
    real(rk8) :: ala, alo, rm, deltalon

    deltalon = lon - pj%p%reflon
    if ( deltalon > +deg180 ) deltalon = deltalon - deg360
    if ( deltalon < -deg180 ) deltalon = deltalon + deg360
    alo = deltalon * degrad_r8
    ala = lat * degrad_r8
    rm = pj%p%rebydx * cos(ala) * &
         pj%p%scale_top/(1.0_rk8 + pj%p%hemi * sin(ala))
    i = real(pj%p%polei + rm * cos(alo),rkx)
    j = real(pj%p%polej + pj%p%hemi * rm * sin(alo),rkx)
  end subroutine llij_ps

  pure subroutine ijll_ps(pj,i,j,lat,lon)
    implicit none
    class(regcm_projection), intent(in) :: pj
    real(rkx), intent(in) :: i, j
    real(rkx), intent(out) :: lat, lon
    real(rk8) :: xx, yy, r2, gi2, arcc

    xx = i - pj%p%polei
    yy = (j - pj%p%polej) * pj%p%hemi
    r2 = xx*xx + yy*yy
    if ( abs(r2) < dlowval_r8 ) then
      lat = real(pj%p%hemi*deg90,rkx)
      lon = real(pj%p%reflon,rkx)
    else
      gi2 = (pj%p%rebydx * pj%p%scale_top)**2
      lat = real(raddeg_r8 * pj%p%hemi * asin((gi2-r2)/(gi2+r2)),rkx)
      arcc = acos(xx/sqrt(r2))
      if ( yy > 0.0_rk8 ) then
        lon = real(pj%p%reflon + raddeg_r8 * arcc,rkx)
      else
        lon = real(pj%p%reflon - raddeg_r8 * arcc,rkx)
      end if
    end if
    if ( lon >  180.0_rkx ) lon = lon - 360.0_rkx
    if ( lon < -180.0_rkx ) lon = lon + 360.0_rkx
  end subroutine ijll_ps

  subroutine setup_mrc(pj,clon,clat,ci,cj,ds)
    implicit none
    type(regcm_projection), intent(inout) :: pj
    real(rk8), intent(in) :: clat, clon, cj, ci, ds
    real(rk8) :: clain

    pj%p%stdlon = clon
    pj%p%stdlat = clat
    clain = cos(clat*degrad_r8)
    pj%p%dlon = ds / (earthrad_r8 * clain)
    pj%p%rsw = 0.0_rk8
    if ( abs(clat) > dlowval_r8 ) then
      pj%p%rsw = (log(tan(0.5_rk8*((clat+deg90)*degrad_r8))))/pj%p%dlon
    end if
    pj%p%polei = ci
    pj%p%polej = cj
  end subroutine setup_mrc

  pure subroutine llij_mc(pj,lat,lon,i,j)
    implicit none
    class(regcm_projection), intent(in) :: pj
    real(rkx), intent(in) :: lat, lon
    real(rkx), intent(out) :: i, j
    real(rk8) :: deltalon

    deltalon = lon - pj%p%stdlon
    if ( deltalon > +deg180 ) deltalon = deltalon - deg360
    if ( deltalon < -deg180 ) deltalon = deltalon + deg360
    i = real(pj%p%polei + (deltalon/(pj%p%dlon*raddeg_r8)),rkx)
    j = real(pj%p%polej + &
         (log(tan(0.5_rk8*((lat+deg90)*degrad_r8)))) / pj%p%dlon - pj%p%rsw,rkx)
  end subroutine llij_mc

  pure subroutine ijll_mc(pj,i,j,lat,lon)
    implicit none
    class(regcm_projection), intent(in) :: pj
    real(rkx), intent(in) :: i, j
    real(rkx), intent(out) :: lat, lon

    lat = real(2.0_rk8 * &
        atan(exp(pj%p%dlon*(pj%p%rsw + j-pj%p%polej)))*raddeg_r8-deg90,rkx)
    lon = real((i-pj%p%polei)*pj%p%dlon*raddeg_r8 + pj%p%stdlon,rkx)
    if ( lon >  180.0_rkx ) lon = lon - 360.0_rkx
    if ( lon < -180.0_rkx ) lon = lon + 360.0_rkx
  end subroutine ijll_mc

  subroutine setup_rmc(pj,clon,clat,ci,cj,ds,plon,plat,luvrot)
    implicit none
    type(regcm_projection), intent(inout) :: pj
    real(rk8), intent(in) :: clat, clon, cj, ci, ds, plon, plat
    logical, intent(in) :: luvrot
    real(rk8) :: plam, pphi, zphipol, pf3
    real(rkx) :: lat, lon, ri, rj
    integer(ik4) :: i, j
    pj%p%dlon = ds*raddeg_r8/earthrad_r8
    pj%p%dlat = ds*raddeg_r8/earthrad_r8
    pj%p%polei = ci
    pj%p%polej = cj
    pphi = deg90 - plat
    plam = plon + deg180
    if ( plam > deg180 ) plam = plam - deg360
    pj%p%zlampol = degrad_r8*plam
    zphipol = degrad_r8*pphi
    pj%p%zsinpol = sin(zphipol)
    pj%p%zcospol = cos(zphipol)
    ! Rotated Mercator coordinates of the domain centre
    pj%p%xoff = 0.0_rk8
    pj%p%yoff = 0.0_rk8
    lat = real(clat,rkx)
    lon = real(clon,rkx)
    call llij_rc(pj,lat,lon,ri,rj)
    pj%p%xoff = (ri-ci)*pj%p%dlon
    pj%p%yoff = (rj-cj)*pj%p%dlat
    if ( luvrot ) then
      if ( plat > 0.0_rk8 ) then
        pphi = (deg90 - plat) * degrad_r8
        pj%p%pollam = plon + deg180
      else
        pphi = (deg90 + plat) * degrad_r8
        pj%p%pollam = plon
      end if
      if ( pj%p%pollam > deg180 ) then
        pj%p%pollam = pj%p%pollam - deg360
      end if
      if ( pj%p%pollam < -deg180 ) then
        pj%p%pollam = deg360 + pj%p%pollam
      end if
      pj%p%polcphi = cos(pphi)
      pj%p%polsphi = sin(pphi)
      call getmem(pj%f1,1,pj%p%nlon,1,pj%p%nlat,'projections:f1')
      call getmem(pj%f2,1,pj%p%nlon,1,pj%p%nlat,'projections:f2')
      call getmem(pj%f3,1,pj%p%nlon,1,pj%p%nlat,'projections:f3')
      do j = 1, pj%p%nlat
        rj = j
        do i = 1, pj%p%nlon
          ri = i
          call ijll_rc(pj,ri,rj,lat,lon)
          call uvrot_rc(pj,lat,lon,pj%f1(i,j),pj%f2(i,j))
          pf3 = acos(pj%f1(i,j))
          pj%f3(i,j) = pf3
        end do
      end do
    end if
  end subroutine setup_rmc

  pure subroutine llij_rc(pj,lat,lon,i,j)
    implicit none
    class(regcm_projection), intent(in) :: pj
    real(rkx), intent(in) :: lat, lon
    real(rkx), intent(out) :: i, j
    real(rk8) :: zarg, zarg1, zarg2, zlam, zphi
    real(rk8) :: lams, phis

    zphi = degrad_r8*lat
    zlam = lon
    if ( zlam > deg180 ) zlam = zlam - deg360
    zlam = degrad_r8*zlam
    zarg = pj%p%zcospol*cos(zphi)*cos(zlam-pj%p%zlampol) + &
           pj%p%zsinpol*sin(zphi)
    phis = asin(zarg)
    phis = log(tan(phis*0.5_rk8+atan(1.0_rk8)))*raddeg_r8
    zarg1 = -sin(zlam-pj%p%zlampol)*cos(zphi)
    zarg2 = -pj%p%zsinpol*cos(zphi)*cos(zlam-pj%p%zlampol) + &
            pj%p%zcospol*sin(zphi)
    if ( abs(zarg2) >= dlowval_r8 ) then
      lams = raddeg_r8*atan2(zarg1,zarg2)
    else if ( abs(zarg1) < dlowval_r8 ) then
      lams = deg00
    else if ( zarg1 > 0.0_rk8 ) then
      lams = deg90
    else
      lams = -deg90
    end if
    i = real(pj%p%polei + (lams-pj%p%xoff)/pj%p%dlon,rkx)
    j = real(pj%p%polej + (phis-pj%p%yoff)/pj%p%dlat,rkx)
  end subroutine llij_rc

  pure subroutine ijll_rc(pj,i,j,lat,lon)
    implicit none
    class(regcm_projection), intent(in) :: pj
    real(rkx), intent(in) :: i, j
    real(rkx), intent(out) :: lat, lon
    real(rk8) :: xr, yr, arg, zarg1, zarg2

    xr = pj%p%xoff + (i-pj%p%polei)*pj%p%dlon
    if ( xr > deg180 ) xr = xr - deg360
    xr = degrad_r8*xr
    yr = pj%p%yoff + (j-pj%p%polej)*pj%p%dlat
    yr = 2.0_rk8*atan(exp(degrad_r8*yr)) - atan(1.0_rk8)*2.0_rk8
    arg = pj%p%zcospol*cos(yr)*cos(xr) + pj%p%zsinpol*sin(yr)
    lat = real(raddeg_r8*asin(arg),rkx)
    zarg1 = sin(pj%p%zlampol)*(-pj%p%zsinpol*cos(xr)*cos(yr)+ &
            pj%p%zcospol*sin(yr))-cos(pj%p%zlampol)*sin(xr)*cos(yr)
    zarg2 = cos(pj%p%zlampol)*(-pj%p%zsinpol*cos(xr)*cos(yr)+ &
            pj%p%zcospol*sin(yr))+sin(pj%p%zlampol)*sin(xr)*cos(yr)
    if ( abs(zarg2) >= dlowval_r8 ) then
      lon = real(raddeg_r8*atan2(zarg1,zarg2),rkx)
      if ( lon >  180.0_rkx ) lon = lon - 360.0_rkx
      if ( lon < -180.0_rkx ) lon = lon + 360.0_rkx
    else if ( abs(zarg1) < dlowval_r8 ) then
      lon = deg00
    else if ( zarg1 > 0.0_rk8 ) then
      lon = deg90
    else
      lon = -deg90
    end if
  end subroutine ijll_rc

  real(rkx) function rounder(xval,ltop)
    implicit none
    real(rkx), intent(in) :: xval
    logical, intent(in) :: ltop
    integer(ik4) :: tmpval
    if ( ltop ) then
      tmpval = ceiling(xval*1000.0_rkx)
    else
      tmpval = floor(xval*1000.0_rkx)
    end if
    rounder = tmpval/1000.0_rkx
  end function rounder

  ! Arguments in radiants

  subroutine get_equator(plat,plon,elat,elon)
    implicit none
    real(rk8), intent(in) :: plat, plon
    real(rk8), intent(out) :: elat, elon

    if ( plon < 0.0_rk8 ) then
      elon = plon + mathpi_r8
    else if ( plon > 0.0_rk8 ) then
      elon = plon - mathpi_r8
    else
      elon = 0.0_rk8
    end if
    if ( plat > 0.0_rk8 ) then
      elat = halfpi_r8-plat
    else if ( plat < 0.0_rk8 ) then
      elat = -(halfpi_r8+plat)
    else
      elat = 0.0_rk8
    end if
  end subroutine get_equator

  pure elemental real(rkx) function fac_rl(pj,xlat,xlon) result(xmap)
    implicit none
    type(regcm_projection), intent(in) :: pj
    real(rkx), intent(in) :: xlat, xlon
    real(rkx) :: ri, rj
    call llij_rl(pj,xlat,xlon,ri,rj)
    xmap = real(1.0_rk8/cos(degrad_r8*(pj%p%rlat0+(rj-1.0_rk8)*pj%p%dlat)),rkx)
  end function fac_rl

  pure elemental real(rkx) function fac_lc(pj,lat) result(xmap)
    implicit none
    type(regcm_projection), intent(in) :: pj
    real(rkx), intent(in) :: lat
    real(rk8) :: colat
    colat = degrad_r8*(deg90-lat)
    if ( .not. pj%p%lamtan ) then
      xmap = real(sin(pj%p%colat2)/sin(colat) * &
          (tan(colat*0.5_rk8)/tan(pj%p%colat2*0.5_rk8))**pj%p%nfac,rkx)
    else
      xmap = real(sin(pj%p%colat1)/sin(colat) * &
          (tan(colat*0.5_rk8)/tan(pj%p%colat1*0.5_rk8))**cos(pj%p%colat1),rkx)
    endif
  end function fac_lc

  pure elemental real(rkx) function fac_ps(pj,lat) result(xmap)
    implicit none
    type(regcm_projection), intent(in) :: pj
    real(rkx), intent(in) :: lat
    xmap = real(pj%p%scale_top/(1.0_rk8 + pj%p%hemi * sin(lat*degrad_r8)),rkx)
  end function fac_ps

  pure elemental real(rkx) function fac_mc(pj,lat) result(xmap)
    implicit none
    type(regcm_projection), intent(in) :: pj
    real(rkx), intent(in) :: lat
    xmap = real(cos(pj%p%stdlat*degrad_r8)/cos(lat*degrad_r8),rkx)
  end function fac_mc

  pure elemental real(rkx) function fac_rc(pj,xlat,xlon) result(xmap)
    implicit none
    type(regcm_projection), intent(in) :: pj
    real(rkx), intent(in) :: xlat, xlon
    real(rkx) :: ri, rj
    real(rk8) :: yr
    call llij_rc(pj,xlat,xlon,ri,rj)
    yr = pj%p%yoff + (rj-pj%p%polej)*pj%p%dlat
    yr = 2.0_rk8*atan(exp(degrad_r8*yr)) - halfpi_r8
    xmap = real(1.0_rk8/cos(yr),rkx)
  end function fac_rc

  pure elemental real(rk8) function uvrot_lc(pj,lon) result(alpha)
    implicit none
    type(regcm_projection), intent(in) :: pj
    real(rkx), intent(in) :: lon
    real(rk8) :: deltalon
    deltalon = pj%p%stdlon - lon
    if ( deltalon > +deg180 ) deltalon = deltalon - deg360
    if ( deltalon < -deg180 ) deltalon = deltalon + deg360
    alpha = deltalon*degrad_r8*pj%p%conefac
  end function uvrot_lc

  pure elemental real(rk8) function uvrot_ps(pj,lon) result(alpha)
    implicit none
    type(regcm_projection), intent(in) :: pj
    real(rkx), intent(in) :: lon
    real(rk8) :: deltalon
    deltalon = pj%p%stdlon - lon
    if ( deltalon > +deg180 ) deltalon = deltalon - deg360
    if ( deltalon < -deg180 ) deltalon = deltalon + deg360
    alpha = deltalon*degrad_r8
  end function uvrot_ps

  pure subroutine uvrot_rc(pj,lat,lon,cosdel,sindel)
    implicit none
    type(regcm_projection), intent(in) :: pj
    real(rkx), intent(in) :: lon, lat
    real(rk8), intent(out) :: cosdel, sindel
    real(rk8) :: zphi, zrla, zrlap, zarg1, zarg2, znorm
    zrla = lon
    if ( lat > deg90-0.00001_rk8 ) zrla = 0.0_rk8
    zrlap = (pj%p%pollam - zrla)*degrad_r8
    zphi = lat*degrad_r8
    zarg1 = pj%p%polcphi*sin(zrlap)
    zarg2 = pj%p%polsphi*cos(zphi) - pj%p%polcphi*sin(zphi)*cos(zrlap)
    znorm = 1.0_rk8/sqrt(zarg1*zarg1 + zarg2*zarg2)
    sindel = zarg1*znorm
    cosdel = zarg2*znorm
  end subroutine uvrot_rc

  subroutine rotate2_lc(pj,u,v)
    implicit none
    class(regcm_projection), intent(in) :: pj
    real(rkx), pointer, contiguous, dimension(:,:), intent(inout) :: u, v
    integer(ik4) :: i1, i2, j1, j2, i, j
    real(rk8) :: tmp
    i1 = lbound(u,1)
    i2 = ubound(u,1)
    j1 = lbound(u,2)
    j2 = ubound(u,2)
    if ( pj%p%stdlat >= 0.0_rk8 ) then
      do j = j1, j2
        do i = i1, i2
          tmp = u(i,j)*pj%f1(i,j) + v(i,j)*pj%f2(i,j)
          v(i,j) = real(-u(i,j)*pj%f2(i,j) + v(i,j)*pj%f1(i,j),rkx)
          u(i,j) = real(tmp,rkx)
        end do
      end do
    else
      do j = j1, j2
        do i = i1, i2
          tmp = u(i,j)*pj%f1(i,j) - v(i,j)*pj%f2(i,j)
          v(i,j) = real(v(i,j)*pj%f1(i,j) + u(i,j)*pj%f2(i,j),rkx)
          u(i,j) = real(tmp,rkx)
        end do
      end do
    end if
  end subroutine rotate2_lc

  subroutine backrotate2_lc(pj,u,v)
    implicit none
    class(regcm_projection), intent(in) :: pj
    real(rkx), pointer, contiguous, dimension(:,:), intent(inout) :: u, v
    integer(ik4) :: i1, i2, j1, j2, i, j
    real(rk8) :: tmp
    i1 = lbound(u,1)
    i2 = ubound(u,1)
    j1 = lbound(u,2)
    j2 = ubound(u,2)
    if ( pj%p%stdlat >= 0.0_rk8 ) then
      do j = j1, j2
        do i = i1, i2
          tmp = u(i,j)*pj%f1(i,j) - v(i,j)*pj%f2(i,j)
          v(i,j) = real(u(i,j)*pj%f2(i,j) + v(i,j)*pj%f1(i,j),rkx)
          u(i,j) = real(tmp,rkx)
        end do
      end do
    else
      do j = j1, j2
        do i = i1, i2
          tmp = u(i,j)*pj%f1(i,j) + v(i,j)*pj%f2(i,j)
          v(i,j) = real(v(i,j)*pj%f1(i,j) - u(i,j)*pj%f2(i,j),rkx)
          u(i,j) = real(tmp,rkx)
        end do
      end do
    end if
  end subroutine backrotate2_lc

  subroutine rotate3_lc(pj,u,v)
    implicit none
    class(regcm_projection), intent(in) :: pj
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: u, v
    integer(ik4) :: i1, i2, j1, j2, k1, k2, i, j, k
    real(rk8) :: tmp
    i1 = lbound(u,1)
    i2 = ubound(u,1)
    j1 = lbound(u,2)
    j2 = ubound(u,2)
    k1 = lbound(u,3)
    k2 = ubound(u,3)
    if ( pj%p%stdlat >= 0.0_rk8 ) then
      do k = k1, k2
        do j = j1, j2
          do i = i1, i2
            tmp = u(i,j,k)*pj%f1(i,j) + v(i,j,k)*pj%f2(i,j)
            v(i,j,k) = real(-u(i,j,k)*pj%f2(i,j) + v(i,j,k)*pj%f1(i,j),rkx)
            u(i,j,k) = real(tmp,rkx)
          end do
        end do
      end do
    else
      do k = k1, k2
        do j = j1, j2
          do i = i1, i2
            tmp = u(i,j,k)*pj%f1(i,j) - v(i,j,k)*pj%f2(i,j)
            v(i,j,k) = real(v(i,j,k)*pj%f1(i,j) + u(i,j,k)*pj%f2(i,j),rkx)
            u(i,j,k) = real(tmp,rkx)
          end do
        end do
      end do
    end if
  end subroutine rotate3_lc

  subroutine backrotate3_lc(pj,u,v)
    implicit none
    class(regcm_projection), intent(in) :: pj
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: u, v
    integer(ik4) :: i1, i2, j1, j2, k1, k2, i, j, k
    real(rk8) :: tmp
    i1 = lbound(u,1)
    i2 = ubound(u,1)
    j1 = lbound(u,2)
    j2 = ubound(u,2)
    k1 = lbound(u,3)
    k2 = ubound(u,3)
    if ( pj%p%stdlat >= 0.0_rk8 ) then
      do k = k1, k2
        do j = j1, j2
          do i = i1, i2
            tmp = u(i,j,k)*pj%f1(i,j) - v(i,j,k)*pj%f2(i,j)
            v(i,j,k) = real(u(i,j,k)*pj%f2(i,j) + v(i,j,k)*pj%f1(i,j),rkx)
            u(i,j,k) = real(tmp,rkx)
          end do
        end do
      end do
    else
      do k = k1, k2
        do j = j1, j2
          do i = i1, i2
            tmp = u(i,j,k)*pj%f1(i,j) + v(i,j,k)*pj%f2(i,j)
            v(i,j,k) = real(v(i,j,k)*pj%f1(i,j) - u(i,j,k)*pj%f2(i,j),rkx)
            u(i,j,k) = real(tmp,rkx)
          end do
        end do
      end do
    end if
  end subroutine backrotate3_lc

  subroutine rotate2_ps(pj,u,v)
    implicit none
    class(regcm_projection), intent(in) :: pj
    real(rkx), pointer, contiguous, dimension(:,:), intent(inout) :: u, v
    integer(ik4) :: i1, i2, j1, j2, i, j
    real(rk8) :: tmp
    i1 = lbound(u,1)
    i2 = ubound(u,1)
    j1 = lbound(u,2)
    j2 = ubound(u,2)
    if ( pj%p%stdlat >= 0.0_rk8 ) then
      do j = j1, j2
        do i = i1, i2
          tmp = u(i,j)*pj%f1(i,j) + v(i,j)*pj%f2(i,j)
          v(i,j) = real(-u(i,j)*pj%f2(i,j) + v(i,j)*pj%f1(i,j),rkx)
          u(i,j) = real(tmp,rkx)
        end do
      end do
    else
      do j = j1, j2
        do i = i1, i2
          tmp = u(i,j)*pj%f1(i,j) - v(i,j)*pj%f2(i,j)
          v(i,j) = real(v(i,j)*pj%f1(i,j) + u(i,j)*pj%f2(i,j),rkx)
          u(i,j) = real(tmp,rkx)
        end do
      end do
    end if
  end subroutine rotate2_ps

  subroutine backrotate2_ps(pj,u,v)
    implicit none
    class(regcm_projection), intent(in) :: pj
    real(rkx), pointer, contiguous, dimension(:,:), intent(inout) :: u, v
    integer(ik4) :: i1, i2, j1, j2, i, j
    real(rk8) :: tmp
    i1 = lbound(u,1)
    i2 = ubound(u,1)
    j1 = lbound(u,2)
    j2 = ubound(u,2)
    if ( pj%p%stdlat >= 0.0_rk8 ) then
      do j = j1, j2
        do i = i1, i2
          tmp = u(i,j)*pj%f1(i,j) - v(i,j)*pj%f2(i,j)
          v(i,j) = real(u(i,j)*pj%f2(i,j) + v(i,j)*pj%f1(i,j),rkx)
          u(i,j) = real(tmp,rkx)
        end do
      end do
    else
      do j = j1, j2
        do i = i1, i2
          tmp = u(i,j)*pj%f1(i,j) + v(i,j)*pj%f2(i,j)
          v(i,j) = real(v(i,j)*pj%f1(i,j) - u(i,j)*pj%f2(i,j),rkx)
          u(i,j) = real(tmp,rkx)
        end do
      end do
    end if
  end subroutine backrotate2_ps

  subroutine rotate3_ps(pj,u,v)
    implicit none
    class(regcm_projection), intent(in) :: pj
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: u, v
    integer(ik4) :: i1, i2, j1, j2, k1, k2, i, j, k
    real(rk8) :: tmp
    i1 = lbound(u,1)
    i2 = ubound(u,1)
    j1 = lbound(u,2)
    j2 = ubound(u,2)
    k1 = lbound(u,3)
    k2 = ubound(u,3)
    if ( pj%p%stdlat >= 0.0_rk8 ) then
      do k = k1, k2
        do j = j1, j2
          do i = i1, i2
            tmp = u(i,j,k)*pj%f1(i,j) + v(i,j,k)*pj%f2(i,j)
            v(i,j,k) = real(-u(i,j,k)*pj%f2(i,j) + v(i,j,k)*pj%f1(i,j),rkx)
            u(i,j,k) = real(tmp,rkx)
          end do
        end do
      end do
    else
      do k = k1, k2
        do j = j1, j2
          do i = i1, i2
            tmp = u(i,j,k)*pj%f1(i,j) - v(i,j,k)*pj%f2(i,j)
            v(i,j,k) = real(v(i,j,k)*pj%f1(i,j) + u(i,j,k)*pj%f2(i,j),rkx)
            u(i,j,k) = real(tmp,rkx)
          end do
        end do
      end do
    end if
  end subroutine rotate3_ps

  subroutine backrotate3_ps(pj,u,v)
    implicit none
    class(regcm_projection), intent(in) :: pj
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: u, v
    integer(ik4) :: i1, i2, j1, j2, k1, k2, i, j, k
    real(rk8) :: tmp
    i1 = lbound(u,1)
    i2 = ubound(u,1)
    j1 = lbound(u,2)
    j2 = ubound(u,2)
    k1 = lbound(u,3)
    k2 = ubound(u,3)
    if ( pj%p%stdlat >= 0.0_rk8 ) then
      do k = k1, k2
        do j = j1, j2
          do i = i1, i2
            tmp = u(i,j,k)*pj%f1(i,j) - v(i,j,k)*pj%f2(i,j)
            v(i,j,k) = real(u(i,j,k)*pj%f2(i,j) + v(i,j,k)*pj%f1(i,j),rkx)
            u(i,j,k) = real(tmp,rkx)
          end do
        end do
      end do
    else
      do k = k1, k2
        do j = j1, j2
          do i = i1, i2
            tmp = u(i,j,k)*pj%f1(i,j) + v(i,j,k)*pj%f2(i,j)
            v(i,j,k) = real(v(i,j,k)*pj%f1(i,j) - u(i,j,k)*pj%f2(i,j),rkx)
            u(i,j,k) = real(tmp,rkx)
          end do
        end do
      end do
    end if
  end subroutine backrotate3_ps

  subroutine rotate2_mc(pj,u,v)
    implicit none
    class(regcm_projection), intent(in) :: pj
    real(rkx), pointer, contiguous, dimension(:,:), intent(inout) :: u, v
    return
  end subroutine rotate2_mc

  subroutine rotate3_mc(pj,u,v)
    implicit none
    class(regcm_projection), intent(in) :: pj
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: u, v
    return
  end subroutine rotate3_mc

  subroutine rotate2_rc(pj,u,v)
    implicit none
    class(regcm_projection), intent(in) :: pj
    real(rkx), pointer, contiguous, dimension(:,:), intent(inout) :: u, v
    integer(ik4) :: i1, i2, j1, j2, i, j
    real(rk8) :: tmp
    i1 = lbound(u,1)
    i2 = ubound(u,1)
    j1 = lbound(u,2)
    j2 = ubound(u,2)
    do j = j1, j2
      do i = i1, i2
        tmp = u(i,j)*pj%f1(i,j) - v(i,j)*pj%f2(i,j)
        v(i,j) = real(v(i,j)*pj%f1(i,j) + u(i,j)*pj%f2(i,j),rkx)
        u(i,j) = real(tmp,rkx)
      end do
    end do
  end subroutine rotate2_rc

  subroutine backrotate2_rc(pj,u,v)
    implicit none
    class(regcm_projection), intent(in) :: pj
    real(rkx), pointer, contiguous, dimension(:,:), intent(inout) :: u, v
    integer(ik4) :: i1, i2, j1, j2, i, j
    real(rk8) :: tmp
    i1 = lbound(u,1)
    i2 = ubound(u,1)
    j1 = lbound(u,2)
    j2 = ubound(u,2)
    do j = j1, j2
      do i = i1, i2
        tmp = u(i,j)*pj%f1(i,j) + v(i,j)*pj%f2(i,j)
        v(i,j) = real(v(i,j)*pj%f1(i,j) - u(i,j)*pj%f2(i,j),rkx)
        u(i,j) = real(tmp,rkx)
      end do
    end do
  end subroutine backrotate2_rc

  subroutine rotate3_rc(pj,u,v)
    implicit none
    class(regcm_projection), intent(in) :: pj
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: u, v
    integer(ik4) :: i1, i2, j1, j2, k1, k2, i, j, k
    real(rk8) :: tmp
    i1 = lbound(u,1)
    i2 = ubound(u,1)
    j1 = lbound(u,2)
    j2 = ubound(u,2)
    k1 = lbound(u,3)
    k2 = ubound(u,3)
    do k = k1, k2
      do j = j1, j2
        do i = i1, i2
          tmp = u(i,j,k)*pj%f1(i,j) - v(i,j,k)*pj%f2(i,j)
          v(i,j,k) = real(v(i,j,k)*pj%f1(i,j) + u(i,j,k)*pj%f2(i,j),rkx)
          u(i,j,k) = real(tmp,rkx)
        end do
      end do
    end do
  end subroutine rotate3_rc

  subroutine backrotate3_rc(pj,u,v)
    implicit none
    class(regcm_projection), intent(in) :: pj
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: u, v
    integer(ik4) :: i1, i2, j1, j2, k1, k2, i, j, k
    real(rk8) :: tmp
    i1 = lbound(u,1)
    i2 = ubound(u,1)
    j1 = lbound(u,2)
    j2 = ubound(u,2)
    k1 = lbound(u,3)
    k2 = ubound(u,3)
    do k = k1, k2
      do j = j1, j2
        do i = i1, i2
          tmp = u(i,j,k)*pj%f1(i,j) + v(i,j,k)*pj%f2(i,j)
          v(i,j,k) = real(v(i,j,k)*pj%f1(i,j) - u(i,j,k)*pj%f2(i,j),rkx)
          u(i,j,k) = real(tmp,rkx)
        end do
      end do
    end do
  end subroutine backrotate3_rc

  subroutine rotate2_rl(pj,u,v)
    implicit none
    class(regcm_projection), intent(in) :: pj
    real(rkx), pointer, contiguous, dimension(:,:), intent(inout) :: u, v
    integer(ik4) :: i1, i2, j1, j2, i, j
    real(rk8) :: tmp
    if ( pj%p%skiprot ) return
    i1 = lbound(u,1)
    i2 = ubound(u,1)
    j1 = lbound(u,2)
    j2 = ubound(u,2)
    do j = j1, j2
      do i = i1, i2
        tmp = u(i,j)*pj%f1(i,j) - v(i,j)*pj%f2(i,j)
        v(i,j) = real(u(i,j)*pj%f2(i,j) + v(i,j)*pj%f1(i,j),rkx)
        u(i,j) = real(tmp,rkx)
      end do
    end do
  end subroutine rotate2_rl

  subroutine backrotate2_rl(pj,u,v)
    implicit none
    class(regcm_projection), intent(in) :: pj
    real(rkx), pointer, contiguous, dimension(:,:), intent(inout) :: u, v
    integer(ik4) :: i1, i2, j1, j2, i, j
    real(rk8) :: tmp
    if ( pj%p%skiprot ) return
    i1 = lbound(u,1)
    i2 = ubound(u,1)
    j1 = lbound(u,2)
    j2 = ubound(u,2)
    do j = j1, j2
      do i = i1, i2
        tmp = u(i,j)*pj%f1(i,j) + v(i,j)*pj%f2(i,j)
        v(i,j) = real(-u(i,j)*pj%f2(i,j) + v(i,j)*pj%f1(i,j),rkx)
        u(i,j) = real(tmp,rkx)
      end do
    end do
  end subroutine backrotate2_rl

  subroutine rotate3_rl(pj,u,v)
    implicit none
    class(regcm_projection), intent(in) :: pj
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: u, v
    integer(ik4) :: i1, i2, j1, j2, k1, k2, i, j, k
    real(rk8) :: tmp
    if ( pj%p%skiprot ) return
    i1 = lbound(u,1)
    i2 = ubound(u,1)
    j1 = lbound(u,2)
    j2 = ubound(u,2)
    k1 = lbound(u,3)
    k2 = ubound(u,3)
    do k = k1, k2
      do j = j1, j2
        do i = i1, i2
          tmp = u(i,j,k)*pj%f1(i,j) - v(i,j,k)*pj%f2(i,j)
          v(i,j,k) = real(u(i,j,k)*pj%f2(i,j) + v(i,j,k)*pj%f1(i,j),rkx)
          u(i,j,k) = real(tmp,rkx)
        end do
      end do
    end do
  end subroutine rotate3_rl

  subroutine backrotate3_rl(pj,u,v)
    implicit none
    class(regcm_projection), intent(in) :: pj
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: u, v
    integer(ik4) :: i1, i2, j1, j2, k1, k2, i, j, k
    real(rk8) :: tmp
    if ( pj%p%skiprot ) return
    i1 = lbound(u,1)
    i2 = ubound(u,1)
    j1 = lbound(u,2)
    j2 = ubound(u,2)
    k1 = lbound(u,3)
    k2 = ubound(u,3)
    do k = k1, k2
      do j = j1, j2
        do i = i1, i2
          tmp = u(i,j,k)*pj%f1(i,j) + v(i,j,k)*pj%f2(i,j)
          v(i,j,k) = real(-u(i,j,k)*pj%f2(i,j) + v(i,j,k)*pj%f1(i,j),rkx)
          u(i,j,k) = real(tmp,rkx)
        end do
      end do
    end do
  end subroutine backrotate3_rl

  pure subroutine mapfac_rl(pj,xlat,xlon,xmap)
    implicit none
    class(regcm_projection), intent(in) :: pj
    real(rkx), pointer, contiguous, dimension(:,:), intent(in) :: xlat, xlon
    real(rkx), pointer, contiguous, dimension(:,:), intent(inout) :: xmap
    xmap = fac_rl(pj,xlat,xlon)
  end subroutine mapfac_rl

  pure subroutine mapfac_lc(pj,xlat,xlon,xmap)
    implicit none
    class(regcm_projection), intent(in) :: pj
    real(rkx), pointer, contiguous, dimension(:,:), intent(in) :: xlat, xlon
    real(rkx), pointer, contiguous, dimension(:,:), intent(inout) :: xmap
    if ( .false. ) then
      xmap(:,:) = xlon(:,:) ! Avoid unused argument warning
    end if
    xmap = fac_lc(pj,xlat)
  end subroutine mapfac_lc

  pure subroutine mapfac_ps(pj,xlat,xlon,xmap)
    implicit none
    class(regcm_projection), intent(in) :: pj
    real(rkx), pointer, contiguous, dimension(:,:), intent(in) :: xlat, xlon
    real(rkx), pointer, contiguous, dimension(:,:), intent(inout) :: xmap
    if ( .false. ) then
      xmap(:,:) = xlon(:,:) ! Avoid unused argument warning
    end if
    xmap = fac_ps(pj,xlat)
  end subroutine mapfac_ps

  pure subroutine mapfac_mc(pj,xlat,xlon,xmap)
    implicit none
    class(regcm_projection), intent(in) :: pj
    real(rkx), pointer, contiguous, dimension(:,:), intent(in) :: xlat, xlon
    real(rkx), pointer, contiguous, dimension(:,:), intent(inout) :: xmap
    if ( .false. ) then
      xmap(:,:) = xlon(:,:) ! Avoid unused argument warning
    end if
    xmap = fac_mc(pj,xlat)
  end subroutine mapfac_mc

  pure subroutine mapfac_ll(pj,xlat,xlon,xmap)
    implicit none
    class(regcm_projection), intent(in) :: pj
    real(rkx), pointer, contiguous, dimension(:,:), intent(in) :: xlat, xlon
    real(rkx), pointer, contiguous, dimension(:,:), intent(inout) :: xmap
    xmap = 1.0_rkx
  end subroutine mapfac_ll

  pure subroutine mapfac_rc(pj,xlat,xlon,xmap)
    implicit none
    class(regcm_projection), intent(in) :: pj
    real(rkx), pointer, contiguous, dimension(:,:), intent(in) :: xlat, xlon
    real(rkx), pointer, contiguous, dimension(:,:), intent(inout) :: xmap
    xmap = fac_rc(pj,xlat,xlon)
  end subroutine mapfac_rc

end module mod_projections
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
