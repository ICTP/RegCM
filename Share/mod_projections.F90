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

module mod_projections

  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_memutil

#ifdef F2008
  use , intrinsic :: iso_fortran_env
#endif

  implicit none

  private

  public :: rounder

  type anyprojparams
    character(len=6) :: pcode
    real(rk8) :: ds
    real(rk8) :: clat , clon
    real(rk8) :: plat , plon
    real(rk8) :: trlat1 , trlat2
    integer(ik4) :: nlat , nlon
    logical :: staggerx , staggery
    logical :: rotparam
  end type anyprojparams

  type regcm_projdata
    character(len=6) :: code
    real(rk8) :: stdlon , stdlat
    real(rk8) :: truelat1 , truelat2 , tl2r , ctl1r
    real(rk8) :: chi1 , chi2 , xct1 , tanchi1h , schi1 , tchi1
    real(rk8) :: colat1 , colat2 , nfac
    real(rk8) :: rsw , rebydx , hemi
    real(rk8) :: reflon , dlon , dlat , scale_top
    real(rk8) :: polei , polej
    real(rk8) :: rlon0 , rlat0 , phi0 , lam0
    real(rk8) :: xoff , yoff
    real(rk8) :: zsinpol , zcospol , zlampol
    real(rk8) :: pollam , polcphi , polsphi
    real(rk8) :: conefac , rconefac
    integer(ik4) :: nlat , nlon
    logical :: lamtan
    logical :: skiprot
  end type regcm_projdata

  type regcm_projection
    private
    type(regcm_projdata) , pointer :: p
    procedure(transform) , pass(pj) , pointer , public :: llij => NULL()
    procedure(transform) , pass(pj) , pointer , public :: ijll => NULL()
    procedure(map_factor) , pass(pj) , pointer , public :: mapfac => NULL()
    procedure(rotate2) , pass(pj) , pointer , public :: uvrotate2 => NULL()
    procedure(rotate2) , pass(pj) , pointer , public :: uvbkrotate2 => NULL()
    procedure(rotate3) , pass(pj) , pointer , public :: uvrotate3 => NULL()
    procedure(rotate3) , pass(pj) , pointer , public :: uvbkrotate3 => NULL()

    real(rk8) , pointer , dimension(:,:) :: f1 , f2 , f3 , f4
    real(rk8) , pointer , dimension(:,:) :: f5 , f6 , f7 , f8

    contains

      procedure , public :: initialize
      procedure , public :: rotation_angle
      procedure , public :: wind_rotate
      procedure , public :: wind2_rotate
      procedure , public :: wind_antirotate
      procedure , public :: wind2_antirotate
      procedure , public :: rl00
      procedure , public :: destruct

  end type regcm_projection

  abstract interface
    pure subroutine transform(pj,a,b,c,d)
      import
      class(regcm_projection) , intent(in) :: pj
#ifdef SINGLE_PRECISION_REAL
#ifdef F2008
      integer , parameter :: rk4  = REAL32
#else
      integer , parameter :: rk4  = kind(1.0)
#endif
      integer , parameter :: rkx = rk4
#else
#ifdef F2008
      integer , parameter :: rk8  = REAL64
#else
      integer , parameter :: rk8  = selected_real_kind(2*precision(1.0))
#endif
      integer , parameter :: rkx = rk8
#endif
      real(rkx) , intent(in) :: a , b
      real(rkx) , intent(out) :: c , d
    end subroutine transform
  end interface

  abstract interface
    subroutine rotate3(pj,u,v)
      import
      class(regcm_projection) , intent(in) :: pj
#ifdef SINGLE_PRECISION_REAL
#ifdef F2008
      integer , parameter :: rk4  = REAL32
#else
      integer , parameter :: rk4  = kind(1.0)
#endif
      integer , parameter :: rkx = rk4
#else
#ifdef F2008
      integer , parameter :: rk8  = REAL64
#else
      integer , parameter :: rk8  = selected_real_kind(2*precision(1.0))
#endif
      integer , parameter :: rkx = rk8
#endif
      real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: u , v
    end subroutine rotate3
  end interface

  abstract interface
    subroutine rotate2(pj,u,v)
      import
      class(regcm_projection) , intent(in) :: pj
#ifdef SINGLE_PRECISION_REAL
#ifdef F2008
      integer , parameter :: rk4  = REAL32
#else
      integer , parameter :: rk4  = kind(1.0)
#endif
      integer , parameter :: rkx = rk4
#else
#ifdef F2008
      integer , parameter :: rk8  = REAL64
#else
      integer , parameter :: rk8  = selected_real_kind(2*precision(1.0))
#endif
      integer , parameter :: rkx = rk8
#endif
      real(rkx) , pointer , dimension(:,:) , intent(inout) :: u , v
    end subroutine rotate2
  end interface

  abstract interface
    pure subroutine map_factor(pj,a,b,c)
      import
      class(regcm_projection) , intent(in) :: pj
#ifdef SINGLE_PRECISION_REAL
#ifdef F2008
      integer , parameter :: rk4  = REAL32
#else
      integer , parameter :: rk4  = kind(1.0)
#endif
      integer , parameter :: rkx = rk4
#else
#ifdef F2008
      integer , parameter :: rk8  = REAL64
#else
      integer , parameter :: rk8  = selected_real_kind(2*precision(1.0))
#endif
      integer , parameter :: rkx = rk8
#endif
      real(rkx) , pointer , dimension(:,:) , intent(in) :: a , b
      real(rkx) , pointer , dimension(:,:) , intent(inout) :: c
    end subroutine map_factor
  end interface

  public :: anyprojparams , regcm_projection

  contains

  subroutine rl00(pj,lat0,lon0)
    implicit none
    class(regcm_projection) , intent(in) :: pj
    real(rkx) , intent(out) :: lat0 , lon0
    lat0 = pj%p%rlat0
    lon0 = pj%p%rlon0
  end subroutine rl00

  subroutine wind_rotate(pj,u,v)
    implicit none
    class(regcm_projection) , intent(in) :: pj
    real(rkx) , dimension(:,:,:) , pointer , intent(inout) :: u , v
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
    class(regcm_projection) , intent(in) :: pj
    real(rkx) , dimension(:,:,:) , pointer , intent(inout) :: u , v
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
    class(regcm_projection) , intent(in) :: pj
    real(rkx) , dimension(:,:) , pointer , intent(inout) :: u , v
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
    class(regcm_projection) , intent(in) :: pj
    real(rkx) , dimension(:,:) , pointer , intent(inout) :: u , v
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
    class(regcm_projection) , intent(in) :: pj
    real(rk8) , dimension(:) , intent(in) :: lon
    real(rk8) , dimension(:) , intent(in) :: lat
    real(rk8) , dimension(:,:) , intent(out) :: p
    integer(ik4) :: i , j , nlon , nlat
    real(rk8) :: f1 , f2
    nlon = size(lon)
    nlat = size(lat)
    if ( size(p,1) /= nlon .or. size(p,2) /= nlat ) then
      return
    end if
    select case (pj%p%code)
      case ('LAMCON')
        do i = 1 , nlat
          do j = 1 , nlon
            p(j,i) = uvrot_lc(pj,lon(j))
          end do
        end do
      case ('ROTMER')
        do i = 1 , nlat
          do j = 1 , nlon
            call uvrot_rc(pj,lat(i),lon(j),f1,f2)
            p(j,i) = acos(f1)
          end do
        end do
      case ('POLSTR')
        do i = 1 , nlat
          do j = 1 , nlon
            p(j,i) = uvrot_ps(pj,lon(j))
          end do
        end do
      case default
        p(:,:) = 1.0_rk8
    end select
  end subroutine rotation_angle

  subroutine destruct(pj)
    implicit none
    class(regcm_projection) , intent(inout) :: pj
    if ( associated(pj%p) ) deallocate(pj%p)
  end subroutine destruct

  subroutine initialize(pj,pjpara)
    implicit none
    class(regcm_projection) , intent(out) :: pj
    type(anyprojparams) , intent(in) :: pjpara
    real(rk8) :: ci , cj
    if ( .not. associated(pj%p) ) then
      allocate(pj%p)
    end if
    ci = real(pjpara%nlon,rk8) * 0.5_rk8
    cj = real(pjpara%nlat,rk8) * 0.5_rk8
    if ( .not. pjpara%staggerx ) ci = ci - 0.5_rk8
    if ( .not. pjpara%staggery ) cj = cj - 0.5_rk8
    pj%p%code = pjpara%pcode
    pj%p%nlon = pjpara%nlon
    pj%p%nlat = pjpara%nlat
    select case (pjpara%pcode)
      case ('LAMCON')
        call setup_lcc(pj,pjpara%clon,pjpara%clat,ci,cj,pjpara%ds, &
                       pjpara%clon,pjpara%trlat1,pjpara%trlat2,pjpara%rotparam)
        pj%llij => llij_lc
        pj%ijll => ijll_lc
        pj%uvrotate2 => rotate2_lc
        pj%uvrotate3 => rotate3_lc
        pj%uvbkrotate2 => backrotate2_lc
        pj%uvbkrotate3 => backrotate3_lc
        pj%mapfac => mapfac_lc
      case ('POLSTR')
          call setup_plr(pj,pjpara%clon,pjpara%clat,ci,cj, &
                       pjpara%ds,pjpara%clon,pjpara%rotparam)
        pj%llij => llij_ps
        pj%ijll => ijll_ps
        pj%uvrotate2 => rotate2_ps
        pj%uvrotate3 => rotate3_ps
        pj%uvbkrotate2 => backrotate2_ps
        pj%uvbkrotate3 => backrotate3_ps
        pj%mapfac => mapfac_ps
      case ('NORMER')
        call setup_mrc(pj,pjpara%clon,pjpara%clat,ci,cj, &
                       pjpara%ds)
        pj%llij => llij_mc
        pj%ijll => ijll_mc
        pj%uvrotate2 => rotate2_mc
        pj%uvrotate3 => rotate3_mc
        pj%uvbkrotate2 => rotate2_mc
        pj%uvbkrotate3 => rotate3_mc
        pj%mapfac => mapfac_mc
      case ('ROTMER')
        call setup_rmc(pj,pjpara%clon,pjpara%clat,ci,cj, &
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
    type(regcm_projection) , intent(inout) :: pj
    real(rk8) , intent(in) :: ci , cj , clat , clon , ds
    pj%p%dlon = raddeg * ds / earthrad
    pj%p%dlat = pj%p%dlon
    pj%p%rlat0 = clat - cj*pj%p%dlat
    pj%p%rlon0 = clon - ci*pj%p%dlon
    if ( pj%p%rlat0 >  deg90 )  pj%p%rlat0 = deg90 - pj%p%rlat0
    if ( pj%p%rlat0 < -deg90 )  pj%p%rlat0 = pj%p%rlat0 + deg90
    if ( pj%p%rlon0 >  deg180 ) pj%p%rlon0 = deg360 - pj%p%rlon0
    if ( pj%p%rlon0 < -deg180 ) pj%p%rlon0 = pj%p%rlon0 + deg360
  end subroutine setup_ll

  pure subroutine ijll_ll(pj,i,j,lat,lon)
    implicit none
    class(regcm_projection) , intent(in) :: pj
    real(rkx) , intent(in) :: i , j
    real(rkx) , intent(out) :: lat , lon
    lat = pj%p%rlat0 + (j-1.0_rk8) * pj%p%dlat
    lon = pj%p%rlon0 + (i-1.0_rk8) * pj%p%dlon
    if ( lat >  90.0_rkx ) lat = 90.0_rkx - lat
    if ( lat < -90.0_rkx ) lat = lat + 90.0_rkx
    if ( lon >  180.0_rkx ) lon = lon - 360.0_rkx
    if ( lon < -180.0_rkx ) lon = lon + 360.0_rkx
  end subroutine ijll_ll

  pure subroutine llij_ll(pj,lat,lon,i,j)
    implicit none
    class(regcm_projection) , intent(in) :: pj
    real(rkx) , intent(in) :: lat , lon
    real(rkx) , intent(out) :: i , j
    i = real(( lon - pj%p%rlon0 ) / pj%p%dlon + 1.0_rk8,rkx)
    j = real(( lat - pj%p%rlat0 ) / pj%p%dlat + 1.0_rk8,rkx)
  end subroutine llij_ll

  subroutine setup_rll(pj,clon,clat,ci,cj,ds,plon,plat,luvrot)
    implicit none
    type(regcm_projection) , intent(inout) :: pj
    real(rk8) , intent(in) :: ci , cj , clat , clon , plon , plat , ds
    logical , intent(in) :: luvrot
    real(rk8) :: phi , lam , dlam
    real(rk8) :: rotlam , rotphi
    real(rkx) :: lon , lat , ri , rj
    integer(ik4) :: i , j
    if ( abs(plat-90.0_rk8) < 0.001 ) pj%p%skiprot = .true.
    pj%p%dlon = raddeg * ds / earthrad
    pj%p%dlat = pj%p%dlon
    phi = degrad*plat
    lam = degrad*plon
    call get_equator(phi,lam,pj%p%phi0,pj%p%lam0)
    phi = clat
    lam = clon
    if ( phi >  deg90 )  phi = deg90 - phi
    if ( phi < -deg90 )  phi = phi + deg90
    if ( lam >  deg180 ) lam = lam - deg360
    if ( lam < -deg180 ) lam = lam + deg360
    phi = degrad * phi
    lam = degrad * lam
    pj%p%rlat0 = asin(-cos(phi)*sin(pj%p%phi0)*cos(lam-pj%p%lam0) + &
                    sin(phi)*cos(pj%p%phi0))
    if ( abs(abs(pj%p%rlat0)-halfpi) > 1.0e-7_rk8 .and. &
         abs(pj%p%phi0) > 1.0e-7_rk8 ) then
      pj%p%rlon0 = (sin(phi)-cos(pj%p%phi0)*sin(pj%p%rlat0)) / &
                 (sin(pj%p%phi0)*cos(pj%p%rlat0))
      if ( pj%p%rlon0 < -1.0_rk8 .and. &
           pj%p%rlon0 > -1.00001_rk8 ) pj%p%rlon0 = -1.0_rk8
      if ( pj%p%rlon0 >  1.0_rk8 .and. &
           pj%p%rlon0 <  1.00001_rk8 ) pj%p%rlon0 =  1.0_rk8
      pj%p%rlon0 = acos(pj%p%rlon0)
      if ( lam < pj%p%lam0 ) then
        pj%p%rlon0 = -pj%p%rlon0
      end if
    else
      pj%p%rlon0 = lam
    end if
    pj%p%rlat0 = raddeg*pj%p%rlat0 - cj*pj%p%dlat
    pj%p%rlon0 = raddeg*pj%p%rlon0 - ci*pj%p%dlon
    if ( pj%p%rlat0 >  deg90 )  pj%p%rlat0 = deg90 - pj%p%rlat0
    if ( pj%p%rlat0 < -deg90 )  pj%p%rlat0 = pj%p%rlat0 + deg90
    if ( pj%p%rlon0 >  deg180 ) pj%p%rlon0 = pj%p%rlon0 - deg360
    if ( pj%p%rlon0 < -deg180 ) pj%p%rlon0 = pj%p%rlon0 + deg360
    if ( luvrot ) then
      call getmem2d(pj%f1,1,pj%p%nlon,1,pj%p%nlat,'projections:f1')
      call getmem2d(pj%f2,1,pj%p%nlon,1,pj%p%nlat,'projections:f2')
      call getmem2d(pj%f3,1,pj%p%nlon,1,pj%p%nlat,'projections:f3')
      call getmem2d(pj%f4,1,pj%p%nlon,1,pj%p%nlat,'projections:f4')
      call getmem2d(pj%f5,1,pj%p%nlon,1,pj%p%nlat,'projections:f5')
      call getmem2d(pj%f6,1,pj%p%nlon,1,pj%p%nlat,'projections:f6')
      call getmem2d(pj%f7,1,pj%p%nlon,1,pj%p%nlat,'projections:f7')
      call getmem2d(pj%f8,1,pj%p%nlon,1,pj%p%nlat,'projections:f8')
      do j = 1 , pj%p%nlat
        rotphi = degrad*(pj%p%rlat0 + (j-1)*pj%p%dlat)
        do i = 1 , pj%p%nlon
          rotlam = degrad*(pj%p%rlon0 + (i-1)*pj%p%dlon)
          ri = i
          rj = j
          call ijll_rl(pj,ri,rj,lat,lon)
          lam = degrad*lon
          phi = degrad*lat
          dlam = lam - pj%p%lam0
          if ( .not. pj%p%skiprot ) then
            if ( abs(rotlam) < 1.0e-7 ) then
              if ( lam > halfpi .or. lam < -halfpi ) then
                pj%f1(i,j) = -1.0_rk8
                pj%f2(i,j) = 0.0_rk8
                pj%f3(i,j) = 0.0_rk8
                pj%f4(i,j) = -1.0_rk8
              else
                pj%f1(i,j) = 1.0_rk8
                pj%f2(i,j) = 0.0_rk8
                pj%f3(i,j) = 0.0_rk8
                pj%f4(i,j) = 1.0_rk8
              end if
            else
              pj%f1(i,j) = sin(dlam)*sin(pj%p%phi0) / cos(rotphi)
              pj%f2(i,j) = (cos(dlam)*sin(phi)*sin(pj%p%phi0) + &
                            cos(pj%p%phi0)*cos(phi)) / &
                            cos(rotphi)
              pj%f3(i,j) = (cos(pj%p%phi0)*cos(rotphi) - &
                            sin(pj%p%phi0)*sin(rotphi)*cos(rotlam)) / &
                           (sin(pj%p%phi0)*sin(rotlam))
              pj%f4(i,j) = -cos(phi) / (sin(pj%p%phi0)*sin(rotlam))
            end if
            if ( abs(cos(phi)) > 1.0e-7_rk8 ) then
              if ( abs(sin(dlam)) > 1.0e-7_rk8 ) then
                pj%f5(i,j) = -sin(pj%p%phi0)*sin(rotlam)/cos(phi)
                pj%f6(i,j) = (cos(pj%p%phi0)*cos(rotphi) - &
                              sin(pj%p%phi0)*sin(rotphi)*cos(rotlam))/cos(phi)
                pj%f7(i,j) = cos(rotphi)/(sin(dlam)*sin(pj%p%phi0))
                pj%f8(i,j) = (cos(dlam)*sin(pj%p%phi0)*sin(phi) + &
                              cos(pj%p%phi0)*cos(phi))/cos(rotphi)
              else
                if ( dlam < -halfpi .or. dlam > halfpi ) then
                  pj%f5(i,j) = -1.0_rk8
                  pj%f6(i,j) = 0.0_rk8
                  pj%f7(i,j) = -1.0_rk8
                  pj%f8(i,j) = 0.0_rk8
                else
                  pj%f5(i,j) = 1.0_rk8
                  pj%f6(i,j) = 0.0_rk8
                  pj%f7(i,j) = 1.0_rk8
                  pj%f8(i,j) = 0.0_rk8
                end if
              end if
            else
              if ( pj%p%phi0 < 0.0 ) then
                pj%f5(i,j) = -lam
                pj%f6(i,j) = 0.0_rk8
                pj%f7(i,j) = 0.0_rk8
                pj%f8(i,j) = 0.0_rk8
              else
                pj%f5(i,j) = lam
                pj%f6(i,j) = 0.0_rk8
                pj%f7(i,j) = 0.0_rk8
                pj%f8(i,j) = 0.0_rk8
              end if
            end if
          end if
        end do
      end do
    end if
  end subroutine setup_rll

  pure subroutine ijll_rl(pj,i,j,lat,lon)
    class(regcm_projection) , intent(in) :: pj
    real(rkx) , intent(in) :: i , j
    real(rkx) , intent(out) :: lat , lon
    real(rk8) :: phi , lam
    phi = pj%p%rlat0 + (j-1.0_rk8) * pj%p%dlat
    lam = pj%p%rlon0 + (i-1.0_rk8) * pj%p%dlon
    if ( phi >  deg90 )  phi = deg90 - phi
    if ( phi < -deg90 )  phi = phi + deg90
    if ( lam >  deg180 ) lam = lam - deg360
    if ( lam < -deg180 ) lam = lam + deg360
    phi = degrad * phi
    lam = degrad * lam
    if ( abs(pj%p%phi0) > 1.0e-7_rk8 ) then
      lat = asin(cos(phi)*sin(pj%p%phi0)*cos(lam) + sin(phi)*cos(pj%p%phi0))
      if ( abs(abs(lat)-halfpi) > 1.0e-7_rk8 ) then
        lon = (-sin(phi)*sin(pj%p%phi0) +  &
               cos(pj%p%phi0)*cos(lam)*cos(phi))/cos(lat)
        if ( lon < -1.0_rk8 .and. lon > -1.00001_rk8 ) lon = -1.0_rk8
        if ( lon >  1.0_rk8 .and. lon <  1.00001_rk8 ) lon =  1.0_rk8
        if ( lam < 0.0_rk8 ) then
          lon = pj%p%lam0 - acos(lon)
        else
          lon = pj%p%lam0 + acos(lon)
        end if
      else
        lon = lam
      end if
    else
      lat = asin(cos(phi)*sin(pj%p%phi0)*cos(lam) + sin(phi)*cos(pj%p%phi0))
      lon = lam
    end if
    lat = raddeg*lat
    lon = raddeg*lon
    if ( lat >  90.0_rkx ) lat = 90.0_rkx - lat
    if ( lat < -90.0_rkx ) lat = lat + 90.0_rkx
    if ( lon >  180.0_rkx ) lon = lon - 360.0_rkx
    if ( lon < -180.0_rkx ) lon = lon + 360.0_rkx
  end subroutine ijll_rl

  pure subroutine llij_rl(pj,lat,lon,i,j)
    implicit none
    class(regcm_projection) , intent(in) :: pj
    real(rkx) , intent(in) :: lat , lon
    real(rkx) , intent(out) :: i , j
    real(rk8) :: rlat , rlon , phi , lam
    phi = lat
    lam = lon
    if ( lam >  deg180 ) lam = lam - deg360
    if ( lam < -deg180 ) lam = lam + deg360
    if ( phi >  deg90 )  phi = deg90 - phi
    if ( phi < -deg90 )  phi = phi + deg90
    phi = degrad * phi
    lam = degrad * lam
    rlat = asin(-cos(phi)*sin(pj%p%phi0)*cos(lam-pj%p%lam0) + &
           sin(phi)*cos(pj%p%phi0))
    if ( abs(abs(rlat)-halfpi) > 1.0e-7_rk8 .and. &
         abs(pj%p%phi0) > 1.0e-7_rk8 ) then
      rlon = (sin(phi)-cos(pj%p%phi0)*sin(rlat))/(sin(pj%p%phi0)*cos(rlat))
      if ( rlon < -1.0_rk8 .and. rlon > -1.00001_rk8 ) rlon = -1.0_rk8
      if ( rlon >  1.0_rk8 .and. rlon <  1.00001_rk8 ) rlon =  1.0_rk8
      rlon = acos(rlon)
      if ( lam < pj%p%lam0 ) then
        rlon = -rlon
      end if
    else
      rlon = lam
    end if
    rlon = raddeg*rlon
    rlat = raddeg*rlat
    i = real(( rlon - pj%p%rlon0 ) / pj%p%dlon + 1.0_rk8,rkx)
    j = real(( rlat - pj%p%rlat0 ) / pj%p%dlat + 1.0_rk8,rkx)
  end subroutine llij_rl

  subroutine setup_lcc(pj,clon,clat,ci,cj,ds,slon,trlat1,trlat2,luvrot)
    implicit none
    type(regcm_projection) , intent(inout) :: pj
    real(rk8) , intent(in) :: ci , cj , slon , clat , clon , ds , &
                              trlat1 , trlat2
    logical , intent(in) :: luvrot
    real(rk8) :: arg , deltalon1 , tl1r , tl2r
    real(rkx) :: ri , rj , lat , lon
    integer(ik4) :: i , j

    pj%p%stdlon = slon
    pj%p%stdlat = clat
    pj%p%truelat1 = trlat1
    pj%p%truelat2 = trlat2
    tl1r = pj%p%truelat1*degrad
    tl2r = pj%p%truelat2*degrad
    pj%p%colat1  = degrad*(deg90 - pj%p%truelat1)
    pj%p%colat2  = degrad*(deg90 - pj%p%truelat2)
    pj%p%nfac = (log(sin(pj%p%colat1))         - log(sin(pj%p%colat2))) / &
                (log(tan(pj%p%colat1*0.5_rk8)) - log(tan(pj%p%colat2*0.5_rk8)))
    if ( pj%p%stdlat > 0.0_rk8 ) then
      pj%p%hemi =  1.0_rk8
    else
      pj%p%hemi = -1.0_rk8
    end if
    pj%p%rebydx = earthrad / ds
    if ( abs(pj%p%truelat1-pj%p%truelat2) > 0.1_rk8 ) then
      pj%p%conefac = log10(cos(tl1r)) - log10(cos(tl2r))
      pj%p%conefac = pj%p%conefac / &
                (log10(tan((deg45-abs(pj%p%truelat1)/2.0_rk8)*degrad)) - &
                 log10(tan((deg45-abs(pj%p%truelat2)/2.0_rk8)*degrad)))
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
    pj%p%xct1 = tan((deg90*pj%p%hemi-pj%p%truelat1)*degrad*0.5_rk8)
    pj%p%rsw = pj%p%rebydx * pj%p%ctl1r*pj%p%rconefac * &
           (tan((deg90*pj%p%hemi-clat)*degrad*0.5_rk8)/pj%p%xct1)**pj%p%conefac
    arg = pj%p%conefac*(deltalon1*degrad)
    pj%p%polei = pj%p%hemi*ci - pj%p%hemi * pj%p%rsw * sin(arg)
    pj%p%polej = pj%p%hemi*cj + pj%p%rsw * cos(arg)
    pj%p%chi1 = (deg90 - pj%p%hemi*pj%p%truelat1)*degrad
    pj%p%chi2 = (deg90 - pj%p%hemi*pj%p%truelat2)*degrad
    pj%p%tanchi1h = tan(pj%p%chi1*0.5_rk8)
    pj%p%tchi1 = tan(pj%p%chi1)
    pj%p%schi1 = sin(pj%p%chi1)
    if ( luvrot ) then
      call getmem2d(pj%f1,1,pj%p%nlon,1,pj%p%nlat,'projections:f1')
      call getmem2d(pj%f2,1,pj%p%nlon,1,pj%p%nlat,'projections:f2')
      call getmem2d(pj%f3,1,pj%p%nlon,1,pj%p%nlat,'projections:f3')
      do j = 1 , pj%p%nlat
        do i = 1 , pj%p%nlon
          ri = i
          rj = j
          call ijll_lc(pj,ri,rj,lat,lon)
          pj%f3(i,j) = uvrot_lc(pj,lon)
          pj%f1(i,j) = cos(pj%f3(i,j))
          pj%f2(i,j) = sin(pj%f3(i,j))
        end do
      end do
    end if
  end subroutine setup_lcc

  pure subroutine ijll_lc(pj,i,j,lat,lon)
    implicit none
    class(regcm_projection) , intent(in) :: pj
    real(rkx) , intent(in) :: i , j
    real(rkx) , intent(out) :: lat , lon
    real(rk8) :: chi
    real(rk8) :: inew , jnew , xx , yy , r2 , r

    inew = pj%p%hemi * i
    jnew = pj%p%hemi * j
    xx = inew - pj%p%polei
    yy = pj%p%polej - jnew
    r2 = (xx*xx + yy*yy)
    r = sqrt(r2)/pj%p%rebydx
    if ( abs(r2) < dlowval ) then
      lat = real(pj%p%hemi*deg90,rkx)
      lon = real(pj%p%stdlon,rkx)
    else
      lon = real(pj%p%stdlon + &
           raddeg * atan2(pj%p%hemi*xx,yy)*pj%p%rconefac,rkx)
      if ( pj%p%lamtan ) then
        chi = 2.0_rk8 * &
          atan(((r/pj%p%tchi1)**pj%p%rconefac)*pj%p%tanchi1h)
      else
        chi = 2.0_rk8 * &
          atan(((r*pj%p%conefac/pj%p%schi1)**pj%p%rconefac)*pj%p%tanchi1h)
      end if
      lat = real((deg90-chi*raddeg)*pj%p%hemi,rkx)
    end if
    if ( lon >  180.0_rkx ) lon = lon - 360.0_rkx
    if ( lon < -180.0_rkx ) lon = lon + 360.0_rkx
  end subroutine ijll_lc

  pure subroutine llij_lc(pj,lat,lon,i,j)
    implicit none
    class(regcm_projection) , intent(in) :: pj
    real(rkx) , intent(in) :: lat , lon
    real(rkx) , intent(out) :: i , j
    real(rk8) :: arg , deltalon , rm

    deltalon = lon - pj%p%stdlon
    if ( deltalon > +deg180 ) deltalon = deltalon - deg360
    if ( deltalon < -deg180 ) deltalon = deltalon + deg360
    rm = pj%p%rebydx * pj%p%ctl1r * pj%p%rconefac * &
           (tan((deg90*pj%p%hemi-lat)*degrad*0.5_rk8)/pj%p%xct1)**pj%p%conefac
    arg = pj%p%conefac*(deltalon*degrad)
    i = real(pj%p%hemi*(pj%p%polei + pj%p%hemi * rm * sin(arg)),rkx)
    j = real(pj%p%hemi*(pj%p%polej - rm * cos(arg)),rkx)
  end subroutine llij_lc

  subroutine setup_plr(pj,clon,clat,ci,cj,ds,slon,luvrot)
    implicit none
    type(regcm_projection) , intent(inout) :: pj
    real(rk8) , intent(in) :: clat , clon , cj , ci , ds , slon
    logical , intent(in) :: luvrot
    real(rk8) :: ala1 , alo1
    real(rkx) :: lat , lon , ri , rj
    integer(ik4) :: i , j

    pj%p%stdlon = slon
    pj%p%stdlat = clat
    if ( pj%p%stdlat > 0.0_rk8 ) then
      pj%p%hemi = 1.0_rk8
    else
      pj%p%hemi = -1.0_rk8
    end if
    pj%p%rebydx = earthrad / ds
    pj%p%reflon = pj%p%stdlon + deg90
    ala1 = clat*degrad
    alo1 = (clon-pj%p%reflon)*degrad
    pj%p%scale_top = 1.0_rk8 + pj%p%hemi * sin(ala1)
    pj%p%rsw = pj%p%rebydx * &
         cos(ala1)*pj%p%scale_top/(1.0_rk8+pj%p%hemi*sin(ala1))
    pj%p%polei = ci - pj%p%rsw * cos(alo1)
    pj%p%polej = cj - pj%p%hemi * pj%p%rsw * sin(alo1)
    if ( luvrot ) then
      call getmem2d(pj%f1,1,pj%p%nlon,1,pj%p%nlat,'projections:f1')
      call getmem2d(pj%f2,1,pj%p%nlon,1,pj%p%nlat,'projections:f2')
      call getmem2d(pj%f3,1,pj%p%nlon,1,pj%p%nlat,'projections:f3')
      do j = 1 , pj%p%nlat
        do i = 1 , pj%p%nlon
          ri = i
          rj = j
          call ijll_ps(pj,ri,rj,lat,lon)
          pj%f3(i,j) = uvrot_ps(pj,lon)
          pj%f1(i,j) = cos(pj%f3(i,j))
          pj%f2(i,j) = sin(pj%f3(i,j))
        end do
      end do
    end if
  end subroutine setup_plr

  pure subroutine llij_ps(pj,lat,lon,i,j)
    implicit none
    class(regcm_projection) , intent(in) :: pj
    real(rkx) , intent(in) :: lat , lon
    real(rkx) , intent(out) :: i , j
    real(rk8) :: ala , alo , rm , deltalon

    deltalon = lon - pj%p%reflon
    if ( deltalon > +deg180 ) deltalon = deltalon - deg360
    if ( deltalon < -deg180 ) deltalon = deltalon + deg360
    alo = deltalon * degrad
    ala = lat * degrad
    rm = pj%p%rebydx * cos(ala) * &
         pj%p%scale_top/(1.0_rk8 + pj%p%hemi * sin(ala))
    i = real(pj%p%polei + rm * cos(alo),rkx)
    j = real(pj%p%polej + pj%p%hemi * rm * sin(alo),rkx)
  end subroutine llij_ps

  pure subroutine ijll_ps(pj,i,j,lat,lon)
    implicit none
    class(regcm_projection) , intent(in) :: pj
    real(rkx) , intent(in) :: i , j
    real(rkx) , intent(out) :: lat , lon
    real(rk8) :: xx , yy , r2 , gi2 , arcc

    xx = i - pj%p%polei
    yy = (j - pj%p%polej) * pj%p%hemi
    r2 = xx*xx + yy*yy
    if ( abs(r2) < dlowval ) then
      lat = real(pj%p%hemi*deg90,rkx)
      lon = real(pj%p%reflon,rkx)
    else
      gi2 = (pj%p%rebydx * pj%p%scale_top)**2
      lat = real(raddeg * pj%p%hemi * asin((gi2-r2)/(gi2+r2)),rkx)
      arcc = acos(xx/sqrt(r2))
      if ( yy > 0.0_rk8 ) then
        lon = real(pj%p%reflon + raddeg * arcc,rkx)
      else
        lon = real(pj%p%reflon - raddeg * arcc,rkx)
      end if
    end if
    if ( lon >  180.0_rkx ) lon = lon - 360.0_rkx
    if ( lon < -180.0_rkx ) lon = lon + 360.0_rkx
  end subroutine ijll_ps

  subroutine setup_mrc(pj,clon,clat,ci,cj,ds)
    implicit none
    type(regcm_projection) , intent(inout) :: pj
    real(rk8) , intent(in) :: clat , clon , cj , ci , ds
    real(rk8) :: clain

    pj%p%stdlon = clon
    clain = cos(clat*degrad)
    pj%p%dlon = ds / (earthrad * clain)
    pj%p%rsw = 0.0_rk8
    if ( abs(clat) > dlowval ) then
      pj%p%rsw = (log(tan(0.5_rk8*((clat+deg90)*degrad))))/pj%p%dlon
    end if
    pj%p%polei = ci
    pj%p%polej = cj
  end subroutine setup_mrc

  pure subroutine llij_mc(pj,lat,lon,i,j)
    implicit none
    class(regcm_projection) , intent(in) :: pj
    real(rkx) , intent(in) :: lat , lon
    real(rkx) , intent(out) :: i , j
    real(rk8) :: deltalon

    deltalon = lon - pj%p%stdlon
    if ( deltalon > +deg180 ) deltalon = deltalon - deg360
    if ( deltalon < -deg180 ) deltalon = deltalon + deg360
    i = real(pj%p%polei + (deltalon/(pj%p%dlon*raddeg)),rkx)
    j = real(pj%p%polej + &
         (log(tan(0.5_rk8*((lat+deg90)*degrad)))) / pj%p%dlon - pj%p%rsw,rkx)
  end subroutine llij_mc

  pure subroutine ijll_mc(pj,i,j,lat,lon)
    implicit none
    class(regcm_projection) , intent(in) :: pj
    real(rkx) , intent(in) :: i , j
    real(rkx) , intent(out) :: lat , lon

    lat = real(2.0_rk8 * &
        atan(exp(pj%p%dlon*(pj%p%rsw + j-pj%p%polej)))*raddeg-deg90,rkx)
    lon = real((i-pj%p%polei)*pj%p%dlon*raddeg + pj%p%stdlon,rkx)
    if ( lon >  180.0_rkx ) lon = lon - 360.0_rkx
    if ( lon < -180.0_rkx ) lon = lon + 360.0_rkx
  end subroutine ijll_mc

  subroutine setup_rmc(pj,clon,clat,ci,cj,ds,plon,plat,luvrot)
    implicit none
    type(regcm_projection) , intent(inout) :: pj
    real(rk8) , intent(in) :: clat , clon , cj , ci , ds , plon , plat
    logical , intent(in) :: luvrot
    real(rk8) :: plam , pphi , zphipol
    real(rkx) :: lat , lon , ri , rj
    integer(ik4) :: i , j
    pj%p%dlon = ds*raddeg/earthrad
    pj%p%dlat = ds*raddeg/earthrad
    pj%p%xoff = clon - plon
    pj%p%yoff = clat - plat
    pj%p%polei = ci
    pj%p%polej = cj
    pphi = deg90 - plat
    plam = plon + deg180
    if ( plam > deg180 ) plam = plam - deg360
    pj%p%zlampol = degrad*plam
    zphipol = degrad*pphi
    pj%p%zsinpol = sin(zphipol)
    pj%p%zcospol = cos(zphipol)
    if ( luvrot ) then
      if ( plat > 0.0_rk8 ) then
        pphi = (deg90 - plat) * degrad
        pj%p%pollam = plon + deg180
      else
        pphi = (deg90 + plat) * degrad
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
      call getmem2d(pj%f1,1,pj%p%nlon,1,pj%p%nlat,'projections:f1')
      call getmem2d(pj%f2,1,pj%p%nlon,1,pj%p%nlat,'projections:f2')
      call getmem2d(pj%f3,1,pj%p%nlon,1,pj%p%nlat,'projections:f3')
      do j = 1 , pj%p%nlat
        rj = j
        do i = 1 , pj%p%nlon
          ri = i
          call ijll_rc(pj,ri,rj,lat,lon)
          call uvrot_rc(pj,lat,lon,pj%f1(i,j),pj%f2(i,j))
          pj%f3(i,j) = acos(pj%f1(i,j))
        end do
      end do
    end if
  end subroutine setup_rmc

  pure subroutine llij_rc(pj,lat,lon,i,j)
    implicit none
    class(regcm_projection) , intent(in) :: pj
    real(rkx) , intent(in) :: lat , lon
    real(rkx) , intent(out) :: i , j
    real(rk8) :: zarg , zarg1 , zarg2 , zlam , zphi
    real(rk8) :: lams , phis

    zphi = degrad*lat
    zlam = lon
    if ( zlam > deg180 ) zlam = zlam - deg360
    zlam = degrad*zlam
    zarg = pj%p%zcospol*cos(zphi)*cos(zlam-pj%p%zlampol) + &
           pj%p%zsinpol*sin(zphi)
    phis = asin(zarg)
    phis = log(tan(phis*0.5_rk8+atan(1.0_rk8)))*raddeg
    zarg1 = -sin(zlam-pj%p%zlampol)*cos(zphi)
    zarg2 = -pj%p%zsinpol*cos(zphi)*cos(zlam-pj%p%zlampol) + &
            pj%p%zcospol*sin(zphi)
    if ( abs(zarg2) >= dlowval ) then
      lams = raddeg*atan2(zarg1,zarg2)
    else if ( abs(zarg1) < dlowval ) then
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
    class(regcm_projection) , intent(in) :: pj
    real(rkx) , intent(in) :: i , j
    real(rkx) , intent(out) :: lat , lon
    real(rk8) :: xr , yr , arg , zarg1 , zarg2

    xr = pj%p%xoff + (i-pj%p%polei)*pj%p%dlon
    if ( xr > deg180 ) xr = xr - deg360
    xr = degrad*xr
    yr = pj%p%yoff + (j-pj%p%polej)*pj%p%dlat
    yr = 2.0_rk8*atan(exp(degrad*yr)) - atan(1.0_rk8)*2.0_rk8
    arg = pj%p%zcospol*cos(yr)*cos(xr) + pj%p%zsinpol*sin(yr)
    lat = real(raddeg*asin(arg),rkx)
    zarg1 = sin(pj%p%zlampol)*(-pj%p%zsinpol*cos(xr)*cos(yr)+ &
            pj%p%zcospol*sin(yr))-cos(pj%p%zlampol)*sin(xr)*cos(yr)
    zarg2 = cos(pj%p%zlampol)*(-pj%p%zsinpol*cos(xr)*cos(yr)+ &
            pj%p%zcospol*sin(yr))+sin(pj%p%zlampol)*sin(xr)*cos(yr)
    if ( abs(zarg2) >= dlowval ) then
      lon = real(raddeg*atan2(zarg1,zarg2),rkx)
      if ( lon >  180.0_rkx ) lon = lon - 360.0_rkx
      if ( lon < -180.0_rkx ) lon = lon + 360.0_rkx
    else if ( abs(zarg1) < dlowval ) then
      lon = deg00
    else if ( zarg1 > 0.0_rk8 ) then
      lon = deg90
    else
      lon = -deg90
    end if
  end subroutine ijll_rc

  real(rkx) function rounder(xval,ltop)
    implicit none
    real(rkx) , intent(in) :: xval
    logical , intent(in) :: ltop
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
    real(rk8) , intent(in) :: plat , plon
    real(rk8) , intent(out) :: elat , elon

    if ( plon < 0.0_rk8 ) then
      elon = plon + mathpi
    else if ( plon > 0.0_rk8 ) then
      elon = plon - mathpi
    else
      elon = 0.0_rk8
    end if
    if ( plat > 0.0_rk8 ) then
      elat = halfpi-plat
    else if ( plat < 0.0_rk8 ) then
      elat = -(halfpi+plat)
    else
      elat = 0.0_rk8
    end if
  end subroutine get_equator

  pure elemental real(rkx) function fac_rl(pj,xlat,xlon) result(xmap)
    implicit none
    type(regcm_projection) , intent(in) :: pj
    real(rkx) , intent(in) :: xlat , xlon
    real(rkx) :: ri , rj
    call llij_rl(pj,xlat,xlon,ri,rj)
    xmap = real(d_one/cos(degrad*(pj%p%rlat0+(rj-1)*pj%p%dlat)),rkx)
  end function fac_rl

  pure elemental real(rkx) function fac_lc(pj,lat) result(xmap)
    implicit none
    type(regcm_projection) , intent(in) :: pj
    real(rkx) , intent(in) :: lat
    real(rk8) :: colat
    colat = degrad*(deg90-lat)
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
    type(regcm_projection) , intent(in) :: pj
    real(rkx) , intent(in) :: lat
    xmap = real(pj%p%scale_top/(1.0_rk8 + pj%p%hemi * sin(lat*degrad)),rkx)
  end function fac_ps

  pure elemental real(rkx) function fac_mc(pj,lat) result(xmap)
    implicit none
    type(regcm_projection) , intent(in) :: pj
    real(rkx) , intent(in) :: lat
    xmap = real(1.0_rk8/cos(lat*degrad),rkx)
  end function fac_mc

  pure elemental real(rkx) function fac_rc(pj,xlat,xlon) result(xmap)
    implicit none
    type(regcm_projection) , intent(in) :: pj
    real(rkx) , intent(in) :: xlat , xlon
    real(rkx) :: ri , rj
    real(rk8) :: yr
    call llij_rc(pj,xlat,xlon,ri,rj)
    yr = pj%p%yoff + (rj-pj%p%polej)*pj%p%dlon
    xmap = real(1.0_rk8/cos(yr*degrad),rkx)
  end function fac_rc

  pure elemental real(rk8) function uvrot_lc(pj,lon) result(alpha)
    implicit none
    type(regcm_projection) , intent(in) :: pj
    real(rkx) , intent(in) :: lon
    real(rk8) :: deltalon
    deltalon = pj%p%stdlon - lon
    if ( deltalon > +deg180 ) deltalon = deltalon - deg360
    if ( deltalon < -deg180 ) deltalon = deltalon + deg360
    alpha = deltalon*degrad*pj%p%conefac
  end function uvrot_lc

  pure elemental real(rk8) function uvrot_ps(pj,lon) result(alpha)
    implicit none
    type(regcm_projection) , intent(in) :: pj
    real(rkx) , intent(in) :: lon
    real(rk8) :: deltalon
    deltalon = pj%p%stdlon - lon
    if ( deltalon > +deg180 ) deltalon = deltalon - deg360
    if ( deltalon < -deg180 ) deltalon = deltalon + deg360
    alpha = real(deltalon*degrad*pj%p%hemi,rkx)
  end function uvrot_ps

  pure subroutine uvrot_rc(pj,lat,lon,cosdel,sindel)
    implicit none
    type(regcm_projection) , intent(in) :: pj
    real(rkx) , intent(in) :: lon , lat
    real(rk8) , intent(out) :: cosdel , sindel
    real(rk8) :: zphi , zrla , zrlap , zarg1 , zarg2 , znorm
    zrla = lon
    if ( lat > deg90-0.00001_rk8 ) zrla = 0.0_rk8
    zrlap = (pj%p%pollam - zrla)*degrad
    zphi = lat*degrad
    zarg1 = pj%p%polcphi*sin(zrlap)
    zarg2 = pj%p%polsphi*cos(zphi) - pj%p%polcphi*sin(zphi)*cos(zrlap)
    znorm = 1.0_rk8/sqrt(zarg1*zarg1 + zarg2*zarg2)
    sindel = zarg1*znorm
    cosdel = zarg2*znorm
  end subroutine uvrot_rc

  subroutine rotate2_lc(pj,u,v)
    implicit none
    class(regcm_projection) , intent(in) :: pj
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: u , v
    integer(ik4) :: i1 , i2 , j1 , j2 , i , j
    real(rk8) :: tmp
    i1 = lbound(u,1)
    i2 = ubound(u,1)
    j1 = lbound(u,2)
    j2 = ubound(u,2)
    if ( pj%p%stdlat >= d_zero ) then
      do j = j1 , j2
        do i = i1 , i2
          tmp = u(i,j)*pj%f1(i,j) + v(i,j)*pj%f2(i,j)
          v(i,j) = real(-u(i,j)*pj%f2(i,j) + v(i,j)*pj%f1(i,j),rkx)
          u(i,j) = real(tmp,rkx)
        end do
      end do
    else
      do j = j1 , j2
        do i = i1 , i2
          tmp = u(i,j)*pj%f1(i,j) - v(i,j)*pj%f2(i,j)
          v(i,j) = real(v(i,j)*pj%f1(i,j) + u(i,j)*pj%f2(i,j),rkx)
          u(i,j) = real(tmp,rkx)
        end do
      end do
    end if
  end subroutine rotate2_lc

  subroutine backrotate2_lc(pj,u,v)
    implicit none
    class(regcm_projection) , intent(in) :: pj
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: u , v
    integer(ik4) :: i1 , i2 , j1 , j2 , i , j
    real(rk8) :: tmp
    i1 = lbound(u,1)
    i2 = ubound(u,1)
    j1 = lbound(u,2)
    j2 = ubound(u,2)
    if ( pj%p%stdlat >= d_zero ) then
      do j = j1 , j2
        do i = i1 , i2
          tmp = u(i,j)*pj%f1(i,j) - v(i,j)*pj%f2(i,j)
          v(i,j) = real(u(i,j)*pj%f2(i,j) + v(i,j)*pj%f1(i,j),rkx)
          u(i,j) = real(tmp,rkx)
        end do
      end do
    else
      do j = j1 , j2
        do i = i1 , i2
          tmp = u(i,j)*pj%f1(i,j) + v(i,j)*pj%f2(i,j)
          v(i,j) = real(v(i,j)*pj%f1(i,j) - u(i,j)*pj%f2(i,j),rkx)
          u(i,j) = real(tmp,rkx)
        end do
      end do
    end if
  end subroutine backrotate2_lc

  subroutine rotate3_lc(pj,u,v)
    implicit none
    class(regcm_projection) , intent(in) :: pj
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: u , v
    integer(ik4) :: i1 , i2 , j1 , j2 , k1 , k2 , i , j , k
    real(rk8) :: tmp
    i1 = lbound(u,1)
    i2 = ubound(u,1)
    j1 = lbound(u,2)
    j2 = ubound(u,2)
    k1 = lbound(u,3)
    k2 = ubound(u,3)
    if ( pj%p%stdlat >= d_zero ) then
      do k = k1 , k2
        do j = j1 , j2
          do i = i1 , i2
            tmp = u(i,j,k)*pj%f1(i,j) + v(i,j,k)*pj%f2(i,j)
            v(i,j,k) = real(-u(i,j,k)*pj%f2(i,j) + v(i,j,k)*pj%f1(i,j),rkx)
            u(i,j,k) = real(tmp,rkx)
          end do
        end do
      end do
    else
      do k = k1 , k2
        do j = j1 , j2
          do i = i1 , i2
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
    class(regcm_projection) , intent(in) :: pj
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: u , v
    integer(ik4) :: i1 , i2 , j1 , j2 , k1 , k2 , i , j , k
    real(rk8) :: tmp
    i1 = lbound(u,1)
    i2 = ubound(u,1)
    j1 = lbound(u,2)
    j2 = ubound(u,2)
    k1 = lbound(u,3)
    k2 = ubound(u,3)
    if ( pj%p%stdlat >= d_zero ) then
      do k = k1 , k2
        do j = j1 , j2
          do i = i1 , i2
            tmp = u(i,j,k)*pj%f1(i,j) - v(i,j,k)*pj%f2(i,j)
            v(i,j,k) = real(u(i,j,k)*pj%f2(i,j) + v(i,j,k)*pj%f1(i,j),rkx)
            u(i,j,k) = real(tmp,rkx)
          end do
        end do
      end do
    else
      do k = k1 , k2
        do j = j1 , j2
          do i = i1 , i2
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
    class(regcm_projection) , intent(in) :: pj
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: u , v
    integer(ik4) :: i1 , i2 , j1 , j2 , i , j
    real(rk8) :: tmp
    i1 = lbound(u,1)
    i2 = ubound(u,1)
    j1 = lbound(u,2)
    j2 = ubound(u,2)
    if ( pj%p%stdlat >= d_zero ) then
      do j = j1 , j2
        do i = i1 , i2
          tmp = u(i,j)*pj%f1(i,j) + v(i,j)*pj%f2(i,j)
          v(i,j) = real(-u(i,j)*pj%f2(i,j) + v(i,j)*pj%f1(i,j),rkx)
          u(i,j) = real(tmp,rkx)
        end do
      end do
    else
      do j = j1 , j2
        do i = i1 , i2
          tmp = u(i,j)*pj%f1(i,j) - v(i,j)*pj%f2(i,j)
          v(i,j) = real(v(i,j)*pj%f1(i,j) + u(i,j)*pj%f2(i,j),rkx)
          u(i,j) = real(tmp,rkx)
        end do
      end do
    end if
  end subroutine rotate2_ps

  subroutine backrotate2_ps(pj,u,v)
    implicit none
    class(regcm_projection) , intent(in) :: pj
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: u , v
    integer(ik4) :: i1 , i2 , j1 , j2 , i , j
    real(rk8) :: tmp
    i1 = lbound(u,1)
    i2 = ubound(u,1)
    j1 = lbound(u,2)
    j2 = ubound(u,2)
    if ( pj%p%stdlat >= d_zero ) then
      do j = j1 , j2
        do i = i1 , i2
          tmp = u(i,j)*pj%f1(i,j) - v(i,j)*pj%f2(i,j)
          v(i,j) = real(u(i,j)*pj%f2(i,j) + v(i,j)*pj%f1(i,j),rkx)
          u(i,j) = real(tmp,rkx)
        end do
      end do
    else
      do j = j1 , j2
        do i = i1 , i2
          tmp = u(i,j)*pj%f1(i,j) + v(i,j)*pj%f2(i,j)
          v(i,j) = real(v(i,j)*pj%f1(i,j) - u(i,j)*pj%f2(i,j),rkx)
          u(i,j) = real(tmp,rkx)
        end do
      end do
    end if
  end subroutine backrotate2_ps

  subroutine rotate3_ps(pj,u,v)
    implicit none
    class(regcm_projection) , intent(in) :: pj
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: u , v
    integer(ik4) :: i1 , i2 , j1 , j2 , k1 , k2 , i , j , k
    real(rk8) :: tmp
    i1 = lbound(u,1)
    i2 = ubound(u,1)
    j1 = lbound(u,2)
    j2 = ubound(u,2)
    k1 = lbound(u,3)
    k2 = ubound(u,3)
    if ( pj%p%stdlat >= d_zero ) then
      do k = k1 , k2
        do j = j1 , j2
          do i = i1 , i2
            tmp = u(i,j,k)*pj%f1(i,j) + v(i,j,k)*pj%f2(i,j)
            v(i,j,k) = real(-u(i,j,k)*pj%f2(i,j) + v(i,j,k)*pj%f1(i,j),rkx)
            u(i,j,k) = real(tmp,rkx)
          end do
        end do
      end do
    else
      do k = k1 , k2
        do j = j1 , j2
          do i = i1 , i2
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
    class(regcm_projection) , intent(in) :: pj
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: u , v
    integer(ik4) :: i1 , i2 , j1 , j2 , k1 , k2 , i , j , k
    real(rk8) :: tmp
    i1 = lbound(u,1)
    i2 = ubound(u,1)
    j1 = lbound(u,2)
    j2 = ubound(u,2)
    k1 = lbound(u,3)
    k2 = ubound(u,3)
    if ( pj%p%stdlat >= d_zero ) then
      do k = k1 , k2
        do j = j1 , j2
          do i = i1 , i2
            tmp = u(i,j,k)*pj%f1(i,j) - v(i,j,k)*pj%f2(i,j)
            v(i,j,k) = real(u(i,j,k)*pj%f2(i,j) + v(i,j,k)*pj%f1(i,j),rkx)
            u(i,j,k) = real(tmp,rkx)
          end do
        end do
      end do
    else
      do k = k1 , k2
        do j = j1 , j2
          do i = i1 , i2
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
    class(regcm_projection) , intent(in) :: pj
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: u , v
    return
  end subroutine rotate2_mc

  subroutine rotate3_mc(pj,u,v)
    implicit none
    class(regcm_projection) , intent(in) :: pj
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: u , v
    return
  end subroutine rotate3_mc

  subroutine rotate2_rc(pj,u,v)
    implicit none
    class(regcm_projection) , intent(in) :: pj
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: u , v
    integer(ik4) :: i1 , i2 , j1 , j2 , i , j
    real(rk8) :: tmp
    i1 = lbound(u,1)
    i2 = ubound(u,1)
    j1 = lbound(u,2)
    j2 = ubound(u,2)
    do j = j1 , j2
      do i = i1 , i2
        tmp = u(i,j)*pj%f1(i,j) - v(i,j)*pj%f2(i,j)
        v(i,j) = real(v(i,j)*pj%f1(i,j) + u(i,j)*pj%f2(i,j),rkx)
        u(i,j) = real(tmp,rkx)
      end do
    end do
  end subroutine rotate2_rc

  subroutine backrotate2_rc(pj,u,v)
    implicit none
    class(regcm_projection) , intent(in) :: pj
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: u , v
    integer(ik4) :: i1 , i2 , j1 , j2 , i , j
    real(rk8) :: tmp
    i1 = lbound(u,1)
    i2 = ubound(u,1)
    j1 = lbound(u,2)
    j2 = ubound(u,2)
    do j = j1 , j2
      do i = i1 , i2
        tmp = u(i,j)*pj%f1(i,j) + v(i,j)*pj%f2(i,j)
        v(i,j) = real(v(i,j)*pj%f1(i,j) - u(i,j)*pj%f2(i,j),rkx)
        u(i,j) = real(tmp,rkx)
      end do
    end do
  end subroutine backrotate2_rc

  subroutine rotate3_rc(pj,u,v)
    implicit none
    class(regcm_projection) , intent(in) :: pj
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: u , v
    integer(ik4) :: i1 , i2 , j1 , j2 , k1 , k2 , i , j , k
    real(rk8) :: tmp
    i1 = lbound(u,1)
    i2 = ubound(u,1)
    j1 = lbound(u,2)
    j2 = ubound(u,2)
    k1 = lbound(u,3)
    k2 = ubound(u,3)
    do k = k1 , k2
      do j = j1 , j2
        do i = i1 , i2
          tmp = u(i,j,k)*pj%f1(i,j) - v(i,j,k)*pj%f2(i,j)
          v(i,j,k) = real(v(i,j,k)*pj%f1(i,j) + u(i,j,k)*pj%f2(i,j),rkx)
          u(i,j,k) = real(tmp,rkx)
        end do
      end do
    end do
  end subroutine rotate3_rc

  subroutine backrotate3_rc(pj,u,v)
    implicit none
    class(regcm_projection) , intent(in) :: pj
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: u , v
    integer(ik4) :: i1 , i2 , j1 , j2 , k1 , k2 , i , j , k
    real(rk8) :: tmp
    i1 = lbound(u,1)
    i2 = ubound(u,1)
    j1 = lbound(u,2)
    j2 = ubound(u,2)
    k1 = lbound(u,3)
    k2 = ubound(u,3)
    do k = k1 , k2
      do j = j1 , j2
        do i = i1 , i2
          tmp = u(i,j,k)*pj%f1(i,j) + v(i,j,k)*pj%f2(i,j)
          v(i,j,k) = real(v(i,j,k)*pj%f1(i,j) - u(i,j,k)*pj%f2(i,j),rkx)
          u(i,j,k) = real(tmp,rkx)
        end do
      end do
    end do
  end subroutine backrotate3_rc

  subroutine rotate2_rl(pj,u,v)
    implicit none
    class(regcm_projection) , intent(in) :: pj
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: u , v
    integer(ik4) :: i1 , i2 , j1 , j2 , i , j
    real(rk8) :: tmp1 , tmp2
    if ( pj%p%skiprot ) return
    i1 = lbound(u,1)
    i2 = ubound(u,1)
    j1 = lbound(u,2)
    j2 = ubound(u,2)
    do j = j1 , j2
      do i = i1 , i2
        if ( pj%f2(i,j) /= 0.0_rk8 .and. pj%f3(i,j) /= 0.0_rk8 ) then
          tmp1 = pj%f1(i,j) * u(i,j) + pj%f2(i,j) * v(i,j)
          tmp2 = pj%f3(i,j) * tmp1 + pj%f4(i,j) * v(i,j)
          v(i,j) = real(tmp1,rkx)
          u(i,j) = real(tmp2,rkx)
        else
          v(i,j) = pj%f1(i,j) * v(i,j)
          u(i,j) = pj%f4(i,j) * u(i,j)
        end if
      end do
    end do
  end subroutine rotate2_rl

  subroutine backrotate2_rl(pj,u,v)
    implicit none
    class(regcm_projection) , intent(in) :: pj
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: u , v
    integer(ik4) :: i1 , i2 , j1 , j2 , i , j
    real(rk8) :: tmp1 , tmp2
    if ( pj%p%skiprot ) return
    i1 = lbound(u,1)
    i2 = ubound(u,1)
    j1 = lbound(u,2)
    j2 = ubound(u,2)
    do j = j1 , j2
      do i = i1 , i2
        if ( pj%f8(i,j) /= 0.0_rk8 .and.  pj%f7(i,j) /= 0.0_rk8 ) then
          tmp1 = pj%f5(i,j) * u(i,j) + pj%f6(i,j) * v(i,j)
          tmp2 = pj%f7(i,j) * (v(i,j)-tmp1*pj%f8(i,j))
          u(i,j) = real(tmp2,rkx)
          v(i,j) = real(tmp1,rkx)
        else
          if ( pj%f7(i,j) /= 0.0_rk8 ) then
            u(i,j) = pj%f5(i,j) * u(i,j)
            v(i,j) = pj%f7(i,j) * v(i,j)
          else
            tmp1 = sqrt(u(i,j)*u(i,j)+v(i,j)*v(i,j))
            tmp2 = halfpi + atan2(v(i,j),u(i,j))
            if ( tmp2 > mathpi ) tmp2 = tmp2 - twopi
            if ( tmp2 < mathpi ) tmp2 = tmp2 + twopi
            u(i,j) = -tmp1 * sin(pj%f5(i,j)-tmp2)
            v(i,j) = -tmp1 * cos(pj%f5(i,j)-tmp2)
          end if
        end if
      end do
    end do
  end subroutine backrotate2_rl

  subroutine rotate3_rl(pj,u,v)
    implicit none
    class(regcm_projection) , intent(in) :: pj
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: u , v
    integer(ik4) :: i1 , i2 , j1 , j2 , k1 , k2 , i , j , k
    real(rk8) :: tmp1 , tmp2
    if ( pj%p%skiprot ) return
    i1 = lbound(u,1)
    i2 = ubound(u,1)
    j1 = lbound(u,2)
    j2 = ubound(u,2)
    k1 = lbound(u,3)
    k2 = ubound(u,3)
    do k = k1 , k2
      do j = j1 , j2
        do i = i1 , i2
          if ( pj%f2(i,j) /= 0.0_rk8 .and. pj%f3(i,j) /= 0.0_rk8 ) then
            tmp1 = pj%f1(i,j) * u(i,j,k) + pj%f2(i,j) * v(i,j,k)
            tmp2 = pj%f3(i,j) * tmp1 + pj%f4(i,j) * v(i,j,k)
            v(i,j,k) = real(tmp1,rkx)
            u(i,j,k) = real(tmp2,rkx)
          else
            v(i,j,k) = pj%f1(i,j) * v(i,j,k)
            u(i,j,k) = pj%f4(i,j) * u(i,j,k)
          end if
        end do
      end do
    end do
  end subroutine rotate3_rl

  subroutine backrotate3_rl(pj,u,v)
    implicit none
    class(regcm_projection) , intent(in) :: pj
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: u , v
    integer(ik4) :: i1 , i2 , j1 , j2 , k1 , k2 , i , j , k
    real(rk8) :: tmp1 , tmp2
    if ( pj%p%skiprot ) return
    i1 = lbound(u,1)
    i2 = ubound(u,1)
    j1 = lbound(u,2)
    j2 = ubound(u,2)
    k1 = lbound(u,3)
    k2 = ubound(u,3)
    do k = k1 , k2
      do j = j1 , j2
        do i = i1 , i2
          if ( pj%f8(i,j) /= 0.0_rk8 .and.  pj%f7(i,j) /= 0.0_rk8 ) then
            tmp1 = pj%f5(i,j) * u(i,j,k) + pj%f6(i,j) * v(i,j,k)
            tmp2 = pj%f7(i,j) * (v(i,j,k) - tmp1*pj%f8(i,j))
            u(i,j,k) = real(tmp2,rkx)
            v(i,j,k) = real(tmp1,rkx)
          else
            if ( pj%f7(i,j) /= 0.0_rk8 ) then
              u(i,j,k) = pj%f5(i,j) * u(i,j,k)
              v(i,j,k) = pj%f7(i,j) * v(i,j,k)
            else
              tmp1 = sqrt(u(i,j,k)*u(i,j,k)+v(i,j,k)*v(i,j,k))
              tmp2 = halfpi + atan2(v(i,j,k),u(i,j,k))
              if ( tmp2 > mathpi ) tmp2 = tmp2 - twopi
              if ( tmp2 < mathpi ) tmp2 = tmp2 + twopi
              u(i,j,k) = -tmp1 * sin(pj%f5(i,j)-tmp2)
              v(i,j,k) = -tmp1 * cos(pj%f5(i,j)-tmp2)
            end if
          end if
        end do
      end do
    end do
  end subroutine backrotate3_rl

  pure subroutine mapfac_rl(pj,xlat,xlon,xmap)
    implicit none
    class(regcm_projection) , intent(in) :: pj
    real(rkx) , pointer , dimension(:,:) , intent(in) :: xlat , xlon
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: xmap
    xmap = fac_rl(pj,xlat,xlon)
  end subroutine mapfac_rl

  pure subroutine mapfac_lc(pj,xlat,xlon,xmap)
    implicit none
    class(regcm_projection) , intent(in) :: pj
    real(rkx) , pointer , dimension(:,:) , intent(in) :: xlat , xlon
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: xmap
    xmap = fac_lc(pj,xlat)
  end subroutine mapfac_lc

  pure subroutine mapfac_ps(pj,xlat,xlon,xmap)
    implicit none
    class(regcm_projection) , intent(in) :: pj
    real(rkx) , pointer , dimension(:,:) , intent(in) :: xlat , xlon
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: xmap
    xmap = fac_ps(pj,xlat)
  end subroutine mapfac_ps

  pure subroutine mapfac_mc(pj,xlat,xlon,xmap)
    implicit none
    class(regcm_projection) , intent(in) :: pj
    real(rkx) , pointer , dimension(:,:) , intent(in) :: xlat , xlon
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: xmap
    xmap = fac_mc(pj,xlat)
  end subroutine mapfac_mc

  pure subroutine mapfac_ll(pj,xlat,xlon,xmap)
    implicit none
    class(regcm_projection) , intent(in) :: pj
    real(rkx) , pointer , dimension(:,:) , intent(in) :: xlat , xlon
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: xmap
    xmap = d_one
  end subroutine mapfac_ll

  pure subroutine mapfac_rc(pj,xlat,xlon,xmap)
    implicit none
    class(regcm_projection) , intent(in) :: pj
    real(rkx) , pointer , dimension(:,:) , intent(in) :: xlat , xlon
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: xmap
    xmap = fac_rc(pj,xlat,xlon)
  end subroutine mapfac_rc

end module mod_projections
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
