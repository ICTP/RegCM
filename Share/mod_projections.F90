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

  type regcm_projection
    private
    real(rk8) :: stdlon , stdlat
    real(rk8) :: truelat1 , truelat2 , tl2r , ctl1r
    real(rk8) :: colat1 , colat2 , nfac
    real(rk8) :: rsw , rebydx , hemi
    real(rk8) :: reflon , dlon , dlat , scale_top
    real(rk8) :: polei , polej
    real(rk8) :: rlon0 , rlat0 , phi0 , lam0
    real(rk8) :: xoff , yoff
    real(rk8) :: zsinpol , zcospol , zlampol
    real(rk8) :: conefac
    logical :: lamtan
    integer(ik4) :: nlat , nlon
    procedure(transform) , pass(pj) , pointer , public :: llij => NULL()
    procedure(transform) , pass(pj) , pointer , public :: ijll => NULL()
    procedure(map_factor) , pass(pj) , pointer , public :: mapfac => NULL()
    procedure(rotate3) , pass(pj) , pointer , public :: uvrotate3 => NULL()
    procedure(rotate2) , pass(pj) , pointer , public :: uvrotate2 => NULL()
    procedure(rotate2) , pass(pj) , pointer , public :: uvbkrotate2 => NULL()
    procedure(rotate3) , pass(pj) , pointer , public :: uvbkrotate3 => NULL()
    real(rk8) , pointer , dimension(:,:) :: f1 , f2 , f3 , f4
    real(rk8) , pointer , dimension(:,:) :: f5 , f6 , f7 , f8
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
    pure subroutine rotate3(pj,u,v,ur,vr)
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
      real(rkx) , pointer , dimension(:,:,:) , intent(out) :: ur , vr
    end subroutine rotate3
  end interface

  abstract interface
    pure subroutine rotate2(pj,u,v,ur,vr)
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
      real(rkx) , pointer , dimension(:,:) , intent(out) :: ur , vr
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
      real(rkx) , pointer , dimension(:,:) , intent(out) :: c
    end subroutine map_factor
  end interface

  public :: anyprojparams , regcm_projection , init_projection
  public :: uvrot_inplace , uvantirot_inplace

  contains

  subroutine uvrot_inplace(pj,u,v)
    implicit none
    type(regcm_projection) , intent(in) :: pj
    real(rkx) , dimension(:,:,:) , pointer , intent(inout) :: u , v
    real(rkx) , dimension(:,:,:) , pointer :: ur , vr
    call assignpnt(u,ur)
    call assignpnt(v,vr)
    call pj%uvrotate3(u,v,ur,vr)
  end subroutine uvrot_inplace

  subroutine uvantirot_inplace(pj,u,v)
    implicit none
    type(regcm_projection) , intent(in) :: pj
    real(rkx) , dimension(:,:,:) , pointer , intent(inout) :: u , v
    real(rkx) , dimension(:,:,:) , pointer :: ur , vr
    call assignpnt(u,ur)
    call assignpnt(v,vr)
    call pj%uvbkrotate3(u,v,ur,vr)
  end subroutine uvantirot_inplace

  subroutine init_projection(pjpara,pj)
    implicit none
    type(anyprojparams) , intent(in) :: pjpara
    type(regcm_projection) , intent(out) :: pj
    real(rk8) :: ci , cj
    ci = real(pjpara%nlon,rk8) * 0.5_rk8
    cj = real(pjpara%nlat,rk8) * 0.5_rk8
    if ( .not. pjpara%staggerx ) ci = ci - 0.5_rk8
    if ( .not. pjpara%staggery ) cj = cj - 0.5_rk8
    pj%nlon = pjpara%nlon
    pj%nlat = pjpara%nlat
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
  end subroutine init_projection

  subroutine setup_ll(pj,clon,clat,ci,cj,ds)
    implicit none
    type(regcm_projection) , intent(inout) :: pj
    real(rk8) , intent(in) :: ci , cj , clat , clon , ds
    real(rk8) :: phi , lam
    pj%dlon = raddeg * ds / earthrad
    pj%dlat = pj%dlon
    pj%rlat0 = clat - cj*pj%dlat
    pj%rlon0 = clon - ci*pj%dlon
    if ( pj%rlat0 >  deg90 )  pj%rlat0 = deg90 - pj%rlat0
    if ( pj%rlat0 < -deg90 )  pj%rlat0 = pj%rlat0 + deg90
    if ( pj%rlon0 >  deg180 ) pj%rlon0 = deg360 - pj%rlon0
    if ( pj%rlon0 < -deg180 ) pj%rlon0 = pj%rlon0 + deg360
  end subroutine setup_ll

  pure subroutine ijll_ll(pj,i,j,lat,lon)
    implicit none
    class(regcm_projection) , intent(in) :: pj
    real(rkx) , intent(in) :: i , j
    real(rkx) , intent(out) :: lat , lon
    real(rk8) :: phi , lam
    lat = pj%rlat0 + (j-1.0_rk8) * pj%dlat
    lon = pj%rlon0 + (i-1.0_rk8) * pj%dlon
    if ( lat >  90.0_rkx ) lat = 90.0_rkx - lat
    if ( lat < -90.0_rkx ) lat = lat + 90.0_rkx
    if ( lon >  180.0_rkx ) lon = 360.0_rkx - lon
    if ( lon < -180.0_rkx ) lon = lon + 360.0_rkx
  end subroutine ijll_ll

  pure subroutine llij_ll(pj,lat,lon,i,j)
    implicit none
    class(regcm_projection) , intent(in) :: pj
    real(rkx) , intent(in) :: lat , lon
    real(rkx) , intent(out) :: i , j
    i = real(( lon - pj%rlon0 ) / pj%dlon + 1.0_rk8,rkx)
    j = real(( lat - pj%rlat0 ) / pj%dlat + 1.0_rk8,rkx)
  end subroutine llij_ll

  subroutine setup_rll(pj,clon,clat,ci,cj,ds,plon,plat,luvrot)
    implicit none
    type(regcm_projection) , intent(inout) :: pj
    real(rk8) , intent(in) :: ci , cj , clat , clon , plon , plat , ds
    logical , intent(in) :: luvrot
    real(rk8) :: phi , lam , dlam , lon , lat
    real(rk8) :: rotlam , rotphi , ri , rj
    integer(ik4) :: i , j
    pj%dlon = raddeg * ds / earthrad
    pj%dlat = pj%dlon
    phi = degrad*plat
    lam = degrad*plon
    call get_equator(phi,lam,pj%phi0,pj%lam0)
    phi = clat
    lam = clon
    if ( phi >  deg90 )  phi = deg90 - phi
    if ( phi < -deg90 )  phi = phi + deg90
    if ( lam >  deg180 ) lam = deg360 - lam
    if ( lam < -deg180 ) lam = lam + deg360
    phi = degrad * phi
    lam = degrad * lam
    pj%rlat0 = asin(-cos(phi)*sin(pj%phi0)*cos(lam-pj%lam0) + &
                    sin(phi)*cos(pj%phi0))
    if ( abs(abs(pj%rlat0)-halfpi) > 1.0e-7_rk8 .and. &
         abs(pj%phi0) > 1.0e-7_rk8 ) then
      pj%rlon0 = (sin(phi)-cos(pj%phi0)*sin(pj%rlat0)) / &
                 (sin(pj%phi0)*cos(pj%rlat0))
      if ( pj%rlon0 < -1.0_rk8 .and. &
           pj%rlon0 > -1.00001_rk8 ) pj%rlon0 = -1.0_rk8
      if ( pj%rlon0 >  1.0_rk8 .and. &
           pj%rlon0 <  1.00001_rk8 ) pj%rlon0 =  1.0_rk8
      pj%rlon0 = acos(pj%rlon0)
      if ( lam < pj%lam0 ) then
        pj%rlon0 = -pj%rlon0
      end if
    else
      pj%rlon0 = lam
    end if
    pj%rlat0 = raddeg*pj%rlat0 - cj*pj%dlat
    pj%rlon0 = raddeg*pj%rlon0 - ci*pj%dlon
    if ( pj%rlat0 >  deg90 )  pj%rlat0 = deg90 - pj%rlat0
    if ( pj%rlat0 < -deg90 )  pj%rlat0 = pj%rlat0 + deg90
    if ( pj%rlon0 >  deg180 ) pj%rlon0 = deg360 - pj%rlon0
    if ( pj%rlon0 < -deg180 ) pj%rlon0 = pj%rlon0 + deg360
    if ( luvrot ) then
      call getmem2d(pj%f1,1,pj%nlon,1,pj%nlat,'projections:f1')
      call getmem2d(pj%f2,1,pj%nlon,1,pj%nlat,'projections:f2')
      call getmem2d(pj%f3,1,pj%nlon,1,pj%nlat,'projections:f3')
      call getmem2d(pj%f4,1,pj%nlon,1,pj%nlat,'projections:f4')
      call getmem2d(pj%f5,1,pj%nlon,1,pj%nlat,'projections:f5')
      call getmem2d(pj%f6,1,pj%nlon,1,pj%nlat,'projections:f6')
      call getmem2d(pj%f7,1,pj%nlon,1,pj%nlat,'projections:f7')
      call getmem2d(pj%f8,1,pj%nlon,1,pj%nlat,'projections:f8')
      do j = 1 , pj%nlat
        rotphi = degrad*(pj%rlat0 + (j-1)*pj%dlat)
        do i = 1 , pj%nlon
          rotlam = degrad*(pj%rlon0 + (i-1)*pj%dlon)
          ri = i
          rj = j
          call ijll_rl(pj,ri,rj,lat,lon)
          lam = degrad*lon
          phi = degrad*lat
          dlam = lam - pj%lam0
          if ( abs(rotlam) < 1.0e-4 ) then
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
            pj%f1(i,j) = sin(dlam)*sin(pj%phi0) / cos(rotphi)
            pj%f2(i,j) = (cos(dlam)*sin(phi)*sin(pj%phi0) + &
                          cos(pj%phi0)*cos(phi)) / &
                          cos(rotphi)
            pj%f3(i,j) = (cos(pj%phi0)*cos(rotphi) - &
                          sin(pj%phi0)*sin(rotphi)*cos(rotlam)) / &
                         (sin(pj%phi0)*sin(rotlam))
            pj%f4(i,j) = -cos(phi) / (sin(pj%phi0)*sin(rotlam))
          end if
          if ( abs(cos(phi)) > 1.0e-2_rk8 ) then
            if ( abs(sin(dlam)) > 1.0e-1_rk8 ) then
              pj%f5(i,j) = -sin(pj%phi0)*sin(rotlam)/cos(phi)
              pj%f6(i,j) = (cos(pj%phi0)*cos(rotphi) - &
                            sin(pj%phi0)*sin(rotphi)*cos(rotlam)) / cos(phi)
              pj%f7(i,j) = cos(rotphi)/sin(dlam)/sin(pj%phi0)
              pj%f8(i,j) = (cos(dlam)*sin(pj%phi0)*sin(phi) + &
                            cos(pj%phi0)*cos(phi))/cos(rotphi)
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
            if ( pj%phi0 < 0.0 ) then
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
        end do
      end do
    end if
  end subroutine setup_rll

  pure subroutine ijll_rl(pj,i,j,lat,lon)
    class(regcm_projection) , intent(in) :: pj
    real(rkx) , intent(in) :: i , j
    real(rkx) , intent(out) :: lat , lon
    real(rk8) :: phi , lam
    phi = pj%rlat0 + (j-1.0_rk8) * pj%dlat
    lam = pj%rlon0 + (i-1.0_rk8) * pj%dlon
    if ( phi >  deg90 )  phi = deg90 - phi
    if ( phi < -deg90 )  phi = phi + deg90
    if ( lam >  deg180 ) lam = deg360 - lam
    if ( lam < -deg180 ) lam = lam + deg360
    phi = degrad * phi
    lam = degrad * lam
    if ( abs(pj%phi0) > 1.0e-7_rk8 ) then
      lat = asin(cos(phi)*sin(pj%phi0)*cos(lam) + sin(phi)*cos(pj%phi0))
      if ( abs(abs(lat)-halfpi) > 1.0e-7_rk8 ) then
        lon = (-sin(phi)*sin(pj%phi0) + cos(pj%phi0)*cos(lam)*cos(phi))/cos(lat)
        if ( lon < -1.0_rk8 .and. lon > -1.00001_rk8 ) lon = -1.0_rk8
        if ( lon >  1.0_rk8 .and. lon <  1.00001_rk8 ) lon =  1.0_rk8
        if ( lam < 0.0_rk8 ) then
          lon = pj%lam0 - acos(lon)
        else
          lon = pj%lam0 + acos(lon)
        end if
      else
        lon = lam
      end if
    else
      lat = asin(cos(phi)*sin(pj%phi0)*cos(lam) + sin(phi)*cos(pj%phi0))
      lon = lam
    end if
    lat = raddeg*lat
    lon = raddeg*lon
    if ( lat >  90.0_rkx ) lat = 90.0_rkx - lat
    if ( lat < -90.0_rkx ) lat = lat + 90.0_rkx
    if ( lon >  180.0_rkx ) lon = 360.0_rkx - lon
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
    if ( lam >  deg180 ) lam = deg360 - lam
    if ( lam < -deg180 ) lam = lam + deg360
    if ( phi >  deg90 )  phi = deg90 - phi
    if ( phi < -deg90 )  phi = phi + deg90
    phi = degrad * phi
    lam = degrad * lam
    rlat = asin(-cos(phi)*sin(pj%phi0)*cos(lam-pj%lam0) + sin(phi)*cos(pj%phi0))
    if ( abs(abs(rlat)-halfpi) > 1.0e-7_rk8 .and. &
         abs(pj%phi0) > 1.0e-7_rk8 ) then
      rlon = (sin(phi)-cos(pj%phi0)*sin(rlat))/(sin(pj%phi0)*cos(rlat))
      if ( rlon < -1.0_rk8 .and. rlon > -1.00001_rk8 ) rlon = -1.0_rk8
      if ( rlon >  1.0_rk8 .and. rlon <  1.00001_rk8 ) rlon =  1.0_rk8
      rlon = acos(rlon)
      if ( lam < pj%lam0 ) then
        rlon = -rlon
      end if
    else
      rlon = lam
    end if
    rlon = raddeg*rlon
    rlat = raddeg*rlat
    i = real(( rlon - pj%rlon0 ) / pj%dlon + 1.0_rk8,rkx)
    j = real(( rlat - pj%rlat0 ) / pj%dlat + 1.0_rk8,rkx)
  end subroutine llij_rl

  subroutine setup_lcc(pj,clon,clat,ci,cj,ds,slon,trlat1,trlat2,luvrot)
    implicit none
    type(regcm_projection) , intent(inout) :: pj
    real(rk8) , intent(in) :: ci , cj , slon , clat , clon , ds , &
                              trlat1 , trlat2
    logical , intent(in) :: luvrot
    real(rk8) :: arg , deltalon1 , tl1r , tl2r , ri , rj , angle , lat , lon
    integer(ik4) :: i , j

    pj%stdlon = slon
    pj%stdlat = clat
    pj%truelat1 = trlat1
    pj%truelat2 = trlat2
    tl1r = pj%truelat1*degrad
    tl2r = pj%truelat2*degrad
    pj%colat1  = degrad*(deg90 - pj%truelat1)
    pj%colat2  = degrad*(deg90 - pj%truelat2)
    pj%nfac = (log(sin(pj%colat1))        - log(sin(pj%colat2))) / &
              (log(tan(pj%colat1*0.5_rk8)) - log(tan(pj%colat2*0.5_rk8)))
    if ( pj%truelat1 > 0.0_rk8 ) then
      pj%hemi =  1.0_rk8
    else
      pj%hemi = -1.0_rk8
    end if
    pj%rebydx = earthrad / ds
    if ( abs(pj%truelat1-pj%truelat2) > 0.1_rk8 ) then
      pj%conefac = log10(cos(tl1r)) - log10(cos(tl2r))
      pj%conefac = pj%conefac / &
                (log10(tan((deg45-abs(pj%truelat1)/2.0_rk8)*degrad)) - &
                 log10(tan((deg45-abs(pj%truelat2)/2.0_rk8)*degrad)))
      pj%lamtan = .false.
    else
      pj%conefac = sin(abs(tl1r))
      pj%lamtan = .true.
    end if
    deltalon1 = clon - pj%stdlon
    if ( deltalon1 >  deg180 ) deltalon1 = deltalon1 - deg360
    if ( deltalon1 < -deg180 ) deltalon1 = deltalon1 + deg360
    pj%ctl1r = cos(tl1r)
    pj%rsw = pj%rebydx * pj%ctl1r/pj%conefac * &
               (tan((deg90*pj%hemi-clat)    *degrad*0.5_rk8) / &
                tan((deg90*pj%hemi-pj%truelat1)*degrad*0.5_rk8))**pj%conefac
    arg = pj%conefac*(deltalon1*degrad)
    pj%polei = pj%hemi*ci - pj%hemi * pj%rsw * sin(arg)
    pj%polej = pj%hemi*cj + pj%rsw * cos(arg)
    if ( luvrot ) then
      call getmem2d(pj%f1,1,pj%nlon,1,pj%nlat,'projections:f1')
      call getmem2d(pj%f2,1,pj%nlon,1,pj%nlat,'projections:f2')
      do j = 1 , pj%nlat
        do i = 1 , pj%nlon
          ri = i
          rj = j
          call ijll_lc(pj,ri,rj,lat,lon)
          angle = uvrot_lc(pj,lon)
          pj%f1(i,j) = cos(angle)
          pj%f2(i,j) = sin(angle)
        end do
      end do
    end if
  end subroutine setup_lcc

  pure subroutine ijll_lc(pj,i,j,lat,lon)
    implicit none
    class(regcm_projection) , intent(in) :: pj
    real(rkx) , intent(in) :: i , j
    real(rkx) , intent(out) :: lat , lon
    real(rk8) :: chi1 , chi2 , chi
    real(rk8) :: inew , jnew , xx , yy , r2 , r

    chi1 = (deg90 - pj%hemi*pj%truelat1)*degrad
    chi2 = (deg90 - pj%hemi*pj%truelat2)*degrad
    inew = pj%hemi * i
    jnew = pj%hemi * j
    xx = inew - pj%polei
    yy = pj%polej - jnew
    r2 = (xx*xx + yy*yy)
    r = sqrt(r2)/pj%rebydx
    if ( abs(r2) < dlowval ) then
      lat = real(pj%hemi * deg90,rkx)
      lon = real(pj%stdlon,rkx)
    else
      lon = real(pj%stdlon + raddeg * atan2(pj%hemi*xx,yy)/pj%conefac,rkx)
      lon = mod(lon+360.0_rkx, 360.0_rkx)
      if ( abs(chi1-chi2) < dlowval ) then
        chi = 2.0_rk8 * &
          atan((r/tan(chi1))**(1.0_rk8/pj%conefac)*tan(chi1*0.5_rk8))
      else
        chi = 2.0_rk8*atan((r*pj%conefac/sin(chi1))**(1.0_rk8/pj%conefac) * &
              tan(chi1*0.5_rk8))
      end if
      lat = real((deg90-chi*raddeg)*pj%hemi,rkx)
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

    deltalon = lon - pj%stdlon
    if ( deltalon > +deg180 ) deltalon = deltalon - deg360
    if ( deltalon < -deg180 ) deltalon = deltalon + deg360
    rm = pj%rebydx * pj%ctl1r/pj%conefac * &
           (tan((deg90*pj%hemi-lat)*degrad*0.5_rk8) / &
            tan((deg90*pj%hemi-pj%truelat1)*degrad*0.5_rk8))**pj%conefac
    arg = pj%conefac*(deltalon*degrad)
    i = real(pj%hemi*(pj%polei + pj%hemi * rm * sin(arg)),rkx)
    j = real(pj%hemi*(pj%polej - rm * cos(arg)),rkx)
  end subroutine llij_lc

  subroutine setup_plr(pj,clon,clat,ci,cj,ds,slon,luvrot)
    implicit none
    type(regcm_projection) , intent(inout) :: pj
    real(rk8) , intent(in) :: clat , clon , cj , ci , ds , slon
    logical , intent(in) :: luvrot
    real(rk8) :: ala1 , alo1 , angle , ri , rj , lat , lon
    integer(ik4) :: i , j

    pj%stdlon = slon
    pj%stdlat = clat
    if ( pj%stdlat > 0.0_rk8 ) then
      pj%hemi = 1.0_rk8
    else
      pj%hemi = -1.0_rk8
    end if
    pj%rebydx = earthrad / ds
    pj%reflon = pj%stdlon + deg90
    ala1 = clat*degrad
    alo1 = (clon-pj%reflon)*degrad
    pj%scale_top = 1.0_rk8 + pj%hemi * sin(ala1)
    pj%rsw = pj%rebydx*cos(ala1)*pj%scale_top/(1.0_rk8+pj%hemi*sin(ala1))
    pj%polei = ci - pj%rsw * cos(alo1)
    pj%polej = cj - pj%hemi * pj%rsw * sin(alo1)
    if ( luvrot ) then
      call getmem2d(pj%f1,1,pj%nlon,1,pj%nlat,'projections:f1')
      call getmem2d(pj%f2,1,pj%nlon,1,pj%nlat,'projections:f2')
      do j = 1 , pj%nlat
        do i = 1 , pj%nlon
          ri = i
          rj = j
          call ijll_ps(pj,ri,rj,lat,lon)
          angle = uvrot_ps(pj,lon)
          pj%f1(i,j) = cos(angle)
          pj%f2(i,j) = sin(angle)
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

    deltalon = lon - pj%reflon
    if ( deltalon > +deg180 ) deltalon = deltalon - deg360
    if ( deltalon < -deg180 ) deltalon = deltalon + deg360
    alo = deltalon * degrad
    ala = lat * degrad
    rm = pj%rebydx * cos(ala) * pj%scale_top/(1.0_rk8 + pj%hemi * sin(ala))
    i = real(pj%polei + rm * cos(alo),rkx)
    j = real(pj%polej + pj%hemi * rm * sin(alo),rkx)
  end subroutine llij_ps

  pure subroutine ijll_ps(pj,i,j,lat,lon)
    implicit none
    class(regcm_projection) , intent(in) :: pj
    real(rkx) , intent(in) :: i , j
    real(rkx) , intent(out) :: lat , lon
    real(rk8) :: xx , yy , r2 , gi2 , arcc

    xx = i - pj%polei
    yy = (j - pj%polej) * pj%hemi
    r2 = xx**2.0_rk8 + yy**2.0_rk8
    if ( abs(r2) < dlowval ) then
      lat = real(pj%hemi*deg90,rkx)
      lon = real(pj%reflon,rkx)
    else
      gi2 = (pj%rebydx * pj%scale_top)**2
      lat = real(raddeg * pj%hemi * asin((gi2-r2)/(gi2+r2)),rkx)
      arcc = acos(xx/sqrt(r2))
      if ( yy > 0.0_rk8 ) then
        lon = real(pj%reflon + raddeg * arcc,rkx)
      else
        lon = real(pj%reflon - raddeg * arcc,rkx)
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

    pj%stdlon = clon
    clain = cos(clat*degrad)
    pj%dlon = ds / (earthrad * clain)
    pj%rsw = 0.0_rk8
    if ( abs(clat) > dlowval ) then
      pj%rsw = (log(tan(0.5_rk8*((clat+deg90)*degrad))))/pj%dlon
    end if
    pj%polei = ci
    pj%polej = cj
  end subroutine setup_mrc

  pure subroutine llij_mc(pj,lat,lon,i,j)
    implicit none
    class(regcm_projection) , intent(in) :: pj
    real(rkx) , intent(in) :: lat , lon
    real(rkx) , intent(out) :: i , j
    real(rk8) :: deltalon

    deltalon = lon - pj%stdlon
    if ( deltalon > +deg180 ) deltalon = deltalon - deg360
    if ( deltalon < -deg180 ) deltalon = deltalon + deg360
    i = real(pj%polei + (deltalon/(pj%dlon*raddeg)),rkx)
    j = real(pj%polej + &
         (log(tan(0.5_rk8*((lat+deg90)*degrad)))) / pj%dlon - pj%rsw,rkx)
  end subroutine llij_mc

  pure subroutine ijll_mc(pj,i,j,lat,lon)
    implicit none
    class(regcm_projection) , intent(in) :: pj
    real(rkx) , intent(in) :: i , j
    real(rkx) , intent(out) :: lat , lon

    lat = real(2.0_rk8 * &
        atan(exp(pj%dlon*(pj%rsw + j-pj%polej)))*raddeg-deg90,rkx)
    lon = real((i-pj%polei)*pj%dlon*raddeg + pj%stdlon,rkx)
    if ( lon >  180.0_rkx ) lon = lon - 360.0_rkx
    if ( lon < -180.0_rkx ) lon = lon + 360.0_rkx
  end subroutine ijll_mc

  subroutine setup_rmc(pj,clon,clat,ci,cj,ds,plon,plat,luvrot)
    implicit none
    type(regcm_projection) , intent(inout) :: pj
    real(rk8) , intent(in) :: clat , clon , cj , ci , ds , plon , plat
    logical , intent(in) :: luvrot
    real(rk8) :: plam , pphi , zphipol , angle , ri , rj , lat , lon
    integer(ik4) :: i , j
    pj%dlon = ds*raddeg/earthrad
    pj%xoff = clon - plon
    pj%yoff = clat - plat
    pj%polei = ci
    pj%polej = cj
    pphi = deg90 - plat
    plam = plon + deg180
    if ( plam>deg180 ) plam = plam - deg360
    pj%zlampol = degrad*plam
    zphipol = degrad*pphi
    pj%zsinpol = sin(zphipol)
    pj%zcospol = cos(zphipol)
    if ( luvrot ) then
      call getmem2d(pj%f1,1,pj%nlon,1,pj%nlat,'projections:f1')
      call getmem2d(pj%f2,1,pj%nlon,1,pj%nlat,'projections:f2')
      do j = 1 , pj%nlat
        do i = 1 , pj%nlon
          ri = i
          rj = j
          call ijll_rc(pj,ri,rj,lat,lon)
          angle = uvrot_rc(pj,lon,lat)
          pj%f1(i,j) = cos(angle)
          pj%f2(i,j) = sin(angle)
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
    if ( zlam>deg180 ) zlam = zlam - deg360
    zlam = degrad*zlam
    zarg = pj%zcospol*cos(zphi)*cos(zlam-pj%zlampol) + pj%zsinpol*sin(zphi)
    phis = asin(zarg)
    phis = log(tan(phis*0.5_rk8+atan(1.0_rk8)))*raddeg
    zarg1 = -sin(zlam-pj%zlampol)*cos(zphi)
    zarg2 = -pj%zsinpol*cos(zphi)*cos(zlam-pj%zlampol) + pj%zcospol*sin(zphi)
    if ( abs(zarg2) >= dlowval ) then
      lams = raddeg*atan2(zarg1,zarg2)
    else if ( abs(zarg1) < dlowval ) then
      lams = deg00
    else if ( zarg1 > 0.0_rk8 ) then
      lams = deg90
    else
      lams = -deg90
    end if
    i = real(pj%polei + (lams-pj%xoff)/pj%dlon,rkx)
    j = real(pj%polej + (phis-pj%yoff)/pj%dlon,rkx)
  end subroutine llij_rc

  pure subroutine ijll_rc(pj,i,j,lat,lon)
    implicit none
    class(regcm_projection) , intent(in) :: pj
    real(rkx) , intent(in) :: i , j
    real(rkx) , intent(out) :: lat , lon
    real(rk8) :: xr , yr , arg , zarg1 , zarg2

    xr = pj%xoff + (i-pj%polei)*pj%dlon
    if ( xr > deg180 ) xr = xr - deg360
    xr = degrad*xr
    yr = pj%yoff + (j-pj%polej)*pj%dlon
    yr = 2.0_rk8*atan(exp(degrad*yr)) - atan(1.0_rk8)*2.0_rk8
    arg = pj%zcospol*cos(yr)*cos(xr) + pj%zsinpol*sin(yr)
    lat = real(raddeg*asin(arg),rkx)
    zarg1 = sin(pj%zlampol)*(-pj%zsinpol*cos(xr)*cos(yr)+ &
            pj%zcospol*sin(yr))-cos(pj%zlampol)*sin(xr)*cos(yr)
    zarg2 = cos(pj%zlampol)*(-pj%zsinpol*cos(xr)*cos(yr)+ &
            pj%zcospol*sin(yr))+sin(pj%zlampol)*sin(xr)*cos(yr)
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
      tmpval = ceiling(xval*100.0_rkx)
    else
      tmpval = floor(xval*100.0_rkx)
    end if
    rounder = tmpval/100.0_rkx
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
    real(rk8) :: ri , rj
    call llij_rl(pj,xlat,xlon,ri,rj)
    xmap = real(d_one/cos(degrad*(pj%rlat0+(rj-1)*pj%dlat)),rkx)
  end function fac_rl

  pure elemental real(rkx) function fac_lc(pj,lat) result(xmap)
    implicit none
    type(regcm_projection) , intent(in) :: pj
    real(rkx) , intent(in) :: lat
    real(rk8) :: colat
    colat = degrad*(deg90-lat)
    if ( .not. pj%lamtan ) then
      xmap = real(sin(pj%colat2)/sin(colat) * &
             (tan(colat*0.5_rk8)/tan(pj%colat2*0.5_rk8))**pj%nfac,rkx)
    else
      xmap = real(sin(pj%colat1)/sin(colat) * &
             (tan(colat*0.5_rk8)/tan(pj%colat1*0.5_rk8))**cos(pj%colat1),rkx)
    endif
  end function fac_lc

  pure elemental real(rkx) function fac_ps(pj,lat) result(xmap)
    implicit none
    type(regcm_projection) , intent(in) :: pj
    real(rkx) , intent(in) :: lat
    xmap = real(pj%scale_top/(1.0_rk8 + pj%hemi * sin(lat*degrad)),rkx)
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
    real(rk8) :: ri , rj , yr
    call llij_rc(pj,xlat,xlon,ri,rj)
    yr = pj%yoff + (rj-pj%polej)*pj%dlon
    xmap = real(1.0_rk8/cos(yr*degrad),rkx)
  end function fac_rc

  pure elemental real(rk8) function uvrot_lc(pj,lon) result(alpha)
    implicit none
    type(regcm_projection) , intent(in) :: pj
    real(rkx) , intent(in) :: lon
    real(rk8) :: deltalon
    deltalon = pj%stdlon - lon
    if ( deltalon > +deg180 ) deltalon = deltalon - deg360
    if ( deltalon < -deg180 ) deltalon = deltalon + deg360
    alpha = deltalon*degrad*pj%conefac
  end function uvrot_lc

  pure elemental real(rk8) function uvrot_ps(pj,lon) result(alpha)
    implicit none
    type(regcm_projection) , intent(in) :: pj
    real(rkx) , intent(in) :: lon
    real(rk8) :: deltalon
    deltalon = pj%stdlon - lon
    if ( deltalon > +deg180 ) deltalon = deltalon - deg360
    if ( deltalon < -deg180 ) deltalon = deltalon + deg360
    alpha = real(deltalon*degrad*pj%hemi,rkx)
  end function uvrot_ps

  pure elemental real(rk8) function uvrot_rc(pj,lat,lon) result(alpha)
    implicit none
    type(regcm_projection) , intent(in) :: pj
    real(rkx) , intent(in) :: lon , lat
    real(rk8) :: zphi , zrla , zrlap , zarg1 , zarg2 , znorm
    zphi = lat*degrad
    zrla = lon*degrad
    if ( lat > deg90-0.00001_rk8 ) zrla = 0.0_rk8
    zrlap = pj%zlampol - zrla
    zarg1 = pj%zcospol*sin(zrlap)
    zarg2 = pj%zsinpol*cos(zphi) - pj%zcospol*sin(zphi)*cos(zrlap)
    znorm = 1.0_rk8/sqrt(zarg1**2.0_rk8+zarg2**2.0_rk8)
    alpha = real(acos(zarg2*znorm),rkx)
  end function uvrot_rc

  pure subroutine rotate2_lc(pj,u,v,ur,vr)
    implicit none
    class(regcm_projection) , intent(in) :: pj
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: u , v
    real(rkx) , pointer , dimension(:,:) , intent(out) :: ur , vr
    integer(ik4) :: i1 , i2 , j1 , j2 , i , j
    real(rk8) :: tmp
    i1 = lbound(u,1)
    i2 = ubound(u,1)
    j1 = lbound(u,2)
    j2 = ubound(u,2)
    if ( pj%stdlat >= d_zero ) then
      do j = j1 , j2
        do i = i1 , i2
          tmp = u(i,j)*pj%f1(i,j) + v(i,j)*pj%f2(i,j)
          vr(i,j) = real(-u(i,j)*pj%f2(i,j) + v(i,j)*pj%f1(i,j),rkx)
          ur(i,j) = real(tmp,rkx)
        end do
      end do
    else
      do j = j1 , j2
        do i = i1 , i2
          tmp = u(i,j)*pj%f1(i,j) - v(i,j)*pj%f2(i,j)
          vr(i,j) = real(v(i,j)*pj%f1(i,j) + u(i,j)*pj%f2(i,j),rkx)
          ur(i,j) = real(tmp,rkx)
        end do
      end do
    end if
  end subroutine rotate2_lc

  pure subroutine backrotate2_lc(pj,u,v,ur,vr)
    implicit none
    class(regcm_projection) , intent(in) :: pj
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: u , v
    real(rkx) , pointer , dimension(:,:) , intent(out) :: ur , vr
    integer(ik4) :: i1 , i2 , j1 , j2 , i , j
    real(rk8) :: tmp
    i1 = lbound(u,1)
    i2 = ubound(u,1)
    j1 = lbound(u,2)
    j2 = ubound(u,2)
    if ( pj%stdlat >= d_zero ) then
      do j = j1 , j2
        do i = i1 , i2
          tmp = u(i,j)*pj%f1(i,j) - v(i,j)*pj%f2(i,j)
          vr(i,j) = real(u(i,j)*pj%f2(i,j) + v(i,j)*pj%f1(i,j),rkx)
          ur(i,j) = real(tmp,rkx)
        end do
      end do
    else
      do j = j1 , j2
        do i = i1 , i2
          tmp = u(i,j)*pj%f1(i,j) + v(i,j)*pj%f2(i,j)
          vr(i,j) = real(v(i,j)*pj%f1(i,j) - u(i,j)*pj%f2(i,j),rkx)
          ur(i,j) = real(tmp,rkx)
        end do
      end do
    end if
  end subroutine backrotate2_lc

  pure subroutine rotate3_lc(pj,u,v,ur,vr)
    implicit none
    class(regcm_projection) , intent(in) :: pj
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: u , v
    real(rkx) , pointer , dimension(:,:,:) , intent(out) :: ur , vr
    integer(ik4) :: i1 , i2 , j1 , j2 , k1 , k2 , i , j , k
    real(rk8) :: tmp
    i1 = lbound(u,1)
    i2 = ubound(u,1)
    j1 = lbound(u,2)
    j2 = ubound(u,2)
    k1 = lbound(u,3)
    k2 = ubound(u,3)
    if ( pj%stdlat >= d_zero ) then
      do k = k1 , k2
        do j = j1 , j2
          do i = i1 , i2
            tmp = u(i,j,k)*pj%f1(i,j) + v(i,j,k)*pj%f2(i,j)
            vr(i,j,k) = real(-u(i,j,k)*pj%f2(i,j) + v(i,j,k)*pj%f1(i,j),rkx)
            ur(i,j,k) = real(tmp,rkx)
          end do
        end do
      end do
    else
      do k = k1 , k2
        do j = j1 , j2
          do i = i1 , i2
            tmp = u(i,j,k)*pj%f1(i,j) - v(i,j,k)*pj%f2(i,j)
            vr(i,j,k) = real(v(i,j,k)*pj%f1(i,j) + u(i,j,k)*pj%f2(i,j),rkx)
            ur(i,j,k) = real(tmp,rkx)
          end do
        end do
      end do
    end if
  end subroutine rotate3_lc

  pure subroutine backrotate3_lc(pj,u,v,ur,vr)
    implicit none
    class(regcm_projection) , intent(in) :: pj
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: u , v
    real(rkx) , pointer , dimension(:,:,:) , intent(out) :: ur , vr
    integer(ik4) :: i1 , i2 , j1 , j2 , k1 , k2 , i , j , k
    real(rk8) :: tmp
    i1 = lbound(u,1)
    i2 = ubound(u,1)
    j1 = lbound(u,2)
    j2 = ubound(u,2)
    k1 = lbound(u,3)
    k2 = ubound(u,3)
    if ( pj%stdlat >= d_zero ) then
      do k = k1 , k2
        do j = j1 , j2
          do i = i1 , i2
            tmp = u(i,j,k)*pj%f1(i,j) - v(i,j,k)*pj%f2(i,j)
            vr(i,j,k) = real(u(i,j,k)*pj%f2(i,j) + v(i,j,k)*pj%f1(i,j),rkx)
            ur(i,j,k) = real(tmp,rkx)
          end do
        end do
      end do
    else
      do k = k1 , k2
        do j = j1 , j2
          do i = i1 , i2
            tmp = u(i,j,k)*pj%f1(i,j) + v(i,j,k)*pj%f2(i,j)
            vr(i,j,k) = real(v(i,j,k)*pj%f1(i,j) - u(i,j,k)*pj%f2(i,j),rkx)
            ur(i,j,k) = real(tmp,rkx)
          end do
        end do
      end do
    end if
  end subroutine backrotate3_lc

  pure subroutine rotate2_ps(pj,u,v,ur,vr)
    implicit none
    class(regcm_projection) , intent(in) :: pj
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: u , v
    real(rkx) , pointer , dimension(:,:) , intent(out) :: ur , vr
    integer(ik4) :: i1 , i2 , j1 , j2 , i , j
    real(rk8) :: tmp
    i1 = lbound(u,1)
    i2 = ubound(u,1)
    j1 = lbound(u,2)
    j2 = ubound(u,2)
    if ( pj%stdlat >= d_zero ) then
      do j = j1 , j2
        do i = i1 , i2
          tmp = u(i,j)*pj%f1(i,j) + v(i,j)*pj%f2(i,j)
          vr(i,j) = real(-u(i,j)*pj%f2(i,j) + v(i,j)*pj%f1(i,j),rkx)
          ur(i,j) = real(tmp,rkx)
        end do
      end do
    else
      do j = j1 , j2
        do i = i1 , i2
          tmp = u(i,j)*pj%f1(i,j) - v(i,j)*pj%f2(i,j)
          vr(i,j) = real(v(i,j)*pj%f1(i,j) + u(i,j)*pj%f2(i,j),rkx)
          ur(i,j) = real(tmp,rkx)
        end do
      end do
    end if
  end subroutine rotate2_ps

  pure subroutine backrotate2_ps(pj,u,v,ur,vr)
    implicit none
    class(regcm_projection) , intent(in) :: pj
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: u , v
    real(rkx) , pointer , dimension(:,:) , intent(out) :: ur , vr
    integer(ik4) :: i1 , i2 , j1 , j2 , i , j
    real(rk8) :: tmp
    i1 = lbound(u,1)
    i2 = ubound(u,1)
    j1 = lbound(u,2)
    j2 = ubound(u,2)
    if ( pj%stdlat >= d_zero ) then
      do j = j1 , j2
        do i = i1 , i2
          tmp = u(i,j)*pj%f1(i,j) - v(i,j)*pj%f2(i,j)
          vr(i,j) = real(u(i,j)*pj%f2(i,j) + v(i,j)*pj%f1(i,j),rkx)
          ur(i,j) = real(tmp,rkx)
        end do
      end do
    else
      do j = j1 , j2
        do i = i1 , i2
          tmp = u(i,j)*pj%f1(i,j) + v(i,j)*pj%f2(i,j)
          vr(i,j) = real(v(i,j)*pj%f1(i,j) - u(i,j)*pj%f2(i,j),rkx)
          ur(i,j) = real(tmp,rkx)
        end do
      end do
    end if
  end subroutine backrotate2_ps

  pure subroutine rotate3_ps(pj,u,v,ur,vr)
    implicit none
    class(regcm_projection) , intent(in) :: pj
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: u , v
    real(rkx) , pointer , dimension(:,:,:) , intent(out) :: ur , vr
    integer(ik4) :: i1 , i2 , j1 , j2 , k1 , k2 , i , j , k
    real(rk8) :: tmp
    i1 = lbound(u,1)
    i2 = ubound(u,1)
    j1 = lbound(u,2)
    j2 = ubound(u,2)
    k1 = lbound(u,3)
    k2 = ubound(u,3)
    if ( pj%stdlat >= d_zero ) then
      do k = k1 , k2
        do j = j1 , j2
          do i = i1 , i2
            tmp = u(i,j,k)*pj%f1(i,j) + v(i,j,k)*pj%f2(i,j)
            vr(i,j,k) = real(-u(i,j,k)*pj%f2(i,j) + v(i,j,k)*pj%f1(i,j),rkx)
            ur(i,j,k) = real(tmp,rkx)
          end do
        end do
      end do
    else
      do k = k1 , k2
        do j = j1 , j2
          do i = i1 , i2
            tmp = u(i,j,k)*pj%f1(i,j) - v(i,j,k)*pj%f2(i,j)
            vr(i,j,k) = real(v(i,j,k)*pj%f1(i,j) + u(i,j,k)*pj%f2(i,j),rkx)
            ur(i,j,k) = real(tmp,rkx)
          end do
        end do
      end do
    end if
  end subroutine rotate3_ps

  pure subroutine backrotate3_ps(pj,u,v,ur,vr)
    implicit none
    class(regcm_projection) , intent(in) :: pj
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: u , v
    real(rkx) , pointer , dimension(:,:,:) , intent(out) :: ur , vr
    integer(ik4) :: i1 , i2 , j1 , j2 , k1 , k2 , i , j , k
    real(rk8) :: tmp
    i1 = lbound(u,1)
    i2 = ubound(u,1)
    j1 = lbound(u,2)
    j2 = ubound(u,2)
    k1 = lbound(u,3)
    k2 = ubound(u,3)
    if ( pj%stdlat >= d_zero ) then
      do k = k1 , k2
        do j = j1 , j2
          do i = i1 , i2
            tmp = u(i,j,k)*pj%f1(i,j) - v(i,j,k)*pj%f2(i,j)
            vr(i,j,k) = real(u(i,j,k)*pj%f2(i,j) + v(i,j,k)*pj%f1(i,j),rkx)
            ur(i,j,k) = real(tmp,rkx)
          end do
        end do
      end do
    else
      do k = k1 , k2
        do j = j1 , j2
          do i = i1 , i2
            tmp = u(i,j,k)*pj%f1(i,j) + v(i,j,k)*pj%f2(i,j)
            vr(i,j,k) = real(v(i,j,k)*pj%f1(i,j) - u(i,j,k)*pj%f2(i,j),rkx)
            ur(i,j,k) = real(tmp,rkx)
          end do
        end do
      end do
    end if
  end subroutine backrotate3_ps

  pure subroutine rotate2_mc(pj,u,v,ur,vr)
    implicit none
    class(regcm_projection) , intent(in) :: pj
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: u , v
    real(rkx) , pointer , dimension(:,:) , intent(out) :: ur , vr
    ur(:,:) = u(:,:)
    vr(:,:) = v(:,:)
  end subroutine rotate2_mc

  pure subroutine rotate3_mc(pj,u,v,ur,vr)
    implicit none
    class(regcm_projection) , intent(in) :: pj
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: u , v
    real(rkx) , pointer , dimension(:,:,:) , intent(out) :: ur , vr
    ur(:,:,:) = u(:,:,:)
    vr(:,:,:) = v(:,:,:)
  end subroutine rotate3_mc

  pure subroutine rotate2_rc(pj,u,v,ur,vr)
    implicit none
    class(regcm_projection) , intent(in) :: pj
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: u , v
    real(rkx) , pointer , dimension(:,:) , intent(out) :: ur , vr
    integer(ik4) :: i1 , i2 , j1 , j2 , i , j
    real(rk8) :: tmp
    i1 = lbound(u,1)
    i2 = ubound(u,1)
    j1 = lbound(u,2)
    j2 = ubound(u,2)
    do j = j1 , j2
      do i = i1 , i2
        tmp = u(i,j)*pj%f1(i,j) - v(i,j)*pj%f2(i,j)
        vr(i,j) = real(v(i,j)*pj%f1(i,j) + u(i,j)*pj%f2(i,j),rkx)
        ur(i,j) = real(tmp,rkx)
      end do
    end do
  end subroutine rotate2_rc

  pure subroutine backrotate2_rc(pj,u,v,ur,vr)
    implicit none
    class(regcm_projection) , intent(in) :: pj
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: u , v
    real(rkx) , pointer , dimension(:,:) , intent(out) :: ur , vr
    integer(ik4) :: i1 , i2 , j1 , j2 , i , j
    real(rk8) :: tmp
    i1 = lbound(u,1)
    i2 = ubound(u,1)
    j1 = lbound(u,2)
    j2 = ubound(u,2)
    do j = j1 , j2
      do i = i1 , i2
        tmp = u(i,j)*pj%f1(i,j) + v(i,j)*pj%f2(i,j)
        vr(i,j) = real(v(i,j)*pj%f1(i,j) - u(i,j)*pj%f2(i,j),rkx)
        ur(i,j) = real(tmp,rkx)
      end do
    end do
  end subroutine backrotate2_rc

  pure subroutine rotate3_rc(pj,u,v,ur,vr)
    implicit none
    class(regcm_projection) , intent(in) :: pj
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: u , v
    real(rkx) , pointer , dimension(:,:,:) , intent(out) :: ur , vr
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
          vr(i,j,k) = real(v(i,j,k)*pj%f1(i,j) + u(i,j,k)*pj%f2(i,j),rkx)
          ur(i,j,k) = real(tmp,rkx)
        end do
      end do
    end do
  end subroutine rotate3_rc

  pure subroutine backrotate3_rc(pj,u,v,ur,vr)
    implicit none
    class(regcm_projection) , intent(in) :: pj
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: u , v
    real(rkx) , pointer , dimension(:,:,:) , intent(out) :: ur , vr
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
          vr(i,j,k) = real(v(i,j,k)*pj%f1(i,j) - u(i,j,k)*pj%f2(i,j),rkx)
          ur(i,j,k) = real(tmp,rkx)
        end do
      end do
    end do
  end subroutine backrotate3_rc

  pure subroutine rotate2_rl(pj,u,v,ur,vr)
    implicit none
    class(regcm_projection) , intent(in) :: pj
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: u , v
    real(rkx) , pointer , dimension(:,:) , intent(out) :: ur , vr
    integer(ik4) :: i1 , i2 , j1 , j2 , i , j
    real(rk8) :: tmp1 , tmp2
    i1 = lbound(u,1)
    i2 = ubound(u,1)
    j1 = lbound(u,2)
    j2 = ubound(u,2)
    do j = j1 , j2
      do i = i1 , i2
        if ( pj%f2(i,j) > 0.0_rk8 .and. pj%f3(i,j) > 0.0_rk8 ) then
          tmp1 = pj%f1(i,j) * u(i,j) + pj%f2(i,j) * v(i,j)
          tmp2 = pj%f3(i,j) * tmp1 + pj%f4(i,j) * v(i,j)
          vr(i,j) = real(tmp1,rkx)
          ur(i,j) = real(tmp2,rkx)
        else
          vr(i,j) = pj%f1(i,j) * v(i,j)
          ur(i,j) = pj%f4(i,j) * u(i,j)
        end if
      end do
    end do
  end subroutine rotate2_rl

  pure subroutine backrotate2_rl(pj,u,v,ur,vr)
    implicit none
    class(regcm_projection) , intent(in) :: pj
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: u , v
    real(rkx) , pointer , dimension(:,:) , intent(out) :: ur , vr
    integer(ik4) :: i1 , i2 , j1 , j2 , i , j
    real(rk8) :: tmp1 , tmp2
    i1 = lbound(u,1)
    i2 = ubound(u,1)
    j1 = lbound(u,2)
    j2 = ubound(u,2)
    do j = j1 , j2
      do i = i1 , i2
        if ( pj%f8(i,j) > 0.0_rk8 ) then
          tmp1 = pj%f5(i,j) * u(i,j) + pj%f6(i,j) * v(i,j)
          tmp2 = pj%f7(i,j) * (v(i,j)/tmp1*pj%f8(j,i))
          ur(i,j) = real(tmp2,rkx)
          vr(i,j) = real(tmp1,rkx)
        else
          if ( pj%f7(i,j) > 0.0_rk8 ) then
            ur(i,j) = pj%f5(i,j) * u(i,j)
            vr(i,j) = pj%f7(i,j) * v(i,j)
          else
            tmp1 = sqrt(u(i,j)*u(i,j)+v(i,j)*v(i,j))
            tmp2 = halfpi + atan2(v(i,j),u(i,j))
            if ( tmp2 > mathpi ) tmp2 = tmp2 - twopi
            if ( tmp2 < mathpi ) tmp2 = tmp2 + twopi
            ur(i,j) = -tmp1 * sin(pj%f5(i,j)-tmp2)
            vr(i,j) = -tmp1 * cos(pj%f5(i,j)-tmp2)
          end if
        end if
      end do
    end do
  end subroutine backrotate2_rl

  pure subroutine rotate3_rl(pj,u,v,ur,vr)
    implicit none
    class(regcm_projection) , intent(in) :: pj
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: u , v
    real(rkx) , pointer , dimension(:,:,:) , intent(out) :: ur , vr
    integer(ik4) :: i1 , i2 , j1 , j2 , k1 , k2 , i , j , k
    real(rk8) :: tmp1 , tmp2
    i1 = lbound(u,1)
    i2 = ubound(u,1)
    j1 = lbound(u,2)
    j2 = ubound(u,2)
    k1 = lbound(u,3)
    k2 = ubound(u,3)
    do k = k1 , k2
      do j = j1 , j2
        do i = i1 , i2
          if ( pj%f2(i,j) > 0.0_rk8 .and. pj%f3(i,j) > 0.0_rk8 ) then
            tmp1 = pj%f1(i,j) * u(i,j,k) + pj%f2(i,j) * v(i,j,k)
            tmp2 = pj%f3(i,j) * tmp1 + pj%f4(i,j) * v(i,j,k)
            vr(i,j,k) = real(tmp1,rkx)
            ur(i,j,k) = real(tmp2,rkx)
          else
            vr(i,j,k) = pj%f1(i,j) * v(i,j,k)
            ur(i,j,k) = pj%f4(i,j) * u(i,j,k)
          end if
        end do
      end do
    end do
  end subroutine rotate3_rl

  pure subroutine backrotate3_rl(pj,u,v,ur,vr)
    implicit none
    class(regcm_projection) , intent(in) :: pj
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: u , v
    real(rkx) , pointer , dimension(:,:,:) , intent(out) :: ur , vr
    integer(ik4) :: i1 , i2 , j1 , j2 , k1 , k2 , i , j , k
    real(rk8) :: tmp1 , tmp2
    i1 = lbound(u,1)
    i2 = ubound(u,1)
    j1 = lbound(u,2)
    j2 = ubound(u,2)
    k1 = lbound(u,3)
    k2 = ubound(u,3)
    do k = k1 , k2
      do j = j1 , j2
        do i = i1 , i2
          if ( pj%f8(i,j) > 0.0_rk8 ) then
            tmp1 = pj%f5(i,j) * u(i,j,k) + pj%f6(i,j) * v(i,j,k)
            tmp2 = pj%f7(i,j) * (v(i,j,k)/tmp1*pj%f8(j,i))
            ur(i,j,k) = real(tmp2,rkx)
            vr(i,j,k) = real(tmp1,rkx)
          else
            if ( pj%f7(i,j) > 0.0_rk8 ) then
              ur(i,j,k) = pj%f5(i,j) * u(i,j,k)
              vr(i,j,k) = pj%f7(i,j) * v(i,j,k)
            else
              tmp1 = sqrt(u(i,j,k)*u(i,j,k)+v(i,j,k)*v(i,j,k))
              tmp2 = halfpi + atan2(v(i,j,k),u(i,j,k))
              if ( tmp2 > mathpi ) tmp2 = tmp2 - twopi
              if ( tmp2 < mathpi ) tmp2 = tmp2 + twopi
              ur(i,j,k) = -tmp1 * sin(pj%f5(i,j)-tmp2)
              vr(i,j,k) = -tmp1 * cos(pj%f5(i,j)-tmp2)
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
    real(rkx) , pointer , dimension(:,:) , intent(out) :: xmap
    xmap = fac_rl(pj,xlat,xlon)
  end subroutine mapfac_rl

  pure subroutine mapfac_lc(pj,xlat,xlon,xmap)
    implicit none
    class(regcm_projection) , intent(in) :: pj
    real(rkx) , pointer , dimension(:,:) , intent(in) :: xlat , xlon
    real(rkx) , pointer , dimension(:,:) , intent(out) :: xmap
    xmap = fac_lc(pj,xlat)
  end subroutine mapfac_lc

  pure subroutine mapfac_ps(pj,xlat,xlon,xmap)
    implicit none
    class(regcm_projection) , intent(in) :: pj
    real(rkx) , pointer , dimension(:,:) , intent(in) :: xlat , xlon
    real(rkx) , pointer , dimension(:,:) , intent(out) :: xmap
    xmap = fac_ps(pj,xlat)
  end subroutine mapfac_ps

  pure subroutine mapfac_mc(pj,xlat,xlon,xmap)
    implicit none
    class(regcm_projection) , intent(in) :: pj
    real(rkx) , pointer , dimension(:,:) , intent(in) :: xlat , xlon
    real(rkx) , pointer , dimension(:,:) , intent(out) :: xmap
    xmap = fac_mc(pj,xlat)
  end subroutine mapfac_mc

  pure subroutine mapfac_ll(pj,xlat,xlon,xmap)
    implicit none
    class(regcm_projection) , intent(in) :: pj
    real(rkx) , pointer , dimension(:,:) , intent(in) :: xlat , xlon
    real(rkx) , pointer , dimension(:,:) , intent(out) :: xmap
    xmap = d_one
  end subroutine mapfac_ll

  pure subroutine mapfac_rc(pj,xlat,xlon,xmap)
    implicit none
    class(regcm_projection) , intent(in) :: pj
    real(rkx) , pointer , dimension(:,:) , intent(in) :: xlat , xlon
    real(rkx) , pointer , dimension(:,:) , intent(out) :: xmap
    xmap = fac_rc(pj,xlat,xlon)
  end subroutine mapfac_rc

end module mod_projections
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
