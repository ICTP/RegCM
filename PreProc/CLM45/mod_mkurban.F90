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
module mod_mkurban
  use mod_realkinds
  use mod_intkinds
  use mod_dynparam
  use mod_constants
  use mod_grid
  use mod_rdldtr
  use mod_message

  implicit none

  private

  public :: mkurban_base , mkurban_param
  public :: ip2d , ip3d , ip4d

  character(len=16) , parameter :: densdim = 'density_class'
  character(len=16) , parameter :: regiondim = 'region'
  character(len=16) , parameter :: varname = 'PCT_URBAN'
  character(len=16) , parameter :: regionname = 'REGION_ID'

  integer(ik4) , public , parameter :: npu2d = 14
  integer(ik4) , public , parameter :: npu3d = 6
  integer(ik4) , public , parameter :: npu4d = 4
  integer(ik4) , parameter :: nparam = npu2d + npu3d + npu4d

  character(len=16) , dimension(npu2d) , public , parameter :: parm2d = &
    ['CANYON_HWR      ', 'EM_IMPROAD      ', 'EM_PERROAD      ', &
      'EM_ROOF         ', 'EM_WALL         ', 'HT_ROOF         ', &
      'NLEV_IMPROAD    ', 'THICK_ROOF      ', 'THICK_WALL      ', &
      'T_BUILDING_MAX  ', 'T_BUILDING_MIN  ', 'WIND_HGT_CANYON ', &
      'WTLUNIT_ROOF    ', 'WTROAD_PERV     ']
  character(len=36) , dimension(npu2d) , public , parameter :: lngn2d = &
    ['canyon height to width ratio        ', &
      'emissivity of impervious road       ', &
      'emissivity of pervious road         ', &
      'emissivity of roof                  ', &
      'emissivity of wall                  ', &
      'height of roof                      ', &
      'number of impervious road layers    ', &
      'thickness of roof                   ', &
      'thickness of wall                   ', &
      'maximum intern building temperature ', &
      'minimum intern building temperature ', &
      'height of wind in canyon            ', &
      'fraction of roof                    ', &
      'fraction of pervious road           ']
  character(len=4) , dimension(npu2d) , public , parameter :: unit2d = &
    ['1   ', '1   ', '1   ', '1   ', '1   ', 'm   ', &
      '1   ', 'm   ', 'm   ', 'K   ', 'K   ', 'm   ', &
      '1   ', '1   ']
  character(len=16) , dimension(npu3d) , public , parameter :: parm3d = &
    ['CV_IMPROAD      ', 'CV_ROOF         ', 'CV_WALL         ', &
      'TK_IMPROAD      ', 'TK_ROOF         ', 'TK_WALL         ']
  character(len=36) , dimension(npu3d) , public , parameter :: lngn3d = &
    ['vol heat capacity of impervious road', &
      'vol heat capacity of roof           ', &
      'vol heat capacity of wall           ', &
      'thermal conductivity of imperv road ', &
      'thermal conductivity of roof        ', &
      'thermal conductivity of wall        ']
  character(len=8) , dimension(npu3d) , public , parameter :: unit3d = &
    ['J/m^3*K ' , 'J/m^3*K ', 'J/m^3*K ', &
      'W/m*K   ' , 'W/m*K   ', 'W/m*K   ']
  character(len=16) , dimension(npu4d) , public , parameter :: parm4d = &
    ['ALB_IMPROAD     ', 'ALB_PERROAD     ', 'ALB_ROOF        ', &
      'ALB_WALL        ']
  character(len=36) , dimension(npu4d) , public , parameter :: lngn4d = &
    ['albedo of impervious road           ', &
      'albedo of pervious road             ', &
      'albedo of roof                      ', &
      'albedo of wall                      ']
  character(len=2) , dimension(npu4d) , public , parameter :: unit4d = &
    ['1 ', '1 ', '1 ', '1 ']
  character(len=16) , dimension(nparam) :: parmname

  integer(ik4) , dimension(nparam) , parameter :: parmdim = &
    [3,3,3,3,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,5,5,5,5]

  real(rkx) , parameter :: vcutoff =  1.0_rkx

  contains

  subroutine mkurban_base(urbanfile,mask,urban)
    implicit none
    character(len=*) , intent(in) :: urbanfile
    real(rkx) , dimension(:,:) , intent(in) :: mask
    real(rkx) , dimension(:,:,:) , intent(out) :: urban
    integer(ik4) :: i , j , n , nurban , nmax(1)
    type(globalfile) :: gfile
    character(len=256) :: inpfile
    real(rkx) , dimension(:,:) , allocatable :: usum

    nurban = size(urban,3)
    inpfile = trim(inpglob)//pthsep//'CLM45'// &
                             pthsep//'surface'//pthsep//urbanfile
    call gfopen(gfile,inpfile,xlat,xlon,ds*nsg,roidem,i_band)
    call gfread(gfile,varname,urban,h_missing_value)
    call gfclose(gfile)

    do n  = 1 , nurban
      do i = 1 , iysg
        do j = 1 , jxsg
          if ( mask(j,i) < 0.5_rkx ) then
            urban(j,i,n) = h_missing_value
          else
            if ( urban(j,i,n) < d_zero ) then
              urban(j,i,n) = d_zero
            else
              urban(j,i,n) = min(max(0,nint(urban(j,i,n))),100)
            end if
          end if
        end do
      end do
    end do

    allocate(usum(jxsg,iysg))
    usum = sum(urban,3)
    do n  = 1 , nurban
      do i = 1 , iysg
        do j = 1 , jxsg
          if ( mask(j,i) > 0.5_rkx ) then
            if ( usum(j,i) < vcutoff ) then
              urban(j,i,n) = d_zero
            else if ( usum(j,i) > 100.0_rkx ) then
              nmax = maxloc(urban(j,i,:))
              urban(j,i,nmax(1)) = urban(j,i,nmax(1))-(usum(j,i)-100.0_rkx)
            end if
          end if
        end do
      end do
    end do
    deallocate(usum)
  end subroutine mkurban_base

  subroutine mkurban_param(urbanfile,mask,urban3d,urban4d,urban5d)
    implicit none
    character(len=*) , intent(in) :: urbanfile
    real(rkx) , dimension(:,:) , intent(in) :: mask
    real(rkx) , dimension(:,:,:,:) , intent(out) :: urban3d
    real(rkx) , dimension(:,:,:,:,:) , intent(out) :: urban4d
    real(rkx) , dimension(:,:,:,:,:,:) , intent(out) :: urban5d
    integer(ik4) :: ipt , ip , n , i , j , n1 , n2 , n3
    integer(ik4) :: i4 , i5 , i6
    type(globalfile) :: gfile

    character(len=256) :: inpfile

    ipt = 1
    do ip = 1 , npu2d
      parmname(ipt) = parm2d(ip)
      ipt = ipt + 1
    end do
    do ip = 1 , npu3d
      parmname(ipt) = parm3d(ip)
      ipt = ipt + 1
    end do
    do ip = 1 , npu4d
      parmname(ipt) = parm4d(ip)
      ipt = ipt + 1
    end do

    inpfile = trim(inpglob)//pthsep//'CLM45'// &
                             pthsep//'surface'//pthsep//urbanfile
    call gfopen(gfile,inpfile,xlat,xlon,ds*nsg,roidem,i_band)

    do n = 1 , nparam
      select case (parmdim(n))
        case (3)
          i4 = ip2d(parmname(n))
          call gfread(gfile,parmname(n),regiondim,regionname, &
                        urban3d(:,:,:,i4),.true.,h_missing_value)
          where ( urban3d(:,:,:,i4) < d_zero )
            urban3d(:,:,:,i4) = h_missing_value
          end where
          do n1 = 1 , size(urban3d,3)
            call bestaround(urban3d(:,:,n1,i4),h_missing_value)
            do i = 1 , iysg
              do j = 1 , jxsg
                if ( mask(j,i) < 0.5_rkx ) then
                  urban3d(j,i,n1,i4) = h_missing_value
                end if
              end do
            end do
          end do
        case (4)
          i5 = ip3d(parmname(n))
          call gfread(gfile,parmname(n),regiondim,regionname, &
                        urban4d(:,:,:,:,i5),.true.,h_missing_value)
          where ( urban4d(:,:,:,:,i5) < d_zero )
            urban4d(:,:,:,:,i5) = h_missing_value
          end where
          do n2 = 1 , size(urban4d,4)
            do n1 = 1 , size(urban4d,3)
              call bestaround(urban4d(:,:,n1,n2,i5),h_missing_value)
              do i = 1 , iysg
                do j = 1 , jxsg
                  if ( mask(j,i) < 0.5_rkx ) then
                    urban4d(j,i,n1,n2,i5) = h_missing_value
                  end if
                end do
              end do
            end do
          end do
        case (5)
          i6 = ip4d(parmname(n))
          call gfread(gfile,parmname(n),regiondim,regionname, &
                        urban5d(:,:,:,:,:,i6),.true.,h_missing_value)
          where ( urban5d(:,:,:,:,:,i6) < d_zero )
            urban5d(:,:,:,:,:,i6) = h_missing_value
          end where
          do n3 = 1 , size(urban5d,5)
            do n2 = 1 , size(urban5d,4)
              do n1 = 1 , size(urban5d,3)
                call bestaround(urban5d(:,:,n1,n2,n3,i6),h_missing_value)
                do i = 1 , iysg
                  do j = 1 , jxsg
                    if ( mask(j,i) < 0.5_rkx ) then
                      urban5d(j,i,n1,n2,n3,i6) = h_missing_value
                    end if
                  end do
                end do
              end do
            end do
          end do
        case default
          call die(__FILE__,'Variable dimension not implemented',__LINE__)
      end select
    end do
    call gfclose(gfile)
  end subroutine mkurban_param

  integer(ik4) function ip2d(pname) result(ip)
    implicit none
    character(len=*) :: pname
    do ip = 1 , npu2d
      if ( pname == parm2d(ip) ) then
        return
      end if
    end do
    call die(__FILE__,'Variable '//pname//' NOT FOUND',__LINE__)
  end function ip2d

  integer(ik4) function ip3d(pname) result(ip)
    implicit none
    character(len=*) :: pname
    do ip = 1 , npu3d
      if ( pname == parm3d(ip) ) then
        return
      end if
    end do
    call die(__FILE__,'Variable '//pname//' NOT FOUND',__LINE__)
  end function ip3d

  integer(ik4) function ip4d(pname) result(ip)
    implicit none
    character(len=*) :: pname
    do ip = 1 , npu4d
      if ( pname == parm4d(ip) ) then
        return
      end if
    end do
    call die(__FILE__,'Variable '//pname//' NOT FOUND',__LINE__)
  end function ip4d

end module mod_mkurban
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
