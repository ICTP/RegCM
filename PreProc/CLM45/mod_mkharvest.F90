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
module mod_mkharvest
#ifdef CN
  use mod_realkinds
  use mod_intkinds
  use mod_dynparam
  use mod_grid
  use mod_rdldtr

  implicit none

  private

  public :: mkharvest

  integer , parameter :: nvarc = 6

  character(len=16) , parameter , dimension(nvarc):: varname = &
          (/ 'HARVEST_VH1' , 'HARVEST_VH2' , 'HARVEST_SH1' , &
             'HARVEST_SH2' , 'HARVEST_SH3' , 'GRAZING    '/)
  character(len=16) , parameter :: maskname = 'LANDMASK'

  real(rkx) :: vmisdat = -9999.0_rkx

  contains

  subroutine mkharvest(harvest,iyear)
    implicit none
    real(rkx) , dimension(:,:,:) , intent(out) :: harvest
    integer(ik4) , intent(in) :: iyear
    integer(ik4) :: i , j , n
    real(rkx) , pointer , dimension(:,:) :: mask
    character(len=32) :: p1 , p2
    character(len=4) :: cy
    type(globalfile) :: gfile

    character(len=256) :: inpfile

    p1 = 'dynamic'
    p2 = '.'
    write(cy,'(i0.4)') iyear

    if ( iyear > 2005 ) then
      select case (dattyp(4:5))
        case ('RF')
          continue
        case ('26')
          p2 = 'SCENARIO'//pthsep//'RCP2.6'
        case ('45')
          p2 = 'SCENARIO'//pthsep//'RCP4.5'
        case ('60')
          p2 = 'SCENARIO'//pthsep//'RCP6.0'
        case ('85', '15')
          p2 = 'SCENARIO'//pthsep//'RCP8.5'
        case default
          if ( dattyp /= "EIN15" .and. &
               dattyp(1:4) /= "NNRP" .and. &
               dattyp /= "JRA55" ) then
            call die(__FILE__, &
              'Dynamic landuse only supported for CMIP5',__LINE__)
          end if
      end select
    end if

    inpfile = trim(inpglob)//pthsep//'CLM45'//pthsep//'surface'// &
            pthsep//trim(p1)//pthsep//trim(p2)//pthsep//&
            'mksrf_landuse_'//cy//'.nc'
    allocate(mask(jxsg,iysg))

    call gfopen(gfile,inpfile,xlat,xlon,ds*nsg,i_band)
    call gfread(gfile,maskname,mask)
    do n = 1 , nvarc
      call gfread(gfile,varname(n),harvest(:,:,n))
    end do

    do i = 1 , iysg
      do j = 1 , jxsg
        if ( mask(j,i) < 1.0_rkx ) then
          harvest(j,i,:) = vmisdat
        end if
      end do
    end do
    deallocate(mask)
    call gfclose(gfile)
  end subroutine mkharvest
#endif

end module mod_mkharvest
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
