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

module mod_cmip7

  use mod_intkinds
  use mod_realkinds
  use mod_cmip7_helper
  use mod_message
  use mod_date
  use mod_stdio
  use mod_dynparam
  use mod_memutil
  use mod_grid
  use mod_kdinterp
  use mod_write
  use mod_vectutil
  use mod_mksst
  use mod_humid
  use mod_hgt
  use mod_vertint
  use netcdf
  use mod_cmip7_ecea

  implicit none

  private

  public :: init_cmip7, get_cmip7, conclude_cmip7

  ! Pressure levels to interpolate to if dataset is on model sigma levels.
  integer(ik4), parameter :: nipl = 41
  real(rkx), target, dimension(nipl) :: fplev = &
     [ 1030.0_rkx,1020.0_rkx,1010.0_rkx,1000.0_rkx, 975.0_rkx, 950.0_rkx, &
        925.0_rkx, 900.0_rkx, 875.0_rkx, 850.0_rkx, 825.0_rkx, 800.0_rkx, &
        775.0_rkx, 750.0_rkx, 700.0_rkx, 650.0_rkx, 600.0_rkx, 550.0_rkx, &
        500.0_rkx, 450.0_rkx, 425.0_rkx, 400.0_rkx, 350.0_rkx, 300.0_rkx, &
        250.0_rkx, 225.0_rkx, 200.0_rkx, 175.0_rkx, 150.0_rkx, 125.0_rkx, &
        100.0_rkx,  70.0_rkx,  50.0_rkx,  30.0_rkx,  20.0_rkx,  10.0_rkx, &
          7.0_rkx,   5.0_rkx,   3.0_rkx,   2.0_rkx,   1.0_rkx ]

  type(cmip7_2d_var), pointer :: orog => null( )
  type(cmip7_2d_var), pointer :: ps => null( )
  type(cmip7_3d_var), pointer :: ta => null( )
  type(cmip7_3d_var), pointer :: qa => null( )
  type(cmip7_3d_var), pointer :: ua => null( )
  type(cmip7_3d_var), pointer :: va => null( )
  type(cmip7_3d_var), pointer :: zg => null( )

  real(rkx), dimension(:), pointer, contiguous :: sigmar
  real(rkx) :: pss, pst

  real(rkx), dimension(:,:,:), pointer, contiguous :: pa_in, zp_in
  real(rkx), dimension(:,:,:), pointer, contiguous :: tvar
  real(rkx), dimension(:,:,:), pointer, contiguous :: uvar
  real(rkx), dimension(:,:,:), pointer, contiguous :: vvar
  real(rkx), dimension(:,:,:), pointer, contiguous :: qvar
  real(rkx), dimension(:,:,:), pointer, contiguous :: zvar

  real(rkx), dimension(:,:,:), pointer, contiguous :: tah
  real(rkx), dimension(:,:,:), pointer, contiguous :: uah
  real(rkx), dimension(:,:,:), pointer, contiguous :: vah
  real(rkx), dimension(:,:,:), pointer, contiguous :: qah
  real(rkx), dimension(:,:,:), pointer, contiguous :: zgh

  real(rkx), pointer, contiguous, dimension(:,:,:) :: dv, du
  real(rkx), pointer, contiguous, dimension(:,:,:) :: hv, hu

  integer(ik4) :: nkin

  logical, parameter :: only_coord = .true.

  contains

    subroutine init_cmip7(idate)
      implicit none
      type(rcm_time_and_date), intent(in) :: idate
      integer(ik4) :: k
      select case (cmip7_model)
        case ( 'EC-Earth3-ESM-1-1' )
          allocate(ps,ua,va,ta,qa,zg,orog)
          ps%vname = 'ps'
          ua%vname = 'ua'
          va%vname = 'va'
          ta%vname = 'ta'
          qa%vname = 'hus'
          zg%vname = 'zg'
          orog%vname = 'orog'
          call read_fx_ecea(orog)
          if ( idynamic == 3 ) then
            allocate(orog%hint(3))
            call h_interpolator_create(orog%hint(1),orog%hcoord%lat1d, &
              orog%hcoord%lon1d, xlat, xlon)
            call h_interpolator_create(orog%hint(2),orog%hcoord%lat1d, &
              orog%hcoord%lon1d, ulat, ulon)
            call h_interpolator_create(orog%hint(3),orog%hcoord%lat1d, &
              orog%hcoord%lon1d, vlat, vlon)
          else
            allocate(orog%hint(2))
            call h_interpolator_create(orog%hint(1),orog%hcoord%lat1d, &
              orog%hcoord%lon1d, xlat, xlon)
            call h_interpolator_create(orog%hint(2),orog%hcoord%lat1d, &
              orog%hcoord%lon1d, dlat, dlon)
          end if
          ps%hint => orog%hint
          ta%hint => orog%hint
          ua%hint => orog%hint
          va%hint => orog%hint
          qa%hint => orog%hint
          zg%hint => orog%hint
          ps%hcoord => orog%hcoord
          ta%hcoord => orog%hcoord
          ua%hcoord => orog%hcoord
          va%hcoord => orog%hcoord
          qa%hcoord => orog%hcoord
          zg%hcoord => orog%hcoord
          call read_2d_ecea(idate,ps,only_coord)
          call read_3d_ecea(idate,ta,only_coord)
          ua%vcoord => ta%vcoord
          va%vcoord => ta%vcoord
          qa%vcoord => qa%vcoord
          call read_3d_ecea(idate,ua,only_coord)
          call read_3d_ecea(idate,va,only_coord)
          call read_3d_ecea(idate,qa,only_coord)
          nkin = nipl
          call getmem(sigmar,1,nkin,'cmip7:ecea:sigmar')
          do k = 1, nkin
            sigmar(k) = (fplev(k)-fplev(nkin))/(fplev(1)-fplev(nkin))
          end do
          pss = (fplev(1)-fplev(nkin))/10.0_rkx
          pst = fplev(nkin)/10.0_rkx
          call getmem(pa_in,1,ta%ni,1,ta%nj,1,ta%nk,'cmip7:ecea:pa_in')
          call getmem(zp_in,1,ta%ni,1,ta%nj,1,ta%nk,'cmip7:ecea:zp_in')
          call getmem(tvar,1,ta%ni,1,ta%nj,1,nkin,'cmip7:ecea:tvar')
          call getmem(uvar,1,ua%ni,1,ua%nj,1,nkin,'cmip7:ecea:uvar')
          call getmem(vvar,1,va%ni,1,va%nj,1,nkin,'cmip7:ecea:vvar')
          call getmem(qvar,1,qa%ni,1,qa%nj,1,nkin,'cmip7:ecea:qvar')
          call getmem(zvar,1,ta%ni,1,ta%nj,1,nkin,'cmip7:ecea:zvar')
          call getmem(tah,1,jx,1,iy,1,nkin,'cmip7:ecea:tah')
          call getmem(qah,1,jx,1,iy,1,nkin,'cmip7:ecea:qah')
          call getmem(uah,1,jx,1,iy,1,nkin,'cmip7:ecea:uah')
          call getmem(vah,1,jx,1,iy,1,nkin,'cmip7:ecea:vah')
          call getmem(zgh,1,jx,1,iy,1,nkin,'cmip7:ecea:zgh')
        case default
          call die(__FILE__,'Unsupported cmi76 model. Stop at line ',__LINE__)
      end select

      if ( idynamic == 3 ) then
        call getmem(du,1,jx,1,iy,1,nkin,'cmip7:du')
        call getmem(dv,1,jx,1,iy,1,nkin,'cmip7:dv')
        call getmem(hu,1,jx,1,iy,1,nkin,'cmip7:hu')
        call getmem(hv,1,jx,1,iy,1,nkin,'cmip7:hv')
      end if

      write (stdout,*) 'Read in Static fields OK'
    end subroutine init_cmip7

    subroutine get_cmip7(idate)
      implicit none
      type(rcm_time_and_date), intent(in) :: idate
      integer(ik4) :: i, j, k

      if ( dattyp == 'CMIP7' ) then
        select case (cmip7_model)
          case ( 'EC-Earth3-ESM-1-1' )
!$OMP SECTIONS
!$OMP SECTION
            call read_2d_ecea(idate,ps)
!$OMP SECTION
            call read_3d_ecea(idate,ua)
!$OMP SECTION
            call read_3d_ecea(idate,va)
!$OMP SECTION
            call read_3d_ecea(idate,ta)
!$OMP SECTION
            call read_3d_ecea(idate,qa)
!$OMP END SECTIONS
            call sph2mxr(qa%var,qa%ni,qa%nj,qa%nk)
            do k = 1, ta%nk
              do j = 1, ta%nj
                do i = 1, ta%ni
                  pa_in(i,j,k) = ta%vcoord%ak(k) + ta%vcoord%bk(k) * ps%var(i,j)
                end do
              end do
            end do
            pa_in = pa_in * 0.01_rkx
!$OMP SECTIONS
!$OMP SECTION
            call intlin(uvar,ua%var,pa_in,ua%ni,ua%nj,ua%nk,fplev,nkin)
!$OMP SECTION
            call intlin(vvar,va%var,pa_in,va%ni,va%nj,va%nk,fplev,nkin)
!$OMP SECTION
            call intlog(tvar,ta%var,pa_in,ta%ni,ta%nj,ta%nk,fplev,nkin)
!$OMP SECTION
            call intlin(qvar,qa%var,pa_in,qa%ni,qa%nj,qa%nk,fplev,nkin)
!$OMP END SECTIONS
            ps%var = ps%var * 0.01_rkx
            call htsig(zp_in,ta%var,pa_in,qa%var,ps%var,orog%var)
            call height(zvar,zp_in,ta%var,ps%var,pa_in,orog%var, &
                        ta%ni,ta%nj,ta%nk,fplev,nkin)
          case default
            call die(__FILE__,'Unsupported cmip7 model. Stop at line ',__LINE__)
        end select
      end if

      write (stdout,*) 'Read in fields at Date: ', tochar(idate)

!$OMP SECTIONS
!$OMP SECTION
      call h_interpolate_cont(ta%hint(1),tvar,tah)
!$OMP SECTION
      call h_interpolate_cont(qa%hint(1),qvar,qah)
!$OMP SECTION
      call h_interpolate_cont(ta%hint(1),zvar,zgh)
!$OMP END SECTIONS
      if ( idynamic == 3 ) then
!$OMP SECTIONS
!$OMP SECTION
        call h_interpolate_cont(ua%hint(2),uvar,uah)
!$OMP SECTION
        call h_interpolate_cont(ua%hint(2),vvar,dv)
!$OMP SECTION
        call h_interpolate_cont(ua%hint(3),uvar,du)
!$OMP SECTION
        call h_interpolate_cont(ua%hint(3),vvar,vah)
        call pju%wind_rotate(uah,dv)
        call pjv%wind_rotate(du,vah)
!$OMP END SECTIONS
      else
!$OMP SECTIONS
!$OMP SECTION
        call h_interpolate_cont(ua%hint(2),uvar,uah)
!$OMP SECTION
        call h_interpolate_cont(ua%hint(2),vvar,vah)
        call pjd%wind_rotate(uah,vah)
!$OMP END SECTIONS
      end if

      if ( idynamic == 3 ) then
        call ucrs2dot(hu,zgh,jx,iy,nkin,i_band)
        call vcrs2dot(hv,zgh,jx,iy,nkin,i_crm)
        call intzps(ps4,topogm,tah,zgh,pss,sigmar,pst, &
                    xlat,yeardayfrac(idate),dayspy,jx,iy,nkin)
        call intz3(ts4,tah,zgh,topogm,0.6_rkx,0.5_rkx,0.85_rkx,jx,iy,nkin)
      else
        call intgtb(pa,za,tlayer,topogm,tah,zgh,pss,sigmar,pst,jx,iy,nkin)
        call intpsn(ps4,topogm,pa,za,tlayer,ptop,jx,iy)
        call crs2dot(pd4,ps4,jx,iy,i_band,i_crm)
        call intv3(ts4,tah,ps4,pss,sigmar,ptop,pst,jx,iy,nkin)
      end if

      call readsst(ts4,idate)

      if ( idynamic == 3 ) then
!$OMP SECTIONS
!$OMP SECTION
        call intz1(u4,uah,zetau,hu,topou,jx,iy,kz,nkin,0.6_rkx,0.2_rkx,0.2_rkx)
!$OMP SECTION
        call intz1(v4,vah,zetav,hv,topov,jx,iy,kz,nkin,0.6_rkx,0.2_rkx,0.2_rkx)
!$OMP SECTION
        call intz1(t4,tah,z0,zgh,topogm,jx,iy,kz,nkin,0.6_rkx,0.5_rkx,0.85_rkx)
!$OMP SECTION
        call intz1(q4,qah,z0,zgh,topogm,jx,iy,kz,nkin,0.7_rkx,0.4_rkx,0.7_rkx)
!$OMP END SECTIONS
      else
!$OMP SECTIONS
!$OMP SECTION
        call intv1(u4,uah,pd4,sigmah,pss,sigmar,ptop,pst,jx,iy,kz,nkin,1)
!$OMP SECTION
        call intv1(v4,vah,pd4,sigmah,pss,sigmar,ptop,pst,jx,iy,kz,nkin,1)
!$OMP SECTION
        call intv2(t4,tah,ps4,sigmah,pss,sigmar,ptop,pst,jx,iy,kz,nkin)
!$OMP SECTION
        call intv1(q4,qah,ps4,sigmah,pss,sigmar,ptop,pst,jx,iy,kz,nkin,1)
!$OMP END SECTIONS
      end if
    end subroutine get_cmip7

    subroutine conclude_cmip7( )
      implicit none
      if ( idynamic == 3 ) then
        call h_interpolator_destroy(orog%hint(1))
        call h_interpolator_destroy(orog%hint(2))
        call h_interpolator_destroy(orog%hint(3))
      else
        call h_interpolator_destroy(orog%hint(1))
        call h_interpolator_destroy(orog%hint(2))
      end if
      deallocate(orog%hint)
    end subroutine conclude_cmip7

end module mod_cmip7

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
