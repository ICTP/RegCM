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

module mod_rad_outrad

  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_dynparam
  use mod_mpmessage
  use mod_rad_common
  use mod_runparams
  use mod_outvars
  use mod_regcm_types

  implicit none

  private

  public :: allocate_mod_rad_outrad , radout

  integer(ik4) :: npr

  contains

  subroutine allocate_mod_rad_outrad
    implicit none
    npr = (jci2-jci1+1)*(ici2-ici1+1)
  end subroutine allocate_mod_rad_outrad

  subroutine radout(lout,solin,solout,frsa,clrst,clrss,qrs,lwout,        &
                    frla,clrlt,clrls,qrl,slwd,sols,soll,solsd,solld,     &
                    totcf,totwv,totcl,totci,cld,clwp,abv,sol,aeradfo,    &
                    aeradfos,aerlwfo,aerlwfos,tauxar3d,tauasc3d,gtota3d, &
                    deltaz,o3,outtaucl,outtauci,asaeradfo,asaeradfos,    &
                    asaerlwfo,asaerlwfos,r2m,m2r)
    implicit none
    !
    ! copy radiation output quantities to model buffer
    !
    ! change units of the radiative fluxes from cgs to mks
    !
    ! compute the total radiative heat flux at the surface for
    ! the surface temperature computation
    !
    !     input/output arguments
    !
    ! solin  - instantaneous incident solar
    ! solout - outgoing short wave flux
    ! frsa   - surface absorbed solar flux
    ! clrst  - clear sky total column abs solar flux
    ! clrss  - clear sky surface absorbed solar flux
    ! qrs    - solar heating rate
    ! lwout  - outgoing lw flux
    ! frla   - longwave cooling of surface (up-dwn flx)
    ! clrlt  - clr sky net up flx top of model (up-dwn f
    ! clrls  - clr sky lw cooling of srf (up-dwn flx)
    ! qrl    - longwave cooling rate
    ! slwd   - surface longwave down flux
    ! cld    - cloud fractional cover
    ! clwp   - cloud liquid water path
    ! soll   - Downward solar rad onto surface (lw direct)
    ! solld  - Downward solar rad onto surface (lw diffuse)
    ! sols   - Downward solar rad onto surface (sw direct)
    ! solsd  - Downward solar rad onto surface (sw diffuse)
    !
    logical , intent(in) :: lout ! Preapre data for outfile
    real(rkx) , pointer , dimension(:) , intent(inout) :: clrls , clrlt ,  &
           clrss , clrst , lwout , frla , frsa , solout , slwd , solin ,   &
           soll , solld , sols , solsd , totcf , totcl , totci , totwv ,   &
           abv , sol , aeradfo , aeradfos , aerlwfo , aerlwfos
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: cld , clwp ,   &
           qrl , qrs , deltaz , o3
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: outtaucl ,  &
           outtauci , tauxar3d , tauasc3d , gtota3d
    real(rkx) , pointer , dimension(:) , intent(inout) :: asaeradfo , &
           asaeradfos , asaerlwfo , asaerlwfos
    type(mod_2_rad) , intent(in) :: m2r
    type(rad_2_mod) , intent(inout) :: r2m

    integer(ik4) :: i , j , k , n
    integer(ik4) :: kh1 , kh2 , km1 , km2 , kl1 , kl2
    integer(ik4) :: visband
    real(rkx) :: hif , mif , lof
    !
    ! total heating rate in deg/s
    !
    do k = 1 , kz
      n = 1
      do i = ici1 , ici2
        do j = jci1 , jci2
          r2m%heatrt(j,i,k) = qrs(n,k) + qrl(n,k)
          n = n + 1
        end do
      end do
    end do
    !
    ! surface absorbed solar flux in watts/m2
    ! net up longwave flux at the surface
    !
    n = 1
    do i = ici1 , ici2
      do j = jci1 , jci2
        r2m%totcf(j,i) = totcf(n)
        r2m%solis(j,i) = sol(n)
        r2m%fsw(j,i)   = frsa(n)
        r2m%flw(j,i)   = frla(n)
        r2m%flwd(j,i)  = slwd(n)
        n = n + 1
      end do
    end do
    !
    ! for coupling with bats
    !
    ! for now assume sabveg (solar absorbed by vegetation) is equal
    ! to frsa (solar absorbed by surface). possible problems are
    ! over sparsely vegetated areas in which vegetation and ground
    ! albedo are significantly different
    !
    n = 1
    do i = ici1 , ici2
      do j = jci1 , jci2
        r2m%sabveg(j,i) = abv(n)
        r2m%solvs(j,i)  = sols(n)
        r2m%solvsd(j,i) = solsd(n)
        r2m%solvl(j,i)  = soll(n)
        r2m%solvld(j,i) = solld(n)
        r2m%sinc(j,i)   = soll(n) + sols(n) + solsd(n) + solld(n)
        n = n + 1
      end do
    end do

    if ( ifrad ) then
      rnrad_for_radfrq = rnrad_for_radfrq + 1.0_rkx
    end if
    if ( ifopt ) then
      rnrad_for_optfrq = rnrad_for_optfrq + 1.0_rkx
    end if

    if ( ifrad .and. associated(rad_higcl_out) .and. &
                     associated(rad_midcl_out) .and. &
                     associated(rad_lowcl_out) ) then
      n = 1
      kh1 = 2
      kl2 = kzp1
      do i = ici1 , ici2
        do j = jci1 , jci2
          hif = d_one
          mif = d_one
          lof = d_one
          kh2 = 2
          km1 = 2
          km2 = 2
          kl1 = 2
          do k = 2 , kzm1
            if ( m2r%phatms(j,i,k) <= 44000.0_rkx ) then
              kh2 = k
              km1 = k+1
            else if ( m2r%phatms(j,i,k) > 44000.0_rkx .and. &
                      m2r%phatms(j,i,k) <= 68000.0_rkx ) then
              km2 = k
              kl1 = k+1
            end if
          end do
          do k = kh1 , kh2
            hif = hif*(d_one - max(cld(n,k-1),cld(n,k)))/(d_one-cld(n,k-1))
          end do
          do k = km1 , km2
            mif = mif*(d_one - max(cld(n,k-1),cld(n,k)))/(d_one-cld(n,k-1))
          end do
          do k = kl1 , kl2
            lof = lof*(d_one - max(cld(n,k-1),cld(n,k)))/(d_one-cld(n,k-1))
          end do
          rad_higcl_out(j,i) = rad_higcl_out(j,i) + d_one - hif
          rad_midcl_out(j,i) = rad_midcl_out(j,i) + d_one - mif
          rad_lowcl_out(j,i) = rad_lowcl_out(j,i) + d_one - lof
          n = n + 1
        end do
      end do
    end if

    if ( rcmtimer%start( ) ) return

    if ( ifopt .and. (iaerosol == 1 .or. iclimaaer > 0) ) then
      ! when outputing aerosol properties, back to extinction m-1, ssa, and g
      if ( irrtm == 1 ) then
        visband = 10
        call copy4d_div(tauxar3d,opt_aext8_out,visband,deltaz,kth-kz+1,kth)
        ! Include the top radiation levels in integrated AOD
        ! (strato contribution) outputs
        ! Note that the stratospheric radiation hat is not visible in
        ! the oppt profiles
        call copy2d_integrate_from3(tauxar3d,opt_aod_out,visband,1,kth)
        call copy4d(tauasc3d,opt_assa8_out,visband,kth-kz+1,kth)
        call copy4d(gtota3d,opt_agfu8_out,visband,kth-kz+1,kth)
      else
        visband = 8
        call copy4d_div(tauxar3d,opt_aext8_out,visband,deltaz,1,kz)
        call copy2d_integrate_from3(tauxar3d,opt_aod_out,visband,0,kz)
        call copy4d_div2(tauasc3d,opt_assa8_out,visband,tauxar3d,1,kz)
        call copy4d_div2(gtota3d,opt_agfu8_out,visband,tauasc3d,1,kz)
      endif

      call copy4d2(deltaz,opt_deltaz_out)
      if ( idirect > 0 .or. iclimaaer > 0 ) then
        call copy2d_add(aeradfo,opt_acstoarf_out)
        call copy2d_add(aeradfos,opt_acstsrrf_out)
        if (associated(asaeradfo))  call copy2d_add(asaeradfo,opt_aastoarf_out)
        if (associated(asaeradfos)) call copy2d_add(asaeradfos,opt_aastsrrf_out)
        call copy2d_add(aerlwfo,opt_acstalrf_out)
        call copy2d_add(aerlwfos,opt_acssrlrf_out)
        if (associated(asaerlwfo))  call copy2d_add(asaerlwfo,opt_aastalrf_out)
        if (associated(asaerlwfos)) call copy2d_add(asaerlwfos,opt_aassrlrf_out)
      end if
    end if

    if ( ifrad ) then
      if ( lout ) then
        call copy3d(cld,rad_cld_out)
        call copy3d(clwp,rad_clwp_out)
        rad_clwp_out = 1.0e-3_rkx * rad_clwp_out * rad_cld_out
        if ( idiag > 0 ) then
          call copy3d(qrs,rad_qrs_out)
          call copy3d(qrl,rad_qrl_out)
          call copy3d(o3,rad_o3_out)
        end if
        if ( icosp == 1 ) then
          call copy4d1(outtaucl,rad_taucl_out,4)
          call copy4d1(outtauci,rad_tauci_out,4)
        end if
        call copy2d(frsa,rad_frsa_out)
        call copy2d(frla,rad_frla_out)
        call copy2d(clrst,rad_clrst_out)
        call copy2d(clrss,rad_clrss_out)
        call copy2d(clrlt,rad_clrlt_out)
        call copy2d(clrls,rad_clrls_out)
        call copy2d(solin,rad_solin_out)
        call copy2d(solout,rad_solout_out)
        call copy2d(totwv,rad_totwv_out)
        call copy2d(totcl,rad_totcl_out)
        call copy2d(totci,rad_totci_out)
        call copy2d(lwout,rad_lwout_out)
      end if
    end if

  end subroutine radout

  subroutine copy2d(a,b)
    implicit none
    real(rkx) , pointer , intent(in) , dimension(:) :: a
    real(rkx) , pointer , intent(inout) , dimension(:,:) :: b
    integer(ik4) :: i , j , n
    if ( associated(b) ) then
      n = 1
      do i = ici1 , ici2
        do j = jci1 , jci2
          b(j,i) = a(n)
          n = n + 1
        end do
      end do
    end if
  end subroutine copy2d

  subroutine copy2d_integrate_from3(a,b,l,ki,kl)
    implicit none
    real(rkx) , pointer , intent(in) , dimension(:,:,:) :: a
    real(rkx) , pointer , intent(inout) , dimension(:,:) :: b
    integer(ik4) , intent(in) :: l , kl , ki
    integer(ik4) :: i , j , k , n
    if ( associated(b) ) then
      b(:,:) = d_zero
      do k = ki , kl
        n = 1
        do i = ici1 , ici2
          do j = jci1 , jci2
            b(j,i) = b(j,i) + a(n,k,l)
            n = n + 1
          end do
        end do
      end do
    end if
  end subroutine copy2d_integrate_from3

  subroutine copy3d(a,b)
    implicit none
    real(rkx) , pointer , intent(in) , dimension(:,:) :: a
    real(rkx) , pointer , intent(inout) , dimension(:,:,:) :: b
    call copy3d1(a,b,1,kz)
  end subroutine copy3d

  subroutine copy3d1(a,b,k1,k2)
    implicit none
    real(rkx) , pointer , intent(in) , dimension(:,:) :: a
    real(rkx) , pointer , intent(inout) , dimension(:,:,:) :: b
    integer(ik4) , intent(in) :: k1 , k2
    integer(ik4) :: i , j , k , n
    if ( associated(b) ) then
      do k = k1 , k2
        n = 1
        do i = ici1 , ici2
          do j = jci1 , jci2
            b(j,i,k) = a(n,k)
            n = n + 1
          end do
        end do
      end do
    end if
  end subroutine copy3d1

  subroutine copy4d(a,b,l,ki,kl)
    implicit none
    real(rkx) , pointer , intent(in) , dimension(:,:,:) :: a
    real(rkx) , pointer , intent(inout) , dimension(:,:,:) :: b
    integer(ik4) , intent(in) :: l , ki , kl
    integer(ik4) :: i , j , k , n , kk
    if ( associated(b) ) then
      kk = 1
      do k = ki , kl
        n = 1
        do i = ici1 , ici2
          do j = jci1 , jci2
            b(j,i,kk) = a(n,k,l)
            n = n + 1
          end do
        end do
        kk = kk + 1
      end do
    end if
  end subroutine copy4d

  subroutine copy4d1(a,b,nl)
    implicit none
    real(rkx) , pointer , intent(in) , dimension(:,:,:) :: a
    real(rkx) , pointer , intent(inout) , dimension(:,:,:,:) :: b
    integer(ik4) , intent(in) :: nl
    integer(ik4) :: i , j , l , k , n
    if ( associated(b) ) then
      do l = 1 , nl
        do k = 1 , kz
          n = 1
          do i = ici1 , ici2
            do j = jci1 , jci2
              b(j,i,k,l) = a(n,k,l)
              n = n + 1
            end do
          end do
        end do
      end do
    end if
  end subroutine copy4d1

  subroutine copy4d2(a,b)
    implicit none
    real(rkx) , pointer , intent(in) , dimension(:,:) :: a
    real(rkx) , pointer , intent(inout) , dimension(:,:,:) :: b
    integer(ik4) :: i , j , k , n
    if ( associated(b) ) then
      do k = 1 , kz
        n = 1
        do i = ici1 , ici2
          do j = jci1 , jci2
            b(j,i,k) = a(n,k)
            n = n + 1
          end do
        end do
      end do
    end if
  end subroutine copy4d2

  subroutine copy4d_mult(a,b,l,c)
    implicit none
    real(rkx) , pointer , intent(in) , dimension(:,:,:) :: a
    real(rkx) , pointer , intent(in) , dimension(:,:) :: c
    real(rkx) , pointer , intent(inout) , dimension(:,:,:) :: b
    integer(ik4) , intent(in) :: l
    integer(ik4) :: i , j , k , n
    if ( associated(b) ) then
      do k = 1 , kz
        n = 1
        do i = ici1 , ici2
          do j = jci1 , jci2
            b(j,i,k) = a(n,k,l) * c(n,k)
            n = n + 1
          end do
        end do
      end do
    end if
  end subroutine copy4d_mult

  subroutine copy4d_div(a,b,l,c,ki,kl)
    implicit none
    real(rkx) , pointer , intent(in) , dimension(:,:,:) :: a
    real(rkx) , pointer , intent(in) , dimension(:,:) :: c
    real(rkx) , pointer , intent(inout) , dimension(:,:,:) :: b
    integer(ik4) , intent(in) :: l , ki , kl
    integer(ik4) :: i , j , k , n , kk
    if ( associated(b) ) then
      kk = 1
      do k = ki , kl
        n = 1
        do i = ici1 , ici2
          do j = jci1 , jci2
            b(j,i,kk) = a(n,k,l) / c(n,kk)
            n = n + 1
          end do
        end do
        kk = kk + 1
      end do
    end if
  end subroutine copy4d_div

  subroutine copy4d_div2(a,b,l,c,ki,kl)
    implicit none
    real(rkx) , pointer , intent(in) , dimension(:,:,:) :: a
    real(rkx) , pointer , intent(in) , dimension(:,:,:) :: c
    real(rkx) , pointer , intent(inout) , dimension(:,:,:) :: b
    integer(ik4) , intent(in) :: l , ki , kl
    integer(ik4) :: i , j , k , n , kk
    if ( associated(b) ) then
      kk = 1
      do k = ki , kl
        n = 1
        do i = ici1 , ici2
          do j = jci1 , jci2
            if ( abs(c(n,k,l)) > dlowval ) then
              b(j,i,kk) = a(n,k,l) / c(n,k,l)
            else
              b(j,i,kk) = smissval
            end if
            n = n + 1
          end do
        end do
        kk = kk + 1
      end do
    end if
  end subroutine copy4d_div2

  subroutine copy2d_add(a,b)
    implicit none
    real(rkx) , pointer , intent(in) , dimension(:) :: a
    real(rkx) , pointer , intent(inout) , dimension(:,:) :: b
    integer(ik4) :: i , j , n
    if ( associated(b) ) then
      n = 1
      do i = ici1 , ici2
        do j = jci1 , jci2
          b(j,i) = b(j,i) + a(n)
          n = n + 1
        end do
      end do
    end if
  end subroutine copy2d_add

  subroutine copy3d_add(a,b)
    implicit none
    real(rkx) , pointer , intent(in) , dimension(:,:) :: a
    real(rkx) , pointer , intent(inout) , dimension(:,:,:) :: b
    integer(ik4) :: i , j , k , n
    if ( associated(b) ) then
      do k = 1 , kz
        n = 1
        do i = ici1 , ici2
          do j = jci1 , jci2
            b(j,i,k) = b(j,i,k) + a(n,k)
            n = n + 1
          end do
        end do
      end do
    end if
  end subroutine copy3d_add

  subroutine copy4d_add(a,b,l)
    implicit none
    real(rkx) , pointer , intent(in) , dimension(:,:,:) :: a
    real(rkx) , pointer , intent(inout) , dimension(:,:,:) :: b
    integer(ik4) , intent(in) :: l
    integer(ik4) :: i , j , k , n
    if ( associated(b) ) then
      do k = 1 , kz
        n = 1
        do i = ici1 , ici2
          do j = jci1 , jci2
            b(j,i,k) = b(j,i,k) + a(n,k,l)
            n = n + 1
          end do
        end do
      end do
    end if
  end subroutine copy4d_add

end module mod_rad_outrad
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
