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
!
  subroutine radout(lout,solin,sabtp,frsa,clrst,clrss,qrs,firtp,         &
                    frla,clrlt,clrls,qrl,slwd,sols,soll,solsd,solld,     &
                    totcf,totcl,totci,cld,clwp,abv,sol,aeradfo,aeradfos, &
                    aerlwfo,aerlwfos,tauxar3d,tauasc3d,gtota3d,deltaz,   &
                    outtaucl,outtauci,r2m, & 
                       asaeradfo ,asaeradfos,asaerlwfo,asaerlwfos)

!
! copy radiation output quantities to model buffer
!
! change units of the radiative fluxes from cgs to mks
!
! compute the total radiative heat flux at the surface for
! the surface temperature computation
!
    implicit none
!
!     input/output arguments
!
! solin  - instantaneous incident solar
! sabtp  - total column absorbed solar flux
! frsa   - surface absorbed solar flux
! clrst  - clear sky total column abs solar flux
! clrss  - clear sky surface absorbed solar flux
! qrs    - solar heating rate
! firtp  - net up flux top of model (up-dwn flx)
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
    real(rk8) , pointer , dimension(:) :: clrls , clrlt ,  &
                clrss , clrst , firtp , frla , frsa ,      &
                sabtp , slwd , solin , soll , solld ,      &
                sols , solsd , totcf , totcl , totci , abv , sol
    real(rk8) , pointer , dimension(:,:) :: cld , clwp , qrl , qrs , deltaz
    real(rk8) , pointer , dimension(:,:,:) :: outtaucl , outtauci
    real(rk8) , pointer , dimension(:,:,:) :: tauxar3d , tauasc3d , gtota3d
    real(rk8) , pointer , dimension(:) :: aeradfo , aeradfos
    real(rk8) , optional, pointer , dimension(:) :: asaeradfo ,asaeradfos,asaerlwfo,asaerlwfos
    real(rk8) , pointer , dimension(:) :: aerlwfo , aerlwfos
    intent (in) cld , clrls , clrlt , clrss , clrst ,             &
                clwp , firtp , frla , frsa , qrl , qrs , sabtp ,  &
                slwd , solin , soll , solld , sols , solsd ,      &
                totcf , totcl , totci , aeradfo , aeradfos,asaeradfo , asaeradfos,       &
                aerlwfo , aerlwfos , asaerlwfo , asaerlwfos ,deltaz , outtaucl , outtauci
    type(rad_2_mod) , intent(inout) :: r2m
!
    integer(ik4) :: i , j , k , n
    integer(ik4) :: visband
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
!   surface absorbed solar flux in watts/m2
!   net up longwave flux at the surface
!
    n = 1
    do i = ici1 , ici2
      do j = jci1 , jci2
        r2m%fsw(j,i) = frsa(n)
        r2m%solis(j,i) = sol(n)
        r2m%flw(j,i) = frla(n)
        r2m%flwd(j,i) = slwd(n)
        n = n + 1
      end do
    end do
!
!   for coupling with bats
!
!   for now assume sabveg (solar absorbed by vegetation) is equal
!   to frsa (solar absorbed by surface). possible problems are
!   over sparsely vegetated areas in which vegetation and ground
!   albedo are significantly different
!
    n = 1
    do i = ici1 , ici2
      do j = jci1 , jci2
        r2m%sabveg(j,i) = abv(n)
        r2m%solvs(j,i) = sols(n)
        r2m%solvsd(j,i) = solsd(n)
        r2m%solvl(j,i) = soll(n)
        r2m%solvld(j,i) = solld(n)
        r2m%sinc(j,i) = soll(n) + sols(n) + solsd(n) + solld(n)
        n = n + 1
      end do
    end do

    if ( ktau == 0 ) return

    if ( ifchem .and. iaerosol == 1 ) then
      if(irrtm == 1) then 
        visband = 9
      else
        visband =8
      endif
      call copy4d_div(tauxar3d,opt_aext8_out,visband,deltaz)
      call copy4d_div(tauasc3d,opt_assa8_out,visband,deltaz)
      call copy4d_div(gtota3d,opt_agfu8_out,visband,deltaz)
      call copy2d_integrate_from3(tauxar3d,opt_aod_out,visband)
      if ( idirect > 0 ) then
        call copy2d_add(aeradfo,opt_acstoarf_out)
        call copy2d_add(aeradfos,opt_acstsrrf_out)
        if (present(asaeradfo))  call copy2d_add(asaeradfo,opt_aastoarf_out)
        if (present(asaeradfos)) call copy2d_add(asaeradfos,opt_aastsrrf_out)
        call copy2d_add(aerlwfo,opt_acstalrf_out)
        call copy2d_add(aerlwfos,opt_acssrlrf_out)        
        if (present(asaerlwfo))  call copy2d_add(asaerlwfo,opt_aastalrf_out)
        if (present(asaerlwfos)) call copy2d_add(asaerlwfos,opt_aassrlrf_out)

      end if
    end if
!
    if ( ifrad ) then
      if ( lout ) then
        call copy3d(cld,rad_cld_out)
        call copy3d(clwp,rad_clwp_out)
        call copy3d(qrs,rad_qrs_out)
        call copy3d(qrl,rad_qrl_out)
        call copy4d1(outtaucl,rad_taucl_out,4)
        call copy4d1(outtauci,rad_tauci_out,4)
        call copy2d(frsa,rad_frsa_out)
        call copy2d(frla,rad_frla_out)
        call copy2d(clrst,rad_clrst_out)
        call copy2d(clrss,rad_clrss_out)
        call copy2d(clrlt,rad_clrlt_out)
        call copy2d(clrls,rad_clrls_out)
        call copy2d(solin,rad_solin_out)
        call copy2d(sabtp,rad_sabtp_out)
        call copy2d(totcf,rad_totcf_out)
        call copy2d(totcl,rad_totcl_out)
        call copy2d(totci,rad_totci_out)
        call copy2d(firtp,rad_firtp_out)
      end if
    end if

  end subroutine radout

  subroutine copy2d(a,b)
    implicit none
    real(rk8) , pointer , intent(in) , dimension(:) :: a
    real(rk8) , pointer , intent(inout) , dimension(:,:) :: b
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

  subroutine copy2d_integrate_from3(a,b,l)
    implicit none
    real(rk8) , pointer , intent(in) , dimension(:,:,:) :: a
    real(rk8) , pointer , intent(inout) , dimension(:,:) :: b
    integer(ik4) , intent(in) :: l
    integer(ik4) :: i , j , k , n
    if ( associated(b) ) then
      b(:,:) = d_zero
      do k = 1 , kz
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
    real(rk8) , pointer , intent(in) , dimension(:,:) :: a
    real(rk8) , pointer , intent(inout) , dimension(:,:,:) :: b
    call copy3d1(a,b,1,kz)
  end subroutine copy3d

  subroutine copy3d1(a,b,k1,k2)
    implicit none
    real(rk8) , pointer , intent(in) , dimension(:,:) :: a
    real(rk8) , pointer , intent(inout) , dimension(:,:,:) :: b
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

  subroutine copy4d(a,b,l)
    implicit none
    real(rk8) , pointer , intent(in) , dimension(:,:,:) :: a
    real(rk8) , pointer , intent(inout) , dimension(:,:,:) :: b
    integer(ik4) , intent(in) :: l
    integer(ik4) :: i , j , k , n
    if ( associated(b) ) then
      do k = 1 , kz
        n = 1
        do i = ici1 , ici2
          do j = jci1 , jci2
            b(j,i,k) = a(n,k,l)
            n = n + 1
          end do
        end do
      end do
    end if
  end subroutine copy4d

  subroutine copy4d1(a,b,nl)
    implicit none
    real(rk8) , pointer , intent(in) , dimension(:,:,:) :: a
    real(rk8) , pointer , intent(inout) , dimension(:,:,:,:) :: b
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

  subroutine copy4d_mult(a,b,l,c)
    implicit none
    real(rk8) , pointer , intent(in) , dimension(:,:,:) :: a
    real(rk8) , pointer , intent(in) , dimension(:,:) :: c
    real(rk8) , pointer , intent(inout) , dimension(:,:,:) :: b
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

  subroutine copy4d_div(a,b,l,c)
    implicit none
    real(rk8) , pointer , intent(in) , dimension(:,:,:) :: a
    real(rk8) , pointer , intent(in) , dimension(:,:) :: c
    real(rk8) , pointer , intent(inout) , dimension(:,:,:) :: b
    integer(ik4) , intent(in) :: l
    integer(ik4) :: i , j , k , n
    if ( associated(b) ) then
      do k = 1 , kz
        n = 1
        do i = ici1 , ici2
          do j = jci1 , jci2
            b(j,i,k) = a(n,k,l) / c(n,k)
            n = n + 1
          end do
        end do
      end do
    end if
  end subroutine copy4d_div

  subroutine copy2d_add(a,b)
    implicit none
    real(rk8) , pointer , intent(in) , dimension(:) :: a
    real(rk8) , pointer , intent(inout) , dimension(:,:) :: b
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
    real(rk8) , pointer , intent(in) , dimension(:,:) :: a
    real(rk8) , pointer , intent(inout) , dimension(:,:,:) :: b
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
    real(rk8) , pointer , intent(in) , dimension(:,:,:) :: a
    real(rk8) , pointer , intent(inout) , dimension(:,:,:) :: b
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
