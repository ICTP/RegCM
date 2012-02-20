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

module mod_split
!
! Split explicit time integration
!
  use mod_dynparam
  use mod_mppparam
  use mod_runparams
  use mod_atm_interface
  use mod_vmodes
  use mod_bdycod
  use mod_atm_interface
  use mod_memutil
  use mod_service
  use mod_mppio
!
  private
!
  public :: allocate_mod_split , spinit , splitf
!
  real(dp) , pointer , dimension(:) :: aam
  real(dp) , pointer , dimension(:) :: an
  real(dp) , pointer , dimension(:,:) :: am
  real(dp) , pointer , dimension(:,:,:) :: uuu , vvv
!
  real(dp) , pointer , dimension(:,:,:) :: ddsum
  real(dp) , pointer , dimension(:,:,:) :: dhsum
  real(dp) , pointer , dimension(:,:,:,:) :: deld
  real(dp) , pointer , dimension(:,:,:,:) :: delh
  real(dp) , pointer , dimension(:,:,:) :: work
  real(dp) , pointer , dimension(:,:) :: uu , vv
  real(dp) , pointer , dimension(:,:) :: xdelh
!
  contains 
!
  subroutine allocate_mod_split
    implicit none
    character (len=64) :: subroutine_name='allocate_mod_split'
    integer :: idindx = 0
!
    call time_begin(subroutine_name,idindx)
    call getmem1d(aam,1,nsplit,'split:aam')
    call getmem2d(am,1,kz,1,nsplit,'split:am')
    call getmem1d(an,1,nsplit,'split:naam')
    call getmem3d(ddsum,1,jxp,idot1,idot2,1,nsplit,'split:ddsum')
    call getmem4d(deld,1,jxp,idot1,idot2,1,nsplit,1,3,'split:deld')
    call getmem4d(delh,1,jxp,idot1,idot2,1,nsplit,1,3,'split:delh')
    call getmem2d(xdelh,0,jxp,idot1,idot2,'split:xdelh')
    call getmem3d(dhsum,0,jxp,idot1,idot2,1,nsplit,'split:dhsum')
    call getmem3d(work,1,jxp,idot1,idot2,1,3,'split:work')
    call getmem2d(uu,1,jxp+1,idot1,idot2,'split:uu')
    call getmem2d(vv,1,jxp+1,idot1,idot2,'split:vv')
    call getmem3d(uuu,1,jxp+1,idot1,idot2,1,kz,'split:uuu')
    call getmem3d(vvv,1,jxp+1,idot1,idot2,1,kz,'split:vvv')
    call time_end(subroutine_name,idindx)
  end subroutine allocate_mod_split
!
! Intial computation of vertical modes.
!
  subroutine spinit
#ifndef IBM
    use mpi
#endif
    implicit none
#ifdef IBM
    include 'mpif.h'
#endif
    real(dp) :: eps , eps1 , fac , pdlog
    integer :: i , ijlx , j , k , l , n , ns
    logical :: lstand
    character (len=64) :: subroutine_name='spinit'
    integer :: idindx=0
!
    call time_begin(subroutine_name,idindx)
!
!   compute m.
!
    do ns = 1 , nsplit
      aam(ns) = dnint(dt/dtau(ns))
      if ( ktau /= 0 ) then
        aam(ns) = dnint((dt*d_half)/dtau(ns))
      end if
    end do
    if ( myid == 0 ) print * , 'dt, dtau = ' , dt , dtau
!
!   compute xps and tbarh for use in vmodes.
!
    call allocate_mod_vmodes
!
    xps = d_zero
    ijlx = (jce2-jce1+1)*(ice2-ice1+1)
    do k = 1 , kz
      tbarh(k) = d_zero
    end do
    do i = ice1 , ice2
      do j = jce1 , jce2
        xps = xps + sfs%psa(j,i)/ijlx
      end do
    end do

    do k = 1 , kz
      do i = ice1 , ice2
        do j = jce1 , jce2
          tbarh(k) = tbarh(k) + atm1%t(j,i,k)/(sfs%psa(j,i)*ijlx)
        end do
      end do
    end do
!
!   compute vertical modes.
!   lstand = .true. if standard atmosphere t to be used (ignore input
!   tbarh and xps in that case). otherwise, xps and tbarh must
!   be defined on input.
!   dtau = time steps(in sec)for modes in split explicit is
!   specified in namelist as array dtsplit
!
    lstand = .true.
    if ( ktau /= 0 ) lstand = .true.
    call vmodes(lstand)
!
!   compute am and an.
!
    do n = 1 , nsplit
      an(n) = d_zero
      do l = 1 , kz
        an(n) = an(n) + dsigma(l)*zmatx(l,n)
      end do
      do k = 1 , kz
        am(k,n) = d_zero
        tau(n,k) = d_zero
      end do
      do l = 1 , kz
        do k = 1 , kz
          am(k,n) = am(k,n) + a0(k,l)*zmatx(l,n)
          tau(n,k) = tau(n,k) + rgas*zmatxr(n,l)*hydros(l,k)
        end do
      end do
!
      do k = 1 , kzp1
        varpa1(n,k) = d_zero
      end do
      do l = 1 , kz
        do k = 1 , kzp1
          varpa1(n,k) = varpa1(n,k) + rgas*zmatxr(n,l)*hydroc(l,k)
        end do
      end do
    end do
!
!   multiply am, an and zmatx by factor.
!
    do l = 1 , nsplit
      fac = d_two*dt/(d_two*aam(l)+d_one)
      if ( ktau /= 0 ) then
        fac = dt/(d_two*aam(l)+d_one)
      end if
      if ( myid == 0 ) print * , 'aam, fac = ' , aam(l) , fac
      an(l) = an(l)*fac
      do k = 1 , kz
        zmatx(k,l) = zmatx(k,l)*fac
        am(k,l) = am(k,l)*fac
      end do
    end do
!
    if ( ifrest ) return
!
!=======================================================================
!
!   zero new arrays
!
    dstor(:,:,:) = d_zero
    hstor(:,:,:) = d_zero
!
!   Divergence manipulations (0)
!
!   compute divergence z from u and v
!   ( u must be pstar * u ; similarly for v )
!   ( note: map scale factors have been inverted in model (init) )
!
    do k = 1 , kz
      do i = ide1 , ide2
        do j = jde1 , jde2
          uuu(j,i,k) = atm2%u(j,i,k)*mddom%msfd(j,i)
          vvv(j,i,k) = atm2%v(j,i,k)*mddom%msfd(j,i)
        end do
      end do
    end do
!
    call deco1_exchange_right(uuu,1,ide1,ide2,1,kz)
    call deco1_exchange_right(vvv,1,ide1,ide2,1,kz)
!
    dstor(:,:,:) = d_zero
!
    do l = 1 , nsplit
      do k = 1 , kz
        do i = ice1 , ice2
          do j = jce1 , jce2
            fac = dx2*mddom%msfx(j,i)*mddom%msfx(j,i)
            dstor(j,i,l) = dstor(j,i,l) + zmatxr(l,k) * & 
                 (-uuu(j,i+1,k)+uuu(j+1,i+1,k)-uuu(j,i,k)+uuu(j+1,i,k) + &
                   vvv(j,i+1,k)+vvv(j+1,i+1,k)-vvv(j,i,k)-vvv(j+1,i,k))/fac
          end do
        end do
      end do
    end do
!
!=======================================================================
!
!     Geopotential manipulations
!
    do l = 1 , nsplit
      pdlog = varpa1(l,kzp1)*dlog(sigmah(kzp1)*pd+ptop)
      eps1 = varpa1(l,kzp1)*sigmah(kzp1)/(sigmah(kzp1)*pd+ptop)
      do i = ice1 , ice2
        do j = jce1 , jce2
          eps = eps1*(sfs%psb(j,i)-pd)
          hstor(j,i,l) = pdlog + eps
        end do
      end do

      do k = 1 , kz
        pdlog = varpa1(l,k)*dlog(sigmah(k)*pd+ptop)
        eps1 = varpa1(l,k)*sigmah(k)/(sigmah(k)*pd+ptop)
        do i = ice1 , ice2
          do j = jce1 , jce2
            eps = eps1*(sfs%psb(j,i)-pd)
            hstor(j,i,l) = hstor(j,i,l) + pdlog + &
                           tau(l,k)*atm2%t(j,i,k)/sfs%psb(j,i) + eps
          end do
        end do
      end do
    end do
!
    call time_end(subroutine_name,idindx)
!
  end subroutine spinit
!
! Compute deld, delh, integrate in time and add correction terms appropriately
!
  subroutine splitf
#ifndef IBM
    use mpi
#endif
    implicit none
#ifdef IBM
    include 'mpif.h'
#endif
    real(dp) :: eps , eps1 , fac , gnuam , gnuan , gnuzm , pdlog , x , y
    integer :: i , j , k , l , n
    character (len=64) :: subroutine_name='splitf'
    integer :: idindx=0
!
    call time_begin(subroutine_name,idindx)
!
    deld(:,:,:,:) = d_zero
    delh(:,:,:,:) = d_zero
!
!   compute pressure on dot grid
!=======================================================================
!
!   this routine determines p(.) from p(x) by a 4-point interpolation.
!   on the x-grid, a p(x) point outside the grid domain is assumed to
!   satisfy p(j,0)=p(j,1); p(j,iy)=p(j,iym1); and similarly for the i's.
!
    call deco1_exchange_left(sfs%psa,1,ice1,ice2)
    call psc2psd(sfs%psa,psdot)
!
!=======================================================================
!
!   get deld(0), delh(0) from storage
!
   do n = 1 , nsplit
     do i = ide1 , ide2
       do j = jde1 , jde2
         deld(j,i,n,1) = dstor(j,i,n)
         delh(j,i,n,1) = hstor(j,i,n)
       end do
     end do
   end do
!
!=======================================================================
!
!   Divergence manipulations (f)
!
    do k = 1 , kz
      do i = ide1 , ide2
        do j = jde1 , jde2
          uuu(j,i,k) = atm1%u(j,i,k)*mddom%msfd(j,i)
          vvv(j,i,k) = atm1%v(j,i,k)*mddom%msfd(j,i)
        end do
      end do
    end do
!
    call deco1_exchange_right(uuu,1,1,iy,1,kz)
    call deco1_exchange_right(vvv,1,1,iy,1,kz)
!
    do l = 1 , nsplit
      do i = ide1 , ide2
        do j = jde1 , jde2
          deld(j,i,l,3) = d_zero
        end do
      end do

      do k = 1 , kz
        do i = ice1 , ice2
          do j = jce1 , jce2
            fac = dx2*mddom%msfx(j,i)*mddom%msfx(j,i)
            deld(j,i,l,3) = deld(j,i,l,3) + zmatxr(l,k) * &
               (-uuu(j,i+1,k)+uuu(j+1,i+1,k)-uuu(j,i,k)+uuu(j+1,i,k) + &
                 vvv(j,i+1,k)+vvv(j+1,i+1,k)-vvv(j,i,k)-vvv(j+1,i,k))/fac
          end do
        end do
      end do
    end do

!
!=======================================================================
! 
    do n = 1 , nsplit
      do i = ide1 , ide2
        do j = jde1 , jde2
          deld(j,i,n,3) = deld(j,i,n,3) - deld(j,i,n,1)
        end do
      end do
    end do
!
!=======================================================================
!
!   Divergence manipulations (0)
!
    do k = 1 , kz
      do i = ide1 , ide2
        do j = jde1 , jde2
          uuu(j,i,k) = atm2%u(j,i,k)*mddom%msfd(j,i)
          vvv(j,i,k) = atm2%v(j,i,k)*mddom%msfd(j,i)
        end do
      end do
    end do
!
    call deco1_exchange_right(uuu,1,ide1,ide2,1,kz)
    call deco1_exchange_right(vvv,1,ide1,ide2,1,kz)
!
    do l = 1 , nsplit
      do i = ide1 , ide2
        do j = jde1 , jde2
          deld(j,i,l,2) = d_zero
        end do
      end do
      do k = 1 , kz
        do i = ice1 , ice2
          do j = jce1 , jce2
            fac = dx2*mddom%msfx(j,i)*mddom%msfx(j,i)
            deld(j,i,l,2) = deld(j,i,l,2) + zmatxr(l,k) * &
              (-uuu(j,i+1,k)+uuu(j+1,i+1,k)-uuu(j,i,k)+uuu(j+1,i,k) + &
                vvv(j,i+1,k)+vvv(j+1,i+1,k)-vvv(j,i,k)-vvv(j+1,i,k))/fac
          end do
        end do
      end do
    end do
!
!=======================================================================
!
    do n = 1 , nsplit
      do i = ide1 , ide2
        do j = jde1 , jde2
          deld(j,i,n,1) = deld(j,i,n,1) - deld(j,i,n,2)
        end do
      end do
    end do
!
!=======================================================================
!
!   Geopotential manipulations (f)
!
    do l = 1 , nsplit
      pdlog = varpa1(l,kzp1)*dlog(sigmah(kzp1)*pd+ptop)
      eps1 = varpa1(l,kzp1)*sigmah(kzp1)/(sigmah(kzp1)*pd+ptop)
      do i = ice1 , ice2
        do j = jce1 , jce2
          eps = eps1*(sfs%psa(j,i)-pd)
          delh(j,i,l,3) = pdlog + eps
        end do
      end do
      do k = 1 , kz
        pdlog = varpa1(l,k)*dlog(sigmah(k)*pd+ptop)
        eps1 = varpa1(l,k)*sigmah(k)/(sigmah(k)*pd+ptop)
        do i = ice1 , ice2
          do j = jce1 , jce2
            eps = eps1*(sfs%psa(j,i)-pd)
            delh(j,i,l,3) = delh(j,i,l,3) + pdlog +  &
                    tau(l,k)*atm1%t(j,i,k)/sfs%psa(j,i) + eps
          end do
        end do
      end do
    end do
!
!=======================================================================
! 
    do n = 1 , nsplit
      do i = ide1 , ide2
        do j = jde1 , jde2
          delh(j,i,n,3) = delh(j,i,n,3) - delh(j,i,n,1)
        end do
      end do
    end do
!
!=======================================================================
!
!   Geopotential manipulations (0)
!
    do l = 1 , nsplit
      pdlog = varpa1(l,kzp1)*dlog(sigmah(kzp1)*pd+ptop)
      eps1 = varpa1(l,kzp1)*sigmah(kzp1)/(sigmah(kzp1)*pd+ptop)
      do i = ice1 , ice2
        do j = jce1 , jce2
          eps = eps1*(sfs%psb(j,i)-pd)
          delh(j,i,l,2) = pdlog + eps
        end do
      end do
      do k = 1 , kz
        pdlog = varpa1(l,k)*dlog(sigmah(k)*pd+ptop)
        eps1 = varpa1(l,k)*sigmah(k)/(sigmah(k)*pd+ptop)
        do i = ice1 , ice2
          do j = jce1 , jce2
            eps = eps1*(sfs%psb(j,i)-pd)
            delh(j,i,l,2) = delh(j,i,l,2) + pdlog +  &
                     tau(l,k)*atm2%t(j,i,k)/sfs%psb(j,i) + eps
          end do
        end do
      end do
    end do
!
!=======================================================================
!
    do n = 1 , nsplit
      do i = ide1 , ide2
        do j = jde1 , jde2
          delh(j,i,n,1) = delh(j,i,n,1) - delh(j,i,n,2)
        end do
      end do
    end do
!
!   put deld(0), delh(0) into storage
!
   do n = 1 , nsplit
     do i = ide1 , ide2
       do j = jde1 , jde2
         dstor(j,i,n) = deld(j,i,n,2)
         hstor(j,i,n) = delh(j,i,n,2)
       end do
     end do
   end do
!
!   split explicit time integration
!
    call spstep
!
!=======================================================================
!
!   Add corrections to t and p;  u and v
!
    do l = 1 , nsplit
      gnuan = gnuhf*an(l)
      do i = ici1 , ici2
        do j = jci1 , jci2
          sfs%psa(j,i) = sfs%psa(j,i) - an(l)*ddsum(j,i,l)
          sfs%psb(j,i) = sfs%psb(j,i) - gnuan*ddsum(j,i,l)
        end do
      end do
    end do
    do l = 1 , nsplit
      do k = 1 , kz
        gnuam = gnuhf*am(k,l)
        do i = ici1 , ici2
          do j = jci1 , jci2
            atm1%t(j,i,k) = atm1%t(j,i,k) + am(k,l)*ddsum(j,i,l)
            atm2%t(j,i,k) = atm2%t(j,i,k) + gnuam*ddsum(j,i,l)
          end do
        end do
      end do
    end do

!=======================================================================

    call deco1_exchange_left(dhsum,1,1,iy,1,nsplit)

    do l = 1 , nsplit
      do k = 1 , kz
        gnuzm = gnuhf*zmatx(k,l)
        do i = idi1 , idi2
          do j = jdi1 , jdi2
            fac = psdot(j,i)/(dx2*mddom%msfd(j,i))
            x = fac*(dhsum(j,i,l)+dhsum(j,i-1,l) - &
                     dhsum(j-1,i,l)-dhsum(j-1,i-1,l))
            y = fac*(dhsum(j,i,l)-dhsum(j,i-1,l) + &
                     dhsum(j-1,i,l)-dhsum(j-1,i-1,l))
            atm1%u(j,i,k) = atm1%u(j,i,k) - zmatx(k,l)*x
            atm1%v(j,i,k) = atm1%v(j,i,k) - zmatx(k,l)*y
            atm2%u(j,i,k) = atm2%u(j,i,k) - gnuzm*x
            atm2%v(j,i,k) = atm2%v(j,i,k) - gnuzm*y
          end do
        end do
      end do
    end do
!
!=======================================================================
!
    call time_end(subroutine_name,idindx)
!
  end subroutine splitf
!
  subroutine spstep
    implicit none
!
    real(dp) :: dtau2 , fac
    integer :: i , j , m2 , n , n0 , n1 , n2 , ns , nw
    character (len=64) :: subroutine_name='spstep'
    integer :: idindx=0
!
    call time_begin(subroutine_name,idindx)
!
    do n = 1 , nsplit
      do i = ide1 , ide2
        do j = jde1 , jde2
          ddsum(j,i,n) = d_zero
          dhsum(j,i,n) = d_zero
        end do
      end do
    end do
!
    do ns = 1 , nsplit
!
      n0 = 1
      n1 = 2
      n2 = n0
      m2 = idint(aam(ns))*2
      dtau2 = dtau(ns)*d_two
!
!     below follows Madala (1987)
!
      do i = ice1 , ice2
        do j = jce1 , jce2
!         deld, delh: 1,ilx on cross grid
          ddsum(j,i,ns) = deld(j,i,ns,n0)
          dhsum(j,i,ns) = delh(j,i,ns,n0)
        end do
      end do
!
!     first step, use forward scheme
!
!=======================================================================
!
!     compute gradient of delh;  output = (work1,work2)
!
      xdelh(jdi1:jdi2,idi1:idi2) = delh(jdi1:jdi2,idi1:idi2,ns,n0)
      call deco1_exchange_left(xdelh,1,1,iy)
      do i = idi1 , idi2
        do j = jdi1 , jdi2
          fac = dx2*mddom%msfx(j,i)
          work(j,i,1) = (xdelh(j,i)  +xdelh(j,i-1) - &
                         xdelh(j-1,i)-xdelh(j-1,i-1))/fac
          work(j,i,2) = (xdelh(j,i)  +xdelh(j-1,i) - &
                         xdelh(j,i-1)-xdelh(j-1,i-1))/fac
        end do
      end do
!
!=======================================================================
!
      do nw = 1 , 2
        do i = idi1 , idi2
          do j = jdi1 , jdi2
!           work: 2,ilx on dot grid
            work(j,i,nw) = work(j,i,nw)*psdot(j,i)
          end do
        end do
      end do
!
!=======================================================================
!
!     compute divergence z from u and v
!     ( u must be pstar * u ; similarly for v )
!     ( note: map scale factors have been inverted in model (init) )
!
      do i = idi1 , idi2
        do j = jdi1 , jdi2
          uu(j,i) = work(j,i,1)*mddom%msfd(j,i)
          vv(j,i) = work(j,i,2)*mddom%msfd(j,i)
        end do
      end do
!
      call deco1_exchange_right(uu,1,1,iy)
      call deco1_exchange_right(vv,1,1,iy)
!
      do i = ici1 , ici2
        do j = jci1 , jci2
          fac = dx2*mddom%msfx(j,i)*mddom%msfx(j,i)
          work(j,i,3) = (-uu(j,i+1)+uu(j+1,i+1)-uu(j,i)+uu(j+1,i) &
                         +vv(j,i+1)+vv(j+1,i+1)-vv(j,i)-vv(j+1,i))/fac
        end do
      end do
!
!=======================================================================
!
      do i = ici1 , ici2
        do j = jci1 , jci2
          deld(j,i,ns,n1) = deld(j,i,ns,n0) - dtau(ns)*work(j,i,3) + &
                            deld(j,i,ns,3)/m2
          delh(j,i,ns,n1) = delh(j,i,ns,n0) - dtau(ns)*hbar(ns) * &
                            deld(j,i,ns,n0)/sfs%psa(j,i)+delh(j,i,ns,3)/m2
        end do
      end do
 
!     not in Madala (1987)
      fac = (aam(ns)-d_one)/aam(ns)
      if ( ma%hasleft ) then
        do i = ici1 , ici2
          delh(jce1,i,ns,n1) = delh(jce1,i,ns,n0)*fac
        end do
      end if
      if ( ma%hasright ) then
        do i = ici1 , ici2
          delh(jce2,i,ns,n1) = delh(jce2,i,ns,n0)*fac
        end do
      end if
      if ( ma%hasbottom ) then
        do j = jce1 , jce2
          delh(j,ice1,ns,n1) = delh(j,ice1,ns,n0)*fac
        end do
      end if
      if ( ma%hastop ) then
        do j = jce1 , jce2
          delh(j,ice2,ns,n1) = delh(j,ice2,ns,n0)*fac
        end do
      end if
!
      do i = ice1 , ice2
        do j = jce1 , jce2
          ddsum(j,i,ns) = ddsum(j,i,ns) + deld(j,i,ns,n1)
          dhsum(j,i,ns) = dhsum(j,i,ns) + delh(j,i,ns,n1)
        end do
      end do
!
!     subsequent steps, use leapfrog scheme
!
      do n = 2 , m2
!
!=======================================================================
!
!       compute gradient of delh;  output = (work1,work2)
!
        xdelh(jdi1:jdi2,idi1:idi2) = delh(jdi1:jdi2,idi1:idi2,ns,n1)
        call deco1_exchange_left(xdelh,1,1,iy)
        do i = idi1 , idi2
          do j = jdi1 , jdi2
            fac = dx2*mddom%msfx(j,i)
            work(j,i,1) = (xdelh(j,i)+xdelh(j,i-1)- &
                           xdelh(j-1,i)-xdelh(j-1,i-1))/fac
            work(j,i,2) = (xdelh(j,i)+xdelh(j-1,i)- &
                           xdelh(j,i-1)-xdelh(j-1,i-1))/fac
          end do
        end do
!
!=======================================================================
!
        do nw = 1 , 2
          do i = idi1 , idi2
            do j = jdi1 , jdi2
              work(j,i,nw) = work(j,i,nw)*psdot(j,i)
            end do
          end do
        end do
!
!=======================================================================
!
!       compute divergence z from u and v
!       ( u must be pstar * u ; similarly for v )
!       ( note: map scale factors have been inverted in model (init) )
!
        do i = idi1 , idi2
          do j = jdi1 , jdi2
            uu(j,i) = work(j,i,1)*mddom%msfd(j,i)
            vv(j,i) = work(j,i,2)*mddom%msfd(j,i)
          end do
        end do
!
        call deco1_exchange_right(uu,1,1,iy)
        call deco1_exchange_right(vv,1,1,iy)
!
        do i = ici1 , ici2
          do j = jci1 , jci2
            fac = dx2*mddom%msfx(j,i)*mddom%msfx(j,i)
            work(j,i,3) = (-uu(j,i+1)+uu(j+1,i+1)-uu(j,i)+uu(j+1,i) + &
                            vv(j,i+1)+vv(j+1,i+1)-vv(j,i)-vv(j+1,i))/fac
          end do
        end do
!
!=======================================================================
!
        do i = ici1 , ici2
          do j = jci1 , jci2
            deld(j,i,ns,n2) = deld(j,i,ns,n0) - dtau2*work(j,i,3) + &
                              deld(j,i,ns,3)/aam(ns)
            delh(j,i,ns,n2) = delh(j,i,ns,n0) - dtau2*hbar(ns) * &
                              deld(j,i,ns,n1)/sfs%psa(j,i) +     &
                              delh(j,i,ns,3)/aam(ns)
          end do
        end do
!
!       not in Madala (1987)
!
        if ( ma%hasleft ) then 
          do i = ici1 , ici2
            delh(jce1,i,ns,n2) = d_two*delh(jce1,i,ns,n1)-delh(jce1,i,ns,n0)
          end do
        end if
        if ( ma%hasright ) then
          do i = ici1 , ici2
            delh(jce2,i,ns,n2) = d_two*delh(jce2,i,ns,n1)-delh(jce2,i,ns,n0)
          end do
        end if
        if ( ma%hasbottom ) then
          do j = jce1 , jce2
            delh(j,ice1,ns,n2) = d_two*delh(j,ice1,ns,n1)-delh(j,ice1,ns,n0)
          end do
        end if
        if ( ma%hastop ) then
          do j = jce1 , jce2
            delh(j,ice2,ns,n2) = d_two*delh(j,ice2,ns,n1)-delh(j,ice2,ns,n0)
          end do
        end if
!
        do i = ice1 , ice2
          do j = jce1 , jce2
            ddsum(j,i,ns) = ddsum(j,i,ns) + deld(j,i,ns,n2)
            dhsum(j,i,ns) = dhsum(j,i,ns) + delh(j,i,ns,n2)
          end do
        end do
!
        n0 = n1
        n1 = n2
        n2 = n0
      end do
!
    end do
!
    call time_end(subroutine_name,idindx)
  end subroutine spstep
!
end module mod_split
