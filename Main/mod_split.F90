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
  use mod_runparams
  use mod_atm_interface
  use mod_vmodes
  use mod_bdycod
  use mod_atm_interface
  use mod_savefile
  use mod_memutil
  use mod_service
  use mod_mppio
!
  private
!
  public :: allocate_mod_split , spinit , splitf
!
  real(8) , pointer , dimension(:) :: aam
  real(8) , pointer , dimension(:) :: an
  real(8) , pointer , dimension(:,:) :: am
  real(8) , pointer , dimension(:,:,:) :: uuu , vvv
!
  real(8) , pointer , dimension(:,:,:) :: ddsum
  real(8) , pointer , dimension(:,:,:) :: dhsum
  real(8) , pointer , dimension(:,:,:,:) :: deld
  real(8) , pointer , dimension(:,:,:,:) :: delh
  real(8) , pointer , dimension(:,:) :: psdot
  real(8) , pointer , dimension(:,:,:) :: work
  real(8) , pointer , dimension(:,:) :: uu , vv
  real(8) , pointer , dimension(:) :: trans1 , trans2
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
    call getmem1d(trans1,1,iy,'split:trans1')
    call getmem1d(trans2,1,iy,'split:trans2')
    call getmem2d(am,1,kz,1,nsplit,'split:am')
    call getmem1d(an,1,nsplit,'split:naam')
    call getmem3d(ddsum,1,jxp,1,iy,1,nsplit,'split:ddsum')
    call getmem4d(deld,1,jxp,1,iy,1,nsplit,1,3,'split:deld')
    call getmem4d(delh,0,jxp,1,iy,1,nsplit,1,3,'split:delh')
    call getmem3d(dhsum,0,jxp,1,iy,1,nsplit,'split:dhsum')
    call getmem2d(psdot,0,jxp,1,iy,'split:psdot')
    call getmem3d(work,1,jxp,1,iy,1,3,'split:work')
    call getmem2d(uu,1,jxp+1,1,iy,'split:uu')
    call getmem2d(vv,1,jxp+1,1,iy,'split:vv')
    call getmem3d(uuu,1,jxp+1,1,iy,1,kz,'split:uuu')
    call getmem3d(vvv,1,jxp+1,1,iy,1,kz,'split:vvv')
    call time_end(subroutine_name,idindx)
  end subroutine allocate_mod_split
!
! Intial computation of vertical modes.
!
  subroutine spinit(jstart,jend,istart,iend)
#ifndef IBM
    use mpi
#endif
    implicit none
#ifdef IBM
    include 'mpif.h'
#endif
!
    integer , intent(in) :: jstart , jend , istart , iend
    real(8) :: eps , eps1 , fac , pdlog
    integer :: i , ijlx , j , k , l , n , ns
    logical :: lstand
    integer :: ierr
    integer :: jp1
    character (len=64) :: subroutine_name='spinit'
    integer :: idindx=0
!
    call time_begin(subroutine_name,idindx)
!
!   zero new arrays
!
    dstor(:,:,:) = d_zero
    hstor(:,:,:) = d_zero
!
    ijlx = jendx*iym1
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
    do k = 1 , kz
      tbarh(k) = d_zero
    end do
    do i = 1 , iym1
      do j = 1 , jendx
        xps = xps + sps1%ps(j,i)/ijlx
      end do
    end do

    do k = 1 , kz
      do i = 1 , iym1
        do j = 1 , jendx
          tbarh(k) = tbarh(k) + atm1%t(i,k,j)/(sps1%ps(j,i)*ijlx)
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
    if ( ifrest ) then
      call read_savefile_part2
!
      if ( myid == 0 ) then
        do j = 1 , jx
          do n = 1 , nsplit
            do i = 1 , iy
              sav_0d(i,n,j) = dstor_io(i,j,n)
              sav_0d(i,n+nsplit,j) = hstor_io(i,j,n)
            end do
          end do
        end do
        do j = 1 , jx
          do k = 1 , kz
            sav_6(k,1,j) = ui1_io(k,j)
            sav_6(k,2,j) = ui2_io(k,j)
            sav_6(k,3,j) = uilx_io(k,j)
            sav_6(k,4,j) = uil_io(k,j)
            sav_6(k,5,j) = vi1_io(k,j)
            sav_6(k,6,j) = vi2_io(k,j)
            sav_6(k,7,j) = vilx_io(k,j)
            sav_6(k,8,j) = vil_io(k,j)
          end do
        end do
      end if
      call mpi_scatter(sav_0d,iy*nsplit*2*jxp,mpi_real8,       &
                       sav0d, iy*nsplit*2*jxp,mpi_real8,       &
                       0,mycomm,ierr)
      do j = 1 , jendl
        do n = 1 , nsplit
          do i = 1 , iy
            dstor(j,i,n) = sav0d(i,n,j)
            hstor(j,i,n) = sav0d(i,n+nsplit,j)
          end do
        end do
      end do
      call mpi_scatter(sav_6,kz*8*jxp,mpi_real8, &
                       sav6, kz*8*jxp,mpi_real8, &
                       0,mycomm,ierr)
      do j = 1 , jendl
        do k = 1 , kz
          ui1(k,j) = sav6(k,1,j)
          ui2(k,j) = sav6(k,2,j)
          uilx(k,j) = sav6(k,3,j)
          uil(k,j) = sav6(k,4,j)
          vi1(k,j) = sav6(k,5,j)
          vi2(k,j) = sav6(k,6,j)
          vilx(k,j) = sav6(k,7,j)
          vil(k,j) = sav6(k,8,j)
        end do
      end do
#ifndef BAND
      call mpi_bcast(uj1,iy*kz,mpi_real8,0,mycomm,ierr)
      call mpi_bcast(uj2,iy*kz,mpi_real8,0,mycomm,ierr)
      call mpi_bcast(vj1,iy*kz,mpi_real8,0,mycomm,ierr)
      call mpi_bcast(vj2,iy*kz,mpi_real8,0,mycomm,ierr)
      call mpi_bcast(ujlx,iy*kz,mpi_real8,0,mycomm,ierr)
      call mpi_bcast(ujl,iy*kz,mpi_real8,0,mycomm,ierr)
      call mpi_bcast(vjlx,iy*kz,mpi_real8,0,mycomm,ierr)
      call mpi_bcast(vjl,iy*kz,mpi_real8,0,mycomm,ierr)
      if ( myid /= nproc-1 ) then
#endif
        do k = 1 , kz
          var1snd(k,1) = ui1(k,jxp)
          var1snd(k,2) = vi1(k,jxp)
          var1snd(k,3) = ui2(k,jxp)
          var1snd(k,4) = vi2(k,jxp)
          var1snd(k,5) = uilx(k,jxp)
          var1snd(k,6) = vilx(k,jxp)
          var1snd(k,7) = uil(k,jxp)
          var1snd(k,8) = vil(k,jxp)
        end do
#ifndef BAND
      end if
#endif
      call mpi_sendrecv(var1snd,kz*8,mpi_real8,ieast,1, &
                        var1rcv,kz*8,mpi_real8,iwest,1, &
                        mycomm,mpi_status_ignore,ierr)
#ifndef BAND
      if ( myid /= 0 ) then
#endif
        do k = 1 , kz
          ui1(k,0) = var1rcv(k,1)
          vi1(k,0) = var1rcv(k,2)
          ui2(k,0) = var1rcv(k,3)
          vi2(k,0) = var1rcv(k,4)
          uilx(k,0) = var1rcv(k,5)
          vilx(k,0) = var1rcv(k,6)
          uil(k,0) = var1rcv(k,7)
          vil(k,0) = var1rcv(k,8)
        end do
#ifndef BAND
      end if
      if ( myid /= 0 ) then
#endif
        do k = 1 , kz
          var1snd(k,1) = ui1(k,1)
          var1snd(k,2) = vi1(k,1)
          var1snd(k,3) = ui2(k,1)
          var1snd(k,4) = vi2(k,1)
          var1snd(k,5) = uilx(k,1)
          var1snd(k,6) = vilx(k,1)
          var1snd(k,7) = uil(k,1)
          var1snd(k,8) = vil(k,1)
        end do
#ifndef BAND
      end if
#endif
      call mpi_sendrecv(var1snd,kz*8,mpi_real8,iwest,2, &
                        var1rcv,kz*8,mpi_real8,ieast,2, &
                        mycomm,mpi_status_ignore,ierr)
#ifndef BAND
      if ( myid /= nproc-1 ) then
#endif
        do k = 1 , kz
          ui1(k,jxp+1) = var1rcv(k,1)
          vi1(k,jxp+1) = var1rcv(k,2)
          ui2(k,jxp+1) = var1rcv(k,3)
          vi2(k,jxp+1) = var1rcv(k,4)
          uilx(k,jxp+1) = var1rcv(k,5)
          vilx(k,jxp+1) = var1rcv(k,6)
          uil(k,jxp+1) = var1rcv(k,7)
          vil(k,jxp+1) = var1rcv(k,8)
        end do
#ifndef BAND
      end if
#endif
    else
!
!=======================================================================
!     Divergence manipulations (0)
!
!     compute divergence z from u and v
!     ( u must be pstar * u ; similarly for v )
!     ( note: map scale factors have been inverted in model (init) )
!
      do k = 1 , kz
        do i = 1 , iy
          do j = 1 , jendl
            uuu(j,i,k) = atm2%u(i,k,j)*mddom%msfd(j,i)
            vvv(j,i,k) = atm2%v(i,k,j)*mddom%msfd(j,i)
          end do
        end do
        if ( iwest /= mpi_proc_null ) then
          trans1 = uuu(1,:,k)
        end if
        call mpi_sendrecv(trans1,iy,mpi_real8,iwest,2, &
                          trans2,iy,mpi_real8,ieast,2, &
                          mycomm,mpi_status_ignore,ierr)
        if ( ieast /= mpi_proc_null ) then
          uuu(jxp+1,:,k) = trans2
        end if
        if ( iwest /= mpi_proc_null ) then
          trans1 = vvv(1,:,k)
        end if
        call mpi_sendrecv(trans1,iy,mpi_real8,iwest,2, &
                          trans2,iy,mpi_real8,ieast,2, &
                          mycomm,mpi_status_ignore,ierr)
        if ( ieast /= mpi_proc_null ) then
          vvv(jxp+1,:,k) = trans2
        end if
      end do
!
      dstor(:,:,:) = d_zero
      do l = 1 , nsplit
        do k = 1 , kz
          do i = 1 , iym1
            do j = 1 , jendx
              jp1 = j+1
              fac = dx2*mddom%msfx(j,i)*mddom%msfx(j,i)
              dstor(j,i,l) = dstor(j,i,l) + zmatxr(l,k) * & 
                   (-uuu(j,i+1,k)+uuu(jp1,i+1,k)-uuu(j,i,k)+uuu(jp1,i,k) + &
                     vvv(j,i+1,k)+vvv(jp1,i+1,k)-vvv(j,i,k)-vvv(jp1,i,k))/fac
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
        do i = 1 , iym1
          do j = 1 , jendx
            eps = eps1*(sps2%ps(j,i)-pd)
            hstor(j,i,l) = pdlog + eps
          end do
        end do

        do k = 1 , kz
          pdlog = varpa1(l,k)*dlog(sigmah(k)*pd+ptop)
          eps1 = varpa1(l,k)*sigmah(k)/(sigmah(k)*pd+ptop)
          do i = 1 , iym1
            do j = 1 , jendx
              eps = eps1*(sps2%ps(j,i)-pd)
              hstor(j,i,l) = hstor(j,i,l) + pdlog + &
                          tau(l,k)*atm2%t(i,k,j)/sps2%ps(j,i) + eps
            end do
          end do
        end do
      end do
    end if
!
    call time_end(subroutine_name,idindx)
!
  end subroutine spinit
!
! Compute deld, delh, integrate in time and add correction terms appropriately
!
  subroutine splitf(jstart,jend,istart,iend)
#ifndef IBM
    use mpi
#endif
    implicit none
#ifdef IBM
    include 'mpif.h'
#endif
!
    integer , intent(in) :: jstart , jend , istart , iend
    real(8) :: eps , eps1 , fac , gnuam , gnuan , gnuzm , pdlog , x , y
    integer :: i , j , k , l , n
    integer :: jm1 , jp1
    integer :: ierr , ii
    real(8) , dimension(iy*nsplit) :: wkrecv , wksend
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
!   satisfy p(0,j)=p(1,j); p(iy,j)=p(iym1,j); and similarly for the i's.
!
    if ( ieast /= mpi_proc_null ) then
      trans1 = sps1%ps(jxp,:)
    end if
    call mpi_sendrecv(trans1,iy,mpi_real8,ieast,1,      &
                      trans2,iy,mpi_real8,iwest,1,        &
                      mycomm,mpi_status_ignore,ierr)
    if ( iwest /= mpi_proc_null ) then
      sps1%ps(0,:) = trans2
    end if
    do i = 2 , iym1
      do j = jbegin , jendx
        jm1 = j-1
        psdot(j,i) = (sps1%ps(j,i)  +sps1%ps(j,i-1)+ &
                      sps1%ps(jm1,i)+sps1%ps(jm1,i-1))*d_rfour
      end do
    end do
!
#ifndef BAND
    do i = 2 , iym1
      if ( myid == 0 ) & 
        psdot(1,i) = (sps1%ps(1,i)+sps1%ps(1,i-1))*d_half
      if ( myid == nproc-1 ) &
        psdot(jendx,i) = (sps1%ps(jendx,i)+sps1%ps(jendx,i-1))*d_half
    end do
#endif
!
    do j = jbegin , jendx
      jm1 = j-1
      psdot(j,1) = (sps1%ps(j,1)+sps1%ps(jm1,1))*d_half
      psdot(j,iy) = (sps1%ps(j,iym1)+sps1%ps(jm1,iym1))*d_half
    end do
!
#ifndef BAND
    if ( myid == 0 ) then
      psdot(1,1) = sps1%ps(1,1)
      psdot(1,iy) = sps1%ps(1,iym1)
    end if
    if ( myid == nproc-1 ) then
      psdot(jendl,1) = sps1%ps(jendx,1)
      psdot(jendl,iy) = sps1%ps(jendx,iym1)
    end if
#endif
!
!=======================================================================
!
!   get deld(0), delh(0) from storage
!
   do n = 1 , nsplit
     do i = 1 , iy
       do j = 1 , jendl
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
      do i = 1 , iy
        do j = 1 , jendl
          uuu(j,i,k) = atm1%u(i,k,j)*mddom%msfd(j,i)
          vvv(j,i,k) = atm1%v(i,k,j)*mddom%msfd(j,i)
        end do
      end do
      if ( iwest /= mpi_proc_null ) then
        trans1 = uuu(1,:,k)
      end if
      call mpi_sendrecv(trans1,iy,mpi_real8,iwest,2, &
                        trans2,iy,mpi_real8,ieast,2, &
                        mycomm,mpi_status_ignore,ierr)
      if ( ieast /= mpi_proc_null ) then
        uuu(jxp+1,:,k) = trans2
      end if
      if ( iwest /= mpi_proc_null ) then
        trans1 = vvv(1,:,k)
      end if
      call mpi_sendrecv(trans1,iy,mpi_real8,iwest,2, &
                        trans2,iy,mpi_real8,ieast,2, &
                        mycomm,mpi_status_ignore,ierr)
      if ( ieast /= mpi_proc_null ) then
        vvv(jxp+1,:,k) = trans2
      end if
    end do
    do l = 1 , nsplit
      do i = 1 , iy
        do j = 1 , jendl
          deld(j,i,l,3) = d_zero
        end do
      end do

      do k = 1 , kz
        do i = 1 , iym1
          do j = 1 , jendx
            jp1 = j+1
            fac = dx2*mddom%msfx(j,i)*mddom%msfx(j,i)
            deld(j,i,l,3) = deld(j,i,l,3) + zmatxr(l,k) * &
               (-uuu(j,i+1,k)+uuu(jp1,i+1,k)-uuu(j,i,k)+uuu(jp1,i,k) + &
                 vvv(j,i+1,k)+vvv(jp1,i+1,k)-vvv(j,i,k)-vvv(jp1,i,k))/fac
          end do
        end do
      end do
    end do
!
!=======================================================================
! 
    do n = 1 , nsplit
      do i = 1 , iy
        do j = 1 , jendl
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
      do i = 1 , iy
        do j = 1 , jendl
          uuu(j,i,k) = atm2%u(i,k,j)*mddom%msfd(j,i)
          vvv(j,i,k) = atm2%v(i,k,j)*mddom%msfd(j,i)
        end do
      end do
      if ( iwest /= mpi_proc_null ) then
        trans1 = uuu(1,:,k)
      end if
      call mpi_sendrecv(trans1,iy,mpi_real8,iwest,2, &
                        trans2,iy,mpi_real8,ieast,2, &
                        mycomm,mpi_status_ignore,ierr)
      if ( ieast /= mpi_proc_null ) then
        uuu(jxp+1,:,k) = trans2
      end if
      if ( iwest /= mpi_proc_null ) then
        trans1 = vvv(1,:,k)
      end if
      call mpi_sendrecv(trans1,iy,mpi_real8,iwest,2, &
                        trans2,iy,mpi_real8,ieast,2, &
                        mycomm,mpi_status_ignore,ierr)
      if ( ieast /= mpi_proc_null ) then
        vvv(jxp+1,:,k) = trans2
      end if
    end do
    do l = 1 , nsplit
      do i = 1 , iy
        do j = 1 , jendl
          deld(j,i,l,2) = d_zero
        end do
      end do
      do k = 1 , kz
        do i = 1 , iym1
          do j = 1 , jendx
            jp1 = j+1
            fac = dx2*mddom%msfx(j,i)*mddom%msfx(j,i)
            deld(j,i,l,2) = deld(j,i,l,2) + zmatxr(l,k) * &
              (-uuu(j,i+1,k)+uuu(jp1,i+1,k)-uuu(j,i,k)+uuu(jp1,i,k) + &
                vvv(j,i+1,k)+vvv(jp1,i+1,k)-vvv(j,i,k)-vvv(jp1,i,k))/fac
          end do
        end do
      end do
    end do
!
!=======================================================================
!
    do n = 1 , nsplit
      do i = 1 , iy
        do j = 1 , jendl
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
      do i = 1 , iym1
        do j = 1 , jendx
          eps = eps1*(sps1%ps(j,i)-pd)
          delh(j,i,l,3) = pdlog + eps
        end do
      end do
      do k = 1 , kz
        pdlog = varpa1(l,k)*dlog(sigmah(k)*pd+ptop)
        eps1 = varpa1(l,k)*sigmah(k)/(sigmah(k)*pd+ptop)
        do i = 1 , iym1
          do j = 1 , jendx
            eps = eps1*(sps1%ps(j,i)-pd)
            delh(j,i,l,3) = delh(j,i,l,3) + pdlog +  &
                    tau(l,k)*atm1%t(i,k,j)/sps1%ps(j,i) + eps
          end do
        end do
      end do
    end do
!
!=======================================================================
! 
    do n = 1 , nsplit
      do i = 1 , iy
        do j = 1 , jendl
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
      do i = 1 , iym1
        do j = 1 , jendx
          eps = eps1*(sps2%ps(j,i)-pd)
          delh(j,i,l,2) = pdlog + eps
        end do
      end do
      do k = 1 , kz
        pdlog = varpa1(l,k)*dlog(sigmah(k)*pd+ptop)
        eps1 = varpa1(l,k)*sigmah(k)/(sigmah(k)*pd+ptop)
        do i = 1 , iym1
          do j = 1 , jendx
            eps = eps1*(sps2%ps(j,i)-pd)
            delh(j,i,l,2) = delh(j,i,l,2) + pdlog +  &
                     tau(l,k)*atm2%t(i,k,j)/sps2%ps(j,i) + eps
          end do
        end do
      end do
    end do
!
!=======================================================================
!
    do n = 1 , nsplit
      do i = 1 , iy
        do j = 1 , jendl
          delh(j,i,n,1) = delh(j,i,n,1) - delh(j,i,n,2)
        end do
      end do
    end do
!
!   put deld(0), delh(0) into storage
!
   do n = 1 , nsplit
     do i = 1 , iy
       do j = 1 , jendl
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
      do i = 2 , iym2
        do j = jbegin , jendm
          sps1%ps(j,i) = sps1%ps(j,i) - an(l)*ddsum(j,i,l)
          sps2%ps(j,i) = sps2%ps(j,i) - gnuan*ddsum(j,i,l)
        end do
      end do
    end do
    do l = 1 , nsplit
      do k = 1 , kz
        gnuam = gnuhf*am(k,l)
        do i = 2 , iym2
          do j = jbegin , jendm
            atm1%t(i,k,j) = atm1%t(i,k,j) + am(k,l)*ddsum(j,i,l)
            atm2%t(i,k,j) = atm2%t(i,k,j) + gnuam*ddsum(j,i,l)
          end do
        end do
      end do
    end do
!=======================================================================
    ii = 0
    do l = 1 , nsplit
      do i = 1 , iy
        ii = ii + 1
        wksend(ii) = dhsum(jxp,i,l)
      end do
    end do
    call mpi_sendrecv(wksend,iy*nsplit,mpi_real8,ieast,1, &
                      wkrecv,iy*nsplit,mpi_real8,iwest,1, &
                      mycomm,mpi_status_ignore,ierr)
    ii = 0
    do l = 1 , nsplit
      do i = 1 , iy
        ii = ii + 1
        dhsum(0,i,l) = wkrecv(ii)
      end do
    end do
    do l = 1 , nsplit
      do k = 1 , kz
        gnuzm = gnuhf*zmatx(k,l)
        do i = 2 , iym1
          do j = jbegin , jendx
            jm1 = j-1
            fac = psdot(j,i)/(dx2*mddom%msfd(j,i))
            x = fac*(dhsum(j,i,l)+dhsum(j,i-1,l) - &
                     dhsum(jm1,i,l)-dhsum(jm1,i-1,l))
            y = fac*(dhsum(j,i,l)-dhsum(j,i-1,l) + &
                     dhsum(jm1,i,l)-dhsum(jm1,i-1,l))
            atm1%u(i,k,j) = atm1%u(i,k,j) - zmatx(k,l)*x
            atm1%v(i,k,j) = atm1%v(i,k,j) - zmatx(k,l)*y
            atm2%u(i,k,j) = atm2%u(i,k,j) - gnuzm*x
            atm2%v(i,k,j) = atm2%v(i,k,j) - gnuzm*y
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
!
#ifndef IBM
    use mpi
#else
    include 'mpif.h'
#endif
    implicit none
!
    real(8) :: dtau2 , fac
    integer :: i , j , m2 , n , n0 , n1 , n2 , ns , nw
    integer :: jm1, jp1
    integer :: ierr
    real(8) , dimension(iy*2) :: wkrecv , wksend
    character (len=64) :: subroutine_name='spstep'
    integer :: idindx=0
!
    call time_begin(subroutine_name,idindx)
!--
    
    do n = 1 , nsplit
      do i = 1 , iy
        do j = 1 , jendl
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
      do i = 1 , iym1
        do j = 1 , jendx
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
      if ( ieast /= mpi_proc_null ) then
        trans1 = delh(jxp,:,ns,n0)
      end if
      call mpi_sendrecv(trans1,iy,mpi_real8,ieast,1, &
                        trans2,iy,mpi_real8,iwest,1, &
                        mycomm,mpi_status_ignore,ierr)
      if ( iwest /= mpi_proc_null ) then
        delh(0,:,ns,n0) = trans2
      end if
      do i = 2 , iym1
        do j = jbegin , jendx
          jm1 = j-1
          fac = dx2*mddom%msfx(j,i)
          work(j,i,1) = (delh(j,i,ns,n0)+delh(j,i-1,ns,n0)  &
                      & -delh(jm1,i,ns,n0)-delh(jm1,i-1,ns,n0))/fac
          work(j,i,2) = (delh(j,i,ns,n0)+delh(jm1,i,ns,n0) &
                      & -delh(j,i-1,ns,n0)-delh(jm1,i-1,ns,n0))/fac
        end do
      end do
!
!=======================================================================
!
      do nw = 1 , 2
        do i = 2 , iym1
          do j = jbegin , jendx
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
      do i = 2 , iym1
        do j = jbegin , jendx
          uu(j,i) = work(j,i,1)*mddom%msfd(j,i)
          vv(j,i) = work(j,i,2)*mddom%msfd(j,i)
        end do
      end do
!
      do i = 1 , iy
        wksend(i) = uu(1,i)
        wksend(i+iy) = vv(1,i)
      end do
      call mpi_sendrecv(wksend,2*iy,mpi_real8,iwest,2,             &
                      & wkrecv,2*iy,mpi_real8,ieast,2,             &
                      & mycomm,mpi_status_ignore,ierr)
      do i = 1 , iy
        uu(jxp+1,i) = wkrecv(i)
        vv(jxp+1,i) = wkrecv(i+iy)
      end do
!
      do i = 2 , iym2
        do j = jbegin , jendm
          jp1 = j+1
          fac = dx2*mddom%msfx(j,i)*mddom%msfx(j,i)
          work(j,i,3) = (-uu(j,i+1)+uu(jp1,i+1)-uu(j,i)+uu(jp1,i) &
                         +vv(j,i+1)+vv(jp1,i+1)-vv(j,i)-vv(jp1,i))/fac
        end do
      end do
!
!=======================================================================
!
      do i = 2 , iym2
        do j = jbegin , jendm
!           work3: 2,iym2 on cross grid
          deld(j,i,ns,n1) = deld(j,i,ns,n0) - dtau(ns)*work(j,i,3) + &
                            deld(j,i,ns,3)/m2
          delh(j,i,ns,n1) = delh(j,i,ns,n0) - dtau(ns)*hbar(ns) * &
                            deld(j,i,ns,n0)/sps1%ps(j,i)+delh(j,i,ns,3)/m2
        end do
      end do
 
!     not in Madala (1987)
      fac = (aam(ns)-d_one)/aam(ns)
#ifndef BAND
      do i = 2 , iym2
        if ( myid == 0 ) &
          delh(1,i,ns,n1) = delh(1,i,ns,n0)*fac
        if ( myid == nproc-1 ) &
          delh(jendx,i,ns,n1) = delh(jendx,i,ns,n0)*fac
      end do
#endif
      do j = 1 , jendx
        delh(j,1,ns,n1) = delh(j,1,ns,n0)*fac
        delh(j,iym1,ns,n1) = delh(j,iym1,ns,n0)*fac
      end do
!
      do i = 1 , iym1
        do j = 1 , jendx
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
        if ( ieast /= mpi_proc_null ) then
          trans1 = delh(jxp,:,ns,n1)
        end if
        call mpi_sendrecv(trans1,iy,mpi_real8,ieast,1, &
                          trans2,iy,mpi_real8,iwest,1, &
                          mycomm,mpi_status_ignore,ierr)
        if ( iwest /= mpi_proc_null ) then
          delh(0,:,ns,n1) = trans2
        end if
        do i = 2 , iym1
          do j = jbegin , jendx
            jm1 = j-1
            fac = dx2*mddom%msfx(j,i)
            work(j,i,1) = (delh(j,i,ns,n1)+delh(j,i-1,ns,n1)- &
                           delh(jm1,i,ns,n1)-delh(jm1,i-1,ns,n1))/fac
            work(j,i,2) = (delh(j,i,ns,n1)+delh(jm1,i,ns,n1)- &
                           delh(j,i-1,ns,n1)-delh(jm1,i-1,ns,n1))/fac
          end do
        end do
!
!=======================================================================
!
        do nw = 1 , 2
          do i = 2 , iym1
            do j = jbegin , jendx
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
        do i = 2 , iym1
          do j = jbegin , jendx
            uu(j,i) = work(j,i,1)*mddom%msfd(j,i)
            vv(j,i) = work(j,i,2)*mddom%msfd(j,i)
          end do
        end do
!
        do i = 1 , iy
          wksend(i) = uu(1,i)
          wksend(i+iy) = vv(1,i)
        end do
        call mpi_sendrecv(wksend,2*iy,mpi_real8,iwest,2,           &
                          wkrecv,2*iy,mpi_real8,ieast,2,           &
                          mycomm,mpi_status_ignore,ierr)
        do i = 1 , iy
          uu(jxp+1,i) = wkrecv(i)
          vv(jxp+1,i) = wkrecv(i+iy)
        end do
!
        do i = 2 , iym2
          do j = jbegin , jendm
            jp1 = j+1
            fac = dx2*mddom%msfx(j,i)*mddom%msfx(j,i)
            work(j,i,3) = (-uu(j,i+1)+uu(jp1,i+1)-uu(j,i)+uu(jp1,i) + &
                            vv(j,i+1)+vv(jp1,i+1)-vv(j,i)-vv(jp1,i))/fac
          end do
        end do
!
!=======================================================================
!
        do i = 2 , iym2
          do j = jbegin , jendm
            deld(j,i,ns,n2) = deld(j,i,ns,n0) - dtau2*work(j,i,3) + &
                              deld(j,i,ns,3)/aam(ns)
            delh(j,i,ns,n2) = delh(j,i,ns,n0) - dtau2*hbar(ns) * &
                              deld(j,i,ns,n1)/sps1%ps(j,i) +     &
                              delh(j,i,ns,3)/aam(ns)
          end do
        end do
!
!       not in Madala (1987)
!
#ifndef BAND
        do i = 2 , iym2
          if ( myid == 0 ) &
            delh(1,i,ns,n2) = d_two*delh(1,i,ns,n1)-delh(1,i,ns,n0)
          if ( myid == nproc-1 ) &
            delh(jendx,i,ns,n2) = d_two*delh(jendx,i,ns,n1)-delh(jendx,i,ns,n0)
        end do
#endif
        do j = 1 , jendx
          delh(j,1,ns,n2) = d_two*delh(j,1,ns,n1)-delh(j,1,ns,n0)
          delh(j,iym1,ns,n2) = d_two*delh(j,iym1,ns,n1)-delh(j,iym1,ns,n0)
        end do
!
        do i = 1 , iym1
          do j = 1 , jendx
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
