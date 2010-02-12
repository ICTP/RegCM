!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of RegCM model.
!
!    RegCM model is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    RegCM model is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with RegCM model.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 
      subroutine spinit(ptop,sigma,kv1)
!
!** compute vertical modes.
!
      use mod_regcm_param
      use mod_param1
      use mod_param2
      use mod_iunits
      use mod_bdycod
      use mod_main
      use mod_split
      use mod_bxq
      use mod_tmpsav
#ifdef MPP1
      use mod_mppio
      use mpi
#endif
      implicit none
!
! Dummy arguments
!
      integer :: kv1
      real(8) :: ptop
      real(8) , dimension(kv1) :: sigma
      intent (in) ptop
!
! Local variables
!
      real(8) :: eps , eps1 , fac , pdlog
      integer :: i , ijlx , j , k , l , n , ns
      logical :: lstand
#ifdef MPP1
      integer :: ierr
#endif
!
!     lstand = .true. if standard atmosphere t to be used (ignore input
!     tbarh and ps in that case).  otherwise, ps and tbarh must
!     be defined on input.  note that in either case, pt must
!     also be defined on input (common block named cvert).
!
!
!     ******dtau = time steps(in sec)for modes in split explicit is
!     ******specified in namelist as array dtsplit
!
!**   zero new arrays

      dstor = 0.0
      hstor = 0.0
!
!**   compute m.
      do ns = 1 , nsplit
        m(ns) = nint(dt/dtau(ns))
        if ( jyear.ne.jyear0 .or. ktau.ne.0 ) m(ns)                     &
           & = nint(.5*dt/dtau(ns))
      end do
#ifdef MPP1
      if ( myid.eq.0 ) print * , 'dt, dtau = ' , dt , dtau
#else
      print * , 'dt, dtau = ' , dt , dtau
#endif
!
!**   compute pt, ps and tbarh for use in vmodes.
      ps = 0.
      do k = 1 , kx
        tbarh(k) = 0.
      end do
      pt = ptop
#ifdef MPP1
      ijlx = ixm1*jendx
      do j = 1 , jendx
        do i = 1 , ixm1
          ps = ps + psa(i,j)/ijlx
        end do
      end do
#else
      ijlx = ixm1*jxm1
      do j = 1 , jxm1
        do i = 1 , ixm1
          ps = ps + psa(i,j)/ijlx
        end do
      end do
#endif
      do k = 1 , kx
#ifdef MPP1
        do j = 1 , jendx
          do i = 1 , ixm1
            tbarh(k) = tbarh(k) + ta(i,k,j)/(psa(i,j)*ijlx)
          end do
        end do
#else
        do j = 1 , jxm1
          do i = 1 , ixm1
            tbarh(k) = tbarh(k) + ta(i,k,j)/(psa(i,j)*ijlx)
          end do
        end do
#endif
      end do
!
!**   compute vertical modes.
      lstand = .true.
      if ( jyear.ne.jyear0 .or. ktau.ne.0 ) lstand = .true.
      call vmodes(lstand,sigma,kv1)
!
!**   subract a4 from a for use in computing am.
      do l = 1 , kx
        do k = 1 , kx
          a(k,l) = a(k,l) - a4(k,l)
        end do
      end do
!
!**   compute am and an.
      do n = 1 , nsplit
        an(n) = 0.
        do l = 1 , kx
          an(n) = an(n) + dsigma(l)*zmatx(l,n)
        end do
        do k = 1 , kx
          am(k,n) = 0.
          tau(n,k) = 0.
        end do
        do l = 1 , kx
          do k = 1 , kx
            am(k,n) = am(k,n) + a(k,l)*zmatx(l,n)
            tau(n,k) = tau(n,k) + r*zmatxr(n,l)*hydros(l,k)
          end do
        end do
!
        do k = 1 , kxp1
          varpa1(n,k) = 0.
        end do
        do l = 1 , kx
          do k = 1 , kxp1
            varpa1(n,k) = varpa1(n,k) + r*zmatxr(n,l)*hydroc(l,k)
          end do
        end do
      end do
!
!**   multiply am, an and zmatx by factor.
      do l = 1 , nsplit
        fac = 2.*dt/(2.*dble(m(l))+1.)
        if ( jyear.ne.jyear0 .or. ktau.ne.0 )                           &
           & fac = dt/(2.*dble(m(l))+1.)
#ifdef MPP1
        if ( myid.eq.0 ) print * , 'm, fac = ' , m(l) , fac
#else
        print * , 'm, fac = ' , m(l) , fac
#endif
        an(l) = an(l)*fac
        do k = 1 , kx
          zmatx(k,l) = zmatx(k,l)*fac
          am(k,l) = am(k,l)*fac
        end do
      end do
!
      if ( ifrest ) then
#ifdef MPP1
        if ( myid.eq.0 ) then
          read (iutrs) dstor_io
          read (iutrs) hstor_io
          read (iutrs) uj1 , uj2 , ujlx , ujl
          read (iutrs) ui1_io , ui2_io , uilx_io , uil_io
          read (iutrs) vj1 , vj2 , vjlx , vjl
          read (iutrs) vi1_io , vi2_io , vilx_io , vil_io
          do j = 1 , jx
            do n = 1 , nsplit
              do i = 1 , ix
                sav_0d(i,n,j) = dstor_io(i,j,n)
                sav_0d(i,n+nsplit,j) = hstor_io(i,j,n)
              end do
            end do
          end do
          do j = 1 , jx
            do k = 1 , kx
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
        call mpi_scatter(sav_0d(1,1,1),ix*nsplit*2*jxp,                 &
                       & mpi_double_precision,sav0d(1,1,1),             &
                       & ix*nsplit*2*jxp,mpi_double_precision,0,        &
                       & mpi_comm_world,ierr)
        do j = 1 , jendl
          do n = 1 , nsplit
            do i = 1 , ix
              dstor(i,j,n) = sav0d(i,n,j)
              hstor(i,j,n) = sav0d(i,n+nsplit,j)
            end do
          end do
        end do
        call mpi_scatter(sav_6(1,1,1),kx*8*jxp,mpi_double_precision,    &
                       & sav6(1,1,1),kx*8*jxp,mpi_double_precision,0,   &
                       & mpi_comm_world,ierr)
        do j = 1 , jendl
          do k = 1 , kx
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
        call mpi_bcast(uj1,ix*kx,mpi_double_precision,0,mpi_comm_world, &
                     & ierr)
        call mpi_bcast(uj2,ix*kx,mpi_double_precision,0,mpi_comm_world, &
                     & ierr)
        call mpi_bcast(vj1,ix*kx,mpi_double_precision,0,mpi_comm_world, &
                     & ierr)
        call mpi_bcast(vj2,ix*kx,mpi_double_precision,0,mpi_comm_world, &
                     & ierr)
        call mpi_bcast(ujlx,ix*kx,mpi_double_precision,0,mpi_comm_world,&
                     & ierr)
        call mpi_bcast(ujl,ix*kx,mpi_double_precision,0,mpi_comm_world, &
                     & ierr)
        call mpi_bcast(vjlx,ix*kx,mpi_double_precision,0,mpi_comm_world,&
                     & ierr)
        call mpi_bcast(vjl,ix*kx,mpi_double_precision,0,mpi_comm_world, &
                     & ierr)
        if ( myid.ne.nproc-1 ) then
          do k = 1 , kx
            var1snd(k,1) = ui1(k,jxp)
            var1snd(k,2) = vi1(k,jxp)
            var1snd(k,3) = ui2(k,jxp)
            var1snd(k,4) = vi2(k,jxp)
            var1snd(k,5) = uilx(k,jxp)
            var1snd(k,6) = vilx(k,jxp)
            var1snd(k,7) = uil(k,jxp)
            var1snd(k,8) = vil(k,jxp)
          end do
        end if
        call mpi_sendrecv(var1snd(1,1),kx*8,mpi_double_precision,ieast, &
                        & 1,var1rcv(1,1),kx*8,mpi_double_precision,     &
                        & iwest,1,mpi_comm_world,mpi_status_ignore,ierr)
        if ( myid.ne.0 ) then
          do k = 1 , kx
            ui1(k,0) = var1rcv(k,1)
            vi1(k,0) = var1rcv(k,2)
            ui2(k,0) = var1rcv(k,3)
            vi2(k,0) = var1rcv(k,4)
            uilx(k,0) = var1rcv(k,5)
            vilx(k,0) = var1rcv(k,6)
            uil(k,0) = var1rcv(k,7)
            vil(k,0) = var1rcv(k,8)
          end do
        end if
        if ( myid.ne.0 ) then
          do k = 1 , kx
            var1snd(k,1) = ui1(k,1)
            var1snd(k,2) = vi1(k,1)
            var1snd(k,3) = ui2(k,1)
            var1snd(k,4) = vi2(k,1)
            var1snd(k,5) = uilx(k,1)
            var1snd(k,6) = vilx(k,1)
            var1snd(k,7) = uil(k,1)
            var1snd(k,8) = vil(k,1)
          end do
        end if
        call mpi_sendrecv(var1snd(1,1),kx*8,mpi_double_precision,iwest, &
                        & 2,var1rcv(1,1),kx*8,mpi_double_precision,     &
                        & ieast,2,mpi_comm_world,mpi_status_ignore,ierr)
        if ( myid.ne.nproc-1 ) then
          do k = 1 , kx
            ui1(k,jxp+1) = var1rcv(k,1)
            vi1(k,jxp+1) = var1rcv(k,2)
            ui2(k,jxp+1) = var1rcv(k,3)
            vi2(k,jxp+1) = var1rcv(k,4)
            uilx(k,jxp+1) = var1rcv(k,5)
            vilx(k,jxp+1) = var1rcv(k,6)
            uil(k,jxp+1) = var1rcv(k,7)
            vil(k,jxp+1) = var1rcv(k,8)
          end do
        end if
#else
        read (iutrs) dstor
        read (iutrs) hstor
        read (iutrs) uj1 , uj2 , ujlx , ujl
        read (iutrs) ui1 , ui2 , uilx , uil
        read (iutrs) vj1 , vj2 , vjlx , vjl
        read (iutrs) vi1 , vi2 , vilx , vil
#endif
      else
!
!=======================================================================
!******* divergence manipulations (0)
!
!**     compute divergence z from u and v
!       ( u must be pstar * u ; similarly for v )
!       ( note: map scale factors have been inverted in model (init) )
!
        do k = 1 , kx
#ifdef MPP1
          do j = 1 , jendl
            do i = 1 , ix
              uuu(i,k,j) = ub(i,k,j)*msfd(i,j)
              vvv(i,k,j) = vb(i,k,j)*msfd(i,j)
            end do
          end do
#else
          do j = 1 , jx
            do i = 1 , ix
              uuu(i,k,j) = ub(i,k,j)*msfd(i,j)
              vvv(i,k,j) = vb(i,k,j)*msfd(i,j)
            end do
          end do
#endif
        end do
#ifdef MPP1
        call mpi_sendrecv(uuu(1,1,1),ix*kx,mpi_double_precision,iwest,2,&
                        & uuu(1,1,jxp+1),ix*kx,mpi_double_precision,    &
                        & ieast,2,mpi_comm_world,mpi_status_ignore,ierr)
        call mpi_sendrecv(vvv(1,1,1),ix*kx,mpi_double_precision,iwest,2,&
                        & vvv(1,1,jxp+1),ix*kx,mpi_double_precision,    &
                        & ieast,2,mpi_comm_world,mpi_status_ignore,ierr)
#endif
!
        do l = 1 , nsplit
#ifdef MPP1
          do j = 1 , jendl
            do i = 1 , ix
              dstor(i,j,l) = 0.
            end do
          end do
#else
          do j = 1 , jx
            do i = 1 , ix
              dstor(i,j,l) = 0.
            end do
          end do
#endif
        end do
        do l = 1 , nsplit
          do k = 1 , kx
#ifdef MPP1
            do j = 1 , jendx
              do i = 1 , ixm1
                fac = dx2*msfx(i,j)*msfx(i,j)
                dstor(i,j,l) = dstor(i,j,l) + zmatxr(l,k)               &
                             & *(-uuu(i+1,k,j)+uuu(i+1,k,j+1)-uuu(i,k,j)&
                             & +uuu(i,k,j+1)+vvv(i+1,k,j)+vvv(i+1,k,j+1)&
                             & -vvv(i,k,j)-vvv(i,k,j+1))/fac
              end do
            end do
#else
            do j = 1 , jxm1
              do i = 1 , ixm1
                fac = dx2*msfx(i,j)*msfx(i,j)
                dstor(i,j,l) = dstor(i,j,l) + zmatxr(l,k)               &
                             & *(-uuu(i+1,k,j)+uuu(i+1,k,j+1)-uuu(i,k,j)&
                             & +uuu(i,k,j+1)+vvv(i+1,k,j)+vvv(i+1,k,j+1)&
                             & -vvv(i,k,j)-vvv(i,k,j+1))/fac
              end do
            end do
#endif
          end do
        end do
!
!=======================================================================
!
!******* geopotential manipulations
        do l = 1 , nsplit
          pdlog = varpa1(l,kxp1)*dlog(sigmah(kxp1)*pd+pt)
          eps1 = varpa1(l,kxp1)*sigmah(kxp1)/(sigmah(kxp1)*pd+pt)
#ifdef MPP1
          do j = 1 , jendx
            do i = 1 , ixm1
              eps = eps1*(psb(i,j)-pd)
              hstor(i,j,l) = pdlog + eps
            end do
          end do
#else
          do j = 1 , jxm1
            do i = 1 , ixm1
              eps = eps1*(psb(i,j)-pd)
              hstor(i,j,l) = pdlog + eps
            end do
          end do
#endif
          do k = 1 , kx
            pdlog = varpa1(l,k)*dlog(sigmah(k)*pd+pt)
            eps1 = varpa1(l,k)*sigmah(k)/(sigmah(k)*pd+pt)
#ifdef MPP1
            do j = 1 , jendx
              do i = 1 , ixm1
                eps = eps1*(psb(i,j)-pd)
                hstor(i,j,l) = hstor(i,j,l) + pdlog + tau(l,k)*tb(i,k,j)&
                             & /psb(i,j) + eps
              end do
            end do
#else
            do j = 1 , jxm1
              do i = 1 , ixm1
                eps = eps1*(psb(i,j)-pd)
                hstor(i,j,l) = hstor(i,j,l) + pdlog + tau(l,k)*tb(i,k,j)&
                             & /psb(i,j) + eps
              end do
            end do
#endif
          end do
        end do
      end if
!
      end subroutine spinit
