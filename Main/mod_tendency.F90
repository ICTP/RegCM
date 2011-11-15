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
 
module mod_tendency

  use mod_runparams
  use mod_mpmessage
  use mod_service
  use mod_memutil
  use mod_atm_interface
  use mod_che_interface
  use mod_cu_interface
  use mod_lm_interface
  use mod_rad_interface
  use mod_pbl_interface
  use mod_bdycod
  use mod_precip
  use mod_sun
  use mod_slice
  use mod_diagnosis
  use mod_advection
  use mod_diffusion
#ifdef CLM
  use mod_clm
  use mod_mtrxclm
  use clm_varsur
#endif

  private

  public :: allocate_mod_tend , tend

  real(8) , pointer , dimension(:,:) :: divl
  real(8) , pointer , dimension(:,:,:) :: ttld , xkc , xkcf , td
  real(8) , pointer , dimension(:,:) :: bdyewrcv , bdyewsnd
  real(8) , pointer , dimension(:,:) :: bdynsrcv , bdynssnd
  real(8) , pointer , dimension(:,:,:) :: ps4
  real(8) , pointer , dimension(:,:,:) :: ps_4 
  real(8) , pointer , dimension(:,:) :: var2rcv , var2snd
  real(8) , pointer , dimension(:,:) :: tvar1rcv , tvar1snd
  real(8) , pointer , dimension(:,:) :: qvcs
  real(8) , pointer , dimension(:,:,:) :: wrkkuo1
  real(8) , pointer , dimension(:,:,:) :: wrkkuo2
  real(8) , pointer , dimension(:) :: pstrans1
  real(8) , pointer , dimension(:) :: pstrans2

  integer(8) , parameter :: irep = 50

  contains

  subroutine allocate_mod_tend(lband)
    implicit none
    logical , intent(in) :: lband
    integer :: n1 , n2

    call getmem2d(divl,1,iy,1,kz,'tendency:divl')
    call getmem3d(ttld,1,iy,1,kz,1,jxp,'tendency:ttld')
    call getmem3d(xkc,1,iy,1,kz,1,jxp,'tendency:xkc')
    call getmem3d(xkcf,1,iy,1,kzp1,1,jxp,'tendency:xkcf')
    call getmem3d(td,1,iy,1,kz,1,jxp,'tendency:td')
    if ( .not. lband ) then
      call getmem2d(bdyewrcv,1,iy,1,kz*16+4,'tendency:bdyewrcv')
      call getmem2d(bdyewsnd,1,iy,1,kz*16+4,'tendency:bdyewsnd')
    end if
    call getmem2d(bdynsrcv,1,nspgx,1,kz*16+4,'tendency:bdynsrcv')
    call getmem2d(bdynssnd,1,nspgx,1,kz*16+4,'tendency:bdynssnd')
    call getmem3d(ps4,1,iy,1,4,1,jxp,'tendency:ps4')
    call getmem3d(ps_4,1,iy,1,4,1,jx,'tendency:ps_4')
    n1 = kz*11 + 1
    n2 = kz*10
    if ( ichem == 1 ) then
      n1 = n1 + ntr*kz*2
      n2 = n2 + ntr*kz*2
    end if
    if ( ibltyp == 2 .or. ibltyp == 99 ) then
      n1 = n1 + kz
      n2 = n2 + 2*kz
    end if
    call getmem2d(tvar1rcv,1,iy,1,n1,'tendency:tvar1rcv')
    call getmem2d(tvar1snd,1,iy,1,n1,'tendency:tvar1snd')
    call getmem2d(var2rcv,1,iy,1,n2,'tendency:var2rcv')
    call getmem2d(var2snd,1,iy,1,n2,'tendency:var2snd')
    call getmem2d(qvcs,1,iy,1,kz,'tendency:qvcs')
    call getmem1d(pstrans1,1,iy,'tendency:pstrans1')
    call getmem1d(pstrans2,1,iy,'tendency:pstrans2')
    if ( icup == 1 ) then
      call getmem3d(wrkkuo1,1,iy,0,jxp+1,1,kz,'tendency:wrkkuo1')
      call getmem3d(wrkkuo2,1,iy,0,jxp+1,1,kz,'tendency:wrkkuo2')
    end if
  end subroutine allocate_mod_tend

  subroutine tend(iexec)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!     this subroutine computes the tendencies of the prognostic       c
!     variables p*, u, v, and t.                                      c
!                                                                     c
!     p*u, p*v, p*t ,p*qv, and p*qc stored in main common block.      c
!                                                                     c
!     all the two-dimension arrays stored in main common block.       c
!                                                                     c
!     east/west boundary conditions stored in common block /bdycod/ . c
!                                                                     c
!     north/south boundary conditions stored in common block          c
!              /bdycod/.                                              c
!                                                                     c
!     all the integers stored in common block /param1/.               c
!                                                                     c
!     all the constants stored in common block /param1/.              c
!                                                                     c
!     iexec  : = 1 ; represents this subroutine is called for the     c
!                    first time in this forecast run.                 c
!              > 1 ; represents subsequent calls to this subroutine.  c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
#ifndef IBM
    use mpi
#else
    include 'mpif.h'
#endif
    implicit none
!
    integer :: iexec
    intent (inout) iexec
!
    real(8) :: cell , chias , chibs , dudx , dudy , dvdx , dvdy ,  &
               psabar , psasum , pt2bar , pt2tot , ptnbar ,        &
               ptntot , qcas , qcbs , qvas , qvbs , rovcpm ,       &
               rtbar , sigpsa , tv , tv1 , tv2 , tv3 , tv4 , tva , &
               tvavg , tvb , tvc , xmsf , xtm1 , theta , eccf
    real(8) , pointer , dimension(:,:,:) :: spchiten , spchi , spchia , &
                                            spchib3d
    integer :: i , iptn , itr , j , k , lev , n
    integer :: jm1, jp1
    integer :: ierr , icons_mpi , numrec
    logical :: loutrad
    character (len=32) :: appdat
    character (len=64) :: subroutine_name='tend'
    integer :: idindx=0
!
    call time_begin(subroutine_name,idindx)
!
!----------------------------------------------------------------------
!   fill up the boundary slices:
!
    if ( .not. ifrest .and. iexec == 1 ) then
      call bdyval(xbctime,iexec)
      iexec = 2
    else
      iexec = 2
    end if
!
!   Calculate eccentricity factor for radiation calculations
!
#ifdef CLM
    eccf  = r2ceccf
#else
    calday = yeardayfrac(idatex)
    theta = twopi*calday/dayspy
    eccf = 1.000110D0 + 0.034221D0*dcos(theta) +  &
           0.001280D0 * dsin(theta) + &
           0.000719D0 * dcos(d_two*theta) + &
           0.000077D0 * dsin(d_two*theta)
#endif
!
!----------------------------------------------------------------------
!   multiply ua and va by inverse of mapscale factor at dot point:
!
    do j = 1 , jendl
      do k = 1 , kz
        do i = 1 , iy
          atm1%u(i,k,j) = atm1%u(i,k,j)*mddom%msfd(i,j)
          atm1%v(i,k,j) = atm1%v(i,k,j)*mddom%msfd(i,j)
        end do
      end do
    end do
    if ( ieast /= mpi_proc_null ) then
      pstrans1 = sps1%ps(jxp,:)
    end if
    call mpi_sendrecv(pstrans1,iy,mpi_real8,ieast,1, &
                      pstrans2,iy,mpi_real8,iwest,1, &
                      mycomm,mpi_status_ignore,ierr)
    if ( iwest /= mpi_proc_null ) then
      sps1%ps(0,:) = pstrans2
    end if
    if ( iwest /= mpi_proc_null ) then
      pstrans1 = sps1%ps(1,:)
    end if
    call mpi_sendrecv(pstrans1,iy,mpi_real8,iwest,2, &
                      pstrans2,iy,mpi_real8,ieast,2, &
                      mycomm,mpi_status_ignore,ierr)
    if ( ieast /= mpi_proc_null ) then
      sps1%ps(jxp+1,:) = pstrans2
    end if
!
!   decouple u, v, t, qv, and qc
!
#ifndef BAND
    do j = 1 , jendl
      if ( myid == 0 .and. j == 1 ) then
!       lateral slices:
!       west boundary:
        do k = 1 , kz
          do i = 1 , iy
            atmx%u(i,k,j) = uj1(i,k)
            atmx%v(i,k,j) = vj1(i,k)
          end do
        end do
        if ( iboudy == 3 .or. iboudy == 4 ) then
!         inflow/outflow dependence:
          do k = 1 , kz
            do i = 1 , iy
              if ( atmx%u(i,k,j) < d_zero ) then
                atmx%v(i,k,j) = vj2(i,k)
                atmx%u(i,k,j) = uj2(i,k)
              end if
            end do
          end do
        end if
      else if ( myid == 0 .and. j == 2 ) then
        do k = 1 , kz
          do i = 1 , iy
            atmx%u(i,k,j) = uj2(i,k)
            atmx%v(i,k,j) = vj2(i,k)
          end do
        end do
        if ( iboudy == 3 .or. iboudy == 4 ) then
!         inflow/outflow dependence:
          do k = 1 , kz
!           south boundary:
            if ( atmx%v(1,k,j) < d_zero ) then
              atmx%v(1,k,j) = atmx%v(2,k,j)
              atmx%u(1,k,j) = atmx%u(2,k,j)
            end if
!           north boundary:
            if ( atmx%v(iy,k,j) >= d_zero ) then
              atmx%v(iy,k,j) = atmx%v(iym1,k,j)
              atmx%u(iy,k,j) = atmx%u(iym1,k,j)
            end if
          end do
        end if
      else if ( myid == nproc-1 .and. j == jendl-1 ) then
        do k = 1 , kz
          do i = 1 , iy
            atmx%u(i,k,j) = ujlx(i,k)
            atmx%v(i,k,j) = vjlx(i,k)
          end do
        end do
        if ( iboudy == 3 .or. iboudy == 4 ) then
!         inflow/outflow dependence:
          do k = 1 , kz
!           south boundary:
            if ( atmx%v(1,k,j) < d_zero ) then
              atmx%v(1,k,j) = atmx%v(2,k,j)
              atmx%u(1,k,j) = atmx%u(2,k,j)
            end if
!           north boundary:
            if ( atmx%v(iy,k,j) >= d_zero ) then
              atmx%v(iy,k,j) = atmx%v(iym1,k,j)
              atmx%u(iy,k,j) = atmx%u(iym1,k,j)
            end if
          end do
        end if
      else if ( myid == nproc-1 .and. j == jendl ) then
!       east boundary:
!       no inflow/outflow dependence:
        do k = 1 , kz
          do i = 1 , iy
            atmx%u(i,k,j) = ujl(i,k)
            atmx%v(i,k,j) = vjl(i,k)
          end do
        end do
        if ( iboudy == 3 .or. iboudy == 4 ) then
!         inflow/outflow dependence:
          do k = 1 , kz
            do i = 1 , iy
              if ( atmx%u(i,k,j) >= d_zero ) then
                atmx%v(i,k,j) = vjlx(i,k)
                atmx%u(i,k,j) = ujlx(i,k)
              end if
            end do
          end do
        end if
      end if
    end do
#endif

    do j = 1 , jendl
      jm1 = j-1
!
!     interior slice:
!     interior points:
!
#ifndef BAND
      if ((.not.(myid == 0 .and. (j == 1 .or. j == 2))) .and. &
         (.not.(myid == nproc-1 .and. &
         (j == jendl-1 .or. j == jendl))) ) then
#endif
        do k = 1 , kz
          do i = 3 , iym2
            psabar=(sps1%ps(j,i)+sps1%ps(jm1,i)+ &
                    sps1%ps(j,i-1)+sps1%ps(jm1,i-1))*d_rfour
            xmsf = mddom%msfd(i,j)
            atmx%u(i,k,j) = atm1%u(i,k,j)/(psabar*xmsf)
            atmx%v(i,k,j) = atm1%v(i,k,j)/(psabar*xmsf)
          end do
        end do
!
!       north/south boundary points:
!       no inflow/outflow dependence:
!
        do k = 1 , kz
          atmx%u(2,k,j) = ui2(k,j)
          atmx%u(iym1,k,j) = uilx(k,j)
          atmx%v(2,k,j) = vi2(k,j)
          atmx%v(iym1,k,j) = vilx(k,j)
          atmx%u(1,k,j) = ui1(k,j)
          atmx%u(iy,k,j) = uil(k,j)
          atmx%v(1,k,j) = vi1(k,j)
          atmx%v(iy,k,j) = vil(k,j)
        end do
        if ( iboudy == 3 .or. iboudy == 4 ) then
!         inflow/outflow dependence:
          do k = 1 , kz
!           south boundary:
            if ( atmx%v(1,k,j) < d_zero ) then
              atmx%v(1,k,j) = atmx%v(2,k,j)
              atmx%u(1,k,j) = atmx%u(2,k,j)
            end if
!           north boundary:
            if ( atmx%v(iy,k,j) >= d_zero ) then
              atmx%v(iy,k,j) = atmx%v(iym1,k,j)
              atmx%u(iy,k,j) = atmx%u(iym1,k,j)
            end if
          end do
        end if
#ifndef BAND
      end if
#endif
    end do
!
    do j = 1 , jendx
      do k = 1 , kz
        do i = 1 , iym1
          atmx%t(i,k,j) = atm1%t(i,k,j)/sps1%ps(j,i)
          atmx%qv(i,k,j) = atm1%qv(i,k,j)/sps1%ps(j,i)
          atmx%qc(i,k,j) = atm1%qc(i,k,j)/sps1%ps(j,i)
        end do
      end do
    end do
!
!------------------------------------------------------------------------
!
    if ( ichem == 1 ) then
!
!     call special tracer decoupling routine for multiple (ntr) species
!
      do n = 1 , ntr
        do j = 1 , jendx
          do k = 1 , kz
            do i = 1 , iym1
              chi(i,k,j,n) = chia(i,k,j,n)/sps1%ps(j,i)
            end do
          end do
        end do
      end do
    end if
!
!=======================================================================
!
#ifndef BAND
    if ( myid /= nproc-1 ) then
#endif
      do i = 1 , iy
        tvar1snd(i,1) = sps2%ps(jxp,i)
      end do
      do k = 1 , kz
        do i = 1 , iy
          tvar1snd(i,1+k)       = atm2%t(i,k,jxp)
          tvar1snd(i,1+kz+k)    = atm2%qv(i,k,jxp)
          tvar1snd(i,1+kz*2+k)  = atm2%u(i,k,jxp)
          tvar1snd(i,1+kz*3+k)  = atm2%v(i,k,jxp)
          tvar1snd(i,1+kz*4+k)  = atmx%u(i,k,jxp)
          tvar1snd(i,1+kz*5+k)  = atmx%v(i,k,jxp)
          tvar1snd(i,1+kz*6+k)  = atmx%t(i,k,jxp)
          tvar1snd(i,1+kz*7+k)  = atmx%qv(i,k,jxp)
          tvar1snd(i,1+kz*8+k)  = atmx%qc(i,k,jxp)
          tvar1snd(i,1+kz*9+k)  = atm1%u(i,k,jxp)
          tvar1snd(i,1+kz*10+k) = atm1%v(i,k,jxp)
        end do
      end do
      numrec = kz*11+1
      if ( ichem == 1 ) then
        do n = 1 , ntr
          do k = 1 , kz
            do i = 1 , iy
              tvar1snd(i,numrec+(n-1)*2*kz+k)    = chi(i,k,jxp,n)
              tvar1snd(i,numrec+(n-1)*2*kz+kz+k) = chib(i,k,jxp,n)
            end do
          end do
        end do
        numrec = numrec + ntr*kz*2
      end if
      if ( ibltyp == 2 .or. ibltyp == 99 ) then
        do k = 1 , kz
          do i = 1 , iy
            tvar1snd(i,numrec+k) = atm1%tke(i,k,jxp)
          end do
        end do
        numrec = numrec + kz
      end if
#ifndef BAND
    else
      numrec = kz*11+1
      if ( ichem == 1 ) numrec = numrec + ntr*kz*2
      if ( ibltyp == 2 .or. ibltyp == 99 ) numrec = numrec + kz
    end if
#endif
    call mpi_sendrecv(tvar1snd,iy*numrec,mpi_real8,ieast,1, &
                      tvar1rcv,iy*numrec,mpi_real8,iwest,1, &
                      mycomm,mpi_status_ignore,ierr)
#ifndef BAND
    if ( myid /= 0 ) then
#endif
      do i = 1 , iy
        sps2%ps(0,i) = tvar1rcv(i,1)
      end do
      do k = 1 , kz
        do i = 1 , iy
          atm2%t(i,k,0) = tvar1rcv(i,1+k)
          atm2%qv(i,k,0) = tvar1rcv(i,1+kz+k)
          atm2%u(i,k,0) = tvar1rcv(i,1+kz*2+k)
          atm2%v(i,k,0) = tvar1rcv(i,1+kz*3+k)
          atmx%u(i,k,0) = tvar1rcv(i,1+kz*4+k)
          atmx%v(i,k,0) = tvar1rcv(i,1+kz*5+k)
          atmx%t(i,k,0) = tvar1rcv(i,1+kz*6+k)
          atmx%qv(i,k,0) = tvar1rcv(i,1+kz*7+k)
          atmx%qc(i,k,0) = tvar1rcv(i,1+kz*8+k)
          atm1%u(i,k,0) = tvar1rcv(i,1+kz*9+k)
          atm1%v(i,k,0) = tvar1rcv(i,1+kz*10+k)
        end do
      end do
      numrec = kz*11+1
      if ( ichem == 1 ) then
        do n = 1 , ntr
          do k = 1 , kz
            do i = 1 , iy
              chi(i,k,0,n) = tvar1rcv(i,numrec+(n-1)*2*kz+k)
              chib(i,k,0,n) = tvar1rcv(i,numrec+(n-1)*2*kz+kz+k)
            end do
          end do
        end do
        numrec = numrec + ntr*kz*2
      end if
      if ( ibltyp == 2 .or. ibltyp == 99 ) then
        do k = 1 , kz
          do i = 1 , iy
            atm2%tke(i,k,0) = tvar1rcv(i,numrec+k)
          end do
        end do
      end if
#ifndef BAND
    end if
    if ( myid /= 0 ) then
#endif
!
      do i = 1 , iy
        tvar1snd(i,1) = sps2%ps(1,i)
      end do
      do k = 1 , kz
        do i = 1 , iy
          tvar1snd(i,1+k) = atm2%t(i,k,1)
          tvar1snd(i,1+kz+k) = atm2%qv(i,k,1)
          tvar1snd(i,1+kz*2+k) = atm2%u(i,k,1)
          tvar1snd(i,1+kz*3+k) = atm2%v(i,k,1)
          tvar1snd(i,1+kz*4+k) = atmx%u(i,k,1)
          tvar1snd(i,1+kz*5+k) = atmx%v(i,k,1)
          tvar1snd(i,1+kz*6+k) = atmx%t(i,k,1)
          tvar1snd(i,1+kz*7+k) = atmx%qv(i,k,1)
          tvar1snd(i,1+kz*8+k) = atmx%qc(i,k,1)
          tvar1snd(i,1+kz*9+k) = atm1%u(i,k,1)
          tvar1snd(i,1+kz*10+k) = atm1%v(i,k,1)
        end do
      end do
      numrec = kz*11+1
      if ( ichem == 1 ) then
        do n = 1 , ntr
          do k = 1 , kz
            do i = 1 , iy
              tvar1snd(i,numrec+(n-1)*kz*2+k) = chi(i,k,1,n)
              tvar1snd(i,numrec+(n-1)*kz*2+kz+k) = chib(i,k,1,n)
            end do
          end do
        end do
        numrec = numrec + ntr*kz*2
      end if
      if ( ibltyp == 2 .or. ibltyp == 99 ) then
        do k = 1 , kz
          do i = 1 , iy
            tvar1snd(i,numrec+k) = atm1%tke(i,k,1)
          end do
        end do
        numrec = numrec + kz
      end if
#ifndef BAND
    else
      numrec = kz*11+1
      if ( ichem == 1 ) numrec = numrec + ntr*kz*2
      if ( ibltyp == 2 .or. ibltyp == 99 ) numrec = numrec + kz
    end if
#endif
    call mpi_sendrecv(tvar1snd,iy*numrec,mpi_real8,iwest,2, &
                      tvar1rcv,iy*numrec,mpi_real8,ieast,2, &
                      mycomm,mpi_status_ignore,ierr)
#ifndef BAND
    if ( myid /= nproc-1 ) then
#endif
      do i = 1 , iy
        sps2%ps(jxp+1,i) = tvar1rcv(i,1)
      end do
      do k = 1 , kz
        do i = 1 , iy
          atm2%t(i,k,jxp+1) = tvar1rcv(i,1+k)
          atm2%qv(i,k,jxp+1) = tvar1rcv(i,1+kz+k)
          atm2%u(i,k,jxp+1) = tvar1rcv(i,1+kz*2+k)
          atm2%v(i,k,jxp+1) = tvar1rcv(i,1+kz*3+k)
          atmx%u(i,k,jxp+1) = tvar1rcv(i,1+kz*4+k)
          atmx%v(i,k,jxp+1) = tvar1rcv(i,1+kz*5+k)
          atmx%t(i,k,jxp+1) = tvar1rcv(i,1+kz*6+k)
          atmx%qv(i,k,jxp+1) = tvar1rcv(i,1+kz*7+k)
          atmx%qc(i,k,jxp+1) = tvar1rcv(i,1+kz*8+k)
          atm1%u(i,k,jxp+1) = tvar1rcv(i,1+kz*9+k)
          atm1%v(i,k,jxp+1) = tvar1rcv(i,1+kz*10+k)
        end do
      end do
      numrec = kz*11+1
      if ( ichem == 1 ) then
        do n = 1 , ntr
          do k = 1 , kz
            do i = 1 , iy
              chi(i,k,jxp+1,n) = tvar1rcv(i,numrec+(n-1)*kz*2+k)
              chib(i,k,jxp+1,n) = tvar1rcv(i,numrec+(n-1)*kz*2+kz+k)
            end do
          end do
        end do
        numrec = numrec + ntr*kz*2
      end if
      if ( ibltyp == 2 .or. ibltyp == 99 ) then
        do k = 1 , kz
          do i = 1 , iy
            atm1%tke(i,k,jxp+1) = tvar1rcv(i,numrec+k)
          end do
        end do
      end if
#ifndef BAND
    end if
#endif
!
!=======================================================================
!
!   interior points:
!
    do j = jbegin , jendx
      jm1 = j-1
      do i = 2 , iym1
        sps2%pdot(j,i)=(sps2%ps(j,i)+sps2%ps(j,i-1)+ &
                        sps2%ps(jm1,i)+sps2%ps(jm1,i-1))*d_rfour
      end do
    end do
!
!   east and west boundaries:
!
#ifndef BAND
    do i = 2 , iym1
      if ( myid == 0 )  &
        sps2%pdot(1,i) = (sps2%ps(1,i)+sps2%ps(1,i-1))*d_half
      if ( myid == nproc-1 ) &
        sps2%pdot(jendl,i) = (sps2%ps(jendx,i)+ &
                              sps2%ps(jendx,i-1))*d_half
    end do
#endif
!
!   north and south boundaries:
!
    do j = jbegin , jendx
      jm1 = j-1
      sps2%pdot(j,1)  = (sps2%ps(j,1)+sps2%ps(jm1,1))*d_half
      sps2%pdot(j,iy) = (sps2%ps(j,iym1)+sps2%ps(jm1,iym1))*d_half
    end do
!
!   corner points:
!
#ifndef BAND
    if ( myid == 0 ) then
      sps2%pdot(1,1) = sps2%ps(1,1)
      sps2%pdot(1,iy) = sps2%ps(1,iym1)
    end if
    if ( myid == nproc-1 ) then
      sps2%pdot(jendl,1) = sps2%ps(jendx,1)
      sps2%pdot(jendl,iy) = sps2%ps(jendx,iym1)
    end if
#endif
!
!=======================================================================
!
    call mkslice(1,jendx,1,iym1)

#ifdef CLM
    if ( init_grid ) then
      call initclm(ifrest,idate1,idate2,dx,dtrad,dtsrf)
      init_grid = .false.
    end if
#endif
!
!=======================================================================
!
#ifndef BAND
    if ( myid /= nproc-1 ) then
#endif
      do k = 1 , kz
        do i = 1 , iy
          var2snd(i,+k)     = atms%ubd3d(jxp-1,i,k)
          var2snd(i,kz+k)   = atms%ubd3d(jxp,i,k)
          var2snd(i,kz*2+k) = atms%vbd3d(jxp-1,i,k)
          var2snd(i,kz*3+k) = atms%vbd3d(jxp,i,k)
          var2snd(i,kz*4+k) = atms%tb3d(jxp-1,i,k)
          var2snd(i,kz*5+k) = atms%tb3d(jxp,i,k)
          var2snd(i,kz*6+k) = atms%qvb3d(jxp-1,i,k)
          var2snd(i,kz*7+k) = atms%qvb3d(jxp,i,k)
          var2snd(i,kz*8+k) = atms%qcb3d(jxp-1,i,k)
          var2snd(i,kz*9+k) = atms%qcb3d(jxp,i,k)
        end do
      end do
      numrec = kz*10
      if ( ichem == 1 ) then
        do n = 1 , ntr
          do k = 1 , kz
            do i = 1 , iy
              var2snd(i,numrec+(n-1)*2*kz+k)    = atms%chib3d(jxp-1,i,k,n)
              var2snd(i,numrec+(n-1)*2*kz+kz+k) = atms%chib3d(jxp,i,k,n)
            end do
          end do
        end do
        numrec = numrec + kz*ntr*2
      end if
      if ( ibltyp == 2 .or. ibltyp == 99 ) then
        do k = 1 , kz
          do i = 1 , iy
            var2snd(i,numrec+k)    = atm2%tke(i,k,jxp-1)
            var2snd(i,numrec+kz+k) = atm2%tke(i,k,jxp)
          end do
        end do
        numrec = numrec + kz*2
      end if
#ifndef BAND
    else
      numrec = kz*10
      if ( ichem == 1 ) numrec = numrec + ntr*kz*2
      if ( ibltyp == 2 .or. ibltyp == 99 ) numrec = numrec + kz*2
    end if
#endif
    call mpi_sendrecv(var2snd,iy*numrec,mpi_real8,ieast,1, &
                      var2rcv,iy*numrec,mpi_real8,iwest,1, &
                      mycomm,mpi_status_ignore,ierr)
#ifndef BAND
    if ( myid /= 0 ) then
#endif
      do k = 1 , kz
        do i = 1 , iy
          atms%ubd3d(-1,i,k) = var2rcv(i,+k)
          atms%ubd3d(0,i,k) = var2rcv(i,kz+k)
          atms%vbd3d(-1,i,k) = var2rcv(i,kz*2+k)
          atms%vbd3d(0,i,k) = var2rcv(i,kz*3+k)
          atms%tb3d(-1,i,k) = var2rcv(i,kz*4+k)
          atms%tb3d(0,i,k) = var2rcv(i,kz*5+k)
          atms%qvb3d(-1,i,k) = var2rcv(i,kz*6+k)
          atms%qvb3d(0,i,k) = var2rcv(i,kz*7+k)
          atms%qcb3d(-1,i,k) = var2rcv(i,kz*8+k)
          atms%qcb3d(0,i,k) = var2rcv(i,kz*9+k)
        end do
      end do
      numrec = kz*10
      if ( ichem == 1 ) then
        do n = 1 , ntr
          do k = 1 , kz
            do i = 1 , iy
              atms%chib3d(-1,i,k,n) = var2rcv(i,kz*10+(n-1)*2*kz+k)
              atms%chib3d(0,i,k,n) = var2rcv(i,kz*10+(n-1)*2*kz+kz+k)
            end do
          end do
        end do
        numrec = numrec + ntr*kz*2
      end if
      if ( ibltyp == 2 .or. ibltyp == 99 ) then
        do k = 1 , kz
          do i = 1 , iy
            atm2%tke(i,k,-1) = var2rcv(i,numrec+k)
            atm2%tke(i,k,0)  = var2rcv(i,numrec+kz+k)
          end do
        end do
      end if
#ifndef BAND
    end if
!
    if ( myid /= 0 ) then
#endif
      do k = 1 , kz
        do i = 1 , iy
          var2snd(i,+k) = atms%ubd3d(1,i,k)
          var2snd(i,kz+k) = atms%ubd3d(2,i,k)
          var2snd(i,kz*2+k) = atms%vbd3d(1,i,k)
          var2snd(i,kz*3+k) = atms%vbd3d(2,i,k)
          var2snd(i,kz*4+k) = atms%tb3d(1,i,k)
          var2snd(i,kz*5+k) = atms%tb3d(2,i,k)
          var2snd(i,kz*6+k) = atms%qvb3d(1,i,k)
          var2snd(i,kz*7+k) = atms%qvb3d(2,i,k)
          var2snd(i,kz*8+k) = atms%qcb3d(1,i,k)
          var2snd(i,kz*9+k) = atms%qcb3d(2,i,k)
        end do
      end do
      numrec = kz*10
      if ( ichem == 1 ) then
        do n = 1 , ntr
          do k = 1 , kz
            do i = 1 , iy
              var2snd(i,numrec+(n-1)*2*kz+k) = atms%chib3d(1,i,k,n)
              var2snd(i,numrec+(n-1)*2*kz+kz+k) = atms%chib3d(2,i,k,n)
            end do
          end do
        end do
        numrec = numrec + ntr*kz*2
      end if
      if ( ibltyp == 2 .or. ibltyp == 99 ) then
        do k = 1 , kz
          do i = 1 , iy
            var2snd(i,numrec+k)    = atm2%tke(i,k,1)
            var2snd(i,numrec+kz+k) = atm2%tke(i,k,2)
          end do
        end do
        numrec = numrec + kz*2
      end if
#ifndef BAND
    else
      numrec = kz*10
      if ( ichem == 1 ) numrec = numrec + ntr*kz*2
      if ( ibltyp == 2 .or. ibltyp == 99 ) numrec = numrec + kz*2
    end if
#endif
    call mpi_sendrecv(var2snd,iy*numrec,mpi_real8,iwest,2, &
                      var2rcv,iy*numrec,mpi_real8,ieast,2, &
                      mycomm,mpi_status_ignore,ierr)
#ifndef BAND
    if ( myid /= nproc-1 ) then
#endif
      do k = 1 , kz
        do i = 1 , iy
          atms%ubd3d(jxp+1,i,k) = var2rcv(i,+k)
          atms%ubd3d(jxp+2,i,k) = var2rcv(i,kz+k)
          atms%vbd3d(jxp+1,i,k) = var2rcv(i,kz*2+k)
          atms%vbd3d(jxp+2,i,k) = var2rcv(i,kz*3+k)
          atms%tb3d(jxp+1,i,k) = var2rcv(i,kz*4+k)
          atms%tb3d(jxp+2,i,k) = var2rcv(i,kz*5+k)
          atms%qvb3d(jxp+1,i,k) = var2rcv(i,kz*6+k)
          atms%qvb3d(jxp+2,i,k) = var2rcv(i,kz*7+k)
          atms%qcb3d(jxp+1,i,k) = var2rcv(i,kz*8+k)
          atms%qcb3d(jxp+2,i,k) = var2rcv(i,kz*9+k)
        end do
      end do
      numrec = kz*10
      if ( ichem == 1 ) then
        do n = 1 , ntr
          do k = 1 , kz
            do i = 1 , iy
              atms%chib3d(jxp+1,i,k,n) = var2rcv(i,numrec+(n-1)*2*kz+k)
              atms%chib3d(jxp+2,i,k,n) = var2rcv(i,numrec+(n-1)*2*kz+kz+k)
            end do
          end do
        end do
        numrec = numrec + ntr*kz*2
      end if
      if ( ibltyp == 2 .or. ibltyp == 99 ) then
        do k = 1 , kz
          do i = 1 , iy
            atm2%tke(i,k,jxp+1) = var2rcv(i,numrec+k)
            atm2%tke(i,k,jxp+2) = var2rcv(i,numrec+kz+k)
          end do
        end do
      end if
#ifndef BAND
    end if
#endif
!
!**********************************************************************
!
!  "j" loop begins here:
!
!**********************************************************************
!
#ifndef BAND
    do j = 1 , jendx
      if ( (myid == 0 .and. j == 1) .or.  &
           (myid == nproc-1 .and. j == jendx) ) then
        do k = 1 , kzp1
          do i = 1 , iym1
            qdot(i,k,j) = d_zero
          end do
        end do
        do k = 1 , kz
          do i = 1 , iym1
            omega(i,k,j) = d_zero
          end do
        end do
      end if
    end do
#endif

    do j = 1 , jendx
      jp1 = j+1
!
#ifndef BAND
      if ((.not.(myid == 0 .and. j == 1)) .and. &
         (.not.(myid == nproc-1 .and. j == jendx)) ) then
#endif
!
!----------------------------------------------------------------------
!
!     compute the pressure tendency:
!
      do i = 2 , iym2
         pten(i,j) = d_zero
      end do
      do k = 1 , kz
         do i = 2 , iym2
            divl(i,k) = (atm1%u(i+1,k,jp1)+atm1%u(i,k,jp1)- &
                         atm1%u(i+1,k,j)-atm1%u(i,k,j))+    &
                        (atm1%v(i+1,k,jp1)+atm1%v(i+1,k,j)- &
                         atm1%v(i,k,jp1)-atm1%v(i,k,j))
            pten(i,j) = pten(i,j) - divl(i,k)*dsigma(k)     &
                        /(dx2*mddom%msfx(i,j)*mddom%msfx(i,j))
         end do
      end do
!
!     compute vertical sigma-velocity (qdot):
!
      do k = 1 , kzp1
         do i = 1 , iym1
            qdot(i,k,j) = d_zero
         end do
      end do
      do k = 2 , kz
         do i = 2 , iym2
            qdot(i,k,j) = qdot(i,k-1,j)                               &
                       - (pten(i,j)+divl(i,k-1)/(dx2*mddom%msfx(i,j) &
                       *mddom%msfx(i,j)))*dsigma(k-1)/sps1%ps(j,i)
         end do
      end do
#ifndef BAND
      end if
#endif
    end do
    call mpi_sendrecv(qdot(:,:,jxp),iy*kzp1,mpi_real8,ieast,1, &
                      qdot(:,:,0),  iy*kzp1,mpi_real8,iwest,1, &
                      mycomm,mpi_status_ignore,ierr)
    call mpi_sendrecv(qdot(:,:,1),    iy*kzp1,mpi_real8,iwest,2, &
                      qdot(:,:,jxp+1),iy*kzp1,mpi_real8,ieast,2, &
                      mycomm,mpi_status_ignore,ierr)
!
!   compute omega
!
    do j = 1 , jendx
      jp1 = j+1
      jm1 = j-1
#ifndef BAND
      if ((.not.(myid == 0 .and. j == 1)) .and. &
         (.not.(myid == nproc-1 .and. j == jendx)) ) then
#endif
      do k = 1 , kz
         do i = 2 , iym2
            omega(i,k,j) = d_half*sps1%ps(j,i)* &
                       (qdot(i,k+1,j)+qdot(i,k,j))+a(k)*(pten(i,j)+   &
                       ((atmx%u(i,k,j)+atmx%u(i+1,k,j)+               &
                         atmx%u(i+1,k,jp1)+atmx%u(i,k,jp1))*          &
                       (sps1%ps(jp1,i)-sps1%ps(jm1,i))+               &
                       (atmx%v(i,k,j)+atmx%v(i+1,k,j)+                &
                        atmx%v(i+1,k,jp1)+atmx%v(i,k,jp1))*           &
                       (sps1%ps(j,i+1)-sps1%ps(j,i-1)))/              &
                       (dx8*mddom%msfx(i,j)))
         end do
      end do
#ifndef BAND
      end if
#endif
    end do

#ifndef BAND
    if ( nspgx >= jxp ) then
      do i = 1 , iy
        bdyewsnd(i,1) = xpsb%eb(i,1)
        bdyewsnd(i,2) = xpsb%wb(i,jxp)
        bdyewsnd(i,3) = xpsb%ebt(i,1)
        bdyewsnd(i,4) = xpsb%wbt(i,jxp)
      end do
      do k = 1 , kz
        do i = 1 , iy
          bdyewsnd(i,4+k) = xtb%eb(i,k,1)
          bdyewsnd(i,4+kz+k) = xtb%wb(i,k,jxp)
          bdyewsnd(i,4+kz*2+k) = xtb%ebt(i,k,1)
          bdyewsnd(i,4+kz*3+k) = xtb%wbt(i,k,jxp)
          bdyewsnd(i,4+kz*4+k) = xqb%eb(i,k,1)
          bdyewsnd(i,4+kz*5+k) = xqb%wb(i,k,jxp)
          bdyewsnd(i,4+kz*6+k) = xqb%ebt(i,k,1)
          bdyewsnd(i,4+kz*7+k) = xqb%wbt(i,k,jxp)
          bdyewsnd(i,4+kz*8+k) = xub%eb(i,k,1)
          bdyewsnd(i,4+kz*9+k) = xub%wb(i,k,jxp)
          bdyewsnd(i,4+kz*10+k) = xub%ebt(i,k,1)
          bdyewsnd(i,4+kz*11+k) = xub%wbt(i,k,jxp)
          bdyewsnd(i,4+kz*12+k) = xvb%eb(i,k,1)
          bdyewsnd(i,4+kz*13+k) = xvb%wb(i,k,jxp)
          bdyewsnd(i,4+kz*14+k) = xvb%ebt(i,k,1)
          bdyewsnd(i,4+kz*15+k) = xvb%wbt(i,k,jxp)
        end do
      end do
      call mpi_sendrecv(bdyewsnd,iy*(kz*16+4),mpi_real8,ieast,1, &
                        bdyewrcv,iy*(kz*16+4),mpi_real8,iwest,1, &
                        mycomm,mpi_status_ignore,ierr)
      do i = 1 , iy
        if ( myid == nproc-1 ) then
          xpsb%eb(i,jendl) = bdyewrcv(i,1)
          xpsb%ebt(i,jendl) = bdyewrcv(i,3)
        else
          xpsb%eb(i,jxp+1) = bdyewrcv(i,1)
          xpsb%ebt(i,jxp+1) = bdyewrcv(i,3)
        end if
        xpsb%wb(i,0) = bdyewrcv(i,2)
        xpsb%wbt(i,0) = bdyewrcv(i,4)
      end do
      do k = 1 , kz
        do i = 1 , iy
          if ( myid == nproc-1 ) then
            xtb%eb(i,k,jendl) = bdyewrcv(i,4+k)
            xtb%ebt(i,k,jendl) = bdyewrcv(i,4+kz*2+k)
            xqb%eb(i,k,jendl) = bdyewrcv(i,4+kz*4+k)
            xqb%ebt(i,k,jendl) = bdyewrcv(i,4+kz*6+k)
          else
            xtb%eb(i,k,jxp+1) = bdyewrcv(i,4+k)
            xtb%ebt(i,k,jxp+1) = bdyewrcv(i,4+kz*2+k)
            xqb%eb(i,k,jxp+1) = bdyewrcv(i,4+kz*4+k)
            xqb%ebt(i,k,jxp+1) = bdyewrcv(i,4+kz*6+k)
          end if
          xub%eb(i,k,jxp+1) = bdyewrcv(i,4+kz*8+k)
          xub%ebt(i,k,jxp+1) = bdyewrcv(i,4+kz*10+k)
          xvb%eb(i,k,jxp+1) = bdyewrcv(i,4+kz*12+k)
          xvb%ebt(i,k,jxp+1) = bdyewrcv(i,4+kz*14+k)
          xtb%wb(i,k,0) = bdyewrcv(i,4+kz+k)
          xtb%wbt(i,k,0) = bdyewrcv(i,4+kz*3+k)
          xqb%wb(i,k,0) = bdyewrcv(i,4+kz*5+k)
          xqb%wbt(i,k,0) = bdyewrcv(i,4+kz*7+k)
          xub%wb(i,k,0) = bdyewrcv(i,4+kz*9+k)
          xub%wbt(i,k,0) = bdyewrcv(i,4+kz*11+k)
          xvb%wb(i,k,0) = bdyewrcv(i,4+kz*13+k)
          xvb%wbt(i,k,0) = bdyewrcv(i,4+kz*15+k)
        end do
      end do
      do i = 1 , iy
        if ( myid == nproc-1 ) then
          bdyewsnd(i,1) = xpsb%eb(i,jendx)
          bdyewsnd(i,3) = xpsb%ebt(i,jendx)
        else
          bdyewsnd(i,1) = xpsb%eb(i,jxp)
          bdyewsnd(i,3) = xpsb%ebt(i,jxp)
        end if
        bdyewsnd(i,2) = xpsb%wb(i,1)
        bdyewsnd(i,4) = xpsb%wbt(i,1)
      end do
      do k = 1 , kz
        do i = 1 , iy
          if ( myid == nproc-1 ) then
            bdyewsnd(i,4+k) = xtb%eb(i,k,jendx)
            bdyewsnd(i,4+kz*2+k) = xtb%ebt(i,k,jendx)
            bdyewsnd(i,4+kz*4+k) = xqb%eb(i,k,jendx)
            bdyewsnd(i,4+kz*6+k) = xqb%ebt(i,k,jendx)
          else
            bdyewsnd(i,4+k) = xtb%eb(i,k,jxp)
            bdyewsnd(i,4+kz*2+k) = xtb%ebt(i,k,jxp)
            bdyewsnd(i,4+kz*4+k) = xqb%eb(i,k,jxp)
            bdyewsnd(i,4+kz*6+k) = xqb%ebt(i,k,jxp)
          end if
          bdyewsnd(i,4+kz*8+k) = xub%eb(i,k,jxp)
          bdyewsnd(i,4+kz*10+k) = xub%ebt(i,k,jxp)
          bdyewsnd(i,4+kz*12+k) = xvb%eb(i,k,jxp)
          bdyewsnd(i,4+kz*14+k) = xvb%ebt(i,k,jxp)
          bdyewsnd(i,4+kz+k) = xtb%wb(i,k,1)
          bdyewsnd(i,4+kz*3+k) = xtb%wbt(i,k,1)
          bdyewsnd(i,4+kz*5+k) = xqb%wb(i,k,1)
          bdyewsnd(i,4+kz*7+k) = xqb%wbt(i,k,1)
          bdyewsnd(i,4+kz*9+k) = xub%wb(i,k,1)
          bdyewsnd(i,4+kz*11+k) = xub%wbt(i,k,1)
          bdyewsnd(i,4+kz*13+k) = xvb%wb(i,k,1)
          bdyewsnd(i,4+kz*15+k) = xvb%wbt(i,k,1)
        end do
      end do
      call mpi_sendrecv(bdyewsnd,iy*(kz*16+4),mpi_real8,iwest,2, &
                        bdyewrcv,iy*(kz*16+4),mpi_real8,ieast,2, &
                        mycomm,mpi_status_ignore,ierr)
      do i = 1 , iy
        xpsb%eb(i,0) = bdyewrcv(i,1)
        xpsb%ebt(i,0) = bdyewrcv(i,3)
        xpsb%wb(i,jxp+1) = bdyewrcv(i,2)
        xpsb%wbt(i,jxp+1) = bdyewrcv(i,4)
      end do
      do k = 1 , kz
        do i = 1 , iy
          xtb%eb(i,k,0) = bdyewrcv(i,4+k)
          xtb%wb(i,k,jxp+1) = bdyewrcv(i,4+kz+k)
          xtb%ebt(i,k,0) = bdyewrcv(i,4+kz*2+k)
          xtb%wbt(i,k,jxp+1) = bdyewrcv(i,4+kz*3+k)
          xqb%eb(i,k,0) = bdyewrcv(i,4+kz*4+k)
          xqb%wb(i,k,jxp+1) = bdyewrcv(i,4+kz*5+k)
          xqb%ebt(i,k,0) = bdyewrcv(i,4+kz*6+k)
          xqb%wbt(i,k,jxp+1) = bdyewrcv(i,4+kz*7+k)
          xub%eb(i,k,0) = bdyewrcv(i,4+kz*8+k)
          xub%wb(i,k,jxp+1) = bdyewrcv(i,4+kz*9+k)
          xub%ebt(i,k,0) = bdyewrcv(i,4+kz*10+k)
          xub%wbt(i,k,jxp+1) = bdyewrcv(i,4+kz*11+k)
          xvb%eb(i,k,0) = bdyewrcv(i,4+kz*12+k)
          xvb%wb(i,k,jxp+1) = bdyewrcv(i,4+kz*13+k)
          xvb%ebt(i,k,0) = bdyewrcv(i,4+kz*14+k)
          xvb%wbt(i,k,jxp+1) = bdyewrcv(i,4+kz*15+k)
        end do
      end do
    end if
#endif
!
#ifndef BAND
    if ( myid /= nproc-1 ) then
#endif
      do i = 1 , nspgx
        bdynssnd(i,1) = xpsb%nb(i,jxp)
        bdynssnd(i,2) = xpsb%nbt(i,jxp)
        bdynssnd(i,3) = xpsb%sb(i,jxp)
        bdynssnd(i,4) = xpsb%sbt(i,jxp)
      end do
      do k = 1 , kz
        do i = 1 , nspgx
          bdynssnd(i,4+k) = xtb%nb(i,k,jxp)
          bdynssnd(i,4+kz+k) = xtb%nbt(i,k,jxp)
          bdynssnd(i,4+kz*2+k) = xtb%sb(i,k,jxp)
          bdynssnd(i,4+kz*3+k) = xtb%sbt(i,k,jxp)
          bdynssnd(i,4+kz*4+k) = xqb%nb(i,k,jxp)
          bdynssnd(i,4+kz*5+k) = xqb%nbt(i,k,jxp)
          bdynssnd(i,4+kz*6+k) = xqb%sb(i,k,jxp)
          bdynssnd(i,4+kz*7+k) = xqb%sbt(i,k,jxp)
          bdynssnd(i,4+kz*8+k) = xub%nb(i,k,jxp)
          bdynssnd(i,4+kz*9+k) = xub%nbt(i,k,jxp)
          bdynssnd(i,4+kz*10+k) = xub%sb(i,k,jxp)
          bdynssnd(i,4+kz*11+k) = xub%sbt(i,k,jxp)
          bdynssnd(i,4+kz*12+k) = xvb%nb(i,k,jxp)
          bdynssnd(i,4+kz*13+k) = xvb%nbt(i,k,jxp)
          bdynssnd(i,4+kz*14+k) = xvb%sb(i,k,jxp)
          bdynssnd(i,4+kz*15+k) = xvb%sbt(i,k,jxp)
        end do
      end do
#ifndef BAND
    end if
#endif
    call mpi_sendrecv(bdynssnd,nspgx*(kz*16+4),mpi_real8,ieast,1, &
                      bdynsrcv,nspgx*(kz*16+4),mpi_real8,iwest,1, &
                      mycomm,mpi_status_ignore,ierr)
#ifndef BAND
    if ( myid /= 0 ) then
#endif
      do i = 1 , nspgx
        xpsb%nb(i,0) = bdynsrcv(i,1)
        xpsb%nbt(i,0) = bdynsrcv(i,2)
        xpsb%sb(i,0) = bdynsrcv(i,3)
        xpsb%sbt(i,0) = bdynsrcv(i,4)
      end do
      do k = 1 , kz
        do i = 1 , nspgx
          xtb%nb(i,k,0) = bdynsrcv(i,4+k)
          xtb%nbt(i,k,0) = bdynsrcv(i,4+kz+k)
          xtb%sb(i,k,0) = bdynsrcv(i,4+kz*2+k)
          xtb%sbt(i,k,0) = bdynsrcv(i,4+kz*3+k)
          xqb%nb(i,k,0) = bdynsrcv(i,4+kz*4+k)
          xqb%nbt(i,k,0) = bdynsrcv(i,4+kz*5+k)
          xqb%sb(i,k,0) = bdynsrcv(i,4+kz*6+k)
          xqb%sbt(i,k,0) = bdynsrcv(i,4+kz*7+k)
          xub%nb(i,k,0) = bdynsrcv(i,4+kz*8+k)
          xub%nbt(i,k,0) = bdynsrcv(i,4+kz*9+k)
          xub%sb(i,k,0) = bdynsrcv(i,4+kz*10+k)
          xub%sbt(i,k,0) = bdynsrcv(i,4+kz*11+k)
          xvb%nb(i,k,0) = bdynsrcv(i,4+kz*12+k)
          xvb%nbt(i,k,0) = bdynsrcv(i,4+kz*13+k)
          xvb%sb(i,k,0) = bdynsrcv(i,4+kz*14+k)
          xvb%sbt(i,k,0) = bdynsrcv(i,4+kz*15+k)
        end do
      end do
#ifndef BAND
    end if
    if ( myid /= 0 ) then
#endif
      do i = 1 , nspgx
        bdynssnd(i,1) = xpsb%nb(i,1)
        bdynssnd(i,2) = xpsb%nbt(i,1)
        bdynssnd(i,3) = xpsb%sb(i,1)
        bdynssnd(i,4) = xpsb%sbt(i,1)
      end do
      do k = 1 , kz
        do i = 1 , nspgx
          bdynssnd(i,4+k) = xtb%nb(i,k,1)
          bdynssnd(i,4+kz+k) = xtb%nbt(i,k,1)
          bdynssnd(i,4+kz*2+k) = xtb%sb(i,k,1)
          bdynssnd(i,4+kz*3+k) = xtb%sbt(i,k,1)
          bdynssnd(i,4+kz*4+k) = xqb%nb(i,k,1)
          bdynssnd(i,4+kz*5+k) = xqb%nbt(i,k,1)
          bdynssnd(i,4+kz*6+k) = xqb%sb(i,k,1)
          bdynssnd(i,4+kz*7+k) = xqb%sbt(i,k,1)
          bdynssnd(i,4+kz*8+k) = xub%nb(i,k,1)
          bdynssnd(i,4+kz*9+k) = xub%nbt(i,k,1)
          bdynssnd(i,4+kz*10+k) = xub%sb(i,k,1)
          bdynssnd(i,4+kz*11+k) = xub%sbt(i,k,1)
          bdynssnd(i,4+kz*12+k) = xvb%nb(i,k,1)
          bdynssnd(i,4+kz*13+k) = xvb%nbt(i,k,1)
          bdynssnd(i,4+kz*14+k) = xvb%sb(i,k,1)
          bdynssnd(i,4+kz*15+k) = xvb%sbt(i,k,1)
        end do
      end do
#ifndef BAND
    end if
#endif
    call mpi_sendrecv(bdynssnd,nspgx*(kz*16+4),mpi_real8,iwest,2, &
                      bdynsrcv,nspgx*(kz*16+4),mpi_real8,ieast,2, &
                      mycomm,mpi_status_ignore,ierr)
#ifndef BAND
    if ( myid /= nproc-1 ) then
#endif
      do i = 1 , nspgx
        xpsb%nb(i,jxp+1) = bdynsrcv(i,1)
        xpsb%nbt(i,jxp+1) = bdynsrcv(i,2)
        xpsb%sb(i,jxp+1) = bdynsrcv(i,3)
        xpsb%sbt(i,jxp+1) = bdynsrcv(i,4)
      end do
      do k = 1 , kz
        do i = 1 , nspgx
          xtb%nb(i,k,jxp+1) = bdynsrcv(i,4+k)
          xtb%nbt(i,k,jxp+1) = bdynsrcv(i,4+kz+k)
          xtb%sb(i,k,jxp+1) = bdynsrcv(i,4+kz*2+k)
          xtb%sbt(i,k,jxp+1) = bdynsrcv(i,4+kz*3+k)
          xqb%nb(i,k,jxp+1) = bdynsrcv(i,4+kz*4+k)
          xqb%nbt(i,k,jxp+1) = bdynsrcv(i,4+kz*5+k)
          xqb%sb(i,k,jxp+1) = bdynsrcv(i,4+kz*6+k)
          xqb%sbt(i,k,jxp+1) = bdynsrcv(i,4+kz*7+k)
          xub%nb(i,k,jxp+1) = bdynsrcv(i,4+kz*8+k)
          xub%nbt(i,k,jxp+1) = bdynsrcv(i,4+kz*9+k)
          xub%sb(i,k,jxp+1) = bdynsrcv(i,4+kz*10+k)
          xub%sbt(i,k,jxp+1) = bdynsrcv(i,4+kz*11+k)
          xvb%nb(i,k,jxp+1) = bdynsrcv(i,4+kz*12+k)
          xvb%nbt(i,k,jxp+1) = bdynsrcv(i,4+kz*13+k)
          xvb%sb(i,k,jxp+1) = bdynsrcv(i,4+kz*14+k)
          xvb%sbt(i,k,jxp+1) = bdynsrcv(i,4+kz*15+k)
        end do
      end do
#ifndef BAND
    end if
#endif

    do j = jbegin , jendx
#ifndef BAND
      if ( myid /= nproc-1 .or. j /= jendx ) then
#endif
!
        if ( iboudy == 4 ) then
!         apply sponge boundary conditions to pten:
          call sponge(.false.,ispgx,wgtx,pten,j,1,xpsb)
!       apply the nudging boundary conditions:
        else if ( iboudy == 1 .or. iboudy == 5 ) then
          xtm1 = xbctime - dtsec
          if ( nbdytime == 0 .and. ktau /= 0 ) xtm1 = -dtsec
          call nudge(.false.,ispgx,fnudge,gnudge, &
                     xtm1,sps2%ps,pten,j,1,iboudy,xpsb)
        end if
#ifndef BAND
      end if !end if (j /= jxm1) test
#endif
    end do

#ifndef BAND
    do j = 1 , jendx
      if ( myid == 0 .and. j == 1 ) then
        do i = 1 , iym1
          psc(i,j) = sps2%ps(j,i) + dt*xpsb%wbt(i,j)
          psd(i,j) = sps1%ps(j,i)
        end do
      else if ( myid == nproc-1 .and. j == jendx ) then
        do i = 1 , iym1
          psd(i,j) = sps1%ps(j,i)
        end do
      end if
    end do
#endif

    do j = 1 , jendx
#ifndef BAND
      if ((.not.(myid == 0 .and. j == 1)) .and. &
         (.not.(myid == nproc-1 .and. j == jendx)) ) then
#endif
!
!      forecast pressure:
!
       do i = 2 , iym2
         psc(i,j) = sps2%ps(j,i) + pten(i,j)*dt
       end do
!
!      weighted p* (psd)
!
       do i = 2 , iym2
         psd(i,j) = sps1%ps(j,i)
       end do
!
       psc(1,j) = sps2%ps(j,1) + dt*xpsb%sbt(1,j)
       psc(iym1,j) = sps2%ps(j,iym1) + dt*xpsb%nbt(1,j)
       psd(1,j) = sps1%ps(j,1)
       psd(iym1,j) = sps1%ps(j,iym1)
#ifndef BAND
      end if
#endif
    end do
    call mpi_sendrecv(psd(1,jxp),iy,mpi_real8,ieast,1,                &
                      psd(1,0),iy,mpi_real8,iwest,1,                  &
                      mycomm,mpi_status_ignore,ierr)
!
!   compute bleck (1977) noise parameters:
!
    do j = 1 , jendl
      do i = 1 , iy
        ps4(i,1,j) = pten(i,j)
        ps4(i,2,j) = psc(i,j)
        ps4(i,3,j) = sps2%ps(j,i)
        ps4(i,4,j) = sps1%ps(j,i)
      end do
    end do
    call mpi_gather(ps4, iy*4*jxp,mpi_real8, &
                    ps_4,iy*4*jxp,mpi_real8, &
                    0,mycomm,ierr)
    iptn = 0
    ptntot = d_zero
    pt2tot = d_zero
    if ( myid == 0 ) then
#ifdef BAND
      do j = 1 , jx
#else
      do j = 2 , jxm2
#endif
        if ( ktau /= 0 ) then
          do i = 2 , iym2
            iptn = iptn + 1
            ptntot = ptntot + dabs(ps_4(i,1,j))
            pt2tot = pt2tot +                       &
                     dabs((ps_4(i,2,j)+ps_4(i,3,j)- &
                           d_two*ps_4(i,4,j))/(dt*dt*d_rfour))
          end do
        end if
      end do
    end if
    call mpi_bcast(iptn,1,mpi_integer,0,mycomm,ierr)
    call mpi_bcast(ptntot,1,mpi_real8,0,mycomm,ierr)
    call mpi_bcast(pt2tot,1,mpi_real8,0,mycomm,ierr)
!
!   calculate solar zenith angle
    if ( ktau == 0 .or. ichem == 1 .or. &
         mod(ktau+1,ntsrf) == 0 .or. mod(ktau+1,ntrad) == 0 ) then
      call zenitm(coszrs,jbegin,jendx,1,iy)
    end if
 
    do j = jbegin , jendx
      jp1 = j+1
!
!     compute the horizontal diffusion coefficient and stored in xkc:
!     the values are calculated at cross points, but they also used
!     for dot-point variables.
 
      do k = 1 , kz
        do i = 2 , iym1
          dudx = atm2%u(i,k,jp1) + atm2%u(i+1,k,jp1) - &
                 atm2%u(i,k,j)   - atm2%u(i+1,k,j)
          dvdx = atm2%v(i,k,jp1) + atm2%v(i+1,k,jp1) - &
                 atm2%v(i,k,j)   - atm2%v(i+1,k,j)
          dudy = atm2%u(i+1,k,j) + atm2%u(i+1,k,jp1) - &
                 atm2%u(i,k,j)   - atm2%u(i,k,jp1)
          dvdy = atm2%v(i+1,k,j) + atm2%v(i+1,k,jp1) - &
                 atm2%v(i,k,j)   - atm2%v(i,k,jp1)
!fil      cell=(xkhz*hgfact(i,j)/5.+c200*dsqrt((dudx-dvdy)*(dudx-dvdy)
          cell = (xkhz*hgfact(i,j)                         &
                 +c200*dsqrt((dudx-dvdy)*(dudx-dvdy)+(dvdx+dudy) &
                 *(dvdx+dudy)))
          xkc(i,k,j) = dmin1(cell,xkhmax)
          ! TAO: Interpolate the diffusion coefficients to full levels
          ! for use in diffusion of TKE
          if ( k /= 1 ) then
            cell = twt(k,1)*xkc(i,k,j) + twt(k,2)*xkc(i,k-1,j)
            ! TAO: Multiply the horizontal diffusion coefficient by
            ! nuk for TKE.  Without this multiplication, it appears
            ! that TKE does not diffuse fast enough, and stabilities
            ! appear in the TKE field.  While this is different from
            ! Bretherton's treatment, it is consistent with the
            ! scaling of the vertical TKE diffusivity.
            xkcf(i,k,j) = dmin1(nuk*cell,xkhmax)
          end if
        end do
        xkcf(i,1,j) = nuk*xkc(i,1,j)
        ! No diffusion of TKE on lower boundary
        xkcf(i,kzp1,j) = 0.0d0
      end do
    end do
!

!
!---------------------------------------------------------------------
!
!   compute the temperature tendency:
!
    do k = 1 , kz
      do i = 2 , iym2
        do j = jbegin , jendx
          aten%t(i,k,j) = d_zero
          aten%qv(i,k,j) = d_zero
          aten%qc(i,k,j) = d_zero
        end do
      end do
    end do
 
!
!   compute the horizontal advection term:
!
    call hadv(.false.,aten%t,atmx%t,jbegin,jendx,1)
!
!   compute the vertical advection term:
!
    if ( ibltyp /= 2 .and. ibltyp /= 99 ) then
      call vadv(aten%t,atm1%t,jbegin,jendx,1)
    else
      if ( iuwvadv == 1 ) then
        call vadv(aten%t,atm1%t,jbegin,jendx,6)
      else
        call vadv(aten%t,atm1%t,jbegin,jendx,1)
      end if
    end if
!
!   compute the adiabatic term:
!
    do k = 1 , kz
      do i = 2 , iym2
        do j = jbegin , jendx
          rovcpm = rgas/(cpd*(d_one+0.8D0*(atmx%qv(i,k,j))))
          tv = atmx%t(i,k,j)*(d_one+ep1*(atmx%qv(i,k,j)))
          aten%t(i,k,j) = aten%t(i,k,j) + (omega(i,k,j)*rovcpm*tv) / &
                          (ptop/sps1%ps(j,i)+a(k))
        end do
      end do
    end do
!
!   compute the diffusion term for t and store in difft:
!
    do k = 1 , kz
      do i = 1 , iym1
        do j = jbegin , jendx
          adf%difft(i,k,j) = d_zero
          adf%diffq(i,k,j) = d_zero
        end do
      end do
    end do
!
    call diffu_x(jbegin,jendx,2,iym2,adf%difft,atms%tb3d,sps2%ps,xkc,kz)
!
!   compute the moisture tendencies:
!
!   icup = 1 : kuo-anthes cumulus parameterizaion scheme
!   icup = 2 : grell cumulus paramterization scheme
!   icup = 3 : betts-miller (1986)
!   icup = 4 : emanuel (1991)
!   icup = 5 : tiedtke (1986)
!   icup = 99: grell over land, emanuel over ocean
!   icup = 98: emanuel over land, grell over ocean
!
    call hadv(.false.,aten%qv,atmx%qv,jbegin,jendx,1)
!
    if ( icup /= 1 ) then
      if ( ibltyp /= 2 .and. ibltyp /= 99 ) then
        call vadv(aten%qv,atm1%qv,jbegin,jendx,2)
      else
        if ( iuwvadv == 1 ) then
          call vadv(aten%qv,atm1%qv,jbegin,jendx,6)
        else
          call vadv(aten%qv,atm1%qv,jbegin,jendx,2)
        end if
      end if
    end if
 
    if ( icup == 1 ) then
      call cupara(jbegin,jendx,2,iym2,ktau)
    end if
    if ( icup == 2 .or. icup == 99 .or. icup == 98 ) then
      call cuparan(jbegin,jendx,2,iym2,ktau)
    end if
    if ( icup == 3 ) then
      call bmpara(jbegin,jendx,2,iym2,ktau)
    end if
    if ( icup == 4 .or. icup == 99 .or. icup == 98 ) then
      call cupemandrv(jbegin,jendx,2,iym2,ktau)
    end if
    if ( icup == 5 ) then
      call tiedtkedrv(jbegin,jendx,2,iym2,ktau)
    end if

    if ( ipptls == 1 ) then
      call hadv(.false.,aten%qc,atmx%qc,jbegin,jendx,1)
      if ( ibltyp /= 2 .and. ibltyp /= 99 ) then
        call vadv(aten%qc,atm1%qc,jbegin,jendx,5)
      else
        if ( iuwvadv == 1 ) then
          call vadv(aten%qc,atm1%qc,jbegin,jendx,6)
        else
          call vadv(aten%qc,atm1%qc,jbegin,jendx,5)
        end if
      end if
      call pcp(jbegin,jendx,2,iym1)
      call cldfrac(jbegin,jendx,2,iym1)
!
!     need also to set diffq to 0 here before calling diffut
!
      do k = 1 , kz
        do i = 1 , iym1
          do j = jbegin , jendx
            adf%diffq(i,k,j) = d_zero
          end do
        end do
      end do
 
!     compute the diffusion terms:
!     the diffusion term for qv is stored in diffq. before
!     completing aten%qv computation, do not use diffq for other
!     purpose.
!
      call diffu_x(jbegin,jendx,2,iym2,adf%diffq,atms%qvb3d,sps2%ps,xkc,kz)
      call diffu_x(jbegin,jendx,2,iym2,aten%qc,atms%qcb3d,sps2%ps,xkc,kz)
    end if
!
    if ( ichem == 1 ) then
      !
      ! TRANSPORT OF TRACERS : initialize tracer tendencies
      !
      do j = jbegin , jendx
        do itr = 1 , ntr
          do k = 1 , kz
            do i = 2 , iym1
              chiten(i,k,j,itr) = d_zero
            end do
          end do
        end do
      end do
      !
      ! horizontal and vertical advection
      !
      do itr = 1 , ntr
        ! Here assignpnt does not work with gfortran with a sliced array.
        ! Doing explicit work on bounds.
        spchiten                      => chiten(:,:,:,itr)
        spchi(1:,1:,lbound(chi,3):)   => chi(:,:,:,itr)
        spchia(1:,1:,lbound(chia,3):) => chia(:,:,:,itr)
        spchib3d(1:,1:,lbound(chia,3):) => chib3d(:,:,:,itr)
        call hadv(.false.,spchiten,spchi,jbegin,jendx,2)
        call vadv(spchiten,spchia,jbegin,jendx,5)
        ! horizontal diffusion: initialize scratch vars to 0.
        ! need to compute tracer tendencies due to diffusion
        call diffu_x(jbegin,jendx,2,iym2,spchiten,spchib3d,cpsb,xkc,kz)
      end do ! end tracer loop
      !
      ! Compute tendencies
      !
      call tractend2(jbegin,jendm,2,iym1,ktau,xmonth,xkc)
      !
    end if ! ichem
!
!----------------------------------------------------------------------
!  compute the pbl fluxes:
!  the diffusion and pbl tendencies of t and qv are stored in
!  difft and diffq.
!
   do j = jbegin , jendx
      do k = 1 , kz
        do i = 2 , iym1
          aten%u(i,k,j) = d_zero
          aten%v(i,k,j) = d_zero
        end do
      end do
    end do
 
!     calculate albedo
    if ( ktau == 0 .or. mod(ktau+1,ntrad) == 0 ) then
#ifdef CLM
      call albedoclm(xmonth,jbegin,jendx,2,iym1)
#else
      call albedov(xmonth,jbegin,jendx,2,iym1)
#endif
    end if
 
!   call ccm3 radiative transfer package
    if ( ktau == 0 .or. mod(ktau+1,ntrad) == 0 ) then
      loutrad = (ktau == 0 .or. mod(ktau+1,krad) == 0)
      if (irrtm ==1 ) then
        call rrtmg_driver(jbegin,jendx,2,iym1,ktau,xyear,eccf,loutrad)
      else
        call colmod3(jbegin,jendx,2,iym1,ktau,xyear,eccf,loutrad)
      end if
    end if
 
#ifndef CLM
!   call mtrxbats for surface physics calculations
    if ( ktau == 0 .or. mod(ktau+1,ntsrf) == 0 ) then
      dtbat = dt*d_half*dble(ntsrf)
      if ( ktau == 0 ) dtbat = dt
      call mtrxbats(jbegin,jendx,2,iym1,ktau)
    end if
#endif
 
#ifdef CLM
!   call mtrxclm for surface physics calculations
    if ( ktau == 0 .or. mod(ktau+1,ntrad) == 0 ) then
      r2cdoalb = .true.
    else
      r2cdoalb = .false.
    end if
    if ( ktau == 0 .or. mod(ktau+1,ntsrf) == 0 ) then
      ! Timestep used is the same as for bats
      if ( ktau == 0 ) then
        r2cnstep = 0
      else
        r2cnstep = (ktau+1)/ntsrf
      end if
      dtbat = dt*d_half*ntsrf
      ! CLM j loop is in mtrxclm
      call mtrxclm(ktau)
    end if
#endif

    if ( icup == 1 ) then
      do k = 1 , kz
        do j = 1 , jendl
          do i = 1 , iy
            wrkkuo1(i,j,k) = rsheat(j,i,k)
            wrkkuo2(i,j,k) = rswat(j,i,k)
          end do
        end do
        call mpi_sendrecv(wrkkuo1(:,jxp,k),iy,mpi_real8,ieast,1, &
                          wrkkuo1(:,0,k),iy,mpi_real8,iwest,1,   &
                          mycomm,mpi_status_ignore,ierr)
        call mpi_sendrecv(wrkkuo1(:,1,k),iy,mpi_real8,iwest,2,     &
                          wrkkuo1(:,jxp+1,k),iy,mpi_real8,ieast,2, &
                          mycomm,mpi_status_ignore,ierr)
        call mpi_sendrecv(wrkkuo2(:,jxp,k),iy,mpi_real8,ieast,1, &
                          wrkkuo2(:,0,k),iy,mpi_real8,iwest,1,   &
                          mycomm,mpi_status_ignore,ierr)
        call mpi_sendrecv(wrkkuo2(:,1,k),iy,mpi_real8,iwest,2,     &
                          wrkkuo2(:,jxp+1,k),iy,mpi_real8,ieast,2, &
                          mycomm,mpi_status_ignore,ierr)
      end do
      call htdiff(wrkkuo1,wrkkuo2,dxsq,akht1,jbegin,jendx,2,iym2)
    end if
#ifndef BAND
      ! diagnostic on total evaporation
      if (debug_level > 2) call conqeva
#endif
!
!   Call medium resolution PBL
!
    if ( ibltyp == 2 .or. ibltyp == 99 ) then
      ! Call the Grenier and Bretherton (2001) / Bretherton (2004) TCM
      call uwtcm
      call get_data_from_tcm(uwstateb,uwten,aten,atm1,atm2,.true.)
    end if
    if ( ibltyp == 1 .or. ibltyp == 99 ) then
      ! exchange internal ghost points
      call mpi_sendrecv(sfsta%uvdrag(1,jxp),iy,mpi_real8,ieast,1,  &
                        sfsta%uvdrag(1,0),  iy,mpi_real8,iwest,1,  &
                        mycomm,mpi_status_ignore,ierr)
      ! Call the Holtslag PBL
      call holtbl(jbegin,jendx,2,iym1)
    end if

    if ( ibltyp == 99 ) then
      call check_conserve_qt(holtten%qv,holtten%qc,uwten,uwstateb,kz)
      adf%diffq = adf%diffq + holtten%qv
      aten%qc = aten%qc + holtten%qc
    end if
!
    do j = jbegin , jendx
#ifndef BAND
      if ( myid /= nproc-1 .or. j /= jendx ) then
#endif
!       add ccm radiative transfer package-calculated heating rates to
!       temperature tendency
        do k = 1 , kz
          do i = 2 , iym2
            ! heating rate in deg/sec
            aten%t(i,k,j) = aten%t(i,k,j) + sps2%ps(j,i)*heatrt(j,i,k)
          end do
        end do
!
!       add horizontal diffusion and pbl tendencies for t and qv to aten%t
!       and aten%qv for calculating condensational term in subroutine
!       "condtq".
!
        do k = 1 , kz
          do i = 2 , iym2
            aten%t(i,k,j) = aten%t(i,k,j) + adf%difft(i,k,j)
          end do
        end do
!
        do k = 1 , kz
          do i = 2 , iym2
            aten%qv(i,k,j) = aten%qv(i,k,j) + adf%diffq(i,k,j)
          end do
        end do
      end if
    end do
!
!   compute the condensation and precipitation terms for explicit
!   moisture scheme:
!
    call condtq(jbegin,jendm,2,iym1,psc,qvcs)
!
!   subtract horizontal diffusion and pbl tendencies from aten%t and
!   aten%qv for appling the sponge boundary conditions on t and qv:
!
    do j = jbegin , jendx
#ifndef BAND
      if ( myid /= nproc-1 .or. j /= jendx ) then
#endif
        if ( iboudy == 4 ) then
          do k = 1 , kz
            do i = 2 , iym2
              aten%t(i,k,j) = aten%t(i,k,j) - adf%difft(i,k,j)
            end do
          end do
          call sponge(.false.,ispgx,wgtx,aten%t,j,kz,xtb)
          do k = 1 , kz
            do i = 2 , iym2
              aten%t(i,k,j) = aten%t(i,k,j) + adf%difft(i,k,j)
            end do
          end do
          do k = 1 , kz
            do i = 2 , iym2
              aten%qv(i,k,j) = aten%qv(i,k,j) - adf%diffq(i,k,j)
            end do
          end do
          call sponge(.false.,ispgx,wgtx,aten%qv,j,kz,xqb)
          do k = 1 , kz
            do i = 2 , iym2
              aten%qv(i,k,j) = aten%qv(i,k,j) + adf%diffq(i,k,j)
            end do
          end do
        end if
!
!       apply the nudging boundary conditions:
!
        if ( iboudy == 1 .or. iboudy == 5 ) then
          xtm1 = xbctime - dtsec
          if ( nbdytime == 0 .and. ktau /= 0 ) xtm1 = -dtsec
          call nudge(.false.,ispgx,fnudge,gnudge, &
                     xtm1,atm2%t,aten%t,j,kz,iboudy,xtb)
          call nudge(.false.,ispgx,fnudge,gnudge, &
                     xtm1,atm2%qv,aten%qv,j,kz,iboudy,xqb)
        end if
!
!       forecast t, qv, and qc at tau+1:
!
        do k = 1 , kz
          do i = 2 , iym2
            atmc%qv(i,k,j) = atm2%qv(i,k,j) + dt*aten%qv(i,k,j)
          end do
        end do
!
        do k = 1 , kz
          do i = 2 , iym2
            atmc%qc(i,k,j) = atm2%qc(i,k,j) + dt*aten%qc(i,k,j)
          end do
        end do
!
        do k = 1 , kz
          do i = 2 , iym2
            atmc%t(i,k,j) = atm2%t(i,k,j) + dt*aten%t(i,k,j)
          end do
        end do
!
!       forecast tracer chi at at tau+1:
        if ( ichem == 1 ) then
!
          do itr = 1 , ntr
            do k = 1 , kz
              do i = 2 , iym2
                chic(i,k,j,itr) = chib(i,k,j,itr) + dt*chiten(i,k,j,itr)
              end do
            end do
          end do
        end if
#ifndef BAND
      end if !end if (j /= jxm1),else test
#endif
    end do
!
#ifndef BAND
    do j = 1 , jendx
      if ( myid == 0 .and. j == 1 ) then
        if ( ipgf == 1 ) then
          do k = 1 , kz
            do i = 1 , iym1
              td(i,k,j) = atm1%t(i,k,j)*(d_one+ep1*(atmx%qv(i,k,j)))
              ttld(i,k,j) = td(i,k,j) - sps1%ps(j,i) * &
                        t00pg*((a(k)*sps1%ps(j,i)+ptop)/p00pg)**pgfaa1
            end do
          end do
        else if ( ipgf == 0 ) then
          do k = 1 , kz
            do i = 1 , iym1
              td(i,k,j) = atm1%t(i,k,j)*(d_one+ep1*(atmx%qv(i,k,j)))
            end do
          end do
        end if
!
      else if ( myid == nproc-1 .and. j == jendx ) then
!
!       set td and psd at j=jlx equal to ta and sps1%ps:
!
        if ( ipgf == 1 ) then
          do k = 1 , kz
            do i = 1 , iym1
              td(i,k,j) = atm1%t(i,k,j)*(d_one+ep1*(atmx%qv(i,k,j)))
              ttld(i,k,j) = td(i,k,j) - sps1%ps(j,i) * &
                     t00pg*((a(k)*sps1%ps(j,i)+ptop)/p00pg)**pgfaa1
            end do
          end do
        else if ( ipgf == 0 ) then
          do k = 1 , kz
            do i = 1 , iym1
              td(i,k,j) = atm1%t(i,k,j)*(d_one+ep1*(atmx%qv(i,k,j)))
            end do
          end do
        end if
!
!       compute weighted p*t (td) for use in ssi:
!
      end if
    end do
#endif

    do j = 1 , jendx
#ifndef BAND
      if ((.not.(myid == 0 .and. j == 1)) .and. &
         (.not.(myid == nproc-1 .and. j == jendx)) ) then
#endif
!
!     compute weighted p*t (td) for use in ssi:
!
      if ( ipgf == 1 ) then
        do k = 1 , kz
          do i = 2 , iym2
            tvc = atmc%t(i,k,j)*(d_one+ep1*(atmc%qv(i,k,j))/psc(i,j))
            tva = atm1%t(i,k,j)*(d_one+ep1*(atmx%qv(i,k,j)))
            tvb = atm2%t(i,k,j)*(d_one+ep1* &
                                 (atm2%qv(i,k,j))/sps2%ps(j,i))
            td(i,k,j) = alpha*(tvc+tvb) + beta*tva
            ttld(i,k,j) = td(i,k,j) - psd(i,j) * &
                      t00pg*((a(k)*psd(i,j)+ptop)/p00pg)**pgfaa1
          end do
        end do
        do k = 1 , kz
          td(1,k,j) = atm1%t(1,k,j)*(d_one+ep1*(atmx%qv(1,k,j)))
          ttld(1,k,j) = td(1,k,j) - sps1%ps(j,1) * &
                   t00pg*((a(k)*sps1%ps(j,1)+ptop)/p00pg)**pgfaa1
          td(iym1,k,j) = atm1%t(iym1,k,j)* &
                        (d_one+ep1*(atmx%qv(iym1,k,j)))
          ttld(iym1,k,j) = td(iym1,k,j) - sps1%ps(j,iym1) *           &
                   t00pg*((a(k)*sps1%ps(j,iym1)+ptop)/p00pg)**pgfaa1
        end do
!
      else if ( ipgf == 0 ) then
!
        do k = 1 , kz
          do i = 2 , iym2
            tvc = atmc%t(i,k,j)*(d_one+ep1*(atmc%qv(i,k,j))/psc(i,j))
            tva = atm1%t(i,k,j)*(d_one+ep1*(atmx%qv(i,k,j)))
            tvb = atm2%t(i,k,j)*(d_one+ep1* &
                 (atm2%qv(i,k,j))/sps2%ps(j,i))
            td(i,k,j) = alpha*(tvc+tvb) + beta*tva
          end do
        end do
        do k = 1 , kz
          td(1,k,j) = atm1%t(1,k,j)*(d_one+ep1*(atmx%qv(1,k,j)))
          td(iym1,k,j) = atm1%t(iym1,k,j)* &
                         (d_one+ep1*(atmx%qv(iym1,k,j)))
        end do
 
      end if
#ifndef BAND
      end if
#endif
    end do
!
!----------------------------------------------------------------------
!   compute the u and v tendencies:
!
    do j = jbegin , jendx
!
!     compute the diffusion terms:
!     put diffusion and pbl tendencies of u and v in difuu and difuv.
!
      do k = 1 , kz
        do i = 2 , iym1
          adf%difuu(i,k,j) = aten%u(i,k,j)
          adf%difuv(i,k,j) = aten%v(i,k,j)
        end do
      end do
    end do
!
    call diffu_d(jbegin,jendx,2,iym1,adf%difuu,atms%ubd3d, &
                 sps2%pdot,mddom%msfd,xkc,1)
    call diffu_d(jbegin,jendx,2,iym1,adf%difuv,atms%vbd3d, &
                 sps2%pdot,mddom%msfd,xkc,1)
!
!   compute the horizontal advection terms for u and v:
!
    do j = jbegin , jendx
      do k = 1 , kz
        do i = 2 , iym1
          aten%u(i,k,j) = d_zero
          aten%v(i,k,j) = d_zero
        end do
      end do
!
      call hadv(.true.,aten%u,atmx%u,j,j,3)
      call hadv(.true.,aten%v,atmx%v,j,j,3)
!
!     compute coriolis terms:
!
      do k = 1 , kz
        do i = 2 , iym1
          aten%u(i,k,j) = aten%u(i,k,j) + &
                       mddom%coriol(i,j)*atm1%v(i,k,j)/mddom%msfd(i,j)
          aten%v(i,k,j) = aten%v(i,k,j) - &
                       mddom%coriol(i,j)*atm1%u(i,k,j)/mddom%msfd(i,j)
        end do
      end do
    end do
!
    do j = jbegin , jendx
      jm1 = j-1
!
!     compute pressure gradient terms:
!
      if ( ipgf == 1 ) then
        do k = 1 , kz
          do i = 2 , iym1
            psasum = psd(i,j) + psd(i-1,j) + psd(i,jm1) + psd(i-1,jm1)
            sigpsa = psasum
            tv1 = atmx%t(i-1,k,jm1)*(d_one+ep1*(atmx%qv(i-1,k,jm1)))
            tv2 = atmx%t(i,k,jm1)*(d_one+ep1*(atmx%qv(i,k,jm1)))
            tv3 = atmx%t(i-1,k,j)*(d_one+ep1*(atmx%qv(i-1,k,j)))
            tv4 = atmx%t(i,k,j)*(d_one+ep1*(atmx%qv(i,k,j)))
            rtbar = tv1 + tv2 + tv3 + tv4 - d_four*t00pg*             &
                    ((a(k)*psasum*d_rfour+ptop)/p00pg)**pgfaa1
            rtbar = rgas*rtbar*sigpsa/16.0D0
            aten%u(i,k,j) = aten%u(i,k,j) - rtbar * &
                  (dlog(d_half*(psd(i,j)+psd(i-1,j))*a(k)+ptop) -     &
                   dlog(d_half*(psd(i,jm1)+psd(i-1,jm1))*a(k)+ptop))/ &
                   (dx*mddom%msfd(i,j))
            aten%v(i,k,j) = aten%v(i,k,j) - rtbar * &
                  (dlog(d_half*(psd(i,j)+psd(i,jm1))*a(k)+ptop) -     &
                   dlog(d_half*(psd(i-1,jm1)+psd(i-1,j))*a(k)+ptop))/ &
                   (dx*mddom%msfd(i,j))
          end do
        end do
      else if ( ipgf == 0 ) then
        do k = 1 , kz
          do i = 2 , iym1
            psasum = psd(i,j) + psd(i-1,j) + psd(i,jm1) + psd(i-1,jm1)
            sigpsa = psasum
            tv1 = atmx%t(i-1,k,jm1)*(d_one+ep1*(atmx%qv(i-1,k,jm1)))
            tv2 = atmx%t(i,k,jm1)*(d_one+ep1*(atmx%qv(i,k,jm1)))
            tv3 = atmx%t(i-1,k,j)*(d_one+ep1*(atmx%qv(i-1,k,j)))
            tv4 = atmx%t(i,k,j)*(d_one+ep1*(atmx%qv(i,k,j)))
            rtbar = rgas*(tv1+tv2+tv3+tv4)*sigpsa/16.0D0
            aten%u(i,k,j) = aten%u(i,k,j) - rtbar * &
                   (dlog(d_half*(psd(i,j)+psd(i-1,j))*a(k)+ptop) -    &
                    dlog(d_half*(psd(i,jm1)+psd(i-1,jm1))*a(k)+ptop))/&
                    (dx*mddom%msfd(i,j))
            aten%v(i,k,j) = aten%v(i,k,j) - rtbar *                   &
                   (dlog(d_half*(psd(i,j)+psd(i,jm1))*a(k)+ptop) -    &
                    dlog(d_half*(psd(i-1,jm1)+psd(i-1,j))*a(k)+ptop))/&
                    (dx*mddom%msfd(i,j))
          end do
        end do
      else   ! ipgf if block
      end if
    end do
!
    do j = 1 , jendx
!
!     compute geopotential height at half-k levels, cross points:
!
      if ( ipgf == 1 ) then
 
        do i = 1 , iym1
          tv = (ttld(i,kz,j)/psd(i,j))/(d_one+atmx%qc(i,kz,j)/ &
                                       (d_one+atmx%qv(i,kz,j)))
          phi(i,kz,j) = mddom%ht(i,j) + &
                   rgas*t00pg/pgfaa1*((psd(i,j)+ptop)/p00pg)**pgfaa1
          phi(i,kz,j) = phi(i,kz,j) - rgas * &
                  tv*dlog((a(kz)+ptop/psd(i,j))/(d_one+ptop/psd(i,j)))
        end do
 
        do k = 1 , kzm1
          lev = kz - k
          do i = 1 , iym1
            tvavg = ((ttld(i,lev,j)*dsigma(lev)+ttld(i,lev+1,j)* &
                    dsigma(lev+1))/(psd(i,j)*(dsigma(lev)+       &
                    dsigma(lev+1))))/(d_one+atmx%qc(i,lev,j)/    &
                    (d_one+atmx%qv(i,lev,j)))
            phi(i,lev,j) = phi(i,lev+1,j) - rgas *    &
                   tvavg*dlog((a(lev)+ptop/psd(i,j))/ &
                             (a(lev+1)+ptop/psd(i,j)))
          end do
        end do
 
      else if ( ipgf == 0 ) then
 
        do i = 1 , iym1
          tv = (td(i,kz,j)/psd(i,j))/(d_one+atmx%qc(i,kz,j)/  &
               (d_one+atmx%qv(i,kz,j)))
          phi(i,kz,j) = mddom%ht(i,j) - rgas * &
               tv*dlog((a(kz)+ptop/psd(i,j))/(d_one+ptop/psd(i,j)))
        end do
 
        do k = 1 , kzm1
          lev = kz - k
          do i = 1 , iym1
            tvavg = ((td(i,lev,j)*dsigma(lev)+td(i,lev+1,j)*   &
                    dsigma(lev+1))/(psd(i,j)*(dsigma(lev)+     &
                    dsigma(lev+1))))/(d_one+atmx%qc(i,lev,j)/  &
                    (d_one+atmx%qv(i,lev,j)))
            phi(i,lev,j) = phi(i,lev+1,j) - rgas *    &
                   tvavg*dlog((a(lev)+ptop/psd(i,j))  &
                           /(a(lev+1)+ptop/psd(i,j)))
          end do
        end do
 
      else   ! ipgf if block
      end if
    end do
    call mpi_sendrecv(phi(1,1,jxp),iy*kz,mpi_real8,ieast,1, &
                      phi(1,1,0),  iy*kz,mpi_real8,iwest,1, &
                      mycomm,mpi_status_ignore,ierr)
!
    do j = jbegin , jendx
      jm1 = j-1
!
!     compute the geopotential gradient terms:
!
      do k = 1 , kz
        do i = 2 , iym1
          aten%u(i,k,j) = aten%u(i,k,j) -                              &
               (psd(i-1,jm1)+psd(i,jm1)+psd(i-1,j)+psd(i,j)) *         &
               (phi(i,k,j)+phi(i-1,k,j)-phi(i,k,jm1)-phi(i-1,k,jm1)) / &
               (dx8*mddom%msfd(i,j))
          aten%v(i,k,j) = aten%v(i,k,j) -                              &
               (psd(i-1,jm1)+psd(i,jm1)+psd(i-1,j)+psd(i,j)) *         &
               (phi(i,k,j)+phi(i,k,jm1)-phi(i-1,k,j)-phi(i-1,k,jm1)) / &
               (dx8*mddom%msfd(i,j))
        end do
      end do
    end do
!
    do j = jbegin , jendx
!
!     compute the vertical advection terms:
!
      call vadv(aten%u,atm1%u,j,j,4)
      call vadv(aten%v,atm1%v,j,j,4)
!
!     apply the sponge boundary condition on u and v:
!
      if ( iboudy == 4 ) then
        call sponge(.true.,ispgd,wgtd,aten%u,j,kz,xub)
        call sponge(.true.,ispgd,wgtd,aten%v,j,kz,xvb)
      end if
!
!     apply the nudging boundary conditions:
!
      if ( iboudy == 1 .or. iboudy == 5 ) then
        call nudge(.true.,ispgd,fnudge,gnudge,xtm1, &
                   atm2%u,aten%u,j,kz,iboudy,xub)
        call nudge(.true.,ispgd,fnudge,gnudge,xtm1, &
                   atm2%v,aten%v,j,kz,iboudy,xvb)
      end if
!
!     add the diffusion and pbl tendencies to aten%u and aten%v:
!
      do k = 1 , kz
        do i = 2 , iym1
          aten%u(i,k,j) = aten%u(i,k,j) + adf%difuu(i,k,j)
          aten%v(i,k,j) = aten%v(i,k,j) + adf%difuv(i,k,j)
        end do
      end do
!
!     forecast p*u and p*v at tau+1:
!
      do k = 1 , kz
        do i = 2 , iym1
          atmc%u(i,k,j) = atm2%u(i,k,j) + dt*aten%u(i,k,j)
          atmc%v(i,k,j) = atm2%v(i,k,j) + dt*aten%v(i,k,j)
        end do
      end do
!
      ! Couple TKE to ps for use in vertical advection
      if ( ibltyp == 2 .or. ibltyp == 99 ) then
        do i = 1 , iym1
          do k = 1 , kz
            uwstatea%tkeps(j,i,k) = atm1%tke(i,k,j)*sps1%ps(j,i)
            uwstatea%advtke(i,k,j) = d_zero
          end do
          uwstatea%tkeps(j,i,kz+1) = atm1%tke(i,kz+1,j)*sps1%ps(j,i)
        end do
      end if
    end do ! end of j loop.

    if ( ibltyp == 2 .or. ibltyp == 99 ) then
      do j = jbegin , jendx
        ! Don't work with TKE on boundary grid-cells
        if ( (.not.(myid == 0 .and. j == 1)) .and. &
             (.not.(myid == (nproc-1) .and. j == jx)) ) then
          ! Calculate the horizontal advective tendency for TKE
          call hadvtke(uwstatea,atm1,twt,dx4,j)
          ! Calculate the vertical advective tendency for TKE
          call vadvtke(uwstatea,qdot,j,2)
          ! Calculate the horizontal, diffusive tendency for TKE
        end if
      end do
      call diffu_x(jbegin,jendx,2,iym2,uwstatea%advtke, &
                   atms%tkeb3d,sps2%ps,xkcf,kzp1)
    end if
!
!**********************************************************************
!
!---------------------------------------------------------------------
!
!   store the xxa variables in xxb and xxc in xxa:
!   perform time smoothing operations.
!
    do j = jbegin , jendx
      do k = 1 , kz
        do i = 2 , iym1
          atm2%u(i,k,j) = omuhf*atm1%u(i,k,j)/mddom%msfd(i,j)  &
                      + gnuhf*(atm2%u(i,k,j)+atmc%u(i,k,j))
          atm2%v(i,k,j) = omuhf*atm1%v(i,k,j)/mddom%msfd(i,j)  &
                      + gnuhf*(atm2%v(i,k,j)+atmc%v(i,k,j))
          atm1%u(i,k,j) = atmc%u(i,k,j)
          atm1%v(i,k,j) = atmc%v(i,k,j)
          ! TAO: Once the full loop above is completed, update the TKE
          ! tendency if the UW PBL is running.  NOTE!!! Do not try to
          ! combine these loops with the above loop Advection MUST be
          ! done in a loop separate from the updates.  (I lost 3 days
          ! of working to disocover that this is a problem because I
          ! thought it would be clever to combine loops--TAO)
          if ( ibltyp == 2 .or. ibltyp == 99 ) then
            ! Add the advective tendency to the TKE tendency calculated
            ! by the UW TKE routine
             aten%tke(i,k,j) = aten%tke(i,k,j) + &
                               uwstatea%advtke(i,k,j)/sps1%ps(j,i)
             ! Do a filtered time integration
             atmc%tke(i,k,j) = max(TKEMIN,atm2%tke(i,k,j) + &
                               dttke*aten%tke(i,k,j))
             atm2%tke(i,k,j) = max(TKEMIN,omuhf*atm1%tke(i,k,j) + &
                               gnuhf*(atm2%tke(i,k,j) + atmc%tke(i,k,j)))
             atm1%tke(i,k,j) = atmc%tke(i,k,j)
          end if ! TKE tendency update
        end do
      end do
    end do
!
    do j = jbegin , jendm
      do k = 1 , kz
        do i = 2 , iym2
          atm2%t(i,k,j) = omuhf*atm1%t(i,k,j) + &
                          gnuhf*(atm2%t(i,k,j)+atmc%t(i,k,j))
          atm1%t(i,k,j) = atmc%t(i,k,j)
          qvas = atmc%qv(i,k,j)
          if ( qvas < dlowval ) qvas = d_zero
          qvbs = omuhf*atm1%qv(i,k,j) + &
                 gnuhf*(atm2%qv(i,k,j)+atmc%qv(i,k,j))
          if ( qvbs < dlowval ) qvbs = d_zero
          atm2%qv(i,k,j) = qvbs
          atm1%qv(i,k,j) = qvas
          qcas = atmc%qc(i,k,j)
          if ( qcas < dlowval ) qcas = d_zero
          qcbs = omu*atm1%qc(i,k,j) + &
                 gnu*(atm2%qc(i,k,j)+atmc%qc(i,k,j))
          if ( qcbs < dlowval ) qcbs = d_zero
          atm2%qc(i,k,j) = qcbs
          atm1%qc(i,k,j) = qcas
        end do
        if ( ichem == 1 ) then
          do itr = 1 , ntr
            do i = 2 , iym2
              chias = chic(i,k,j,itr)
              if ( chias < dlowval ) chias = d_zero
              chibs = omu*chia(i,k,j,itr)                             &
                      + gnu*(chib(i,k,j,itr)+chic(i,k,j,itr))
              if ( chibs < dlowval ) chibs = d_zero
              chib(i,k,j,itr) = chibs
              chia(i,k,j,itr) = chias
            end do
          end do
        end if
      end do
      do i = 2 , iym2
        sps2%ps(j,i) = omuhf*sps1%ps(j,i) + &
                       gnuhf*(sps2%ps(j,i)+psc(i,j))
        sps1%ps(j,i) = psc(i,j)
      end do
    end do
!
!----------------------------------------------------------------------
!
!   increment elapsed forecast time:
!
    ktau = ktau + 1
    xbctime = xbctime + dtsec
    nbdytime = nbdytime + ntsec
    idatex = idatex + intmdl

    if ( mod(ktau,khour) == 0 ) then
      call split_idate(idatex,xyear,xmonth,xday,xhour)
    end if
    if ( mod(ktau,kbdy) == 0 ) then
      nbdytime = 0
      xbctime = d_zero
    end if

    if ( iexec == 2 ) then
      dt = dt2
      dtcum = dt2
      dtche = dt2
      dtpbl = dt2
      rdtpbl = d_one/dt2
      dttke = dt2
      iexec = 3
    end if
!
!     compute the amounts advected through the lateral boundaries:
!     *** note *** we must calculate the amounts advected through
!     the lateral boundaries before updating the values
!     at boundary slices.
!
#ifndef BAND
    if (debug_level > 2) then
      call conadv
      if ( ichem == 1 ) call tracdiag(xkc)
    end if
#endif
! 
!   fill up the boundary values for xxb and xxa variables:
!
    call bdyval(xbctime,iexec)
!
!   compute the nonconvective precipitation:
!
!   do cumulus transport of tracers
    if ( ichem == 1 .and. ichcumtra == 1 ) call cumtran
! 
!   trace the mass conservation of dry air and water substance:
!
#ifndef BAND
    if (debug_level > 2) call conmas
#endif
!
!   budgets for tracers
!
    if ( ichem == 1 ) then
      call tracbud
#ifndef BAND
      if (debug_level > 2) call contrac
#endif
    end if
!
!   print out noise parameter:
!
    if ( ktau > 1 ) then
      ptnbar = ptntot/dble(iptn)
      pt2bar = pt2tot/dble(iptn)
      icons_mpi = 0
      call mpi_allreduce(total_precip_points,icons_mpi,1,mpi_integer, &
                         mpi_sum,mycomm,ierr)
      ! Added a check for nan... The following inequality is wanted.
      if ((ptnbar /= ptnbar) .or. &
         ((ptnbar > d_zero) .eqv. (ptnbar <= d_zero))) then
        if ( myid == 0 ) then
          write (*,*) 'WHUUUUBBBASAAAGASDDWD!!!!!!!!!!!!!!!!'
          write (*,*) 'No more atmosphere here....'
          write (*,*) 'CFL violation detected, so model STOP'
          write (*,*) '#####################################'
          write (*,*) '#            DECREASE DT !!!!       #'
          write (*,*) '#####################################'
          call fatal(__FILE__,__LINE__,'CFL VIOLATION')
        end if
      end if
      if ( myid == 0 ) then
        if ( mod(ktau,irep) == 0 ) then
          appdat = tochar(idatex)
          write(6,99001) appdat , ktau , ptnbar , pt2bar , icons_mpi
        end if
      end if

99001 format (a23,', ktau = ',i10, ' :  1st, 2nd time deriv of ps = ',2E12.5, &
             ',  no. of points w/convection = ',i7)
    end if
!
!----------------------------------------------------------------------
!
!   recalculate solar declination angle if forecast time larger than
!   24 hours:
!
    if ( nbdytime == 0 ) then
      if (myid == 0) then
        write (6,*) 'Recalculate solar declination angle at ',toint10(idatex)
      end if
      call solar1
    end if
!
    call time_end(subroutine_name,idindx)
  end subroutine tend
!
end module mod_tendency
