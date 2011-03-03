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

      module mod_diagnosis
!
! Diagnostic printout subroutines
!
#ifndef BAND

      use mod_constants
      use mod_dynparam
      use mod_runparams
      use mod_main
      use mod_mainchem
      use mod_date
      use mod_trachem
#ifdef MPP1
      use mod_mppio
#endif
!
      private
!
      public :: allocate_mod_diagnosis , initdiag
      public :: restdiag , restchemdiag
      public :: savediag , savechemdiag
      public :: conadv , conmas , conqeva
      public :: tracdiag , contrac
#ifdef MPP1
      public :: mpidiag
#endif
!
      real(8) :: tdadv , tdini , tqadv , tqeva , tqini , tqrai
!
      real(8) , allocatable , dimension(:) :: tchiad , tchie , tchitb
      real(8) , allocatable , dimension(:,:) :: tremcvc , tremdrd ,     &
               & tremlsc , trxsaq1 , trxsaq2 , trxsg , ttrace
      contains

      subroutine allocate_mod_diagnosis
      implicit none
        if ( ichem == 1 ) then
          allocate(tchie(ntr))
          allocate(tchiad(ntr))
          allocate(tchitb(ntr))
          allocate(tremcvc(ntr,2))
          allocate(tremdrd(ntr,2))
          allocate(tremlsc(ntr,2))
          allocate(trxsaq1(ntr,2))
          allocate(trxsaq2(ntr,2))
          allocate(trxsg(ntr,2))
          allocate(ttrace(ntr,2))
          tchiad = d_zero
          tchie = d_zero
          tchitb = d_zero
          tremcvc = d_zero
          tremdrd = d_zero
          tremlsc = d_zero
          trxsaq1 = d_zero
          trxsaq2 = d_zero
          trxsg = d_zero
          ttrace = d_zero
        end if
      end subroutine allocate_mod_diagnosis
!
      subroutine initdiag
#ifdef MPP1
#ifndef IBM
        use mpi
#else 
        include 'mpif.h'
#endif 
#endif
        implicit none
#ifdef MPP1
        integer :: ierr
#endif
        real(8) :: tvmass , tcmass , tttmp
        integer :: i , j , k

        tvmass = d_zero
        tcmass = d_zero
        tttmp = d_zero
        tdini = d_zero
        tdadv = d_zero
        tqini = d_zero
        tqadv = d_zero
        tqeva = d_zero
        tqrai = d_zero
        if ( ichem == 1 ) then
          ttrace(:,:) = d_zero
          tchie(:) = d_zero
          tchiad(:) = d_zero
          tchitb(:) = d_zero
        end if

#ifdef MPP1
!=======================================================================
!
!-----dry air (unit = kg):
!
        call mpi_gather(sps1%ps, iy*jxp,mpi_real8, &
                      & psa_io,  iy*jxp,mpi_real8, &
                      & 0,mpi_comm_world,ierr)
        if ( myid == 0 ) then
          do k = 1 , kz
            tttmp = d_zero
            do j = 1 , jxm1
              do i = 1 , iym1
                tttmp = tttmp + psa_io(i,j)
              end do
            end do
            tdini = tdini + tttmp*dsigma(k)
          end do
          tdini = tdini*dx*dx*d_1000*rgti
        end if

        call mpi_bcast(tdini,1,mpi_real8,0,mpi_comm_world,ierr)
!
!-----water substance (unit = kg):
!
        call mpi_gather(atm1%qv,   iy*kz*jxp,mpi_real8, &
                      & atm1_io%qv,iy*kz*jxp,mpi_real8, &
                      & 0,mpi_comm_world,ierr)
        if ( myid == 0 ) then
          do k = 1 , kz
            tttmp = d_zero
            do j = 1 , jxm1
              do i = 1 , iym1
                tttmp = tttmp + atm1_io%qv(i,k,j)
              end do
            end do
            tvmass = tvmass + tttmp*dsigma(k)
          end do
          tvmass = tvmass*dx*dx*d_1000*rgti
        end if

        call mpi_bcast(tvmass,1,mpi_real8,0,mpi_comm_world,ierr)
!
        call mpi_gather(atm1%qc,   iy*kz*jxp,mpi_real8,              &
                      & atm1_io%qc,iy*kz*jxp,mpi_real8,              &
                      & 0,mpi_comm_world,ierr)
        if ( myid == 0 ) then
          do k = 1 , kz
            tttmp = d_zero
            do j = 1 , jxm1
              do i = 1 , iym1
                tttmp = tttmp + atm1_io%qc(i,k,j)
              end do
            end do
            tcmass = tcmass + tttmp*dsigma(k)
          end do
          tcmass = tcmass*dx*dx*d_1000*rgti
        end if

        call mpi_bcast(tcmass,1,mpi_real8,0,mpi_comm_world,ierr)

        tqini = tvmass + tcmass
!
!=======================================================================
!
        if ( myid == 0 ) write(6,99003) tdini , tqini
#else
!=======================================================================
!
!-----dry air (unit = kg):
!
        do k = 1 , kz
          tttmp = d_zero
          do j = 1 , jxm1
            do i = 1 , iym1
              tttmp = tttmp + sps1%ps(i,j)
            end do
          end do
          tdini = tdini + tttmp*dsigma(k)
        end do
        tdini = tdini*dx*dx*d_1000*rgti
!
!-----water substance (unit = kg):
!
        do k = 1 , kz
          tttmp = d_zero
          do j = 1 , jxm1
            do i = 1 , iym1
              tttmp = tttmp + atm1%qv(i,k,j)
            end do
          end do
          tvmass = tvmass + tttmp*dsigma(k)
        end do
        tvmass = tvmass*dx*dx*d_1000*rgti
!
        do k = 1 , kz
          tttmp = d_zero
          do j = 1 , jxm1
            do i = 1 , iym1
              tttmp = tttmp + atm1%qc(i,k,j)
            end do
          end do
          tcmass = tcmass + tttmp*dsigma(k)
        end do
        tcmass = tcmass*dx*dx*d_1000*rgti
        tqini = tvmass + tcmass
!=======================================================================
        write(6,99003) tdini , tqini
#endif
99003 format (' *** initial total air = ',e12.5,' kg, total water = ',  &
            & e12.5,' kg in large domain.')
!
      end subroutine initdiag
!
!
#ifdef MPP1
      subroutine mpidiag
#ifndef IBM
        use mpi
#else 
        include 'mpif.h'
#endif 
        implicit none
        integer :: ierr
        call mpi_bcast(tdini,1,mpi_real8,0,mpi_comm_world,ierr)
        call mpi_bcast(tdadv,1,mpi_real8,0,mpi_comm_world,ierr)
        call mpi_bcast(tqini,1,mpi_real8,0,mpi_comm_world,ierr)
        call mpi_bcast(tqadv,1,mpi_real8,0,mpi_comm_world,ierr)
        call mpi_bcast(tqeva,1,mpi_real8,0,mpi_comm_world,ierr)
        call mpi_bcast(tqrai,1,mpi_real8,0,mpi_comm_world,ierr)
        if ( ichem == 1 ) then
          call mpi_bcast(tchiad,ntr,mpi_real8,0,mpi_comm_world,ierr)
          call mpi_bcast(tchitb,ntr,mpi_real8,0,mpi_comm_world,ierr)
          call mpi_bcast(tchie,ntr,mpi_real8,0,mpi_comm_world,ierr)
        end if
      end subroutine mpidiag
#endif
!
!
      subroutine restdiag(iunit)
        implicit none
        integer , intent(in) :: iunit
        read (iunit) tdini , tdadv , tqini , tqadv , tqeva , tqrai
      end subroutine restdiag
!
      subroutine restchemdiag(iunit)
        implicit none
        integer , intent(in) :: iunit
        if ( ichem == 1 ) then
          read (iunit) tchiad
          read (iunit) tchitb
        end if
      end subroutine restchemdiag
!
      subroutine savediag(iunit)
        implicit none
        integer , intent(in) :: iunit
        write (iunit) tdini , tdadv , tqini , tqadv , tqeva , tqrai
      end subroutine savediag
!
      subroutine savechemdiag(iunit)
        implicit none
        integer , intent(in) :: iunit
        if ( ichem == 1 ) then
          write (iunit) tchiad
          write (iunit) tchitb
          write (iunit) tchie
        end if
      end subroutine savechemdiag
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!     this subroutine computes the amounts of dry air and water       c
!     substance advected through the lateral boundaries.              c
!                                                                     c
!     ---the unit from advection is converted to "kg".                c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      subroutine conadv
!
#ifdef MPP1
#ifndef IBM
      use mpi
#else 
      include 'mpif.h'
#endif 
#endif
!
      implicit none
!
#ifdef MPP1
      integer :: ierr
      real(8) , dimension(jxp) :: psa01 , psailx
      real(8) , dimension(jx) :: psa01_g , psailx_g
      real(8) , dimension(kz,jxp) :: qca01 , qcailx , qva01 , qvailx ,  &
                                   & va01 , vaix
      real(8) , dimension(kz,jx) :: qca01_g , qcailx_g , qva01_g ,      &
                                   & qvailx_g , va01_g , vaix_g
#endif
      real(8) , dimension(iym1,kz) :: worka , workb
      integer :: i ,  j , k
!
!----------------------------------------------------------------------
!-----advection of dry air through the lateral boundaries:
!
!.....advection through east-west boundaries:
!
#ifdef MPP1
      if ( myid == nproc-1 ) then
        do k = 1 , kz
          do i = 1 , iym1
            worka(i,k) = (atm1%u(i+1,k,jendl)+atm1%u(i,k,jendl)) /      &
                       &  (mddom%msfx(i,jendx)*mddom%msfx(i,jendx))
          end do
        end do
      end if
      if ( myid == 0 ) then
        do k = 1 , kz
          do i = 1 , iym1
            workb(i,k) = (atm1%u(i+1,k,1)+atm1%u(i,k,1)) / &
                          (mddom%msfx(i,1)*mddom%msfx(i,1))
          end do
        end do
      end if
      call mpi_bcast(worka,iym1*kz,mpi_real8,nproc-1,                   &
                   & mpi_comm_world,ierr)
      call mpi_bcast(workb,iym1*kz,mpi_real8,0,                         &
                   & mpi_comm_world,ierr)
#else
      do k = 1 , kz
        do i = 1 , iym1
          worka(i,k) = (atm1%u(i+1,k,jx)+atm1%u(i,k,jx)) / &
                        (mddom%msfx(i,jxm1)*mddom%msfx(i,jxm1))
          workb(i,k) = (atm1%u(i+1,k,1)+atm1%u(i,k,1)) / &
                        (mddom%msfx(i,1)*mddom%msfx(i,1))
        end do
      end do
#endif
      do k = 1 , kz
        do i = 1 , iym1
          tdadv = tdadv - dtmin*3.0D4*dsigma(k)                          &
                & *dx*(worka(i,k)-workb(i,k))*rgti
        end do
      end do
!
!.....advection through north-south boundaries:
!
#ifdef MPP1
      do j = 1 , jendl
        do k = 1 , kz
          vaix(k,j) = atm1%v(iy,k,j)
          va01(k,j) = atm1%v(1,k,j)
        end do
      end do
      call mpi_gather(vaix,  kz*jxp,mpi_real8,                     &
                    & vaix_g,kz*jxp,mpi_real8,                     &
                    & 0,mpi_comm_world,ierr)
      call mpi_gather(va01,  kz*jxp,mpi_real8,                     &
                    & va01_g,kz*jxp,mpi_real8,                     &
                    & 0,mpi_comm_world,ierr)
      if ( myid == 0 ) then
        do k = 1 , kz
          do j = 1 , jxm1
            tdadv = tdadv - dtmin*3.0D4*dsigma(k)                   &
                  & *dx*((vaix_g(k,j+1)+vaix_g(k,j))               &
                  & /(mddom_io%msfx(iym1,j)*mddom_io%msfx(iym1,j)) &
                  & -(va01_g(k,j+1)+va01_g(k,j))                   &
                  & /(mddom_io%msfx(1,j)*mddom_io%msfx(1,j)))*rgti
          end do
        end do
      end if
      call mpi_bcast(tdadv,1,mpi_real8,0,mpi_comm_world,ierr)
#else
      do k = 1 , kz
        do j = 1 , jxm1
          tdadv = tdadv - dtmin*3.0D4*dsigma(k)                &
                  *dx*((atm1%v(iy,k,j+1)+atm1%v(iy,k,j)) /    &
                  (mddom%msfx(iym1,j)*mddom%msfx(iym1,j)) -   &
                  (atm1%v(1,k,j+1)+atm1%v(1,k,j)) /           &
                  (mddom%msfx(1,j)*mddom%msfx(1,j)))*rgti
        end do
      end do
#endif
!
!----------------------------------------------------------------------
!-----advection of water vapor through the lateral boundaries:
!
!.....advection through east-west boundaries:
!
#ifdef MPP1
      if ( myid == nproc-1 ) then
        do k = 1 , kz
          do i = 1 , iym1
            worka(i,k) = (atm1%u(i+1,k,jendl)+atm1%u(i,k,jendl))    &
                       & *(atm1%qv(i,k,jendx)/sps1%ps(i,jendx))     &
                       & /(mddom%msfx(i,jendx)*mddom%msfx(i,jendx))
          end do
        end do
      end if
      if ( myid == 0 ) then
        do k = 1 , kz
          do i = 1 , iym1
            workb(i,k) = (atm1%u(i+1,k,1)+atm1%u(i,k,1)) * &
                          (atm1%qv(i,k,1)/sps1%ps(i,1)) / &
                          (mddom%msfx(i,1)*mddom%msfx(i,1))
          end do
        end do
      end if
      call mpi_bcast(worka,iym1*kz,mpi_real8,nproc-1,                   &
                   & mpi_comm_world,ierr)
      call mpi_bcast(workb,iym1*kz,mpi_real8,0,                         &
                   & mpi_comm_world,ierr)
#else
      do k = 1 , kz
        do i = 1 , iym1
          worka(i,k) = (atm1%u(i+1,k,jx)+atm1%u(i,k,jx))     &
                     & *(atm1%qv(i,k,jxm1)/sps1%ps(i,jxm1))  &
                     & /(mddom%msfx(i,jxm1)*mddom%msfx(i,jxm1))
          workb(i,k) = (atm1%u(i+1,k,1)+atm1%u(i,k,1))*(atm1%qv(i,k,1)/ &
                       sps1%ps(i,1))/(mddom%msfx(i,1)*mddom%msfx(i,1))
        end do
      end do
#endif
      do k = 1 , kz
        do i = 1 , iym1
          tqadv = tqadv - dtmin*3.0D4*dsigma(k)                          &
                & *dx*(worka(i,k)-workb(i,k))*rgti
        end do
      end do
!
!....advection through north-south boundaries:
!
#ifdef MPP1
      do j = 1 , jendl
        do k = 1 , kz
          qvailx(k,j) = atm1%qv(iym1,k,j)
          qva01(k,j) = atm1%qv(1,k,j)
        end do
        psailx(j) = sps1%ps(iym1,j)
        psa01(j) = sps1%ps(1,j)
      end do
      call mpi_gather(qvailx,  kz*jxp,mpi_real8,                   &
                    & qvailx_g,kz*jxp,mpi_real8,                   &
                    & 0,mpi_comm_world,ierr)
      call mpi_gather(qva01,  kz*jxp,mpi_real8,                    &
                    & qva01_g,kz*jxp,mpi_real8,                    &
                    & 0,mpi_comm_world,ierr)
      call mpi_gather(psailx,  jxp,mpi_real8,                        &
                    & psailx_g,jxp,mpi_real8,0,mpi_comm_world,ierr)
      call mpi_gather(psa01,  jxp,mpi_real8,                         &
                      psa01_g,jxp,mpi_real8,0,mpi_comm_world,ierr)
      if ( myid == 0 ) then
        do k = 1 , kz
          do j = 1 , jxm1
            tqadv = tqadv - dtmin*3.0D4*dsigma(k)                        &
                  & *dx*((vaix_g(k,j+1)+vaix_g(k,j))                    &
                  & *(qvailx_g(k,j)/psailx_g(j))                        &
                  & /(mddom_io%msfx(iym1,j)*mddom_io%msfx(iym1,j))      &
                  & -(va01_g(k,j+1)+va01_g(k,j))                        &
                  & *(qva01_g(k,j)/psa01_g(j))                          &
                  & /(mddom_io%msfx(1,j)*mddom_io%msfx(1,j)))*rgti
          end do
        end do
      end if
      call mpi_bcast(tqadv,1,mpi_real8,0,mpi_comm_world,ierr)
#else
      do k = 1 , kz
        do j = 1 , jxm1
          tqadv = tqadv - dtmin*3.0D4*dsigma(k)*            &
               dx*((atm1%v(iy,k,j+1)+atm1%v(iy,k,j))*      &
               (atm1%qv(iym1,k,j)/sps1%ps(iym1,j)) /       &
               (mddom%msfx(iym1,j)*mddom%msfx(iym1,j)) -   &
               (atm1%v(1,k,j+1)+atm1%v(1,k,j))*(atm1%qv(1,k,j) /   &
               sps1%ps(1,j))/(mddom%msfx(1,j)*mddom%msfx(1,j)))*rgti
        end do
      end do
#endif
!
!-----advection of cloud water and rainwater through lateral boundaries:
!
!.....advection through east-west boundaries:
!
#ifdef MPP1
      if ( myid == nproc-1 ) then
        do k = 1 , kz
          do i = 1 , iym1
            worka(i,k) = (atm1%u(i+1,k,jendl)+atm1%u(i,k,jendl))        &
                       & *(atm1%qc(i,k,jendx)/sps1%ps(i,jendx))         &
                       & /(mddom%msfx(i,jendx)*mddom%msfx(i,jendx))
          end do
        end do
      end if
      if ( myid == 0 ) then
        do k = 1 , kz
          do i = 1 , iym1
            workb(i,k) = (atm1%u(i+1,k,1)+atm1%u(i,k,1))* &
                         (atm1%qc(i,k,1)/sps1%ps(i,1)) /  &
                         (mddom%msfx(i,1)*mddom%msfx(i,1))
          end do
        end do
      end if
      call mpi_bcast(worka,iym1*kz,mpi_real8,nproc-1,                   &
                   & mpi_comm_world,ierr)
      call mpi_bcast(workb,iym1*kz,mpi_real8,0,                         &
                   & mpi_comm_world,ierr)
#else
      do k = 1 , kz
        do i = 1 , iym1
          worka(i,k) = (atm1%u(i+1,k,jx)+atm1%u(i,k,jx))      &
                     & *(atm1%qc(i,k,jxm1)/sps1%ps(i,jxm1))   &
                     & /(mddom%msfx(i,jxm1)*mddom%msfx(i,jxm1))
          workb(i,k) = (atm1%u(i+1,k,1)+atm1%u(i,k,1))*(atm1%qc(i,k,1) /&
                     sps1%ps(i,1))/(mddom%msfx(i,1)*mddom%msfx(i,1))
        end do
      end do
#endif
      do k = 1 , kz
        do i = 1 , iym1
          tqadv = tqadv - dtmin*3.0D4*dsigma(k)                          &
                & *dx*(worka(i,k)-workb(i,k))*rgti
        end do
      end do
!
!....advection through north-south boundaries:
!
#ifdef MPP1
      do j = 1 , jendl
        do k = 1 , kz
          qcailx(k,j) = atm1%qc(iym1,k,j)
          qca01(k,j) = atm1%qc(1,k,j)
        end do
      end do
      call mpi_gather(qcailx,  kz*jxp,mpi_real8,                   &
                    & qcailx_g,kz*jxp,mpi_real8,                   &
                    & 0,mpi_comm_world,ierr)
      call mpi_gather(qca01,  kz*jxp,mpi_real8,                    &
                    & qca01_g,kz*jxp,mpi_real8,                    &
                    & 0,mpi_comm_world,ierr)
      if ( myid == 0 ) then
        do k = 1 , kz
          do j = 1 , jxm1
            tqadv = tqadv - dtmin*3.0D4*dsigma(k)                        &
                  & *dx*((vaix_g(k,j+1)+vaix_g(k,j))                    &
                  & *(qcailx_g(k,j)/psailx_g(j))                        &
                  & /(mddom_io%msfx(iym1,j)*mddom_io%msfx(iym1,j))      &
                  & -(va01_g(k,j+1)+va01_g(k,j))                        &
                  & *(qca01_g(k,j)/psa01_g(j))                          &
                  & /(mddom_io%msfx(1,j)*mddom_io%msfx(1,j)))*rgti
          end do
        end do
      end if
      call mpi_bcast(tqadv,1,mpi_real8,0,mpi_comm_world,ierr)
#else
      do k = 1 , kz
        do j = 1 , jxm1
          tqadv = tqadv - dtmin*3.0D4*dsigma(k)*            &
                dx*((atm1%v(iy,k,j+1)+atm1%v(iy,k,j))*     &
                (atm1%qc(iym1,k,j)/sps1%ps(iym1,j)) /      &
                (mddom%msfx(iym1,j)*mddom%msfx(iym1,j))-   &
                (atm1%v(1,k,j+1)+atm1%v(1,k,j))*(atm1%qc(1,k,j) /  &
                sps1%ps(1,j))/(mddom%msfx(1,j)*mddom%msfx(1,j)))*rgti
        end do
      end do
#endif
!
      end subroutine conadv
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!     this subroutine computes the total dry air and water substance  c
!     within the domain and compares with the initial values.         c
!                                                                     c
!     ---the unit used in all the calculation is "kg".                c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      subroutine conmas
!
#ifdef MPP1
#ifndef IBM
      use mpi
#else 
      include 'mpif.h'
#endif 
#endif
      implicit none
!
! Local variables
!
      real(8) :: error1 , error2 , tcmass , tcrai , tdrym , tncrai ,    &
               & tqmass , tttmp , tvmass , xh
      integer :: i , j , k
#ifdef MPP1
      integer :: ierr
#endif
!
!----------------------------------------------------------------------
!
      error1 = d_zero
      error2 = d_zero
!
!-----compute the total dry air and water substance in the model at
!     this time:
!
!=======================================================================
!
!-----dry air (unit = kg):
!
      tdrym = d_zero
#ifdef MPP1
      call mpi_gather(sps1%ps,   iy*jxp,mpi_real8,                     &
                    & psa_io,iy*jxp,mpi_real8,                     &
                    & 0,mpi_comm_world,ierr)
      if ( myid == 0 ) then
        do k = 1 , kz
          tttmp = d_zero
          do j = 1 , jxm1
            do i = 1 , iym1
              tttmp = tttmp + psa_io(i,j)
            end do
          end do
          tdrym = tdrym + tttmp*dsigma(k)
        end do
        tdrym = tdrym*dx*dx*d_1000*rgti
      end if
      call mpi_bcast(tdrym,1,mpi_real8,0,mpi_comm_world,ierr)
#else
      do k = 1 , kz
        tttmp = d_zero
        do j = 1 , jxm1
          do i = 1 , iym1
            tttmp = tttmp + sps1%ps(i,j)
          end do
        end do
        tdrym = tdrym + tttmp*dsigma(k)
      end do
      tdrym = tdrym*dx*dx*d_1000*rgti
#endif
!
!-----water substance (unit = kg):
!
      tvmass = d_zero
#ifdef MPP1
      call mpi_gather(atm1%qv,   iy*kz*jxp,mpi_real8,                &
                    & atm1_io%qv,iy*kz*jxp,mpi_real8,                &
                    & 0,mpi_comm_world,ierr)
      if ( myid == 0 ) then
        do k = 1 , kz
          tttmp = d_zero
          do j = 1 , jxm1
            do i = 1 , iym1
              tttmp = tttmp + atm1_io%qv(i,k,j)
            end do
          end do
          tvmass = tvmass + tttmp*dsigma(k)
        end do
        tvmass = tvmass*dx*dx*d_1000*rgti
      end if
      call mpi_bcast(tvmass,1,mpi_real8,0,mpi_comm_world,               &
                   & ierr)
#else
      do k = 1 , kz
        tttmp = d_zero
        do j = 1 , jxm1
          do i = 1 , iym1
            tttmp = tttmp + atm1%qv(i,k,j)
          end do
        end do
        tvmass = tvmass + tttmp*dsigma(k)
      end do
      tvmass = tvmass*dx*dx*d_1000*rgti
#endif
!
      tcmass = d_zero

#ifdef MPP1
      call mpi_gather(atm1%qc,   iy*kz*jxp,mpi_real8,                &
                    & atm1_io%qc,iy*kz*jxp,mpi_real8,                &
                    & 0,mpi_comm_world,ierr)
      if ( myid == 0 ) then
        do k = 1 , kz
          tttmp = d_zero
          do j = 1 , jxm1
            do i = 1 , iym1
              tttmp = tttmp + atm1_io%qc(i,k,j)
            end do
          end do
          tcmass = tcmass + tttmp*dsigma(k)
        end do
        tcmass = tcmass*dx*dx*d_1000*rgti
      end if
      call mpi_bcast(tcmass,1,mpi_real8,0,mpi_comm_world,               &
                   & ierr)
#else
      do k = 1 , kz
        tttmp = d_zero
        do j = 1 , jxm1
          do i = 1 , iym1
            tttmp = tttmp + atm1%qc(i,k,j)
          end do
        end do
        tcmass = tcmass + tttmp*dsigma(k)
      end do
      tcmass = tcmass*dx*dx*d_1000*rgti
#endif

      tqmass = tvmass + tcmass

!=======================================================================
!
!-----conservation of dry air:
!
      tdrym = tdrym - tdadv
      error1 = (tdrym-tdini)/tdini*d_100
!
!-----conservation of water substance:
!
!-----total raifall at this time:
!
#ifdef MPP1
      call mpi_gather(sfsta%rainc,   iy*jxp,mpi_real8,                   &
                    & rainc_io,iy*jxp,mpi_real8,                   &
                    & 0,mpi_comm_world,ierr)
      call mpi_gather(sfsta%rainnc,   iy*jxp,mpi_real8,                  &
                    & rainnc_io,iy*jxp,mpi_real8,                  &
                    & 0,mpi_comm_world,ierr)
      if ( myid == 0 ) then
        tcrai = d_zero
        tncrai = d_zero
        do j = 1 , jxm1
          do i = 1 , iym1
            tcrai = tcrai + rainc_io(i,j)*dxsq
            tncrai = tncrai + rainnc_io(i,j)*dxsq
          end do
        end do
        tqrai = tcrai + tncrai
      end if
      call mpi_bcast(tqrai,1,mpi_real8,0,mpi_comm_world,ierr)
#else
      tcrai = d_zero
      tncrai = d_zero
      do j = 1 , jxm1
        do i = 1 , iym1
          tcrai = tcrai + sfsta%rainc(i,j)*dxsq
          tncrai = tncrai + sfsta%rainnc(i,j)*dxsq
        end do
      end do
      tqrai = tcrai + tncrai
#endif
!
      tqmass = tqmass + tqrai - tqeva - tqadv
      error2 = (tqmass-tqini)/tqini*d_100
!
!-----print out the information:
!
#ifdef MPP1
      if ( myid == 0 ) then
#endif
        if ( debug_level > 3 .and. mod(ntime,ndbgfrq) == 0 ) then
          xh = xtime/minpd
          write(6,*)  '***** day = ' , ldatez + xh , ' *****'
          write(6,99001) tdrym , error1
          write(6,99002) tdadv
          write(6,99003) tqmass , error2
          write(6,99004) tvmass
          write(6,99005) tcmass
          write(6,99006) tqadv
          write(6,99007) tcrai
          write(6,99008) tncrai
          write(6,99009) tqeva
        end if
#ifdef MPP1
      end if
#endif

99001 format ('   total air =',e12.5,' kg, error = ',e12.5)
99002 format ('   horizontal advection = ',e12.5,' kg.')
99003 format ('   total water =',e12.5,' kg, error = ',e12.5)
99004 format ('   qv                      = ',e12.5,' kg.')
99005 format ('   qc                      = ',e12.5,' kg.')
99006 format ('   horizontal advection    = ',e12.5,' kg.')
99007 format ('   convective railfall     = ',e12.5,' kg.')
99008 format ('   nonconvective rainfall  = ',e12.5,' kg.')
99009 format ('   evaporation from ground = ',e12.5,' kg.')
 
      end subroutine conmas
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
! Tracer diagnostic
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      subroutine tracdiag(xkc)
#ifdef MPP1
#ifndef IBM
      use mpi
#else 
      include 'mpif.h'
#endif 
#endif
      implicit none
!
#ifdef MPP1
      real(8) , dimension(iy,kz,jxp) :: xkc
#else
      real(8) , dimension(iy,kz,jx) :: xkc
#endif
      intent (in) xkc
!
#ifdef MPP1
      integer :: ierr
      real(8) , dimension(kz,ntr,jxp) :: chia01 , chia02 , chiaill ,    &
           & chiaill1
      real(8) , dimension(kz,ntr,jx) :: chia01_g , chia02_g ,           &
           & chiaill1_g , chiaill_g
      real(8) , dimension(jxp) :: psa01 , psa02 , psaill , psaill1
      real(8) , dimension(jx) :: psa01_g , psa02_g , psaill1_g ,        &
                                & psaill_g
      real(8) , dimension(kz,jxp) :: va02 , vaill , xkc02 , xkcill1
      real(8) , dimension(kz,jx) :: va02_g , vaill_g , xkc02_g ,        &
                                   & xkcill1_g
      real(8) :: chid1 , chid2
#endif
      real(8) :: fact1 , fact2 , fx1 , fx2 , uavg1 , uavg2 , vavg1 ,    &
           &  vavg2
      integer :: i , j , k , n
      real(8) , dimension(iym1,kz,ntr) :: worka , workb
!
!ccccccccccccccccccccccccccccccccccccccccccccccc
 
!-------------------------
!     1  ADVECTION budgets
!-----------------------
 
!-----advection of tracer through lateral boundaries
 
!
!.....advection through east-west boundaries:
!
!     the 'relaxed' upstream scheme
      fact1 = 0.6D0
      fact2 = d_one - fact1
 
!     inflow/outflow
#ifdef MPP1
      do n = 1 , ntr
        do k = 1 , kz
          do i = 2 , iym2
            if ( myid == nproc-1 ) then
              uavg2 = d_half*(atm1%u(i+1,k,jendx)+atm1%u(i,k,jendx))
              if ( uavg2 < d_zero ) then
                worka(i,k,n) =  &
                    -uavg2*(fact1*chia(i,k,jendx,n)/sps1%ps(i,jendx)/ &
                    (mddom%msfx(i,jendx)*mddom%msfx(i,jendx))+      &
                    fact2*chia(i,k,jendm,n)/sps1%ps(i,jendm)/         &
                    (mddom%msfx(i,jendm)*mddom%msfx(i,jendm)))
              else
                worka(i,k,n) = &
                    -uavg2*(fact1*chia(i,k,jendm,n)/sps1%ps(i,jendm)/ &
                    (mddom%msfx(i,jendm)*mddom%msfx(i,jendm))+      &
                    fact2*chia(i,k,jendx,n)/sps1%ps(i,jendx) /        &
                    (mddom%msfx(i,jendx)*mddom%msfx(i,jendx)))
              end if
            end if
            if ( myid == 0 ) then
              uavg1 = d_half*(atm1%u(i+1,k,1+1)+atm1%u(i,k,1+1))
              if ( uavg1 > d_zero ) then
                workb(i,k,n) = &
                    -uavg1*(fact1*chia(i,k,1,n)/sps1%ps(i,1)/ &
                    (mddom%msfx(i,1)*mddom%msfx(i,1)) +     &
                    fact2*chia(i,k,1+1,n)/sps1%ps(i,1+1) /    &
                    (mddom%msfx(i,1+1)*mddom%msfx(i,1+1)))
              else
                workb(i,k,n) = & 
                    -uavg1*(fact1*chia(i,k,1+1,n)/sps1%ps(i,1+1) / &
                    (mddom%msfx(i,1+1)*mddom%msfx(i,1+1)) +      &
                    fact2*chia(i,k,1,n)/sps1%ps(i,1) /             &
                    (mddom%msfx(i,1)*mddom%msfx(i,1)))
              end if
            end if
          end do
        end do
      end do
      call mpi_bcast(worka,iym1*kz*ntr,mpi_real8,nproc-1,               &
                   & mpi_comm_world,ierr)
#else
      do n = 1 , ntr
        do k = 1 , kz
          do i = 2 , iym2
            uavg2 = d_half*(atm1%u(i+1,k,jxm1)+atm1%u(i,k,jxm1))
            if ( uavg2 < d_zero ) then
              worka(i,k,n) =  &
                  -uavg2*(fact1*chia(i,k,jxm1,n)/sps1%ps(i,jxm1) / &
                  (mddom%msfx(i,jxm1)*mddom%msfx(i,jxm1)) +      &
                  fact2*chia(i,k,jxm2,n)/sps1%ps(i,jxm2) /         &
                  (mddom%msfx(i,jxm2)*mddom%msfx(i,jxm2)))
            else
              worka(i,k,n) = &
                  -uavg2*(fact1*chia(i,k,jxm2,n)/sps1%ps(i,jxm2) / & 
                  (mddom%msfx(i,jxm2)*mddom%msfx(i,jxm2)) +      &
                  fact2*chia(i,k,jxm1,n)/sps1%ps(i,jxm1) /         &
                  (mddom%msfx(i,jxm1)*mddom%msfx(i,jxm1)))
            end if
 
            uavg1 = d_half*(atm1%u(i+1,k,1+1)+atm1%u(i,k,1+1))
            if ( uavg1 > d_zero ) then
              workb(i,k,n) = &
                  -uavg1*(fact1*chia(i,k,1,n)/sps1%ps(i,1) / &
                  (mddom%msfx(i,1)*mddom%msfx(i,1)) +      &
                  fact2*chia(i,k,1+1,n)/sps1%ps(i,1+1) /     &
                  (mddom%msfx(i,1+1)*mddom%msfx(i,1+1)))
            else
              workb(i,k,n) = &
                  -uavg1*(fact1*chia(i,k,1+1,n)/sps1%ps(i,1+1) / &
                  (mddom%msfx(i,1+1)*mddom%msfx(i,1+1)) +      &
                  fact2*chia(i,k,1,n)/sps1%ps(i,1) /             &
                  (mddom%msfx(i,1)*mddom%msfx(i,1)))
            end if
          end do
        end do
      end do
#endif 

#ifdef MPP1
      do j = 1 , jendl
        do k = 1 , kz
          vaill(k,j) = atm1%v(iym1,k,j)
          va02(k,j) = atm1%v(2,k,j)
          xkcill1(k,j) = xkc(iym2,k,j)
          xkc02(k,j) = xkc(2,k,j)
          do n = 1 , ntr
            chiaill(k,n,j) = chia(iym1,k,j,n)
            chiaill1(k,n,j) = chia(iym2,k,j,n)
            chia01(k,n,j) = chia(1,k,j,n)
            chia02(k,n,j) = chia(2,k,j,n)
          end do
        end do
        psaill(j) = sps1%ps(iym1,j)
        psaill1(j) = sps1%ps(iym2,j)
        psa01(j) = sps1%ps(1,j)
        psa02(j) = sps1%ps(2,j)
      end do
      call mpi_gather(vaill,  kz*jxp,mpi_real8,                    &
                    & vaill_g,kz*jxp,mpi_real8,                    &
                    & 0,mpi_comm_world,ierr)
      call mpi_gather(va02,  kz*jxp,mpi_real8,                     &
                    & va02_g,kz*jxp,mpi_real8,                     &
                    & 0,mpi_comm_world,ierr)
      call mpi_gather(xkcill1,  kz*jxp,mpi_real8,                  &
                    & xkcill1_g,kz*jxp,mpi_real8,                  &
                    & 0,mpi_comm_world,ierr)
      call mpi_gather(xkc02,  kz*jxp,mpi_real8,                    &
                    & xkc02_g,kz*jxp,mpi_real8,                    &
                    & 0,mpi_comm_world,ierr)
      call mpi_gather(chiaill,  kz*ntr*jxp,mpi_real8,            &
                    & chiaill_g,kz*ntr*jxp,mpi_real8,            &
                    & 0,mpi_comm_world,ierr)
      call mpi_gather(chiaill1,  kz*ntr*jxp,mpi_real8,           &
                    & chiaill1_g,kz*ntr*jxp,mpi_real8,           &
                    & 0,mpi_comm_world,ierr)
      call mpi_gather(chia01,  kz*ntr*jxp,mpi_real8,             &
                    & chia01_g,kz*ntr*jxp,mpi_real8,             &
                    & 0,mpi_comm_world,ierr)
      call mpi_gather(chia02,  kz*ntr*jxp,mpi_real8,             &
                    & chia02_g,kz*ntr*jxp,mpi_real8,             &
                    & 0,mpi_comm_world,ierr)
      call mpi_gather(psaill,  jxp,mpi_real8,                        &
                    & psaill_g,jxp,mpi_real8,0,mpi_comm_world,ierr)
      call mpi_gather(psaill1,  jxp,mpi_real8,                       &
                    & psaill1_g,jxp,mpi_real8,0,mpi_comm_world,ierr)
      call mpi_gather(psa01,  jxp,mpi_real8,                         &
                    & psa01_g,jxp,mpi_real8,0,mpi_comm_world,ierr)
      call mpi_gather(psa02,  jxp,mpi_real8,                         &
                    & psa02_g,jxp,mpi_real8,0,mpi_comm_world,ierr)
      if ( myid == 0 ) then
        do n = 1 , ntr
          do k = 1 , kz
            do i = 2 , iym2
              tchiad(n) = tchiad(n) + dtmin*6.0D4*dsigma(k)     &
                        & *dx*(worka(i,k,n)-workb(i,k,n))*rgti
            end do
          end do
!.....
!.....advection through north-south boundaries:
!
          do k = 1 , kz
            do j = 2 , jxm2
!hy           inflow/outflow
              vavg2 = d_half*(vaill_g(k,j+1)+vaill_g(k,j))
              if ( vavg2 < d_zero ) then
                fx2 = -vavg2*(fact1*chiaill_g(k,n,j)/psaill_g(j)        &
                    & /(mddom_io%msfx(iym1,j)*mddom_io%msfx(iym1,j))    &
                    & +fact2*chiaill1_g(k,n,j)/psaill1_g(j)             &
                    & /(mddom_io%msfx(iym2,j)*mddom_io%msfx(iym2,j)))
              else
                fx2 = -vavg2*(fact1*chiaill1_g(k,n,j)/psaill1_g(j)      &
                    & /(mddom_io%msfx(iym2,j)*mddom_io%msfx(iym2,j))    &
                    & +fact2*chiaill_g(k,n,j)/psaill_g(j)               &
                    & /(mddom_io%msfx(iym1,j)*mddom_io%msfx(iym1,j)))
              end if
 
              vavg1 = d_half*(va02_g(k,j+1)+va02_g(k,j))
              if ( vavg1 > d_zero ) then
                fx1 = -vavg1*(fact1*chia01_g(k,n,j)/psa01_g(j)          &
                    & /(mddom_io%msfx(1,j)*mddom_io%msfx(1,j))+         &
                    & fact2*chia02_g(k,n,j)/psa02_g(j)                  &
                    & /(mddom_io%msfx(2,j)*mddom_io%msfx(2,j)))
              else
                fx1 = -vavg1*(fact1*chia02_g(k,n,j)/psa02_g(j)     &
                    & /(mddom_io%msfx(2,j)*mddom_io%msfx(2,j))+    &
                    & fact2*chia01_g(k,n,j)/psa01_g(j)             &
                    & /(mddom_io%msfx(1,j)*mddom_io%msfx(1,j)))
              end if
              tchiad(n) = tchiad(n) + dtmin*6.0D4*dsigma(k)*dx*(fx2-fx1) &
                        & *rgti
            end do
          end do
        end do
      end if
      call mpi_bcast(tchiad,ntr,mpi_real8,0,mpi_comm_world,ierr)
#else
      do n = 1 , ntr
        do k = 1 , kz
          do i = 2 , iym2
            tchiad(n) = tchiad(n) + dtmin*6.0D4*dsigma(k)                &
                      & *dx*(worka(i,k,n)-workb(i,k,n))*rgti
          end do
        end do
      end do
!.....
!.....advection through north-south boundaries:
!
      do n = 1 , ntr
        do k = 1 , kz
          do j = 2 , jxm2
!hy         inflow/outflow
            vavg2 = d_half*(atm1%v(iym1,k,j+1)+atm1%v(iym1,k,j))
            if ( vavg2 < d_zero ) then
              fx2 = -vavg2*(fact1*chia(iym1,k,j,n)/sps1%ps(iym1,j) / &
                      (mddom%msfx(iym1,j)*mddom%msfx(iym1,j)) +    &
                       fact2*chia(iym2,k,j,n)/sps1%ps(iym2,j) /      &
                       (mddom%msfx(iym2,j)*mddom%msfx(iym2,j)))
            else
              fx2 = -vavg2*(fact1*chia(iym2,k,j,n)/sps1%ps(iym2,j) &
                  & /(mddom%msfx(iym2,j)*mddom%msfx(iym2,j))     &
                  & +fact2*chia(iym1,k,j,n)/sps1%ps(iym1,j)        &
                  & /(mddom%msfx(iym1,j)*mddom%msfx(iym1,j)))
            end if
 
            vavg1 = d_half*(atm1%v(1+1,k,j+1)+atm1%v(1+1,k,j))
            if ( vavg1 > d_zero ) then
              fx1 = -vavg1*(fact1*chia(1,k,j,n)/sps1%ps(1,j) / &
                      (mddom%msfx(1,j)*mddom%msfx(1,j)) +    &
                       fact2*chia(1+1,k,j,n)/sps1%ps(1+1,j) /  &
                       (mddom%msfx(1+1,j)*mddom%msfx(1+1,j)))
            else
              fx1 = -vavg1*(fact1*chia(1+1,k,j,n)/sps1%ps(1+1,j) / &
                      (mddom%msfx(1+1,j)*mddom%msfx(1+1,j)) +    &
                       fact2*chia(1,k,j,n)/sps1%ps(1,j) /          &
                       (mddom%msfx(1,j)*mddom%msfx(1,j)))
            end if
 
            tchiad(n) = tchiad(n) + dtmin*6.0D4*dsigma(k)*dx*            &
                  & (fx2-fx1)*rgti
 
          end do
        end do
      end do
#endif
 
!
!..... diffusion through east-west boundaries:
!
#ifdef MPP1
      do n = 1 , ntr
        do k = 1 , kz
          do i = 2 , iym2
            if ( myid == nproc-1 )  &
              worka(i,k,n) = xkc(i,k,jendm)*sps1%ps(i,jendm) * &
                  (chia(i,k,jendm,n)/sps1%ps(i,jendm)-         &
                   chia(i,k,jendx,n)/sps1%ps(i,jendx))
            if ( myid == 0 ) &
              workb(i,k,n) = xkc(i,k,2)*sps1%ps(i,2) *  &
                  (chia(i,k,2,n)/sps1%ps(i,2) -         &
                   chia(i,k,1,n)/sps1%ps(i,1))
          end do
        end do
      end do
      call mpi_bcast(worka,iym1*kz*ntr,mpi_real8,nproc-1,               &
                   & mpi_comm_world,ierr)
#else
      do n = 1 , ntr
        do k = 1 , kz
          do i = 2 , iym2
            worka(i,k,n) = xkc(i,k,jxm2)*sps1%ps(i,jxm2)       &
                         & *(chia(i,k,jxm2,n)/sps1%ps(i,jxm2)  &
                         & -chia(i,k,jxm1,n)/sps1%ps(i,jxm1))
            workb(i,k,n) = xkc(i,k,2)*sps1%ps(i,2)             &
                         & *(chia(i,k,2,n)/sps1%ps(i,2)-chia(i,k,1,n) &
                         & /sps1%ps(i,1))
          end do
        end do
      end do
#endif

#ifdef MPP1
      if ( myid == 0 ) then
        do n = 1 , ntr
          do k = 1 , kz
            do i = 2 , iym2
              tchitb(n) = tchitb(n) - dtmin*6.0D4*dsigma(k)              &
                        & *(workb(i,k,n)+worka(i,k,n))*rgti
            end do
          end do
 
!.....  diffusion through north-south boundaries:
 
          do k = 1 , kz
            do j = 2 , jxm2
              chid1 = xkcill1_g(k,j)*psaill1_g(j)                       &
                    & *(chiaill1_g(k,n,j)/psaill1_g(j)-chiaill_g(k,n,j) &
                    & /psaill_g(j))
              chid2 = xkc02_g(k,j)*psa02_g(j)                           &
                    & *(chia02_g(k,n,j)/psa02_g(j)-chia01_g(k,n,j)      &
                    & /psa01_g(j))
              tchitb(n) = tchitb(n) - dtmin*6.0D4*dsigma(k)*(chid2+chid1)&
                        & *rgti
            end do
          end do
        end do
      end if
      call mpi_bcast(tchitb,ntr,mpi_real8,0,mpi_comm_world,ierr)
#else
      do n = 1 , ntr
        do k = 1 , kz
          do i = 2 , iym2
            tchitb(n) = tchitb(n) - dtmin*6.0D4*dsigma(k)                &
                      & *(workb(i,k,n)+worka(i,k,n))*rgti
          end do
        end do
      end do
#endif

      end subroutine tracdiag
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
! Tracer budget diagnostic
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      subroutine contrac
#ifdef MPP1
#ifndef IBM
      use mpi
#else 
      include 'mpif.h'
#endif 
#endif 
      implicit none
!
      integer :: i , itr , j , k
      real(8) :: tttmp
#ifdef MPP1
      integer :: ierr , l
#endif
!
!---- add tracers,  tremlsc, tremcvc, tremdrd are total amount mass
!     and chemical conversion term
!     removed by large scale, convective precipation and dry deposition
!     respectively
 
      do itr = 1 , ntr
        ttrace(itr,1) = d_zero
        tremlsc(itr,1) = d_zero
        tremcvc(itr,1) = d_zero
        trxsg(itr,1) = d_zero
        trxsaq1(itr,1) = d_zero
        trxsaq2(itr,1) = d_zero
        tremdrd(itr,1) = d_zero
      end do
!
#ifdef MPP1
      do itr = 1 , ntr
        call mpi_gather(chia(1,1,1,itr),   iy*kz*jxp,mpi_real8,         &
                      & chia_io(1,1,1,itr),iy*kz*jxp,mpi_real8,         &
                      & 0,mpi_comm_world,ierr)
        call mpi_gather(remlsc(1,1,1,itr),   iy*kz*jxp,mpi_real8,       &
                      & remlsc_io(1,1,1,itr),iy*kz*jxp,mpi_real8,       &
                      & 0,mpi_comm_world,ierr)
        call mpi_gather(remcvc(1,1,1,itr),   iy*kz*jxp,mpi_real8,       &
                      & remcvc_io(1,1,1,itr),iy*kz*jxp,mpi_real8,       &
                      & 0,mpi_comm_world,ierr)
        call mpi_gather(rxsg(1,1,1,itr),   iy*kz*jxp,mpi_real8,         &
                      & rxsg_io(1,1,1,itr),iy*kz*jxp,mpi_real8,         &
                      & 0,mpi_comm_world,ierr)
        call mpi_gather(rxsaq1(1,1,1,itr),   iy*kz*jxp,mpi_real8,       &
                      & rxsaq1_io(1,1,1,itr),iy*kz*jxp,mpi_real8,       &
                      & 0,mpi_comm_world,ierr)
        call mpi_gather(rxsaq2(1,1,1,itr),   iy*kz*jxp,mpi_real8,       &
                      & rxsaq2_io(1,1,1,itr),iy*kz*jxp,mpi_real8,       &
                      & 0,mpi_comm_world,ierr)
        call mpi_gather(remdrd(1,1,itr),   iy*jxp,mpi_real8,            &
                      & remdrd_io(1,1,itr),iy*jxp,mpi_real8,            &
                      & 0,mpi_comm_world,ierr)
        do l = 1 , mpy
          call mpi_gather(chemsrc(1,1,l,itr),   iy*jxp,mpi_real8,       &
                        & chemsrc_io(1,1,l,itr),iy*jxp,mpi_real8,       &
                        & 0,mpi_comm_world,ierr)
        end do
      end do
      if ( myid == 0 ) then
        do itr = 1 , ntr
          do k = 1 , kz
            tttmp = d_zero
            do j = 2 , jxm2
              do i = 2 , iym2
                tttmp = tttmp + chia_io(i,k,j,itr)
              end do
            end do
            ttrace(itr,1) = ttrace(itr,1) + tttmp*dsigma(k)
            tttmp = d_zero
            do j = 2 , jxm2
              do i = 2 , iym2
                tttmp = tttmp + remlsc_io(i,k,j,itr)
              end do
            end do
            tremlsc(itr,1) = tremlsc(itr,1) + tttmp*dsigma(k)
            tttmp = d_zero
            do j = 2 , jxm2
              do i = 2 , iym2
                tttmp = tttmp + remcvc_io(i,k,j,itr)
              end do
            end do
            tremcvc(itr,1) = tremcvc(itr,1) + tttmp*dsigma(k)
            tttmp = d_zero
            do j = 2 , jxm2
              do i = 2 , iym2
                tttmp = tttmp + rxsg_io(i,k,j,itr)
              end do
            end do
            trxsg(itr,1) = trxsg(itr,1) + tttmp*dsigma(k)
            tttmp = d_zero
            do j = 2 , jxm2
              do i = 2 , iym2
                tttmp = tttmp + rxsaq1_io(i,k,j,itr)
              end do
            end do
            trxsaq1(itr,1) = trxsaq1(itr,1) + tttmp*dsigma(k)
            tttmp = d_zero
            do j = 2 , jxm2
              do i = 2 , iym2
                tttmp = tttmp + rxsaq2_io(i,k,j,itr)
              end do
            end do
            trxsaq2(itr,1) = trxsaq2(itr,1) + tttmp*dsigma(k)
          end do
        end do
 
        do itr = 1 , ntr
          ttrace(itr,1) = ttrace(itr,1)*dx*dx
          tremlsc(itr,1) = tremlsc(itr,1)*dx*dx
          tremcvc(itr,1) = tremcvc(itr,1)*dx*dx
          trxsg(itr,1) = trxsg(itr,1)*dx*dx
          trxsaq1(itr,1) = trxsaq1(itr,1)*dx*dx
          trxsaq2(itr,1) = trxsaq2(itr,1)*dx*dx
        end do
 
        do itr = 1 , ntr
          tttmp = d_zero
          do j = 2 , jxm2
            do i = 2 , iym2
              tttmp = tttmp + remdrd_io(i,j,itr)
            end do
          end do
          tremdrd(itr,1) = tremdrd(itr,1) + tttmp*dx*dx*dsigma(kz)
 
!         emissions
          tttmp = d_zero
          do j = 2 , jxm2
            do i = 2 , iym2
              tttmp = tttmp + chemsrc_io(i,j,lmonth,itr)*dtmin*60.*dx*dx
            end do
          end do
          tchie(itr) = tchie(itr) + tttmp
        end do
      end if

      call mpi_bcast(ttrace(1,1),ntr,mpi_real8,0,mpi_comm_world,ierr)
      call mpi_bcast(tremlsc(1,1),ntr,mpi_real8,0,mpi_comm_world,ierr)
      call mpi_bcast(tremcvc(1,1),ntr,mpi_real8,0,mpi_comm_world,ierr)
      call mpi_bcast(trxsg(1,1),ntr,mpi_real8,0,mpi_comm_world,ierr)
      call mpi_bcast(trxsaq1(1,1),ntr,mpi_real8,0,mpi_comm_world,ierr)
      call mpi_bcast(trxsaq2(1,1),ntr,mpi_real8,0,mpi_comm_world,ierr)
      call mpi_bcast(tremdrd(1,1),ntr,mpi_real8,0,mpi_comm_world,ierr)
      call mpi_bcast(tchie(1),ntr,mpi_real8,0,mpi_comm_world,ierr)
#else
      do itr = 1 , ntr
        do k = 1 , kz
          tttmp = d_zero
          do j = 2 , jxm2
            do i = 2 , iym2
              tttmp = tttmp + chia(i,k,j,itr)
            end do
          end do
          ttrace(itr,1) = ttrace(itr,1) + tttmp*dsigma(k)
          tttmp = d_zero
          do j = 2 , jxm2
            do i = 2 , iym2
              tttmp = tttmp + remlsc(i,k,j,itr)
            end do
          end do
          tremlsc(itr,1) = tremlsc(itr,1) + tttmp*dsigma(k)
          tttmp = d_zero
          do j = 2 , jxm2
            do i = 2 , iym2
              tttmp = tttmp + remcvc(i,k,j,itr)
            end do
          end do
          tremcvc(itr,1) = tremcvc(itr,1) + tttmp*dsigma(k)
          tttmp = d_zero
          do j = 2 , jxm2
            do i = 2 , iym2
              tttmp = tttmp + rxsg(i,k,j,itr)
            end do
          end do
          trxsg(itr,1) = trxsg(itr,1) + tttmp*dsigma(k)
          tttmp = d_zero
          do j = 2 , jxm2
            do i = 2 , iym2
              tttmp = tttmp + rxsaq1(i,k,j,itr)
            end do
          end do
          trxsaq1(itr,1) = trxsaq1(itr,1) + tttmp*dsigma(k)
          tttmp = d_zero
          do j = 2 , jxm2
            do i = 2 , iym2
              tttmp = tttmp + rxsaq2(i,k,j,itr)
            end do
          end do
          trxsaq2(itr,1) = trxsaq2(itr,1) + tttmp*dsigma(k)
        end do
      end do

      do itr = 1 , ntr
        ttrace(itr,1) = ttrace(itr,1)*dx*dx
        tremlsc(itr,1) = tremlsc(itr,1)*dx*dx
        tremcvc(itr,1) = tremcvc(itr,1)*dx*dx
        trxsg(itr,1) = trxsg(itr,1)*dx*dx
        trxsaq1(itr,1) = trxsaq1(itr,1)*dx*dx
        trxsaq2(itr,1) = trxsaq2(itr,1)*dx*dx
      end do

      do itr = 1 , ntr
        tttmp = d_zero
        do j = 2 , jxm2
          do i = 2 , iym2
            tttmp = tttmp + remdrd(i,j,itr)
          end do
        end do
        tremdrd(itr,1) = tremdrd(itr,1) + tttmp*dx*dx*dsigma(kz)

!       emissions
        tttmp = d_zero
        do j = 2 , jxm2
          do i = 2 , iym2
            tttmp = tttmp + chemsrc(i,j,lmonth,itr)*dtmin*60.*dx*dx
          end do
        end do
        tchie(itr) = tchie(itr) + tttmp
      end do
#endif
 
      do itr = 1 , ntr
        ttrace(itr,1) = ttrace(itr,1)*d_1000*rgti
        tremlsc(itr,1) = tremlsc(itr,1)*d_1000*rgti
        tremcvc(itr,1) = tremcvc(itr,1)*d_1000*rgti
        tremdrd(itr,1) = tremdrd(itr,1)*d_1000*rgti
        trxsg(itr,1) = trxsg(itr,1)*d_1000*rgti
        trxsaq1(itr,1) = trxsaq1(itr,1)*d_1000*rgti
        trxsaq2(itr,1) = trxsaq2(itr,1)*d_1000*rgti
      end do

 
!-----print out the information:
#ifdef MPP1
      if ( myid == 0 ) then
#endif
 
        if ( debug_level > 3 .and. mod(ntime,ndbgfrq) == 0 ) then
 
!----     tracers
 
          write(6,*)  '************************************************'
          write(6,*)  ' Budgets for tracers (intergrated quantitites)'
          write(6,*)  ' day = ' , ldatez , ' *****'
          write(6,*)  '************************************************'
 
          do itr = 1 , ntr
 
!           errort= ttrace(itr,1)-
!           &     (tchie(itr) - (tchiad(itr)-tchitb(itr)+
!           &     tremdrd(itr,1)+ tremlsc(itr,1)+tremcvc(itr,1)))
 
            write(6,*)  '****************************'
            write(6,99001) itr , tchie(itr)
            write(6,99002) itr , tchiad(itr)
            write(6,99003) itr , tchitb(itr)
            write(6,99004) itr , tremdrd(itr,1)
            write(6,99005) itr , tremlsc(itr,1)
            write(6,99006) itr , tremcvc(itr,1)
            write(6,99007) itr , ttrace(itr,1)
!           write(6,99008) itr , errort
 
          end do
 
        end if
#ifdef MPP1
      end if
#endif
! 
99001 format ('total emission tracer',i1,':',e12.5,' kg')
99002 format ('total advected tracer',i1,':',e12.5,' kg')
99003 format ('total diffused tracer',i1,':',e12.5,' kg')
99004 format ('total dry deposition tracer',i1,':',e12.5,' kg')
99005 format ('total large scale wet dep.tracer',i1,':',e12.5,' kg')
99006 format ('total convective wet dep.tracer',i1,':',e12.5,' kg')
99007 format ('total mass (at ktau), tracer',i1,':',e12.5,' kg')
!99008 format('error trac',I1,':',e12.5,' % ')
      end subroutine contrac
!
! diagnostic on total evaporation
!
      subroutine conqeva
#ifdef MPP1
#ifndef IBM
      use mpi
#else 
      include 'mpif.h'
#endif 
#endif
      implicit none
#ifdef MPP1
      integer :: ierr
#endif
      integer :: i , j
#ifdef MPP1
      call mpi_gather(sfsta%qfx,iy*jxp,mpi_real8,  &
                    & qfx_io,   iy*jxp,mpi_real8,  &
                    & 0,mpi_comm_world,ierr)
      if ( myid == 0 ) then
        do j = 2 , jxm2
          do i = 2 , iym2
            tqeva = tqeva + qfx_io(i,j)*dx*dx*dtmin*minph
          end do
        end do
      end if
      call mpi_bcast(tqeva,1,mpi_real8,0,mpi_comm_world,ierr)
#else
      do j = 2 , jxm2
        do i = 2 , iym2
          tqeva = tqeva + sfsta%qfx(i,j)*dx*dx*dtmin*minph
        end do
      end do
#endif
      end subroutine conqeva
!
#endif
!
      end module mod_diagnosis
