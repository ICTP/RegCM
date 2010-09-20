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
        allocate(tchiad(ntr))
        allocate(tchie(ntr))
        allocate(tchitb(ntr))
        allocate(tremcvc(ntr,2))
        allocate(tremdrd(ntr,2))
        allocate(tremlsc(ntr,2))
        allocate(trxsaq1(ntr,2))
        allocate(trxsaq2(ntr,2))
        allocate(trxsg(ntr,2))
        allocate(ttrace(ntr,2))
        tchiad = 0.0D0
        tchie = 0.0D0
        tchitb = 0.0D0
        tremcvc = 0.0D0
        tremdrd = 0.0D0
        tremlsc = 0.0D0
        trxsaq1 = 0.0D0
        trxsaq2 = 0.0D0
        trxsg = 0.0D0
        ttrace = 0.0D0
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

        tvmass = 0.0
        tcmass = 0.0
        tttmp = 0.0
        tdini = 0.0
        tdadv = 0.0
        tqini = 0.0
        tqadv = 0.0
        tqeva = 0.0
        tqrai = 0.0
        ttrace = 0.0
        tchie = 0.0
        tchiad = 0.0
        tchitb = 0.0

#ifdef MPP1
!=======================================================================
!
!-----dry air (unit = kg):
!
        call mpi_gather(psa,   iy*jxp,mpi_real8, &
                      & psa_io,iy*jxp,mpi_real8, &
                      & 0,mpi_comm_world,ierr)
        if ( myid.eq.0 ) then
          do k = 1 , kz
            tttmp = 0.
            do j = 1 , jxm1
              do i = 1 , iym1
                tttmp = tttmp + psa_io(i,j)
              end do
            end do
            tdini = tdini + tttmp*dsigma(k)
          end do
          tdini = tdini*dx*dx*1000.*rgti
        end if
        call mpi_bcast(tdini,1,mpi_real8,0,mpi_comm_world,ierr)
!
!-----water substance (unit = kg):
!
        call mpi_gather(qva,   iy*kz*jxp,mpi_real8, &
                      & qva_io,iy*kz*jxp,mpi_real8, &
                      & 0,mpi_comm_world,ierr)
        if ( myid.eq.0 ) then
          do k = 1 , kz
            tttmp = 0.
            do j = 1 , jxm1
              do i = 1 , iym1
                tttmp = tttmp + qva_io(i,k,j)
              end do
            end do
            tvmass = tvmass + tttmp*dsigma(k)
          end do
          tvmass = tvmass*dx*dx*1000.*rgti
        end if
        call mpi_bcast(tvmass,1,mpi_real8,0,mpi_comm_world,ierr)
!
        call mpi_gather(qca,   iy*kz*jxp,mpi_real8,              &
                      & qca_io,iy*kz*jxp,mpi_real8,              &
                      & 0,mpi_comm_world,ierr)
        if ( myid.eq.0 ) then
          do k = 1 , kz
            tttmp = 0.
            do j = 1 , jxm1
              do i = 1 , iym1
                tttmp = tttmp + qca_io(i,k,j)
              end do
            end do
            tcmass = tcmass + tttmp*dsigma(k)
          end do
          tcmass = tcmass*dx*dx*1000.*rgti
        end if
        call mpi_bcast(tcmass,1,mpi_real8,0,mpi_comm_world,ierr)
        tqini = tvmass + tcmass
!=======================================================================
        if ( myid.eq.0 ) print 99003 , tdini , tqini
#else
!=======================================================================
!
!-----dry air (unit = kg):
!
        do k = 1 , kz
          tttmp = 0.
          do j = 1 , jxm1
            do i = 1 , iym1
              tttmp = tttmp + psa(i,j)
            end do
          end do
          tdini = tdini + tttmp*dsigma(k)
        end do
        tdini = tdini*dx*dx*1000.*rgti
!
!-----water substance (unit = kg):
!
        do k = 1 , kz
          tttmp = 0.
          do j = 1 , jxm1
            do i = 1 , iym1
              tttmp = tttmp + qva(i,k,j)
            end do
          end do
          tvmass = tvmass + tttmp*dsigma(k)
        end do
        tvmass = tvmass*dx*dx*1000.*rgti
!
        do k = 1 , kz
          tttmp = 0.
          do j = 1 , jxm1
            do i = 1 , iym1
              tttmp = tttmp + qca(i,k,j)
            end do
          end do
          tcmass = tcmass + tttmp*dsigma(k)
        end do
        tcmass = tcmass*dx*dx*1000.*rgti
        tqini = tvmass + tcmass
!=======================================================================
        print 99003 , tdini , tqini
#endif
99003 format (' *** initial total air = ',e12.5,' kg, total water = ',  &
            & e12.5,' kg in large domain.')
      end subroutine initdiag
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
        if ( ichem.eq.1 ) then
          call mpi_bcast(tchiad,ntr,mpi_real8,0,mpi_comm_world,ierr)
          call mpi_bcast(tchitb,ntr,mpi_real8,0,mpi_comm_world,ierr)
          call mpi_bcast(tchie,ntr,mpi_real8,0,mpi_comm_world,ierr)
        end if
      end subroutine mpidiag
#endif
!
      subroutine restdiag(iunit)
#ifdef MPP1
#ifndef IBM
        use mpi
#else 
        include 'mpif.h'
#endif 
#endif
        implicit none
        integer , intent(in) :: iunit
#ifdef MPP1
        integer :: ierr
#endif
        read (iunit) tdini , tdadv , tqini , tqadv , tqeva , tqrai
#ifdef MPP1
        call mpi_bcast(tdini,1,mpi_real8,0,mpi_comm_world,ierr)
        call mpi_bcast(tdadv,1,mpi_real8,0,mpi_comm_world,ierr)
        call mpi_bcast(tqini,1,mpi_real8,0,mpi_comm_world,ierr)
        call mpi_bcast(tqadv,1,mpi_real8,0,mpi_comm_world,ierr)
        call mpi_bcast(tqeva,1,mpi_real8,0,mpi_comm_world,ierr)
        call mpi_bcast(tqrai,1,mpi_real8,0,mpi_comm_world,ierr)
#endif
      end subroutine restdiag
!
      subroutine restchemdiag(iunit)
#ifdef MPP1
#ifndef IBM
        use mpi
#else 
        include 'mpif.h'
#endif 
#endif
        implicit none
        integer , intent(in) :: iunit
#ifdef MPP1
        integer :: ierr
#endif
        read (iunit) tchiad
        read (iunit) tchitb
#ifdef MPP1
        call mpi_bcast(tchiad,ntr,mpi_real8,0,mpi_comm_world,ierr)
        call mpi_bcast(tchitb,ntr,mpi_real8,0,mpi_comm_world,ierr)
        call mpi_bcast(tchie,ntr,mpi_real8,0,mpi_comm_world,ierr)
#endif
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
        write (iunit) tchiad
        write (iunit) tchitb
        write (iunit) tchie
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
      if ( myid.eq.nproc-1 ) then
        do k = 1 , kz
          do i = 1 , iym1
            worka(i,k) = (ua(i+1,k,jendl)+ua(i,k,jendl))                &
                       & /(msfx(i,jendx)*msfx(i,jendx))
          end do
        end do
      end if
      if ( myid.eq.0 ) then
        do k = 1 , kz
          do i = 1 , iym1
            workb(i,k) = (ua(i+1,k,1)+ua(i,k,1))/(msfx(i,1)*msfx(i,1))
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
          worka(i,k) = (ua(i+1,k,jx)+ua(i,k,jx))                        &
                     & /(msfx(i,jxm1)*msfx(i,jxm1))
          workb(i,k) = (ua(i+1,k,1)+ua(i,k,1))/(msfx(i,1)*msfx(i,1))
        end do
      end do
#endif
      do k = 1 , kz
        do i = 1 , iym1
          tdadv = tdadv - dtmin*3.E4*dsigma(k)                          &
                & *dx*(worka(i,k)-workb(i,k))*rgti
        end do
      end do
!
!.....advection through north-south boundaries:
!
#ifdef MPP1
      do j = 1 , jendl
        do k = 1 , kz
          vaix(k,j) = va(iy,k,j)
          va01(k,j) = va(1,k,j)
        end do
      end do
      call mpi_gather(vaix,  kz*jxp,mpi_real8,                     &
                    & vaix_g,kz*jxp,mpi_real8,                     &
                    & 0,mpi_comm_world,ierr)
      call mpi_gather(va01,  kz*jxp,mpi_real8,                     &
                    & va01_g,kz*jxp,mpi_real8,                     &
                    & 0,mpi_comm_world,ierr)
      if ( myid.eq.0 ) then
        do k = 1 , kz
          do j = 1 , jxm1
            tdadv = tdadv - dtmin*3.E4*dsigma(k)                        &
                  & *dx*((vaix_g(k,j+1)+vaix_g(k,j))                    &
                  & /(msfx_io(iym1,j)*msfx_io(iym1,j))                  &
                  & -(va01_g(k,j+1)+va01_g(k,j))                        &
                  & /(msfx_io(1,j)*msfx_io(1,j)))*rgti
          end do
        end do
      end if
      call mpi_bcast(tdadv,1,mpi_real8,0,mpi_comm_world,ierr)
#else
      do k = 1 , kz
        do j = 1 , jxm1
          tdadv = tdadv - dtmin*3.E4*dsigma(k)                          &
                & *dx*((va(iy,k,j+1)+va(iy,k,j))                        &
                & /(msfx(iym1,j)*msfx(iym1,j))-(va(1,k,j+1)+va(1,k,j))  &
                & /(msfx(1,j)*msfx(1,j)))*rgti
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
      if ( myid.eq.nproc-1 ) then
        do k = 1 , kz
          do i = 1 , iym1
            worka(i,k) = (ua(i+1,k,jendl)+ua(i,k,jendl))                &
                       & *(qva(i,k,jendx)/psa(i,jendx))                 &
                       & /(msfx(i,jendx)*msfx(i,jendx))
          end do
        end do
      end if
      if ( myid.eq.0 ) then
        do k = 1 , kz
          do i = 1 , iym1
            workb(i,k) = (ua(i+1,k,1)+ua(i,k,1))*(qva(i,k,1)/psa(i,1))  &
                       & /(msfx(i,1)*msfx(i,1))
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
          worka(i,k) = (ua(i+1,k,jx)+ua(i,k,jx))                        &
                     & *(qva(i,k,jxm1)/psa(i,jxm1))                     &
                     & /(msfx(i,jxm1)*msfx(i,jxm1))
          workb(i,k) = (ua(i+1,k,1)+ua(i,k,1))*(qva(i,k,1)/psa(i,1))    &
                     & /(msfx(i,1)*msfx(i,1))
        end do
      end do
#endif
      do k = 1 , kz
        do i = 1 , iym1
          tqadv = tqadv - dtmin*3.E4*dsigma(k)                          &
                & *dx*(worka(i,k)-workb(i,k))*rgti
        end do
      end do
!
!....advection through north-south boundaries:
!
#ifdef MPP1
      do j = 1 , jendl
        do k = 1 , kz
          qvailx(k,j) = qva(iym1,k,j)
          qva01(k,j) = qva(1,k,j)
        end do
        psailx(j) = psa(iym1,j)
        psa01(j) = psa(1,j)
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
      if ( myid.eq.0 ) then
        do k = 1 , kz
          do j = 1 , jxm1
            tqadv = tqadv - dtmin*3.E4*dsigma(k)                        &
                  & *dx*((vaix_g(k,j+1)+vaix_g(k,j))                    &
                  & *(qvailx_g(k,j)/psailx_g(j))                        &
                  & /(msfx_io(iym1,j)*msfx_io(iym1,j))                  &
                  & -(va01_g(k,j+1)+va01_g(k,j))                        &
                  & *(qva01_g(k,j)/psa01_g(j))                          &
                  & /(msfx_io(1,j)*msfx_io(1,j)))*rgti
          end do
        end do
      end if
      call mpi_bcast(tqadv,1,mpi_real8,0,mpi_comm_world,ierr)
#else
      do k = 1 , kz
        do j = 1 , jxm1
          tqadv = tqadv - dtmin*3.E4*dsigma(k)                          &
                & *dx*((va(iy,k,j+1)+va(iy,k,j))                        &
                & *(qva(iym1,k,j)/psa(iym1,j))/(msfx(iym1,j)            &
                & *msfx(iym1,j))-(va(1,k,j+1)+va(1,k,j))                &
                & *(qva(1,k,j)/psa(1,j))/(msfx(1,j)*msfx(1,j)))*rgti
        end do
      end do
#endif
!
!-----advection of cloud water and rainwater through lateral boundaries:
!
!.....advection through east-west boundaries:
!
#ifdef MPP1
      if ( myid.eq.nproc-1 ) then
        do k = 1 , kz
          do i = 1 , iym1
            worka(i,k) = (ua(i+1,k,jendl)+ua(i,k,jendl))                &
                       & *(qca(i,k,jendx)/psa(i,jendx))                 &
                       & /(msfx(i,jendx)*msfx(i,jendx))
          end do
        end do
      end if
      if ( myid.eq.0 ) then
        do k = 1 , kz
          do i = 1 , iym1
            workb(i,k) = (ua(i+1,k,1)+ua(i,k,1))*(qca(i,k,1)/psa(i,1))  &
                       & /(msfx(i,1)*msfx(i,1))
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
          worka(i,k) = (ua(i+1,k,jx)+ua(i,k,jx))                        &
                     & *(qca(i,k,jxm1)/psa(i,jxm1))                     &
                     & /(msfx(i,jxm1)*msfx(i,jxm1))
          workb(i,k) = (ua(i+1,k,1)+ua(i,k,1))*(qca(i,k,1)/psa(i,1))    &
                     & /(msfx(i,1)*msfx(i,1))
        end do
      end do
#endif
      do k = 1 , kz
        do i = 1 , iym1
          tqadv = tqadv - dtmin*3.E4*dsigma(k)                          &
                & *dx*(worka(i,k)-workb(i,k))*rgti
        end do
      end do
!
!....advection through north-south boundaries:
!
#ifdef MPP1
      do j = 1 , jendl
        do k = 1 , kz
          qcailx(k,j) = qca(iym1,k,j)
          qca01(k,j) = qca(1,k,j)
        end do
      end do
      call mpi_gather(qcailx,  kz*jxp,mpi_real8,                   &
                    & qcailx_g,kz*jxp,mpi_real8,                   &
                    & 0,mpi_comm_world,ierr)
      call mpi_gather(qca01,  kz*jxp,mpi_real8,                    &
                    & qca01_g,kz*jxp,mpi_real8,                    &
                    & 0,mpi_comm_world,ierr)
      if ( myid.eq.0 ) then
        do k = 1 , kz
          do j = 1 , jxm1
            tqadv = tqadv - dtmin*3.E4*dsigma(k)                        &
                  & *dx*((vaix_g(k,j+1)+vaix_g(k,j))                    &
                  & *(qcailx_g(k,j)/psailx_g(j))                        &
                  & /(msfx_io(iym1,j)*msfx_io(iym1,j))                  &
                  & -(va01_g(k,j+1)+va01_g(k,j))                        &
                  & *(qca01_g(k,j)/psa01_g(j))                          &
                  & /(msfx_io(1,j)*msfx_io(1,j)))*rgti
          end do
        end do
      end if
      call mpi_bcast(tqadv,1,mpi_real8,0,mpi_comm_world,ierr)
#else
      do k = 1 , kz
        do j = 1 , jxm1
          tqadv = tqadv - dtmin*3.E4*dsigma(k)                          &
                & *dx*((va(iy,k,j+1)+va(iy,k,j))                        &
                & *(qca(iym1,k,j)/psa(iym1,j))/(msfx(iym1,j)            &
                & *msfx(iym1,j))-(va(1,k,j+1)+va(1,k,j))                &
                & *(qca(1,k,j)/psa(1,j))/(msfx(1,j)*msfx(1,j)))*rgti
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
      error1 = 0.
      error2 = 0.
!
!-----compute the total dry air and water substance in the model at
!     this time:
!
!=======================================================================
!
!-----dry air (unit = kg):
!
      tdrym = 0.
#ifdef MPP1
      call mpi_gather(psa,   iy*jxp,mpi_real8,                     &
                    & psa_io,iy*jxp,mpi_real8,                     &
                    & 0,mpi_comm_world,ierr)
      if ( myid.eq.0 ) then
        do k = 1 , kz
          tttmp = 0.
          do j = 1 , jxm1
            do i = 1 , iym1
              tttmp = tttmp + psa_io(i,j)
            end do
          end do
          tdrym = tdrym + tttmp*dsigma(k)
        end do
        tdrym = tdrym*dx*dx*1000.*rgti
      end if
      call mpi_bcast(tdrym,1,mpi_real8,0,mpi_comm_world,ierr)
#else
      do k = 1 , kz
        tttmp = 0.
        do j = 1 , jxm1
          do i = 1 , iym1
            tttmp = tttmp + psa(i,j)
          end do
        end do
        tdrym = tdrym + tttmp*dsigma(k)
      end do
      tdrym = tdrym*dx*dx*1000.*rgti
#endif
!
!-----water substance (unit = kg):
!
      tvmass = 0.
#ifdef MPP1
      call mpi_gather(qva,   iy*kz*jxp,mpi_real8,                &
                    & qva_io,iy*kz*jxp,mpi_real8,                &
                    & 0,mpi_comm_world,ierr)
      if ( myid.eq.0 ) then
        do k = 1 , kz
          tttmp = 0.
          do j = 1 , jxm1
            do i = 1 , iym1
              tttmp = tttmp + qva_io(i,k,j)
            end do
          end do
          tvmass = tvmass + tttmp*dsigma(k)
        end do
        tvmass = tvmass*dx*dx*1000.*rgti
      end if
      call mpi_bcast(tvmass,1,mpi_real8,0,mpi_comm_world,               &
                   & ierr)
#else
      do k = 1 , kz
        tttmp = 0.
        do j = 1 , jxm1
          do i = 1 , iym1
            tttmp = tttmp + qva(i,k,j)
          end do
        end do
        tvmass = tvmass + tttmp*dsigma(k)
      end do
      tvmass = tvmass*dx*dx*1000.*rgti
#endif
!
      tcmass = 0.

#ifdef MPP1
      call mpi_gather(qca,   iy*kz*jxp,mpi_real8,                &
                    & qca_io,iy*kz*jxp,mpi_real8,                &
                    & 0,mpi_comm_world,ierr)
      if ( myid.eq.0 ) then
        do k = 1 , kz
          tttmp = 0.
          do j = 1 , jxm1
            do i = 1 , iym1
              tttmp = tttmp + qca_io(i,k,j)
            end do
          end do
          tcmass = tcmass + tttmp*dsigma(k)
        end do
        tcmass = tcmass*dx*dx*1000.*rgti
      end if
      call mpi_bcast(tcmass,1,mpi_real8,0,mpi_comm_world,               &
                   & ierr)
#else
      do k = 1 , kz
        tttmp = 0.
        do j = 1 , jxm1
          do i = 1 , iym1
            tttmp = tttmp + qca(i,k,j)
          end do
        end do
        tcmass = tcmass + tttmp*dsigma(k)
      end do
      tcmass = tcmass*dx*dx*1000.*rgti
#endif

      tqmass = tvmass + tcmass

!=======================================================================
!
!-----conservation of dry air:
!
      tdrym = tdrym - tdadv
      error1 = (tdrym-tdini)/tdini*100.
!
!-----conservation of water substance:
!
!-----total raifall at this time:
!
#ifdef MPP1
      call mpi_gather(rainc,   iy*jxp,mpi_real8,                   &
                    & rainc_io,iy*jxp,mpi_real8,                   &
                    & 0,mpi_comm_world,ierr)
      call mpi_gather(rainnc,   iy*jxp,mpi_real8,                  &
                    & rainnc_io,iy*jxp,mpi_real8,                  &
                    & 0,mpi_comm_world,ierr)
      if ( myid.eq.0 ) then
        tcrai = 0.
        tncrai = 0.
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
      tcrai = 0.
      tncrai = 0.
      do j = 1 , jxm1
        do i = 1 , iym1
          tcrai = tcrai + rainc(i,j)*dxsq
          tncrai = tncrai + rainnc(i,j)*dxsq
        end do
      end do
      tqrai = tcrai + tncrai
#endif
!
      tqmass = tqmass + tqrai - tqeva - tqadv
      error2 = (tqmass-tqini)/tqini*100.
!
!-----print out the information:
!
#ifdef MPP1
      if ( myid.eq.0 ) then
        if ( debug_level > 3 .and. mod(ntime,ndbgfrq).eq.0 ) then
          xh = xtime/1440.
          print * , '***** day = ' , ldatez + xh , ' *****'
          print 99001 , tdrym , error1
          print 99002 , tdadv
          print 99003 , tqmass , error2
          print 99004 , tvmass
          print 99005 , tcmass
          print 99006 , tqadv
          print 99007 , tcrai
          print 99008 , tncrai
          print 99009 , tqeva
        end if
      end if
#else
      if ( debug_level > 3 .and. mod(ntime,ndbgfrq).eq.0 ) then
        xh = xtime/1440.
        print * , '***** day = ' , ldatez + xh , ' *****'
        print 99001 , tdrym , error1
        print 99002 , tdadv
        print 99003 , tqmass , error2
        print 99004 , tvmass
        print 99005 , tcmass
        print 99006 , tqadv
        print 99007 , tcrai
        print 99008 , tncrai
        print 99009 , tqeva
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
!     real(kind=8)  chixp1,chiym1,chiyp1,chiym1,chi00,chidx,chidy
!     real(kind=8)  chixp2,chiym2,chiyp2,chiym2
 
!ccccccccccccccccccccccccccccccccccccccccccccccc
 
!-------------------------
!     1  ADVECTION budgets
!-----------------------
 
!-----advection of tracer through lateral boundaries
 
!
!.....advection through east-west boundaries:
!
!     the 'relaxed' upstream scheme
      fact1 = 0.6
      fact2 = 1 - fact1
 
!     inflow/outflow
#ifdef MPP1
      do n = 1 , ntr
        do k = 1 , kz
          do i = 2 , iym2
            if ( myid.eq.nproc-1 ) then
              uavg2 = 0.5*(ua(i+1,k,jendx)+ua(i,k,jendx))
              if ( uavg2.lt.0. ) then
                worka(i,k,n) = -uavg2*(fact1*chia(i,k,jendx,n)/psa(i,   &
                             & jendx)/(msfx(i,jendx)*msfx(i,jendx))     &
                             & +fact2*chia(i,k,jendm,n)/psa(i,jendm)    &
                             & /(msfx(i,jendm)*msfx(i,jendm)))
              else
                worka(i,k,n) = -uavg2*(fact1*chia(i,k,jendm,n)/psa(i,   &
                             & jendm)/(msfx(i,jendm)*msfx(i,jendm))     &
                             & +fact2*chia(i,k,jendx,n)/psa(i,jendx)    &
                             & /(msfx(i,jendx)*msfx(i,jendx)))
              end if
            end if
            if ( myid.eq.0 ) then
              uavg1 = 0.5*(ua(i+1,k,1+1)+ua(i,k,1+1))
              if ( uavg1.gt.0. ) then
                workb(i,k,n) = -uavg1*(fact1*chia(i,k,1,n)/psa(i,1)/(   &
                             & msfx(i,1)*msfx(i,1))                     &
                             & +fact2*chia(i,k,1+1,n)/psa(i,1+1)        &
                             & /(msfx(i,1+1)*msfx(i,1+1)))
              else
                workb(i,k,n) = -uavg1*(fact1*chia(i,k,1+1,n)/psa(i,1+1) &
                             & /(msfx(i,1+1)*msfx(i,1+1))               &
                             & +fact2*chia(i,k,1,n)/psa(i,1)            &
                             & /(msfx(i,1)*msfx(i,1)))
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
            uavg2 = 0.5*(ua(i+1,k,jxm1)+ua(i,k,jxm1))
            if ( uavg2.lt.0. ) then
              worka(i,k,n) = -uavg2*(fact1*chia(i,k,jxm1,n)/psa(i,jxm1) &
                           & /(msfx(i,jxm1)*msfx(i,jxm1))               &
                           & +fact2*chia(i,k,jxm2,n)/psa(i,jxm2)        &
                           & /(msfx(i,jxm2)*msfx(i,jxm2)))
            else
              worka(i,k,n) = -uavg2*(fact1*chia(i,k,jxm2,n)/psa(i,jxm2  &
                           & )/(msfx(i,jxm2)*msfx(i,jxm2))              &
                           & +fact2*chia(i,k,jxm1,n)/psa(i,jxm1)        &
                           & /(msfx(i,jxm1)*msfx(i,jxm1)))
            end if
 
            uavg1 = 0.5*(ua(i+1,k,1+1)+ua(i,k,1+1))
            if ( uavg1.gt.0. ) then
              workb(i,k,n) = -uavg1*(fact1*chia(i,k,1,n)/psa(i,1)/(msfx(&
                           & i,1)*msfx(i,1))+fact2*chia(i,k,1+1,n)      &
                           & /psa(i,1+1)/(msfx(i,1+1)*msfx(i,1+1)))
            else
              workb(i,k,n) = -uavg1*(fact1*chia(i,k,1+1,n)/psa(i,1+1)   &
                           & /(msfx(i,1+1)*msfx(i,1+1))                 &
                           & +fact2*chia(i,k,1,n)/psa(i,1)              &
                           & /(msfx(i,1)*msfx(i,1)))
            end if
          end do
        end do
      end do
#endif 

#ifdef MPP1
      do j = 1 , jendl
        do k = 1 , kz
          vaill(k,j) = va(iym1,k,j)
          va02(k,j) = va(2,k,j)
          xkcill1(k,j) = xkc(iym2,k,j)
          xkc02(k,j) = xkc(2,k,j)
          do n = 1 , ntr
            chiaill(k,n,j) = chia(iym1,k,j,n)
            chiaill1(k,n,j) = chia(iym2,k,j,n)
            chia01(k,n,j) = chia(1,k,j,n)
            chia02(k,n,j) = chia(2,k,j,n)
          end do
        end do
        psaill(j) = psa(iym1,j)
        psaill1(j) = psa(iym2,j)
        psa01(j) = psa(1,j)
        psa02(j) = psa(2,j)
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
      if ( myid.eq.0 ) then
        do n = 1 , ntr
          do k = 1 , kz
            do i = 2 , iym2
              tchiad(n) = tchiad(n) + dtmin*6.E4*dsigma(k)              &
                        & *dx*(worka(i,k,n)-workb(i,k,n))*rgti
            end do
          end do
!.....
!.....advection through north-south boundaries:
!
          do k = 1 , kz
            do j = 2 , jxm2
!hy           inflow/outflow
              vavg2 = 0.5*(vaill_g(k,j+1)+vaill_g(k,j))
              if ( vavg2.lt.0. ) then
                fx2 = -vavg2*(fact1*chiaill_g(k,n,j)/psaill_g(j)        &
                    & /(msfx_io(iym1,j)*msfx_io(iym1,j))                &
                    & +fact2*chiaill1_g(k,n,j)/psaill1_g(j)             &
                    & /(msfx_io(iym2,j)*msfx_io(iym2,j)))
              else
                fx2 = -vavg2*(fact1*chiaill1_g(k,n,j)/psaill1_g(j)      &
                    & /(msfx_io(iym2,j)*msfx_io(iym2,j))                &
                    & +fact2*chiaill_g(k,n,j)/psaill_g(j)               &
                    & /(msfx_io(iym1,j)*msfx_io(iym1,j)))
              end if
 
              vavg1 = 0.5*(va02_g(k,j+1)+va02_g(k,j))
              if ( vavg1.gt.0. ) then
                fx1 = -vavg1*(fact1*chia01_g(k,n,j)/psa01_g(j)          &
                    & /(msfx_io(1,j)*msfx_io(1,j))+fact2*chia02_g(k,n,j)&
                    & /psa02_g(j)/(msfx_io(2,j)*msfx_io(2,j)))
              else
                fx1 = -vavg1*(fact1*chia02_g(k,n,j)/psa02_g(j)          &
                    & /(msfx_io(2,j)*msfx_io(2,j))+fact2*chia01_g(k,n,j)&
                    & /psa01_g(j)/(msfx_io(1,j)*msfx_io(1,j)))
              end if
              tchiad(n) = tchiad(n) + dtmin*6.E4*dsigma(k)*dx*(fx2-fx1) &
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
            tchiad(n) = tchiad(n) + dtmin*6.E4*dsigma(k)                &
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
            vavg2 = 0.5*(va(iym1,k,j+1)+va(iym1,k,j))
            if ( vavg2.lt.0. ) then
              fx2 = -vavg2*(fact1*chia(iym1,k,j,n)/psa(iym1,j)          &
                  & /(msfx(iym1,j)*msfx(iym1,j))+fact2*chia(iym2,k,j,n) &
                  & /psa(iym2,j)/(msfx(iym2,j)*msfx(iym2,j)))
            else
              fx2 = -vavg2*(fact1*chia(iym2,k,j,n)/psa(iym2,j)          &
                  & /(msfx(iym2,j)*msfx(iym2,j))                        &
                  & +fact2*chia(iym1,k,j,n)/psa(iym1,j)                 &
                  & /(msfx(iym1,j)*msfx(iym1,j)))
            end if
 
            vavg1 = 0.5*(va(1+1,k,j+1)+va(1+1,k,j))
            if ( vavg1.gt.0. ) then
              fx1 = -vavg1*(fact1*chia(1,k,j,n)/psa(1,j)                &
                  & /(msfx(1,j)*msfx(1,j))+fact2*chia(1+1,k,j,n)        &
                  & /psa(1+1,j)/(msfx(1+1,j)*msfx(1+1,j)))
            else
              fx1 = -vavg1*(fact1*chia(1+1,k,j,n)/psa(1+1,j)            &
                  & /(msfx(1+1,j)*msfx(1+1,j))+fact2*chia(1,k,j,n)      &
                  & /psa(1,j)/(msfx(1,j)*msfx(1,j)))
            end if
 
            tchiad(n) = tchiad(n) + dtmin*6.E4*dsigma(k)*dx*            &
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
            if ( myid.eq.nproc-1 ) worka(i,k,n) = xkc(i,k,jendm)        &
               & *psa(i,jendm)                                          &
               & *(chia(i,k,jendm,n)/psa(i,jendm)-chia(i,k,jendx,n)     &
               & /psa(i,jendx))
            if ( myid.eq.0 ) workb(i,k,n) = xkc(i,k,2)*psa(i,2)         &
               & *(chia(i,k,2,n)/psa(i,2)-chia(i,k,1,n)/psa(i,1))
          end do
        end do
      end do
      call mpi_bcast(worka,iym1*kz*ntr,mpi_real8,nproc-1,               &
                   & mpi_comm_world,ierr)
#else
      do n = 1 , ntr
        do k = 1 , kz
          do i = 2 , iym2
            worka(i,k,n) = xkc(i,k,jxm2)*psa(i,jxm2)                    &
                         & *(chia(i,k,jxm2,n)/psa(i,jxm2)               &
                         & -chia(i,k,jxm1,n)/psa(i,jxm1))
            workb(i,k,n) = xkc(i,k,2)*psa(i,2)                          &
                         & *(chia(i,k,2,n)/psa(i,2)-chia(i,k,1,n)       &
                         & /psa(i,1))
          end do
        end do
      end do
#endif

#ifdef MPP1
      if ( myid.eq.0 ) then
        do n = 1 , ntr
          do k = 1 , kz
            do i = 2 , iym2
              tchitb(n) = tchitb(n) - dtmin*6.E4*dsigma(k)              &
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
              tchitb(n) = tchitb(n) - dtmin*6.E4*dsigma(k)*(chid2+chid1)&
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
            tchitb(n) = tchitb(n) - dtmin*6.E4*dsigma(k)                &
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
        ttrace(itr,1) = 0.
        tremlsc(itr,1) = 0.
        tremcvc(itr,1) = 0.
        trxsg(itr,1) = 0.
        trxsaq1(itr,1) = 0.
        trxsaq2(itr,1) = 0.
        tremdrd(itr,1) = 0.
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
      if ( myid.eq.0 ) then
        do itr = 1 , ntr
          do k = 1 , kz
            tttmp = 0.
            do j = 2 , jxm2
              do i = 2 , iym2
                tttmp = tttmp + chia_io(i,k,j,itr)
              end do
            end do
            ttrace(itr,1) = ttrace(itr,1) + tttmp*dsigma(k)
            tttmp = 0.
            do j = 2 , jxm2
              do i = 2 , iym2
                tttmp = tttmp + remlsc_io(i,k,j,itr)
              end do
            end do
            tremlsc(itr,1) = tremlsc(itr,1) + tttmp*dsigma(k)
            tttmp = 0.
            do j = 2 , jxm2
              do i = 2 , iym2
                tttmp = tttmp + remcvc_io(i,k,j,itr)
              end do
            end do
            tremcvc(itr,1) = tremcvc(itr,1) + tttmp*dsigma(k)
            tttmp = 0.
            do j = 2 , jxm2
              do i = 2 , iym2
                tttmp = tttmp + rxsg_io(i,k,j,itr)
              end do
            end do
            trxsg(itr,1) = trxsg(itr,1) + tttmp*dsigma(k)
            tttmp = 0.
            do j = 2 , jxm2
              do i = 2 , iym2
                tttmp = tttmp + rxsaq1_io(i,k,j,itr)
              end do
            end do
            trxsaq1(itr,1) = trxsaq1(itr,1) + tttmp*dsigma(k)
            tttmp = 0.
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
          tttmp = 0.
          do j = 2 , jxm2
            do i = 2 , iym2
              tttmp = tttmp + remdrd_io(i,j,itr)
            end do
          end do
          tremdrd(itr,1) = tremdrd(itr,1) + tttmp*dx*dx*dsigma(kz)
 
!         emissions
          tttmp = 0.
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
          tttmp = 0.
          do j = 2 , jxm2
            do i = 2 , iym2
              tttmp = tttmp + chia(i,k,j,itr)
            end do
          end do
          ttrace(itr,1) = ttrace(itr,1) + tttmp*dsigma(k)
          tttmp = 0.
          do j = 2 , jxm2
            do i = 2 , iym2
              tttmp = tttmp + remlsc(i,k,j,itr)
            end do
          end do
          tremlsc(itr,1) = tremlsc(itr,1) + tttmp*dsigma(k)
          tttmp = 0.
          do j = 2 , jxm2
            do i = 2 , iym2
              tttmp = tttmp + remcvc(i,k,j,itr)
            end do
          end do
          tremcvc(itr,1) = tremcvc(itr,1) + tttmp*dsigma(k)
          tttmp = 0.
          do j = 2 , jxm2
            do i = 2 , iym2
              tttmp = tttmp + rxsg(i,k,j,itr)
            end do
          end do
          trxsg(itr,1) = trxsg(itr,1) + tttmp*dsigma(k)
          tttmp = 0.
          do j = 2 , jxm2
            do i = 2 , iym2
              tttmp = tttmp + rxsaq1(i,k,j,itr)
            end do
          end do
          trxsaq1(itr,1) = trxsaq1(itr,1) + tttmp*dsigma(k)
          tttmp = 0.
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
        tttmp = 0.
        do j = 2 , jxm2
          do i = 2 , iym2
            tttmp = tttmp + remdrd(i,j,itr)
          end do
        end do
        tremdrd(itr,1) = tremdrd(itr,1) + tttmp*dx*dx*dsigma(kz)

!       emissions
        tttmp = 0.
        do j = 2 , jxm2
          do i = 2 , iym2
            tttmp = tttmp + chemsrc(i,j,lmonth,itr)*dtmin*60.*dx*dx
          end do
        end do
        tchie(itr) = tchie(itr) + tttmp
      end do
#endif
 
      do itr = 1 , ntr
        ttrace(itr,1) = ttrace(itr,1)*1000.*rgti
        tremlsc(itr,1) = tremlsc(itr,1)*1000.*rgti
        tremcvc(itr,1) = tremcvc(itr,1)*1000.*rgti
        tremdrd(itr,1) = tremdrd(itr,1)*1000.*rgti
        trxsg(itr,1) = trxsg(itr,1)*1000.*rgti
        trxsaq1(itr,1) = trxsaq1(itr,1)*1000.*rgti
        trxsaq2(itr,1) = trxsaq2(itr,1)*1000.*rgti
      end do

 
!-----print out the information:
#ifdef MPP1
      if ( myid.eq.0 ) then
 
        if ( debug_level > 3 .and. mod(ntime,ndbgfrq).eq.0 ) then
 
!----     tracers
 
          print * , '************************************************'
          print * , ' Budgets for tracers (intergrated quantitites)'
          print * , ' day = ' , ldatez , ' *****'
          print * , '************************************************'
 
          do itr = 1 , ntr
 
!           errort= ttrace(itr,1)-
!           &     (tchie(itr) - (tchiad(itr)-tchitb(itr)+
!           &     tremdrd(itr,1)+ tremlsc(itr,1)+tremcvc(itr,1)))
 
            print * , '****************************'
            print 99001 , itr , tchie(itr)
            print 99002 , itr , tchiad(itr)
            print 99003 , itr , tchitb(itr)
            print 99004 , itr , tremdrd(itr,1)
            print 99005 , itr , tremlsc(itr,1)
            print 99006 , itr , tremcvc(itr,1)
            print 99007 , itr , ttrace(itr,1)
!           print 99008 , itr , errort
 
          end do
 
        end if
      end if
#else
      if ( debug_level > 3 .and. mod(ntime,ndbgfrq).eq.0 ) then

!----   tracers

        print * , '************************************************'
        print * , ' Budgets for tracers (intergrated quantitites)'
        print * , ' day = ' , ldatez , ' *****'
        print * , '************************************************'

        do itr = 1 , ntr

!         errort= ttrace(itr,1)-
!         &     (tchie(itr) - (tchiad(itr)-tchitb(itr)+
!         &     tremdrd(itr,1)+ tremlsc(itr,1)+tremcvc(itr,1)))

          print * , '****************************'
          print 99001 , itr , tchie(itr)
          print 99002 , itr , tchiad(itr)
          print 99003 , itr , tchitb(itr)
          print 99004 , itr , tremdrd(itr,1)
          print 99005 , itr , tremlsc(itr,1)
          print 99006 , itr , tremcvc(itr,1)
          print 99007 , itr , ttrace(itr,1)
!         print 99008 , itr , errort
        end do
      end if
#endif
! 
99001       format ('total emission tracer',i1,':',e12.5,' kg')
99002       format ('total advected tracer',i1,':',e12.5,' kg')
99003       format ('total diffused tracer',i1,':',e12.5,' kg')
99004       format ('total dry deposition tracer',i1,':',e12.5,' kg')
99005       format ('total large scale wet dep.tracer',i1,':',e12.5,    &
                   &' kg')
99006       format ('total convective wet dep.tracer',i1,':',e12.5,     &
                  & ' kg')
99007       format ('total mass (at ktau), tracer',i1,':',e12.5,' kg')
!99008       format('error trac',I1,':',e12.5,' % ')
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
      call mpi_gather(qfx,   iy*jxp,mpi_real8,  &
                    & qfx_io,iy*jxp,mpi_real8,  &
                    & 0,mpi_comm_world,ierr)
      if ( myid.eq.0 ) then
        do j = 2 , jxm2
          do i = 2 , iym2
            tqeva = tqeva + qfx_io(i,j)*dx*dx*dtmin*60.
          end do
        end do
      end if
      call mpi_bcast(tqeva,1,mpi_real8,0,mpi_comm_world,ierr)
#else
      do j = 2 , jxm2
        do i = 2 , iym2
          tqeva = tqeva + qfx(i,j)*dx*dx*dtmin*60.
        end do
      end do
#endif
      end subroutine conqeva
!
#endif
!
      end module mod_diagnosis
