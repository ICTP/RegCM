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
 
      subroutine tracbud
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use mod_regcm_param
      use mod_param1
      use mod_param2
      use mod_param3 , only : dsigma
      use mod_main
      use mod_mainchem
      use mod_trachem
      use mod_date
      use mod_constants , only : rgti
#ifdef DIAG
      use mod_diagnosis
#endif
#ifdef MPP1
      use mod_mppio
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
      integer :: i , itr , j , k
#ifdef DIAG
#ifdef MPP1
      integer :: ierr , l
#endif
      real(8) :: tttmp
#endif
!
!---- add tracers,  tremlsc, tremcvc, tremdrd are total amount mass
!     and chemical conversion term
!     removed by large scale, convective precipation and dry deposition
!     respectively
 
      do itr = 1 , ntr
!
#ifdef DIAG
        ttrace(itr,1) = 0.
        tremlsc(itr,1) = 0.
        tremcvc(itr,1) = 0.
        trxsg(itr,1) = 0.
        trxsaq1(itr,1) = 0.
        trxsaq2(itr,1) = 0.
        tremdrd(itr,1) = 0.
#endif
!
#ifdef MPP1
        do j = 1 , jendx
#else
        do j = 1 , jxm1
#endif
          do i = 1 , ixm1
            dtrace(i,j,itr) = 0.0
            wdlsc(i,j,itr) = 0.0
            wdcvc(i,j,itr) = 0.0
            wxsg(i,j,itr) = 0.0
            wxaq(i,j,itr) = 0.0
            ddsfc(i,j,itr) = 0.0
          end do
        end do
      end do
 
 
!-----tracers (unit = kg):
      do itr = 1 , ntr
#ifdef MPP1
        do j = 1 , jendx
#else
        do j = 1 , jxm1
#endif
          do i = 1 , ixm1
            do k = 1 , kx
              dtrace(i,j,itr) = dtrace(i,j,itr) + chia(i,k,j,itr)       &
                              & *dsigma(k)
              wdlsc(i,j,itr) = wdlsc(i,j,itr) + remlsc(i,k,j,itr)       &
                             & *dsigma(k)
              wdcvc(i,j,itr) = wdcvc(i,j,itr) + remcvc(i,k,j,itr)       &
                             & *dsigma(k)
              wxsg(i,j,itr) = wxsg(i,j,itr) + rxsg(i,k,j,itr)*dsigma(k)
!             sum ls and conv contribution
              wxaq(i,j,itr) = wxaq(i,j,itr)                             &
                            & + (rxsaq1(i,k,j,itr)+rxsaq2(i,k,j,itr))   &
                            & *dsigma(k)
            end do
            ddsfc(i,j,itr) = ddsfc(i,j,itr) + remdrd(i,j,itr)*dsigma(kx)
!           Source cumulated diag(care the unit are alredy .m-2)
            cemtrac(i,j,itr) = cemtr(i,j,itr)
          end do
        end do
      end do
 
#ifdef DIAG
#ifdef MPP1
      do itr = 1 , ntr
        call mpi_gather(chia(1,1,1,itr),   ix*kx*jxp,mpi_real8,         &
                      & chia_io(1,1,1,itr),ix*kx*jxp,mpi_real8,         &
                      & 0,mpi_comm_world,ierr)
        call mpi_gather(remlsc(1,1,1,itr),   ix*kx*jxp,mpi_real8,       &
                      & remlsc_io(1,1,1,itr),ix*kx*jxp,mpi_real8,       &
                      & 0,mpi_comm_world,ierr)
        call mpi_gather(remcvc(1,1,1,itr),   ix*kx*jxp,mpi_real8,       &
                      & remcvc_io(1,1,1,itr),ix*kx*jxp,mpi_real8,       &
                      & 0,mpi_comm_world,ierr)
        call mpi_gather(rxsg(1,1,1,itr),   ix*kx*jxp,mpi_real8,         &
                      & rxsg_io(1,1,1,itr),ix*kx*jxp,mpi_real8,         &
                      & 0,mpi_comm_world,ierr)
        call mpi_gather(rxsaq1(1,1,1,itr),   ix*kx*jxp,mpi_real8,       &
                      & rxsaq1_io(1,1,1,itr),ix*kx*jxp,mpi_real8,       &
                      & 0,mpi_comm_world,ierr)
        call mpi_gather(rxsaq2(1,1,1,itr),   ix*kx*jxp,mpi_real8,       &
                      & rxsaq2_io(1,1,1,itr),ix*kx*jxp,mpi_real8,       &
                      & 0,mpi_comm_world,ierr)
        call mpi_gather(remdrd(1,1,itr),   ix*jxp,mpi_real8,            &
                      & remdrd_io(1,1,itr),ix*jxp,mpi_real8,            &
                      & 0,mpi_comm_world,ierr)
        do l = 1 , 12
          call mpi_gather(chemsrc(1,1,l,itr),   ix*jxp,mpi_real8,       &
                        & chemsrc_io(1,1,l,itr),ix*jxp,mpi_real8,       &
                        & 0,mpi_comm_world,ierr)
        end do
      end do
      if ( myid.eq.0 ) then
        do itr = 1 , ntr
          do k = 1 , kx
            tttmp = 0.
            do j = 2 , jxm2
              do i = 2 , ixm2
                tttmp = tttmp + chia_io(i,k,j,itr)
              end do
            end do
            ttrace(itr,1) = ttrace(itr,1) + tttmp*dsigma(k)
            tttmp = 0.
            do j = 2 , jxm2
              do i = 2 , ixm2
                tttmp = tttmp + remlsc_io(i,k,j,itr)
              end do
            end do
            tremlsc(itr,1) = tremlsc(itr,1) + tttmp*dsigma(k)
            tttmp = 0.
            do j = 2 , jxm2
              do i = 2 , ixm2
                tttmp = tttmp + remcvc_io(i,k,j,itr)
              end do
            end do
            tremcvc(itr,1) = tremcvc(itr,1) + tttmp*dsigma(k)
            tttmp = 0.
            do j = 2 , jxm2
              do i = 2 , ixm2
                tttmp = tttmp + rxsg_io(i,k,j,itr)
              end do
            end do
            trxsg(itr,1) = trxsg(itr,1) + tttmp*dsigma(k)
            tttmp = 0.
            do j = 2 , jxm2
              do i = 2 , ixm2
                tttmp = tttmp + rxsaq1_io(i,k,j,itr)
              end do
            end do
            trxsaq1(itr,1) = trxsaq1(itr,1) + tttmp*dsigma(k)
            tttmp = 0.
            do j = 2 , jxm2
              do i = 2 , ixm2
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
            do i = 2 , ixm2
              tttmp = tttmp + remdrd_io(i,j,itr)
            end do
          end do
          tremdrd(itr,1) = tremdrd(itr,1) + tttmp*dx*dx*dsigma(kx)
 
!         emissions
          tttmp = 0.
          do j = 2 , jxm2
            do i = 2 , ixm2
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
        do k = 1 , kx
          tttmp = 0.
          do j = 2 , jxm2
            do i = 2 , ixm2
              tttmp = tttmp + chia(i,k,j,itr)
            end do
          end do
          ttrace(itr,1) = ttrace(itr,1) + tttmp*dsigma(k)
          tttmp = 0.
          do j = 2 , jxm2
            do i = 2 , ixm2
              tttmp = tttmp + remlsc(i,k,j,itr)
            end do
          end do
          tremlsc(itr,1) = tremlsc(itr,1) + tttmp*dsigma(k)
          tttmp = 0.
          do j = 2 , jxm2
            do i = 2 , ixm2
              tttmp = tttmp + remcvc(i,k,j,itr)
            end do
          end do
          tremcvc(itr,1) = tremcvc(itr,1) + tttmp*dsigma(k)
          tttmp = 0.
          do j = 2 , jxm2
            do i = 2 , ixm2
              tttmp = tttmp + rxsg(i,k,j,itr)
            end do
          end do
          trxsg(itr,1) = trxsg(itr,1) + tttmp*dsigma(k)
          tttmp = 0.
          do j = 2 , jxm2
            do i = 2 , ixm2
              tttmp = tttmp + rxsaq1(i,k,j,itr)
            end do
          end do
          trxsaq1(itr,1) = trxsaq1(itr,1) + tttmp*dsigma(k)
          tttmp = 0.
          do j = 2 , jxm2
            do i = 2 , ixm2
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
          do i = 2 , ixm2
            tttmp = tttmp + remdrd(i,j,itr)
          end do
        end do
        tremdrd(itr,1) = tremdrd(itr,1) + tttmp*dx*dx*dsigma(kx)

!       emissions
        tttmp = 0.
        do j = 2 , jxm2
          do i = 2 , ixm2
            tttmp = tttmp + chemsrc(i,j,lmonth,itr)*dtmin*60.*dx*dx
          end do
        end do
        tchie(itr) = tchie(itr) + tttmp
      end do
#endif
#endif
 
      do itr = 1 , ntr
#ifdef DIAG
        ttrace(itr,1) = ttrace(itr,1)*1000.*rgti
        tremlsc(itr,1) = tremlsc(itr,1)*1000.*rgti
        tremcvc(itr,1) = tremcvc(itr,1)*1000.*rgti
        tremdrd(itr,1) = tremdrd(itr,1)*1000.*rgti
        trxsg(itr,1) = trxsg(itr,1)*1000.*rgti
        trxsaq1(itr,1) = trxsaq1(itr,1)*1000.*rgti
        trxsaq2(itr,1) = trxsaq2(itr,1)*1000.*rgti
#endif
!
#ifdef MPP1
        do j = 1 , jendx
#else
        do j = 1 , jxm1
#endif
          do i = 1 , ixm1
            dtrace(i,j,itr) = 1.E6*dtrace(i,j,itr)*1000.*rgti
                                                        ! unit: mg/m2
            wdlsc(i,j,itr) = 1.E6*wdlsc(i,j,itr)*1000.*rgti
            wdcvc(i,j,itr) = 1.E6*wdcvc(i,j,itr)*1000.*rgti
            ddsfc(i,j,itr) = 1.E6*ddsfc(i,j,itr)*1000.*rgti
            wxsg(i,j,itr) = 1.E6*wxsg(i,j,itr)*1000.*rgti
            wxaq(i,j,itr) = 1.E6*wxaq(i,j,itr)*1000.*rgti
!           emtrac isbuilt from chsurfem so just need the 1e6*dt/2
!           factor to to pass im mg/m2
            cemtrac(i,j,itr) = 1.E6*cemtrac(i,j,itr)
          end do
        end do
      end do
 
!-----print out the information:
#ifdef DIAG
#ifdef MPP1
      if ( myid.eq.0 ) then
 
        if ( ifprt .and. mod(ntime,nprtfrq).eq.0 ) then
 
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
      if ( ifprt .and. mod(ntime,nprtfrq).eq.0 ) then

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
#endif

      end subroutine tracbud
