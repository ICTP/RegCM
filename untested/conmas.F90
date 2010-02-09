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
 
      subroutine conmas

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!     this subroutine computes the total dry air and water substance  c
!     within the domain and compares with the initial values.         c
!                                                                     c
!     ---the unit used in all the calculation is "kg".                c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use regcm_param
      use main
      use param1
      use diagnosis
      use param2
      use param3
      use date
      use main
#ifdef MPP1
      use mpi
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
      call mpi_gather(psa(1,1),ix*jxp,mpi_double_precision,psa_io(1,1), &
                    & ix*jxp,mpi_double_precision,0,mpi_comm_world,ierr)
      if ( myid.eq.0 ) then
        do k = 1 , kx
          tttmp = 0.
          do j = 1 , mjx - 1
            do i = 1 , ix - 1
              tttmp = tttmp + psa_io(i,j)
            end do
          end do
          tdrym = tdrym + tttmp*dsigma(k)
        end do
        tdrym = tdrym*dx*dx*1000./g
      end if
      call mpi_bcast(tdrym,1,mpi_double_precision,0,mpi_comm_world,ierr)
#else
      do k = 1 , kx
        tttmp = 0.
        do j = 1 , jx - 1
          do i = 1 , ix - 1
            tttmp = tttmp + psa(i,j)
          end do
        end do
        tdrym = tdrym + tttmp*dsigma(k)
      end do
      tdrym = tdrym*dx*dx*1000./g
#endif
!
!-----water substance (unit = kg):
!
      tvmass = 0.
#ifdef MPP1
      call mpi_gather(qva(1,1,1),ix*kx*jxp,mpi_double_precision,        &
                    & qva_io(1,1,1),ix*kx*jxp,mpi_double_precision,0,   &
                    & mpi_comm_world,ierr)
      if ( myid.eq.0 ) then
        do k = 1 , kx
          tttmp = 0.
          do j = 1 , mjx - 1
            do i = 1 , ix - 1
              tttmp = tttmp + qva_io(i,k,j)
            end do
          end do
          tvmass = tvmass + tttmp*dsigma(k)
        end do
        tvmass = tvmass*dx*dx*1000./g
      end if
      call mpi_bcast(tvmass,1,mpi_double_precision,0,mpi_comm_world,    &
                   & ierr)
#else
      do k = 1 , kx
        tttmp = 0.
        do j = 1 , jx - 1
          do i = 1 , ix - 1
            tttmp = tttmp + qva(i,k,j)
          end do
        end do
        tvmass = tvmass + tttmp*dsigma(k)
      end do
      tvmass = tvmass*dx*dx*1000./g
#endif
!
      tcmass = 0.

#ifdef MPP1
      call mpi_gather(qca(1,1,1),ix*kx*jxp,mpi_double_precision,        &
                    & qca_io(1,1,1),ix*kx*jxp,mpi_double_precision,0,   &
                    & mpi_comm_world,ierr)
      if ( myid.eq.0 ) then
        do k = 1 , kx
          tttmp = 0.
          do j = 1 , mjx - 1
            do i = 1 , ix - 1
              tttmp = tttmp + qca_io(i,k,j)
            end do
          end do
          tcmass = tcmass + tttmp*dsigma(k)
        end do
        tcmass = tcmass*dx*dx*1000./g
      end if
      call mpi_bcast(tcmass,1,mpi_double_precision,0,mpi_comm_world,    &
                   & ierr)
#else
      do k = 1 , kx
        tttmp = 0.
        do j = 1 , jx - 1
          do i = 1 , ix - 1
            tttmp = tttmp + qca(i,k,j)
          end do
        end do
        tcmass = tcmass + tttmp*dsigma(k)
      end do
      tcmass = tcmass*dx*dx*1000./g
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
      call mpi_gather(rainc(1,1),ix*jxp,mpi_double_precision,           &
                    & rainc_io(1,1),ix*jxp,mpi_double_precision,0,      &
                    & mpi_comm_world,ierr)
      call mpi_gather(rainnc(1,1),ix*jxp,mpi_double_precision,          &
                    & rainnc_io(1,1),ix*jxp,mpi_double_precision,0,     &
                    & mpi_comm_world,ierr)
      if ( myid.eq.0 ) then
        tcrai = 0.
        tncrai = 0.
        do j = 1 , mjx - 1
          do i = 1 , ilx
            tcrai = tcrai + rainc_io(i,j)*dxsq
            tncrai = tncrai + rainnc_io(i,j)*dxsq
          end do
        end do
        tqrai = tcrai + tncrai
      end if
      call mpi_bcast(tqrai,1,mpi_double_precision,0,mpi_comm_world,ierr)
#else
      tcrai = 0.
      tncrai = 0.
      do j = 1 , jlx
        do i = 1 , ilx
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
        if ( ifprt .and. mod(ntime,nprtfrq).eq.0 ) then
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
      if ( ifprt .and. mod(ntime,nprtfrq).eq.0 ) then
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
