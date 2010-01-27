C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C    This file is part of RegCM model.
C
C    RegCM model is free software: you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation, either version 3 of the License, or
C    (at your option) any later version.
C
C    RegCM model is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License
C    along with RegCM model.  If not, see <http://www.gnu.org/licenses/>.
C
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      subroutine conmas
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                     c
c     this subroutine computes the total dry air and water substance  c
c     within the domain and compares with the initial values.         c
c                                                                     c
c     ---the unit used in all the calculation is "kg".                c
c                                                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
#     include "parameter.inc"
#     include "include/main.cb"
#ifdef MPP1
#     include "mpiregcm.inc"
      integer ierr
#     include "include/io.cb"
#endif
#     include "include/param1.cb"
#     include "include/diagnosis.cb"
#     include "include/param2.cb"
#     include "include/param3.cb"
#     include "include/date.cb"
c
      real(kind=8)  error1,error2,tttmp,tdrym,tqmass,tvmass,tcmass
      real(kind=8)  tcrai,tncrai,xh
      integer i,j,k
c
c----------------------------------------------------------------------
c
      error1 = 0.
      error2 = 0.
c
c-----compute the total dry air and water substance in the model at
c     this time:
c
C=======================================================================
c
c-----dry air (unit = kg):
c
      tdrym = 0.
#ifdef MPP1
      call MPI_GATHER(psa(1,1),ix*jxp,MPI_DOUBLE_PRECISION,
     &             psa_io(1,1),ix*jxp,MPI_DOUBLE_PRECISION,
     &                              0,MPI_COMM_WORLD,ierr)
      if(myid.eq.0) then
        do k=1,kx
          tttmp= 0.
          do j=1,mjx-1
            do i=1,ix-1
              tttmp=tttmp+psa_io(i,j)
            end do
          end do
          tdrym=tdrym+tttmp*dsigma(k)
        end do
        tdrym=tdrym*dx*dx*1000./g
      end if
      call MPI_BCAST(tdrym,1,MPI_DOUBLE_PRECISION,
     &                     0,MPI_COMM_WORLD,ierr)
#else
      do k=1,kx
        tttmp= 0.
        do j=1,jx-1
          do i=1,ix-1
            tttmp=tttmp+psa(i,j)
          end do
        end do
        tdrym=tdrym+tttmp*dsigma(k)
      end do
      tdrym=tdrym*dx*dx*1000./g
#endif
c
c-----water substance (unit = kg):
c
      tvmass = 0.
#ifdef MPP1
      call MPI_GATHER(qva(1,1,1),ix*kx*jxp,MPI_DOUBLE_PRECISION,
     &             qva_io(1,1,1),ix*kx*jxp,MPI_DOUBLE_PRECISION,
     &                                   0,MPI_COMM_WORLD,ierr)
      if(myid.eq.0) then
        do k=1,kx
          tttmp= 0.
          do j=1,mjx-1
            do i=1,ix-1
              tttmp=tttmp+qva_io(i,k,j)
            end do
          end do
          tvmass=tvmass+tttmp*dsigma(k)
        end do
        tvmass=tvmass*dx*dx*1000./g
      end if
      call MPI_BCAST(tvmass,1,MPI_DOUBLE_PRECISION,
     &                      0,MPI_COMM_WORLD,ierr)
#else
      do k=1,kx
        tttmp= 0.
        do j=1,jx-1
          do i=1,ix-1
            tttmp=tttmp+qva(i,k,j)
          end do
        end do
        tvmass=tvmass+tttmp*dsigma(k)
      end do
      tvmass=tvmass*dx*dx*1000./g
#endif
c
      tcmass = 0.
#ifdef MPP1
      call MPI_GATHER(qca(1,1,1),ix*kx*jxp,MPI_DOUBLE_PRECISION,
     &             qca_io(1,1,1),ix*kx*jxp,MPI_DOUBLE_PRECISION,
     &                                   0,MPI_COMM_WORLD,ierr)
      if(myid.eq.0) then
        do k=1,kx
          tttmp= 0.
          do j=1,mjx-1
            do i=1,ix-1
              tttmp=tttmp+qca_io(i,k,j)
            end do
          end do
          tcmass=tcmass+tttmp*dsigma(k)
        end do
        tcmass=tcmass*dx*dx*1000./g
      end if
      call MPI_BCAST(tcmass,1,MPI_DOUBLE_PRECISION,
     &                      0,MPI_COMM_WORLD,ierr)
#else
      do k=1,kx
        tttmp= 0.
        do j=1,jx-1
          do i=1,ix-1
            tttmp=tttmp+qca(i,k,j)
          end do
        end do
        tcmass=tcmass+tttmp*dsigma(k)
      end do
      tcmass=tcmass*dx*dx*1000./g
#endif
      tqmass=tvmass+tcmass
C=======================================================================
c
c-----conservation of dry air:
c
      tdrym=tdrym-tdadv
      error1=(tdrym-tdini)/tdini*100.
c
c-----conservation of water substance:
c
c-----total raifall at this time:
c
#ifdef MPP1
      call MPI_GATHER(rainc(1,1),ix*jxp,MPI_DOUBLE_PRECISION,
     &             rainc_io(1,1),ix*jxp,MPI_DOUBLE_PRECISION,
     &                                0,MPI_COMM_WORLD,ierr)
      call MPI_GATHER(rainnc(1,1),ix*jxp,MPI_DOUBLE_PRECISION,
     &             rainnc_io(1,1),ix*jxp,MPI_DOUBLE_PRECISION,
     &                                 0,MPI_COMM_WORLD,ierr)
      if(myid.eq.0) then
        tcrai = 0.
        tncrai = 0.
        do j=1,mjx-1
        do i=1,ilx
         tcrai =tcrai + rainc_io(i,j)*dxsq
         tncrai=tncrai+rainnc_io(i,j)*dxsq
        end do
        end do
        tqrai=tcrai+tncrai
      end if
      call MPI_BCAST(tqrai,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
#else
      tcrai = 0.
      tncrai = 0.
      do j=1,jlx
      do i=1,ilx
       tcrai=tcrai+rainc(i,j)*dxsq
       tncrai=tncrai+rainnc(i,j)*dxsq
      end do
      end do
      tqrai=tcrai+tncrai
#endif
c
      tqmass=tqmass+tqrai-tqeva-tqadv
      error2=(tqmass-tqini)/tqini*100.
c
c-----print out the information:
c
#ifdef MPP1
      if (myid .eq. 0) then
#endif
      if (ifprt .and. mod(ntime,nprtfrq).eq.0) then
        xh = xtime / 1440.
        print *,'***** day = ',ldatez+xh,' *****'
        print 31,tdrym,error1
  31    format('   total air =',e12.5,' kg, error = ',e12.5)
        print 32,tdadv
  32    format('   horizontal advection = ',e12.5,' kg.')
        print 41,tqmass,error2
  41    format('   total water =',e12.5,' kg, error = ',e12.5)
        print 42,tvmass
  42    format('   qv                      = ',e12.5,' kg.')
        print 43,tcmass
  43    format('   qc                      = ',e12.5,' kg.')
        print 45,tqadv
  45    format('   horizontal advection    = ',e12.5,' kg.')
        print 46,tcrai
  46    format('   convective railfall     = ',e12.5,' kg.')
        print 47,tncrai
  47    format('   nonconvective rainfall  = ',e12.5,' kg.')
        print 48,tqeva
  48    format('   evaporation from ground = ',e12.5,' kg.')
      end if
#ifdef MPP1
      end if
#endif

      return
      end subroutine conmas
