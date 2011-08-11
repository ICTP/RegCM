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

  use mod_runparams
  use mod_main
  use mod_che_main
  use mod_che_trac
  use mod_mppio
  use mod_memutil
!
  private
!
  public :: allocate_mod_diagnosis , initdiag
  public :: restdiag , restchemdiag
  public :: savediag , savechemdiag
  public :: conadv , conmas , conqeva
  public :: tracdiag , contrac
  public :: mpidiag
!
  real(8) :: tdadv , tdini , tqadv , tqeva , tqini , tqrai
!
  real(8) , pointer , dimension(:) :: tchiad , tchie , tchitb
  real(8) , pointer , dimension(:,:) :: tremcvc , tremdrd ,     &
              tremlsc , trxsaq1 , trxsaq2 , trxsg , ttrace
  contains

  subroutine allocate_mod_diagnosis
  implicit none
    if ( ichem == 1 ) then
      call getmem1d(tchie,1,ntr,'diagnosis:tchie')
      call getmem1d(tchiad,1,ntr,'diagnosis:tchiad')
      call getmem1d(tchitb,1,ntr,'diagnosis:tchitb')
      call getmem2d(tremcvc,1,ntr,1,2,'diagnosis:tremcvc')
      call getmem2d(tremdrd,1,ntr,1,2,'diagnosis:tremdrd')
      call getmem2d(tremlsc,1,ntr,1,2,'diagnosis:tremlsc')
      call getmem2d(trxsaq1,1,ntr,1,2,'diagnosis:trxsaq1')
      call getmem2d(trxsaq2,1,ntr,1,2,'diagnosis:trxsaq2')
      call getmem2d(trxsg,1,ntr,1,2,'diagnosis:trxsg')
      call getmem2d(ttrace,1,ntr,1,2,'diagnosis:ttrace')
    end if
  end subroutine allocate_mod_diagnosis
!
  subroutine initdiag
#ifndef IBM
    use mpi
#else 
    include 'mpif.h'
#endif 
    implicit none
    integer :: ierr
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
!
!   dry air (unit = kg):
!
    call mpi_gather(sps1%ps, iy*jxp,mpi_real8, &
                    psa_io,  iy*jxp,mpi_real8, &
                    0,mpi_comm_world,ierr)
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
      tdini = tdini*dx*dx*d_1000*regrav
    end if

    call mpi_bcast(tdini,1,mpi_real8,0,mpi_comm_world,ierr)
!
!   water substance (unit = kg):
!
    call mpi_gather(atm1%qv,   iy*kz*jxp,mpi_real8, &
                    atm1_io%qv,iy*kz*jxp,mpi_real8, &
                    0,mpi_comm_world,ierr)
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
      tvmass = tvmass*dx*dx*d_1000*regrav
    end if

    call mpi_bcast(tvmass,1,mpi_real8,0,mpi_comm_world,ierr)
!
    call mpi_gather(atm1%qc,   iy*kz*jxp,mpi_real8,              &
                    atm1_io%qc,iy*kz*jxp,mpi_real8,              &
                    0,mpi_comm_world,ierr)
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
      tcmass = tcmass*dx*dx*d_1000*regrav
    end if

    call mpi_bcast(tcmass,1,mpi_real8,0,mpi_comm_world,ierr)

    tqini = tvmass + tcmass
!
!=======================================================================
!
    if ( myid == 0 ) write(6,99003) tdini , tqini
99003 format (' *** initial total air = ',e12.5,' kg, total water = ',  &
          e12.5,' kg in large domain.')
!
  end subroutine initdiag
!
!
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
      read (iunit) tchie
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
!        the unit from advection is converted to "kg".                c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine conadv
!
#ifndef IBM
    use mpi
#else 
    include 'mpif.h'
#endif 
!
    implicit none
!
    integer :: ierr
    real(8) , dimension(jxp) :: psa01 , psailx
    real(8) , dimension(jx) :: psa01_g , psailx_g
    real(8) , dimension(kz,jxp) :: qca01 , qcailx , qva01 , qvailx ,  &
                                   va01 , vaix
    real(8) , dimension(kz,jx) :: qca01_g , qcailx_g , qva01_g ,      &
                                  qvailx_g , va01_g , vaix_g
    real(8) , dimension(iym1,kz) :: worka , workb
    integer :: i ,  j , k
!
!----------------------------------------------------------------------
!
!   advection of dry air through the lateral boundaries:
!
!   advection through east-west boundaries:
!
    if ( myid == nproc-1 ) then
      do k = 1 , kz
        do i = 1 , iym1
          worka(i,k) = (atm1%u(i+1,k,jendl)+atm1%u(i,k,jendl)) / &
                       (mddom%msfx(i,jendx)*mddom%msfx(i,jendx))
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
    call mpi_bcast(worka,iym1*kz,mpi_real8,nproc-1,mpi_comm_world,ierr)
    call mpi_bcast(workb,iym1*kz,mpi_real8,0,mpi_comm_world,ierr)
    do k = 1 , kz
      do i = 1 , iym1
        tdadv = tdadv - dtsec*5.0D2*dsigma(k)*dx*(worka(i,k)-workb(i,k))*regrav
      end do
    end do
!
!   advection through north-south boundaries:
!
    do j = 1 , jendl
      do k = 1 , kz
        vaix(k,j) = atm1%v(iy,k,j)
        va01(k,j) = atm1%v(1,k,j)
      end do
    end do
    call mpi_gather(vaix,  kz*jxp,mpi_real8,vaix_g,kz*jxp,mpi_real8,  &
                    0,mpi_comm_world,ierr)
    call mpi_gather(va01,  kz*jxp,mpi_real8,va01_g,kz*jxp,mpi_real8,  &
                    0,mpi_comm_world,ierr)
    if ( myid == 0 ) then
      do k = 1 , kz
        do j = 1 , jxm1
          tdadv = tdadv - dtsec*5.0D2*dsigma(k)*dx*             &
               ((vaix_g(k,j+1)+vaix_g(k,j)) /                   &
                (mddom_io%msfx(iym1,j)*mddom_io%msfx(iym1,j)) - &
                (va01_g(k,j+1)+va01_g(k,j)) /                   &
                (mddom_io%msfx(1,j)*mddom_io%msfx(1,j)))*regrav
        end do
      end do
    end if
    call mpi_bcast(tdadv,1,mpi_real8,0,mpi_comm_world,ierr)
!
!----------------------------------------------------------------------
!
!   advection of water vapor through the lateral boundaries:
!
!   advection through east-west boundaries:
!
    if ( myid == nproc-1 ) then
      do k = 1 , kz
        do i = 1 , iym1
          worka(i,k) = (atm1%u(i+1,k,jendl)+atm1%u(i,k,jendl)) * &
                       (atm1%qv(i,k,jendx)/sps1%ps(i,jendx)) /   &
                       (mddom%msfx(i,jendx)*mddom%msfx(i,jendx))
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
    call mpi_bcast(worka,iym1*kz,mpi_real8,nproc-1,mpi_comm_world,ierr)
    call mpi_bcast(workb,iym1*kz,mpi_real8,0,mpi_comm_world,ierr)
    do k = 1 , kz
      do i = 1 , iym1
        tqadv = tqadv - dtsec*5.0D2*dsigma(k)*dx*(worka(i,k)-workb(i,k))*regrav
      end do
    end do
!
!   advection through north-south boundaries:
!
    do j = 1 , jendl
      do k = 1 , kz
        qvailx(k,j) = atm1%qv(iym1,k,j)
        qva01(k,j) = atm1%qv(1,k,j)
      end do
      psailx(j) = sps1%ps(iym1,j)
      psa01(j) = sps1%ps(1,j)
    end do
    call mpi_gather(qvailx,kz*jxp,mpi_real8,qvailx_g,kz*jxp,mpi_real8, &
                    0,mpi_comm_world,ierr)
    call mpi_gather(qva01,kz*jxp,mpi_real8,qva01_g,kz*jxp,mpi_real8,   &
                    0,mpi_comm_world,ierr)
    call mpi_gather(psailx,jxp,mpi_real8,psailx_g,jxp,mpi_real8, &
                    0,mpi_comm_world,ierr)
    call mpi_gather(psa01,jxp,mpi_real8,psa01_g,jxp,mpi_real8,   &
                    0,mpi_comm_world,ierr)
    if ( myid == 0 ) then
      do k = 1 , kz
        do j = 1 , jxm1
          tqadv = tqadv - dtsec*5.0D2*dsigma(k)*dx *                        &
                 ((vaix_g(k,j+1)+vaix_g(k,j))*(qvailx_g(k,j)/psailx_g(j)) / &
                  (mddom_io%msfx(iym1,j)*mddom_io%msfx(iym1,j)) -           &
                  (va01_g(k,j+1)+va01_g(k,j))*(qva01_g(k,j)/psa01_g(j)) /   &
                  (mddom_io%msfx(1,j)*mddom_io%msfx(1,j)))*regrav
        end do
      end do
    end if
    call mpi_bcast(tqadv,1,mpi_real8,0,mpi_comm_world,ierr)
!
!   advection of cloud water and rainwater through lateral boundaries:
!
!   advection through east-west boundaries:
!
    if ( myid == nproc-1 ) then
      do k = 1 , kz
        do i = 1 , iym1
          worka(i,k) = (atm1%u(i+1,k,jendl)+atm1%u(i,k,jendl)) * &
                       (atm1%qc(i,k,jendx)/sps1%ps(i,jendx))  /  &
                       (mddom%msfx(i,jendx)*mddom%msfx(i,jendx))
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
    call mpi_bcast(worka,iym1*kz,mpi_real8,nproc-1,mpi_comm_world,ierr)
    call mpi_bcast(workb,iym1*kz,mpi_real8,0,mpi_comm_world,ierr)
    do k = 1 , kz
      do i = 1 , iym1
        tqadv = tqadv - dtsec*5.0D2*dsigma(k)*dx*(worka(i,k)-workb(i,k))*regrav
      end do
    end do
!
!   advection through north-south boundaries:
!
    do j = 1 , jendl
      do k = 1 , kz
        qcailx(k,j) = atm1%qc(iym1,k,j)
        qca01(k,j) = atm1%qc(1,k,j)
      end do
    end do
    call mpi_gather(qcailx,kz*jxp,mpi_real8,qcailx_g,kz*jxp,mpi_real8, &
                    0,mpi_comm_world,ierr)
    call mpi_gather(qca01,kz*jxp,mpi_real8,qca01_g,kz*jxp,mpi_real8,   &
                    0,mpi_comm_world,ierr)
    if ( myid == 0 ) then
      do k = 1 , kz
        do j = 1 , jxm1
          tqadv = tqadv - dtsec*5.0D2*dsigma(k)*dx *                        &
                 ((vaix_g(k,j+1)+vaix_g(k,j))*(qcailx_g(k,j)/psailx_g(j)) / &
                  (mddom_io%msfx(iym1,j)*mddom_io%msfx(iym1,j)) -           &
                  (va01_g(k,j+1)+va01_g(k,j))*(qca01_g(k,j)/psa01_g(j)) /   &
                  (mddom_io%msfx(1,j)*mddom_io%msfx(1,j)))*regrav
        end do
      end do
    end if
    call mpi_bcast(tqadv,1,mpi_real8,0,mpi_comm_world,ierr)
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
#ifndef IBM
    use mpi
#else 
    include 'mpif.h'
#endif 
    implicit none
!
    real(8) :: error1 , error2 , tcmass , tcrai , tdrym , tncrai ,    &
               tqmass , tttmp , tvmass
    integer :: i , j , k
    integer :: ierr
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
    call mpi_gather(sps1%ps,iy*jxp,mpi_real8,psa_io,iy*jxp,mpi_real8, &
                    0,mpi_comm_world,ierr)
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
      tdrym = tdrym*dx*dx*d_1000*regrav
    end if
    call mpi_bcast(tdrym,1,mpi_real8,0,mpi_comm_world,ierr)
!
!-----water substance (unit = kg):
!
    tvmass = d_zero
    call mpi_gather(atm1%qv,iy*kz*jxp,mpi_real8,atm1_io%qv,iy*kz*jxp,mpi_real8, &
                    0,mpi_comm_world,ierr)
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
      tvmass = tvmass*dx*dx*d_1000*regrav
    end if
    call mpi_bcast(tvmass,1,mpi_real8,0,mpi_comm_world,ierr)
!
    tcmass = d_zero
    call mpi_gather(atm1%qc,iy*kz*jxp,mpi_real8,atm1_io%qc,iy*kz*jxp,mpi_real8, &
                    0,mpi_comm_world,ierr)
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
      tcmass = tcmass*dx*dx*d_1000*regrav
    end if
    call mpi_bcast(tcmass,1,mpi_real8,0,mpi_comm_world,ierr)

    tqmass = tvmass + tcmass

!=======================================================================
!
!   conservation of dry air:
!
    tdrym = tdrym - tdadv
    error1 = (tdrym-tdini)/tdini*d_100
!
!   conservation of water substance:
!
!   total raifall at this time:
!
    call mpi_gather(sfsta%rainc,iy*jxp,mpi_real8,rainc_io,iy*jxp,mpi_real8,   &
                    0,mpi_comm_world,ierr)
    call mpi_gather(sfsta%rainnc,iy*jxp,mpi_real8,rainnc_io,iy*jxp,mpi_real8, &
                    0,mpi_comm_world,ierr)
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
!
    tqmass = tqmass + tqrai - tqeva - tqadv
    error2 = (tqmass-tqini)/tqini*d_100
!
!   print out the information:
!
    if ( myid == 0 ) then
      if ( debug_level > 3 .and. mod(ntime,ndbgfrq) == 0 ) then
        write(6,*)  '***** ' , idatex%tostring(), ' *****'
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
    end if

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
#ifndef IBM
    use mpi
#else 
    include 'mpif.h'
#endif 
    implicit none
!
    real(8) , dimension(iy,kz,jxp) :: xkc
    intent (in) xkc
!
    integer :: ierr
    real(8) , dimension(kz,ntr,jxp) :: chia01 , chia02 , chiaill , chiaill1
    real(8) , dimension(kz,ntr,jx) :: chia01_g , chia02_g , chiaill1_g , chiaill_g
    real(8) , dimension(jxp) :: psa01 , psa02 , psaill , psaill1
    real(8) , dimension(jx) :: psa01_g , psa02_g , psaill1_g , psaill_g
    real(8) , dimension(kz,jxp) :: va02 , vaill , xkc02 , xkcill1
    real(8) , dimension(kz,jx) :: va02_g , vaill_g , xkc02_g , xkcill1_g
    real(8) :: chid1 , chid2
    real(8) :: fact1 , fact2 , fx1 , fx2 , uavg1 , uavg2 , vavg1 , vavg2
    integer :: i , j , k , n
    real(8) , dimension(iym1,kz,ntr) :: worka , workb
!
!
!   1  ADVECTION budgets
!
!   advection of tracer through lateral boundaries
!
!   advection through east-west boundaries: the 'relaxed' upstream scheme
!
    fact1 = 0.6D0
    fact2 = d_one - fact1
   
!   inflow/outflow
    do n = 1 , ntr
      do k = 1 , kz
        do i = 2 , iym2
          if ( myid == nproc-1 ) then
            uavg2 = (atm1%u(i+1,k,jendx)+atm1%u(i,k,jendx))*d_half
            if ( uavg2 < d_zero ) then
              worka(i,k,n) =  &
                  -uavg2*(fact1*chia(i,k,jendx,n)/sps1%ps(i,jendx)/ &
                  (mddom%msfx(i,jendx)*mddom%msfx(i,jendx))+        &
                  fact2*chia(i,k,jendm,n)/sps1%ps(i,jendm)/         &
                  (mddom%msfx(i,jendm)*mddom%msfx(i,jendm)))
            else
              worka(i,k,n) = &
                  -uavg2*(fact1*chia(i,k,jendm,n)/sps1%ps(i,jendm)/ &
                  (mddom%msfx(i,jendm)*mddom%msfx(i,jendm))+        &
                  fact2*chia(i,k,jendx,n)/sps1%ps(i,jendx) /        &
                  (mddom%msfx(i,jendx)*mddom%msfx(i,jendx)))
            end if
          end if
          if ( myid == 0 ) then
            uavg1 = (atm1%u(i+1,k,2)+atm1%u(i,k,2))*d_half
            if ( uavg1 > d_zero ) then
              workb(i,k,n) = &
                  -uavg1*(fact1*chia(i,k,1,n)/sps1%ps(i,1)/ &
                  (mddom%msfx(i,1)*mddom%msfx(i,1)) +       &
                  fact2*chia(i,k,2,n)/sps1%ps(i,2) /        &
                  (mddom%msfx(i,2)*mddom%msfx(i,2)))
            else
              workb(i,k,n) = & 
                  -uavg1*(fact1*chia(i,k,2,n)/sps1%ps(i,2) / &
                  (mddom%msfx(i,2)*mddom%msfx(i,2)) +        &
                  fact2*chia(i,k,1,n)/sps1%ps(i,1) /         &
                  (mddom%msfx(i,1)*mddom%msfx(i,1)))
            end if
          end if
        end do
      end do
    end do
    call mpi_bcast(worka,iym1*kz*ntr,mpi_real8,nproc-1,mpi_comm_world,ierr)

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
    call mpi_gather(vaill,  kz*jxp,mpi_real8,                      &
                    vaill_g,kz*jxp,mpi_real8,0,mpi_comm_world,ierr)
    call mpi_gather(va02,  kz*jxp,mpi_real8,                       &
                    va02_g,kz*jxp,mpi_real8,0,mpi_comm_world,ierr)
    call mpi_gather(xkcill1,  kz*jxp,mpi_real8,                    &
                    xkcill1_g,kz*jxp,mpi_real8,0,mpi_comm_world,ierr)
    call mpi_gather(xkc02,  kz*jxp,mpi_real8,                      &
                    xkc02_g,kz*jxp,mpi_real8,0,mpi_comm_world,ierr)
    call mpi_gather(chiaill,  kz*ntr*jxp,mpi_real8,                &
                    chiaill_g,kz*ntr*jxp,mpi_real8,0,mpi_comm_world,ierr)
    call mpi_gather(chiaill1,  kz*ntr*jxp,mpi_real8,               &
                    chiaill1_g,kz*ntr*jxp,mpi_real8,0,mpi_comm_world,ierr)
    call mpi_gather(chia01,  kz*ntr*jxp,mpi_real8,                 &
                    chia01_g,kz*ntr*jxp,mpi_real8,0,mpi_comm_world,ierr)
    call mpi_gather(chia02,  kz*ntr*jxp,mpi_real8,                 &
                    chia02_g,kz*ntr*jxp,mpi_real8,0,mpi_comm_world,ierr)
    call mpi_gather(psaill,  jxp,mpi_real8,                        &
                    psaill_g,jxp,mpi_real8,0,mpi_comm_world,ierr)
    call mpi_gather(psaill1,  jxp,mpi_real8,                       &
                    psaill1_g,jxp,mpi_real8,0,mpi_comm_world,ierr)
    call mpi_gather(psa01,  jxp,mpi_real8,                         &
                    psa01_g,jxp,mpi_real8,0,mpi_comm_world,ierr)
    call mpi_gather(psa02,  jxp,mpi_real8,                         &
                    psa02_g,jxp,mpi_real8,0,mpi_comm_world,ierr)
    if ( myid == 0 ) then
      do n = 1 , ntr
        do k = 1 , kz
          do i = 2 , iym2
            tchiad(n) = tchiad(n) + dtsec*5.0D2*dsigma(k)*dx * &
                        (worka(i,k,n)-workb(i,k,n))*regrav
          end do
        end do
!
!       advection through north-south boundaries:
!
        do k = 1 , kz
          do j = 2 , jxm2
!hy         inflow/outflow
            vavg2 = (vaill_g(k,j+1)+vaill_g(k,j))*d_half
            if ( vavg2 < d_zero ) then
              fx2 = -vavg2*(fact1*chiaill_g(k,n,j)/psaill_g(j) /      &
                     (mddom_io%msfx(iym1,j)*mddom_io%msfx(iym1,j)) +  &
                     fact2*chiaill1_g(k,n,j)/psaill1_g(j) /           &
                     (mddom_io%msfx(iym2,j)*mddom_io%msfx(iym2,j)))
            else
              fx2 = -vavg2*(fact1*chiaill1_g(k,n,j)/psaill1_g(j) /    &
                     (mddom_io%msfx(iym2,j)*mddom_io%msfx(iym2,j)) +  &
                     fact2*chiaill_g(k,n,j)/psaill_g(j) /             &
                     (mddom_io%msfx(iym1,j)*mddom_io%msfx(iym1,j)))
            end if
   
            vavg1 = (va02_g(k,j+1)+va02_g(k,j))*d_half
            if ( vavg1 > d_zero ) then
              fx1 = -vavg1*(fact1*chia01_g(k,n,j)/psa01_g(j) /        &
                     (mddom_io%msfx(1,j)*mddom_io%msfx(1,j)) +        &
                     fact2*chia02_g(k,n,j)/psa02_g(j) /               &
                     (mddom_io%msfx(2,j)*mddom_io%msfx(2,j)))
            else
              fx1 = -vavg1*(fact1*chia02_g(k,n,j)/psa02_g(j) /   &
                     (mddom_io%msfx(2,j)*mddom_io%msfx(2,j)) +   &
                     fact2*chia01_g(k,n,j)/psa01_g(j) /          &
                     (mddom_io%msfx(1,j)*mddom_io%msfx(1,j)))
            end if
            tchiad(n) = tchiad(n) + dtsec*5.0D2*dsigma(k)*dx*(fx2-fx1)*regrav
          end do
        end do
      end do
    end if
    call mpi_bcast(tchiad,ntr,mpi_real8,0,mpi_comm_world,ierr)
!
!   diffusion through east-west boundaries:
!
    do n = 1 , ntr
      do k = 1 , kz
        do i = 2 , iym2
          if ( myid == nproc-1 )  &
            worka(i,k,n) = xkc(i,k,jendm)*sps1%ps(i,jendm) * &
               (chia(i,k,jendm,n)/sps1%ps(i,jendm)-chia(i,k,jendx,n)/sps1%ps(i,jendx))
          if ( myid == 0 ) &
            workb(i,k,n) = xkc(i,k,2)*sps1%ps(i,2) *  &
               (chia(i,k,2,n)/sps1%ps(i,2)-chia(i,k,1,n)/sps1%ps(i,1))
        end do
      end do
    end do
    call mpi_bcast(worka,iym1*kz*ntr,mpi_real8,nproc-1,mpi_comm_world,ierr)

    if ( myid == 0 ) then
      do n = 1 , ntr
        do k = 1 , kz
          do i = 2 , iym2
            tchitb(n) = tchitb(n) - dtsec*5.0D2*dsigma(k) * &
                        (workb(i,k,n)+worka(i,k,n))*regrav
          end do
        end do
   
!       diffusion through north-south boundaries:
   
        do k = 1 , kz
          do j = 2 , jxm2
            chid1 = xkcill1_g(k,j)*psaill1_g(j)                       &
                    *(chiaill1_g(k,n,j)/psaill1_g(j)-chiaill_g(k,n,j) &
                    /psaill_g(j))
            chid2 = xkc02_g(k,j)*psa02_g(j)                           &
                    *(chia02_g(k,n,j)/psa02_g(j)-chia01_g(k,n,j)      &
                    /psa01_g(j))
            tchitb(n) = tchitb(n) - dtsec*5.0D2*dsigma(k)*(chid2+chid1)&
                        *regrav
          end do
        end do
      end do
    end if
    call mpi_bcast(tchitb,ntr,mpi_real8,0,mpi_comm_world,ierr)

  end subroutine tracdiag
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
! Tracer budget diagnostic
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine contrac
#ifndef IBM
    use mpi
#else 
    include 'mpif.h'
#endif 
    implicit none
!
    integer :: i , itr , j , k
    real(8) :: tttmp
    integer :: ierr , l
!
!   add tracers,  tremlsc, tremcvc, tremdrd are total amount mass
!   and chemical conversion term
!   removed by large scale, convective precipation and dry deposition
!   respectively
   
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
    do itr = 1 , ntr
      call mpi_gather(chia(1,1,1,itr),   iy*kz*jxp,mpi_real8,         &
                      chia_io(1,1,1,itr),iy*kz*jxp,mpi_real8,         &
                      0,mpi_comm_world,ierr)
      call mpi_gather(remlsc(1,1,1,itr),   iy*kz*jxp,mpi_real8,       &
                      remlsc_io(1,1,1,itr),iy*kz*jxp,mpi_real8,       &
                      0,mpi_comm_world,ierr)
      call mpi_gather(remcvc(1,1,1,itr),   iy*kz*jxp,mpi_real8,       &
                      remcvc_io(1,1,1,itr),iy*kz*jxp,mpi_real8,       &
                      0,mpi_comm_world,ierr)
      call mpi_gather(rxsg(1,1,1,itr),   iy*kz*jxp,mpi_real8,         &
                      rxsg_io(1,1,1,itr),iy*kz*jxp,mpi_real8,         &
                      0,mpi_comm_world,ierr)
      call mpi_gather(rxsaq1(1,1,1,itr),   iy*kz*jxp,mpi_real8,       &
                      rxsaq1_io(1,1,1,itr),iy*kz*jxp,mpi_real8,       &
                      0,mpi_comm_world,ierr)
      call mpi_gather(rxsaq2(1,1,1,itr),   iy*kz*jxp,mpi_real8,       &
                      rxsaq2_io(1,1,1,itr),iy*kz*jxp,mpi_real8,       &
                      0,mpi_comm_world,ierr)
      call mpi_gather(remdrd(1,1,itr),   iy*jxp,mpi_real8,            &
                      remdrd_io(1,1,itr),iy*jxp,mpi_real8,            &
                      0,mpi_comm_world,ierr)
      do l = 1 , mpy
        call mpi_gather(chemsrc(1,1,l,itr),   iy*jxp,mpi_real8,       &
                        chemsrc_io(1,1,l,itr),iy*jxp,mpi_real8,       &
                        0,mpi_comm_world,ierr)
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
   
!       emissions
        tttmp = d_zero
        do j = 2 , jxm2
          do i = 2 , iym2
            tttmp = tttmp + chemsrc_io(i,j,idatex%month,itr)*dtsec*dx*dx
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
   
    do itr = 1 , ntr
      ttrace(itr,1) = ttrace(itr,1)*d_1000*regrav
      tremlsc(itr,1) = tremlsc(itr,1)*d_1000*regrav
      tremcvc(itr,1) = tremcvc(itr,1)*d_1000*regrav
      tremdrd(itr,1) = tremdrd(itr,1)*d_1000*regrav
      trxsg(itr,1) = trxsg(itr,1)*d_1000*regrav
      trxsaq1(itr,1) = trxsaq1(itr,1)*d_1000*regrav
      trxsaq2(itr,1) = trxsaq2(itr,1)*d_1000*regrav
    end do

   
!   print out the information:
    if ( myid == 0 ) then
      if ( debug_level > 3 .and. mod(ntime,ndbgfrq) == 0 ) then
   
!       tracers
   
        write(6,*)  '************************************************'
        write(6,*)  ' Budgets for tracers (intergrated quantitites)'
        write(6,*)  ' day = ' , idatex%tostring() , ' *****'
        write(6,*)  '************************************************'
   
        do itr = 1 , ntr
   
!         errort= ttrace(itr,1)- &
!                 (tchie(itr) - (tchiad(itr)-tchitb(itr)+ &
!                 tremdrd(itr,1)+ tremlsc(itr,1)+tremcvc(itr,1)))
   
          write(6,*)  '****************************'
          write(6,99001) itr , tchie(itr)
          write(6,99002) itr , tchiad(itr)
          write(6,99003) itr , tchitb(itr)
          write(6,99004) itr , tremdrd(itr,1)
          write(6,99005) itr , tremlsc(itr,1)
          write(6,99006) itr , tremcvc(itr,1)
          write(6,99007) itr , ttrace(itr,1)
!         write(6,99008) itr , errort
   
        end do
   
      end if
    end if
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
#ifndef IBM
    use mpi
#else 
    include 'mpif.h'
#endif 
    implicit none
    integer :: ierr
    integer :: i , j
!
    call mpi_gather(sfsta%qfx,iy*jxp,mpi_real8,  &
                    qfx_io,   iy*jxp,mpi_real8,  &
                    0,mpi_comm_world,ierr)
    if ( myid == 0 ) then
      do j = 2 , jxm2
        do i = 2 , iym2
          tqeva = tqeva + qfx_io(i,j)*dx*dx*dtsec
        end do
      end do
    end if
    call mpi_bcast(tqeva,1,mpi_real8,0,mpi_comm_world,ierr)
  end subroutine conqeva
!
#endif
!
end module mod_diagnosis
