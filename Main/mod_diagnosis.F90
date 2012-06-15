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
  use mod_runparams
  use mod_atm_interface
  use mod_che_interface
  use mod_domain
  use mod_memutil
  use mod_mppparam
  use mod_mppio
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
  real(dp) :: tdadv , tdini , tqadv , tqeva , tqini , tqrai
!
  real(dp) , pointer , dimension(:) :: tchiad , tchie , tchitb
  real(dp) , pointer , dimension(:,:) :: tremcvc , tremdrd ,     &
              tremlsc , trxsaq1 , trxsaq2 , trxsg , ttrace
  real(dp) , pointer , dimension(:) :: psa01 , psailx , psa01_g , psailx_g , &
              psa02 , psaill , psaill1 , psa02_g , psaill_g , psaill1_g
  real(dp) , pointer , dimension(:,:) :: qca01 , qcailx , qva01 , qvailx ,  &
              va01 , vaix , va02 , vaill , xkc02 , xkcill1
  real(dp) , pointer , dimension(:,:) :: qca01_g , qcailx_g , qva01_g ,  &
              qvailx_g , va01_g , vaix_g , va02_g , vaill_g , xkc02_g , &
              xkcill1_g
  real(dp) , pointer , dimension(:,:,:) :: chia01 , chia02 , chiaill , chiaill1
  real(dp) , pointer , dimension(:,:,:) :: chia01_g , chia02_g ,  &
              chiaill1_g , chiaill_g
  real(dp) , pointer , dimension(:,:,:,:) :: gatmp4
  real(dp) , pointer , dimension(:,:,:) :: gatmp3

  contains

  subroutine allocate_mod_diagnosis
  implicit none
    if ( debug_level > 2 ) then
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
        call getmem1d(psa02,jce1,jce2,'diagnosis:psa02')
        call getmem1d(psaill,jce1,jce2,'diagnosis:psaill')
        call getmem1d(psaill1,jce1,jce2,'diagnosis:psaill1')
        call getmem1d(psa02_g,jcross1,jcross2,'diagnosis:psa02_g')
        call getmem1d(psaill_g,jcross1,jcross2,'diagnosis:psaill_g')
        call getmem1d(psaill1_g,jcross1,jcross2,'diagnosis:psaill1_g')
        call getmem2d(va02,jce1,jce2,1,kz,'diagnosis:va02')
        call getmem2d(vaill,jce1,jce2,1,kz,'diagnosis:vaill')
        call getmem2d(xkc02,jce1,jce2,1,kz,'diagnosis:xkc02')
        call getmem2d(xkcill1,jce1,jce2,1,kz,'diagnosis:xkcill1')
        call getmem2d(va02_g,jcross1,jcross2,1,kz,'diagnosis:va02_g')
        call getmem2d(vaill_g,jcross1,jcross2,1,kz,'diagnosis:vaill_g')
        call getmem2d(xkc02_g,jcross1,jcross2,1,kz,'diagnosis:xkc02_g')
        call getmem2d(xkcill1_g,jcross1,jcross2,1,kz,'diagnosis:xkcill1_g')
        call getmem3d(chia01,jce1,jce2,1,kz,1,ntr,'diagnosis:chia01')
        call getmem3d(chia02,jce1,jce2,1,kz,1,ntr,'diagnosis:chia02')
        call getmem3d(chiaill,jce1,jce2,1,kz,1,ntr,'diagnosis:chiaill')
        call getmem3d(chiaill1,jce1,jce2,1,kz,1,ntr,'diagnosis:chiaill1')
        call getmem3d(chia01_g,jcross1,jcross2,1,kz,1,ntr,'diagnosis:chia01_g')
        call getmem3d(chia02_g,jcross1,jcross2,1,kz,1,ntr,'diagnosis:chia02_g')
        call getmem3d(chiaill_g,jcross1,jcross2,1,kz,1,ntr,'diagnosis:chiaill_g')
        call getmem3d(chiaill1_g,jcross1,jcross2,1,kz,1,ntr,'diagnosis:chiaill1_g')
        call getmem4d(gatmp4,jcross1,jcross2,icross1,icross2, &
                      1,kz,1,ntr,'diagnosis:gatmp4')
        call getmem3d(gatmp3,jcross1,jcross2,icross1,icross2, &
                      1,ntr,'diagnosis:gatmp3')
      end if
      call getmem1d(psa01,jce1,jce2,'diagnosis:psa01')
      call getmem1d(psailx,jce1,jce2,'diagnosis:psailx')
      call getmem1d(psa01_g,jcross1,jcross2,'diagnosis:psa01_g')
      call getmem1d(psailx_g,jcross1,jcross2,'diagnosis:psailx_g')
      call getmem2d(qca01,jce1,jce2,1,kz,'diagnosis:qca01')
      call getmem2d(qcailx,jce1,jce2,1,kz,'diagnosis:qcailx')
      call getmem2d(qva01,jce1,jce2,1,kz,'diagnosis:qva01')
      call getmem2d(qvailx,jce1,jce2,1,kz,'diagnosis:qvailx')
      call getmem2d(va01,jde1,jde2,1,kz,'diagnosis:va01')
      call getmem2d(vaix,jde1,jde2,1,kz,'diagnosis:vaix')
      call getmem2d(qca01_g,jcross1,jcross2,1,kz,'diagnosis:qca01_g')
      call getmem2d(qcailx_g,jcross1,jcross2,1,kz,'diagnosis:qcailx_g')
      call getmem2d(qva01_g,jcross1,jcross2,1,kz,'diagnosis:qva01_g')
      call getmem2d(qvailx_g,jcross1,jcross2,1,kz,'diagnosis:qvailx_g')
      call getmem2d(va01_g,jdot1,jdot2,1,kz,'diagnosis:va01_g')
      call getmem2d(vaix_g,jdot1,jdot2,1,kz,'diagnosis:vaix_g')
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
    real(dp) :: tvmass , tcmass , tttmp
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
    call deco1_gather(sfs%psa,sfs_io%psa,jcross1,jcross2,icross1,icross2)
    if ( myid == 0 ) then
      do k = 1 , kz
        tttmp = d_zero
        do i = ice1 , ice2
          do j = jce1 , jce2
            tttmp = tttmp + sfs_io%psa(j,i)
          end do
        end do
        tdini = tdini + tttmp*dsigma(k)
      end do
      tdini = tdini*dx*dx*d_1000*regrav
    end if

    call mpi_bcast(tdini,1,mpi_real8,0,mycomm,ierr)
!
!   water substance (unit = kg):
!
    call deco1_gather(atm1%qv,atm1_io%qv,jcross1,jcross2,icross1,icross2,1,kz)
    if ( myid == 0 ) then
      do k = 1 , kz
        tttmp = d_zero
        do i = ice1 , ice2
          do j = jce1 , jce2
            tttmp = tttmp + atm1_io%qv(j,i,k)
          end do
        end do
        tvmass = tvmass + tttmp*dsigma(k)
      end do
      tvmass = tvmass*dx*dx*d_1000*regrav
    end if

    call mpi_bcast(tvmass,1,mpi_real8,0,mycomm,ierr)
!
    call deco1_gather(atm1%qc,atm1_io%qc,jcross1,jcross2,icross1,icross2,1,kz)
    if ( myid == 0 ) then
      do k = 1 , kz
        tttmp = d_zero
        do i = ice1 , ice2
          do j = jce1 , jce2
            tttmp = tttmp + atm1_io%qc(j,i,k)
          end do
        end do
        tcmass = tcmass + tttmp*dsigma(k)
      end do
      tcmass = tcmass*dx*dx*d_1000*regrav
    end if

    call mpi_bcast(tcmass,1,mpi_real8,0,mycomm,ierr)

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
    call mpi_bcast(tdini,1,mpi_real8,0,mycomm,ierr)
    call mpi_bcast(tdadv,1,mpi_real8,0,mycomm,ierr)
    call mpi_bcast(tqini,1,mpi_real8,0,mycomm,ierr)
    call mpi_bcast(tqadv,1,mpi_real8,0,mycomm,ierr)
    call mpi_bcast(tqeva,1,mpi_real8,0,mycomm,ierr)
    call mpi_bcast(tqrai,1,mpi_real8,0,mycomm,ierr)
    if ( ichem == 1 ) then
      call mpi_bcast(tchiad,ntr,mpi_real8,0,mycomm,ierr)
      call mpi_bcast(tchitb,ntr,mpi_real8,0,mycomm,ierr)
      call mpi_bcast(tchie,ntr,mpi_real8,0,mycomm,ierr)
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
    integer :: i ,  j , k
    real(dp) , dimension(nicross,kz) :: worka , workb
!
!----------------------------------------------------------------------
!
!   advection of dry air through the lateral boundaries:
!
!   advection through east-west boundaries:
!
    if ( ma%hasright ) then
      do k = 1 , kz
        do i = ice1 , ice2
          worka(i,k) = (atm1%u(jdi2,i+1,k)+atm1%u(jdi2,i,k)) / &
                       (mddom%msfx(jci2,i)*mddom%msfx(jci2,i))
        end do
      end do
    end if
    if ( ma%hasleft ) then
      do k = 1 , kz
        do i = ice1 , ice2
          workb(i,k) = (atm1%u(jdi1,i+1,k)+atm1%u(jdi1,i,k)) / &
                       (mddom%msfx(jci1,i)*mddom%msfx(jci1,i))
        end do
      end do
    end if
    call mpi_bcast(worka,nicross*kz,mpi_real8,nproc-1,mycomm,ierr)
    call mpi_bcast(workb,nicross*kz,mpi_real8,0,mycomm,ierr)
    do k = 1 , kz
      do i = ice1 , ice2
        tdadv = tdadv - dtsec*5.0D2*dsigma(k)*dx*(worka(i,k)-workb(i,k))*regrav
      end do
    end do
!
!   advection through north-south boundaries:
!
    do k = 1 , kz
      do j = jde1 , jde2
        vaix(j,k) = atm1%v(j,ide2,k)
        va01(j,k) = atm1%v(j,ide1,k)
      end do
    end do
    call deco1_gather(vaix,vaix_g,jdot1,jdot2,1,kz)
    call deco1_gather(va01,va01_g,jdot1,jdot2,1,kz)
    if ( myid == 0 ) then
      do k = 1 , kz
        do j = jce1 , jce2
          tdadv = tdadv - dtsec*5.0D2*dsigma(k)*dx*             &
               ((vaix_g(j+1,k)+vaix_g(j,k)) /                   &
                (mddom_io%msfx(j,ice2)*mddom_io%msfx(j,ice2)) - &
                (va01_g(j+1,k)+va01_g(j,k)) /                   &
                (mddom_io%msfx(j,1)*mddom_io%msfx(j,1)))*regrav
        end do
      end do
    end if
    call mpi_bcast(tdadv,1,mpi_real8,0,mycomm,ierr)
!
!----------------------------------------------------------------------
!
!   advection of water vapor through the lateral boundaries:
!
!   advection through east-west boundaries:
!
    if ( ma%hasright ) then
      do k = 1 , kz
        do i = ice1 , ice2
          worka(i,k) = (atm1%u(jde2,i+1,k)+atm1%u(jde2,i,k)) * &
                       (atm1%qv(jce2,i,k)/sfs%psa(jce2,i)) /   &
                       (mddom%msfx(jce2,i)*mddom%msfx(jce2,i))
        end do
      end do
    end if
    if ( ma%hasleft ) then
      do k = 1 , kz
        do i = ice1 , ice2
          workb(i,k) = (atm1%u(jde1,i+1,k)+atm1%u(jde1,i,k)) * &
                        (atm1%qv(jce1,i,k)/sfs%psa(jce1,i)) / &
                        (mddom%msfx(jce1,i)*mddom%msfx(jce1,i))
        end do
      end do
    end if
    call mpi_bcast(worka,nicross*kz,mpi_real8,nproc-1,mycomm,ierr)
    call mpi_bcast(workb,nicross*kz,mpi_real8,0,mycomm,ierr)
    do k = 1 , kz
      do i = ice1 , ice2
        tqadv = tqadv - dtsec*5.0D2*dsigma(k)*dx*(worka(i,k)-workb(i,k))*regrav
      end do
    end do
!
!   advection through north-south boundaries:
!
    do k = 1 , kz
      do j = jce1 , jce2
        qvailx(j,k) = atm1%qv(j,ice2,k)
        qva01(j,k) = atm1%qv(j,ice1,k)
      end do
    end do
    do j = jce1 , jce2
      psailx(j) = sfs%psa(j,ice2)
      psa01(j) = sfs%psa(j,ice1)
    end do
    call deco1_gather(qvailx,qvailx_g,jcross1,jcross2,1,kz)
    call deco1_gather(qva01,qva01_g,jcross1,jcross2,1,kz)
    call deco1_gather(psailx,psailx_g)
    call deco1_gather(psa01,psa01_g)
    if ( myid == 0 ) then
      do k = 1 , kz
        do j = jce1 , jce2
          tqadv = tqadv - dtsec*5.0D2*dsigma(k)*dx *                        &
                 ((vaix_g(j+1,k)+vaix_g(j,k))*(qvailx_g(j,k)/psailx_g(j)) / &
                  (mddom_io%msfx(j,ice2)*mddom_io%msfx(j,ice2)) -           &
                  (va01_g(j+1,k)+va01_g(j,k))*(qva01_g(j,k)/psa01_g(j)) /   &
                  (mddom_io%msfx(j,ice1)*mddom_io%msfx(j,ice1)))*regrav
        end do
      end do
    end if
    call mpi_bcast(tqadv,1,mpi_real8,0,mycomm,ierr)
!
!   advection of cloud water and rainwater through lateral boundaries:
!
!   advection through east-west boundaries:
!
    if ( ma%hasright ) then
      do k = 1 , kz
        do i = ice1 , ice2
          worka(i,k) = (atm1%u(jde2,i+1,k)+atm1%u(jde2,i,k)) * &
                       (atm1%qc(jce2,i,k)/sfs%psa(jce2,i))  /  &
                       (mddom%msfx(jce2,i)*mddom%msfx(jce2,i))
        end do
      end do
    end if
    if ( ma%hasleft ) then
      do k = 1 , kz
        do i = ice1 , ice2
          workb(i,k) = (atm1%u(jde1,i+1,k)+atm1%u(jde1,i,k))* &
                       (atm1%qc(jce1,i,k)/sfs%psa(jce1,i)) /  &
                       (mddom%msfx(jce1,i)*mddom%msfx(jce1,i))
        end do
      end do
    end if
    call mpi_bcast(worka,nicross*kz,mpi_real8,nproc-1,mycomm,ierr)
    call mpi_bcast(workb,nicross*kz,mpi_real8,0,mycomm,ierr)
    do k = 1 , kz
      do i = ice1 , ice2
        tqadv = tqadv - dtsec*5.0D2*dsigma(k)*dx*(worka(i,k)-workb(i,k))*regrav
      end do
    end do
!
!   advection through north-south boundaries:
!
    do k = 1 , kz
      do j = jce1 , jce2
        qcailx(j,k) = atm1%qc(j,ice2,k)
        qca01(j,k)  = atm1%qc(j,ice1,k)
      end do
    end do
    call deco1_gather(qcailx,qcailx_g,jcross1,jcross2,1,kz)
    call deco1_gather(qca01,qca01_g,jcross1,jcross2,1,kz)
    if ( myid == 0 ) then
      do k = 1 , kz
        do j = jce1 , jce2
          tqadv = tqadv - dtsec*5.0D2*dsigma(k)*dx *                        &
                 ((vaix_g(j+1,k)+vaix_g(j,k))*(qcailx_g(j,k)/psailx_g(j)) / &
                  (mddom_io%msfx(j,ice2)*mddom_io%msfx(j,ice2)) -           &
                  (va01_g(j+1,k)+va01_g(j,k))*(qca01_g(j,k)/psa01_g(j)) /   &
                  (mddom_io%msfx(j,ice1)*mddom_io%msfx(j,ice1)))*regrav
        end do
      end do
    end if
    call mpi_bcast(tqadv,1,mpi_real8,0,mycomm,ierr)
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
    real(dp) :: error1 , error2 , tcmass , tcrai , tdrym , tncrai ,    &
               tqmass , tttmp , tvmass
    character(len=32) :: appdat
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
    call deco1_gather(sfs%psa,sfs_io%psa,jcross1,jcross2,icross1,icross2)
    if ( myid == 0 ) then
      do k = 1 , kz
        tttmp = d_zero
        do i = ice1 , ice2
          do j = jce1 , jce2
            tttmp = tttmp + sfs_io%psa(j,i)
          end do
        end do
        tdrym = tdrym + tttmp*dsigma(k)
      end do
      tdrym = tdrym*dx*dx*d_1000*regrav
    end if
    call mpi_bcast(tdrym,1,mpi_real8,0,mycomm,ierr)
!
!-----water substance (unit = kg):
!
    tvmass = d_zero
    call deco1_gather(atm1%qv,atm1_io%qv,jcross1,jcross2,icross1,icross2,1,kz)
    if ( myid == 0 ) then
      do k = 1 , kz
        tttmp = d_zero
        do i = ice1 , ice2
          do j = jce1 , jce2
            tttmp = tttmp + atm1_io%qv(j,i,k)
          end do
        end do
        tvmass = tvmass + tttmp*dsigma(k)
      end do
      tvmass = tvmass*dx*dx*d_1000*regrav
    end if
    call mpi_bcast(tvmass,1,mpi_real8,0,mycomm,ierr)
!
    tcmass = d_zero
    call deco1_gather(atm1%qc,atm1_io%qc,jcross1,jcross2,icross1,icross2,1,kz)
    if ( myid == 0 ) then
      do k = 1 , kz
        tttmp = d_zero
        do i = ice1 , ice2
          do j = jce1 , jce2
            tttmp = tttmp + atm1_io%qc(j,i,k)
          end do
        end do
        tcmass = tcmass + tttmp*dsigma(k)
      end do
      tcmass = tcmass*dx*dx*d_1000*regrav
    end if
    call mpi_bcast(tcmass,1,mpi_real8,0,mycomm,ierr)

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
    call deco1_gather(sfs%rainc,sfs_io%rainc,jcross1,jcross2,icross1,icross2)
    call deco1_gather(sfs%rainnc,sfs_io%rainnc,jcross1,jcross2,icross1,icross2)
    if ( myid == 0 ) then
      tcrai = d_zero
      tncrai = d_zero
      do i = ice1 , ice2
        do j = jce1 , jce2
          tcrai = tcrai + sfs_io%rainc(j,i)*dxsq
          tncrai = tncrai + sfs_io%rainnc(j,i)*dxsq
        end do
      end do
      tqrai = tcrai + tncrai
    end if
    call mpi_bcast(tqrai,1,mpi_real8,0,mycomm,ierr)
!
    tqmass = tqmass + tqrai - tqeva - tqadv
    error2 = (tqmass-tqini)/tqini*d_100
!
!   print out the information:
!
    if ( myid == 0 ) then
      if ( debug_level > 3 .and. mod(ktau,kdbg) == 0 ) then
        appdat = tochar(idatex)
        write(6,*)  '***** ' , appdat, ' *****'
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
    real(dp) , pointer , dimension(:,:,:) , intent(in) :: xkc
!
    integer :: ierr
    integer :: i , j , k , n
    real(dp) :: chid1 , chid2
    real(dp) :: fact1 , fact2 , fx1 , fx2 , uavg1 , uavg2 , vavg1 , vavg2
    real(dp) , dimension(nicross,kz,ntr) :: worka , workb
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
        do i = ici1 , ici2
          if ( ma%hasright ) then
            uavg2 = (atm1%u(jdi2,i+1,k)+atm1%u(jdi2,i,k))*d_half
            if ( uavg2 < d_zero ) then
              worka(i,k,n) =  &
                  -uavg2*(fact1*chia(jce2,i,k,n)/sfs%psa(jce2,i)/ &
                  (mddom%msfx(jce2,i)*mddom%msfx(jce2,i))+        &
                  fact2*chia(jci2,i,k,n)/sfs%psa(jci2,i)/         &
                  (mddom%msfx(jci2,i)*mddom%msfx(jci2,i)))
            else
              worka(i,k,n) = &
                  -uavg2*(fact1*chia(jci2,i,k,n)/sfs%psa(jci2,i)/ &
                  (mddom%msfx(jci2,i)*mddom%msfx(jci2,i))+        &
                  fact2*chia(jce2,i,k,n)/sfs%psa(jce2,i) /        &
                  (mddom%msfx(jce2,i)*mddom%msfx(jce2,i)))
            end if
          end if
          if ( ma%hasleft ) then
            uavg1 = (atm1%u(jdi1,i+1,k)+atm1%u(jdi1,i,k))*d_half
            if ( uavg1 > d_zero ) then
              workb(i,k,n) = &
                  -uavg1*(fact1*chia(jce1,i,k,n)/sfs%psa(jce1,i)/ &
                  (mddom%msfx(jce1,i)*mddom%msfx(jce1,i)) +       &
                  fact2*chia(jci1,i,k,n)/sfs%psa(jci1,i) /        &
                  (mddom%msfx(jci1,i)*mddom%msfx(jci1,i)))
            else
              workb(i,k,n) = & 
                  -uavg1*(fact1*chia(jci1,i,k,n)/sfs%psa(jci1,i) / &
                  (mddom%msfx(jci1,i)*mddom%msfx(jci1,i)) +        &
                  fact2*chia(jce1,i,k,n)/sfs%psa(jce1,i) /         &
                  (mddom%msfx(jce1,i)*mddom%msfx(jce1,i)))
            end if
          end if
        end do
      end do
    end do
    call mpi_bcast(worka,nicross*kz*ntr,mpi_real8,nproc-1,mycomm,ierr)

    do k = 1 , kz
      do j = jce1 , jce2
        vaill(j,k) = atm1%v(j,idi2,k)
        va02(j,k)  = atm1%v(j,idi1,k)
        xkcill1(j,k) = xkc(j,ici2,k)
        xkc02(j,k)   = xkc(j,ici1,k)
        do n = 1 , ntr
          chiaill(j,k,n)  = chia(j,ice2,k,n)
          chiaill1(j,k,n) = chia(j,ici2,k,n)
          chia01(j,k,n) = chia(j,ice1,k,n)
          chia02(j,k,n) = chia(j,ici1,k,n)
        end do
      end do
    end do
    do j = jce1 , jce2
      psaill(j)  = sfs%psa(j,ice2)
      psaill1(j) = sfs%psa(j,ici2)
      psa01(j) = sfs%psa(j,ice1)
      psa02(j) = sfs%psa(j,ici1)
    end do
    call deco1_gather(vaill,vaill_g,jcross1,jcross2,1,kz)
    call deco1_gather(va02,va02_g,jcross1,jcross2,1,kz)
    call deco1_gather(xkcill1,xkcill1_g,jcross1,jcross2,1,kz)
    call deco1_gather(xkc02,xkc02_g,jcross1,jcross2,1,kz)
    call deco1_gather(chiaill,chiaill_g,jcross1,jcross2,1,kz,1,ntr)
    call deco1_gather(chiaill1,chiaill1_g,jcross1,jcross2,1,kz,1,ntr)
    call deco1_gather(chia01,chia01_g,jcross1,jcross2,1,kz,1,ntr)
    call deco1_gather(chia02,chia02_g,jcross1,jcross2,1,kz,1,ntr)
    call deco1_gather(psaill,psaill_g)
    call deco1_gather(psaill1,psaill1_g)
    call deco1_gather(psa01,psa01_g)
    call deco1_gather(psa02,psa02_g)
    if ( myid == 0 ) then
      do n = 1 , ntr
        do k = 1 , kz
          do i = iout1 , iout2
            tchiad(n) = tchiad(n) + dtsec*d_1000*dsigma(k)*dx * &
                        (worka(i,k,n)-workb(i,k,n))*regrav
          end do
        end do
!
!       advection through north-south boundaries:
!
        do k = 1 , kz
          do j = jout1 , jout2
!hy         inflow/outflow
            vavg2 = (vaill_g(j+1,k)+vaill_g(j,k))*d_half
            if ( vavg2 < d_zero ) then
              fx2 = -vavg2*(fact1*chiaill_g(j,k,n)/psaill_g(j) /      &
                     (mddom_io%msfx(j,icross2)*mddom_io%msfx(j,icross2)) +  &
                     fact2*chiaill1_g(j,k,n)/psaill1_g(j) /           &
                     (mddom_io%msfx(j,iout2)*mddom_io%msfx(j,iout2)))
            else
              fx2 = -vavg2*(fact1*chiaill1_g(j,k,n)/psaill1_g(j) /    &
                     (mddom_io%msfx(j,iout2)*mddom_io%msfx(j,iout2)) +  &
                     fact2*chiaill_g(j,k,n)/psaill_g(j) /             &
                     (mddom_io%msfx(j,icross2)*mddom_io%msfx(j,icross2)))
            end if
   
            vavg1 = (va02_g(j+1,k)+va02_g(j,k))*d_half
            if ( vavg1 > d_zero ) then
              fx1 = -vavg1*(fact1*chia01_g(j,k,n)/psa01_g(j) /        &
                     (mddom_io%msfx(j,1)*mddom_io%msfx(j,1)) +        &
                     fact2*chia02_g(j,k,n)/psa02_g(j) /               &
                     (mddom_io%msfx(j,2)*mddom_io%msfx(j,2)))
            else
              fx1 = -vavg1*(fact1*chia02_g(j,k,n)/psa02_g(j) /   &
                     (mddom_io%msfx(j,2)*mddom_io%msfx(j,2)) +   &
                     fact2*chia01_g(j,k,n)/psa01_g(j) /          &
                     (mddom_io%msfx(j,1)*mddom_io%msfx(j,1)))
            end if
            tchiad(n) = tchiad(n) + dtsec*d_1000*dsigma(k)*dx*(fx2-fx1)*regrav
          end do
        end do
      end do
    end if
    call mpi_bcast(tchiad,ntr,mpi_real8,0,mycomm,ierr)
!
!   diffusion through east-west boundaries:
!
    do n = 1 , ntr
      do k = 1 , kz
        do i = iout1 , iout2
          if ( myid == nproc-1 )  &
            worka(i,k,n) = xkc(jci2,i,k)*sfs%psa(jci2,i) * &
               (chia(jci2,i,k,n)/sfs%psa(jci2,i) - &
                chia(jce2,i,k,n)/sfs%psa(jce2,i))
          if ( myid == 0 ) &
            workb(i,k,n) = xkc(2,i,k)*sfs%psa(2,i) *  &
               (chia(2,i,k,n)/sfs%psa(2,i)-chia(1,i,k,n)/sfs%psa(1,i))
        end do
      end do
    end do
    call mpi_bcast(worka,nicross*kz*ntr,mpi_real8,nproc-1,mycomm,ierr)

    if ( myid == 0 ) then
      do n = 1 , ntr
        do k = 1 , kz
          do i = iout1 , iout2
            tchitb(n) = tchitb(n) - dtsec*d_1000*dsigma(k) * &
                        (workb(i,k,n)+worka(i,k,n))*regrav
          end do
        end do
   
!       diffusion through north-south boundaries:
   
        do k = 1 , kz
          do j = jout1 , jout2
            chid1 = xkcill1_g(j,k)*psaill1_g(j)                       &
                    *(chiaill1_g(j,k,n)/psaill1_g(j)-chiaill_g(j,k,n) &
                    /psaill_g(j))
            chid2 = xkc02_g(j,k)*psa02_g(j)                           &
                    *(chia02_g(j,k,n)/psa02_g(j)-chia01_g(j,k,n)      &
                    /psa01_g(j))
            tchitb(n) = tchitb(n) - dtsec*d_1000*dsigma(k)*(chid2+chid1)&
                        *regrav
          end do
        end do
      end do
    end if
    call mpi_bcast(tchitb,ntr,mpi_real8,0,mycomm,ierr)

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
    real(dp) :: tttmp
    character(len=32) :: appdat
    integer :: ierr
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
    call deco1_gather(chia,gatmp4,jcross1,jcross2,icross1,icross2,1,kz,1,ntr)
    if ( myid == 0 ) then
      do itr = 1 , ntr
        do k = 1 , kz
          tttmp = d_zero
          do j = jout1 , jout2
            do i = iout1 , iout2
              tttmp = tttmp + gatmp4(j,i,k,itr)
            end do
          end do
          ttrace(itr,1) = ttrace(itr,1) + tttmp*dsigma(k)
        end do
      end do
    end if
    call deco1_gather(remlsc,gatmp4,jcross1,jcross2,icross1,icross2,1,kz,1,ntr)
    if ( myid == 0 ) then
      do itr = 1 , ntr
        do k = 1 , kz
          tttmp = d_zero
          do j = jout1 , jout2
            do i = iout1 , iout2
              tttmp = tttmp + gatmp4(j,i,k,itr)
            end do
          end do
          tremlsc(itr,1) = tremlsc(itr,1) + tttmp*dsigma(k)
        end do
      end do
    end if
    call deco1_gather(remcvc,gatmp4,jcross1,jcross2,icross1,icross2,1,kz,1,ntr)
    if ( myid == 0 ) then
      do itr = 1 , ntr
        do k = 1 , kz
          tttmp = d_zero
          do j = jout1 , jout2
            do i = iout1 , iout2
              tttmp = tttmp + gatmp4(j,i,k,itr)
            end do
          end do
          tremcvc(itr,1) = tremcvc(itr,1) + tttmp*dsigma(k)
        end do
      end do
    end if
    call deco1_gather(rxsg,gatmp4,jcross1,jcross2,icross1,icross2,1,kz,1,ntr)
    if ( myid == 0 ) then
      do itr = 1 , ntr
        do k = 1 , kz
          tttmp = d_zero
          do j = jout1 , jout2
            do i = iout1 , iout2
              tttmp = tttmp + gatmp4(j,i,k,itr)
            end do
          end do
          trxsg(itr,1) = trxsg(itr,1) + tttmp*dsigma(k)
        end do
      end do
    end if
    call deco1_gather(rxsaq1,gatmp4,jcross1,jcross2,icross1,icross2,1,kz,1,ntr)
    if ( myid == 0 ) then
      do itr = 1 , ntr
        do k = 1 , kz
          tttmp = d_zero
          do j = jout1 , jout2
            do i = iout1 , iout2
              tttmp = tttmp + gatmp4(j,i,k,itr)
            end do
          end do
          trxsaq1(itr,1) = trxsaq1(itr,1) + tttmp*dsigma(k)
        end do
      end do
    end if
    call deco1_gather(rxsaq2,gatmp4,jcross1,jcross2,icross1,icross2,1,kz,1,ntr)
    if ( myid == 0 ) then
      do itr = 1 , ntr
        do k = 1 , kz
          tttmp = d_zero
          do j = jout1 , jout2
            do i = iout1 , iout2
              tttmp = tttmp + gatmp4(j,i,k,itr)
            end do
          end do
          trxsaq2(itr,1) = trxsaq2(itr,1) + tttmp*dsigma(k)
        end do
      end do
    end if
   
    if ( myid == 0 ) then
      do itr = 1 , ntr
        ttrace(itr,1) = ttrace(itr,1)*dx*dx
        tremlsc(itr,1) = tremlsc(itr,1)*dx*dx
        tremcvc(itr,1) = tremcvc(itr,1)*dx*dx
        trxsg(itr,1) = trxsg(itr,1)*dx*dx
        trxsaq1(itr,1) = trxsaq1(itr,1)*dx*dx
        trxsaq2(itr,1) = trxsaq2(itr,1)*dx*dx
      end do
    end if
   
    call deco1_gather(remdrd,gatmp3,jcross1,jcross2,icross1,icross2,1,ntr)
    if ( myid == 0 ) then
      do itr = 1 , ntr
        tttmp = d_zero
        do j = jout1 , jout2
          do i = iout1 , iout2
            tttmp = tttmp + gatmp3(j,i,itr)
          end do
        end do
        tremdrd(itr,1) = tremdrd(itr,1) + tttmp*dx*dx*dsigma(kz)
      end do
    end if
   
    call deco1_gather(chemsrc,chemsrc_io, &
                      jcross1,jcross2,icross1,icross2,1,ntr)
    if ( myid == 0 ) then
      do itr = 1 , ntr
!       emissions
        tttmp = d_zero
        do j = jout1 , jout2
          do i = iout1 , iout2
            tttmp = tttmp + chemsrc_io(j,i,itr)*dtsec*dx*dx
          end do
        end do
        tchie(itr) = tchie(itr) + tttmp
      end do
    end if

    call mpi_bcast(ttrace,ntr,mpi_real8,0,mycomm,ierr)
    call mpi_bcast(tremlsc,ntr,mpi_real8,0,mycomm,ierr)
    call mpi_bcast(tremcvc,ntr,mpi_real8,0,mycomm,ierr)
    call mpi_bcast(trxsg,ntr,mpi_real8,0,mycomm,ierr)
    call mpi_bcast(trxsaq1,ntr,mpi_real8,0,mycomm,ierr)
    call mpi_bcast(trxsaq2,ntr,mpi_real8,0,mycomm,ierr)
    call mpi_bcast(tremdrd,ntr,mpi_real8,0,mycomm,ierr)
    call mpi_bcast(tchie,ntr,mpi_real8,0,mycomm,ierr)
   
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
      if ( debug_level > 3 .and. mod(ktau,kdbg) == 0 ) then
   
!       tracers
   
        appdat = tochar(idatex)
        write(6,*)  '************************************************'
        write(6,*)  ' Budgets for tracers (intergrated quantitites)'
        write(6,*)  ' day = ' , appdat, ' *****'
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
    call deco1_gather(sfs%qfx,sfs_io%qfx,jcross1,jcross2,icross1,icross2)
    if ( myid == 0 ) then
      do j = jout1 , jout2
        do i = iout1 , iout2
          tqeva = tqeva + sfs_io%qfx(j,i)*dx*dx*dtsec
        end do
      end do
    end if
    call mpi_bcast(tqeva,1,mpi_real8,0,mycomm,ierr)
  end subroutine conqeva
!
end module mod_diagnosis
