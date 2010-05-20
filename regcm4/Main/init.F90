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
 
      subroutine init

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!     this subroutine reads in the initial and boundary conditions.   c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      use mod_dynparam
      use mod_param1 , only : dt , dt2 , dx , ibdyfrq
      use mod_param2 , only : ibltyp , ichem , icup , iemiss , ifrest , &
                   & iocnflx , ipptls , lakemod , icnt
      use mod_param3 , only : a , dsigma , r8pt
      use mod_iunits
      use mod_bats , only : ssw2da , sdeltk2d , sdelqk2d , sfracv2d ,   &
                   & sfracb2d , sfracs2d , svegfrac2d , ht1 , satbrt1 , &
                   & taf2d , tlef2d , ssw2d , srw2d , sol2d , solvd2d , &
                   & solvs2d , flw2d , tgb2d , swt2d , scv2d , gwet2d , &
                   & flwd2d , fsw2d, sabv2d , sinc2d , veg2d1 , sag2d , &
                   & sice2d , dew2d , pptnc , pptc , prca2d , prnca2d , &
                   & ircp2d , text2d , col2d , ocld2d , tg2d , veg2d ,  &
                   & emiss2d , psmn_o , t2mn_o , t2mx_o , tgmn_o ,      &
                   & tgmx_o , w10x_o , albvgl , albvgs
      use mod_pmoist
      use mod_main
      use mod_mainchem
      use mod_bdycod
      use mod_rad
      use mod_trachem
      use mod_message
      use mod_date , only : dectim , mdate , mdate0 , mmrec , ldatez ,  &
                   & idate1 , lyear , lmonth , lday , lhour , ndate0 ,  &
                   & ndate1 , nnnchk , jyear , jyear0, jyearr, ntime ,  &
                   & ktau , ktaur , xtime
      use mod_radbuf
      use mod_tmpsav
      use mod_constants , only : rgti
#ifdef DIAG
      use mod_diagnosis
#endif
#ifdef MPP1
      use mod_mppio
#ifdef CLM
      use mod_clm
      use clm_varsur , only : init_tgb , init_grid
      use mod_bats , only : sols2d , soll2d , solsd2d , solld2d ,       &
                   &        aldirs2d, aldirl2d , aldifs2d , aldifl2d ,  &
                   &        coszrs2d
#endif
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
!----------------------------------------------------------------------
!-----dimension the arrays for parameterizing the sfc. variables.
!     change the variable surface parameters
!
      integer :: depth , freeze , i , ibdydiff , ibdyhr0 , nxxx , nyyy ,&
               & ibdyhr1 , ibin , ilake , im1h , ip1h , ist , jlake ,   &
               & itr , j , jm1h , jp1h , k , kzzz , n 
      real(8) :: eta , hg1 , hg2 , hg3 , hg4 , hgmax , hi , hii , hs ,  &
               & tlp , ts00
      real(4) , dimension(iy,jx) :: io2d
      character(256) :: finm
#ifdef MPP1
      real(8) , dimension(iy,jxp) :: psdot
      integer :: allrec , ierr , l
#else
      real(8) , dimension(iy,jx) :: psdot
#endif
      real(8) , dimension(400) :: tlake
      logical :: existing
!
#ifdef DIAG
      real(8) :: tvmass , tcmass , tttmp
#endif
!
! absnxt  - Nearest layer absorptivities
! abstot  - Non-adjacent layer absorptivites
! emstot  - Total emissivity
!
      data tlp , ts00/50.D0 , 288.D0/
!
      existing = .false.
#ifdef MPP1
      peb  = 0.0
      pwb  = 0.0
      pebt = 0.0
      pwbt = 0.0
      teb  = 0.0
      twb  = 0.0
      tebt = 0.0
      twbt = 0.0
      qeb  = 0.0
      qwb  = 0.0
      qebt = 0.0
      qwbt = 0.0
      ueb  = 0.0
      uwb  = 0.0
      uebt = 0.0
      uwbt = 0.0
      veb  = 0.0
      vwb  = 0.0
      vebt = 0.0
      vwbt = 0.0
      pnb  = 0.0
      pss  = 0.0
      pnbt = 0.0
      psbt = 0.0
      tnb  = 0.0
      tsb  = 0.0
      tnbt = 0.0
      tsbt = 0.0
      qnb  = 0.0
      qsb  = 0.0
      qnbt = 0.0
      qsbt = 0.0
      unb  = 0.0
      usb  = 0.0
      unbt = 0.0
      usbt = 0.0
      vnb  = 0.0
      vsb  = 0.0
      vnbt = 0.0
      vsbt = 0.0
#endif
      tgmx_o = -1.E30
      t2mx_o = -1.E30
      tgmn_o =  1.E30
      t2mn_o =  1.E30
      w10x_o = -1.E30
      psmn_o =  1.E30

#ifdef MPP1
      if ( myid.eq.0 ) then
        write (finm,99001) trim(dirglob),pthsep,trim(domname),'_ICBC',  &
             & ((ndate0/10000)*100+1)*100
        inquire (file=finm,exist=existing)
        if (.not.existing) then
          write (aline,*) 'The following ICBC File does not exist: ' ,  &
              &            trim(finm), 'please check location'
          call say
          call fatal(__FILE__,__LINE__, 'ICBC FILE NOT FOUND')
        else
          open (iutbc,file=finm,form='unformatted',status='old',    &
          & access='direct',recl=iy*jx*ibyte)
        endif  
        mmrec = 0
      end if
#else
      write (finm,99001) trim(dirglob),pthsep,trim(domname),'_ICBC',    &
             & ((ndate0/10000)*100+1)*100
      inquire(file=finm,exist=existing)
        if (.not.existing) then
          write (aline,*) 'The following IBC File does not exist: ' ,   &
              &            trim(finm), 'please check location'
          call say
          call fatal(__FILE__,__LINE__,' ICBC FILE NOT FOUND')
        else 
           open (iutbc,file=finm,form='unformatted',status='old',       &
           &access='direct',recl=iy*jx*ibyte)
           mmrec = 0
        endif 
#endif
!
      if ( .not.ifrest ) then
!-----for initial run--not using restart
!
!------set rainwater and cloud water equal to zero initially.
!
        qca = 0.0
        qcb = 0.0
!
#ifdef DIAG
        tdini = 0.
        tdadv = 0.
        tqini = 0.
        tqadv = 0.
        tqeva = 0.
        tqrai = 0.
#endif
!
!chem2
        if ( ichem.eq.1 ) then
!-----    total tracer concs (initial, emission, advected)
#ifdef DIAG
          do itr = 1 , ntr
            ttrace(itr,1) = 0.
            ttrace(itr,2) = 0.
            tchie(itr) = 0.
            tchiad(itr) = 0.
            tchitb(itr) = 0.
          end do
#endif
!qhy      tchie, tchitb(replace tchidp:deposition)
!         initialize removal terms
!
          remlsc = 0.0
          remcvc = 0.0
          rxsg   = 0.0
          rxsaq1 = 0.0
          rxsaq2 = 0.0
          remdrd = 0.0
          wdlsc  = 0.0
        end if
!chem2_
!------set the variables related to blackadar pbl equal to 0 initially.
!
        if ( ibltyp.ne.0 ) then
          hfx = 0.0
          qfx = 0.0
        end if
!
        if ( icup.eq.1 ) then
          rsheat = 0.0
          rswat  = 0.0
        end if
!
!------read in the initial conditions for large domain:
!       the initial conditions are the output from PREPROC/ICBC.
!
#ifdef MPP1
#ifdef CLM
        if ( .not. allocated(init_tgb) ) allocate(init_tgb(iy,jx))
#endif
        do
          if ( myid.eq.0 ) then
            mmrec = mmrec + 1
            read (iutbc,rec=mmrec) ndate0 , nxxx , nyyy , kzzz
            if ( nyyy.ne.iy .or. nxxx.ne.jx .or. kzzz.ne.kz ) then
              write (aline,*) 'SET IN regcm.param: IY=' , iy , ' JX=' , &
                            & jx , ' KX=' , kz
              call say
              write (aline,*) 'SET IN ICBC: NY=' , nyyy , ' NX=' ,      &
                            & nxxx , ' NZ=' , kzzz
              call fatal(__FILE__,__LINE__,                             &
                        &'IMPROPER DIMENSION SPECIFICATION')
            end if
            print * , 'READING INITAL CONDITIONS' , ndate0
            if ( ndate0.lt.mdatez(nnnchk) ) then
              print * , ndate0 , mdatez(nnnchk) , nnnchk
              print * , 'read in datasets at :' , ndate0
              if ( ehso4 ) then
                if ( lsmtyp.ne.'USGS' ) then
                  mmrec = mmrec + kz*5 + 2
                else
                  mmrec = mmrec + kz*5 + 2 + 13
                end if
              else if ( lsmtyp.ne.'USGS' ) then
                mmrec = mmrec + kz*4 + 2
              else
                mmrec = mmrec + kz*4 + 2 + 13
              end if
              print * , 'Searching for proper date: ' , ndate1 ,        &
                  & mdatez(nnnchk+1)
              print * , ndate0 , mdatez(nnnchk)
              cycle ! Proper date not found
            else if ( ndate0.gt.mdatez(nnnchk) ) then
              write (aline,*) ndate0 , mdatez(nnnchk)
              call say
              call fatal(__FILE__,__LINE__,                             &
                        &'DATE IN ICBC FILE EXCEEDED DATE IN RegCM')
            else
            end if
          end if
          exit ! Found proper date
        end do
!
        if ( myid.eq.0 ) then
          print * , 'U'
          do k = kz , 1 , -1
            mmrec = mmrec + 1
            read (iutbc,rec=mmrec) ((io2d(i,j),j=1,jx),i=1,iy)
            do j = 1 , jx
              do i = 1 , iy
                ub0_io(i,k,j) = dble(io2d(i,j))
              end do
            end do
          end do
          print * , 'V'
          do k = kz , 1 , -1
            mmrec = mmrec + 1
            read (iutbc,rec=mmrec) ((io2d(i,j),j=1,jx),i=1,iy)
            do j = 1 , jx
              do i = 1 , iy
                vb0_io(i,k,j) = dble(io2d(i,j))
              end do
            end do
          end do
          print * , 'TA'
          do k = kz , 1 , -1
            mmrec = mmrec + 1
            read (iutbc,rec=mmrec) ((io2d(i,j),j=1,jx),i=1,iy)
            do j = 1 , jx
              do i = 1 , iy
                tb0_io(i,k,j) = dble(io2d(i,j))
              end do
            end do
          end do
          print * , 'QV'
          do k = kz , 1 , -1
            mmrec = mmrec + 1
            read (iutbc,rec=mmrec) ((io2d(i,j),j=1,jx),i=1,iy)
            do j = 1 , jx
              do i = 1 , iy
                qb0_io(i,k,j) = dble(io2d(i,j))
              end do
            end do
          end do
          print * , 'PS'
          mmrec = mmrec + 1
          read (iutbc,rec=mmrec) ((io2d(i,j),j=1,jx),i=1,iy)
          do j = 1 , jx
            do i = 1 , iy
              ps0_io(i,j) = dble(io2d(i,j))
            end do
          end do
          print * , 'TS'
          mmrec = mmrec + 1
          read (iutbc,rec=mmrec) ((io2d(i,j),j=1,jx),i=1,iy)
          do j = 1 , jx
            do i = 1 , iy
              ts0_io(i,j) = dble(io2d(i,j))
            end do
          end do
          if ( ehso4 ) then
            print * , 'SO0'
            do k = kz , 1 , -1
              mmrec = mmrec + 1
              read (iutbc,rec=mmrec) ((io2d(i,j),j=1,jx),i=1,iy)
              do j = 1 , jx
                do i = 1 , iy
                  so0_io(i,k,j) = dble(io2d(i,j))
                end do
              end do
            end do
          end if
          if ( lsmtyp.eq.'USGS' ) mmrec = mmrec + 13

          do j = 1 , jx
            do k = 1 , kz
              do i = 1 , iy
                sav_0(i,k,j) = ub0_io(i,k,j)
                sav_0(i,kz+k,j) = vb0_io(i,k,j)
                sav_0(i,kz*2+k,j) = qb0_io(i,k,j)
                sav_0(i,kz*3+k,j) = tb0_io(i,k,j)
              end do
            end do
            do i = 1 , iy
              sav_0(i,kz*4+1,j) = ps0_io(i,j)
              sav_0(i,kz*4+2,j) = ts0_io(i,j)
            end do
          end do
          if ( ehso4 ) then
            do j = 1 , jx
              do k = 1 , kz
                do i = 1 , iy
                  sav_0s(i,k,j) = so0_io(i,k,j)
                end do
              end do
            end do
          end if
        end if ! end if myid=0
!
!       Start transmission of data to other processors
!
        call mpi_bcast(ndate0,1,mpi_integer,0,mpi_comm_world,ierr)
        call mpi_bcast(nxxx,1,mpi_integer,0,mpi_comm_world,ierr)
        call mpi_bcast(nyyy,1,mpi_integer,0,mpi_comm_world,ierr)
        call mpi_bcast(kzzz,1,mpi_integer,0,mpi_comm_world,ierr)
        call mpi_scatter(sav_0(1,1,1),iy*(kz*4+2)*jxp,mpi_real8,        &
                       & sav0(1,1,1), iy*(kz*4+2)*jxp,mpi_real8,        &
                       & 0,mpi_comm_world,ierr)
        if ( ehso4 )                                                    &
          & call mpi_scatter(sav_0s(1,1,1),iy*kz*jxp,mpi_real8,         &
          &                  sav0s(1,1,1), iy*kz*jxp,mpi_real8,         &
          &                  0,mpi_comm_world,ierr)
        do j = 1 , jendl
          do k = 1 , kz
            do i = 1 , iy
              ub0(i,k,j) = sav0(i,k,j)
              vb0(i,k,j) = sav0(i,kz+k,j)
              qb0(i,k,j) = sav0(i,kz*2+k,j)
              tb0(i,k,j) = sav0(i,kz*3+k,j)
            end do
          end do
          do i = 1 , iy
            ps0(i,j) = sav0(i,kz*4+1,j)
            ts0(i,j) = sav0(i,kz*4+2,j)
          end do
          if ( ehso4 ) then
            do k = 1 , kz
              do i = 1 , iy
                so0(i,k,j) = sav0s(i,k,j)
              end do
            end do
          end if
        end do
!
!       Convert surface pressure to pstar
!
        do j = 1 , jendl
          do i = 1 , iy
            ps0(i,j) = ps0(i,j) - r8pt
          end do
        end do
!=======================================================================
!
!       this routine determines p(.) from p(x) by a 4-point
!       interpolation. on the x-grid, a p(x) point outside the grid
!       domain is assumed to satisfy p(0,j)=p(1,j); p(iy,j)=p(iym1,j);
!       and similarly for the i's.
!
        call mpi_sendrecv(ps0(1,jxp),iy,mpi_real8,ieast,1,              &
                        & ps0(1,0),iy,mpi_real8,iwest,1,                &
                        & mpi_comm_world,mpi_status_ignore,ierr)
        do j = jbegin , jendx
          do i = 2 , iym1
            psdot(i,j) = 0.25*(ps0(i,j)   + ps0(i-1,j) +                &
                         &     ps0(i,j-1) + ps0(i-1,j-1))
          end do
        end do
!
        do i = 2 , iym1
          if ( myid.eq.0 ) psdot(i,1) = 0.5*(ps0(i,1)+ps0(i-1,1))
          if ( myid.eq.nproc-1 ) psdot(i,jendl)                         &
             & = 0.5*(ps0(i,jendx)+ps0(i-1,jendx))
        end do
!
        do j = jbegin , jendx
          psdot(1,j) = 0.5*(ps0(1,j)+ps0(1,j-1))
          psdot(iy,j) = 0.5*(ps0(iym1,j)+ps0(iym1,j-1))
        end do
!
        if ( myid.eq.0 ) then
          psdot(1,1) = ps0(1,1)
          psdot(iy,1) = ps0(iym1,1)
        end if
        if ( myid.eq.nproc-1 ) then
          psdot(1,jendl) = ps0(1,jendx)
          psdot(iy,jendl) = ps0(iym1,jendx)
        end if
!
!=======================================================================
!       Couple pressure u,v,t,q
        do k = 1 , kz
          do j = 1 , jendl
            do i = 1 , iy
              ub0(i,k,j) = ub0(i,k,j)*psdot(i,j)
              vb0(i,k,j) = vb0(i,k,j)*psdot(i,j)
              tb0(i,k,j) = tb0(i,k,j)*ps0(i,j)
              qb0(i,k,j) = qb0(i,k,j)*ps0(i,j)
            end do
         end do
        end do
        if ( myid.eq.0 ) then
          mmrec = mmrec + 1
          read (iutbc,rec=mmrec) ndate1
          mmrec = mmrec - 1
          ibdyhr0 = ndate0 - (ndate0/100)*100
          ibdyhr1 = ndate1 - (ndate1/100)*100
          if ( ibdyhr1.eq.0 ) ibdyhr1 = 24
          ibdydiff = ibdyhr1 - ibdyhr0
          if ( ibdydiff.ne.ibdyfrq ) then
            write (aline,*) '  ndate0=' , ndate0 , 'ndate1=' , ndate1
            call say
            write (aline,*) '  ibdyfrq=' , ibdyfrq , 'ibdydiff=' ,    &
                          & ibdydiff
            call say
            write (aline,*) '  ibdyhr0=' , ibdyhr0 , 'ibdyhr1=' ,     &
                          & ibdyhr1
            call say
            call fatal(__FILE__,__LINE__,                             &
                      &'BOUNDARY CONDITION FREQUENCY INCOMPATIBILITY')
          end if
        end if
!
        mdate = ndate0
!
!       Initialize variables and convert to double precision
!
        do k = 1 , kz
          do j = 1 , jendl
            do i = 1 , iy
              ua(i,k,j) = ub0(i,k,j)
              ub(i,k,j) = ub0(i,k,j)
              va(i,k,j) = vb0(i,k,j)
              vb(i,k,j) = vb0(i,k,j)
              qva(i,k,j) = qb0(i,k,j)
              qvb(i,k,j) = qb0(i,k,j)
              ta(i,k,j) = tb0(i,k,j)
              tb(i,k,j) = tb0(i,k,j)
            end do
          end do
        end do
        do j = 1 , jendl
          do i = 1 , iy
            psa(i,j) = ps0(i,j)
            psb(i,j) = ps0(i,j)
            tga(i,j) = ts0(i,j)
            tgb(i,j) = ts0(i,j)
          end do
        end do
#ifdef SEAICE
        do j = 1 , jendx
          do i = 1 , iym1
            if ( veg2d(i,j).le.0.00001 ) then
              if ( ts0(i,j).le.271.38 ) then
                tga(i,j) = 271.38
                tgb(i,j) = 271.38
                ts0(i,j) = 271.38
!               write(*,*) 'Sea Ice point:', i,j
                do n = 1, nnsg
                  ocld2d(n,i,j)=2.
                  sice2d(n,i,j)=1000.
                end do
              else
                do n = 1, nnsg
                  ocld2d(n,i,j)=0.
                  sice2d(n,i,j)=0.
                end do
              end if
            end if
          end do
        end do
#endif
        do k = 1 , kz
          do j = 1 , jendl
            do i = 1 , iy
              tbase(i,k,j) = ts00 + tlp*dlog((psa(i,j)*a(k)+r8pt)/100.)
            end do
          end do
        end do
        if ( ehso4 ) then
          do k = 1 , kz
            do j = 1 , jendl
              do i = 1 , iy
                so4(i,k,j) = so0(i,k,j)
              end do
            end do
          end do
        end if
!
        do j = 1 , jendx
          do i = 1 , iym1
            tga(i,j) = ta(i,kz,j)/psa(i,j)
            tgb(i,j) = tb(i,kz,j)/psb(i,j)
            tgbb(i,j) = tb(i,kz,j)/psb(i,j)
            zpbl(i,j) = 500.  ! For Zeng Ocean Flux Scheme
          end do
        end do
        do j = 1 , jendx
          do i = 1 , iym1
            do k = 1 , nnsg
              snowc(k,i,j) = 0.
            end do
          end do
        end do
        if ( ichem.eq.1 ) then
          ssw2da    = 0.0
          sdeltk2d  = 0.0
          sdelqk2d  = 0.0
          sfracv2d  = 0.5
          sfracb2d  = 0.5
          sfracs2d  = 0.0
          svegfrac2d = 0.0
        end if
#else
        do
          mmrec = mmrec + 1
          read (iutbc,rec=mmrec) ndate0 , nxxx , nyyy , kzzz
          if ( nyyy.ne.iy .or. nxxx.ne.jx .or. kzzz.ne.kz ) then
            write (aline,*) 'SET IN regcm.param: IY=' , iy , ' JX=' ,   &
                          & jx , ' KX=' , kz
            call say
            write (aline,*) 'SET IN ICBC: NY=' , nyyy , ' NX=' , nxxx , &
                           &' NZ=' , kzzz
            call fatal(__FILE__,__LINE__,                               &
                      &'IMPROPER DIMENSION SPECIFICATION')
          end if
          print * , 'READING INITAL CONDITIONS' , ndate0
          if ( ndate0.lt.mdatez(nnnchk) ) then
            print * , ndate0 , mdatez(nnnchk) , nnnchk
            print * , 'read in datasets at :' , ndate0
            if ( ehso4 ) then
              if ( lsmtyp.ne.'USGS' ) then
                mmrec = mmrec + kz*5 + 2
              else
                mmrec = mmrec + kz*5 + 2 + 13
              end if
            else if ( lsmtyp.ne.'USGS' ) then
              mmrec = mmrec + kz*4 + 2
            else
              mmrec = mmrec + kz*4 + 2 + 13
            end if
            print * , 'Searching for proper date: ' , ndate1 ,          &
                & mdatez(nnnchk+1)
            print * , ndate0 , mdatez(nnnchk)
            cycle ! Proper date still not found
          else if ( ndate0.gt.mdatez(nnnchk) ) then
            write (aline,*) ndate0 , mdatez(nnnchk)
            call say
            call fatal(__FILE__,__LINE__,                               &
                      &'DATE IN ICBC FILE EXCEEDED DATE IN RegCM')
          else
          end if
          exit ! Found proper date
        end do
!
        print * , 'U'
        do k = kz , 1 , -1
          mmrec = mmrec + 1
          read (iutbc,rec=mmrec) ((io2d(i,j),j=1,jx),i=1,iy)
          do j = 1 , jx
            do i = 1 , iy
              ub0(i,k,j) = dble(io2d(i,j))
            end do
          end do
        end do
        print * , 'V'
        do k = kz , 1 , -1
          mmrec = mmrec + 1
          read (iutbc,rec=mmrec) ((io2d(i,j),j=1,jx),i=1,iy)
          do j = 1 , jx
            do i = 1 , iy
              vb0(i,k,j) = dble(io2d(i,j))
            end do
          end do
        end do
        print * , 'TA'
        do k = kz , 1 , -1
          mmrec = mmrec + 1
          read (iutbc,rec=mmrec) ((io2d(i,j),j=1,jx),i=1,iy)
          do j = 1 , jx
            do i = 1 , iy
              tb0(i,k,j) = dble(io2d(i,j))
            end do
          end do
        end do
        print * , 'QV'
        do k = kz , 1 , -1
          mmrec = mmrec + 1
          read (iutbc,rec=mmrec) ((io2d(i,j),j=1,jx),i=1,iy)
          do j = 1 , jx
            do i = 1 , iy
              qb0(i,k,j) = dble(io2d(i,j))
            end do
          end do
        end do
        print * , 'PS'
        mmrec = mmrec + 1
        read (iutbc,rec=mmrec) ((io2d(i,j),j=1,jx),i=1,iy)
        do j = 1 , jx
          do i = 1 , iy
            ps0(i,j) = dble(io2d(i,j))
          end do
        end do
        print * , 'TS'
        mmrec = mmrec + 1
        read (iutbc,rec=mmrec) ((io2d(i,j),j=1,jx),i=1,iy)
        do j = 1 , jx
          do i = 1 , iy
            ts0(i,j) = dble(io2d(i,j))
          end do
        end do
        if ( ehso4 ) then
          print * , 'SO0'
          do k = kz , 1 , -1
            mmrec = mmrec + 1
            read (iutbc,rec=mmrec) ((io2d(i,j),j=1,jx),i=1,iy)
            do j = 1 , jx
              do i = 1 , iy
                so0(i,k,j) = dble(io2d(i,j))
              end do
            end do
          end do
        end if
        if ( lsmtyp.eq.'USGS' ) mmrec = mmrec + 13
!
!       Convert surface pressure to pstar
!
        do j = 1 , jx
          do i = 1 , iy
            ps0(i,j) = ps0(i,j) - r8pt
          end do
        end do
!=======================================================================
!
!       this routine determines p(.) from p(x) by a 4-point
!       interpolation. on the x-grid, a p(x) point outside the grid
!       domain is assumed to satisfy p(0,j)=p(1,j);
!       p(iy,j)=p(iym1,j); and similarly for the i's.
!
        do j = 2 , jxm1
          do i = 2 , iym1
            psdot(i,j) = 0.25*(ps0(i,j)+ps0(i-1,j)+                     &
                       &       ps0(i,j-1)+ps0(i-1,j-1))
          end do
        end do
!
        do i = 2 , iym1
          psdot(i,1) = 0.5*(ps0(i,1)+ps0(i-1,1))
          psdot(i,jx) = 0.5*(ps0(i,jxm1)+ps0(i-1,jxm1))
        end do
!
        do j = 2 , jxm1
          psdot(1,j) = 0.5*(ps0(1,j)+ps0(1,j-1))
          psdot(iy,j) = 0.5*(ps0(iym1,j)+ps0(iym1,j-1))
        end do
!
        psdot(1,1) = ps0(1,1)
        psdot(iy,1) = ps0(iym1,1)
        psdot(1,jx) = ps0(1,jxm1)
        psdot(iy,jx) = ps0(iym1,jxm1)
!
!=======================================================================
!       Couple pressure u,v,t,q
!
        do k = 1 , kz
          do j = 1 , jx
            do i = 1 , iy
              ub0(i,k,j) = ub0(i,k,j)*psdot(i,j)
              vb0(i,k,j) = vb0(i,k,j)*psdot(i,j)
              tb0(i,k,j) = tb0(i,k,j)*ps0(i,j)
              qb0(i,k,j) = qb0(i,k,j)*ps0(i,j)
            end do
          end do
        end do
        mmrec = mmrec + 1
        read (iutbc,rec=mmrec) ndate1
        mmrec = mmrec - 1
        ibdyhr0 = ndate0 - (ndate0/100)*100
        ibdyhr1 = ndate1 - (ndate1/100)*100
        if ( ibdyhr1.eq.0 ) ibdyhr1 = 24
        ibdydiff = ibdyhr1 - ibdyhr0
        if ( ibdydiff.ne.ibdyfrq ) then
          write (aline,*) '  ndate0=' , ndate0 , 'ndate1=' , ndate1
          call say
          write (aline,*) '  ibdyfrq=' , ibdyfrq , 'ibdydiff=' ,      &
                        & ibdydiff
          call say
          write (aline,*) '  ibdyhr0=' , ibdyhr0 , 'ibdyhr1=' ,       &
                        & ibdyhr1
          call say
          call fatal(__FILE__,__LINE__,                               &
                    &'BOUNDARY CONDITION FREQUENCY INCOMPATIBILITY')
        end if
!
        mdate = ndate0
!
!       Initialize variables and convert to double precision
!
        do k = 1 , kz
          do j = 1 , jx
            do i = 1 , iy
              ua(i,k,j) = ub0(i,k,j)
              ub(i,k,j) = ub0(i,k,j)
              va(i,k,j) = vb0(i,k,j)
              vb(i,k,j) = vb0(i,k,j)
              qva(i,k,j) = qb0(i,k,j)
              qvb(i,k,j) = qb0(i,k,j)
              ta(i,k,j) = tb0(i,k,j)
              tb(i,k,j) = tb0(i,k,j)
            end do
          end do
        end do
        do j = 1 , jx
          do i = 1 , iy
            psa(i,j) = ps0(i,j)
            psb(i,j) = ps0(i,j)
            tga(i,j) = ts0(i,j)
            tgb(i,j) = ts0(i,j)
          end do
        end do
#ifdef SEAICE
        do j = 1 , jxm1
          do i = 1 , iym1
            if ( veg2d(i,j).le.0.00001 ) then
              if ( ts0(i,j).le.271.38 ) then
                tga(i,j) = 271.38
                tgb(i,j) = 271.38
                ts0(i,j) = 271.38
!               write(*,*) 'Sea Ice point:', i,j
                do n = 1, nnsg
                  ocld2d(n,i,j)=2.
                  sice2d(n,i,j)=1000.
                end do
              else
                do n = 1, nnsg
                  ocld2d(n,i,j)=0.
                  sice2d(n,i,j)=0.
                end do
              end if
            end if
          end do
        end do
#endif
        do k = 1 , kz
          do j = 1 , jx
            do i = 1 , iy
              tbase(i,k,j) = ts00 + tlp*dlog((psa(i,j)*a(k)+r8pt)/100.)
            end do
          end do
        end do
        if ( ehso4 ) then
          do k = 1 , kz
            do j = 1 , jx
              do i = 1 , iy
                so4(i,k,j) = so0(i,k,j)
              end do
            end do
          end do
        end if
!
        do j = 1 , jxm1
          do i = 1 , iym1
            tga(i,j) = ta(i,kz,j)/psa(i,j)
            tgb(i,j) = tb(i,kz,j)/psb(i,j)
            tgbb(i,j) = tb(i,kz,j)/psb(i,j)
            zpbl(i,j) = 500.
                       ! For Zeng Ocean Flux Scheme
          end do
        end do
        do j = 1 , jxm1
          do i = 1 , iym1
            do k = 1 , nnsg
              snowc(k,i,j) = 0.
            end do
          end do
        end do
        if ( ichem.eq.1 ) then
          ssw2da = 0.0
          sdeltk2d = 0.0
          sdelqk2d = 0.0
          sfracv2d = 0.5
          sfracb2d = 0.5
          sfracs2d = 0.0
          svegfrac2d = 0.0
        end if
#endif

#ifdef DIAG
#ifdef MPP1
!=======================================================================
!
!-----dry air (unit = kg):
!
        tdini = 0.
        call mpi_gather(psa(1,1),   iy*jxp,mpi_real8,                   &
                      & psa_io(1,1),iy*jxp,mpi_real8,                   &
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
        tvmass = 0.
        call mpi_gather(qva(1,1,1),   iy*kz*jxp,mpi_real8,              &
                      & qva_io(1,1,1),iy*kz*jxp,mpi_real8,              &
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
        tcmass = 0.
        call mpi_gather(qca(1,1,1),   iy*kz*jxp,mpi_real8,              &
                      & qca_io(1,1,1),iy*kz*jxp,mpi_real8,              &
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
        tdini = 0.
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
        tvmass = 0.
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
        tcmass = 0.
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
#endif
!
!chem2
        if ( ichem.eq.1 ) then
!-----set tracer concs to 1 (kg/kg) initially. Must convert this to p*
!-----mixing ratio to compute tendencies:
!US       mass test zero concs init input for advection
!qhy      initial chia is 10ppt
!hy       set the initial tracer concentration 10ppt (1.e-11), 9/4/98
 
          do itr = 1 , ntr
            do k = 1 , kz
#ifdef MPP1
              do j = 1 , jendx
                do i = 1 , iym1
                  chia(i,k,j,itr) = psa(i,j)*0.0D0
                  chib(i,k,j,itr) = psb(i,j)*0.0D0
!                 chia(i,k,j,itr)=psa(i,j)*1.e-11
!                 chib(i,k,j,itr)=psb(i,j)*1.e-11
                end do
              end do
#else
              do j = 1 , jxm1
                do i = 1 , iym1
                  chia(i,k,j,itr) = psa(i,j)*0.0D0
                  chib(i,k,j,itr) = psb(i,j)*0.0D0
!                 chia(i,k,j,itr)=psa(i,j)*1.e-11
!                 chib(i,k,j,itr)=psb(i,j)*1.e-11
                end do
              end do
#endif
            end do
          end do
 
        end if
!chem2_
!
!------set rainc and rainnc equal to 0. initially
!
        rainc  = 0.0
        rainnc = 0.0
 
        if ( icup.eq.4 ) then
          cbmf2d = 0.0
        end if
!
      else ! ifrest=.true.
!-----when ifrest=.true., read in the data saved from previous run
!       for large domain from unit 14.
!
#ifdef MPP1
        if ( myid.eq.0 ) then
          write (finm,99002) trim(dirout),pthsep,'SAV.',                &
            &    ((ndate0/10000)*100+1)*100
          inquire (file=finm,exist=existing)
          if ( .not.existing ) then
            write (aline,*) 'The following SAV File does not exist: ' , &
                &            trim(finm), 'please check location'
            call say
            call fatal(__FILE__,__LINE__, 'SAV FILE NOT FOUND')
          else
            open (iutrs,file=finm,form='unformatted',status='old')
          end if
          do ! Loop while ldatez.ne.idate1
            read (iutrs) mdate0
            jyear0 = mdate0/1000000
            read (iutrs) ktau , xtime , ldatez , lyear , lmonth , lday ,&
                       & lhour , ntime
            jyear = lyear
            jyearr = jyear
            ktaur = ktau
            if ( ehso4 ) then
              read (iutrs) ub0_io , vb0_io , qb0_io , tb0_io , ps0_io , &
                         & ts0_io , so0_io
            else
              read (iutrs) ub0_io , vb0_io , qb0_io , tb0_io , ps0_io , &
                         & ts0_io
            end if
            read (iutrs) ua_io
            read (iutrs) ub_io
            read (iutrs) va_io
            read (iutrs) vb_io
            read (iutrs) ta_io
            read (iutrs) tb_io
            read (iutrs) qva_io
            read (iutrs) qvb_io
            read (iutrs) qca_io
            read (iutrs) qcb_io
            read (iutrs) psa_io , psb_io , satbrt_io , satbrt1_io , f_io
            read (iutrs) ht_io , ht1_io , msfx_io , msfd_io , xlat_io , &
                       & xlong_io
            read (iutrs) tga_io , tgb_io , rainc_io , rainnc_io
            if ( icup.eq.1 ) then
              read (iutrs) rsheat_io , rswat_io
            else if ( icup.eq.3 ) then
              read (iutrs) tbase_io , cldefi_io
            else if ( icup.eq.4 ) then
              read (iutrs) cbmf2d_io
            else
            end if
            read (iutrs) hfx_io , qfx_io , snowc_io , uvdrag_io
#ifdef DIAG
            read (iutrs) tdini , tdadv , tqini , tqadv , tqeva , tqrai
#endif
            read (iutrs) absnxt_io , abstot_io , emstot_io
            if ( ipptls.eq.1 ) read (iutrs) fcc_io
#ifdef CLM
            read (iutrs) sols2d_io
            read (iutrs) soll2d_io
            read (iutrs) solsd2d_io
            read (iutrs) solld2d_io
            read (iutrs) flwd2d_io
            read (iutrs) aldirs2d_io
            read (iutrs) aldirl2d_io
            read (iutrs) aldifs2d_io
            read (iutrs) aldifl2d_io
            read (iutrs) coszrs2d_io
            read (iutrs) ocld2d_io
            read (iutrs) heatrt_io
            read (iutrs) o3prof_io
            read (iutrs) tgbb_io
            read (iutrs) flw2d_io
            read (iutrs) swt2d_io
            read (iutrs) sinc2d_io
            read (iutrs) fsw2d_io
            read (iutrs) taf2d_io
#else
            read (iutrs) sol2d_io , solvd2d_io , solvs2d_io , flw2d_io ,&
                       & flwd2d_io , fsw2d_io , sabv2d_io , sinc2d_io
            read (iutrs) taf2d_io , tlef2d_io , tgbb_io , ssw2d_io ,    &
                       & srw2d_io , tg2d_io , tgb2d_io , swt2d_io ,     &
                       & scv2d_io , gwet2d_io , veg2d_io , veg2d1_io ,  &
                       & sag2d_io , sice2d_io , dew2d_io , ircp2d_io ,  &
                       & text2d_io , col2d_io , ocld2d_io , heatrt_io , &
                       & o3prof_io
#endif
            read (iutrs) pptnc_io , pptc_io , prca2d_io , prnca2d_io
            if ( iocnflx.eq.2 ) read (iutrs) zpbl_io
!chem2---
            if ( ichem.eq.1 ) then
              read (iutrs) chia_io
              read (iutrs) chib_io
!             cumul removal terms (3d, 2d)
              read (iutrs) remlsc_io
              read (iutrs) remcvc_io
              read (iutrs) remdrd_io
              read (iutrs) ssw2da_io
              read (iutrs) sdeltk2d_io
              read (iutrs) sdelqk2d_io
              read (iutrs) sfracv2d_io
              read (iutrs) sfracb2d_io
              read (iutrs) sfracs2d_io
              read (iutrs) svegfrac2d_io
!             cumul ad, dif, emis terms ( scalar)
#ifdef DIAG
              read (iutrs) tchiad
              read (iutrs) tchitb
              read (iutrs) tchie
#endif
            end if
 
!------lake model
            if ( lakemod.eq.1 ) then
              lcount = 0
              iin = 41
              iout = 42
              rewind (iin)
              read (iutrs) numpts
              print * , 'reading lake model restart file. numpts = ' ,  &
                  & numpts
              print * , 'jyear, ktau, xtime = ' , jyear , ktau , xtime
              do n = 1 , numpts
                read (iutrs) ilake , jlake , depth , freeze , hi , hii ,&
                           & hs , eta , (tlake(j),j=1,depth)
                print * , 'reading restart file at i, j = ' , ilake ,   &
                    & jlake
                write (iin) ilake , jlake , depth , freeze , hi , hii , &
                          & hs , eta , (tlake(j),j=1,depth)
              end do
              rewind (iin)
            end if
!
            print * , 'ozone profiles restart'
            do k = 1 , kzp1
              write (6,99004) o3prof_io(3,3,k)
            end do
            print 99005 , xtime , ktau , jyear , finm
!
            if ( ldatez.ne.idate1 ) then
              write (*,*) 'INIT: ldatez, idate1=' , ldatez , idate1
              cycle
            end if

            exit ! We have ldatez.eq.idate1
          end do
        end if
!
!       Start sending data to all processors : surface data
!
        if ( myid.eq.0 ) then
          do j = 1 , jx
            do i = 1 , iy
              inisrf_0(i,1,j) = ht_io(i,j)
              inisrf_0(i,2,j) = satbrt_io(i,j)
              inisrf_0(i,3,j) = xlat_io(i,j)
              inisrf_0(i,4,j) = xlong_io(i,j)
              inisrf_0(i,5,j) = msfx_io(i,j)
              inisrf_0(i,6,j) = msfd_io(i,j)
              inisrf_0(i,7,j) = f_io(i,j)
            end do
            do n = 1 , nnsg
              do i = 1 , iy
                inisrf_0(i,7+n,j) = ht1_io(n,i,j)
                inisrf_0(i,7+nnsg+n,j) = satbrt1_io(n,i,j)
              end do
            end do
          end do
        end if

        call mpi_scatter(inisrf_0(1,1,1),iy*(nnsg*3+8)*jxp,mpi_real8,   &
                       & inisrf0(1,1,1), iy*(nnsg*3+8)*jxp,mpi_real8,   &
                       & 0,mpi_comm_world,ierr)

        do j = 1 , jxp
          do i = 1 , iy
            ht(i,j) = inisrf0(i,1,j)
            satbrt(i,j) = inisrf0(i,2,j)
            xlat(i,j) = inisrf0(i,3,j)
            xlong(i,j) = inisrf0(i,4,j)
            msfx(i,j) = inisrf0(i,5,j)
            msfd(i,j) = inisrf0(i,6,j)
            f(i,j) = inisrf0(i,7,j)
          end do
          do n = 1 , nnsg
            do i = 1 , iy
              ht1(n,i,j) = inisrf0(i,7+n,j)
              satbrt1(n,i,j) = inisrf0(i,7+nnsg+n,j)
            end do
          end do
        end do

        if ( myid.eq.0 ) then
          do j = 1 , jx
            do k = 1 , kz
              do i = 1 , iy
                sav_0(i,k,j) = ub0_io(i,k,j)
                sav_0(i,kz+k,j) = vb0_io(i,k,j)
                sav_0(i,kz*2+k,j) = qb0_io(i,k,j)
                sav_0(i,kz*3+k,j) = tb0_io(i,k,j)
              end do
            end do
            do i = 1 , iy
              sav_0(i,kz*4+1,j) = ps0_io(i,j)
              sav_0(i,kz*4+2,j) = ts0_io(i,j)
            end do
            if ( ehso4 ) then
              do k = 1 , kz
                do i = 1 , iy
                  sav_0s(i,k,j) = so0_io(i,k,j)
                end do
              end do
            end if
          end do
        end if
        call mpi_scatter(sav_0(1,1,1),iy*(kz*4+2)*jxp,mpi_real8,        &
                       & sav0(1,1,1), iy*(kz*4+2)*jxp,mpi_real8,        &
                       & 0,mpi_comm_world,ierr)
        if ( ehso4 )                                                    &
          &  call mpi_scatter(sav_0s(1,1,1),iy*kz*jxp,mpi_real8,        &
          &                   sav0s(1,1,1), iy*kz*jxp,mpi_real8,        &
          &                   0,mpi_comm_world,ierr)
        do j = 1 , jendl
          do k = 1 , kz
            do i = 1 , iy
              ub0(i,k,j) = sav0(i,k,j)
              vb0(i,k,j) = sav0(i,kz+k,j)
              qb0(i,k,j) = sav0(i,kz*2+k,j)
              tb0(i,k,j) = sav0(i,kz*3+k,j)
            end do
          end do
          do i = 1 , iy
            ps0(i,j) = sav0(i,kz*4+1,j)
            ts0(i,j) = sav0(i,kz*4+2,j)
          end do
          if ( ehso4 ) then
            do k = 1 , kz
              do i = 1 , iy
                so0(i,k,j) = sav0s(i,k,j)
              end do
            end do
          end if
        end do

        if ( myid.eq.0 ) then
          do j = 1 , jx
            do k = 1 , kz
              do i = 1 , iy
                sav_0(i,k,j) = ua_io(i,k,j)
                sav_0(i,kz+k,j) = ub_io(i,k,j)
                sav_0(i,kz*2+k,j) = va_io(i,k,j)
                sav_0(i,kz*3+k,j) = vb_io(i,k,j)
              end do
            end do
            do i = 1 , iy
              sav_0(i,kz*4+1,j) = psa_io(i,j)
              sav_0(i,kz*4+2,j) = psb_io(i,j)
            end do
          end do
        end if
        call mpi_scatter(sav_0(1,1,1),iy*(kz*4+2)*jxp,mpi_real8,        &
                       & sav0(1,1,1), iy*(kz*4+2)*jxp,mpi_real8,        &
                       & 0,mpi_comm_world,ierr)
        do j = 1 , jendl
          do k = 1 , kz
            do i = 1 , iy
              ua(i,k,j) = sav0(i,k,j)
              ub(i,k,j) = sav0(i,kz+k,j)
              va(i,k,j) = sav0(i,kz*2+k,j)
              vb(i,k,j) = sav0(i,kz*3+k,j)
            end do
          end do
          do i = 1 , iy
            psa(i,j) = sav0(i,kz*4+1,j)
            psb(i,j) = sav0(i,kz*4+2,j)
          end do
        end do
        if ( myid.eq.0 ) then
          do j = 1 , jx
            do k = 1 , kz
              do i = 1 , iy
                sav_0(i,k,j) = ta_io(i,k,j)
                sav_0(i,kz+k,j) = tb_io(i,k,j)
                sav_0(i,kz*2+k,j) = qva_io(i,k,j)
                sav_0(i,kz*3+k,j) = qvb_io(i,k,j)
              end do
            end do
            do i = 1 , iy
              sav_0(i,kz*4+1,j) = tga_io(i,j)
              sav_0(i,kz*4+2,j) = tgb_io(i,j)
            end do
          end do
        end if
        call mpi_scatter(sav_0(1,1,1),iy*(kz*4+2)*jxp,mpi_real8,        &
                       & sav0(1,1,1), iy*(kz*4+2)*jxp,mpi_real8,        &
                       & 0,mpi_comm_world,ierr)
        do j = 1 , jendl
          do k = 1 , kz
            do i = 1 , iy
              ta(i,k,j) = sav0(i,k,j)
              tb(i,k,j) = sav0(i,kz+k,j)
              qva(i,k,j) = sav0(i,kz*2+k,j)
              qvb(i,k,j) = sav0(i,kz*3+k,j)
            end do
          end do
          do i = 1 , iy
            tga(i,j) = sav0(i,kz*4+1,j)
            tgb(i,j) = sav0(i,kz*4+2,j)
          end do
        end do
        if ( myid.eq.0 ) then
          do j = 1 , jx
            do k = 1 , kz
              do i = 1 , iy
                sav_0(i,k,j) = qca_io(i,k,j)
                sav_0(i,kz+k,j) = qcb_io(i,k,j)
                sav_0(i,kz*2+k,j) = fcc_io(i,k,j)
              end do
            end do
            do i = 1 , iy
              sav_0(i,kz*4+1,j) = rainc_io(i,j)
              sav_0(i,kz*4+2,j) = rainnc_io(i,j)
            end do
          end do
          do j = 1 , jxm1
            do k = 1 , kz
              do i = 1 , iym1
                sav_0(i,kz*3+k,j) = heatrt_io(i,k,j)
              end do
            end do
          end do
        end if
        call mpi_scatter(sav_0(1,1,1),iy*(kz*4+2)*jxp,mpi_real8,        &
                       & sav0(1,1,1), iy*(kz*4+2)*jxp,mpi_real8,        &
                       & 0,mpi_comm_world,ierr)
        do j = 1 , jendl
          do k = 1 , kz
            do i = 1 , iy
              qca(i,k,j) = sav0(i,k,j)
              qcb(i,k,j) = sav0(i,kz+k,j)
              fcc(i,k,j) = sav0(i,kz*2+k,j)
            end do
          end do
          do i = 1 , iy
            rainc(i,j) = sav0(i,kz*4+1,j)
            rainnc(i,j) = sav0(i,kz*4+2,j)
          end do
        end do
        do j = 1 , jendx
          do k = 1 , kz
            do i = 1 , iym1
              heatrt(i,k,j) = sav0(i,kz*3+k,j)
            end do
          end do
        end do
        if ( myid.eq.0 ) then
          do j = 1 , jx
            do i = 1 , iy
              sav_0a(i,1,j) = hfx_io(i,j)
              sav_0a(i,2,j) = qfx_io(i,j)
              sav_0a(i,3,j) = uvdrag_io(i,j)
              sav_0a(i,4,j) = tgbb_io(i,j)
            end do
            do n = 1 , nnsg
              do i = 1 , iy
                sav_0a(i,4+n,j) = snowc_io(n,i,j)
              end do
            end do
          end do
          do j = 1 , jxm1
            do k = 1 , kzp1
              do i = 1 , iym1
                sav_0a(i,nnsg+4+k,j) = o3prof_io(i,k,j)
              end do
            end do
          end do
        end if
        allrec = kz + 5 + nnsg
        call mpi_scatter(sav_0a(1,1,1),iy*allrec*jxp,mpi_real8,         &
                       & sav0a(1,1,1), iy*allrec*jxp,mpi_real8,         &
                       & 0,mpi_comm_world,ierr)
        do j = 1 , jendl
          do i = 1 , iy
            hfx(i,j) = sav0a(i,1,j)
            qfx(i,j) = sav0a(i,2,j)
            uvdrag(i,j) = sav0a(i,3,j)
            tgbb(i,j) = sav0a(i,4,j)
          end do
          do n = 1 , nnsg
            do i = 1 , iy
              snowc(n,i,j) = sav0a(i,4+n,j)
            end do
          end do
        end do
        do j = 1 , jendx
          do k = 1 , kzp1
            do i = 1 , iym1
              o3prof(i,k,j) = sav0a(i,nnsg+4+k,j)
            end do
          end do
        end do
        if ( iocnflx.eq.2 )                                             &
          & call mpi_scatter(zpbl_io(1,1),iy*jxp,mpi_real8,             &
          &                  zpbl(1,1),   iy*jxp,mpi_real8,             &
          &                  0,mpi_comm_world,ierr)
        if ( icup.eq.1 ) then
          if ( myid.eq.0 ) then
            do j = 1 , jx
              do k = 1 , kz
                do i = 1 , iy
                  sav_0c(i,k,j) = rsheat_io(i,k,j)
                  sav_0c(i,kz+k,j) = rswat_io(i,k,j)
                end do
              end do
            end do
          end if
          call mpi_scatter(sav_0c(1,1,1),iy*kz*2*jxp,mpi_real8,         &
                         & sav0c(1,1,1), iy*kz*2*jxp,mpi_real8,         &
                         & 0,mpi_comm_world,ierr)
          do j = 1 , jendl
            do k = 1 , kz
              do i = 1 , iy
                rsheat(i,k,j) = sav0c(i,k,j)
                rswat(i,k,j) = sav0c(i,kz+k,j)
              end do
            end do
          end do
        else if ( icup.eq.3 ) then
          if ( myid.eq.0 ) then
            do j = 1 , jx
              do k = 1 , kz
                do i = 1 , iy
                  sav_0b(i,k,j) = tbase_io(i,k,j)
                end do
              end do
              do i = 1 , iy
                sav_0b(i,kzp1,j) = cldefi_io(i,j)
              end do
            end do
          end if
          call mpi_scatter(sav_0b(1,1,1),iy*(kzp1)*jxp,mpi_real8,       &
                         & sav0b(1,1,1), iy*(kzp1)*jxp,mpi_real8,       &
                         & 0,mpi_comm_world,ierr)
          do j = 1 , jendl
            do k = 1 , kz
              do i = 1 , iy
                tbase(i,k,j) = sav0b(i,k,j)
              end do
            end do
            do i = 1 , iy
              cldefi(i,j) = sav0b(i,kzp1,j)
            end do
          end do
        else if ( icup.eq.4 ) then
          call mpi_scatter(cbmf2d_io(1,1),iy*jxp,mpi_real8,             &
                         & cbmf2d(1,1),   iy*jxp,mpi_real8,             &
                         & 0,mpi_comm_world,ierr)
        else
        end if
        if ( myid.eq.0 ) then
          do j = 1 , jxm1
            do l = 1 , 4
              do k = 1 , kz
                do i = 1 , iym1
                  sav_1(i,(l-1)*kz+k,j) = absnxt_io(i,k,l,j)
                end do
              end do
            end do
          end do
          allrec = kz*4
          do j = 1 , jxm1
            do l = 1 , kzp1
              do k = 1 , kzp1
                do i = 1 , iym1
                  sav_1(i,allrec+(l-1)*(kzp1)+k,j) = abstot_io(i,k,l,j)
                end do
              end do
            end do
          end do
          allrec = allrec + (kzp1)*(kz+1)
          do j = 1 , jxm1
            do k = 1 , kzp1
              do i = 1 , iym1
                sav_1(i,allrec+k,j) = emstot_io(i,k,j)
              end do
            end do
          end do
          allrec = allrec + kzp1
        end if
        allrec = kz*4 + (kzp1*kzp2)
        call mpi_scatter(sav_1(1,1,1),iym1*allrec*jxp,mpi_real8,        &
                       & sav1(1,1,1), iym1*allrec*jxp,mpi_real8,        &
                       & 0,mpi_comm_world,ierr)
        do j = 1 , jendx
          do l = 1 , 4
            do k = 1 , kz
              do i = 1 , iym1
                absnxt(i,k,l,j) = sav1(i,(l-1)*kz+k,j)
              end do
            end do
          end do
        end do
        allrec = kz*4
        do j = 1 , jendx
          do l = 1 , kzp1
            do k = 1 , kzp1
              do i = 1 , iym1
                abstot(i,k,l,j) = sav1(i,allrec+(l-1)*(kzp1)+k,j)
              end do
            end do
          end do
        end do
        allrec = allrec + (kzp1)*(kz+1)
        do j = 1 , jendx
          do k = 1 , kzp1
            do i = 1 , iym1
              emstot(i,k,j) = sav1(i,allrec+k,j)
            end do
          end do
        end do
        if ( myid.eq.0 ) then
          do j = 1 , jxm1
            do n = 1 , nnsg
              do i = 1 , iym1
                sav_2(i,n,j) = taf2d_io(n,i,j)
                sav_2(i,nnsg+n,j) = tlef2d_io(n,i,j)
                sav_2(i,nnsg*2+n,j) = ssw2d_io(n,i,j)
                sav_2(i,nnsg*3+n,j) = srw2d_io(n,i,j)
              end do
            end do
            do i = 1 , iym1
              sav_2(i,nnsg*4+1,j) = sol2d_io(i,j)
              sav_2(i,nnsg*4+2,j) = solvd2d_io(i,j)
              sav_2(i,nnsg*4+3,j) = solvs2d_io(i,j)
              sav_2(i,nnsg*4+4,j) = flw2d_io(i,j)
            end do
          end do
        end if
        allrec = nnsg*4 + 4
        call mpi_scatter(sav_2(1,1,1),iym1*allrec*jxp,mpi_real8,        &
                       & sav2(1,1,1), iym1*allrec*jxp,mpi_real8,        &
                       & 0,mpi_comm_world,ierr)
        do j = 1 , jendx
          do n = 1 , nnsg
            do i = 1 , iym1
              taf2d(n,i,j) = sav2(i,n,j)
              tlef2d(n,i,j) = sav2(i,nnsg+n,j)
              ssw2d(n,i,j) = sav2(i,nnsg*2+n,j)
              srw2d(n,i,j) = sav2(i,nnsg*3+n,j)
            end do
          end do
          do i = 1 , iym1
            sol2d(i,j) = sav2(i,nnsg*4+1,j)
            solvd2d(i,j) = sav2(i,nnsg*4+2,j)
            solvs2d(i,j) = sav2(i,nnsg*4+3,j)
            flw2d(i,j) = sav2(i,nnsg*4+4,j)
          end do
        end do
        if ( myid.eq.0 ) then
          do j = 1 , jxm1
            do n = 1 , nnsg
              do i = 1 , iym1
                sav_2(i,n,j) = tgb2d_io(n,i,j)
                sav_2(i,nnsg+n,j) = swt2d_io(n,i,j)
                sav_2(i,nnsg*2+n,j) = scv2d_io(n,i,j)
                sav_2(i,nnsg*3+n,j) = gwet2d_io(n,i,j)
              end do
            end do
            do i = 1 , iym1
              sav_2(i,nnsg*4+1,j) = flwd2d_io(i,j)
              sav_2(i,nnsg*4+2,j) = fsw2d_io(i,j)
              sav_2(i,nnsg*4+3,j) = sabv2d_io(i,j)
              sav_2(i,nnsg*4+4,j) = sinc2d_io(i,j)
            end do
          end do
        end if
        allrec = nnsg*4 + 4
        call mpi_scatter(sav_2(1,1,1),iym1*allrec*jxp,mpi_real8,        &
                       & sav2(1,1,1), iym1*allrec*jxp,mpi_real8,        &
                       & 0,mpi_comm_world,ierr)
        do j = 1 , jendx
          do n = 1 , nnsg
            do i = 1 , iym1
              tgb2d(n,i,j) = sav2(i,n,j)
              swt2d(n,i,j) = sav2(i,nnsg+n,j)
              scv2d(n,i,j) = sav2(i,nnsg*2+n,j)
              gwet2d(n,i,j) = sav2(i,nnsg*3+n,j)
            end do
          end do
          do i = 1 , iym1
            flwd2d(i,j) = sav2(i,nnsg*4+1,j)
            fsw2d(i,j) = sav2(i,nnsg*4+2,j)
            sabv2d(i,j) = sav2(i,nnsg*4+3,j)
            sinc2d(i,j) = sav2(i,nnsg*4+4,j)
          end do
        end do
        if ( myid.eq.0 ) then
          do j = 1 , jxm1
            do n = 1 , nnsg
              do i = 1 , iym1
                sav_2(i,n,j) = veg2d1_io(n,i,j)
                sav_2(i,nnsg+n,j) = sag2d_io(n,i,j)
                sav_2(i,nnsg*2+n,j) = sice2d_io(n,i,j)
                sav_2(i,nnsg*3+n,j) = dew2d_io(n,i,j)
              end do
            end do
            do i = 1 , iym1
              sav_2(i,nnsg*4+1,j) = pptnc_io(i,j)
              sav_2(i,nnsg*4+2,j) = pptc_io(i,j)
              sav_2(i,nnsg*4+3,j) = prca2d_io(i,j)
              sav_2(i,nnsg*4+4,j) = prnca2d_io(i,j)
            end do
          end do
        end if
        allrec = nnsg*4 + 4
        call mpi_scatter(sav_2(1,1,1),iym1*allrec*jxp,mpi_real8,        &
                       & sav2(1,1,1), iym1*allrec*jxp,mpi_real8,        &
                       & 0,mpi_comm_world,ierr)
        do j = 1 , jendx
          do n = 1 , nnsg
            do i = 1 , iym1
              veg2d1(n,i,j) = sav2(i,n,j)
              sag2d(n,i,j) = sav2(i,nnsg+n,j)
              sice2d(n,i,j) = sav2(i,nnsg*2+n,j)
              dew2d(n,i,j) = sav2(i,nnsg*3+n,j)
            end do
          end do
          do i = 1 , iym1
            pptnc(i,j) = sav2(i,nnsg*4+1,j)
            pptc(i,j) = sav2(i,nnsg*4+2,j)
            prca2d(i,j) = sav2(i,nnsg*4+3,j)
            prnca2d(i,j) = sav2(i,nnsg*4+4,j)
          end do
        end do
        if ( myid.eq.0 ) then
          do j = 1 , jxm1
            do n = 1 , nnsg
              do i = 1 , iym1
                sav_2a(i,n,j) = ircp2d_io(n,i,j)
                sav_2a(i,nnsg+n,j) = text2d_io(n,i,j)
                sav_2a(i,nnsg*2+n,j) = col2d_io(n,i,j)
                sav_2a(i,nnsg*3+n,j) = ocld2d_io(n,i,j)
                sav_2a(i,nnsg*4+n,j) = tg2d_io(n,i,j)
              end do
            end do
            do i = 1 , iym1
              sav_2a(i,nnsg*5+1,j) = veg2d_io(i,j)
            end do
          end do
        end if
        allrec = nnsg*5 + 1
        call mpi_scatter(sav_2a(1,1,1),iym1*allrec*jxp,mpi_real8,       &
                       & sav2a(1,1,1), iym1*allrec*jxp,mpi_real8,       &
                       & 0,mpi_comm_world,ierr)
        do j = 1 , jendx
          do n = 1 , nnsg
            do i = 1 , iym1
              ircp2d(n,i,j) = sav2a(i,n,j)
              text2d(n,i,j) = sav2a(i,nnsg+n,j)
              col2d(n,i,j) = sav2a(i,nnsg*2+n,j)
              ocld2d(n,i,j) = sav2a(i,nnsg*3+n,j)
              tg2d(n,i,j) = sav2a(i,nnsg*4+n,j)
            end do
          end do
          do i = 1 , iym1
            veg2d(i,j) = sav2a(i,nnsg*5+1,j)
          end do
        end do
        if ( ichem.eq.1 ) then
          if ( myid.eq.0 ) then
            do j = 1 , jx
              do n = 1 , ntr
                do k = 1 , kz
                  do i = 1 , iy
                    sav_4(i,(n-1)*kz+k,j) = chia_io(i,k,j,n)
                    sav_4(i,ntr*kz+(n-1)*kz+k,j) = chib_io(i,k,j,n)
                    sav_4(i,ntr*kz*2+(n-1)*kz+k,j) = remlsc_io(i,k,j,n)
                    sav_4(i,ntr*kz*3+(n-1)*kz+k,j) = remcvc_io(i,k,j,n)
                  end do
                end do
              end do
            end do
            allrec = 4*ntr*kz
            do j = 1 , jx
              do n = 1 , ntr
                do i = 1 , iy
                  sav_4(i,allrec+n,j) = remdrd_io(i,j,n)
                end do
              end do
            end do
            allrec = allrec + ntr
          end if
          allrec = ntr*(kz*4+1)
          call mpi_scatter(sav_4(1,1,1),iy*allrec*jxp,mpi_real8,        &
                         & sav4(1,1,1), iy*allrec*jxp,mpi_real8,        &
                         & 0,mpi_comm_world,ierr)
          do j = 1 , jendl
            do n = 1 , ntr
              do k = 1 , kz
                do i = 1 , iy
                  chia(i,k,j,n) = sav4(i,(n-1)*kz+k,j)
                  chib(i,k,j,n) = sav4(i,ntr*kz+(n-1)*kz+k,j)
                  remlsc(i,k,j,n) = sav4(i,ntr*kz*2+(n-1)*kz+k,j)
                  remcvc(i,k,j,n) = sav4(i,ntr*kz*3+(n-1)*kz+k,j)
                end do
              end do
            end do
          end do
          allrec = 4*ntr*kz
          do j = 1 , jendl
            do n = 1 , ntr
              do i = 1 , iy
                remdrd(i,j,n) = sav4(i,allrec+n,j)
              end do
            end do
          end do
          if ( myid.eq.0 ) then
            do j = 1 , jxm1
              do i = 1 , iym1
                sav_4a(i,1,j) = ssw2da_io(i,j)
                sav_4a(i,2,j) = sdeltk2d_io(i,j)
                sav_4a(i,3,j) = sdelqk2d_io(i,j)
                sav_4a(i,4,j) = sfracv2d_io(i,j)
                sav_4a(i,5,j) = sfracb2d_io(i,j)
                sav_4a(i,6,j) = sfracs2d_io(i,j)
                sav_4a(i,7,j) = svegfrac2d_io(i,j)
              end do
            end do
          end if
          call mpi_scatter(sav_4a,iym1*7*jxp,mpi_real8,                 &
                         & sav4a, iym1*7*jxp,mpi_real8,                 &
                         & 0,mpi_comm_world,ierr)
          do j = 1 , jendx
            do i = 1 , iym1
              ssw2da(i,j) = sav4a(i,1,j)
              sdeltk2d(i,j) = sav4a(i,2,j)
              sdelqk2d(i,j) = sav4a(i,3,j)
              sfracv2d(i,j) = sav4a(i,4,j)
              sfracb2d(i,j) = sav4a(i,5,j)
              sfracs2d(i,j) = sav4a(i,6,j)
              svegfrac2d(i,j) = sav4a(i,7,j)
            end do
          end do
#ifdef CLM
          if ( myid.eq.0 ) then
            do j = 1 , jxm1
              do i = 1 , iym1
                sav_clmout(i,1,j) = sols2d_io(i,j)
                sav_clmout(i,2,j) = soll2d_io(i,j)
                sav_clmout(i,3,j) = solsd2d_io(i,j)
                sav_clmout(i,4,j) = solld2d_io(i,j)
                sav_clmout(i,5,j) = aldirs2d_io(i,j)
                sav_clmout(i,6,j) = aldirl2d_io(i,j)
                sav_clmout(i,7,j) = aldifs2d_io(i,j)
                sav_clmout(i,8,j) = aldifl2d_io(i,j)
                sav_clmout(i,9,j) = coszrs2d_io(i,j)
              end do
            end do
          end if
          call mpi_scatter(sav_clmout,iym1*9*jxp,mpi_real8,             &
                         & sav_clmin, iym1*9*jxp,mpi_real8,             &
                         & 0,mpi_comm_world,ierr)
          do j = 1 , jendx
            do i = 1 , iym1
              sols2d(i,j) = sav_clmin(i,1,j)
              soll2d(i,j) = sav_clmin(i,2,j)
              solsd2d(i,j) = sav_clmin(i,3,j)
              solld2d(i,j) = sav_clmin(i,4,j)
              aldirs2d(i,j) = sav_clmin(i,5,j)
              aldirl2d(i,j) = sav_clmin(i,6,j)
              aldifs2d(i,j) = sav_clmin(i,7,j)
              aldifl2d(i,j) = sav_clmin(i,8,j)
              coszrs2d(i,j) = sav_clmin(i,9,j)
            end do
          end do
#endif
        end if
        call mpi_bcast(mdate0,1,mpi_integer,0,mpi_comm_world,ierr)
        call mpi_bcast(jyear0,1,mpi_integer,0,mpi_comm_world,ierr)
        call mpi_bcast(ktau,1,mpi_integer,0,mpi_comm_world,ierr)
        call mpi_bcast(jyear,1,mpi_integer,0,mpi_comm_world,ierr)
        call mpi_bcast(xtime,1,mpi_real8,0,mpi_comm_world,ierr)
        call mpi_bcast(ldatez,1,mpi_integer,0,mpi_comm_world,ierr)
        call mpi_bcast(lyear,1,mpi_integer,0,mpi_comm_world,ierr)
        call mpi_bcast(lmonth,1,mpi_integer,0,mpi_comm_world,ierr)
        call mpi_bcast(lday,1,mpi_integer,0,mpi_comm_world,ierr)
        call mpi_bcast(lhour,1,mpi_integer,0,mpi_comm_world,ierr)
        call mpi_bcast(ntime,1,mpi_integer,0,mpi_comm_world,ierr)
        call mpi_bcast(jyearr,1,mpi_integer,0,mpi_comm_world,ierr)
        call mpi_bcast(ktaur,1,mpi_integer,0,mpi_comm_world,ierr)
#ifdef DIAG
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
#endif
!------lake model
        if ( lakemod.eq.1 ) then
          call mpi_bcast(ilake,1,mpi_integer,0,mpi_comm_world,ierr)
          call mpi_bcast(jlake,1,mpi_integer,0,mpi_comm_world,ierr)
          call mpi_bcast(depth,1,mpi_integer,0,mpi_comm_world,ierr)
          call mpi_bcast(hs,1,mpi_real8,0,mpi_comm_world,ierr)
          call mpi_bcast(eta,1,mpi_real8,0,mpi_comm_world,ierr)
          call mpi_bcast(tlake,depth,mpi_real8,0,mpi_comm_world,ierr)
        end if
        call mpi_sendrecv(ht(1,jxp),iy,mpi_real8,ieast,1,               &
                        & ht(1,0),iy,mpi_real8,iwest,1,                 &
                        & mpi_comm_world,mpi_status_ignore,ierr)
        call mpi_sendrecv(ht(1,1),iy,mpi_real8,iwest,2,                 &
                        & ht(1,jxp+1),iy,mpi_real8,ieast,2,             &
                        & mpi_comm_world,mpi_status_ignore,ierr)
        call mpi_sendrecv(msfx(1,jxp-1),iy*2,mpi_real8,ieast,           &
                        & 1,msfx(1,-1),iy*2,mpi_real8,iwest,            &
                        & 1,mpi_comm_world,mpi_status_ignore,ierr)
        call mpi_sendrecv(msfx(1,1),iy*2,mpi_real8,iwest,2,             &
                        & msfx(1,jxp+1),iy*2,mpi_real8,ieast,           &
                        & 2,mpi_comm_world,mpi_status_ignore,ierr)
        call mpi_sendrecv(msfd(1,jxp-1),iy*2,mpi_real8,ieast,           &
                        & 1,msfd(1,-1),iy*2,mpi_real8,iwest,            &
                        & 1,mpi_comm_world,mpi_status_ignore,ierr)
        call mpi_sendrecv(msfd(1,1),iy*2,mpi_real8,iwest,2,             &
                        & msfd(1,jxp+1),iy*2,mpi_real8,ieast,           &
                        & 2,mpi_comm_world,mpi_status_ignore,ierr)

        dt = dt2 ! First timestep successfully read in

#else
        write (finm,99002) trim(dirout),pthsep,'SAV.',                  &
          &    ((ndate0/10000)*100+1)*100
        inquire (file=finm,exist=existing)
        if ( .not.existing ) then
          write (aline,*) 'The following SAV File does not exist: ' ,   &
              &            trim(finm), 'please check location'
          call say
          call fatal(__FILE__,__LINE__, 'SAV FILE NOT FOUND')
        else
          open (iutrs,file=finm,form='unformatted',status='old')
        end if
        do ! Loop while ldatez.ne.idate1
!
!-----when ifrest=.true., read in the data saved from previous run
!         for large domain from unit 14.
!
          read (iutrs) mdate0
          jyear0 = mdate0/1000000
          read (iutrs) ktau , xtime , ldatez , lyear , lmonth , lday ,  &
                     & lhour , ntime
          jyear = lyear
          jyearr = jyear
          ktaur = ktau
          if ( ehso4 ) then
            read (iutrs) ub0 , vb0 , qb0 , tb0 , ps0 , ts0 , so0
          else
            read (iutrs) ub0 , vb0 , qb0 , tb0 , ps0 , ts0
          end if
          read (iutrs) ua
          read (iutrs) ub
          read (iutrs) va
          read (iutrs) vb
          read (iutrs) ta
          read (iutrs) tb
          read (iutrs) qva
          read (iutrs) qvb
          read (iutrs) qca
          read (iutrs) qcb
          read (iutrs) psa , psb , satbrt , satbrt1 , f
          read (iutrs) ht , ht1 , msfx , msfd , xlat , xlong
          read (iutrs) tga , tgb , rainc , rainnc
          if ( icup.eq.1 ) then
            read (iutrs) rsheat , rswat
          else if ( icup.eq.3 ) then
            read (iutrs) tbase , cldefi
          else if ( icup.eq.4 ) then
            read (iutrs) cbmf2d
          else
          end if
          read (iutrs) hfx , qfx , snowc , uvdrag
#ifdef    DIAG
          read (iutrs) tdini , tdadv , tqini , tqadv , tqeva , tqrai
#endif
          read (iutrs) absnxt , abstot , emstot
          if ( ipptls.eq.1 ) read (iutrs) fcc
          read (iutrs) sol2d , solvd2d , solvs2d , flw2d , flwd2d ,     &
                     & fsw2d , sabv2d , sinc2d
          read (iutrs) taf2d , tlef2d , tgbb , ssw2d , srw2d , tg2d ,   &
                     & tgb2d , swt2d , scv2d , gwet2d , veg2d , veg2d1 ,&
                     & sag2d , sice2d , dew2d , ircp2d , text2d ,       &
                     & col2d , ocld2d , heatrt , o3prof
          read (iutrs) pptnc , pptc , prca2d , prnca2d
          if ( iocnflx.eq.2 ) read (iutrs) zpbl

!chem2---
          if ( ichem.eq.1 ) then
            read (iutrs) chia
            read (iutrs) chib
!           cumul removal terms (3d, 2d)
            read (iutrs) remlsc
            read (iutrs) remcvc
            read (iutrs) remdrd
            read (iutrs) ssw2da
            read (iutrs) sdeltk2d
            read (iutrs) sdelqk2d
            read (iutrs) sfracv2d
            read (iutrs) sfracb2d
            read (iutrs) sfracs2d
            read (iutrs) svegfrac2d
!           cumul ad, dif, emis terms ( scalar)
#ifdef DIAG
            read (iutrs) tchiad
            read (iutrs) tchitb
            read (iutrs) tchie
#endif
          end if
!chem2_
!
!------lake model
          if ( lakemod.eq.1 ) then
            lcount = 0
            iin = 41
            iout = 42
            rewind (iin)
            read (iutrs) numpts
            print * , 'reading lake model restart file. numpts = ' ,    &
                & numpts
            print * , 'jyear, ktau, xtime = ' , jyear , ktau , xtime
            do n = 1 , numpts
              read (iutrs) ilake , jlake , depth , freeze , hi , hii ,  &
                         & hs , eta , (tlake(j),j=1,depth)
              print * , 'reading restart file at i, j = ' , ilake ,     &
                  & jlake
              write (iin) ilake , jlake , depth , freeze , hi , hii ,   &
                        & hs , eta , (tlake(j),j=1,depth)
            end do
            rewind (iin)
          end if
!
          print * , 'ozone profiles restart'
          do k = 1 , kzp1
            write (6,99004) o3prof(3,3,k)
          end do
          print 99005 , xtime , ktau , jyear , iutrs
!
          if ( ldatez.ne.idate1 ) then
            write (*,*) 'INIT: ldatez, idate1=' , ldatez , idate1
            cycle
          end if

          exit ! We have now ldatez.eq.idate1

        end do

        dt = dt2 ! First timestep successfully read in

#endif
!
!-----end of initial/restart if test
!
      end if
!
!     Move from param.F to fix the reatart problem found by
!     Zhang DongFeng
!
      if ( ipptls.eq.1 ) then
#ifdef MPP1
        do j = 1 , jendx
          do i = 1 , iym1
            if ( satbrt(i,j).gt.13.9 .and. satbrt(i,j).lt.15.1 ) then
              qck1(i,j) = qck1oce  ! OCEAN
              cgul(i,j) = guloce   ! OCEAN
              rh0(i,j) = rh0oce    ! OCEAN
            else
              qck1(i,j) = qck1land ! LAND
              cgul(i,j) = gulland  ! LAND
              rh0(i,j) = rh0land   ! LAND
            end if
          end do
        end do
#else
        do j = 1 , jxm1
          do i = 1 , iym1
            if ( satbrt(i,j).gt.13.9 .and. satbrt(i,j).lt.15.1 ) then
              qck1(i,j) = qck1oce  ! OCEAN
              cgul(i,j) = guloce   ! OCEAN
              rh0(i,j) = rh0oce    ! OCEAN
            else
              qck1(i,j) = qck1land ! LAND
              cgul(i,j) = gulland  ! LAND
              rh0(i,j) = rh0land   ! LAND
            end if
          end do
        end do
#endif
      end if
!chem2
      if ( ichem.eq.1 ) then
        iso2 = 0
        iso4 = 0
        ibchl = 0
        ibchb = 0
        iochl = 0
        iochb = 0
        ibin = 0
        do itr = 1 , ntr
          if ( chtrname(itr).eq.'SO2' ) iso2 = itr
          if ( chtrname(itr).eq.'SO4' ) iso4 = itr
          if ( chtrname(itr).eq.'BC_HL' ) ibchl = itr
          if ( chtrname(itr).eq.'BC_HB' ) ibchb = itr
          if ( chtrname(itr).eq.'OC_HL' ) iochl = itr
          if ( chtrname(itr).eq.'OC_HB' ) iochb = itr
          if ( chtrname(itr).eq.'DUST' ) then
            ibin = ibin + 1
            idust(ibin) = itr
          end if
        end do
      end if
!chem2_
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!     ****** initialize and define constants for vector bats
 
      if ( jyear.eq.jyear0 .and. ktau.eq.0 ) call initb
      if ( iemiss.eq.1 ) then
#ifdef MPP1
        do j = 1 , jendx
          do i = 1 , iym1
            do n = 1 , nnsg
              ist = nint(veg2d1(n,i,j))
              if ( ist.eq.0 ) then
                emiss2d(n,i,j) = 0.955D0
              else if ( ist.eq.8 ) then
                emiss2d(n,i,j) = 0.76D0
              else if ( ist.eq.11 ) then
                emiss2d(n,i,j) = 0.85D0
              else if ( ist.eq.12 ) then
                emiss2d(n,i,j) = 0.97D0
              else
                emiss2d(n,i,j) = 0.99D0 - (albvgs(ist)+albvgl(ist))     &
                               & *0.1D0
              end if
!             emiss2d(n,i,j) = 1.0d0
            end do
          end do
        end do
#else
        do j = 1 , jxm1
          do i = 1 , iym1
            do n = 1 , nnsg
              ist = nint(veg2d1(n,i,j))
              if ( ist.eq.0 ) then
                emiss2d(n,i,j) = 0.955D0
              else if ( ist.eq.8 ) then
                emiss2d(n,i,j) = 0.76D0
              else if ( ist.eq.11 ) then
                emiss2d(n,i,j) = 0.85D0
              else if ( ist.eq.12 ) then
                emiss2d(n,i,j) = 0.97D0
              else
                emiss2d(n,i,j) = 0.99D0 - (albvgs(ist)+albvgl(ist))     &
                               & *0.1D0
              end if
!             emiss2d(n,i,j) = 1.0d0
            end do
          end do
        end do
#endif
      end if
! ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!-----read in the boundary conditions for large domain:
!
!-----compute the solar declination angle:
!
#ifdef CLM
      call solar1clm(xtime)
      if ( ( jyear.eq.jyear0 .and. ktau.eq.0 ) .or.                     &
         & ( ktau.eq.ktaur ) ) then
        init_grid = .true.
      else
        init_grid = .false.
      end if
#else
      call solar1(xtime)
#endif
      call inirad
!
!-----calculating topographical correction to diffusion coefficient
#ifdef MPP1
      do j = 1 , jendl
        do i = 1 , iy
          hgfact(i,j) = 1.
        end do
      end do
      do j = jbegin , jendm
        if ( myid.eq.0 ) then
          jm1h = max0(j-1,2)
        else
          jm1h = j - 1
        end if
        if ( myid.eq.nproc-1 ) then
          jp1h = min0(j+1,jxp-2)
        else
          jp1h = j + 1
        end if
        do i = 2 , iym2
          im1h = max0(i-1,2)
          ip1h = min0(i+1,iym2)
          hg1 = dabs((ht(i,j)-ht(im1h,j))/dx)
          hg2 = dabs((ht(i,j)-ht(ip1h,j))/dx)
          hg3 = dabs((ht(i,j)-ht(i,jm1h))/dx)
          hg4 = dabs((ht(i,j)-ht(i,jp1h))/dx)
          hgmax = dmax1(hg1,hg2,hg3,hg4)*rgti
          hgfact(i,j) = 1./(1.+(hgmax/0.001)**2.)
        end do
      end do
#else
      do j = 1 , jx
        do i = 1 , iy
          hgfact(i,j) = 1.
        end do
      end do
      do j = 2 , jxm2
        jm1h = max0(j-1,2)
        jp1h = min0(j+1,jxm2)
        do i = 2 , iym2
          im1h = max0(i-1,2)
          ip1h = min0(i+1,iym2)
          hg1 = dabs((ht(i,j)-ht(im1h,j))/dx)
          hg2 = dabs((ht(i,j)-ht(ip1h,j))/dx)
          hg3 = dabs((ht(i,j)-ht(i,jm1h))/dx)
          hg4 = dabs((ht(i,j)-ht(i,jp1h))/dx)
          hgmax = dmax1(hg1,hg2,hg3,hg4)*rgti
          hgfact(i,j) = 1./(1.+(hgmax/0.001)**2.)
        end do
      end do
#endif
!
!-----set up output time:
!
      icnt = 11        ! set counter for safety-file moves on cycad
      dectim = anint(xtime+dectim)
      write (aline, *) 'dectim = ' , dectim
      call say

99001 format (a,a,a,a,i0.10)
99002 format (a,a,a,i0.10)
#ifdef DIAG
99003 format (' *** initial total air = ',e12.5,' kg, total water = ',  &
            & e12.5,' kg in large domain.')
#endif
99004 format (1x,7E12.4)
99005 format (' ***** restart file for large domain at time = ',f8.0,   &
             &' minutes, ktau = ',i7,' in year = ',i4,                  &
             &'  read in from ',a)
!
      end subroutine init
