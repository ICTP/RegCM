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
 
      subroutine bdyin

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!     this subroutine reads in the boundary conditions.               c
!                                                                     c
!        iunit : is the unit number from which the data are read in.  c
!                                                                     c
!        xtime : is the time in minutes into the forecast.            c
!                                                                     c
!        bdytim : is the time in minutes after which the boundary     c
!                 conditions are needed.                              c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      use regcm_param
      use iunits
      use param1
      use param3
      use main
      use bdycod
      use mod_bats , only : veg2d
      use date
      use tmpsav
#ifdef MPP1
      use mpi
#endif
      implicit none
!
! Local variables
!
      real(8) :: dtbdys
      character(8) :: finm
      integer :: i , ierr1 , j , k , nn , nnb
#ifdef MPP1
      integer :: ierr , ndeb , ndwb , nkk , nxeb , nxwb
      integer , dimension(mpi_status_size) :: status
      real(4) , dimension(ix,mjx) :: io2d
      real(8) , dimension(ix,jxp) :: psdot , tdum
#else
      real(4) , dimension(ix,jx) :: io2d
      real(8) , dimension(ix,jx) :: psdot , tdum
#endif
!
      if ( dabs(xtime).gt.0.0001 ) return

#ifdef MPP1
      do
!
        dtbdys = ibdyfrq*60.*60.
        if ( myid.eq.0 ) then
          do
            if ( ehso4 ) then
              do k = 1 , kx
                do j = 1 , jendl
                  do i = 1 , ix
                    so4(i,k,j) = so0(i,k,j)
                  end do
                end do
              end do
            end if
            mmrec = mmrec + 1
            read (iutbc,rec=mmrec,iostat=ierr1) ndate1
            if ( ierr1.ne.0 ) then
              close (iutbc)
              iutbc = iutbc + 1
              write (finm,99001) iutbc
!_sgi         101    format('fort.',I2)
              open (iutbc,file=finm,form='unformatted',status='old',    &
                  & access='direct',recl=ix*mjx*ibyte)
              mmrec = 0
              print * , 'CHANGING BDY UNIT NUMBER:  iutbc=' , iutbc
              if ( iutbc.gt.999 )                                       &
                  & call fatal(__FILE__,__LINE__,'BDY UNIT MAX EXCEEDED'&
                 & )
              cycle
            end if
            if ( ndate1.lt.mdatez(nnnchk+1) ) then
              if ( ndate1.lt.mdatez(nnnchk+1) ) then
                print * , 'Searching for proper date: ' , ndate1 ,      &
                    & mdatez(nnnchk+1)
                print * , 'read in datasets at :' , ndate0
                if ( ehso4 ) then
                  if ( lsmtyp.ne.'USGS' ) then
                    mmrec = mmrec + kx*5 + 2
                  else
                    mmrec = mmrec + kx*5 + 2 + 13
                  end if
                else if ( lsmtyp.ne.'USGS' ) then
                  mmrec = mmrec + kx*4 + 2
                else
                  mmrec = mmrec + kx*4 + 2 + 13
                end if
                cycle
              end if
            else if ( ndate1.gt.mdatez(nnnchk+1) ) then
              print * , 'DATE IN BC FILE EXCEEDED DATE IN RegCM'
              print * , ndate1 , mdatez(nnnchk+1) , nnnchk + 1
              call fatal(__FILE__,__LINE__,'ICBC date')
            else
            end if
!           print*,'UB1'
            do k = kx , 1 , -1
              mmrec = mmrec + 1
              read (iutbc,rec=mmrec) ((io2d(i,j),j=1,mjx),i=1,ix)
              do j = 1 , mjx
                do i = 1 , ix
                  ub1_io(i,k,j) = dble(io2d(i,j))
                end do
              end do
            end do
!           print*,'VB1'
            do k = kx , 1 , -1
              mmrec = mmrec + 1
              read (iutbc,rec=mmrec) ((io2d(i,j),j=1,mjx),i=1,ix)
              do j = 1 , mjx
                do i = 1 , ix
                  vb1_io(i,k,j) = dble(io2d(i,j))
                end do
              end do
            end do
!           print*,'TB1'
            do k = kx , 1 , -1
              mmrec = mmrec + 1
              read (iutbc,rec=mmrec) ((io2d(i,j),j=1,mjx),i=1,ix)
              do j = 1 , mjx
                do i = 1 , ix
                  tb1_io(i,k,j) = dble(io2d(i,j))
                end do
              end do
            end do
!           print*,'QB1'
            do k = kx , 1 , -1
              mmrec = mmrec + 1
              read (iutbc,rec=mmrec) ((io2d(i,j),j=1,mjx),i=1,ix)
              do j = 1 , mjx
                do i = 1 , ix
                  qb1_io(i,k,j) = dble(io2d(i,j))
                end do
              end do
            end do
!           print*,'PS1'
            mmrec = mmrec + 1
            read (iutbc,rec=mmrec) ((io2d(i,j),j=1,mjx),i=1,ix)
            do j = 1 , mjx
              do i = 1 , ix
                ps1_io(i,j) = dble(io2d(i,j))
              end do
            end do
!           print*,'TS1'
            mmrec = mmrec + 1
            read (iutbc,rec=mmrec) ((io2d(i,j),j=1,mjx),i=1,ix)
            do j = 1 , mjx
              do i = 1 , ix
                ts1_io(i,j) = dble(io2d(i,j))
              end do
            end do
            if ( ehso4 ) then
!             print*,'SO1'
              do k = kx , 1 , -1
                mmrec = mmrec + 1
                read (iutbc,rec=mmrec) ((io2d(i,j),j=1,mjx),i=1,ix)
                do j = 1 , mjx
                  do i = 1 , ix
                    so1_io(i,k,j) = dble(io2d(i,j))
                  end do
                end do
              end do
            end if
            if ( lsmtyp.eq.'USGS' ) mmrec = mmrec + 13
            do j = 1 , mjx
              do k = 1 , kx
                do i = 1 , ix
                  sav_0(i,k,j) = ub1_io(i,k,j)
                  sav_0(i,kx+k,j) = vb1_io(i,k,j)
                  sav_0(i,kx*2+k,j) = qb1_io(i,k,j)
                  sav_0(i,kx*3+k,j) = tb1_io(i,k,j)
                end do
              end do
              do i = 1 , ix
                sav_0(i,kx*4+1,j) = ps1_io(i,j)
                sav_0(i,kx*4+2,j) = ts1_io(i,j)
              end do
            end do
            if ( ehso4 ) then
              do j = 1 , mjx
                do k = 1 , kx
                  do i = 1 , ix
                    sav_0s(i,k,j) = so1_io(i,k,j)
                  end do
                end do
              end do
            end if
            exit
          end do
        end if
        call mpi_bcast(ndate1,1,mpi_integer,0,mpi_comm_world,ierr)
        call mpi_scatter(sav_0(1,1,1),ix*(kx*4+2)*jxp,                  &
                       & mpi_double_precision,sav0(1,1,1),ix*(kx*4+2)   &
                       & *jxp,mpi_double_precision,0,mpi_comm_world,    &
                       & ierr)
        do j = 1 , jendl
          do k = 1 , kx
            do i = 1 , ix
              ub1(i,k,j) = sav0(i,k,j)
              vb1(i,k,j) = sav0(i,kx+k,j)
              qb1(i,k,j) = sav0(i,kx*2+k,j)
              tb1(i,k,j) = sav0(i,kx*3+k,j)
            end do
          end do
          do i = 1 , ix
            ps1(i,j) = sav0(i,kx*4+1,j)
            ts1(i,j) = sav0(i,kx*4+2,j)
          end do
        end do
        if ( ehso4 ) then
          call mpi_scatter(sav_0s(1,1,1),ix*kx*jxp,mpi_double_precision,&
                         & sav0s(1,1,1),ix*kx*jxp,mpi_double_precision, &
                         & 0,mpi_comm_world,ierr)
          do j = 1 , jendl
            do k = 1 , kx
              do i = 1 , ix
                so1(i,k,j) = sav0s(i,k,j)
              end do
            end do
          end do
        end if
!       Convert surface pressure to pstar
        do j = 1 , jendl
          do i = 1 , ix
            ps1(i,j) = ps1(i,j) - ptop
          end do
        end do
!=======================================================================
!
!       this routine determines p(.) from p(x) by a 4-point
!       interpolation. on the x-grid, a p(x) point outside the grid
!       domain is assumed to satisfy p(0,j)=p(1,j); p(ix,j)=p(ix-1,j);
!       and similarly for the i's.
        call mpi_sendrecv(ps1(1,jxp),ix,mpi_double_precision,ieast,1,   &
                        & ps1(1,0),ix,mpi_double_precision,iwest,1,     &
                        & mpi_comm_world,status,ierr)
        do j = jbegin , jendx
          do i = 2 , ilx
            psdot(i,j) = 0.25*(ps1(i,j)+ps1(i-1,j)+ps1(i,j-1)+ps1(i-1,j-&
                       & 1))
          end do
        end do
!
        do i = 2 , ilx
          if ( myid.eq.0 ) psdot(i,1) = 0.5*(ps1(i,1)+ps1(i-1,1))
          if ( myid.eq.nproc-1 ) psdot(i,jendl)                         &
             & = 0.5*(ps1(i,jendx)+ps1(i-1,jendx))
        end do
!
        do j = jbegin , jendx
          psdot(1,j) = 0.5*(ps1(1,j)+ps1(1,j-1))
          psdot(ix,j) = 0.5*(ps1(ilx,j)+ps1(ilx,j-1))
        end do
!
        if ( myid.eq.0 ) then
          psdot(1,1) = ps1(1,1)
          psdot(ix,1) = ps1(ilx,1)
        end if
        if ( myid.eq.nproc-1 ) then
          psdot(1,jendl) = ps1(1,jendx)
          psdot(ix,jendl) = ps1(ilx,jendx)
        end if
!
!=======================================================================
!       Couple pressure u,v,t,q
        do k = 1 , kx
          do j = 1 , jendl
            do i = 1 , ix
              ub1(i,k,j) = ub1(i,k,j)*psdot(i,j)
              vb1(i,k,j) = vb1(i,k,j)*psdot(i,j)
              tb1(i,k,j) = tb1(i,k,j)*ps1(i,j)
              qb1(i,k,j) = qb1(i,k,j)*ps1(i,j)
            end do
          end do
        end do
 
        mdate = ndate0
        nnnchk = nnnchk + 1
 
!       print*,'read in datasets at :',ndate1
!
!-----compute boundary conditions for p*:
!
 
        if ( nspgx.le.jxp ) then
          nxwb = nspgx
        else
          nkk = nspgx/jxp
          if ( nspgx.eq.nkk*jxp ) then
            nxwb = jxp
          else
            nxwb = nspgx - nkk*jxp
          end if
        end if
        if ( nxwb+myid*jxp.gt.nspgx ) then
          nxwb = 0
        else if ( nxwb+myid*jxp.lt.nspgx ) then
          nxwb = jxp
        else
        end if
!
        if ( nspgx.le.jxp-1 ) then
          nxeb = nspgx
        else
          nkk = (nspgx-jxp+1)/jxp
          if ( (nspgx-jxp+1).eq.nkk*jxp ) then
            nxeb = jxp
          else
            nxeb = (nspgx-jxp+1) - nkk*jxp
          end if
        end if
        if ( mjx-1-(myid*jxp+jxp-nxeb).gt.nspgx ) then
          nxeb = 0
        else if ( mjx-1-(myid*jxp+jxp-nxeb).lt.nspgx ) then
          nxeb = min(jendx,jxp)
        else
        end if
        do nn = 1 , nxwb
          do i = 1 , ilx
            pwb(i,nn) = ps0(i,nn)
            pwbt(i,nn) = (ps1(i,nn)-ps0(i,nn))/dtbdys
          end do
        end do
        do nn = 1 , nxeb
          nnb = min(jendx,jxp) - nn + 1
          do i = 1 , ilx
            peb(i,nn) = ps0(i,nnb)
            pebt(i,nn) = (ps1(i,nnb)-ps0(i,nnb))/dtbdys
          end do
        end do
        do nn = 1 , nspgx
          nnb = ilx - nn + 1
          do j = 1 , jendx
            pnb(nn,j) = ps0(nnb,j)
            pss(nn,j) = ps0(nn,j)
            pnbt(nn,j) = (ps1(nnb,j)-ps0(nnb,j))/dtbdys
            psbt(nn,j) = (ps1(nn,j)-ps0(nn,j))/dtbdys
          end do
        end do
!
!-----compute boundary conditions for p*u and p*v:
!
        if ( nspgd.le.jxp ) then
          ndwb = nspgd
        else
          nkk = nspgd/jxp
          if ( nspgd.eq.nkk*jxp ) then
            ndwb = jxp
          else
            ndwb = nspgd - nkk*jxp
          end if
        end if
        if ( ndwb+myid*jxp.gt.nspgd ) then
          ndwb = 0
        else if ( ndwb+myid*jxp.lt.nspgd ) then
          ndwb = jxp
        else
        end if
!
        if ( nspgd.le.jendl ) then
          ndeb = nspgd
        else
          nkk = nspgd/jxp
          if ( nspgd.eq.nkk*jxp ) then
            ndeb = jxp
          else
            ndeb = nspgd - nkk*jxp
          end if
        end if
        if ( mjx-(myid*jxp+jxp-ndeb).gt.nspgd ) then
          ndeb = 0
        else if ( mjx-(myid*jxp+jxp-ndeb).lt.nspgd ) then
          ndeb = jxp
        else
        end if
        do nn = 1 , ndwb
          do k = 1 , kx
            do i = 1 , ix
              uwb(i,k,nn) = ub0(i,k,nn)
              vwb(i,k,nn) = vb0(i,k,nn)
              uwbt(i,k,nn) = (ub1(i,k,nn)-ub0(i,k,nn))/dtbdys
              vwbt(i,k,nn) = (vb1(i,k,nn)-vb0(i,k,nn))/dtbdys
            end do
          end do
        end do
        do nn = 1 , ndeb
          nnb = min(jendl,jxp) - nn + 1
          do k = 1 , kx
            do i = 1 , ix
              ueb(i,k,nn) = ub0(i,k,nnb)
              veb(i,k,nn) = vb0(i,k,nnb)
              uebt(i,k,nn) = (ub1(i,k,nnb)-ub0(i,k,nnb))/dtbdys
              vebt(i,k,nn) = (vb1(i,k,nnb)-vb0(i,k,nnb))/dtbdys
            end do
          end do
        end do
        do nn = 1 , nspgd
          nnb = ix - nn + 1
          do k = 1 , kx
            do j = 1 , jendl
              unb(nn,k,j) = ub0(nnb,k,j)
              usb(nn,k,j) = ub0(nn,k,j)
              vnb(nn,k,j) = vb0(nnb,k,j)
              vsb(nn,k,j) = vb0(nn,k,j)
              unbt(nn,k,j) = (ub1(nnb,k,j)-ub0(nnb,k,j))/dtbdys
              usbt(nn,k,j) = (ub1(nn,k,j)-ub0(nn,k,j))/dtbdys
              vnbt(nn,k,j) = (vb1(nnb,k,j)-vb0(nnb,k,j))/dtbdys
              vsbt(nn,k,j) = (vb1(nn,k,j)-vb0(nn,k,j))/dtbdys
            end do
          end do
        end do
!
!-----compute boundary conditions for p*t and p*qv:
!
        do nn = 1 , nxwb
          do k = 1 , kx
            do i = 1 , ilx
              twb(i,k,nn) = tb0(i,k,nn)
              qwb(i,k,nn) = qb0(i,k,nn)
              twbt(i,k,nn) = (tb1(i,k,nn)-tb0(i,k,nn))/dtbdys
              qwbt(i,k,nn) = (qb1(i,k,nn)-qb0(i,k,nn))/dtbdys
            end do
          end do
        end do
        do nn = 1 , nxeb
          nnb = min(jendx,jxp) - nn + 1
          do k = 1 , kx
            do i = 1 , ilx
              teb(i,k,nn) = tb0(i,k,nnb)
              qeb(i,k,nn) = qb0(i,k,nnb)
              tebt(i,k,nn) = (tb1(i,k,nnb)-tb0(i,k,nnb))/dtbdys
              qebt(i,k,nn) = (qb1(i,k,nnb)-qb0(i,k,nnb))/dtbdys
            end do
          end do
        end do
        do nn = 1 , nspgx
          nnb = ilx - nn + 1
          do k = 1 , kx
            do j = 1 , jendx
              tnb(nn,k,j) = tb0(nnb,k,j)
              tsb(nn,k,j) = tb0(nn,k,j)
              qnb(nn,k,j) = qb0(nnb,k,j)
              qsb(nn,k,j) = qb0(nn,k,j)
              tnbt(nn,k,j) = (tb1(nnb,k,j)-tb0(nnb,k,j))/dtbdys
              tsbt(nn,k,j) = (tb1(nn,k,j)-tb0(nn,k,j))/dtbdys
              qnbt(nn,k,j) = (qb1(nnb,k,j)-qb0(nnb,k,j))/dtbdys
              qsbt(nn,k,j) = (qb1(nn,k,j)-qb0(nn,k,j))/dtbdys
            end do
          end do
        end do
        if ( myid.eq.0 ) print * , 'BCs are ready from ' , ndate0 ,     &
                              &'  to ' , ndate1
        idatex = ndate0
        ndate0 = ndate1
        do j = 1 , jendx
          do i = 1 , ilx
            tdum(i,j) = ts1(i,j)
          end do
        end do
        do k = 1 , kx
          do j = 1 , jendl
            do i = 1 , ix
              ub0(i,k,j) = ub1(i,k,j)
              vb0(i,k,j) = vb1(i,k,j)
              qb0(i,k,j) = qb1(i,k,j)
              tb0(i,k,j) = tb1(i,k,j)
            end do
          end do
        end do
        do j = 1 , jendl
          do i = 1 , ix
            ps0(i,j) = ps1(i,j)
            ts0(i,j) = ts1(i,j)
          end do
        end do
        if ( ehso4 ) then
          do k = 1 , kx
            do j = 1 , jendl
              do i = 1 , ix
                so0(i,k,j) = so1(i,k,j)
              end do
            end do
          end do
        end if
!bxqOCT2001_
 
        nnbase = nnnnnn
        nyear = mdate/1000000
        nmonth = (mdate-nyear*1000000)/10000
 
!-----------------------------------------------------------------------
        if ( ldatez.lt.ndate1 ) then
 
          do j = 1 , jendx
            do i = 1 , ix - 1
              if ( veg2d(i,j).le.0.00001 ) then
                tga(i,j) = tdum(i,j)
                tgb(i,j) = tdum(i,j)
!elguindi       if (tdum(i,j).le.271.35) then
!               if (tdum(i,j).le.271.38) then
!               print *,'Setting ocld2d to ice at i=',i,' j=',j,'
!               t=',tdum(i,j) do n=1,NNSG
!               ocld2d(n,i,j)=2.
!               sice2d(n,i,j)=1000.
!               end do
!               else
!               do n=1,NNSG
!               ocld2d(n,i,j)=0.
!               sice2d(n,i,j)=0.
!               end do
!               end if
              end if
            end do
          end do
          exit
        end if
      end do
#else
 100  continue
      dtbdys = ibdyfrq*60.*60.
      do
        if ( ehso4 ) then
          do k = 1 , kx
            do j = 1 , jx
              do i = 1 , ix
                so4(i,k,j) = so0(i,k,j)
              end do
            end do
          end do
        end if
        mmrec = mmrec + 1
        read (iutbc,rec=mmrec,iostat=ierr1) ndate1
        if ( ierr1.ne.0 ) then
          close (iutbc)
          iutbc = iutbc + 1
          write (finm,99001) iutbc
!_sgi     101    format('fort.',I2)
          open (iutbc,file=finm,form='unformatted',status='old',        &
               &access='direct',recl=ix*jx*ibyte)
          mmrec = 0
          print * , 'CHANGING BDY UNIT NUMBER:  iutbc=' , iutbc
          if ( iutbc.gt.999 )                                           &
              & call fatal(__FILE__,__LINE__,'BDY UNIT MAX EXCEEDED')
          cycle
        end if
        if ( ndate1.lt.mdatez(nnnchk+1) ) then
          if ( ndate1.lt.mdatez(nnnchk+1) ) then
            print * , 'Searching for proper date: ' , ndate1 ,          &
                & mdatez(nnnchk+1)
            print * , 'read in datasets at :' , ndate0
            if ( ehso4 ) then
              if ( lsmtyp.ne.'USGS' ) then
                mmrec = mmrec + kx*5 + 2
              else
                mmrec = mmrec + kx*5 + 2 + 13
              end if
            else if ( lsmtyp.ne.'USGS' ) then
              mmrec = mmrec + kx*4 + 2
            else
              mmrec = mmrec + kx*4 + 2 + 13
            end if
            cycle
          end if
        else if ( ndate1.gt.mdatez(nnnchk+1) ) then
          print * , 'DATE IN BC FILE EXCEEDED DATE IN RegCM'
          print * , ndate1 , mdatez(nnnchk+1) , nnnchk + 1
          call fatal(__FILE__,__LINE__,'ICBC date')
        else
        end if
!       print*,'UB1'
        do k = kx , 1 , -1
          mmrec = mmrec + 1
          read (iutbc,rec=mmrec) ((io2d(i,j),j=1,jx),i=1,ix)
          do j = 1 , jx
            do i = 1 , ix
              ub1(i,k,j) = dble(io2d(i,j))
            end do
          end do
        end do
!       print*,'VB1'
        do k = kx , 1 , -1
          mmrec = mmrec + 1
          read (iutbc,rec=mmrec) ((io2d(i,j),j=1,jx),i=1,ix)
          do j = 1 , jx
            do i = 1 , ix
              vb1(i,k,j) = dble(io2d(i,j))
            end do
          end do
        end do
!       print*,'TB1'
        do k = kx , 1 , -1
          mmrec = mmrec + 1
          read (iutbc,rec=mmrec) ((io2d(i,j),j=1,jx),i=1,ix)
          do j = 1 , jx
            do i = 1 , ix
              tb1(i,k,j) = dble(io2d(i,j))
            end do
          end do
        end do
!       print*,'QB1'
        do k = kx , 1 , -1
          mmrec = mmrec + 1
          read (iutbc,rec=mmrec) ((io2d(i,j),j=1,jx),i=1,ix)
          do j = 1 , jx
            do i = 1 , ix
              qb1(i,k,j) = dble(io2d(i,j))
            end do
          end do
        end do
!       print*,'PS1'
        mmrec = mmrec + 1
        read (iutbc,rec=mmrec) ((io2d(i,j),j=1,jx),i=1,ix)
        do j = 1 , jx
          do i = 1 , ix
            ps1(i,j) = dble(io2d(i,j))
          end do
        end do
!       print*,'TS1'
        mmrec = mmrec + 1
        read (iutbc,rec=mmrec) ((io2d(i,j),j=1,jx),i=1,ix)
        do j = 1 , jx
          do i = 1 , ix
            ts1(i,j) = dble(io2d(i,j))
          end do
        end do
        if ( ehso4 ) then
!         print*,'SO1'
          do k = kx , 1 , -1
            mmrec = mmrec + 1
            read (iutbc,rec=mmrec) ((io2d(i,j),j=1,jx),i=1,ix)
            do j = 1 , jx
              do i = 1 , ix
                so1(i,k,j) = dble(io2d(i,j))
              end do
            end do
          end do
        end if
        if ( lsmtyp.eq.'USGS' ) mmrec = mmrec + 13
!       Convert surface pressure to pstar
        do j = 1 , jx
          do i = 1 , ix
            ps1(i,j) = ps1(i,j) - ptop
          end do
        end do
!=======================================================================
!
!       this routine determines p(.) from p(x) by a 4-point
!       interpolation. on the x-grid, a p(x) point outside the grid
!       domain is assumed to satisfy p(0,j)=p(1,j); p(ix,j)=p(ix-1,j);
!       and similarly for the i's.
        do j = 2 , jlx
          do i = 2 , ilx
            psdot(i,j) = 0.25*(ps1(i,j)+ps1(i-1,j)+ps1(i,j-1)+ps1(i-1,j-&
                       & 1))
          end do
        end do
!
        do i = 2 , ilx
          psdot(i,1) = 0.5*(ps1(i,1)+ps1(i-1,1))
          psdot(i,jx) = 0.5*(ps1(i,jlx)+ps1(i-1,jlx))
        end do
!
        do j = 2 , jlx
          psdot(1,j) = 0.5*(ps1(1,j)+ps1(1,j-1))
          psdot(ix,j) = 0.5*(ps1(ilx,j)+ps1(ilx,j-1))
        end do
!
        psdot(1,1) = ps1(1,1)
        psdot(ix,1) = ps1(ilx,1)
        psdot(1,jx) = ps1(1,jlx)
        psdot(ix,jx) = ps1(ilx,jlx)
!
!=======================================================================
!       Couple pressure u,v,t,q
        do k = 1 , kx
          do j = 1 , jx
            do i = 1 , ix
              ub1(i,k,j) = ub1(i,k,j)*psdot(i,j)
              vb1(i,k,j) = vb1(i,k,j)*psdot(i,j)
              tb1(i,k,j) = tb1(i,k,j)*ps1(i,j)
              qb1(i,k,j) = qb1(i,k,j)*ps1(i,j)
            end do
          end do
        end do
 
        mdate = ndate0
        nnnchk = nnnchk + 1
 
!       print*,'read in datasets at :',ndate1
!
!-----compute boundary conditions for p*:
!
 
        do nn = 1 , nspgx
          do i = 1 , ilx
            pwb(i,nn) = ps0(i,nn)
            pwbt(i,nn) = (ps1(i,nn)-ps0(i,nn))/dtbdys
          end do
        end do
        do nn = 1 , nspgx
          nnb = jlx - nn + 1
          do i = 1 , ilx
            peb(i,nn) = ps0(i,nnb)
            pebt(i,nn) = (ps1(i,nnb)-ps0(i,nnb))/dtbdys
          end do
        end do
        do nn = 1 , nspgx
          nnb = ilx - nn + 1
          do j = 1 , jlx
            pnb(nn,j) = ps0(nnb,j)
            pss(nn,j) = ps0(nn,j)
            pnbt(nn,j) = (ps1(nnb,j)-ps0(nnb,j))/dtbdys
            psbt(nn,j) = (ps1(nn,j)-ps0(nn,j))/dtbdys
          end do
        end do
!
!-----compute boundary conditions for p*u and p*v:
!
        do nn = 1 , nspgd
          do k = 1 , kx
            do i = 1 , ix
              uwb(i,k,nn) = ub0(i,k,nn)
              vwb(i,k,nn) = vb0(i,k,nn)
              uwbt(i,k,nn) = (ub1(i,k,nn)-ub0(i,k,nn))/dtbdys
              vwbt(i,k,nn) = (vb1(i,k,nn)-vb0(i,k,nn))/dtbdys
            end do
          end do
        end do
        do nn = 1 , nspgd
          nnb = jx - nn + 1
          do k = 1 , kx
            do i = 1 , ix
              ueb(i,k,nn) = ub0(i,k,nnb)
              veb(i,k,nn) = vb0(i,k,nnb)
              uebt(i,k,nn) = (ub1(i,k,nnb)-ub0(i,k,nnb))/dtbdys
              vebt(i,k,nn) = (vb1(i,k,nnb)-vb0(i,k,nnb))/dtbdys
            end do
          end do
        end do
        do nn = 1 , nspgd
          nnb = ix - nn + 1
          do k = 1 , kx
            do j = 1 , jx
              unb(nn,k,j) = ub0(nnb,k,j)
              usb(nn,k,j) = ub0(nn,k,j)
              vnb(nn,k,j) = vb0(nnb,k,j)
              vsb(nn,k,j) = vb0(nn,k,j)
              unbt(nn,k,j) = (ub1(nnb,k,j)-ub0(nnb,k,j))/dtbdys
              usbt(nn,k,j) = (ub1(nn,k,j)-ub0(nn,k,j))/dtbdys
              vnbt(nn,k,j) = (vb1(nnb,k,j)-vb0(nnb,k,j))/dtbdys
              vsbt(nn,k,j) = (vb1(nn,k,j)-vb0(nn,k,j))/dtbdys
            end do
          end do
        end do
!
!-----compute boundary conditions for p*t and p*qv:
!
        do nn = 1 , nspgx
          do k = 1 , kx
            do i = 1 , ilx
              twb(i,k,nn) = tb0(i,k,nn)
              qwb(i,k,nn) = qb0(i,k,nn)
              twbt(i,k,nn) = (tb1(i,k,nn)-tb0(i,k,nn))/dtbdys
              qwbt(i,k,nn) = (qb1(i,k,nn)-qb0(i,k,nn))/dtbdys
            end do
          end do
        end do
        do nn = 1 , nspgx
          nnb = jlx - nn + 1
          do k = 1 , kx
            do i = 1 , ilx
              teb(i,k,nn) = tb0(i,k,nnb)
              qeb(i,k,nn) = qb0(i,k,nnb)
              tebt(i,k,nn) = (tb1(i,k,nnb)-tb0(i,k,nnb))/dtbdys
              qebt(i,k,nn) = (qb1(i,k,nnb)-qb0(i,k,nnb))/dtbdys
            end do
          end do
        end do
        do nn = 1 , nspgx
          nnb = ilx - nn + 1
          do k = 1 , kx
            do j = 1 , jlx
              tnb(nn,k,j) = tb0(nnb,k,j)
              tsb(nn,k,j) = tb0(nn,k,j)
              qnb(nn,k,j) = qb0(nnb,k,j)
              qsb(nn,k,j) = qb0(nn,k,j)
              tnbt(nn,k,j) = (tb1(nnb,k,j)-tb0(nnb,k,j))/dtbdys
              tsbt(nn,k,j) = (tb1(nn,k,j)-tb0(nn,k,j))/dtbdys
              qnbt(nn,k,j) = (qb1(nnb,k,j)-qb0(nnb,k,j))/dtbdys
              qsbt(nn,k,j) = (qb1(nn,k,j)-qb0(nn,k,j))/dtbdys
            end do
          end do
        end do
        print * , 'BCs are ready from ' , ndate0 , '  to ' , ndate1
        idatex = ndate0
        ndate0 = ndate1
        do j = 1 , jlx
          do i = 1 , ilx
            tdum(i,j) = ts1(i,j)
          end do
        end do
        do k = 1 , kx
          do j = 1 , jx
            do i = 1 , ix
              ub0(i,k,j) = ub1(i,k,j)
              vb0(i,k,j) = vb1(i,k,j)
              qb0(i,k,j) = qb1(i,k,j)
              tb0(i,k,j) = tb1(i,k,j)
            end do
          end do
        end do
        do j = 1 , jx
          do i = 1 , ix
            ps0(i,j) = ps1(i,j)
            ts0(i,j) = ts1(i,j)
          end do
        end do
        if ( ehso4 ) then
          do k = 1 , kx
            do j = 1 , jx
              do i = 1 , ix
                so0(i,k,j) = so1(i,k,j)
              end do
            end do
          end do
        end if
!bxqOCT2001_
 
        nnbase = nnnnnn
        nyear = mdate/1000000
        nmonth = (mdate-nyear*1000000)/10000
 
!-----------------------------------------------------------------------
        if ( ldatez.ge.ndate1 ) go to 100
 
        do j = 1 , jx - 1
          do i = 1 , ix - 1
            if ( veg2d(i,j).le.0.00001 ) then
              tga(i,j) = tdum(i,j)
              tgb(i,j) = tdum(i,j)
!elguindi     if (tdum(i,j).le.271.35) then
!             if (tdum(i,j).le.271.38) then
!             print *,'Setting ocld2d to ice at i=',i,' j=',j,'
!             t=',tdum(i,j) do n=1,NNSG
!             ocld2d(n,i,j)=2.
!             sice2d(n,i,j)=1000.
!             end do
!             else
!             do n=1,NNSG
!             ocld2d(n,i,j)=0.
!             sice2d(n,i,j)=0.
!             end do
!             end if
            end if
          end do
        end do
        exit
      end do
#endif

99001 format ('fort.',i3)
 
      end subroutine bdyin
