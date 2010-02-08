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
 
      subroutine chsrfem

      use regcm_param
      use mainchem
      use trachem
      use dust
      use iunits
      use aerosol
      use message
#ifdef MPP1
      use mpi
#endif
      implicit none
!
! Local variables
!
      character(5) :: aerctl
      integer :: i , ierr , itr , j , k , l , m , n
      logical :: rd_tex , there
#ifdef MPP1
      real(4) , dimension(ix,mjx) :: toto
#else
      real(4) , dimension(ix,jx) :: toto
#endif
!
      rd_tex = .false.
#ifdef MPP1
      if ( myid.eq.0 ) then
        do itr = 1 , ntr
          aerctl = chtrname(itr)
          if ( aerctl(1:4).eq.'DUST' ) then
            rd_tex = .true.
            exit
          end if
        end do
      end if
#else
      do itr = 1 , ntr
        aerctl = chtrname(itr)
        if ( aerctl(1:4).eq.'DUST' ) then
          rd_tex = .true.
          exit
        end if
      end do
#endif

#ifdef MPP1
      call mpi_bcast(rd_tex,1,mpi_logical,0,mpi_comm_world,ierr)
      if ( myid.eq.0 ) then
        do itr = 1 , ntr
          aerctl = chtrname(itr)
          print * , itr , aerctl
          if ( aerctl(1:4).ne.'DUST' .and. aertyp(4:5).ne.'00' ) then
            open (unit=iutchsrc,file='AERO.dat',status='old',           &
                 &form='unformatted',access='direct',recl=ix*mjx*ibyte)
            if ( aerctl(1:3).eq.'SO2' ) then
              if ( aertyp(4:4).eq.'1' ) then
                read (iutchsrc,rec=1) ((toto(i,j),j=1,mjx),i=1,ix)
                do m = 1 , 12
                  do j = 1 , mjx
                    do i = 1 , ix
                      chemsrc_io(i,j,m,itr) = toto(i,j)
                    end do
                  end do
                end do
              else
                do m = 1 , 12
                  do j = 1 , mjx
                    do i = 1 , ix
                      chemsrc_io(i,j,m,itr) = 0.0D0
                    end do
                  end do
                end do
              end if
              if ( aertyp(5:5).eq.'1' ) then
                do m = 1 , 12
                  read (iutchsrc,rec=3+m) ((toto(i,j),j=1,mjx),i=1,ix)
                  do j = 1 , mjx
                    do i = 1 , ix
                      chemsrc_io(i,j,m,itr) = chemsrc_io(i,j,m,itr)     &
                      & + toto(i,j)
                    end do
                  end do
                end do
              end if
            else if ( aerctl(1:2).eq.'BC' ) then
              if ( aertyp(4:4).eq.'1' ) then
                read (iutchsrc,rec=2) ((toto(i,j),j=1,mjx),i=1,ix)
                do m = 1 , 12
                  do j = 1 , mjx
                    do i = 1 , ix
                      chemsrc_io(i,j,m,itr) = toto(i,j)
                    end do
                  end do
                end do
              else
                do m = 1 , 12
                  do j = 1 , mjx
                    do i = 1 , ix
                      chemsrc_io(i,j,m,itr) = 0.0D0
                    end do
                  end do
                end do
              end if
              if ( aertyp(5:5).eq.'1' ) then
                do m = 1 , 12
                  read (iutchsrc,rec=15+m) ((toto(i,j),j=1,mjx),i=1,ix)
                  do j = 1 , mjx
                    do i = 1 , ix
                      chemsrc_io(i,j,m,itr) = chemsrc_io(i,j,m,itr)     &
                      & + toto(i,j)
                    end do
                  end do
                end do
              end if
            else if ( aerctl(1:2).eq.'OC' ) then
              if ( aertyp(4:4).eq.'1' ) then
                read (iutchsrc,rec=3) ((toto(i,j),j=1,mjx),i=1,ix)
                do m = 1 , 12
                  do j = 1 , mjx
                    do i = 1 , ix
                      chemsrc_io(i,j,m,itr) = toto(i,j)
                    end do
                  end do
                end do
              else
                do m = 1 , 12
                  do j = 1 , mjx
                    do i = 1 , ix
                      chemsrc_io(i,j,m,itr) = 0.0D0
                    end do
                  end do
                end do
              end if
              if ( aertyp(5:5).eq.'1' ) then
                do m = 1 , 12
                  read (iutchsrc,rec=27+m) ((toto(i,j),j=1,mjx),i=1,ix)
                  do j = 1 , mjx
                    do i = 1 , ix
                      chemsrc_io(i,j,m,itr) = chemsrc_io(i,j,m,itr)     &
                      & + toto(i,j)
                    end do
                  end do
                end do
              end if
            else
            end if
          end if
        end do
        do j = 1 , mjx
          do itr = 1 , ntr
            do m = 1 , 12
              do i = 1 , ix
                src_0(i,m,itr,j) = chemsrc_io(i,j,m,itr)
              end do
            end do
          end do
        end do
      end if
!     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call mpi_scatter(src_0(1,1,1,1),ix*12*ntr*jxp,                    &
                     & mpi_double_precision,src0(1,1,1,1),ix*12*ntr*jxp,&
                     & mpi_double_precision,0,mpi_comm_world,ierr)
      do j = 1 , jendl
        do itr = 1 , ntr
          do m = 1 , 12
            do i = 1 , ix
              chemsrc(i,j,m,itr) = src0(i,m,itr,j)
            end do
          end do
        end do
      end do
#else
      do itr = 1 , ntr
        aerctl = chtrname(itr)
        print * , itr , aerctl
        if ( aerctl(1:4).ne.'DUST' .and. aertyp(4:5).ne.'00' ) then
          open (unit=iutchsrc,file='AERO.dat',status='old',             &
               &form='unformatted',access='direct',recl=ix*jx*ibyte)
          if ( aerctl(1:3).eq.'SO2' ) then
            if ( aertyp(4:4).eq.'1' ) then
              read (iutchsrc,rec=1) ((toto(i,j),j=1,jx),i=1,ix)
              do m = 1 , 12
                do j = 1 , jx
                  do i = 1 , ix
                    chemsrc(i,j,m,itr) = toto(i,j)
                  end do
                end do
              end do
            else
              do m = 1 , 12
                do j = 1 , jx
                  do i = 1 , ix
                    chemsrc(i,j,m,itr) = 0.0D0
                  end do
                end do
              end do
            end if
            if ( aertyp(5:5).eq.'1' ) then
              do m = 1 , 12
                read (iutchsrc,rec=3+m) ((toto(i,j),j=1,jx),i=1,ix)
                do j = 1 , jx
                  do i = 1 , ix
                    chemsrc(i,j,m,itr) = chemsrc(i,j,m,itr) + toto(i,j)
                  end do
                end do
              end do
            end if
          else if ( aerctl(1:2).eq.'BC' ) then
            if ( aertyp(4:4).eq.'1' ) then
              read (iutchsrc,rec=2) ((toto(i,j),j=1,jx),i=1,ix)
              do m = 1 , 12
                do j = 1 , jx
                  do i = 1 , ix
                    chemsrc(i,j,m,itr) = toto(i,j)
                  end do
                end do
              end do
            else
              do m = 1 , 12
                do j = 1 , jx
                  do i = 1 , ix
                    chemsrc(i,j,m,itr) = 0.0D0
                  end do
                end do
              end do
            end if
            if ( aertyp(5:5).eq.'1' ) then
              do m = 1 , 12
                read (iutchsrc,rec=15+m) ((toto(i,j),j=1,jx),i=1,ix)
                do j = 1 , jx
                  do i = 1 , ix
                    chemsrc(i,j,m,itr) = chemsrc(i,j,m,itr) + toto(i,j)
                  end do
                end do
              end do
            end if
          else if ( aerctl(1:2).eq.'OC' ) then
            if ( aertyp(4:4).eq.'1' ) then
              read (iutchsrc,rec=3) ((toto(i,j),j=1,jx),i=1,ix)
              do m = 1 , 12
                do j = 1 , jx
                  do i = 1 , ix
                    chemsrc(i,j,m,itr) = toto(i,j)
                  end do
                end do
              end do
            else
              do m = 1 , 12
                do j = 1 , jx
                  do i = 1 , ix
                    chemsrc(i,j,m,itr) = 0.0D0
                  end do
                end do
              end do
            end if
            if ( aertyp(5:5).eq.'1' ) then
              do m = 1 , 12
                read (iutchsrc,rec=27+m) ((toto(i,j),j=1,jx),i=1,ix)
                do j = 1 , jx
                  do i = 1 , ix
                    chemsrc(i,j,m,itr) = chemsrc(i,j,m,itr) + toto(i,j)
                  end do
                end do
              end do
            end if
          else
          end if
        end if
      end do
#endif
 
!     modification dust : read the soil texture type

#ifdef MPP1
      if ( myid.eq.0 ) then
        if ( rd_tex ) then
          if ( lsmtyp.eq.'BATS' ) then
            read (iutin,rec=14) ((toto(i,j),j=1,mjx),i=1,ix)
          else if ( lsmtyp.eq.'USGS' ) then
            read (iutin,rec=43) ((toto(i,j),j=1,mjx),i=1,ix)
          else
          end if
          close (iutin)
          do j = 1 , mjx
            do i = 1 , ix
              dustsotex_io(i,j) = dble(toto(i,j))
            end do
          end do
        end if
      end if
!     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call mpi_scatter(dustsotex_io(1,1),ix*jxp,mpi_double_precision,   &
                     & dustsotex(1,1),ix*jxp,mpi_double_precision,0,    &
                     & mpi_comm_world,ierr)
#else
      if ( rd_tex ) then
        if ( lsmtyp.eq.'BATS' ) then
          read (iutin,rec=14) ((toto(i,j),j=1,jx),i=1,ix)
        else if ( lsmtyp.eq.'USGS' ) then
          read (iutin,rec=43) ((toto(i,j),j=1,jx),i=1,ix)
        else
        end if
        close (iutin)
        do j = 1 , jx
          do i = 1 , ix
            dustsotex(i,j) = dble(toto(i,j))
          end do
        end do
      end if
#endif
      if ( rd_tex ) call inidust
 
!     sulfates sources
      do m = 1 , 12
#ifdef MPP1
        do j = 1 , jendl
#else
        do j = 1 , jx
#endif
          do i = 1 , ix
            if ( iso4.gt.0 ) chemsrc(i,j,m,iso4)                        &
               & = 0.02*chemsrc(i,j,m,iso2)
            if ( iso2.gt.0 ) chemsrc(i,j,m,iso2)                        &
               & = 0.98*chemsrc(i,j,m,iso2)
 
!           partition hydrophilic hydrophonic ( cooke et al.1999)
!           BC
            if ( ibchb.gt.0 .and. ibchl.gt.0 ) then
              chemsrc(i,j,m,ibchl) = 0.2*chemsrc(i,j,m,ibchb)
              chemsrc(i,j,m,ibchb) = 0.8*chemsrc(i,j,m,ibchb)
            end if
!           OC
            if ( iochb.gt.0 .and. iochl.gt.0 ) then
              chemsrc(i,j,m,iochl) = 0.5*chemsrc(i,j,m,iochb)
              chemsrc(i,j,m,iochb) = 0.5*chemsrc(i,j,m,iochb)
            end if
          end do
        end do
      end do
 
!     OPtical properties / internal mixing
      if ( mixtype.eq.2 ) then
#ifdef MPP1
        if ( myid.eq.0 ) then
          inquire (file='optdat.bin',exist=there)
          if ( .not.there ) then
            write (*,*) 'For mixtype=2, optdat.bin is required'
            write (*,*) 'ln -s ../Main/Commons/optdat.bin optdat.bin'
            call fatal(__FILE__,__LINE__,'optdat.bin is required')
          end if
          open (iutopt,file='optdat.bin',form='unformatted',            &
              & recl=4*19*11*11*11*11*ibyte,access='direct')
          read (iutopt,rec=1) ((((((dextmix(i,j,k,l,m,n),i=1,4),j=1,19),&
                            & k=1,11),l=1,11),m=1,11),n=1,11)
          read (iutopt,rec=2) ((((((dssamix(i,j,k,l,m,n),i=1,4),j=1,19),&
                            & k=1,11),l=1,11),m=1,11),n=1,11)
          read (iutopt,rec=3) ((((((dgmix(i,j,k,l,m,n),i=1,4),j=1,19),k=&
                            & 1,11),l=1,11),m=1,11),n=1,11)
          close (iutopt)
        end if
        call mpi_bcast(dextmix,4*19*11*11*11*11,mpi_real,0,             &
                     & mpi_comm_world,ierr)
        call mpi_bcast(dssamix,4*19*11*11*11*11,mpi_real,0,             &
                     & mpi_comm_world,ierr)
        call mpi_bcast(dgmix,4*19*11*11*11*11,mpi_real,0,mpi_comm_world,&
                     & ierr)
#else
        inquire (file='optdat.bin',exist=there)
        if ( .not.there ) then
          write (*,*) 'For mixtype=2, optdat.bin is required'
          write (*,*) 'ln -s ../Main/Commons/optdat.bin optdat.bin'
          call fatal(__FILE__,__LINE__,'optdat.bin is required')
        end if
        open (iutopt,file='optdat.bin',form='unformatted',              &
            & recl=4*19*11*11*11*11*ibyte,access='direct')
        read (iutopt,rec=1) ((((((dextmix(i,j,k,l,m,n),i=1,4),j=1,19),k=&
                          & 1,11),l=1,11),m=1,11),n=1,11)
        read (iutopt,rec=2) ((((((dssamix(i,j,k,l,m,n),i=1,4),j=1,19),k=&
                          & 1,11),l=1,11),m=1,11),n=1,11)
        read (iutopt,rec=3) ((((((dgmix(i,j,k,l,m,n),i=1,4),j=1,19),k=1,&
                          & 11),l=1,11),m=1,11),n=1,11)
        close (iutopt)
#endif
 
!       Check !
 
        do k = 1 , 11
          do l = 1 , 11
            do m = 1 , 11
              do n = 1 , 11
 
                if ( k+l+m+n.eq.14 ) then
                  do i = 1 , 4
                    do j = 1 , 19
 
                      if ( (dextmix(i,j,k,l,m,n).lt.0.) .or.            &
                         & (dextmix(i,j,k,l,m,n).gt.20.) ) then
                        write (aline,*) 'problem in dextmix ' ,         &
                                      & dextmix(i,j,k,l,m,n)
                        call say
                        call fatal(__FILE__,__LINE__,'DETMIX ERROR')
                      end if
 
                      if ( (dssamix(i,j,k,l,m,n).lt.0.) .or.            &
                         & (dssamix(i,j,k,l,m,n).gt.1.) ) then
                        write (aline,*) 'problem in dssamix ' ,         &
                                      & dssamix(i,j,k,l,m,n)
                        call say
                        call fatal(__FILE__,__LINE__,'DSSAMIX ERROR')
                      end if
 
                      if ( (dgmix(i,j,k,l,m,n).lt.0.) .or.              &
                         & (dgmix(i,j,k,l,m,n).gt.1.) ) then
                        write (aline,*) 'problem in dgmix ' ,           &
                                      & dgmix(i,j,k,l,m,n)
                        call say
                        call fatal(__FILE__,__LINE__,'DGMIX ERROR')
                      end if
 
                    end do
                  end do
                end if
 
              end do
            end do
          end do
        end do
 
#ifdef MPP1
        if ( myid.eq.0 ) write (*,*) '! OPDATA CHECKED !'
#else
        write (*,*) '! OPDATA CHECKED !'
#endif
      end if
      end subroutine chsrfem
