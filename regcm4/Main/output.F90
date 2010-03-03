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
 
      subroutine output(iexec)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!     this subroutine handles all of the output                       c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use mod_regcm_param
      use mod_param1
      use mod_param2
      use mod_iunits
      use mod_main
      use mod_mainchem
      use mod_bats
      use mod_date
      use mod_message
      use mod_bdycod
      use mod_pmoist
      use mod_rad
      use mod_trachem
      use mod_cvaria
      use mod_outrad
      use mod_tmpsav
      use mod_radbuf
      use mod_split
#ifdef MPP1
      use mpi
      use mod_mppio
#endif
      implicit none
!
! Dummy arguments
!
      integer :: iexec
      intent (in) iexec
!
! Local variables
!
      integer :: i , iresult , j
      logical :: there
      character(3) :: itype
      character(14) :: newfil
      character(17) :: tmpfil
#ifdef MPP1
      integer :: allrec , idum , ierr , l , k , n
      real(8) , dimension(ix,kx*6+3+nnsg*4,jxp) :: atm0
      real(8) , dimension(ix,kx*6+3+nnsg*4,jx) :: atm_0
      real(4) , dimension(ixm2,numbat,jxp) :: bat0
      real(4) , dimension(ixm2,numbat,jx) :: bat_0
      real(8) , dimension(ix,3,jxp) :: out0
      real(8) , dimension(ix,3,jx) :: out_0
      real(4) , dimension(ixm2,nrad3d*kx+nrad2d,jxp) :: rad0
      real(4) , dimension(ixm2,nrad3d*kx+nrad2d,jx) :: rad_0
      real(4) , dimension(ixm2,nnsg,numsub,jxp) :: sub0
      real(8) , dimension(ixm2,nnsg,numsub,jx) :: sub_0
      real(8) , dimension(ix,ntr*kx+kx*3+ntr*7+3,jxp) :: chem0
      real(8) , dimension(ix,ntr*kx+kx*3+ntr*7+3,jx) :: chem_0
#endif
!
!----------------------------------------------------------------------
!
      if ( jyear.ne.jyear0 .or. ktau.ne.0 ) then
        if ( mod(nint(xtime),60).lt.mod(nint(xtime-dtmin),60) )         &
           & idatex = idatex + 1
        if ( dabs(xtime).lt.0.00001 ) idatex = ldatez
      end if
 
#ifdef MPP1
      if ( myid.eq.0 ) then
        if ( jyear.eq.jyearr .and. ktau.eq.ktaur ) then
          if ( iotyp.eq.1 ) then
            print * , 'Writing output files in direct access format'
          else if ( iotyp.eq.2 ) then
            print * , 'Writing output files in sequential format'
          else
            write (aline,*) 'iotyp = ' , iotyp
            call say
            call fatal(__FILE__,__LINE__,'Output format does not exist')
          end if
          nrcout = 0
          call mkfile
 
        end if
      end if
 
      if ( jyear.eq.jyear0 .and. ktau.eq.0 ) then
!=======================================================================
!       among    ht,htsd,veg2d,satbrt,xlat,xlong,msfx,msfd,f
!       we just need gather                      msfx,msfd
        do j = 1 , jendl
          do i = 1 , ixm1
            out0(i,1,j) = veg2d(i,j)
          end do
        end do
        do j = 1 , jendl
          do i = 1 , ix
            out0(i,2,j) = msfx(i,j)
            out0(i,3,j) = msfd(i,j)
          end do
        end do
        call mpi_gather(out0(1,1,1), ix*3*jxp,mpi_real8,                &
                      & out_0(1,1,1),ix*3*jxp,mpi_real8,                &
                      & 0,mpi_comm_world,ierr)
        if ( myid.eq.0 ) then
          do j = 1 , jxm1
            do i = 1 , ixm1
              veg2d_io(i,j) = out_0(i,1,j)
            end do
          end do
          do j = 1 , jx
            do i = 1 , ix
              msfx_io(i,j) = out_0(i,2,j)
              msfd_io(i,j) = out_0(i,3,j)
            end do
          end do
          call outtap0
          call gradsctl('OUT_HEAD.CTL')
        end if
      end if
 
!
!-----output for dataflow analyses:
!
      if ( iftape ) then
        if ( (jyear.eq.jyear0 .and. ktau.eq.0) .or.                     &
           & (mod(ntime,ntapfrq).eq.0 .and.                             &
           & (.not.(jyear.eq.jyearr.and.ktau.eq.ktaur))) ) then
!=======================================================================
!         gather  ua,va,ta,qva,qca,rainc,rainnc,tgb2d,swt2d,olcd2d,rno2d
          do j = 1 , jendl
            do k = 1 , kx
              do i = 1 , ix
                atm0(i,k,j) = ua(i,k,j)
                atm0(i,k+kx,j) = va(i,k,j)
                atm0(i,k+kx*2,j) = omega(i,k,j)
                atm0(i,k+kx*3,j) = ta(i,k,j)
                atm0(i,k+kx*4,j) = qva(i,k,j)
                atm0(i,k+kx*5,j) = qca(i,k,j)
              end do
            end do
            do i = 1 , ix
              atm0(i,1+kx*6,j) = psa(i,j)
              atm0(i,2+kx*6,j) = rainc(i,j)
              atm0(i,3+kx*6,j) = rainnc(i,j)
            end do
          end do
          do j = 1 , jendx
            do n = 1 , nnsg
              do i = 1 , ixm1
                atm0(i,3+kx*6+n,j) = ocld2d(n,i,j)
                atm0(i,3+kx*6+n+nnsg,j) = tgb2d(n,i,j)
                atm0(i,3+kx*6+n+nnsg*2,j) = swt2d(n,i,j)
                atm0(i,3+kx*6+n+nnsg*3,j) = rno2d(n,i,j)
              end do
            end do
          end do
          call mpi_gather(atm0(1,1,1), ix*(kx*6+3+nnsg*4)*jxp,mpi_real8,&
                        & atm_0(1,1,1),ix*(kx*6+3+nnsg*4)*jxp,mpi_real8,&
                        & 0,mpi_comm_world,ierr)
          if ( myid.eq.0 ) then
            do j = 1 , jx
              do k = 1 , kx
                do i = 1 , ix
                  ua_io(i,k,j) = atm_0(i,k,j)
                  va_io(i,k,j) = atm_0(i,k+kx,j)
                  omega_io(i,k,j) = atm_0(i,k+kx*2,j)
                  ta_io(i,k,j) = atm_0(i,k+kx*3,j)
                  qva_io(i,k,j) = atm_0(i,k+kx*4,j)
                  qca_io(i,k,j) = atm_0(i,k+kx*5,j)
                end do
              end do
              do i = 1 , ix
                psa_io(i,j) = atm_0(i,1+kx*6,j)
                rainc_io(i,j) = atm_0(i,2+kx*6,j)
                rainnc_io(i,j) = atm_0(i,3+kx*6,j)
              end do
            end do
            do j = 1 , jxm1
              do n = 1 , nnsg
                do i = 1 , ixm1
                  ocld2d_io(n,i,j) = atm_0(i,3+kx*6+n,j)
                  tgb2d_io(n,i,j) = atm_0(i,3+kx*6+n+nnsg,j)
                  swt2d_io(n,i,j) = atm_0(i,3+kx*6+n+nnsg*2,j)
                  rno2d_io(n,i,j) = atm_0(i,3+kx*6+n+nnsg*3,j)
                end do
              end do
            end do
            call outtap
          end if
          do j = 1 , jendx
            do i = 1 , ixm1
              do n = 1 , nnsg
                rno2d(n,i,j) = 0.
              end do
              rainc(i,j) = 0.
              rainnc(i,j) = 0.
            end do
          end do
        end if
      end if
 
!     Call surface output
 
      if ( ifbat ) then
        if ( (mod(ntime,kbats).eq.0 .and. (.not.(jyear.eq.jyearr.and.   &
           & ktau.eq.ktaur))) .or. (jyear.eq.jyear0 .and. ktau.eq.1) )  &
           & then

          call fillbat
          do j = 1 , jendx
            do l = 1 , numbat
              do i = 1 , ixm2
                bat0(i,l,j) = fbat(j,i,l)
              end do
            end do
          end do
          call mpi_gather(bat0(1,1,1), ixm2*numbat*jxp,mpi_real4,       &
                        & bat_0(1,1,1),ixm2*numbat*jxp,mpi_real4,       &
                        & 0,mpi_comm_world,ierr)
          if ( myid.eq.0 ) then
            do l = 1 , numbat
              do i = 1 , ixm2
                do j = 1 , jxm2
                  fbat_io(j,i,l) = bat_0(i,l,j+1)
                end do
              end do
            end do
            call outsrf
          end if

          do i = 1 , ixm2
            do j = 1 , jxp
              tgmx_o(j,i) = -1.E30
              t2mx_o(j,i) = -1.E30
              tgmn_o(j,i) = 1.E30
              t2mn_o(j,i) = 1.E30
              w10x_o(j,i) = -1.E30
              psmn_o(j,i) = 1.E30
            end do
          end do

          if ( ifsub .and. nsg.gt.1 ) then

            call fillsub
            do j = 1 , jxp
              do l = 1 , numsub
                do n = 1 , nnsg
                  do i = 1 , ixm2
                    sub0(i,n,l,j) = fsub(n,j,i,l)
                  end do
                end do
              end do
            end do
            call mpi_gather(sub0(1,1,1,1), ixm2*nnsg*numsub*jxp,        &
                          & mpi_real4,                                  &
                          & sub_0(1,1,1,1),ixm2*nnsg*numsub*jxp,        &
                          & mpi_real4,                                  &
                          & 0,mpi_comm_world,ierr)

            if ( myid.eq.0 ) then
              do l = 1 , numsub
                do j = 1 , jxm2
                  do n = 1 , nnsg
                    do i = 1 , ixm2
                      fsub_io(n,j,i,l) = sub_0(i,n,l,j+1)
                    end do
                  end do
                end do
              end do
            end if

            call outsub
          end if

        end if
      end if
 
!     Call radiation output
      if ( ifrad ) then
        if ( (mod(ntime,nradisp).eq.0 .and. (.not.(jyear.eq.jyearr.and. &
           & ktau.eq.ktaur))) .or. (jyear.eq.jyear0 .and. ktau.eq.1) )  &
           & then
!=======================================================================
!         frad2d, frad3d
          do n = 1 , nrad2d
            do j = 1 , jxp
              do i = 1 , ixm2
                rad0(i,n,j) = frad2d(j,i,n)
              end do
            end do
          end do
          do n = 1 , nrad3d
            do k = 1 , kx
              do j = 1 , jxp
                do i = 1 , ixm2
                  rad0(i,nrad2d+(n-1)*kx+k,j) = frad3d(j,i,k,n)
                end do
              end do
            end do
          end do
          call mpi_gather(rad0(1,1,1),ixm2*(nrad3d*kx+nrad2d)*jxp,      &
                        & mpi_real4,rad_0(1,1,1),(ixm2)                 &
                        & *(nrad3d*kx+nrad2d)*jxp,mpi_real4,0,          &
                        & mpi_comm_world,ierr)
          if ( myid.eq.0 ) then
            do n = 1 , nrad2d
              do j = 1 , jxm2
                do i = 1 , ixm2
                  frad2d_io(j,i,n) = rad_0(i,n,j+1)
                end do
              end do
            end do
            do n = 1 , nrad3d
              do k = 1 , kx
                do j = 1 , jxm2
                  do i = 1 , ixm2
                    frad3d_io(j,i,k,n) = rad_0(i,nrad2d+(n-1)*kx+k,j+1)
                  end do
                end do
              end do
            end do
            call radtap
          end if
        end if
      end if
 
!chem2
!     Call chem output
      if ( ifchem ) then
        if ( (jyear.eq.jyear0 .and. ktau.eq.1) .or.                     &
           & (mod(ntime,kchem).eq.0 .and.                               &
           & (.not.(jyear.eq.jyearr.and.ktau.eq.ktaur))) ) then
          do j = 1 , jendl
            do n = 1 , ntr
              do k = 1 , kx
                do i = 1 , ix
                  chem0(i,(n-1)*kx+k,j) = chia(i,k,j,n)
                end do
              end do
            end do
          end do
          do j = 1 , jendx
            do k = 1 , kx
              do i = 1 , ixm1
                chem0(i,ntr*kx+k,j) = aerext(i,k,j)
                chem0(i,ntr*kx+kx+k,j) = aerssa(i,k,j)
                chem0(i,ntr*kx+kx*2+k,j) = aerasp(i,k,j)
              end do
            end do
          end do
          do j = 1 , jendl
            do n = 1 , ntr
              do i = 1 , ix
                chem0(i,(ntr+3)*kx+n,j) = dtrace(i,j,n)
                chem0(i,(ntr+3)*kx+ntr+n,j) = wdlsc(i,j,n)
                chem0(i,(ntr+3)*kx+ntr*2+n,j) = wdcvc(i,j,n)
                chem0(i,(ntr+3)*kx+ntr*3+n,j) = ddsfc(i,j,n)
                chem0(i,(ntr+3)*kx+ntr*4+n,j) = wxsg(i,j,n)
                chem0(i,(ntr+3)*kx+ntr*5+n,j) = wxaq(i,j,n)
                chem0(i,(ntr+3)*kx+ntr*6+n,j) = cemtrac(i,j,n)
              end do
            end do
          end do
          do j = 1 , jendx
            do i = 1 , ixm1
              chem0(i,(ntr+3)*kx+ntr*7+1,j) = aertarf(i,j)
              chem0(i,(ntr+3)*kx+ntr*7+2,j) = aersrrf(i,j)
            end do
          end do
          do j = 1 , jendl
            do i = 1 , ix
              chem0(i,(ntr+3)*kx+ntr*7+3,j) = psa(i,j)
            end do
          end do
          call mpi_gather(chem0(1,1,1),ix*((ntr+3)*kx+ntr*7+3)*jxp,     &
                        & mpi_real8,chem_0(1,1,1),                      &
                        & ix*((ntr+3)*kx+ntr*7+3)*jxp,                  &
                        & mpi_real8,0,mpi_comm_world,ierr)
          if ( myid.eq.0 ) then
            do j = 1 , jx
              do n = 1 , ntr
                do k = 1 , kx
                  do i = 1 , ix
                    chia_io(i,k,j,n) = chem_0(i,(n-1)*kx+k,j)
                  end do
                end do
              end do
            end do
            do j = 1 , jxm1
              do k = 1 , kx
                do i = 1 , ixm1
                  aerext_io(i,k,j) = chem_0(i,ntr*kx+k,j+1)
                  aerssa_io(i,k,j) = chem_0(i,ntr*kx+kx+k,j+1)
                  aerasp_io(i,k,j) = chem_0(i,ntr*kx+kx*2+k,j+1)
                end do
              end do
            end do
            do j = 1 , jx
              do n = 1 , ntr
                do i = 1 , ix
                  dtrace_io(i,j,n) = chem_0(i,(ntr+3)*kx+n,j)
                  wdlsc_io(i,j,n) = chem_0(i,(ntr+3)*kx+ntr+n,j)
                  wdcvc_io(i,j,n) = chem_0(i,(ntr+3)*kx+ntr*2+n,j)
                  ddsfc_io(i,j,n) = chem_0(i,(ntr+3)*kx+ntr*3+n,j)
                  wxsg_io(i,j,n) = chem_0(i,(ntr+3)*kx+ntr*4+n,j)
                  wxaq_io(i,j,n) = chem_0(i,(ntr+3)*kx+ntr*5+n,j)
                  cemtrac_io(i,j,n) = chem_0(i,(ntr+3)*kx+ntr*6+n,j)
                end do
              end do
            end do
            do j = 1 , jxm1
              do i = 1 , ixm1
                aertarf_io(i,j) = chem_0(i,(ntr+3)*kx+ntr*7+1,j+1)
                aersrrf_io(i,j) = chem_0(i,(ntr+3)*kx+ntr*7+2,j+1)
              end do
            end do
            do j = 1 , jx
              do i = 1 , ix
                psa_io(i,j) = chem_0(i,(ntr+3)*kx+ntr*7+3,j)
              end do
            end do
            call chemtap
          end if
          do n = 1 , ntr
            do j = 1 , jendl
              do k = 1 , kx
                do i = 1 , ix
                  remlsc(i,k,j,n) = 0.
                  remcvc(i,k,j,n) = 0.
                  rxsg(i,k,j,n) = 0.
                  rxsaq1(i,k,j,n) = 0.
                  rxsaq2(i,k,j,n) = 0.
                end do
              end do
            end do
          end do
          do n = 1 , ntr
            do j = 1 , jendl
              do i = 1 , ix
                cemtr(i,j,n) = 0.
                remdrd(i,j,n) = 0.
                wdlsc(i,j,n) = 0.
                wdcvc(i,j,n) = 0.
                ddsfc(i,j,n) = 0.
                wxsg(i,j,n) = 0.
                wxaq(i,j,n) = 0.
                cemtrac(i,j,n) = 0.
              end do
            end do
          end do
          do j = 1 , jendl
            do i = 1 , ixm1
              aertarf(i,j) = 0.
              aersrrf(i,j) = 0.
            end do
          end do
        end if
      end if
!chem2
!
!-----output for restart:
!
      if ( ifsave ) then
        if ( ((lday.eq.1 .and. lhour.eq.0 .and. dabs(xtime).lt.0.00001) &
           & .and. ldatez.ne.idate1) .or. nnnnnn.eq.nnnend ) then
          do j = 1 , jendl
            do k = 1 , kx
              do i = 1 , ix
                sav0(i,k,j) = ub0(i,k,j)
                sav0(i,kx+k,j) = vb0(i,k,j)
                sav0(i,kx*2+k,j) = qb0(i,k,j)
                sav0(i,kx*3+k,j) = tb0(i,k,j)
              end do
            end do
            do i = 1 , ix
              sav0(i,kx*4+1,j) = ps0(i,j)
              sav0(i,kx*4+2,j) = ts0(i,j)
            end do
          end do
          allrec = kx*4 + 2
          call mpi_gather(sav0(1,1,1), ix*allrec*jxp,mpi_real8,         &
                        & sav_0(1,1,1),ix*allrec*jxp,mpi_real8,         &
                        & 0,mpi_comm_world,ierr)
          if ( myid.eq.0 ) then
            do j = 1 , jx
              do k = 1 , kx
                do i = 1 , ix
                  ub0_io(i,k,j) = sav_0(i,k,j)
                  vb0_io(i,k,j) = sav_0(i,kx+k,j)
                  qb0_io(i,k,j) = sav_0(i,kx*2+k,j)
                  tb0_io(i,k,j) = sav_0(i,kx*3+k,j)
                end do
              end do
              do i = 1 , ix
                ps0_io(i,j) = sav_0(i,kx*4+1,j)
                ts0_io(i,j) = sav_0(i,kx*4+2,j)
              end do
            end do
          end if
          if ( ehso4 ) then
            do j = 1 , jendl
              do k = 1 , kx
                do i = 1 , ix
                  sav0s(i,k,j) = so0(i,k,j)
                end do
              end do
            end do
            call mpi_gather(sav0s(1,1,1), ix*kx*jxp,mpi_real8,          &
                          & sav_0s(1,1,1),ix*kx*jxp,mpi_real8,          &
                          & 0,mpi_comm_world,ierr)
            if ( myid.eq.0 ) then
              do j = 1 , jx
                do k = 1 , kx
                  do i = 1 , ix
                    so0_io(i,k,j) = sav_0s(i,k,j)
                  end do
                end do
              end do
            end if
          end if
          do j = 1 , jendl
            do k = 1 , kx
              do i = 1 , ix
                sav0(i,k,j) = ua(i,k,j)
                sav0(i,kx+k,j) = ub(i,k,j)
                sav0(i,kx*2+k,j) = va(i,k,j)
                sav0(i,kx*3+k,j) = vb(i,k,j)
              end do
            end do
            do i = 1 , ix
              sav0(i,kx*4+1,j) = psa(i,j)
              sav0(i,kx*4+2,j) = psb(i,j)
            end do
          end do
          allrec = kx*4 + 2
          call mpi_gather(sav0(1,1,1), ix*allrec*jxp,mpi_real8,         &
                        & sav_0(1,1,1),ix*allrec*jxp,mpi_real8,         &
                        & 0,mpi_comm_world,ierr)
          if ( myid.eq.0 ) then
            do j = 1 , jx
              do k = 1 , kx
                do i = 1 , ix
                  ua_io(i,k,j) = sav_0(i,k,j)
                  ub_io(i,k,j) = sav_0(i,kx+k,j)
                  va_io(i,k,j) = sav_0(i,kx*2+k,j)
                  vb_io(i,k,j) = sav_0(i,kx*3+k,j)
                end do
              end do
              do i = 1 , ix
                psa_io(i,j) = sav_0(i,kx*4+1,j)
                psb_io(i,j) = sav_0(i,kx*4+2,j)
              end do
            end do
          end if
          do j = 1 , jendl
            do k = 1 , kx
              do i = 1 , ix
                sav0(i,k,j) = ta(i,k,j)
                sav0(i,kx+k,j) = tb(i,k,j)
                sav0(i,kx*2+k,j) = qva(i,k,j)
                sav0(i,kx*3+k,j) = qvb(i,k,j)
              end do
            end do
            do i = 1 , ix
              sav0(i,kx*4+1,j) = tga(i,j)
              sav0(i,kx*4+2,j) = tgb(i,j)
            end do
          end do
          allrec = kx*4 + 2
          call mpi_gather(sav0(1,1,1), ix*allrec*jxp,mpi_real8,         &
                        & sav_0(1,1,1),ix*allrec*jxp,mpi_real8,         &
                        & 0,mpi_comm_world,ierr)
          if ( myid.eq.0 ) then
            do j = 1 , jx
              do k = 1 , kx
                do i = 1 , ix
                  ta_io(i,k,j) = sav_0(i,k,j)
                  tb_io(i,k,j) = sav_0(i,kx+k,j)
                  qva_io(i,k,j) = sav_0(i,kx*2+k,j)
                  qvb_io(i,k,j) = sav_0(i,kx*3+k,j)
                end do
              end do
              do i = 1 , ix
                tga_io(i,j) = sav_0(i,kx*4+1,j)
                tgb_io(i,j) = sav_0(i,kx*4+2,j)
              end do
            end do
          end if
          do j = 1 , jendl
            do k = 1 , kx
              do i = 1 , ix
                sav0(i,k,j) = qca(i,k,j)
                sav0(i,kx+k,j) = qcb(i,k,j)
                sav0(i,kx*2+k,j) = fcc(i,k,j)
              end do
            end do
            do i = 1 , ix
              sav0(i,kx*4+1,j) = rainc(i,j)
              sav0(i,kx*4+2,j) = rainnc(i,j)
            end do
          end do
          do j = 1 , jendx
            do k = 1 , kx
              do i = 1 , ixm1
                sav0(i,kx*3+k,j) = heatrt(i,k,j)
              end do
            end do
          end do
          allrec = kx*4 + 2
          call mpi_gather(sav0(1,1,1), ix*allrec*jxp,mpi_real8,         &
                        & sav_0(1,1,1),ix*allrec*jxp,mpi_real8,         &
                        & 0,mpi_comm_world,ierr)
          if ( myid.eq.0 ) then
            do j = 1 , jx
              do k = 1 , kx
                do i = 1 , ix
                  qca_io(i,k,j) = sav_0(i,k,j)
                  qcb_io(i,k,j) = sav_0(i,kx+k,j)
                  fcc_io(i,k,j) = sav_0(i,kx*2+k,j)
                end do
              end do
              do i = 1 , ix
                rainc_io(i,j) = sav_0(i,kx*4+1,j)
                rainnc_io(i,j) = sav_0(i,kx*4+2,j)
              end do
            end do
            do j = 1 , jxm1
              do k = 1 , kx
                do i = 1 , ixm1
                  heatrt_io(i,k,j) = sav_0(i,kx*3+k,j)
                end do
              end do
            end do
          end if
          do j = 1 , jendl
            do i = 1 , ix
              sav0a(i,1,j) = hfx(i,j)
              sav0a(i,2,j) = qfx(i,j)
              sav0a(i,3,j) = uvdrag(i,j)
              sav0a(i,4,j) = tgbb(i,j)
            end do
            do n = 1 , nnsg
              do i = 1 , ix
                sav0a(i,4+n,j) = snowc(n,i,j)
              end do
            end do
          end do
          do j = 1 , jendx
            do k = 1 , kxp1
              do i = 1 , ixm1
                sav0a(i,nnsg+4+k,j) = o3prof(i,k,j)
              end do
            end do
          end do
          allrec = 5 + nnsg + kx
          call mpi_gather(sav0a(1,1,1), ix*allrec*jxp,mpi_real8,        &
                        & sav_0a(1,1,1),ix*allrec*jxp,mpi_real8,        &
                        & 0,mpi_comm_world,ierr)
          if ( myid.eq.0 ) then
            do j = 1 , jx
              do i = 1 , ix
                hfx_io(i,j) = sav_0a(i,1,j)
                qfx_io(i,j) = sav_0a(i,2,j)
                uvdrag_io(i,j) = sav_0a(i,3,j)
                tgbb_io(i,j) = sav_0a(i,4,j)
              end do
              do n = 1 , nnsg
                do i = 1 , ix
                  snowc_io(n,i,j) = sav_0a(i,4+n,j)
                end do
              end do
            end do
            do j = 1 , jxm1
              do k = 1 , kxp1
                do i = 1 , ixm1
                  o3prof_io(i,k,j) = sav_0a(i,4+nnsg+k,j)
                end do
              end do
            end do
          end if
          if ( iocnflx.eq.2 )                                           &
             & call mpi_gather(zpbl(1,1),   ix*jxp,mpi_real8,           &
             &                 zpbl_io(1,1),ix*jxp,mpi_real8,           &
             &                 0,mpi_comm_world,ierr)
          if ( icup.eq.1 ) then
            do j = 1 , jendl
              do k = 1 , kx
                do i = 1 , ix
                  sav0c(i,k,j) = rsheat(i,k,j)
                  sav0c(i,kx+k,j) = rswat(i,k,j)
                end do
              end do
            end do
            allrec = kx*2
            call mpi_gather(sav0c(1,1,1), ix*allrec*jxp,mpi_real8,      &
                          & sav_0c(1,1,1),ix*allrec*jxp,mpi_real8,      &
                          & 0,mpi_comm_world,ierr)
            if ( myid.eq.0 ) then
              do j = 1 , jx
                do k = 1 , kx
                  do i = 1 , ix
                    rsheat_io(i,k,j) = sav_0c(i,k,j)
                    rswat_io(i,k,j) = sav_0c(i,kx+k,j)
                  end do
                end do
              end do
            end if
          else if ( icup.eq.3 ) then
            do j = 1 , jendl
              do k = 1 , kx
                do i = 1 , ix
                  sav0b(i,k,j) = tbase(i,k,j)
                end do
              end do
              do i = 1 , ix
                sav0b(i,kxp1,j) = cldefi(i,j)
              end do
            end do
            allrec = kxp1
            call mpi_gather(sav0b(1,1,1), ix*allrec*jxp,mpi_real8,      &
                          & sav_0b(1,1,1),ix*allrec*jxp,mpi_real8,      &
                          & 0,mpi_comm_world,ierr)
            if ( myid.eq.0 ) then
              do j = 1 , jx
                do k = 1 , kx
                  do i = 1 , ix
                    tbase_io(i,k,j) = sav_0b(i,k,j)
                  end do
                end do
                do i = 1 , ix
                  cldefi_io(i,j) = sav_0b(i,kxp1,j)
                end do
              end do
            end if
          else if ( icup.eq.4 ) then
            call mpi_gather(cbmf2d(1,1),   ix*jxp,mpi_real8,            &
                          & cbmf2d_io(1,1),ix*jxp,mpi_real8,            &
                          & 0,mpi_comm_world,ierr)
          else
          end if
          do j = 1 , jendx
            do l = 1 , 4
              do k = 1 , kx
                do i = 1 , ixm1
                  sav1(i,(l-1)*kx+k,j) = absnxt(i,k,l,j)
                end do
              end do
            end do
          end do
          allrec = kx*4
          do j = 1 , jendx
            do l = 1 , kxp1
              do k = 1 , kxp1
                do i = 1 , ixm1
                  sav1(i,allrec+(l-1)*(kxp1)+k,j) = abstot(i,k,l,j)
                end do
              end do
            end do
          end do
          allrec = allrec + (kxp1)*(kx+1)
          do j = 1 , jendx
            do k = 1 , kxp1
              do i = 1 , ixm1
                sav1(i,allrec+k,j) = emstot(i,k,j)
              end do
            end do
          end do
          allrec = kx*4+(kxp1*kxp2)
          call mpi_gather(sav1(1,1,1), ixm1*allrec*jxp,mpi_real8,       &
                        & sav_1(1,1,1),ixm1*allrec*jxp,mpi_real8,       &
                        & 0,mpi_comm_world,ierr)
          if ( myid.eq.0 ) then
            do j = 1 , jxm1
              do l = 1 , 4
                do k = 1 , kx
                  do i = 1 , ixm1
                    absnxt_io(i,k,l,j) = sav_1(i,(l-1)*kx+k,j)
                  end do
                end do
              end do
            end do
            allrec = kx*4
            do j = 1 , jxm1
              do l = 1 , kxp1
                do k = 1 , kxp1
                  do i = 1 , ixm1
                    abstot_io(i,k,l,j)                                  &
                    & = sav_1(i,allrec+(l-1)*(kxp1)+k,j)
                  end do
                end do
              end do
            end do
            allrec = allrec + (kxp1)*(kx+1)
            do j = 1 , jxm1
              do k = 1 , kxp1
                do i = 1 , ixm1
                  emstot_io(i,k,j) = sav_1(i,allrec+k,j)
                end do
              end do
            end do
          end if
          do j = 1 , jendx
            do n = 1 , nnsg
              do i = 1 , ixm1
                sav2(i,n,j) = taf2d(n,i,j)
                sav2(i,nnsg+n,j) = tlef2d(n,i,j)
                sav2(i,nnsg*2+n,j) = ssw2d(n,i,j)
                sav2(i,nnsg*3+n,j) = srw2d(n,i,j)
              end do
            end do
            do i = 1 , ixm1
              sav2(i,nnsg*4+1,j) = sol2d(i,j)
              sav2(i,nnsg*4+2,j) = solvd2d(i,j)
              sav2(i,nnsg*4+3,j) = solvs2d(i,j)
              sav2(i,nnsg*4+4,j) = flw2d(i,j)
            end do
          end do
          allrec = nnsg*4 + 4
          call mpi_gather(sav2(1,1,1), ixm1*allrec*jxp,mpi_real8,       &
                        & sav_2(1,1,1),ixm1*allrec*jxp,mpi_real8,       &
                        & 0,mpi_comm_world,ierr)
          if ( myid.eq.0 ) then
            do j = 1 , jxm1
              do n = 1 , nnsg
                do i = 1 , ixm1
                  taf2d_io(n,i,j) = sav_2(i,n,j)
                  tlef2d_io(n,i,j) = sav_2(i,nnsg+n,j)
                  ssw2d_io(n,i,j) = sav_2(i,nnsg*2+n,j)
                  srw2d_io(n,i,j) = sav_2(i,nnsg*3+n,j)
                end do
              end do
              do i = 1 , ixm1
                sol2d_io(i,j) = sav_2(i,nnsg*4+1,j)
                solvd2d_io(i,j) = sav_2(i,nnsg*4+2,j)
                solvs2d_io(i,j) = sav_2(i,nnsg*4+3,j)
                flw2d_io(i,j) = sav_2(i,nnsg*4+4,j)
              end do
            end do
          end if
          do j = 1 , jendx
            do n = 1 , nnsg
              do i = 1 , ixm1
                sav2(i,n,j) = tgb2d(n,i,j)
                sav2(i,nnsg+n,j) = swt2d(n,i,j)
                sav2(i,nnsg*2+n,j) = scv2d(n,i,j)
                sav2(i,nnsg*3+n,j) = gwet2d(n,i,j)
              end do
            end do
            do i = 1 , ixm1
              sav2(i,nnsg*4+1,j) = flwd2d(i,j)
              sav2(i,nnsg*4+2,j) = fsw2d(i,j)
              sav2(i,nnsg*4+3,j) = sabv2d(i,j)
              sav2(i,nnsg*4+4,j) = sinc2d(i,j)
            end do
          end do
          allrec = nnsg*4 + 4
          call mpi_gather(sav2(1,1,1), ixm1*allrec*jxp,mpi_real8,       &
                        & sav_2(1,1,1),ixm1*allrec*jxp,mpi_real8,       &
                        & 0,mpi_comm_world,ierr)
          if ( myid.eq.0 ) then
            do j = 1 , jxm1
              do n = 1 , nnsg
                do i = 1 , ixm1
                  tgb2d_io(n,i,j) = sav_2(i,n,j)
                  swt2d_io(n,i,j) = sav_2(i,nnsg+n,j)
                  scv2d_io(n,i,j) = sav_2(i,nnsg*2+n,j)
                  gwet2d_io(n,i,j) = sav_2(i,nnsg*3+n,j)
                end do
              end do
              do i = 1 , ixm1
                flwd2d_io(i,j) = sav_2(i,nnsg*4+1,j)
                fsw2d_io(i,j) = sav_2(i,nnsg*4+2,j)
                sabv2d_io(i,j) = sav_2(i,nnsg*4+3,j)
                sinc2d_io(i,j) = sav_2(i,nnsg*4+4,j)
              end do
            end do
          end if
          do j = 1 , jendx
            do n = 1 , nnsg
              do i = 1 , ixm1
                sav2(i,n,j) = veg2d1(n,i,j)
                sav2(i,nnsg+n,j) = sag2d(n,i,j)
                sav2(i,nnsg*2+n,j) = sice2d(n,i,j)
                sav2(i,nnsg*3+n,j) = dew2d(n,i,j)
              end do
            end do
            do i = 1 , ixm1
              sav2(i,nnsg*4+1,j) = pptnc(i,j)
              sav2(i,nnsg*4+2,j) = pptc(i,j)
              sav2(i,nnsg*4+3,j) = prca2d(i,j)
              sav2(i,nnsg*4+4,j) = prnca2d(i,j)
            end do
          end do
          allrec = nnsg*4 + 4
          call mpi_gather(sav2(1,1,1), ixm1*allrec*jxp,mpi_real8,       &
                        & sav_2(1,1,1),ixm1*allrec*jxp,mpi_real8,       &
                        & 0,mpi_comm_world,ierr)
          if ( myid.eq.0 ) then
            do j = 1 , jxm1
              do n = 1 , nnsg
                do i = 1 , ixm1
                  veg2d1_io(n,i,j) = sav_2(i,n,j)
                  sag2d_io(n,i,j) = sav_2(i,nnsg+n,j)
                  sice2d_io(n,i,j) = sav_2(i,nnsg*2+n,j)
                  dew2d_io(n,i,j) = sav_2(i,nnsg*3+n,j)
                end do
              end do
              do i = 1 , ixm1
                pptnc_io(i,j) = sav_2(i,nnsg*4+1,j)
                pptc_io(i,j) = sav_2(i,nnsg*4+2,j)
                prca2d_io(i,j) = sav_2(i,nnsg*4+3,j)
                prnca2d_io(i,j) = sav_2(i,nnsg*4+4,j)
              end do
            end do
          end if
          do j = 1 , jendx
            do n = 1 , nnsg
              do i = 1 , ixm1
                sav2a(i,n,j) = ircp2d(n,i,j)
                sav2a(i,nnsg+n,j) = text2d(n,i,j)
                sav2a(i,nnsg*2+n,j) = col2d(n,i,j)
                sav2a(i,nnsg*3+n,j) = ocld2d(n,i,j)
                sav2a(i,nnsg*4+n,j) = tg2d(n,i,j)
              end do
            end do
            do i = 1 , ixm1
              sav2a(i,nnsg*5+1,j) = veg2d(i,j)
            end do
          end do
          allrec = nnsg*5 + 1
          call mpi_gather(sav2a(1,1,1), ixm1*allrec*jxp,mpi_real8,      &
                        & sav_2a(1,1,1),ixm1*allrec*jxp,mpi_real8,      &
                        & 0,mpi_comm_world,ierr)
          if ( myid.eq.0 ) then
            do j = 1 , jxm1
              do n = 1 , nnsg
                do i = 1 , ixm1
                  ircp2d_io(n,i,j) = sav_2a(i,n,j)
                  text2d_io(n,i,j) = sav_2a(i,nnsg+n,j)
                  col2d_io(n,i,j) = sav_2a(i,nnsg*2+n,j)
                  ocld2d_io(n,i,j) = sav_2a(i,nnsg*3+n,j)
                  tg2d_io(n,i,j) = sav_2a(i,nnsg*4+n,j)
                end do
              end do
              do i = 1 , ixm1
                veg2d_io(i,j) = sav_2a(i,nnsg*5+1,j)
              end do
            end do
          end if
 
          if ( ichem.eq.1 ) then
            do j = 1 , jendl
              do n = 1 , ntr
                do k = 1 , kx
                  do i = 1 , ix
                    sav4(i,(n-1)*kx+k,j) = chia(i,k,j,n)
                    sav4(i,ntr*kx+(n-1)*kx+k,j) = chib(i,k,j,n)
                    sav4(i,ntr*kx*2+(n-1)*kx+k,j) = remlsc(i,k,j,n)
                    sav4(i,ntr*kx*3+(n-1)*kx+k,j) = remcvc(i,k,j,n)
                  end do
                end do
              end do
            end do
            allrec = 4*ntr*kx
            do j = 1 , jendl
              do n = 1 , ntr
                do i = 1 , ix
                  sav4(i,allrec+n,j) = remdrd(i,j,n)
                end do
              end do
            end do
            allrec = ntr*(kx*4+1)
            call mpi_gather(sav4(1,1,1), ix*allrec*jxp,mpi_real8,       &
                          & sav_4(1,1,1),ix*allrec*jxp,mpi_real8,       &
                          & 0,mpi_comm_world,ierr)
            if ( myid.eq.0 ) then
              do j = 1 , jx
                do n = 1 , ntr
                  do k = 1 , kx
                    do i = 1 , ix
                      chia_io(i,k,j,n) = sav_4(i,(n-1)*kx+k,j)
                      chib_io(i,k,j,n) = sav_4(i,ntr*kx+(n-1)*kx+k,j)
                      remlsc_io(i,k,j,n)                                &
                      & = sav_4(i,ntr*kx*2+(n-1)*kx+k,j)
                      remcvc_io(i,k,j,n)                                &
                      & = sav_4(i,ntr*kx*3+(n-1)*kx+k,j)
                    end do
                  end do
                end do
              end do
              allrec = 4*ntr*kx
              do j = 1 , jx
                do n = 1 , ntr
                  do i = 1 , ix
                    remdrd_io(i,j,n) = sav_4(i,allrec+n,j)
                  end do
                end do
              end do
            end if
            do j = 1 , jendx
              do i = 1 , ixm1
                sav4a(i,1,j) = ssw2da(i,j)
                sav4a(i,2,j) = sdeltk2d(i,j)
                sav4a(i,3,j) = sdelqk2d(i,j)
                sav4a(i,4,j) = sfracv2d(i,j)
                sav4a(i,5,j) = sfracb2d(i,j)
                sav4a(i,6,j) = sfracs2d(i,j)
                sav4a(i,7,j) = svegfrac2d(i,j)
              end do
            end do
            call mpi_gather(sav4a, ixm1*7*jxp,mpi_real8,                &
                          & sav_4a,ixm1*7*jxp,mpi_real8,                &
                          & 0,mpi_comm_world,ierr)
            if ( myid.eq.0 ) then
              do j = 1 , jxm1
                do i = 1 , ixm1
                  ssw2da_io(i,j) = sav_4a(i,1,j)
                  sdeltk2d_io(i,j) = sav_4a(i,2,j)
                  sdelqk2d_io(i,j) = sav_4a(i,3,j)
                  sfracv2d_io(i,j) = sav_4a(i,4,j)
                  sfracb2d_io(i,j) = sav_4a(i,5,j)
                  sfracs2d_io(i,j) = sav_4a(i,6,j)
                  svegfrac2d_io(i,j) = sav_4a(i,7,j)
                end do
              end do
            end if
          end if
          do j = 1 , jendl
            do n = 1 , nsplit
              do i = 1 , ix
                sav0d(i,n,j) = dstor(i,j,n)
                sav0d(i,n+nsplit,j) = hstor(i,j,n)
              end do
            end do
          end do
          call mpi_gather(sav0d(1,1,1), ix*nsplit*2*jxp,mpi_real8,      &
                        & sav_0d(1,1,1),ix*nsplit*2*jxp,mpi_real8,      &
                        & 0,mpi_comm_world,ierr)
          if ( myid.eq.0 ) then
            do j = 1 , jx
              do n = 1 , nsplit
                do i = 1 , ix
                  dstor_io(i,j,n) = sav_0d(i,n,j)
                  hstor_io(i,j,n) = sav_0d(i,n+nsplit,j)
                end do
              end do
            end do
          end if
          do j = 1 , jendl
            do k = 1 , kx
              sav6(k,1,j) = ui1(k,j)
              sav6(k,2,j) = ui2(k,j)
              sav6(k,3,j) = uilx(k,j)
              sav6(k,4,j) = uil(k,j)
              sav6(k,5,j) = vi1(k,j)
              sav6(k,6,j) = vi2(k,j)
              sav6(k,7,j) = vilx(k,j)
              sav6(k,8,j) = vil(k,j)
            end do
          end do
          call mpi_gather(sav6(1,1,1), kx*8*jxp,mpi_real8,              &
                        & sav_6(1,1,1),kx*8*jxp,mpi_real8,              &
                        & 0,mpi_comm_world,ierr)
          if ( myid.eq.0 ) then
            do j = 1 , jx
              do k = 1 , kx
                ui1_io(k,j) = sav_6(k,1,j)
                ui2_io(k,j) = sav_6(k,2,j)
                uilx_io(k,j) = sav_6(k,3,j)
                uil_io(k,j) = sav_6(k,4,j)
                vi1_io(k,j) = sav_6(k,5,j)
                vi2_io(k,j) = sav_6(k,6,j)
                vilx_io(k,j) = sav_6(k,7,j)
                vil_io(k,j) = sav_6(k,8,j)
              end do
            end do
          end if
          call mpi_bcast(ujlx(1,1),ix*kx,mpi_real8,nproc-1,             &
                       & mpi_comm_world,ierr)
          call mpi_bcast(ujl(1,1),ix*kx,mpi_real8,nproc-1,              &
                       & mpi_comm_world,ierr)
          call mpi_bcast(vjlx(1,1),ix*kx,mpi_real8,nproc-1,             &
                       & mpi_comm_world,ierr)
          call mpi_bcast(vjl(1,1),ix*kx,mpi_real8,nproc-1,              &
                       & mpi_comm_world,ierr)
          if ( myid.eq.0 ) then
            close (iutsav)
            itype = 'SAV'
            write (newfil,99001) itype , idatex
            open (iutsav,file='output/'//newfil,status='replace',       &
                 &form='unformatted')
            print * , 'OPENING NEW SAV FILE: ' , newfil
            call outsav(iutsav)
            print * , 'restart written date = ' , ldatez + xtime/1440.
            close (iutsav)
          end if
        else if ( mod(ntime,nsavfrq).eq.0 .and.                         &
                & (.not.(jyear.eq.jyearr .and. ktau.eq.ktaur)) ) then
          if ( myid.eq.0 ) then
            close (iutsav)
            itype = 'SAV'
            write (tmpfil,99002) itype , idatex
            open (iutsav,file=tmpfil,status='unknown',                  &
                & form='unformatted')
            call outsav(iutsav)
            close (iutsav)
            print * , 'SAVTMP RESTART WRITTEN: idatex=' , idatex ,      &
                 &'ktau=' , ktau
            if ( oldsav(1:3).eq.'SAV' ) then
              inquire (file=oldsav,exist=there)
              if ( there ) then
                call unlink(oldsav, iresult)
              end if
            end if
            oldsav = tmpfil
          end if
        else
        end if
      end if
!
!-----printer output:
!
      if ( myid.eq.0 ) then
        if ( ifprt ) then
          if ( (jyear.eq.jyear0 .and. ktau.eq.0) .or. mod(ntime,nprtfrq)&
             & .eq.0 ) then
            write (*,*) 'outprt is not available for parallel code'
            idum = iexec
          end if
        end if
 
        if ( lday.eq.1 .and. lhour.eq.0 .and. nint(xtime).eq.0 .and.    &
           & (.not.(jyear.eq.jyearr .and. ktau.eq.ktaur)) .and.         &
           & nnnnnn.ne.nnnend ) call mkfile
      end if

#else

      if ( jyear.eq.jyearr .and. ktau.eq.ktaur ) then
        if ( iotyp.eq.1 ) then
          print * , 'Writing output files in direct access format'
        else if ( iotyp.eq.2 ) then
          print * , 'Writing output files in sequential format'
        else
          write (aline,*) 'iotyp = ' , iotyp
          call say
          call fatal(__FILE__,__LINE__,'Output format does not exist')
        end if
        nrcout = 0
        call mkfile
      end if
 
      if ( jyear.eq.jyear0 .and. ktau.eq.0 ) then
        call outtap0
        call gradsctl('OUT_HEAD.CTL')
      end if
 
!
!-----output for dataflow analyses:
!
      if ( iftape ) then
        if ( (jyear.eq.jyear0 .and. ktau.eq.0) .or.                     &
           & (mod(ntime,ntapfrq).eq.0 .and.                             &
           & (.not.(jyear.eq.jyearr.and.ktau.eq.ktaur))) ) call outtap
      end if
 
 
!     Call surface output
      if ( ifbat ) then
        if ( (mod(ntime,kbats).eq.0 .and. (.not.(jyear.eq.jyearr.and.   &
           & ktau.eq.ktaur))) .or. (jyear.eq.jyear0 .and. ktau.eq.1) )  &
           & then
          call fillbat
          call outsrf
          do i = 1 , ixm2
            do j = 1 , jxm2
              tgmx_o(j,i) = -1.E30
              t2mx_o(j,i) = -1.E30
              tgmn_o(j,i) = 1.E30
              t2mn_o(j,i) = 1.E30
              w10x_o(j,i) = -1.E30
              psmn_o(j,i) = 1.E30
            end do
          end do
          if ( ifsub .and. nsg.gt.1 ) then
            call fillsub
            call outsub
          end if
        end if
      end if
 
!     Call radiation output
      if ( ifrad ) then
        if ( (mod(ntime,nradisp).eq.0 .and. (.not.(jyear.eq.jyearr.and. &
           & ktau.eq.ktaur))) .or. (jyear.eq.jyear0 .and. ktau.eq.1) )  &
           & call radtap
      end if
 
!chem2
!     Call chem output
      if ( ifchem ) then
        if ( (jyear.eq.jyear0 .and. ktau.eq.1) .or.                     &
           & (mod(ntime,kchem).eq.0 .and.                               &
           & (.not.(jyear.eq.jyearr.and.ktau.eq.ktaur))) ) call chemtap
      end if
!chem2
!
!-----output for restart:
!
      if ( ifsave ) then
        if ( ((lday.eq.1 .and. lhour.eq.0 .and. dabs(xtime).lt.0.00001) &
           & .and. ldatez.ne.idate1) .or. nnnnnn.eq.nnnend ) then
          close (iutsav)
          itype = 'SAV'
          write (newfil,99001) itype , idatex
          open (iutsav,file='output/'//newfil,status='replace',         &
               &form='unformatted')
          print * , 'OPENING NEW SAV FILE: ' , newfil
          call outsav(iutsav)
          print * , 'restart written date = ' , ldatez + xtime/1440.
          close (iutsav)
        else if ( mod(ntime,nsavfrq).eq.0 .and.                         &
                & (.not.(jyear.eq.jyearr .and. ktau.eq.ktaur)) ) then
          close (iutsav)
          itype = 'SAV'
          write (tmpfil,99002) itype , idatex
          open (iutsav,file=tmpfil,status='unknown',form='unformatted')
          call outsav(iutsav)
          close (iutsav)
          print * , 'SAVTMP RESTART WRITTEN: idatex=' , idatex ,        &
              & 'ktau=' , ktau
          if ( oldsav(1:3).eq.'SAV' ) then
            inquire (file=oldsav,exist=there)
            if ( there ) then
              call unlink(oldsav, iresult)
            end if
          end if
          oldsav = tmpfil
        else
        end if
      end if
!
!-----printer output:
!
      if ( ifprt ) then
        if ( (jyear.eq.jyear0 .and. ktau.eq.0) .or. mod(ntime,nprtfrq)  &
           & .eq.0 ) call outprt(iexec)
      end if
 
      if ( lday.eq.1 .and. lhour.eq.0 .and. nint(xtime).eq.0 .and.      &
         & (.not.(jyear.eq.jyearr .and. ktau.eq.ktaur)) .and.           &
         & nnnnnn.ne.nnnend ) call mkfile
#endif

99001 format (a3,'.',i10)
99002 format (a3,'TMP.',i10)
 
      end subroutine output
