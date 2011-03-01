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
 
      module mod_output

      use mod_dynparam
      use mod_runparams
      use mod_ncio
      use mod_date
      use mod_lake
      use mod_main
      use mod_mainchem
      use mod_bats
      use mod_message
      use mod_bdycod
      use mod_pmoist
      use mod_rad
      use mod_trachem
      use mod_cvaria
      use mod_outrad
      use mod_radiation
      use mod_split
      use mod_savefile
      use mod_service
      use mod_cu_bm
#ifdef MPP1
      use mod_mppio
#ifdef CLM
      use mod_clm
#endif
#else
      use mod_lake
#endif

      private

      integer :: iolak
      logical :: lskipsrf , lskiprad , lskipche

      public :: output , mkfile

      data iolak /0/
      data lskipsrf /.false./
      data lskiprad /.false./
      data lskipche /.false./

      contains

      subroutine output

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!     this subroutine handles all of the output                       c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#ifdef MPP1
      use mpi
#endif
      implicit none
!
      integer :: i , j
#ifdef MPP1
      integer :: allrec , ierr , l , k , n
#endif
      logical :: ldoatm , ldosrf , ldorad , ldoche , ldosav , ldotmp
      logical :: lstartup
      character (len=50) :: subroutine_name='output'
      integer :: idindx=0
!
!
      call time_begin(subroutine_name,idindx)
!
!----------------------------------------------------------------------
!
      if ( jyear /= jyear0 .or. ktau /= 0 ) then
        if ( mod(idnint(xtime),60) < mod(idnint(xtime-dtmin),60) )   &
           & idatex = idatex + 1
        if ( dabs(xtime) < 0.00001D0 ) idatex = ldatez
      end if
 
      lstartup = .false.
#ifdef MPP1
      if ( myid == 0 ) then
#endif        
        if ( (jyear == jyear0 .and. ktau == 0) .or. &
             (ifrest .and. .not. done_restart) ) then
          call mkfile
          lstartup = .true.
        end if
#ifdef MPP1
      end if
#endif
!
      ldoatm = .false.
      ldosrf = .false.
      ldorad = .false.
      ldoche = .false.
      ldosav = .false.
      ldotmp = .false.

      if ( mod(ntime,nsavfrq) == 0 .and. ldatez /= idate1 ) then
        ldotmp = .true.
      end if
      if ( ((lday==1 .and. lhour==0 .and. dabs(xtime)<0.00001) .and. &
            ldatez /= idate1) .or. nnnnnn == nnnend ) then
        ldosav = .true.
        ldotmp = .false.
      end if
      if ( (jyear == jyear0 .and. ktau == 0) .or. &
           mod(ntime,ntapfrq) == 0) then
        ldoatm = .true.
      end if
      if ( (jyear == jyear0 .and. ktau == 0) .or. & 
           mod(ntime,kbats) == 0) then
        ldosrf = .true.
      end if
      if ( (jyear == jyear0 .and. ktau == 0) .or. &
           mod(ntime,nradisp) == 0) then
        ldorad = .true.
      end if
      if ( (jyear == jyear0 .and. ktau == 0) .or. &
           mod(ntime,kchem) == 0) then
        ldoche = .true.
      end if

      if ( ifrest .and. .not. done_restart ) then
        ldoatm = .false.
        ldosrf = .false.
        ldorad = .false.
        ldoche = .false.
      end if

      if ( lskipsrf ) then
        lskipsrf = .false.
        ldosrf = .true.
      end if
      if ( lskiprad ) then
        lskiprad = .false.
        ldorad = .true.
      end if
      if ( lskipche ) then
        lskipche = .false.
        ldoche = .true.
      end if
!
      if ( jyear == jyear0.and.ktau == 0 ) then
        ldosrf = .false.
        ldorad = .false.
        ldoche = .false.
        lskipsrf = .true.
        lskiprad = .true.
        lskipche = .true.
      end if
!
#ifdef MPP1
!
!-----output for dataflow analyses:
!
      if ( iftape ) then
        if ( ldoatm ) then
!=======================================================================
!         gather  ua,va,ta,qva,qca,rainc,rainnc,tgb2d,swt2d,olcd2d,rno2d
          do j = 1 , jendl
            do k = 1 , kz
              do i = 1 , iy
                atm0(i,k,j) = atm1%u(i,k,j)
                atm0(i,k+kz,j) = atm1%v(i,k,j)
                atm0(i,k+kz*2,j) = omega(i,k,j)
                atm0(i,k+kz*3,j) = atm1%t(i,k,j)
                atm0(i,k+kz*4,j) = atm1%qv(i,k,j)
                atm0(i,k+kz*5,j) = atm1%qc(i,k,j)
              end do
            end do
            do i = 1 , iy
              atm0(i,1+kz*6,j) = sps1%ps(i,j)
              atm0(i,2+kz*6,j) = sfsta%rainc(i,j)
              atm0(i,3+kz*6,j) = sfsta%rainnc(i,j)
            end do
          end do
          do j = 1 , jendx
            do n = 1 , nnsg
              do i = 1 , iym1
                atm0(i,3+kz*6+n,j) = ocld2d(n,i,j)
                atm0(i,3+kz*6+n+nnsg,j) = tgb2d(n,i,j)
                atm0(i,3+kz*6+n+nnsg*2,j) = swt2d(n,i,j)
                atm0(i,3+kz*6+n+nnsg*3,j) = rno2d(n,i,j)
              end do
            end do
          end do
          call mpi_gather(atm0, iy*(kz*6+3+nnsg*4)*jxp,mpi_real8,&
                        & atm_0,iy*(kz*6+3+nnsg*4)*jxp,mpi_real8,&
                        & 0,mpi_comm_world,ierr)
          if ( myid == 0 ) then
            rainc_io  = 0.0D0
            rainnc_io = 0.0D0
            do j = 1 , jx
              do k = 1 , kz
                do i = 1 , iy
                  atm1_io%u(i,k,j) = atm_0(i,k,j)
                  atm1_io%v(i,k,j) = atm_0(i,k+kz,j)
                  omega_io(i,k,j) = atm_0(i,k+kz*2,j)
                  atm1_io%t(i,k,j) = atm_0(i,k+kz*3,j)
                  atm1_io%qv(i,k,j) = atm_0(i,k+kz*4,j)
                  atm1_io%qc(i,k,j) = atm_0(i,k+kz*5,j)
                end do
              end do
              do i = 1 , iy
                psa_io(i,j) = atm_0(i,1+kz*6,j)
                if (atm_0(i,2+kz*6,j) > 1D-30) &
                  rainc_io(i,j) = atm_0(i,2+kz*6,j)
                if (atm_0(i,3+kz*6,j) > 1D-30) &
                  rainnc_io(i,j) = atm_0(i,3+kz*6,j)
              end do
            end do
#ifdef BAND
            do j = 1 , jx
#else
            do j = 1 , jxm1
#endif
              do n = 1 , nnsg
                do i = 1 , iym1
                  ocld2d_io(n,i,j) = atm_0(i,3+kz*6+n,j)
                  tgb2d_io(n,i,j) = atm_0(i,3+kz*6+n+nnsg,j)
                  swt2d_io(n,i,j) = atm_0(i,3+kz*6+n+nnsg*2,j)
                  rno2d_io(n,i,j) = atm_0(i,3+kz*6+n+nnsg*3,j)
                end do
              end do
            end do
            call outatm
          end if
          do j = 1 , jendx
            do i = 1 , iym1
              do n = 1 , nnsg
                rno2d(n,i,j) = 0.0D0
              end do
              sfsta%rainc(i,j) = 0.0D0
              sfsta%rainnc(i,j) = 0.0D0
            end do
          end do
        end if
      end if
 
!     Call surface output
 
      if ( ifbat ) then
        if ( ldosrf ) then
          if ( lakemod == 1 .and. iflak .and. mod(iolak,klak) == 0) then
           call lakegather
          end if
          if ( iseaice == 1 ) then
            do j = 1 , jendx
              do n = 1 , nnsg
                do i = 1 , iym1
                  var2d0(i,n,j) = ocld2d(n,i,j)
                end do
              end do
            end do
            call mpi_gather(var2d0, iy*nnsg*jxp,mpi_real8, &
                          & var2d_0,iy*nnsg*jxp,mpi_real8, &
                          & 0,mpi_comm_world,ierr)
            if (myid == 0) then
#ifdef BAND
              do j = 1 , jx
#else
              do j = 1 , jxm1
#endif
                do n = 1 , nnsg
                  do i = 1 , iym1
                    ocld2d_io(n,i,j) = var2d_0(i,n,j)
                  end do
                end do
              end do
            end if
          end if
          do j = 1 , jendx
            do l = 1 , numbat
              do i = 1 , iym2
                bat0(i,l,j) = fbat(j,i,l)
              end do
            end do
          end do
          call mpi_gather(bat0, iym2*numbat*jxp,mpi_real4,       &
                        & bat_0,iym2*numbat*jxp,mpi_real4,       &
                        & 0,mpi_comm_world,ierr)
          if ( myid == 0 ) then
            do l = 1 , numbat
              do i = 1 , iym2
#ifdef BAND
                do j = 1 , jx
                  fbat_io(j,i,l) = bat_0(i,l,j)
#else
                do j = 1 , jxm2
                  fbat_io(j,i,l) = bat_0(i,l,j+1)
#endif
                end do
              end do
            end do
            call outsrf
          end if
          iolak = iolak + 1

          do i = 1 , iym2
            do j = 1 , jxp
              tgmx_o(j,i) = -1.E30
              t2mx_o(j,i) = -1.E30
              tgmn_o(j,i) =  1.E30
              t2mn_o(j,i) =  1.E30
              w10x_o(j,i) = -1.E30
              psmn_o(j,i) =  1.E30
            end do
          end do

          if ( ifsub .and. nsg > 1 ) then

            do j = 1 , jxp
              do l = 1 , numsub
                do n = 1 , nnsg
                  do i = 1 , iym2
                    sub0(i,n,l,j) = fsub(n,j,i,l)
                  end do
                end do
              end do
            end do
            call mpi_gather(sub0, iym2*nnsg*numsub*jxp,mpi_real4, &
                          & sub_0,iym2*nnsg*numsub*jxp,mpi_real4, &
                          & 0,mpi_comm_world,ierr)

            if ( myid == 0 ) then
              do l = 1 , numsub
#ifdef BAND
                do j = 1 , jx
                  do n = 1 , nnsg
                    do i = 1 , iym2
                      fsub_io(n,j,i,l) = sub_0(i,n,l,j)
#else
                do j = 1 , jxm2
                  do n = 1 , nnsg
                    do i = 1 , iym2
                      fsub_io(n,j,i,l) = sub_0(i,n,l,j+1)
#endif
                    end do
                  end do
                end do
              end do
              call outsub
            end if
          end if

        end if
      end if
 
!     Call radiation output
      if ( ifrad ) then
        if ( ldorad ) then
!=======================================================================
!         frad2d, frad3d , psa
          do n = 1 , nrad2d
            do j = 1 , jxp
              do i = 1 , iym2
                rad0(i,n,j) = frad2d(j,i,n)
              end do
            end do
          end do
          do n = 1 , nrad3d
            do k = 1 , kz
              do j = 1 , jxp
                do i = 1 , iym2
                  rad0(i,nrad2d+(n-1)*kz+k,j) = frad3d(j,i,k,n)
                end do
              end do
            end do
          end do
          call mpi_gather(rad0,iym2*(nrad3d*kz+nrad2d)*jxp,      &
                        & mpi_real4,rad_0,(iym2)                 &
                        & *(nrad3d*kz+nrad2d)*jxp,mpi_real4,0,   &
                        & mpi_comm_world,ierr)
         call mpi_gather(sps1%ps(:,1:jxp), iy*jxp,mpi_real8,     &
                        & psa_io,iy*jxp,mpi_real8,               &
                        & 0,mpi_comm_world,ierr)
          if ( myid == 0 ) then
            do n = 1 , nrad2d
#ifdef BAND
              do j = 1 , jx
                do i = 1 , iym2
                  frad2d_io(j,i,n) = rad_0(i,n,j)
#else
              do j = 1 , jxm2
                do i = 1 , iym2
                  frad2d_io(j,i,n) = rad_0(i,n,j+1)
#endif
                end do
              end do
            end do
            do n = 1 , nrad3d
              do k = 1 , kz
#ifdef BAND
                do j = 1 , jx
                  do i = 1 , iym2
                    frad3d_io(j,i,k,n) = rad_0(i,nrad2d+(n-1)*kz+k,j)
#else
                do j = 1 , jxm2
                  do i = 1 , iym2
                    frad3d_io(j,i,k,n) = rad_0(i,nrad2d+(n-1)*kz+k,j+1)
#endif
                  end do
                end do
              end do
            end do
            call outrad
          end if
        end if
      end if
 
!chem2
!     Call chem output
      if ( ifchem ) then
        if ( ldoche ) then
          do j = 1 , jendl
            do n = 1 , ntr
              do k = 1 , kz
                do i = 1 , iy
                  chem0(i,(n-1)*kz+k,j) = chia(i,k,j,n)
                end do
              end do
            end do
          end do
          do j = 1 , jendx
            do k = 1 , kz
              do i = 1 , iym1
                chem0(i,ntr*kz+k,j) = aerext(i,k,j)
                chem0(i,ntr*kz+kz+k,j) = aerssa(i,k,j)
                chem0(i,ntr*kz+kz*2+k,j) = aerasp(i,k,j)
              end do
            end do
          end do
          do j = 1 , jendl
            do n = 1 , ntr
              do i = 1 , iy
                chem0(i,(ntr+3)*kz+n,j) = dtrace(i,j,n)
                chem0(i,(ntr+3)*kz+ntr+n,j) = wdlsc(i,j,n)
                chem0(i,(ntr+3)*kz+ntr*2+n,j) = wdcvc(i,j,n)
                chem0(i,(ntr+3)*kz+ntr*3+n,j) = ddsfc(i,j,n)
                chem0(i,(ntr+3)*kz+ntr*4+n,j) = wxsg(i,j,n)
                chem0(i,(ntr+3)*kz+ntr*5+n,j) = wxaq(i,j,n)
                chem0(i,(ntr+3)*kz+ntr*6+n,j) = cemtrac(i,j,n)
              end do
            end do
          end do
          do j = 1 , jendx
            do i = 1 , iym1
              chem0(i,(ntr+3)*kz+ntr*7+1,j) = aertarf(i,j)
              chem0(i,(ntr+3)*kz+ntr*7+2,j) = aersrrf(i,j)
              chem0(i,(ntr+3)*kz+ntr*7+3,j) = aertalwrf(i,j)
              chem0(i,(ntr+3)*kz+ntr*7+4,j) = aersrlwrf(i,j)             

            end do
          end do
          do j = 1 , jendl
            do i = 1 , iy
              chem0(i,(ntr+3)*kz+ntr*7+5,j) = sps1%ps(i,j)
            end do
          end do
          call mpi_gather(chem0,iy*((ntr+3)*kz+ntr*7+5)*jxp,            &
                        & mpi_real8,chem_0,iy*((ntr+3)*kz+ntr*7+5)*jxp, &
                        & mpi_real8,0,mpi_comm_world,ierr)
          if ( myid == 0 ) then
            do j = 1 , jx
              do n = 1 , ntr
                do k = 1 , kz
                  do i = 1 , iy
                    chia_io(i,k,j,n) = chem_0(i,(n-1)*kz+k,j)
                  end do
                end do
              end do
            end do
#ifdef BAND
            do j = 1 , jx
              do k = 1 , kz
                do i = 1 , iym1
                  aerext_io(i,k,j) = chem_0(i,ntr*kz+k,j)
                  aerssa_io(i,k,j) = chem_0(i,ntr*kz+kz+k,j)
                  aerasp_io(i,k,j) = chem_0(i,ntr*kz+kz*2+k,j)
#else
            do j = 1 , jxm1
              do k = 1 , kz
                do i = 1 , iym1
                  aerext_io(i,k,j) = chem_0(i,ntr*kz+k,j+1)
                  aerssa_io(i,k,j) = chem_0(i,ntr*kz+kz+k,j+1)
                  aerasp_io(i,k,j) = chem_0(i,ntr*kz+kz*2+k,j+1)
#endif
                end do
              end do
            end do
            do j = 1 , jx
              do n = 1 , ntr
                do i = 1 , iy
                  dtrace_io(i,j,n) = chem_0(i,(ntr+3)*kz+n,j)
                  wdlsc_io(i,j,n) = chem_0(i,(ntr+3)*kz+ntr+n,j)
                  wdcvc_io(i,j,n) = chem_0(i,(ntr+3)*kz+ntr*2+n,j)
                  ddsfc_io(i,j,n) = chem_0(i,(ntr+3)*kz+ntr*3+n,j)
                  wxsg_io(i,j,n) = chem_0(i,(ntr+3)*kz+ntr*4+n,j)
                  wxaq_io(i,j,n) = chem_0(i,(ntr+3)*kz+ntr*5+n,j)
                  cemtrac_io(i,j,n) = chem_0(i,(ntr+3)*kz+ntr*6+n,j)
                end do
              end do
            end do
#ifdef BAND
            do j = 1 , jx
              do i = 1 , iym1
                aertarf_io(i,j) = chem_0(i,(ntr+3)*kz+ntr*7+1,j)
                aersrrf_io(i,j) = chem_0(i,(ntr+3)*kz+ntr*7+2,j)
                aertalwrf_io(i,j) = chem_0(i,(ntr+3)*kz+ntr*7+3,j)
                aersrlwrf_io(i,j) = chem_0(i,(ntr+3)*kz+ntr*7+4,j)
#else
            do j = 1 , jxm1
              do i = 1 , iym1
                aertarf_io(i,j) = chem_0(i,(ntr+3)*kz+ntr*7+1,j+1)
                aersrrf_io(i,j) = chem_0(i,(ntr+3)*kz+ntr*7+2,j+1)
                aertalwrf_io(i,j) = chem_0(i,(ntr+3)*kz+ntr*7+3,j+1)
                aersrlwrf_io(i,j) = chem_0(i,(ntr+3)*kz+ntr*7+4,j+1)
#endif
              end do
            end do
            do j = 1 , jx
              do i = 1 , iy
                psa_io(i,j) = chem_0(i,(ntr+3)*kz+ntr*7+5,j)
              end do
            end do
            call outche
            remlsc_io  = 0.0D0
            remcvc_io  = 0.0D0
            rxsg_io    = 0.0D0
            rxsaq1_io  = 0.0D0
            rxsaq2_io  = 0.0D0
            cemtr_io   = 0.0D0
            remdrd_io  = 0.0D0
            wdlsc_io   = 0.0D0
            wdcvc_io   = 0.0D0
            ddsfc_io   = 0.0D0
            wxsg_io    = 0.0D0
            wxaq_io    = 0.0D0
            cemtrac_io = 0.0D0
            aertarf_io = 0.0D0
            aersrrf_io = 0.0D0
            aersrlwrf_io=0.0D0
            aertalwrf_io=0.0D0
          end if
          do n = 1 , ntr
            do j = 1 , jendl
              do k = 1 , kz
                do i = 1 , iy
                  remlsc(i,k,j,n) = 0.0D0
                  remcvc(i,k,j,n) = 0.0D0
                  rxsg(i,k,j,n) = 0.0D0
                  rxsaq1(i,k,j,n) = 0.0D0
                  rxsaq2(i,k,j,n) = 0.0D0
                end do
              end do
            end do
          end do
          do n = 1 , ntr
            do j = 1 , jendl
              do i = 1 , iy
                cemtr(i,j,n) = 0.0D0
                remdrd(i,j,n) = 0.0D0
                wdlsc(i,j,n) = 0.0D0
                wdcvc(i,j,n) = 0.0D0
                ddsfc(i,j,n) = 0.0D0
                wxsg(i,j,n) = 0.0D0
                wxaq(i,j,n) = 0.0D0
                cemtrac(i,j,n) = 0.0D0
              end do
            end do
          end do
          do j = 1 , jendl
            do i = 1 , iym1
              aertarf(i,j) = 0.0D0
              aersrrf(i,j) = 0.0D0
              aertalwrf(i,j) = 0.0D0              
              aersrlwrf(i,j) = 0.0D0
            end do
          end do
        end if
      end if
!
!-----output for restart:
!
      if ( ifsave ) then
        if ( ldosav .or. ldotmp ) then
          if ( lakemod == 1 ) call lakegather
          do j = 1 , jendl
            do k = 1 , kz
              do i = 1 , iy
                sav0(i,k,j) = ub0(i,k,j)
                sav0(i,kz+k,j) = vb0(i,k,j)
                sav0(i,kz*2+k,j) = qb0(i,k,j)
                sav0(i,kz*3+k,j) = tb0(i,k,j)
              end do
            end do
            do i = 1 , iy
              sav0(i,kz*4+1,j) = ps0(i,j)
              sav0(i,kz*4+2,j) = ts0(i,j)
            end do
          end do
          allrec = kz*4 + 2
          call mpi_gather(sav0, iy*allrec*jxp,mpi_real8,         &
                        & sav_0,iy*allrec*jxp,mpi_real8,         &
                        & 0,mpi_comm_world,ierr)
          if ( myid == 0 ) then
            do j = 1 , jx
              do k = 1 , kz
                do i = 1 , iy
                  ub0_io(i,k,j) = sav_0(i,k,j)
                  vb0_io(i,k,j) = sav_0(i,kz+k,j)
                  qb0_io(i,k,j) = sav_0(i,kz*2+k,j)
                  tb0_io(i,k,j) = sav_0(i,kz*3+k,j)
                end do
              end do
              do i = 1 , iy
                ps0_io(i,j) = sav_0(i,kz*4+1,j)
                ts0_io(i,j) = sav_0(i,kz*4+2,j)
              end do
            end do
          end if
          if ( ehso4 ) then
            do j = 1 , jendl
              do k = 1 , kz
                do i = 1 , iy
                  sav0s(i,k,j) = so0(i,k,j)
                end do
              end do
            end do
            call mpi_gather(sav0s, iy*kz*jxp,mpi_real8,          &
                          & sav_0s,iy*kz*jxp,mpi_real8,          &
                          & 0,mpi_comm_world,ierr)
            if ( myid == 0 ) then
              do j = 1 , jx
                do k = 1 , kz
                  do i = 1 , iy
                    so0_io(i,k,j) = sav_0s(i,k,j)
                  end do
                end do
              end do
            end if
          end if
          do j = 1 , jendl
            do k = 1 , kz
              do i = 1 , iy
                sav0(i,k,j) = atm1%u(i,k,j)
                sav0(i,kz+k,j) = atm2%u(i,k,j)
                sav0(i,kz*2+k,j) = atm1%v(i,k,j)
                sav0(i,kz*3+k,j) = atm2%v(i,k,j)
              end do
            end do
            do i = 1 , iy
              sav0(i,kz*4+1,j) = sps1%ps(i,j)
              sav0(i,kz*4+2,j) = sps2%ps(i,j)
            end do
          end do
          allrec = kz*4 + 2
          call mpi_gather(sav0, iy*allrec*jxp,mpi_real8,         &
                        & sav_0,iy*allrec*jxp,mpi_real8,         &
                        & 0,mpi_comm_world,ierr)
          if ( myid == 0 ) then
            do j = 1 , jx
              do k = 1 , kz
                do i = 1 , iy
                  atm1_io%u(i,k,j) = sav_0(i,k,j)
                  atm2_io%u(i,k,j) = sav_0(i,kz+k,j)
                  atm1_io%v(i,k,j) = sav_0(i,kz*2+k,j)
                  atm2_io%v(i,k,j) = sav_0(i,kz*3+k,j)
                end do
              end do
              do i = 1 , iy
                psa_io(i,j) = sav_0(i,kz*4+1,j)
                psb_io(i,j) = sav_0(i,kz*4+2,j)
              end do
            end do
          end if
          do j = 1 , jendl
            do k = 1 , kz
              do i = 1 , iy
                sav0(i,k,j) = atm1%t(i,k,j)
                sav0(i,kz+k,j) = atm2%t(i,k,j)
                sav0(i,kz*2+k,j) = atm1%qv(i,k,j)
                sav0(i,kz*3+k,j) = atm2%qv(i,k,j)
              end do
            end do
            do i = 1 , iy
              sav0(i,kz*4+1,j) = sts1%tg(i,j)
              sav0(i,kz*4+2,j) = sts2%tg(i,j)
            end do
          end do
          allrec = kz*4 + 2
          call mpi_gather(sav0, iy*allrec*jxp,mpi_real8,         &
                        & sav_0,iy*allrec*jxp,mpi_real8,         &
                        & 0,mpi_comm_world,ierr)
          if ( myid == 0 ) then
            do j = 1 , jx
              do k = 1 , kz
                do i = 1 , iy
                  atm1_io%t(i,k,j) = sav_0(i,k,j)
                  atm2_io%t(i,k,j) = sav_0(i,kz+k,j)
                  atm1_io%qv(i,k,j) = sav_0(i,kz*2+k,j)
                  atm2_io%qv(i,k,j) = sav_0(i,kz*3+k,j)
                end do
              end do
              do i = 1 , iy
                tga_io(i,j) = sav_0(i,kz*4+1,j)
                tgb_io(i,j) = sav_0(i,kz*4+2,j)
              end do
            end do
          end if
          do j = 1 , jendl
            do k = 1 , kz
              do i = 1 , iy
                sav0(i,k,j) = atm1%qc(i,k,j)
                sav0(i,kz+k,j) = atm2%qc(i,k,j)
                sav0(i,kz*2+k,j) = fcc(i,k,j)
              end do
            end do
            do i = 1 , iy
              sav0(i,kz*4+1,j) = sfsta%rainc(i,j)
              sav0(i,kz*4+2,j) = sfsta%rainnc(i,j)
            end do
          end do
          do j = 1 , jendx
            do k = 1 , kz
              do i = 1 , iym1
                sav0(i,kz*3+k,j) = heatrt(i,k,j)
              end do
            end do
          end do
          allrec = kz*4 + 2
          call mpi_gather(sav0, iy*allrec*jxp,mpi_real8,         &
                        & sav_0,iy*allrec*jxp,mpi_real8,         &
                        & 0,mpi_comm_world,ierr)
          if ( myid == 0 ) then
            rainc_io  = 0.0D0
            rainnc_io = 0.0D0
            do j = 1 , jx
              do k = 1 , kz
                do i = 1 , iy
                  atm1_io%qc(i,k,j) = sav_0(i,k,j)
                  atm2_io%qc(i,k,j) = sav_0(i,kz+k,j)
                  fcc_io(i,k,j) = sav_0(i,kz*2+k,j)
                end do
              end do
              do i = 1 , iy
                if (sav_0(i,kz*4+1,j) > 1D-30) &
                  rainc_io(i,j) = sav_0(i,kz*4+1,j)
                if (sav_0(i,kz*4+2,j) > 1D-30) &
                  rainnc_io(i,j) = sav_0(i,kz*4+2,j)
              end do
            end do
#ifdef BAND
            do j = 1 , jx
#else
            do j = 1 , jxm1
#endif
              do k = 1 , kz
                do i = 1 , iym1
                  heatrt_io(i,k,j) = sav_0(i,kz*3+k,j)
                end do
              end do
            end do
          end if
          do j = 1 , jendl
            do i = 1 , iy
              sav0a(i,1,j) = sfsta%hfx(i,j)
              sav0a(i,2,j) = sfsta%qfx(i,j)
              sav0a(i,3,j) = sfsta%uvdrag(i,j)
              sav0a(i,4,j) = sfsta%tgbb(i,j)
            end do
            do n = 1 , nnsg
              do i = 1 , iy
                sav0a(i,4+n,j) = snowc(n,i,j)
              end do
            end do
          end do
          do j = 1 , jendx
            do k = 1 , kzp1
              do i = 1 , iym1
                sav0a(i,nnsg+4+k,j) = o3prof(i,k,j)
              end do
            end do
          end do
          allrec = 4 + nnsg + kzp1
          call mpi_gather(sav0a, iy*allrec*jxp,mpi_real8,        &
                        & sav_0a,iy*allrec*jxp,mpi_real8,        &
                        & 0,mpi_comm_world,ierr)
          if ( myid == 0 ) then
            do j = 1 , jx
              do i = 1 , iy
                hfx_io(i,j) = sav_0a(i,1,j)
                qfx_io(i,j) = sav_0a(i,2,j)
                uvdrag_io(i,j) = sav_0a(i,3,j)
                tgbb_io(i,j) = sav_0a(i,4,j)
              end do
              do n = 1 , nnsg
                do i = 1 , iy
                  snowc_io(n,i,j) = sav_0a(i,4+n,j)
                end do
              end do
            end do
#ifdef BAND
            do j = 1 , jx
#else
            do j = 1 , jxm1
#endif
              do k = 1 , kzp1
                do i = 1 , iym1
                  o3prof_io(i,k,j) = sav_0a(i,4+nnsg+k,j)
                end do
              end do
            end do
          end if
          if ( iocnflx == 2 )                                      &
             & call mpi_gather(sfsta%zpbl,   iy*jxp,mpi_real8,     &
             &                 zpbl_io,iy*jxp,mpi_real8,           &
             &                 0,mpi_comm_world,ierr)
          if ( icup == 1 ) then
            do j = 1 , jendl
              do k = 1 , kz
                do i = 1 , iy
                  sav0c(i,k,j) = rsheat(i,k,j)
                  sav0c(i,kz+k,j) = rswat(i,k,j)
                end do
              end do
            end do
            allrec = kz*2
            call mpi_gather(sav0c, iy*allrec*jxp,mpi_real8,      &
                          & sav_0c,iy*allrec*jxp,mpi_real8,      &
                          & 0,mpi_comm_world,ierr)
            if ( myid == 0 ) then
              do j = 1 , jx
                do k = 1 , kz
                  do i = 1 , iy
                    rsheat_io(i,k,j) = sav_0c(i,k,j)
                    rswat_io(i,k,j) = sav_0c(i,kz+k,j)
                  end do
                end do
              end do
            end if
          end if
          if ( icup == 3 ) then
            do j = 1 , jendl
              do k = 1 , kz
                do i = 1 , iy
                  sav0b(i,k,j) = tbase(i,k,j)
                end do
              end do
              do i = 1 , iy
                sav0b(i,kzp1,j) = cldefi(i,j)
              end do
            end do
            allrec = kzp1
            call mpi_gather(sav0b, iy*allrec*jxp,mpi_real8,      &
                          & sav_0b,iy*allrec*jxp,mpi_real8,      &
                          & 0,mpi_comm_world,ierr)
            if ( myid == 0 ) then
              do j = 1 , jx
                do k = 1 , kz
                  do i = 1 , iy
                    tbase_io(i,k,j) = sav_0b(i,k,j)
                  end do
                end do
                do i = 1 , iy
                  cldefi_io(i,j) = sav_0b(i,kzp1,j)
                end do
              end do
            end if
          end if
          if ( icup==4 .or. icup==99 .or. icup==98 ) then
            call mpi_gather(cbmf2d,   iy*jxp,mpi_real8,            &
                          & cbmf2d_io,iy*jxp,mpi_real8,            &
                          & 0,mpi_comm_world,ierr)
          end if
          do j = 1 , jendx
            do l = 1 , 4
              do k = 1 , kz
                do i = 1 , iym1
                  sav1(i,(l-1)*kz+k,j) = absnxt(i,k,l,j)
                end do
              end do
            end do
          end do
          allrec = kz*4
          do j = 1 , jendx
            do l = 1 , kzp1
              do k = 1 , kzp1
                do i = 1 , iym1
                  sav1(i,allrec+(l-1)*(kzp1)+k,j) = abstot(i,k,l,j)
                end do
              end do
            end do
          end do
          allrec = allrec + (kzp1)*(kz+1)
          do j = 1 , jendx
            do k = 1 , kzp1
              do i = 1 , iym1
                sav1(i,allrec+k,j) = emstot(i,k,j)
              end do
            end do
          end do
          allrec = kz*4+(kzp1*kzp2)
          call mpi_gather(sav1, iym1*allrec*jxp,mpi_real8,       &
                        & sav_1,iym1*allrec*jxp,mpi_real8,       &
                        & 0,mpi_comm_world,ierr)
          if ( myid == 0 ) then
#ifdef BAND
            do j = 1 , jx
#else
            do j = 1 , jxm1
#endif
              do l = 1 , 4
                do k = 1 , kz
                  do i = 1 , iym1
                    absnxt_io(i,k,l,j) = sav_1(i,(l-1)*kz+k,j)
                  end do
                end do
              end do
            end do
            allrec = kz*4
#ifdef BAND
            do j = 1 , jx
#else
            do j = 1 , jxm1
#endif
              do l = 1 , kzp1
                do k = 1 , kzp1
                  do i = 1 , iym1
                    abstot_io(i,k,l,j)                                  &
                    & = sav_1(i,allrec+(l-1)*(kzp1)+k,j)
                  end do
                end do
              end do
            end do
            allrec = allrec + (kzp1)*(kz+1)
#ifdef BAND
            do j = 1 , jx
#else
            do j = 1 , jxm1
#endif
              do k = 1 , kzp1
                do i = 1 , iym1
                  emstot_io(i,k,j) = sav_1(i,allrec+k,j)
                end do
              end do
            end do
          end if
          do j = 1 , jendx
            do n = 1 , nnsg
              do i = 1 , iym1
                sav2(i,n,j) = taf2d(n,i,j)
                sav2(i,nnsg+n,j) = tlef2d(n,i,j)
                sav2(i,nnsg*2+n,j) = ssw2d(n,i,j)
                sav2(i,nnsg*3+n,j) = srw2d(n,i,j)
                sav2(i,nnsg*4+n,j) = col2d(n,i,j)
              end do
            end do
            do i = 1 , iym1
              sav2(i,nnsg*5+1,j) = sol2d(i,j)
              sav2(i,nnsg*5+2,j) = solvd2d(i,j)
              sav2(i,nnsg*5+3,j) = solvs2d(i,j)
              sav2(i,nnsg*5+4,j) = flw2d(i,j)
            end do
          end do
          allrec = nnsg*5 + 4
          call mpi_gather(sav2, iym1*allrec*jxp,mpi_real8,       &
                        & sav_2,iym1*allrec*jxp,mpi_real8,       &
                        & 0,mpi_comm_world,ierr)
          if ( myid == 0 ) then
#ifdef BAND
            do j = 1 , jx
#else
            do j = 1 , jxm1
#endif
              do n = 1 , nnsg
                do i = 1 , iym1
                  taf2d_io(n,i,j) = sav_2(i,n,j)
                  tlef2d_io(n,i,j) = sav_2(i,nnsg+n,j)
                  ssw2d_io(n,i,j) = sav_2(i,nnsg*2+n,j)
                  srw2d_io(n,i,j) = sav_2(i,nnsg*3+n,j)
                  col2d_io(n,i,j) = sav_2(i,nnsg*4+n,j)
                end do
              end do
              do i = 1 , iym1
                sol2d_io(i,j) = sav_2(i,nnsg*5+1,j)
                solvd2d_io(i,j) = sav_2(i,nnsg*5+2,j)
                solvs2d_io(i,j) = sav_2(i,nnsg*5+3,j)
                flw2d_io(i,j) = sav_2(i,nnsg*5+4,j)
              end do
            end do
          end if
#ifdef CLM
          do j = 1 , jendx
            do i = 1 , iym1
              sav_clmin(i,1,j)  = sols2d(i,j)
              sav_clmin(i,2,j)  = soll2d(i,j)
              sav_clmin(i,3,j)  = solsd2d(i,j)
              sav_clmin(i,4,j)  = solld2d(i,j)
              sav_clmin(i,5,j)  = aldirs2d(i,j)
              sav_clmin(i,6,j)  = aldirl2d(i,j)
              sav_clmin(i,7,j)  = aldifs2d(i,j)
              sav_clmin(i,8,j)  = aldifl2d(i,j)
            end do
          end do
          call mpi_gather(sav_clmin, iym1*8*jxp,mpi_real8,       &
                        & sav_clmout,iym1*8*jxp,mpi_real8,       &
                        & 0,mpi_comm_world,ierr)
          if ( myid == 0 ) then
#ifdef BAND
            do j = 1 , jx
#else
            do j = 1 , jxm1
#endif
              do i = 1 , iym1
                sols2d_io(i,j)   = sav_clmout(i,1,j)
                soll2d_io(i,j)   = sav_clmout(i,2,j)
                solsd2d_io(i,j)  = sav_clmout(i,3,j)
                solld2d_io(i,j)  = sav_clmout(i,4,j)
                aldirs2d_io(i,j) = sav_clmout(i,5,j)
                aldirl2d_io(i,j) = sav_clmout(i,6,j)
                aldifs2d_io(i,j) = sav_clmout(i,7,j)
                aldifl2d_io(i,j) = sav_clmout(i,8,j)
              end do
            end do
          end if
          call mpi_gather(satbrt2d,   iy*jxp,mpi_real8, &
                        & satbrt2d_io,iy*jxp,mpi_real8, &
                        & 0,mpi_comm_world,ierr)
#endif
          do j = 1 , jendx
            do n = 1 , nnsg
              do i = 1 , iym1
                sav2(i,n,j) = tgb2d(n,i,j)
                sav2(i,nnsg+n,j) = swt2d(n,i,j)
                sav2(i,nnsg*2+n,j) = scv2d(n,i,j)
                sav2(i,nnsg*3+n,j) = gwet2d(n,i,j)
                sav2(i,nnsg*4+n,j) = tg2d(n,i,j)
              end do
            end do
            do i = 1 , iym1
              sav2(i,nnsg*5+1,j) = flwd2d(i,j)
              sav2(i,nnsg*5+2,j) = fsw2d(i,j)
              sav2(i,nnsg*5+3,j) = sabv2d(i,j)
              sav2(i,nnsg*5+4,j) = sinc2d(i,j)
            end do
          end do
          allrec = nnsg*5 + 4
          call mpi_gather(sav2, iym1*allrec*jxp,mpi_real8,       &
                        & sav_2,iym1*allrec*jxp,mpi_real8,       &
                        & 0,mpi_comm_world,ierr)
          if ( myid == 0 ) then
#ifdef BAND
            do j = 1 , jx
#else
            do j = 1 , jxm1
#endif
              do n = 1 , nnsg
                do i = 1 , iym1
                  tgb2d_io(n,i,j) = sav_2(i,n,j)
                  swt2d_io(n,i,j) = sav_2(i,nnsg+n,j)
                  scv2d_io(n,i,j) = sav_2(i,nnsg*2+n,j)
                  gwet2d_io(n,i,j) = sav_2(i,nnsg*3+n,j)
                  tg2d_io(n,i,j) = sav_2(i,nnsg*4+n,j)
                end do
              end do
              do i = 1 , iym1
                flwd2d_io(i,j) = sav_2(i,nnsg*5+1,j)
                fsw2d_io(i,j) = sav_2(i,nnsg*5+2,j)
                sabv2d_io(i,j) = sav_2(i,nnsg*5+3,j)
                sinc2d_io(i,j) = sav_2(i,nnsg*5+4,j)
              end do
            end do
          end if
          do j = 1 , jendx
            do n = 1 , nnsg
              do i = 1 , iym1
                sav2(i,n,j) = veg2d1(n,i,j)
                sav2(i,nnsg+n,j) = sag2d(n,i,j)
                sav2(i,nnsg*2+n,j) = sice2d(n,i,j)
                sav2(i,nnsg*3+n,j) = dew2d(n,i,j)
                sav2(i,nnsg*4+n,j) = ocld2d(n,i,j)
              end do
            end do
            do i = 1 , iym1
              sav2(i,nnsg*5+1,j) = pptnc(i,j)
              sav2(i,nnsg*5+2,j) = pptc(i,j)
              sav2(i,nnsg*5+3,j) = prca2d(i,j)
              sav2(i,nnsg*5+4,j) = prnca2d(i,j)
            end do
          end do
          allrec = nnsg*5 + 4
          call mpi_gather(sav2, iym1*allrec*jxp,mpi_real8,       &
                        & sav_2,iym1*allrec*jxp,mpi_real8,       &
                        & 0,mpi_comm_world,ierr)
          if ( myid == 0 ) then
#ifdef BAND
            do j = 1 , jx
#else
            do j = 1 , jxm1
#endif
              do n = 1 , nnsg
                do i = 1 , iym1
                  veg2d1_io(n,i,j) = sav_2(i,n,j)
                  sag2d_io(n,i,j) = sav_2(i,nnsg+n,j)
                  sice2d_io(n,i,j) = sav_2(i,nnsg*2+n,j)
                  dew2d_io(n,i,j) = sav_2(i,nnsg*3+n,j)
                  ocld2d_io(n,i,j) = sav_2(i,nnsg*4+n,j)
                end do
              end do
              do i = 1 , iym1
                pptnc_io(i,j) = sav_2(i,nnsg*5+1,j)
                pptc_io(i,j) = sav_2(i,nnsg*5+2,j)
                prca2d_io(i,j) = sav_2(i,nnsg*5+3,j)
                prnca2d_io(i,j) = sav_2(i,nnsg*5+4,j)
              end do
            end do
          end if
          do j = 1 , jendx
            do n = 1 , nnsg
              do i = 1 , iym1
                sav2a(i,n,j) = ircp2d(n,i,j)
                sav2a(i,nnsg+n,j) = text2d(n,i,j)
              end do
            end do
            do i = 1 , iym1
              sav2a(i,nnsg*2+1,j) = veg2d(i,j)
            end do
          end do
          allrec = nnsg*2 + 1
          call mpi_gather(sav2a, iym1*allrec*jxp,mpi_real8,      &
                        & sav_2a,iym1*allrec*jxp,mpi_real8,      &
                        & 0,mpi_comm_world,ierr)
          if ( myid == 0 ) then
#ifdef BAND
            do j = 1 , jx
#else
            do j = 1 , jxm1
#endif
              do n = 1 , nnsg
                do i = 1 , iym1
                  ircp2d_io(n,i,j) = sav_2a(i,n,j)
                  text2d_io(n,i,j) = sav_2a(i,nnsg+n,j)
                end do
              end do
              do i = 1 , iym1
                veg2d_io(i,j) = sav_2a(i,nnsg*2+1,j)
              end do
            end do
          end if
 
          if ( ichem == 1 ) then
            do j = 1 , jendl
              do n = 1 , ntr
                do k = 1 , kz
                  do i = 1 , iy
                    sav4(i,(n-1)*kz+k,j) = chia(i,k,j,n)
                    sav4(i,ntr*kz+(n-1)*kz+k,j) = chib(i,k,j,n)
                    sav4(i,ntr*kz*2+(n-1)*kz+k,j) = remlsc(i,k,j,n)
                    sav4(i,ntr*kz*3+(n-1)*kz+k,j) = remcvc(i,k,j,n)
                  end do
                end do
              end do
            end do
            allrec = 4*ntr*kz
            do j = 1 , jendl
              do n = 1 , ntr
                do i = 1 , iy
                  sav4(i,allrec+n,j) = remdrd(i,j,n)
                end do
              end do
            end do
            allrec = ntr*(kz*4+1)
            call mpi_gather(sav4, iy*allrec*jxp,mpi_real8,       &
                          & sav_4,iy*allrec*jxp,mpi_real8,       &
                          & 0,mpi_comm_world,ierr)
            if ( myid == 0 ) then
              do j = 1 , jx
                do n = 1 , ntr
                  do k = 1 , kz
                    do i = 1 , iy
                      chia_io(i,k,j,n) = sav_4(i,(n-1)*kz+k,j)
                      chib_io(i,k,j,n) = sav_4(i,ntr*kz+(n-1)*kz+k,j)
                      remlsc_io(i,k,j,n)                                &
                      & = sav_4(i,ntr*kz*2+(n-1)*kz+k,j)
                      remcvc_io(i,k,j,n)                                &
                      & = sav_4(i,ntr*kz*3+(n-1)*kz+k,j)
                    end do
                  end do
                end do
              end do
              allrec = 4*ntr*kz
              do j = 1 , jx
                do n = 1 , ntr
                  do i = 1 , iy
                    remdrd_io(i,j,n) = sav_4(i,allrec+n,j)
                  end do
                end do
              end do
            end if
            do j = 1 , jendx
              do i = 1 , iym1
                sav4a(i,1,j) = ssw2da(i,j)
                sav4a(i,2,j) = sdeltk2d(i,j)
                sav4a(i,3,j) = sdelqk2d(i,j)
                sav4a(i,4,j) = sfracv2d(i,j)
                sav4a(i,5,j) = sfracb2d(i,j)
                sav4a(i,6,j) = sfracs2d(i,j)
                sav4a(i,7,j) = svegfrac2d(i,j)
              end do
            end do
            call mpi_gather(sav4a, iym1*7*jxp,mpi_real8,                &
                          & sav_4a,iym1*7*jxp,mpi_real8,                &
                          & 0,mpi_comm_world,ierr)
            if ( myid == 0 ) then
#ifdef BAND
              do j = 1 , jx
#else
              do j = 1 , jxm1
#endif
                do i = 1 , iym1
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
              do i = 1 , iy
                sav0d(i,n,j) = spsav%dstor(i,j,n)
                sav0d(i,n+nsplit,j) = spsav%hstor(i,j,n)
              end do
            end do
          end do
          call mpi_gather(sav0d, iy*nsplit*2*jxp,mpi_real8,      &
                        & sav_0d,iy*nsplit*2*jxp,mpi_real8,      &
                        & 0,mpi_comm_world,ierr)
          if ( myid == 0 ) then
            do j = 1 , jx
              do n = 1 , nsplit
                do i = 1 , iy
                  dstor_io(i,j,n) = sav_0d(i,n,j)
                  hstor_io(i,j,n) = sav_0d(i,n+nsplit,j)
                end do
              end do
            end do
          end if
          do j = 1 , jendl
            do k = 1 , kz
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
          call mpi_gather(sav6, kz*8*jxp,mpi_real8,              &
                        & sav_6,kz*8*jxp,mpi_real8,              &
                        & 0,mpi_comm_world,ierr)
          if ( myid == 0 ) then
            do j = 1 , jx
              do k = 1 , kz
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
#ifndef BAND
          call mpi_bcast(ujlx(1,1),iy*kz,mpi_real8,nproc-1,             &
                       & mpi_comm_world,ierr)
          call mpi_bcast(ujl(1,1),iy*kz,mpi_real8,nproc-1,              &
                       & mpi_comm_world,ierr)
          call mpi_bcast(vjlx(1,1),iy*kz,mpi_real8,nproc-1,             &
                       & mpi_comm_world,ierr)
          call mpi_bcast(vjl(1,1),iy*kz,mpi_real8,nproc-1,              &
                       & mpi_comm_world,ierr)
#endif
          if ( ldosav ) then
            call write_savefile(idatex, .false.)
          else
            call write_savefile(idatex, .true.)
          end if
        end if
      end if

#else
!
!-----output for dataflow analyses:
!
      if ( iftape ) then
        if ( ldoatm ) then
          call outatm
        end if
      end if
!
!
!     Call surface output
      if ( ifbat ) then
        if ( ldosrf ) then
          call outsrf
          do i = 1 , iym2
#ifdef BAND
            do j = 1 , jx
#else
            do j = 1 , jxm2
#endif
              tgmx_o(j,i) = -1.E30
              t2mx_o(j,i) = -1.E30
              tgmn_o(j,i) =  1.E30
              t2mn_o(j,i) =  1.E30
              w10x_o(j,i) = -1.E30
              psmn_o(j,i) =  1.E30
            end do
          end do
          if ( ifsub .and. nsg > 1 ) then
            call outsub
          end if
        end if
      end if
! 
! 
!     Call radiation output
      if ( ifrad ) then
        if ( ldorad ) then
          call outrad
        end if
      end if
!
!
!     Call chem output
      if ( ifchem ) then
        if ( ldoche ) then
          call outche
          remlsc    = 0.0D0
          remcvc    = 0.0D0
          rxsg      = 0.0D0
          rxsaq1    = 0.0D0
          rxsaq2    = 0.0D0
          cemtr     = 0.0D0
          remdrd    = 0.0D0
          wdlsc     = 0.0D0
          wdcvc     = 0.0D0
          ddsfc     = 0.0D0
          wxsg      = 0.0D0
          wxaq      = 0.0D0
          cemtrac   = 0.0D0
          aertarf   = 0.0D0
          aersrrf   = 0.0D0
          aersrlwrf = 0.0D0
          aertalwrf = 0.0D0
        end if
      end if
!
!-----output for restart:
!
      if ( ifsave ) then
        if (ldosav) then
          call write_savefile(idatex,.false.)
        else
          call write_savefile(idatex,.true.)
        end if
      end if
!
#endif

#ifdef MPP1
      if ( myid == 0 ) then
#endif        
        if ( lday == 1 .and. lhour == 0 .and. dabs(xtime)<0.00001 ) then
          if ( .not. lstartup .and. idatex /= idate2 ) then
            call mkfile
          end if
        end if
#ifdef MPP1
      end if
#endif
!
      call time_end(subroutine_name,idindx) 

      end subroutine output
!
      subroutine mkfile
 
      implicit none
!
      if (myid /= 0) return

      print * , ' '
      print * , '******* OPENING NEW OUTPUT FILES:' , idatex
      print * , ' '

      if ( iftape ) then
        call prepare_common_out(idatex,'ATM')
      end if
 
      if ( ifbat ) then
        call prepare_common_out(idatex,'SRF')
        if (lakemod == 1 .and. iflak) then
          call prepare_common_out(idatex,'LAK')
        end if
      end if
 
      if ( nsg > 1 .and. ifsub ) then
        call prepare_common_out(idatex,'SUB')
      end if
 
      if ( ifrad ) then
        call prepare_common_out(idatex,'RAD')
      end if
 
      if ( ichem == 1 ) then
        if ( ifchem ) then
          call prepare_common_out(idatex,'CHE')
        end if
      end if

      end subroutine mkfile
!
      subroutine outatm

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!     this subroutine writes the model output to tape or disk for use c
!     in dataflow analyses.                                           c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none
      integer :: jjx , iiy
#ifdef BAND
      jjx = jx
      iiy = iym1
#else
      jjx = jxm1
      iiy = iym1
#endif

#ifdef MPP1
      call writerec_atm(jx,iy,jjx,iiy,kz,nnsg,atm1_io%u,atm1_io%v,  &
              omega_io,atm1_io%t,atm1_io%qv,atm1_io%qc,psa_io,      &
              rainc_io,rainnc_io,tgb2d_io,swt2d_io,rno2d_io,        &
              ocld2d_io,idatex)
#else
      call writerec_atm(jx,iy,jjx,iiy,kz,nnsg,atm1%u,atm1%v,omega, &
                        atm1%t,atm1%qv,atm1%qc,sps1%ps,sfsta%rainc,&
                        sfsta%rainnc,tgb2d,swt2d,rno2d,ocld2d,idatex)
#endif
 
      write (*,*) 'ATM variables written at ' , idatex , xtime
 
      end subroutine outatm
!
      subroutine outsrf

      implicit none
!
! Local variables
!
      integer :: i , j
!
#ifdef BAND
      i = iy-2
      j = jx
#else
      i = iy-2
      j = jx-2
#endif

#ifdef MPP1
      call writerec_srf(j,i,numbat,fbat_io,ocld2d_io,idatex)
#else
      call writerec_srf(j,i,numbat,fbat,ocld2d,idatex)
#endif
      write (*,*) 'SRF variables written at ' , idatex , xtime
 
      if (lakemod == 1 .and. iflak .and. mod(iolak,klak) == 0) then
#ifdef MPP1
        call writerec_lak(j,i,numbat,fbat_io,evl2d_io,aveice2d_io, &
                          hsnow2d_io,tlak3d_io,idatex)
#else
        call writerec_lak(j,i,numbat,fbat,evpa2d,aveice2d,hsnow2d, &
                          tlak3d,idatex)
#endif
        write (*,*) 'LAK variables written at ' , idatex , xtime
      end if

      end subroutine outsrf
!
      subroutine outsub

      implicit none

      integer :: i , j

#ifdef BAND
      i = iym2sg
      j = jxsg
#else
      i = iym2sg
      j = jxm2sg
#endif

#ifdef MPP1
      call writerec_sub(j,i,nsg,numsub,fsub_io,idatex)
#else
      call writerec_sub(j,i,nsg,numsub,fsub,idatex)
#endif

      write (*,*) 'SUB variables written at ' , idatex , xtime

      end subroutine outsub
!
      subroutine outrad
!
      implicit none
!
      integer :: i , j , imax , jmax , istart, jstart
!
!      character (len=50) :: subroutine_name='outrad'
!      integer :: idindx=0
!
!      call time_begin(subroutine_name,idindx)
#ifdef BAND
      imax = iym2
      jmax = jx
      istart = 1
      jstart = 0
#else
      imax = iym2
      jmax = jxm2
      istart = 1
      jstart = 1
#endif

      do i = 1 , imax
        do j = 1 , jmax
#ifdef MPP1
          radpsa_io(j,i) = (psa_io(i+istart,j+jstart)+r8pt)*10.0D0
#else
          radpsa(j,i) = (sps1%ps(i+istart,j+jstart)+r8pt)*10.0D0
#endif
        end do
      end do

#ifdef MPP1
      call writerec_rad(jmax, imax, kz, 4, 9, &
                        frad3d_io(:,:,:,2:5), frad2d_io(:,:,1:10), &
                        radpsa_io, idatex)
#else
      call writerec_rad(jmax, imax, kz, 4, 9, &
                        frad3d(:,:,:,2:5), frad2d(:,:,1:10), &
                        radpsa, idatex)
#endif

      print * , 'RAD variables written at ' , idatex , xtime
 
!      call time_end(subroutine_name,idindx)
      end subroutine outrad
!
      subroutine outche

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!     this subroutine writes the model chem                           c
!                                                                     c
!     iutl : is the output unit number for large-domain variables.    c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none
!
      integer :: ni , itr , nj , nk , is , ie , js , je

!      character (len=50) :: subroutine_name='outche'
!      integer :: idindx=0
!
!      call time_begin(subroutine_name,idindx)
#ifdef BAND
      ni = iym2
      nj = jx
      nk = kz
      itr = ntr
      is = 2
      js = 1
      ie = iym1
      je = jx
#else
      ni = iym2
      nj = jxm2
      nk = kz
      itr = ntr
      is = 2
      js = 2
      ie = iym1
      je = jxm1
#endif

#ifdef MPP1
     call writerec_che(nj, ni, je, ie, nk, itr, chia_io,     &
                aerext_io, aerssa_io, aerasp_io, dtrace_io,  &
                wdlsc_io, wdcvc_io, ddsfc_io, wxsg_io,       &
                wxaq_io, cemtrac_io, aertarf_io, aersrrf_io, &
                aertalwrf_io, aersrlwrf_io, psa_io, idatex)
#else
     call writerec_che(nj, ni, je , ie , nk, itr, chia, aerext,   &
                aerssa, aerasp, dtrace, wdlsc, wdcvc, ddsfc,      &
                wxsg, wxaq, cemtrac, aertarf, aersrrf, aertalwrf, &
                aersrlwrf, sps1%ps, idatex)
#endif
      write (*,*) 'CHE variables written at ' , idatex , xtime

      end subroutine outche
!
      end module mod_output
