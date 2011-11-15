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

  use mod_runparams
  use mod_mpmessage
  use mod_service
  use mod_atm_interface
  use mod_che_interface
  use mod_lm_interface
  use mod_rad_interface
  use mod_cu_interface
  use mod_pbl_interface
  use mod_ncio
  use mod_bdycod
  use mod_precip
  use mod_split
  use mod_savefile
  use mod_mppio
#ifdef CLM
  use mod_clm
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
  use mpi
  implicit none
!
  integer :: i , j
  integer :: allrec , ierr , l , k , n
  logical :: ldoatm , ldosrf , ldorad , ldoche , ldosav , ldotmp
  logical :: lstartup
  character (len=64) :: subroutine_name='output'
  integer :: idindx=0
!
!
  call time_begin(subroutine_name,idindx)
!
!----------------------------------------------------------------------
!
  lstartup = .false.
  if ( myid == 0 ) then
    if ( ktau == 0 .or. doing_restart ) then
      call mkfile
      lstartup = .true.
    end if
  end if
!
  ldoatm = .false.
  ldosrf = .false.
  ldorad = .false.
  ldoche = .false.
  ldosav = .false.
  ldotmp = .false.

  if ( mod(ktau,ksav) == 0 .and. ktau > 0 ) then
    ldotmp = .true.
  end if
  if ( ktau > 0 .and. ( idatex == idate2 .or. &
      (lfdomonth(idatex) .and. lmidnight(idatex))) ) then
    ldosav = .true.
    ldotmp = .false.
  end if
  if ( ktau == 0 .or. mod(ktau,katm) == 0 ) then
    ldoatm = .true.
  end if
  if ( ktau == 0 .or. mod(ktau,ksrf) == 0 ) then
    ldosrf = .true.
  end if
  if ( ktau == 0 .or. mod(ktau,krad) == 0 ) then
    ldorad = .true.
  end if
  if ( ktau == 0 .or. mod(ktau,kche) == 0 ) then
    ldoche = .true.
  end if

  if ( doing_restart ) then
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
  if ( ktau == 0 ) then
    ldosrf = .false.
    ldorad = .false.
    ldoche = .false.
    lskipsrf = .true.
    lskiprad = .true.
    lskipche = .true.
  end if
!
!-----output for dataflow analyses:
!
  if ( ifatm ) then
    if ( ldoatm ) then
!=======================================================================
!     gather  ua,va,ta,qva,qca,rainc,rainnc,tgb2d,swt2d,olcd2d,rno2d
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
          atm0(i,1+kz*6,j) = sps1%ps(j,i)
          atm0(i,2+kz*6,j) = sfsta%rainc(j,i)
          atm0(i,3+kz*6,j) = sfsta%rainnc(j,i)
        end do
      end do
      do j = 1 , jendx
        do n = 1 , nnsg
          do i = 1 , iym1
            atm0(i,3+kz*6+n,j)        = tgbrd(n,j,i)
            atm0(i,3+kz*6+n+nnsg,j)   = tsw(n,j,i)
            atm0(i,3+kz*6+n+nnsg*2,j) = runoff(n,j,i)
          end do
        end do
      end do
      call mpi_gather(atm0, iy*(kz*6+3+nnsg*3)*jxp,mpi_real8,&
                    & atm_0,iy*(kz*6+3+nnsg*3)*jxp,mpi_real8,&
                    & 0,mycomm,ierr)
      if ( myid == 0 ) then
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
            psa_io(i,j)    = atm_0(i,1+kz*6,j)
            rainc_io(i,j)  = atm_0(i,2+kz*6,j)
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
              tgb2d_io(n,i,j) = atm_0(i,3+kz*6+n,j)
              swt2d_io(n,i,j) = atm_0(i,3+kz*6+n+nnsg,j)
              rno2d_io(n,i,j) = atm_0(i,3+kz*6+n+nnsg*2,j)
            end do
          end do
        end do
      end if
!=======================================================================
!     gather UW Scheme variables
      if ( ibltyp == 2 .or. ibltyp == 99 ) then
        do j = 1 , jendl
          do k = 1 , kz
            do i = 1 , iy
              uw0(i,k,j)      = atm1%tke(i,k,j)
              uw0(i,k+kz,j)   = uwstateb%kth(j,i,k)
              uw0(i,k+kz*2,j) = uwstateb%kzm(j,i,k)
            end do
          end do
        end do
        call mpi_gather(uw0, iy*(kz*3)*jxp,mpi_real8,&
                      & uw_0,iy*(kz*3)*jxp,mpi_real8,&
                      & 0,mycomm,ierr)
        if ( myid == 0 ) then
          do j = 1 , jx
            do k = 1 , kz
              do i = 1 , iy
                atm1_io%tke(i,k,j)     = uw_0(i,k,j)
                tcmstate_io%kth(i,k,j) = uw_0(i,k+kz,j)
                tcmstate_io%kzm(i,k,j) = uw_0(i,k+kz*2,j)
              end do
            end do
          end do
        end if
      end if
!
      do j = 1 , jendx
        do n = 1 , nnsg
          do i = 1 , iym1
            var2d1(i,n,j) = ocld2d(n,j,i)
          end do
        end do
      end do
      call mpi_gather(var2d1, iy*nnsg*jxp,mpi_integer, &
                    & var2d_1,iy*nnsg*jxp,mpi_integer, &
                    & 0,mycomm,ierr)
      if (myid == 0) then
#ifdef BAND
        do j = 1 , jx
#else
        do j = 1 , jxm1
#endif
          do n = 1 , nnsg
            do i = 1 , iym1
              ocld2d_io(n,i,j) = var2d_1(i,n,j)
            end do
          end do
        end do
      end if
!
      if ( myid == 0 ) then
        call outatm
      end if
      do j = 1 , jendx
        do i = 1 , iym1
          do n = 1 , nnsg
            runoff(n,j,i) = d_zero
          end do
          sfsta%rainc(j,i)  = d_zero
          sfsta%rainnc(j,i) = d_zero
        end do
      end do
    end if
  end if
 
!     Call surface output
 
  if ( ifsrf ) then
    if ( ldosrf ) then
#ifndef CLM
      if ( lakemod == 1 .and. iflak .and. mod(iolak,klak) == 0) then
       call lakegather
      end if
#endif
      if ( iseaice == 1 .or. lakemod == 1 ) then
        do j = 1 , jendx
          do i = 1 , iym1
            var2d0(i,j) = ldmsk(j,i)
          end do
        end do
        call mpi_gather(var2d0, iy*jxp,mpi_integer, &
                      & var2d_0,iy*jxp,mpi_integer, &
                      & 0,mycomm,ierr)
        if (myid == 0) then
#ifdef BAND
          do j = 1 , jx
#else
          do j = 1 , jxm1
#endif
            do i = 1 , iym1
              ldmsk_io(i,j) = var2d_0(i,j)
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
                    & 0,mycomm,ierr)
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
                      & 0,mycomm,ierr)

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
      call mpi_gather(rad0, iym2*(nrad3d*kz+nrad2d)*jxp,mpi_real4, &
                      rad_0,iym2*(nrad3d*kz+nrad2d)*jxp,mpi_real4, &
                      0,mycomm,ierr)
      swapv = transpose(sps1%ps(1:jxp,:))
      call mpi_gather(swapv, iy*jxp,mpi_real8, &
                      psa_io,iy*jxp,mpi_real8, &
                      0,mycomm,ierr)
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
          chem0(i,(ntr+3)*kz+ntr*7+5,j) = sps1%ps(j,i)
        end do
      end do
      call mpi_gather(chem0,iy*((ntr+3)*kz+ntr*7+5)*jxp,            &
                    & mpi_real8,chem_0,iy*((ntr+3)*kz+ntr*7+5)*jxp, &
                    & mpi_real8,0,mycomm,ierr)
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
        remlsc_io  = d_zero
        remcvc_io  = d_zero
        rxsg_io    = d_zero
        rxsaq1_io  = d_zero
        rxsaq2_io  = d_zero
        cemtr_io   = d_zero
        remdrd_io  = d_zero
        wdlsc_io   = d_zero
        wdcvc_io   = d_zero
        ddsfc_io   = d_zero
        wxsg_io    = d_zero
        wxaq_io    = d_zero
        cemtrac_io = d_zero
        aertarf_io = d_zero
        aersrrf_io = d_zero
        aersrlwrf_io=d_zero
        aertalwrf_io=d_zero
      end if
      do n = 1 , ntr
        do j = 1 , jendl
          do k = 1 , kz
            do i = 1 , iy
              remlsc(i,k,j,n) = d_zero
              remcvc(i,k,j,n) = d_zero
              rxsg(i,k,j,n) = d_zero
              rxsaq1(i,k,j,n) = d_zero
              rxsaq2(i,k,j,n) = d_zero
            end do
          end do
        end do
      end do
      do n = 1 , ntr
        do j = 1 , jendl
          do i = 1 , iy
            cemtr(i,j,n) = d_zero
            remdrd(i,j,n) = d_zero
            wdlsc(i,j,n) = d_zero
            wdcvc(i,j,n) = d_zero
            ddsfc(i,j,n) = d_zero
            wxsg(i,j,n) = d_zero
            wxaq(i,j,n) = d_zero
            cemtrac(i,j,n) = d_zero
          end do
        end do
      end do
      do j = 1 , jendl
        do i = 1 , iym1
          aertarf(i,j) = d_zero
          aersrrf(i,j) = d_zero
          aertalwrf(i,j) = d_zero              
          aersrlwrf(i,j) = d_zero
        end do
      end do
    end if
  end if
!
!-----output for restart:
!
  if ( ifsave ) then
    if ( ldosav .or. ldotmp ) then
#ifndef CLM
      if ( lakemod == 1 ) call lakegather
#endif
      do j = 1 , jendl
        do k = 1 , kz
          do i = 1 , iy
            sav0(i,k,j) = xub%b0(i,k,j)
            sav0(i,kz+k,j) = xvb%b0(i,k,j)
            sav0(i,kz*2+k,j) = xqb%b0(i,k,j)
            sav0(i,kz*3+k,j) = xtb%b0(i,k,j)
          end do
        end do
        do i = 1 , iy
          sav0(i,kz*4+1,j) = xpsb%b0(i,j)
          sav0(i,kz*4+2,j) = ts0(i,j)
        end do
      end do
      allrec = kz*4 + 2
      call mpi_gather(sav0, iy*allrec*jxp,mpi_real8,         &
                    & sav_0,iy*allrec*jxp,mpi_real8,         &
                    & 0,mycomm,ierr)
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
                      & 0,mycomm,ierr)
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
          sav0(i,kz*4+1,j) = sps1%ps(j,i)
          sav0(i,kz*4+2,j) = sps2%ps(j,i)
        end do
      end do
      allrec = kz*4 + 2
      call mpi_gather(sav0, iy*allrec*jxp,mpi_real8,         &
                    & sav_0,iy*allrec*jxp,mpi_real8,         &
                    & 0,mycomm,ierr)
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
          sav0(i,kz*4+1,j) = sts1%tg(j,i)
          sav0(i,kz*4+2,j) = sts2%tg(j,i)
        end do
      end do
      allrec = kz*4 + 2
      call mpi_gather(sav0, iy*allrec*jxp,mpi_real8,         &
                    & sav_0,iy*allrec*jxp,mpi_real8,         &
                    & 0,mycomm,ierr)
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
          sav0(i,kz*4+1,j) = sfsta%rainc(j,i)
          sav0(i,kz*4+2,j) = sfsta%rainnc(j,i)
        end do
      end do
      do j = 1 , jendx
        do k = 1 , kz
          do i = 1 , iym1
            sav0(i,kz*3+k,j) = heatrt(j,i,k)
          end do
        end do
      end do
      allrec = kz*4 + 2
      call mpi_gather(sav0, iy*allrec*jxp,mpi_real8,         &
                    & sav_0,iy*allrec*jxp,mpi_real8,         &
                    & 0,mycomm,ierr)
      if ( myid == 0 ) then
        do j = 1 , jx
          do k = 1 , kz
            do i = 1 , iy
              atm1_io%qc(i,k,j) = sav_0(i,k,j)
              atm2_io%qc(i,k,j) = sav_0(i,kz+k,j)
              fcc_io(i,k,j) = sav_0(i,kz*2+k,j)
            end do
          end do
          do i = 1 , iy
            rainc_io(i,j)  = sav_0(i,kz*4+1,j)
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

!!!!!!!!!!!!!!!!!!!!
!Begin UW variable gather
      if ( ibltyp == 2 .or. ibltyp == 99 ) then
        do j = 1 , jendx
          do k = 1 , kzp1
            do i = 1 , iy
              sav0b(i,k,j) = atm1%tke(i,k,j)
            end do
          end do
        end do
        allrec = kzp1
        call mpi_gather(sav0b, iy*allrec*jxp,mpi_real8,    &
                        sav_0b,iy*allrec*jxp,mpi_real8,    &
                        0,mycomm,ierr)
        if ( myid == 0 ) then
#ifdef BAND
          do j = 1 , jx
#else
          do j = 1 , jxm1
#endif
            do k = 1 , kzp1
              do i = 1 , iy
                atm1_io%tke(i,k,j) = sav_0b(i,k,j)
              end do
            end do
          end do
        end if
        do j = 1 , jendx
          do k = 1 , kzp1
            do i = 1 , iy
              sav0b(i,k,j) = atm2%tke(i,k,j)
            end do
          end do
        end do
        allrec = kzp1
        call mpi_gather(sav0b, iy*allrec*jxp,mpi_real8,    &
                        sav_0b,iy*allrec*jxp,mpi_real8,    &
                        0,mycomm,ierr)
        if ( myid == 0 ) then
#ifdef BAND
          do j = 1 , jx
#else
          do j = 1 , jxm1
#endif
            do k = 1 , kzp1
              do i = 1 , iy
                atm2_io%tke(i,k,j) = sav_0b(i,k,j)
              end do
            end do
          end do
        end if
        do j = 1 , jendx
          do i = 1 , iym1
            var2d0(i,j) = kpbl(j,i)
          end do
        end do
        call mpi_gather(var2d0, iy*jxp,mpi_integer, &
                      & var2d_0,iy*jxp,mpi_integer, &
                      & 0,mycomm,ierr)
        if (myid == 0) then
#ifdef BAND
          do j = 1 , jx
#else
          do j = 1 , jxm1
#endif
            do i = 1 , iym1
              kpbl_io(i,j) = var2d_0(i,j)
            end do
          end do
        end if
      end if
!end UW variable gather
!!!!!!!!!!!!!!!!!!!!

      do j = 1 , jendl
        do i = 1 , iy
          sav0a(i,1,j) = sfsta%hfx(j,i)
          sav0a(i,2,j) = sfsta%qfx(j,i)
          sav0a(i,3,j) = sfsta%uvdrag(j,i)
          sav0a(i,4,j) = sfsta%tgbb(j,i)
        end do
      end do
      do j = 1 , jendx
        do k = 1 , kzp1
          do i = 1 , iym1
            sav0a(i,4+k,j) = o3prof(j,i,k)
          end do
        end do
      end do
      allrec = 4 + kzp1
      call mpi_gather(sav0a, iy*allrec*jxp,mpi_real8,        &
                    & sav_0a,iy*allrec*jxp,mpi_real8,        &
                    & 0,mycomm,ierr)
      if ( myid == 0 ) then
        do j = 1 , jx
          do i = 1 , iy
            hfx_io(i,j) = sav_0a(i,1,j)
            qfx_io(i,j) = sav_0a(i,2,j)
            uvdrag_io(i,j) = sav_0a(i,3,j)
            tgbb_io(i,j) = sav_0a(i,4,j)
          end do
        end do
#ifdef BAND
        do j = 1 , jx
#else
        do j = 1 , jxm1
#endif
          do k = 1 , kzp1
            do i = 1 , iym1
              o3prof_io(i,k,j) = sav_0a(i,4+k,j)
            end do
          end do
        end do
      end if
      if ( iocnflx == 2 ) then
        swapv = transpose(zpbl)
        call mpi_gather(swapv,  iy*jxp,mpi_real8, &
                        zpbl_io,iy*jxp,mpi_real8, &
                        0,mycomm,ierr)
      end if
      if ( icup == 1 ) then
        do j = 1 , jendl
          do k = 1 , kz
            do i = 1 , iy
              sav0c(i,k,j) = rsheat(j,i,k)
              sav0c(i,kz+k,j) = rswat(j,i,k)
            end do
          end do
        end do
        allrec = kz*2
        call mpi_gather(sav0c, iy*allrec*jxp,mpi_real8,      &
                      & sav_0c,iy*allrec*jxp,mpi_real8,      &
                      & 0,mycomm,ierr)
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
                      & 0,mycomm,ierr)
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
                      & 0,mycomm,ierr)
      end if
      do j = 1 , jendx
        do l = 1 , 4
          do k = 1 , kz
            do i = 1 , iym1
              sav1(i,(l-1)*kz+k,j) = absnxt(j,i,k,l)
            end do
          end do
        end do
      end do
      allrec = kz*4
      do j = 1 , jendx
        do l = 1 , kzp1
          do k = 1 , kzp1
            do i = 1 , iym1
              sav1(i,allrec+(l-1)*(kzp1)+k,j) = abstot(j,i,k,l)
            end do
          end do
        end do
      end do
      allrec = allrec + (kzp1)*(kz+1)
      do j = 1 , jendx
        do k = 1 , kzp1
          do i = 1 , iym1
            sav1(i,allrec+k,j) = emstot(j,i,k)
          end do
        end do
      end do
      allrec = kz*4+(kzp1*kzp2)
      call mpi_gather(sav1, iym1*allrec*jxp,mpi_real8,       &
                    & sav_1,iym1*allrec*jxp,mpi_real8,       &
                    & 0,mycomm,ierr)
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
                abstot_io(i,k,l,j) = sav_1(i,allrec+(l-1)*(kzp1)+k,j)
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
            sav2(i,n,j) = taf(n,j,i)
            sav2(i,nnsg+n,j) = tlef(n,j,i)
            sav2(i,nnsg*2+n,j) = ssw(n,j,i)
            sav2(i,nnsg*3+n,j) = rsw(n,j,i)
          end do
        end do
        do i = 1 , iym1
          sav2(i,nnsg*5+1,j) = solis(j,i)
          sav2(i,nnsg*5+2,j) = solvd(j,i)
          sav2(i,nnsg*5+3,j) = solvs(j,i)
          sav2(i,nnsg*5+4,j) = flw2d(j,i)
        end do
      end do
      allrec = nnsg*5 + 4
      call mpi_gather(sav2, iym1*allrec*jxp,mpi_real8,       &
                    & sav_2,iym1*allrec*jxp,mpi_real8,       &
                    & 0,mycomm,ierr)
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
            end do
          end do
          do i = 1 , iym1
            solis_io(i,j) = sav_2(i,nnsg*5+1,j)
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
                    & 0,mycomm,ierr)
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
      call mpi_gather(lndcat2d,   iy*jxp,mpi_real8, &
                    & lndcat2d_io,iy*jxp,mpi_real8, &
                    & 0,mycomm,ierr)
#endif
      do j = 1 , jendx
        do n = 1 , nnsg
          do i = 1 , iym1
            sav2(i,n,j) = tgbrd(n,j,i)
            sav2(i,nnsg+n,j) = tsw(n,j,i)
            sav2(i,nnsg*2+n,j) = sncv(n,j,i)
            sav2(i,nnsg*3+n,j) = gwet(n,j,i)
            sav2(i,nnsg*4+n,j) = tgrd(n,j,i)
          end do
        end do
        do i = 1 , iym1
          sav2(i,nnsg*5+1,j) = flwd2d(j,i)
          sav2(i,nnsg*5+2,j) = fsw2d(j,i)
          sav2(i,nnsg*5+3,j) = sabveg(j,i)
          sav2(i,nnsg*5+4,j) = sinc(j,i)
        end do
      end do
      allrec = nnsg*5 + 4
      call mpi_gather(sav2, iym1*allrec*jxp,mpi_real8,       &
                    & sav_2,iym1*allrec*jxp,mpi_real8,       &
                    & 0,mycomm,ierr)
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
            sabveg_io(i,j) = sav_2(i,nnsg*5+3,j)
            sinc2d_io(i,j) = sav_2(i,nnsg*5+4,j)
          end do
        end do
      end if
      do j = 1 , jendx
        do n = 1 , nnsg
          do i = 1 , iym1
            sav2(i,n,j)        = ircp(n,j,i)
            sav2(i,nnsg+n,j)   = snag(n,j,i)
            sav2(i,nnsg*2+n,j) = sfice(n,j,i)
            sav2(i,nnsg*3+n,j) = dew2d(n,j,i)
            sav2(i,nnsg*4+n,j) = emiss(n,j,i)
          end do
        end do
        do i = 1 , iym1
          sav2(i,nnsg*5+1,j) = pptnc(j,i)
          sav2(i,nnsg*5+2,j) = pptc(j,i)
          sav2(i,nnsg*5+3,j) = prca2d(j,i)
          sav2(i,nnsg*5+4,j) = prnca2d(j,i)
        end do
      end do
      allrec = nnsg*5 + 4
      call mpi_gather(sav2, iym1*allrec*jxp,mpi_real8,       &
                    & sav_2,iym1*allrec*jxp,mpi_real8,       &
                    & 0,mycomm,ierr)
      if ( myid == 0 ) then
#ifdef BAND
        do j = 1 , jx
#else
        do j = 1 , jxm1
#endif
          do n = 1 , nnsg
            do i = 1 , iym1
              ircp2d_io(n,i,j) = sav_2(i,n,j)
              sag2d_io(n,i,j) = sav_2(i,nnsg+n,j)
              sice2d_io(n,i,j) = sav_2(i,nnsg*2+n,j)
              dew2d_io(n,i,j) = sav_2(i,nnsg*3+n,j)
              emiss2d_io(n,i,j) = sav_2(i,nnsg*4+n,j)
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
            sav2a(i,n,j)      = veg2d1(n,j,i)
            sav2a(i,nnsg+n,j) = ocld2d(n,j,i)
          end do
        end do
        do i = 1 , iym1
          sav2a(i,nnsg*2+1,j) = veg2d(j,i)
          sav2a(i,nnsg*2+2,j) = ldmsk(j,i)
        end do
      end do
      allrec = nnsg*2 + 2
      call mpi_gather(sav2a, iym1*allrec*jxp,mpi_integer, &
                    & sav_2a,iym1*allrec*jxp,mpi_integer, &
                    & 0,mycomm,ierr)
      if ( myid == 0 ) then
#ifdef BAND
        do j = 1 , jx
#else
        do j = 1 , jxm1
#endif
          do n = 1 , nnsg
            do i = 1 , iym1
              veg2d1_io(n,i,j) = sav_2a(i,n,j)
              ocld2d_io(n,i,j) = sav_2a(i,nnsg+n,j)
            end do
          end do
          do i = 1 , iym1
            veg2d_io(i,j) = sav_2a(i,nnsg*2+1,j)
            ldmsk_io(i,j) = sav_2a(i,nnsg*2+2,j)
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
                      & 0,mycomm,ierr)
        if ( myid == 0 ) then
          do j = 1 , jx
            do n = 1 , ntr
              do k = 1 , kz
                do i = 1 , iy
                  chia_io(i,k,j,n)   = sav_4(i,(n-1)*kz+k,j)
                  chib_io(i,k,j,n)   = sav_4(i,ntr*kz+(n-1)*kz+k,j)
                  remlsc_io(i,k,j,n) = sav_4(i,ntr*kz*2+(n-1)*kz+k,j)
                  remcvc_io(i,k,j,n) = sav_4(i,ntr*kz*3+(n-1)*kz+k,j)
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
            sav4a(i,1,j) = ssw2da(j,i)
            sav4a(i,2,j) = sdeltk2d(j,i)
            sav4a(i,3,j) = sdelqk2d(j,i)
            sav4a(i,4,j) = sfracv2d(j,i)
            sav4a(i,5,j) = sfracb2d(j,i)
            sav4a(i,6,j) = sfracs2d(j,i)
            sav4a(i,7,j) = svegfrac2d(j,i)
          end do
        end do
        call mpi_gather(sav4a, iym1*7*jxp,mpi_real8,                &
                      & sav_4a,iym1*7*jxp,mpi_real8,                &
                      & 0,mycomm,ierr)
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
            sav0d(i,n,j)        = dstor(i,j,n)
            sav0d(i,n+nsplit,j) = hstor(i,j,n)
          end do
        end do
      end do
      call mpi_gather(sav0d, iy*nsplit*2*jxp,mpi_real8,      &
                    & sav_0d,iy*nsplit*2*jxp,mpi_real8,      &
                    & 0,mycomm,ierr)
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
                    & 0,mycomm,ierr)
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
                   & mycomm,ierr)
      call mpi_bcast(ujl(1,1),iy*kz,mpi_real8,nproc-1,              &
                   & mycomm,ierr)
      call mpi_bcast(vjlx(1,1),iy*kz,mpi_real8,nproc-1,             &
                   & mycomm,ierr)
      call mpi_bcast(vjl(1,1),iy*kz,mpi_real8,nproc-1,              &
                   & mycomm,ierr)
#endif
      if ( ldosav ) then
        call write_savefile(idatex, .false.)
      else
        call write_savefile(idatex, .true.)
      end if
    end if
  end if

  if ( myid == 0 ) then
    if ( lfdomonth(idatex) .and. lmidnight(idatex) ) then
      if ( .not. lstartup .and. idatex /= idate2 ) then
        call mkfile
      end if
    end if
  end if
!
  call time_end(subroutine_name,idindx) 

  end subroutine output
!
  subroutine mkfile
 
  implicit none
!
  if (myid /= 0) return

  print * , ' '
  print * , '******* OPENING NEW OUTPUT FILES : ' , tochar(idatex)
  print * , ' '

  if ( ifatm ) then
    call prepare_common_out(idatex,'ATM')
  end if
 
  if ( ifsrf ) then
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

  call writerec_atm(jx,iy,jjx,iiy,kz,nnsg,atm1_io%u,atm1_io%v,  &
          omega_io,atm1_io%t,atm1_io%qv,atm1_io%qc,atm1_io%tke, &
          tcmstate_io%kth,tcmstate_io%kzm,psa_io,rainc_io,      &
          rainnc_io,tgb2d_io,swt2d_io,rno2d_io,ocld2d_io,idatex)
 
  print *, 'ATM variables written at ' , tochar(idatex)
 
  end subroutine outatm
!
  subroutine outsrf

  implicit none
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

  call writerec_srf(j,i,numbat,fbat_io,ldmsk_io,idatex)
  print *, 'SRF variables written at ' , tochar(idatex)
 
#ifndef CLM
  if (lakemod == 1 .and. iflak .and. mod(iolak,klak) == 0) then
    call writerec_lak(j,i,numbat,fbat_io,evl2d_io,aveice2d_io, &
                      hsnow2d_io,tlak3d_io,idatex)
    print *, 'LAK variables written at ' , tochar(idatex)
  end if
#endif

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

  call writerec_sub(j,i,nsg,numsub,fsub_io,idatex)

  print *, 'SUB variables written at ' , tochar(idatex)

  end subroutine outsub
!
  subroutine outrad
!
  implicit none
!
  integer :: i , j , imax , jmax , istart, jstart
!
!      character (len=64) :: subroutine_name='outrad'
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
      radpsa_io(j,i) = real((psa_io(i+istart,j+jstart)+ptop)*d_10)
    end do
  end do

  call writerec_rad(jmax, imax, kz, 4, 9, &
                    frad3d_io(:,:,:,2:5), frad2d_io(:,:,1:10), &
                    radpsa_io, idatex)

  print * , 'RAD variables written at ' , tochar(idatex)
 
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
  integer :: ni , itr , nj , nk , ie , je

!      character (len=64) :: subroutine_name='outche'
!      integer :: idindx=0
!
!      call time_begin(subroutine_name,idindx)
#ifdef BAND
  ni = iym2
  nj = jx
  nk = kz
  itr = ntr
  ie = iym1
  je = jx
#else
  ni = iym2
  nj = jxm2
  nk = kz
  itr = ntr
  ie = iym1
  je = jxm1
#endif

     call writerec_che(nj, ni, je, ie, nk, itr, chia_io,     &
            aerext_io, aerssa_io, aerasp_io, dtrace_io,  &
            wdlsc_io, wdcvc_io, ddsfc_io, wxsg_io,       &
            wxaq_io, cemtrac_io, aertarf_io, aersrrf_io, &
            aertalwrf_io, aersrlwrf_io, psa_io, idatex)

  print *, 'CHE variables written at ' , tochar(idatex)

  end subroutine outche
!
end module mod_output
