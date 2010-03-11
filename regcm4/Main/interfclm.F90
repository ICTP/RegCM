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

#ifdef MPP1

      subroutine interfclm(inout,nstep)

!=======================================================================
!l  built for clm version 3.0
!=======================================================================
! inout = 1 : regcm -> clm
! inout = 2 : clm -> regcm
!
      use mod_regcm_param
      use clm_varsur,    only : landmask, landfrac
      use clmtype
      use clm_varsur,    only : c2r_allout,omap_i,omap_j
      use mpi
      use mod_clm
      use mod_date
      use mod_main
      use mod_param1
      use mod_param2
      use mod_param3
      use mod_slice
      use mod_pbldim
      use mod_bats
      use mod_constants
      use mod_interfaces
      implicit none
!
! Dummy arguments
!
      integer :: inout , nstep
      intent (in) inout
!
! Local variables
!
      real(8) :: amxtem , facb , facs , fact , factuv , facv , fracb ,  &
                    & fracs , fracv , mmpd , sfac , solvt , wpm2
      integer :: ci , cj , counter , i , icount , ii , iii , j , je ,   &
               & ji , jj , jo , js , kk , locid , n , nn1 , nnn , nout ,&
               & p , x
      real(8) , dimension(jxp,ix) :: r2cflwd , r2cpsb , r2cqb ,    &
                & r2crnc , r2crnnc , r2csoll , r2csolld , r2csols ,     &
                & r2csolsd , r2ctb , r2cuxb , r2cvxb , r2czga
      real(4) :: real_4
      real(8) , dimension(jxp*ix*13) :: workin
      real(8) , dimension(jx*ix*13) :: workout
      integer :: ierr
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!     About the dimension ordering:
!     regcm: ix=lat,jx=lon, arrays are lat by lon
!     clm: i=lon, j=lat, arrays are lon by lat
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
      if ( inout==1 ) then
 
        locid = 1
!       clm3 currently works on all ix,jx instead of 2-ilx and 2-jlx so
!       copy neighboring values for now
        do j = 1 , jxp
          do i = 1 , ix
 
!jlb        10/05 treat all variables identically. Copy i=2 to i=1 and
!           i=ilx to i=ix. For myid = 0, copy j=2 to j=1. For myid =
!           nproc-1, copy j=jendx to j=jxp.
            if ( myid==0 .and. j==1 ) then
              cj = 2
            else if ( myid==(nproc-1) .and. j==jxp ) then
              cj = jxp - 1
            else
              cj = j
            end if
            if ( i==1 ) then
              ci = 2
            else if ( i==ix ) then
              ci = ix - 1
            else
              ci = i
            end if
 
!           T(K) at bottom layer
            r2ctb(j,i) = tb3d(ci,kx,cj)
!           Specific Humidity ?
            r2cqb(j,i) = qvb3d(ci,kx,cj)/(1+qvb3d(ci,kx,cj))
!           Reference Height (m)
            r2czga(j,i) = za(ci,kx,cj)
!           Surface winds
            r2cuxb(j,i) = ubx3d(ci,kx,cj)
            r2cvxb(j,i) = vbx3d(ci,kx,cj)
!           Surface Pressure in Pa from hPa
            r2cpsb(j,i) = (psb(ci,cj)+ptop)*1000.
!           Rainfall
            r2crnc(j,i) = pptc(ci,cj)
            r2crnnc(j,i) = pptnc(ci,cj)
!           Incident Solar Radiation
            r2csols(j,i) = sols2d(ci,cj)
            r2csoll(j,i) = soll2d(ci,cj)
            r2csolsd(j,i) = solsd2d(ci,cj)
            r2csolld(j,i) = solld2d(ci,cj)
            r2cflwd(j,i) = flwd2d(ci,cj)
 
          end do
        end do
 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c      cc
!c      1. Copy 2d (jxp,ix) arrays to 1d work_in (jx*ix) array.      cc
!c      2. Gather jxp values of each nproc work_in array and fill     cc
!c      work_out(jx*ix) array.                                   cc
!c      3. Copy 1d work_out array to 2d (jx,ix) array for passing    cc
!c      to clm.                                                   cc
!c      abt updated below 1/09                                        cc
!c      UPDATE:  copy all r2c vars to one large array; this allows    cc
!c      for one MPI_ALLGATHER call instead of several        cc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
        ii = 1
        do j = 1 , jxp
          do i = 1 , ix
            workin(ii) = r2ctb(j,i)
            workin(ii+(jxp*ix)) = r2cqb(j,i)
            workin(ii+(2*jxp*ix)) = r2czga(j,i)
            workin(ii+(3*jxp*ix)) = r2cpsb(j,i)
            workin(ii+(4*jxp*ix)) = r2cuxb(j,i)
            workin(ii+(5*jxp*ix)) = r2cvxb(j,i)
            workin(ii+(6*jxp*ix)) = r2crnc(j,i)
            workin(ii+(7*jxp*ix)) = r2crnnc(j,i)
            workin(ii+(8*jxp*ix)) = r2csols(j,i)
            workin(ii+(9*jxp*ix)) = r2csoll(j,i)
            workin(ii+(10*jxp*ix)) = r2csolsd(j,i)
            workin(ii+(11*jxp*ix)) = r2csolld(j,i)
            workin(ii+(12*jxp*ix)) = r2cflwd(j,i)
            ii = ii + 1
          end do
        end do
        call mpi_allgather(workin,13*jxp*ix,mpi_double_precision,       &
                         & workout,13*jxp*ix,mpi_double_precision,      &
                         & mpi_comm_world,ierr)
 
        ii = 1
        kk = 1
        counter = 1
        do j = 1 , jx
          do i = 1 , ix
            r2ctb_all(j,i) = workout(ii)
            r2cqb_all(j,i) = workout(ii+(jxp*ix))
            r2czga_all(j,i) = workout(ii+(2*jxp*ix))
            r2cpsb_all(j,i) = workout(ii+(3*jxp*ix))
            r2cuxb_all(j,i) = workout(ii+(4*jxp*ix))
            r2cvxb_all(j,i) = workout(ii+(5*jxp*ix))
            r2crnc_all(j,i) = workout(ii+(6*jxp*ix))
            r2crnnc_all(j,i) = workout(ii+(7*jxp*ix))
            r2csols_all(j,i) = workout(ii+(8*jxp*ix))
            r2csoll_all(j,i) = workout(ii+(9*jxp*ix))
            r2csolsd_all(j,i) = workout(ii+(10*jxp*ix))
            r2csolld_all(j,i) = workout(ii+(11*jxp*ix))
            r2cflwd_all(j,i) = workout(ii+(12*jxp*ix))
            ii = ii + 1
            counter = counter + 1
          end do
          if ( counter>jxp*ix ) then
            kk = kk + 1
            counter = 1
            ii = jxp*(kk-1)*13*ix + 1
          end if
        end do

      else if ( inout==2 ) then ! end of inout = 1

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c
!c 1. Copy the parts of 2d (jx,ix) clm arrays that contain
!c     data to 1d work_in (jx*ix) array.
!c 2. Gather jxp values of each nproc work_in array and fill
!c     work_out (jx*ix) array.
!c 3. Copy 1d work_out array to 2d (jx,ix) array for passing
!c     to nproc 2d (jxp,ix) arrays.
!c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
!JLB    3/06: NEED TO CREATE MASK/FILTER FOR EACH CPU. CLM DOES NOT
!       ASSIGN GRID CELLS IN QUITE THE ROUND ROBIN FASHION INDICATED.
!       TGB
!       3/06 -> New routine utilizing c2rprocmap in case of non-round
!       robin cpu assignment in clm.
!abt    updated below
!       1/09 -> update use MPI_ALLGATHER call located in clm_atmlnd.F90
!       module This module saves each land surface variable regcm needs
!       into a large array (c2r_all).  gathers that array (c2r_allout)
!       and the corresponding grid point on the i/j grid to omap_i/j
!       Finally, below c2r_allout vars are placed into the proper place
!       on the grid for each c2r variable
!       NOTE: CLM refers to i = lon while REGCM refers to j = lon
!
        iii = 0
        jj = 1
        if ( aertyp/='AER00D0' ) then
          nout = 22
        else
          nout = 20
        end if
        do nn1 = 1 , nproc
          do ii = 1 , c2rngc(nn1)
            kk = c2rngc(nn1)
            j = omap_i(jj)
            i = omap_j(jj)
 
            c2rtgb(j,i) = c2r_allout(ii+iii)
            c2rsnowc(j,i) = c2r_allout(ii+kk+iii)
            c2rsenht(j,i) = c2r_allout(ii+(2*kk)+iii)
            c2rlatht(j,i) = c2r_allout(ii+(3*kk)+iii)
            c2ruvdrag(j,i) = c2r_allout(ii+(4*kk)+iii)
            c2ralbdirs(j,i) = c2r_allout(ii+(5*kk)+iii)
            c2ralbdirl(j,i) = c2r_allout(ii+(6*kk)+iii)
            c2ralbdifs(j,i) = c2r_allout(ii+(7*kk)+iii)
            c2ralbdifl(j,i) = c2r_allout(ii+(8*kk)+iii)
            c2rtgbb(j,i) = c2r_allout(ii+(9*kk)+iii)
            c2r2mt(j,i) = c2r_allout(ii+(10*kk)+iii)
            c2r2mq(j,i) = c2r_allout(ii+(11*kk)+iii)
            c2ru10(j,i) = c2r_allout(ii+(12*kk)+iii)
            c2rtlef(j,i) = c2r_allout(ii+(13*kk)+iii)
            c2rsm10cm(j,i) = c2r_allout(ii+(14*kk)+iii)
            c2rsm1m(j,i) = c2r_allout(ii+(15*kk)+iii)
            c2rsmtot(j,i) = c2r_allout(ii+(16*kk)+iii)
            c2rinfl(j,i) = c2r_allout(ii+(17*kk)+iii)
            c2rro_sur(j,i) = c2r_allout(ii+(18*kk)+iii)
            c2rro_sub(j,i) = c2r_allout(ii+(19*kk)+iii)
            if ( aertyp/='AER00D0' ) then
              c2rfracsno(j,i) = c2r_allout(ii+(20*kk)+iii)
              c2rfvegnosno(j,i) = c2r_allout(ii+(21*kk)+iii)
            end if
            jj = jj + 1
          end do
          iii = iii + c2rngc(nn1)*nout
        end do
 
        deallocate(c2r_allout)
 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!c      cc
!c      Fill nproc 2d (jxp,ix) arrays from full 2d (jx,ix) clm data. cc
!c      cc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
        if ( jyear==jyear0 .and. ktau<=1 ) then
          mmpd = 86400./dtbat
          wpm2 = 1./dtbat
        else if ( jyear==jyear0 .and. dble(ktau*dtmin)<=batfrq*60.+     &
                & 0.01 ) then
          mmpd = 24./(batfrq-dtmin/60.)
          wpm2 = 1./((batfrq-dtmin/60.)*3600.)
        else
          mmpd = 24./batfrq
          wpm2 = 1./(batfrq*3600.)
        end if
 
 
        jj = 0
        do j = jbegin , jendx
          jj = (jxp*myid) + j
 
          call interf(1 , j , kx , 2 , ixm1 , nnsg)
          if ( iocnflx==2 ) call zengocndrv(j, nnsg , 2 , ixm1 , kx)
 
          do i = 2 , ixm1
            ci = i
            uvdrag(i,j) = 0.0
            hfx(i,j) = 0.0
            qfx(i,j) = 0.0
            tgb(i,j) = 0.0
            tga(i,j) = 0.0
            tgbb(i,j) = 0.0
!chem2
            ssw2da(i,j) = 0.0
            sdeltk2d(i,j) = 0.0
            sdelqk2d(i,j) = 0.0
            sfracv2d(i,j) = 0.0
            sfracb2d(i,j) = 0.0
            sfracs2d(i,j) = 0.0
!chem2_
            if ( landmask(jj,ci)==1 ) then
              tgb(i,j) = c2rtgb(jj,ci)
              tga(i,j) = c2rtgb(jj,ci)
              hfx(i,j) = c2rsenht(jj,ci)
              qfx(i,j) = c2rlatht(jj,ci)
              uvdrag(i,j) = c2ruvdrag(jj,ci)
              tgbb(i,j) = c2rtgbb(jj,ci)
 
              if ( r2cdoalb ) coszrs2d(i,j) = c2rcosz(jj,ci)
              if ( i<=ix-1 ) then
                aldirs2d(i,j) = c2ralbdirs(jj,ci)
                aldirl2d(i,j) = c2ralbdirl(jj,ci)
                aldifs2d(i,j) = c2ralbdifs(jj,ci)
                aldifl2d(i,j) = c2ralbdifl(jj,ci)
              end if
 
              do n = 1 , nnsg
                snowc(n,i,j) = c2rsnowc(jj,ci)
                tg2d(n,i,j) = c2rtgb(jj,ci)
                tgb2d(n,i,j) = c2rtgb(jj,ci)
                !supposed to be lower soil layer temp not tgrnd
                taf2d(n,i,j) = c2r2mt(jj,ci)
                tlef2d(n,i,j) = c2rtlef(jj,ci)
                swt2d(n,i,j) = c2rsmtot(jj,ci)
                srw2d(n,i,j) = c2rsm1m(jj,ci)
                ssw2d(n,i,j) = c2rsm10cm(jj,ci)
                dew2d(n,i,j) = ldew1d(n,i)
                sag2d(n,i,j) = sag1d(n,i)    !snow age
                scv2d(n,i,j) = c2rsnowc(jj,ci)
                sice2d(n,i,j) = sice1d(n,i)  ! sea ice
                gwet2d(n,i,j) = gwet1d(n,i)
                ircp2d(n,i,j) = ircp1d(n,i)
                evpa2d(n,i,j) = evpa2d(n,i,j) + dtbat*qfx(i,j)
                sena2d(n,i,j) = sena2d(n,i,j) + dtbat*hfx(i,j)
                rnos2d(n,i,j) = c2rro_sur(jj,ci)*dtbat
                rno2d(n,i,j) = (c2rro_sub(jj,ci)+c2rro_sur(jj,ci))*dtbat
 
                ssw2da(i,j) = ssw2da(i,j) + ssw2d(n,i,j)
                sdeltk2d(i,j) = sdeltk2d(i,j) + delt1d(n,i)
                sdelqk2d(i,j) = sdelqk2d(i,j) + delq1d(n,i)
                sfracv2d(i,j) = sfracv2d(i,j) + c2rfvegnosno(jj,ci)
                sfracb2d(i,j) = sfracb2d(i,j)                           &
                              & + 1 - (c2rfvegnosno(jj,ci)+             &
                              & c2rfracsno(jj,ci))
                sfracs2d(i,j) = sfracs2d(i,j) + c2rfracsno(jj,ci)
              end do
 
              !abt added for 2m humidity when landmask = 1 or 3
              q2d(i,j) = c2r2mq(jj,ci)
!
!             quantities stored on 2d surface array for bats use only
!
              prca2d(i,j) = prca2d(i,j) + dtbat*pptc(i,j)
              prnca2d(i,j) = prnca2d(i,j) + dtbat*pptnc(i,j)
              flwa2d(i,j) = flwa2d(i,j) + dtbat*flw1d(i)
              flwda2d(i,j) = flwda2d(i,j) + dtbat*flwd2d(i,j)
              fswa2d(i,j) = fswa2d(i,j) + dtbat*fsw1d(i)
              svga2d(i,j) = svga2d(i,j) + dtbat*sabveg(i)
              sina2d(i,j) = sina2d(i,j) + dtbat*sinc2d(i,j)
              pptnc(i,j) = 0.
              pptc(i,j) = 0.
!chem2
              ssw2da(i,j) = ssw2da(i,j)/float(nnsg)
              sdeltk2d(i,j) = sdeltk2d(i,j)/float(nnsg)
              sdelqk2d(i,j) = sdelqk2d(i,j)/float(nnsg)
              sfracv2d(i,j) = sfracv2d(i,j)/float(nnsg)
              sfracb2d(i,j) = sfracb2d(i,j)/float(nnsg)
              sfracs2d(i,j) = sfracs2d(i,j)/float(nnsg)
!             svegfrac2d(i,j)= svegfrac2d(i,j)/float(NNSG)
!chem2_
            else if ( landmask(jj,ci)==0 ) then !ocean
 
              do n = 1 , nnsg
                uvdrag(i,j) = uvdrag(i,j) + drag1d(n,i)
                hfx(i,j) = hfx(i,j) + sent1d(n,i)
                qfx(i,j) = qfx(i,j) + evpr1d(n,i)
                tgb(i,j) = tgb(i,j) + tg1d(n,i)
                tga(i,j) = tga(i,j) + tg1d(n,i)
!chem2
                ssw2da(i,j) = ssw2da(i,j) + ssw1d(n,i)
                sdeltk2d(i,j) = sdeltk2d(i,j) + delt1d(n,i)
                sdelqk2d(i,j) = sdelqk2d(i,j) + delq1d(n,i)
                sfracv2d(i,j) = sfracv2d(i,j) + sigf(n,i)
                sfracb2d(i,j) = sfracb2d(i,j) + (1.-sigf(n,i))          &
                              & *(1.-scvk(n,i))
                sfracs2d(i,j) = sfracs2d(i,j) + sigf(n,i)*wt(n,i)       &
                              & + (1.-sigf(n,i))*scvk(n,i)
!               svegfrac2d(i,j)= svegfrac2d(i,j)+veg1d(n,i)
!               abt *** svegfrac2d is assumed to not change over time
!               *** therefore svegfrac2d is set in init_clmpara/ser
!               *** beware of this if attempting to run CLM with
!chem2_         Dynamic Veg
 
                if ( iocnflx==1 .or.                                    &
                   & (iocnflx==2 .and. ocld2d(n,i,j)>=0.5) ) then
                  tgbb(i,j) = tgbb(i,j)                                 &
                            & + ((1.-veg1d(n,i))*tg1d(n,i)**4+veg1d(n,i)&
                            & *tlef1d(n,i)**4)**0.25
                else
                  tgbb(i,j) = tgbb(i,j) + tg1d(n,i)
                end if
                ssw1d(n,i) = -1.E34
                rsw1d(n,i) = -1.E34
                tsw1d(n,i) = -1.E34
                rno1d(n,i) = -1.E34
                rnos1d(n,i) = -1.E34
                scv1d(n,i) = -1.E34
              end do
 
              uvdrag(i,j) = uvdrag(i,j)/float(nnsg)
              hfx(i,j) = hfx(i,j)/float(nnsg)
              qfx(i,j) = qfx(i,j)/float(nnsg)
              tgb(i,j) = tgb(i,j)/float(nnsg)
              tga(i,j) = tga(i,j)/float(nnsg)
              tgbb(i,j) = tgbb(i,j)/float(nnsg)
!chem2
              ssw2da(i,j) = ssw2da(i,j)/float(nnsg)
              sdeltk2d(i,j) = sdeltk2d(i,j)/float(nnsg)
              sdelqk2d(i,j) = sdelqk2d(i,j)/float(nnsg)
              sfracv2d(i,j) = sfracv2d(i,j)/float(nnsg)
              sfracb2d(i,j) = sfracb2d(i,j)/float(nnsg)
              sfracs2d(i,j) = sfracs2d(i,j)/float(nnsg)
!             svegfrac2d(i,j)= svegfrac2d(i,j)/float(NNSG)
!chem2_
              do n = 1 , nnsg
                snowc(n,i,j) = scv1d(n,i)
!fix
                tg2d(n,i,j) = tg1d(n,i)
!fix_
                tgb2d(n,i,j) = tgb1d(n,i)
!               taf2d(n,i,j)=taf1d(n,i)
                taf2d(n,i,j) = t2m_1d(n,i)
                !note taf2d is not temp in canopy but 2m temp
                tlef2d(n,i,j) = tlef1d(n,i)
                swt2d(n,i,j) = tsw1d(n,i)
                srw2d(n,i,j) = rsw1d(n,i)
                ssw2d(n,i,j) = ssw1d(n,i)
                dew2d(n,i,j) = ldew1d(n,i)
                sag2d(n,i,j) = sag1d(n,i)
                scv2d(n,i,j) = scv1d(n,i)
                sice2d(n,i,j) = sice1d(n,i)
                gwet2d(n,i,j) = gwet1d(n,i)
                ircp2d(n,i,j) = ircp1d(n,i)
                evpa2d(n,i,j) = evpa2d(n,i,j) + dtbat*evpr1d(n,i)
                sena2d(n,i,j) = sena2d(n,i,j) + dtbat*sent1d(n,i)
                if ( rnos2d(n,i,j)>-1.E10 .and. rnos1d(n,i)>-1.E10 )    &
                   & then
                  rnos2d(n,i,j) = rnos2d(n,i,j) + rnos1d(n,i)/tau1*dtbat
                else
                  rnos2d(n,i,j) = -1.E34
                end if
                if ( rno2d(n,i,j)>-1.E10 .and. rnos1d(n,i)>-1.E10 .and. &
                   & rno1d(n,i)>-1.E10 ) then
                  rno2d(n,i,j) = rno2d(n,i,j) + (rno1d(n,i)-rnos1d(n,i))&
                               & /tau1*dtbat
                else
                  rno2d(n,i,j) = -1.E34
                end if
              end do
!
!             quantities stored on 2d surface array for bats use only
!
              prca2d(i,j) = prca2d(i,j) + dtbat*pptc(i,j)
              prnca2d(i,j) = prnca2d(i,j) + dtbat*pptnc(i,j)
              flwa2d(i,j) = flwa2d(i,j) + dtbat*flw1d(i)
              flwda2d(i,j) = flwda2d(i,j) + dtbat*flwd2d(i,j)
              fswa2d(i,j) = fswa2d(i,j) + dtbat*fsw1d(i)
              svga2d(i,j) = svga2d(i,j) + dtbat*sabveg(i)
              sina2d(i,j) = sina2d(i,j) + dtbat*sinc2d(i,j)
              pptnc(i,j) = 0.
              pptc(i,j) = 0.
 
            else if ( landmask(jj,ci)==3 ) then
            !gridcell with some % land and ocean
 
              do n = 1 , nnsg
                uvdrag(i,j) = uvdrag(i,j) + drag1d(n,i)
                hfx(i,j) = hfx(i,j) + sent1d(n,i)
                qfx(i,j) = qfx(i,j) + evpr1d(n,i)
                tgb(i,j) = tgb(i,j) + tg1d(n,i)
                tga(i,j) = tga(i,j) + tg1d(n,i)
!chem2
                ssw2da(i,j) = ssw2da(i,j) + ssw2d(n,i,j)
                sdeltk2d(i,j) = sdeltk2d(i,j) + delt1d(n,i)
                sdelqk2d(i,j) = sdelqk2d(i,j) + delq1d(n,i)
                sfracv2d(i,j) = sfracv2d(i,j) + c2rfvegnosno(jj,ci)
                sfracb2d(i,j) = sfracb2d(i,j)                           &
                              & + 1 - (c2rfvegnosno(jj,ci)+             &
                              & c2rfracsno(jj,ci))
                sfracs2d(i,j) = sfracs2d(i,j) + c2rfracsno(jj,ci)
 
                ssw2da(i,j) = ssw2da(i,j)*landfrac(jj,ci)               &
                            & + (1-landfrac(jj,ci))*ssw1d(n,i)
                sdeltk2d(i,j) = sdeltk2d(i,j)*landfrac(jj,ci)           &
                              & + (1-landfrac(jj,ci))*delt1d(n,i)
                sdelqk2d(i,j) = sdelqk2d(i,j)*landfrac(jj,ci)           &
                              & + (1-landfrac(jj,ci))*delq1d(n,i)
                sfracv2d(i,j) = sfracv2d(i,j)*landfrac(jj,ci)           &
                              & + (1-landfrac(jj,ci))*sigf(n,i)
                sfracb2d(i,j) = sfracb2d(i,j)*landfrac(jj,ci)           &
                              & + (1-landfrac(jj,ci))*(1.-sigf(n,i))    &
                              & *(1.-scvk(n,i))
                sfracs2d(i,j) = sfracs2d(i,j)*landfrac(jj,ci)           &
                              & + (1-landfrac(jj,ci))                   &
                              & *(sigf(n,i)*wt(n,i)+(1.-sigf(n,i))      &
                              & *scvk(n,i))
!chem2_
 
                if ( iocnflx==1 .or.                                    &
                   & (iocnflx==2 .and. ocld2d(n,i,j)>=0.5) ) then
                  tgbb(i,j) = tgbb(i,j)                                 &
                            & + ((1.-veg1d(n,i))*tg1d(n,i)**4+veg1d(n,i)&
                            & *tlef1d(n,i)**4)**0.25
                else
                  tgbb(i,j) = tgbb(i,j) + tg1d(n,i)
                end if
!               ssw1d(n,i)=-1.e34
!               rsw1d(n,i)=-1.e34
!               tsw1d(n,i)=-1.e34
!               rno1d(n,i)=-1.e34
!               rnos1d(n,i)=-1.e34
!               scv1d(n,i)=-1.e34
              end do
 
              uvdrag(i,j) = uvdrag(i,j)*(1-landfrac(jj,ci))             &
                          & + c2ruvdrag(jj,ci)*landfrac(jj,ci)
              hfx(i,j) = hfx(i,j)*(1-landfrac(jj,ci)) + c2rsenht(jj,ci) &
                       & *landfrac(jj,ci)
              qfx(i,j) = qfx(i,j)*(1-landfrac(jj,ci)) + c2rlatht(jj,ci) &
                       & *landfrac(jj,ci)
              tgb(i,j) = tgb(i,j)*(1-landfrac(jj,ci)) + c2rtgb(jj,ci)   &
                       & *landfrac(jj,ci)
              tgbb(i,j) = tgbb(i,j)*(1-landfrac(jj,ci)) + c2rtgbb(jj,ci)&
                        & *landfrac(jj,ci)
              tga(i,j) = tga(i,j)*(1-landfrac(jj,ci)) + c2rtgb(jj,ci)   &
                       & *landfrac(jj,ci)
 
              uvdrag(i,j) = uvdrag(i,j)/float(nnsg)
              hfx(i,j) = hfx(i,j)/float(nnsg)
              qfx(i,j) = qfx(i,j)/float(nnsg)
              tgb(i,j) = tgb(i,j)/float(nnsg)
              tga(i,j) = tga(i,j)/float(nnsg)
              tgbb(i,j) = tgbb(i,j)/float(nnsg)
!chem2
              ssw2da(i,j) = ssw2da(i,j)/float(nnsg)
              sdeltk2d(i,j) = sdeltk2d(i,j)/float(nnsg)
              sdelqk2d(i,j) = sdelqk2d(i,j)/float(nnsg)
              sfracv2d(i,j) = sfracv2d(i,j)/float(nnsg)
              sfracb2d(i,j) = sfracb2d(i,j)/float(nnsg)
              sfracs2d(i,j) = sfracs2d(i,j)/float(nnsg)
!             svegfrac2d(i,j)= svegfrac2d(i,j)/float(NNSG)
!chem2_
              do n = 1 , nnsg
                dew2d(n,i,j) = ldew1d(n,i)
                sag2d(n,i,j) = sag1d(n,i)
                scv2d(n,i,j) = scv1d(n,i)
                sice2d(n,i,j) = sice1d(n,i)
                gwet2d(n,i,j) = gwet1d(n,i)
                ircp2d(n,i,j) = ircp1d(n,i)
!abt            added below for the landfraction method
                snowc(n,i,j) = c2rsnowc(jj,ci)*landfrac(jj,ci)          &
                             & + scv1d(n,i)*(1-landfrac(jj,ci))
                tg2d(n,i,j) = c2rtgb(jj,ci)*landfrac(jj,ci) + tg1d(n,i) &
                            & *(1-landfrac(jj,ci))
                tgb2d(n,i,j) = c2rtgb(jj,ci)*landfrac(jj,ci)            &
                             & + tgb1d(n,i)*(1-landfrac(jj,ci))
                taf2d(n,i,j) = c2r2mt(jj,ci)*landfrac(jj,ci)            &
                             & + t2m_1d(n,i)*(1-landfrac(jj,ci))
                !note taf2d is 2m temp not temp in foilage
                tlef2d(n,i,j) = c2rtlef(jj,ci)*landfrac(jj,ci)          &
                              & + tlef1d(n,i)*(1-landfrac(jj,ci))
                swt2d(n,i,j) = c2rsmtot(jj,ci)*landfrac(jj,ci)          &
                             & + tsw1d(n,i)*(1-landfrac(jj,ci))
                srw2d(n,i,j) = c2rsm1m(jj,ci)*landfrac(jj,ci)           &
                             & + rsw1d(n,i)*(1-landfrac(jj,ci))
                ssw2d(n,i,j) = c2rsm10cm(jj,ci)*landfrac(jj,ci)         &
                             & + ssw1d(n,i)*(1-landfrac(jj,ci))
                q2d(i,j) = c2r2mq(jj,ci)*landfrac(jj,ci) + q2m_1d(n,i)  &
                         & *(1-landfrac(jj,ci))
 
 
                evpa2d(n,i,j) = evpa2d(n,i,j) + dtbat*qfx(i,j)
                sena2d(n,i,j) = sena2d(n,i,j) + dtbat*hfx(i,j)
                rnos2d(n,i,j) = c2rro_sur(jj,ci)*dtbat
                rno2d(n,i,j) = c2rro_sub(jj,ci)*dtbat + c2rro_sur(jj,ci)&
                             & *dtbat
!abt            above
              end do
!
!             quantities stored on 2d surface array for bats use only
!
              prca2d(i,j) = prca2d(i,j) + dtbat*pptc(i,j)
              prnca2d(i,j) = prnca2d(i,j) + dtbat*pptnc(i,j)
              flwa2d(i,j) = flwa2d(i,j) + dtbat*flw1d(i)
              flwda2d(i,j) = flwda2d(i,j) + dtbat*flwd2d(i,j)
              fswa2d(i,j) = fswa2d(i,j) + dtbat*fsw1d(i)
              svga2d(i,j) = svga2d(i,j) + dtbat*sabveg(i)
              sina2d(i,j) = sina2d(i,j) + dtbat*sinc2d(i,j)
              pptnc(i,j) = 0.
              pptc(i,j) = 0.
 
            else
               !landmask
            end if
          end do !i loop
 
!!!!!!!!!!!!!!!! addition from new RegCM !!!!!!!!!!!!!!!!!!!
 
          do i = 2 , ixm1
            ci = i
 
            u10m_o(j,i-1) = 0.0
            v10m_o(j,i-1) = 0.0
            tg_o(j,i-1) = 0.0
            t2m_o(j,i-1) = 0.0
 
            do n = 1 , nnsg
              if ( ocld2d(n,i,j)>0.5 ) then
                u10m_s(n,j,i-1) = ubx3d(i,kx,j)
                v10m_s(n,j,i-1) = vbx3d(i,kx,j)
                tg_s(n,j,i-1) = tg2d(n,i,j)
                t2m_s(n,j,i-1) = taf2d(n,i,j)
!               abt            u10m_o(j,i-1)= u10m_o(j,i-1)+ u10m1d(n,i)
!               abt            v10m_o(j,i-1)= v10m_o(j,i-1)+ v10m1d(n,i)
                u10m_o(j,i-1) = u10m_o(j,i-1) + ubx3d(i,kx,j)
                v10m_o(j,i-1) = v10m_o(j,i-1) + vbx3d(i,kx,j)
                t2m_o(j,i-1) = t2m_o(j,i-1) + taf2d(n,i,j)
                tg_o(j,i-1) = tg_o(j,i-1) + tg2d(n,i,j)
              else if ( ocld2d(n,i,j)<0.5 ) then
                tg_s(n,j,i-1) = tg1d(n,i)
                u10m_s(n,j,i-1) = u10m1d(n,i)
                v10m_s(n,j,i-1) = v10m1d(n,i)
                t2m_s(n,j,i-1) = t2m_1d(n,i)
 
                u10m_o(j,i-1) = u10m_o(j,i-1) + u10m1d(n,i)
                v10m_o(j,i-1) = v10m_o(j,i-1) + v10m1d(n,i)
                t2m_o(j,i-1) = t2m_o(j,i-1) + t2m_1d(n,i)
                tg_o(j,i-1) = tg_o(j,i-1) + tg1d(n,i)
              else
              end if
            end do
 
            u10m_o(j,i-1) = u10m_o(j,i-1)/float(nnsg)
            v10m_o(j,i-1) = v10m_o(j,i-1)/float(nnsg)
            t2m_o(j,i-1) = t2m_o(j,i-1)/float(nnsg)
            tg_o(j,i-1) = tg_o(j,i-1)/float(nnsg)
 
            tgmx_o(j,i-1) = amax1(tgmx_o(j,i-1),tg_o(j,i-1))
            tgmn_o(j,i-1) = amin1(tgmn_o(j,i-1),tg_o(j,i-1))
            t2mx_o(j,i-1) = amax1(t2mx_o(j,i-1),t2m_o(j,i-1))
            t2mn_o(j,i-1) = amin1(t2mn_o(j,i-1),t2m_o(j,i-1))
            w10x_o(j,i-1) = amax1(w10x_o(j,i-1),sqrt(u10m_o(j,i-1)**2+  &
                          & v10m_o(j,i-1)**2))
            real_4 = (psb(i,j)+ptop)*10.
            psmn_o(j,i-1) = amin1(psmn_o(j,i-1),real_4)
 
          end do
             !i loop
 
          if ( mod(ntime+nint(dtmin*60.),kbats)==0 .or.                 &
             & (jyear==jyearr .and. ktau==ktaur) ) then
            if ( jyear==jyear0 .and. ktau<=1 ) then
              mmpd = 86400./dtbat
              wpm2 = 1./dtbat
            else if ( jyear==jyear0 .and. dble(ktau*dtmin)<=batfrq*60.+ &
                    & 0.01 ) then
              mmpd = 24./(batfrq-dtmin/60.)
              wpm2 = 1./((batfrq-dtmin/60.)*3600.)
            else
              mmpd = 24./batfrq
              wpm2 = 1./(batfrq*3600.)
            end if
 
            do i = 2 , ixm1
              ci = i
 
              drag_o(j,i-1) = 0.0
              q2m_o(j,i-1) = 0.0
              evpa_o(j,i-1) = 0.0
              sena_o(j,i-1) = 0.0
              do n = 1 , nnsg
                if ( ocld2d(n,i,j)>=0.5 ) then
                  q2m_s(n,j,i-1) = q2d(i,j)
                  drag_s(n,j,i-1) = uvdrag(i,j)
                  evpa_s(n,j,i-1) = evpa2d(n,ci,j)*mmpd
                  sena_s(n,j,i-1) = sena2d(n,ci,j)*wpm2
                  tpr_s(n,j,i-1) = (prnca2d(ci,j)+prca2d(ci,j))*mmpd
                  prcv_s(n,j,i-1) = prca2d(ci,j)*mmpd
                  ps_s(n,j,i-1) = p1d(n,i)*0.01
 
                  q2m_o(j,i-1) = q2m_o(j,i-1) + q2d(i,j)
                  drag_o(j,i-1) = drag_o(j,i-1) + uvdrag(i,j)
                  evpa_o(j,i-1) = evpa_o(j,i-1) + evpa2d(n,ci,j)
                  sena_o(j,i-1) = sena_o(j,i-1) + sena2d(n,ci,j)
                else if ( ocld2d(n,i,j)<=0.5 ) then
                  q2m_s(n,j,i-1) = q2m_1d(n,i)
                  drag_s(n,j,i-1) = drag1d(n,i)
                  evpa_s(n,j,i-1) = evpa2d(n,i,j)*mmpd
                  sena_s(n,j,i-1) = sena2d(n,i,j)*wpm2
                  tpr_s(n,j,i-1) = (prnca2d(i,j)+prca2d(i,j))*mmpd
                  prcv_s(n,j,i-1) = prca2d(i,j)*mmpd
                  ps_s(n,j,i-1) = p1d(n,i)*0.01
 
                  q2m_o(j,i-1) = q2m_o(j,i-1) + q2m_1d(n,i)
                  drag_o(j,i-1) = drag_o(j,i-1) + drag1d(n,i)
                  evpa_o(j,i-1) = evpa_o(j,i-1) + evpa2d(n,i,j)
                  sena_o(j,i-1) = sena_o(j,i-1) + sena2d(n,i,j)
                else
                end if
              end do
              tpr_o(j,i-1) = (prnca2d(ci,j)+prca2d(ci,j))*mmpd
              q2m_o(j,i-1) = q2m_o(j,i-1)/float(nnsg)
              drag_o(j,i-1) = drag_o(j,i-1)/float(nnsg)
              evpa_o(j,i-1) = evpa_o(j,i-1)/float(nnsg)*mmpd
              sena_o(j,i-1) = sena_o(j,i-1)/float(nnsg)*wpm2
              flwa_o(j,i-1) = flwa2d(ci,j)*wpm2
              fswa_o(j,i-1) = fswa2d(ci,j)*wpm2
              flwd_o(j,i-1) = flwda2d(ci,j)*wpm2
              sina_o(j,i-1) = sina2d(ci,j)*wpm2
              prcv_o(j,i-1) = prca2d(ci,j)*mmpd
              ps_o(j,i-1) = (psb(i,j)+ptop)*10.
              zpbl_o(j,i-1) = zpbl(i,j)
 
              tlef_o(j,i-1) = 0.0
              ssw_o(j,i-1) = 0.0
              rsw_o(j,i-1) = 0.0
              rnos_o(j,i-1) = 0.0
              scv_o(j,i-1) = 0.0
              nnn = 0
              do n = 1 , nnsg
!abt           if(ocld2d(n,ci,j).ge.0.5) then
                if ( ocld2d(n,ci,j)>=0.5 .and. landmask(jj,ci)/=3 ) then
                  tlef_o(j,i-1) = tlef_o(j,i-1) + c2rtlef(jj,ci)
                  ssw_o(j,i-1) = ssw_o(j,i-1) + c2rsm10cm(jj,ci)
                  rsw_o(j,i-1) = rsw_o(j,i-1) + c2rsm1m(jj,ci)
                  rnos_o(j,i-1) = rnos_o(j,i-1) + rnos2d(n,ci,j)
                  scv_o(j,i-1) = scv_o(j,i-1) + c2rsnowc(jj,ci)
                  tlef_s(n,j,i-1) = c2rtlef(jj,ci)
                  ssw_s(n,j,i-1) = c2rsm10cm(jj,ci)
                  rsw_s(n,j,i-1) = c2rsm1m(jj,ci)
                  rnos_s(n,j,i-1) = rnos2d(n,ci,j)*mmpd
                  scv_s(n,j,i-1) = c2rsnowc(jj,ci)
                  nnn = nnn + 1
                else
                  tlef_s(n,j,i-1) = -1.E34
                  ssw_s(n,j,i-1) = -1.E34
                  rsw_s(n,j,i-1) = -1.E34
                  rnos_s(n,j,i-1) = -1.E34
                  scv_s(n,j,i-1) = -1.E34
                end if
              end do
              if ( nnn>=max0(nnsg/2,1) ) then
                tlef_o(j,i-1) = tlef_o(j,i-1)/float(nnn)
                ssw_o(j,i-1) = ssw_o(j,i-1)/float(nnn)
                rsw_o(j,i-1) = rsw_o(j,i-1)/float(nnn)
                rnos_o(j,i-1) = rnos_o(j,i-1)/float(nnn)*mmpd
                scv_o(j,i-1) = scv_o(j,i-1)/float(nnn)
              else
                tlef_o(j,i-1) = -1.E34
                ssw_o(j,i-1) = -1.E34
                rsw_o(j,i-1) = -1.E34
                rnos_o(j,i-1) = -1.E34
                scv_o(j,i-1) = -1.E34
              end if
!             ******    reset accumulation arrays to zero
              do n = 1 , nnsg
                evpa2d(n,ci,j) = 0.
                rnos2d(n,ci,j) = 0.
                sena2d(n,ci,j) = 0.
              end do
              prnca2d(ci,j) = 0.
              prca2d(ci,j) = 0.
              flwa2d(ci,j) = 0.
              flwda2d(ci,j) = 0.
              fswa2d(ci,j) = 0.
              svga2d(ci,j) = 0.
              sina2d(ci,j) = 0.
 
            end do  ! end of i loop
          end if    ! end if jyear eq jyearr
        end do      ! end of j loop
      else          ! end if inout = 2
      end if
 
      end subroutine interfclm

#endif
