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

      subroutine output(iexec)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                     c
c     this subroutine handles all of the output                       c
c                                                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer iexec
#     include "parameter.inc"
#ifdef MPP1
#     include "mpiregcm.inc"
      integer ierr
      integer idum
#endif
#     include "include/param1.cb"
#     include "include/param2.cb"
#     include "include/param3.cb"
#     include "include/iunits.cb"
#     include "include/main.cb"
#     include "include/mainchem.cb"
#     include "include/bats.cb2"
#ifdef MPP1
#     include "include/bdycod.cb"
#     include "include/pmoist.cb"
#     include "include/rad.cb"
#     include "include/trachem.cb"
#     include "include/cvaria.cb"
#     include "include/io.cb"
#endif
#     include "include/date.cb"
#     include "include/message.cb"
c**********************
      integer nrcout
      COMMON /GENREC/ nrcout
#ifdef MPP1
      integer m,n
      real(kind=8)  am
      real(kind=8)  an
      real(kind=8)  dstor
      real(kind=8)  hstor
      common /dpass/ m(nsplit)
     a             , am(kx,nsplit), an(nsplit)
     b             , dstor(ix,0:jxp+1,nsplit), hstor(ix,0:jxp+1,nsplit)
#     include "include/parrad.cb"
      real(kind=8)  absnxt,              ! Nearest layer absorptivities
     $     abstot,              ! Non-adjacent layer absorptivites
     $     emstot               ! Total emissivity
      common /radbuf/ absnxt(plond,plev,4,jxp)
     &      , abstot(plond,plevp,plevp,jxp),emstot(plond,plevp,jxp)
      real(kind=8)  absnxt_io
      real(kind=8)  abstot_io
      real(kind=8)  emstot_io
      common /radbuf_io/absnxt_io(plond,plev,4,mjx-1)
     &                 ,abstot_io(plond,plevp,plevp,mjx-1)
     &                 ,emstot_io(plond,plevp,mjx-1)
      integer numbat
      parameter( numbat=21+6 )
      real(kind=4)  fbat(jxp,ix-2,numbat)
      equivalence (fbat(1,1,1),u10m_o(1,1))
      real(kind=4)  fbat_io(mjx-2,ix-2,numbat)
      COMMON /BATOUTio/ fbat_io
      integer numsub
      parameter( numsub=16 )
      real(kind=4)  fsub(NNSG,jxp,ix-2,numsub)
      equivalence (fsub(1,1,1,1),u10m_s(1,1,1))
      real(kind=4)  fsub_io(NNSG,mjx-2,ix-2,numsub)
      COMMON /SUBOUTio/ fsub_io
      integer nrad2d,nrad3d
      parameter(nrad2d=21,nrad3d=5)
      real(kind=4)  frad2d(jxp,ix-2,nrad2d),frad3d(jxp,ix-2,kx,nrad3d)
      COMMON /OUTRAD/ frad2d,frad3d
      real(kind=4)  frad2d_io(mjx-2,ix-2,nrad2d)
      real(kind=4)  frad3d_io(mjx-2,ix-2,kx,nrad3d)
      COMMON /OUTRADio/ frad2d_io,frad3d_io
      real(kind=8)  out0(ix,3,jxp),out_0(ix,3,mjx)
      real(kind=8)  atm0(ix,kx*6+3+NNSG*4,jxp)
      real(kind=8)  atm_0(ix,kx*6+3+NNSG*4,mjx)
      real(kind=4)  bat0(ix-2,numbat,jxp),bat_0(ix-2,numbat,mjx)
      real(kind=4)  sub0(ix-2,NNSG,numsub,jxp)
      real(kind=8)  sub_0(ix-2,NNSG,numsub,mjx)
      real(kind=4)  rad0(ix-2,nrad3d*kx+nrad2d,jxp)
      real(kind=4)  rad_0(ix-2,nrad3d*kx+nrad2d,mjx)
      COMMON /TMPLOC/ out0,out_0,atm0,atm_0,bat0,bat_0,sub0,sub_0
     &               ,rad0,rad_0
      real(kind=8)   chem0(ix,ntr*kx+kx*3+ntr*7+3,jxp)
      real(kind=8)  chem_0(ix,ntr*kx+kx*3+ntr*7+3,mjx)
      COMMON /TMPLOC1/ chem0,chem_0
      real(kind=8)   sav0(ix,kx*4+2,jxp),sav_0(ix,kx*4+2,mjx)
      real(kind=8)   sav0a(ix,kx+NNSG+5,jxp),sav_0a(ix,kx+NNSG+5,mjx)
      real(kind=8)   sav0c(ix,kx*2,jxp),sav_0c(ix,kx*2,mjx)
      real(kind=8)   sav0b(ix,kx+1,jxp),sav_0b(ix,kx+1,mjx)
      real(kind=8)   sav0s(ix,kx,jxp),sav_0s(ix,kx,mjx)
      COMMON /TMPSAV0/ sav0,sav_0,sav0a,sav_0a,sav0c,sav_0c
     &               ,sav0b,sav_0b,sav0s,sav_0s
      real(kind=8)   sav0d(ix,nsplit*2,jxp),sav_0d(ix,nsplit*2,mjx)
      COMMON /TMPSAV0d/ sav0d,sav_0d
      real(kind=8)   sav1(ilx,kx*4+(kx+1)*(kx+2),jxp)
      real(kind=8)  sav_1(ilx,kx*4+(kx+1)*(kx+2),mjx)
      COMMON /TMPSAV1/ sav1,sav_1
      real(kind=8)   sav2(ilx,NNSG*4+4,jxp),sav_2(ilx,NNSG*4+4,mjx)
      real(kind=8)   sav2a(ilx,NNSG*5+1,jxp),sav_2a(ilx,NNSG*5+1,mjx)
      COMMON /TMPSAV2/ sav2,sav_2,sav2a,sav_2a
      real(kind=8)   sav4(ix,ntr*(kx*4+1),jxp),sav4a(ilx,7,jxp)
      real(kind=8)  sav_4(ix,ntr*(kx*4+1),mjx),sav_4a(ilx,7,mjx)
      COMMON /TMPSAV4/ sav4,sav_4,sav4a,sav_4a
      real(kind=8)   sav6(kx,8,jxp),sav_6(kx,8,mjx)
      COMMON /TMPSAV6/ sav6,sav_6
      integer allrec,k,l
#else
      integer nrad2d,nrad3d
      parameter(nrad2d=21,nrad3d=5)
      real(kind=4)  frad2d(jx-2,ix-2,nrad2d),frad3d(jx-2,ix-2,kx,nrad3d)
      COMMON /OUTRAD/ frad2d,frad3d
#endif
      integer i,j
      integer nrcrad
      COMMON /RADREC/ nrcrad
c
Cchem2
      integer nrcchem
      COMMON /CHEREC/ nrcchem
Cchem2_
      character itype*3, newfil*14, tmpfil*17, rmsav*70
      logical there
c
      integer isystm,system
      external system
c----------------------------------------------------------------------

      if (.not.(jyear.eq.jyear0.and.ktau.eq.0)) then
        if (mod(nint(xtime),60).lt.mod(nint(xtime-dtmin),60)) then
          idatex = idatex + 1
        end if
        if(dabs(xtime).lt.0.00001) idatex=ldatez
      end if

#ifdef MPP1
      if (myid .eq. 0) then
#endif
      if (jyear.eq.jyearr.and.ktau.eq.ktaur) then
        if (iotyp.eq.1) then
          print*,'Writing output files in direct access format'
        else if (iotyp.eq.2) then
          print*,'Writing output files in sequential format'
        else
          write(aline, *) 'iotyp = ', iotyp
          call say
          call fatal(__FILE__,__LINE__, 'Output format does not exist')
        end if
        inquire(file='output',exist=there)
        if (.not.there) then
          isystm=system('mkdir output/')
        end if 
 199    format(a3,'.',i10)
        nrcout=0
        call mkfile

      end if
#ifdef MPP1
      end if
#endif

      if(jyear.eq.jyear0.and.ktau.eq.0) then
#ifdef MPP1
!=======================================================================
!     among    ht,htsd,veg2d,satbrt,xlat,xlong,msfx,msfd,f
!     we just need gather                      msfx,msfd
        do j=1,jendl
        do i=1,ilx
          out0(i,1,j) =veg2d(i,j)
        end do
        end do
        do j=1,jendl
        do i=1,ix
          out0(i,2,j) = msfx(i,j)
          out0(i,3,j) = msfd(i,j)
        end do
        end do
        call MPI_GATHER(out0(1,1,1),ix*3*jxp,MPI_DOUBLE_PRECISION,
     &                 out_0(1,1,1),ix*3*jxp,MPI_DOUBLE_PRECISION,
     &                                     0,MPI_COMM_WORLD,ierr)
        if (myid .eq. 0) then
          do j=1,mjx-1
          do i=1,ilx
           veg2d_io(i,j) = out_0(i,1,j)
          end do
          end do
          do j=1,mjx
          do i=1,ix
            msfx_io(i,j) = out_0(i,2,j)
            msfd_io(i,j) = out_0(i,3,j)
          end do
          end do
#endif
          call outtap0
          call gradsctl('OUT_HEAD.CTL')
#ifdef MPP1
        end if
#endif
      end if

c
c-----output for dataflow analyses:
c
      if(iftape) then
        if((jyear.eq.jyear0.and.ktau.eq.0) .or.
     &     (mod(ntime,ntapfrq).eq.0.and.
     &      (.not.(jyear.eq.jyearr.and.ktau.eq.ktaur)))) then
#ifdef MPP1
!=======================================================================
!   gather  ua,va,ta,qva,qca,rainc,rainnc,tgb2d,swt2d,olcd2d,rno2d
          do j=1,jendl
            do k=1,kx
            do i=1,ix
              atm0(i,k     ,j) =    ua(i,k,j)
              atm0(i,k+kx  ,j) =    va(i,k,j)
              atm0(i,k+kx*2,j) = omega(i,k,j)
              atm0(i,k+kx*3,j) =    ta(i,k,j)
              atm0(i,k+kx*4,j) =   qva(i,k,j)
              atm0(i,k+kx*5,j) =   qca(i,k,j)
            end do
            end do
            do i=1,ix
              atm0(i,1+kx*6,j) =   psa(i,j)
              atm0(i,2+kx*6,j) = rainc(i,j)
              atm0(i,3+kx*6,j) =rainnc(i,j)
            end do
          end do
          do j=1,jendx
            do n=1,NNSG
            do i=1,ilx
              atm0(i,3+kx*6+n     ,j) = ocld2d(n,i,j)
              atm0(i,3+kx*6+n+NNSG  ,j) =  tgb2d(n,i,j)
              atm0(i,3+kx*6+n+NNSG*2,j) =  swt2d(n,i,j)
              atm0(i,3+kx*6+n+NNSG*3,j) =  rno2d(n,i,j)
            end do
            end do
          end do
          call MPI_GATHER(atm0(1,1,1),ix*(kx*6+3+NNSG*4)*jxp
     &        ,MPI_DOUBLE_PRECISION,atm_0(1,1,1),ix*(kx*6+3+NNSG*4)*jxp
     &        ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
          if (myid .eq. 0) then
            do j=1,mjx
              do k=1,kx
              do i=1,ix
                ua_io(i,k,j) = atm_0(i,k     ,j)
                va_io(i,k,j) = atm_0(i,k+kx  ,j)
             omega_io(i,k,j) = atm_0(i,k+kx*2,j)
                ta_io(i,k,j) = atm_0(i,k+kx*3,j)
                qva_io(i,k,j) =atm_0(i,k+kx*4,j)
                qca_io(i,k,j) =atm_0(i,k+kx*5,j)
              end do
              end do
              do i=1,ix
                psa_io(i,j)    = atm_0(i,1+kx*6,j)
                rainc_io(i,j)  = atm_0(i,2+kx*6,j)
                rainnc_io(i,j) = atm_0(i,3+kx*6,j)
              end do
            end do
            do j=1,mjx-1
              do n=1,NNSG
              do i=1,ilx
                ocld2d_io(n,i,j) = atm_0(i,3+kx*6+n       ,j)
                tgb2d_io(n,i,j)  = atm_0(i,3+kx*6+n+NNSG  ,j)
                swt2d_io(n,i,j)  = atm_0(i,3+kx*6+n+NNSG*2,j)
                rno2d_io(n,i,j)  = atm_0(i,3+kx*6+n+NNSG*3,j)
              end do
              end do
            end do
#endif
          call outtap
#ifdef MPP1
          end if
          do j=1,jendx
          do i=1,ilx
            do n=1,NNSG
              rno2d(n,i,j) = 0.
            end do
            rainc(i,j) = 0.
            rainnc(i,j) = 0.
          end do
          end do
#endif
        end if
      end if


c   Call surface output
      if(ifbat) then
        if((mod(ntime,kbats).eq.0.and.
     &      (.not.(jyear.eq.jyearr.and.ktau.eq.ktaur)))
     &      .or. (jyear.eq.jyear0.and.ktau.eq.1)) then
#ifdef MPP1
          do j=1,jendx
            do l=1,numbat
            do i=1,ix-2
              bat0(i,l,j) = fbat(j,i,l)
            end do
            end do
          end do
          call MPI_GATHER(bat0(1,1,1),(ix-2)*numbat*jxp,MPI_REAL,
     &                   bat_0(1,1,1),(ix-2)*numbat*jxp,MPI_REAL,
     &                       0,MPI_COMM_WORLD,ierr)
          if(myid.eq.0) then
            do l=1,numbat
              do i=1,ix-2
              do j=1,mjx-2
                fbat_io(j,i,l)=bat_0(i,l,j+1)
              end do
              end do
            end do
#endif
          call outsrf
#ifdef MPP1
          end if
          do i=1,ix-2
          do j=1,jxp
            tgmx_o(j,i) = -1.e30
            t2mx_o(j,i) = -1.e30
            tgmn_o(j,i) =  1.e30
            t2mn_o(j,i) =  1.e30
            w10x_o(j,i) = -1.e30
            psmn_o(j,i) =  1.e30
          end do
          end do
#else
          do i=1,ix-2
          do j=1,jx-2
            tgmx_o(j,i) = -1.e30
            t2mx_o(j,i) = -1.e30
            tgmn_o(j,i) =  1.e30
            t2mn_o(j,i) =  1.e30
            w10x_o(j,i) = -1.e30
            psmn_o(j,i) =  1.e30
          end do
          end do
#endif
#ifdef MPP1
          do j=1,jendx
            do l=1,numsub
              do n=1,NNSG
              do i=1,ix-2
                sub0(i,n,l,j) = fsub(n,j,i,l)
              end do
              end do
            end do
          end do
          call MPI_GATHER(sub0(1,1,1,1),(ix-2)*NNSG*numsub*jxp,MPI_REAL,
     &                   sub_0(1,1,1,1),(ix-2)*NNSG*numsub*jxp,MPI_REAL,
     &                       0,MPI_COMM_WORLD,ierr)
          if(myid.eq.0) then
            do l=1,numsub
              do j=1,mjx-2
              do n=1,NNSG
              do i=1,ix-2
                fsub_io(n,j,i,l)=sub_0(i,n,l,j+1)
              end do
              end do
              end do
            end do
#endif
          if(ifsub.and.NSG.gt.1) call outsub
#ifdef MPP1
          end if
#endif
        end if
      end if

c   Call radiation output
      if (ifrad) then
        if ((mod(ntime,nradisp).eq.0.and.
     &       (.not.(jyear.eq.jyearr.and.ktau.eq.ktaur)))
     &      .or. (jyear.eq.jyear0.and.ktau.eq.1)) then
#ifdef MPP1
!=======================================================================
!      frad2d, frad3d
          do n=1,nrad2d
            do j=1,jxp
            do i=1,ix-2
              rad0(i,n,j)=frad2d(j,i,n)
            end do
            end do
          end do
          do n=1,nrad3d
            do k=1,kx
              do j=1,jxp
              do i=1,ix-2
                rad0(i,nrad2d+(n-1)*kx+k,j) = frad3d(j,i,k,n)
              end do
              end do
            end do
          end do
          call MPI_GATHER(rad0(1,1,1),(ix-2)*(nrad3d*kx+nrad2d)*jxp,
     &         MPI_REAL, rad_0(1,1,1),(ix-2)*(nrad3d*kx+nrad2d)*jxp,
     &         MPI_REAL, 0,MPI_COMM_WORLD,ierr)
          if(myid.eq.0) then
            do n=1,nrad2d
              do j=1,mjx-2
              do i=1,ix-2
                frad2d_io(j,i,n)=rad_0(i,n,j+1)
              end do
              end do
            end do
            do n=1,nrad3d
              do k=1,kx
                do j=1,mjx-2
                do i=1,ix-2
                  frad3d_io(j,i,k,n) = rad_0(i,nrad2d+(n-1)*kx+k,j+1)
                end do
                end do
              end do
            end do
#endif
          call radtap
#ifdef MPP1
          end if
#endif
        end if
      end if

Cchem2
c   Call chem output
      if(ifchem) then
       if((jyear.eq.jyear0.and.ktau.eq.1).or.
     &    (mod(ntime,kchem).eq.0.and.
     &     (.not.(jyear.eq.jyearr.and.ktau.eq.ktaur)))) then
#ifdef MPP1
         do j=1,jendl
           do n=1,ntr
             do k=1,kx
               do i=1,ix
                 chem0(i,(n-1)*kx+k,j) = chia(i,k,j,n)
               end do
             end do
           end do
         end do
         do j=1,jendx
           do k=1,kx
             do i=1,ilx
               chem0(i,ntr*kx+k     ,j) = aerext(i,k,j)
               chem0(i,ntr*kx+kx+k  ,j) = aerssa(i,k,j)
               chem0(i,ntr*kx+kx*2+k,j) = aerasp(i,k,j)
             end do
           end do
         end do
         do j=1,jendl
           do n=1,ntr
             do i=1,ix
               chem0(i,(ntr+3)*kx+n      ,j) =  dtrace(i,j,n)
               chem0(i,(ntr+3)*kx+ntr+n  ,j) =   wdlsc(i,j,n)
               chem0(i,(ntr+3)*kx+ntr*2+n,j) =   wdcvc(i,j,n)
               chem0(i,(ntr+3)*kx+ntr*3+n,j) =   ddsfc(i,j,n)
               chem0(i,(ntr+3)*kx+ntr*4+n,j) =    wxsg(i,j,n)
               chem0(i,(ntr+3)*kx+ntr*5+n,j) =    wxaq(i,j,n)
               chem0(i,(ntr+3)*kx+ntr*6+n,j) = cemtrac(i,j,n)
             end do
           end do
         end do
         do j=1,jendx
           do i=1,ilx
             chem0(i,(ntr+3)*kx+ntr*7+1,j) = aertarf(i,j)
             chem0(i,(ntr+3)*kx+ntr*7+2,j) = aersrrf(i,j)
           end do
         end do
         do j=1,jendl
           do i=1,ix
             chem0(i,(ntr+3)*kx+ntr*7+3,j) = psa(i,j)
           end do
         end do
         call MPI_GATHER(chem0(1,1,1),ix*((ntr+3)*kx+ntr*7+3)*jxp,
     &    MPI_DOUBLE_PRECISION,chem_0(1,1,1),ix*((ntr+3)*kx+ntr*7+3)*jxp
     &   ,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
         if(myid.eq.0) then
           do j=1,mjx
             do n=1,ntr
               do k=1,kx
                 do i=1,ix
                   chia_io(i,k,j,n) = chem_0(i,(n-1)*kx+k,j)
                 end do
               end do
             end do
           end do
           do j=1,mjx-1
             do k=1,kx
               do i=1,ilx
                 aerext_io(i,k,j) = chem_0(i,ntr*kx+k     ,j+1)
                 aerssa_io(i,k,j) = chem_0(i,ntr*kx+kx+k  ,j+1)
                 aerasp_io(i,k,j) = chem_0(i,ntr*kx+kx*2+k,j+1)
               end do
             end do
           end do
           do j=1,mjx
             do n=1,ntr
               do i=1,ix
                 dtrace_io(i,j,n)  = chem_0(i,(ntr+3)*kx+n      ,j)
                 wdlsc_io(i,j,n)   = chem_0(i,(ntr+3)*kx+ntr+n  ,j)
                 wdcvc_io(i,j,n)   = chem_0(i,(ntr+3)*kx+ntr*2+n,j)
                 ddsfc_io(i,j,n)   = chem_0(i,(ntr+3)*kx+ntr*3+n,j)
                 wxsg_io(i,j,n)    = chem_0(i,(ntr+3)*kx+ntr*4+n,j)
                 wxaq_io(i,j,n)    = chem_0(i,(ntr+3)*kx+ntr*5+n,j)
                 cemtrac_io(i,j,n) = chem_0(i,(ntr+3)*kx+ntr*6+n,j)
               end do
             end do
           end do
           do j=1,mjx-1
             do i=1,ilx
               aertarf_io(i,j) = chem_0(i,(ntr+3)*kx+ntr*7+1,j+1)
               aersrrf_io(i,j) = chem_0(i,(ntr+3)*kx+ntr*7+2,j+1)
             end do
           end do
           do j=1,mjx
             do i=1,ix
               psa_io(i,j) = chem_0(i,(ntr+3)*kx+ntr*7+3,j)
             end do
           end do
#endif
         call chemtap
#ifdef MPP1
         end if
         do n=1,ntr
         do j=1,jendl
         do k=1,kx
         do i=1,ix
           remlsc(i,k,j,n) = 0.
           remcvc(i,k,j,n) = 0.
           rxsg(i,k,j,n) = 0.
           rxsaq1(i,k,j,n) = 0.
           rxsaq2(i,k,j,n) = 0.
         end do
         end do
         end do
         end do
         do n=1,ntr
         do j=1,jendl
         do i=1,ix
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
         do j=1,jendl
         do i=1,ilx
           aertarf(i,j)=0.
           aersrrf(i,j)=0.
         end do
         end do
#endif
       end if
      end if
Cchem2
c
c-----output for restart:
c
      if (ifsave) then
        if(((lday.eq.1.and.lhour.eq.0.and.dabs(xtime).lt.0.00001)
     &     .and.ldatez.ne.IDATE1) .or. nnnnnn.eq.nnnend) then
#ifdef MPP1
          do j=1,jendl
            do k=1,kx
              do i=1,ix
                sav0(i,      k,j) = UB0(i,k,j)
                sav0(i,kx   +k,j) = VB0(i,k,j)
                sav0(i,kx*2 +k,j) = QB0(i,k,j)
                sav0(i,kx*3 +k,j) = TB0(i,k,j)
              end do
            end do
            do i=1,ix
              sav0(i,kx*4+1 ,j) =    PS0(i,j)
              sav0(i,kx*4+2 ,j) =    TS0(i,j)
            end do
          end do
          allrec = kx*4+2
          call MPI_GATHER(
     &             sav0(1,1,1),ix*allrec*jxp,MPI_DOUBLE_PRECISION,
     &            sav_0(1,1,1),ix*allrec*jxp,MPI_DOUBLE_PRECISION,
     &                                     0,MPI_COMM_WORLD,ierr)
          if(myid.eq.0) then
            do j=1,mjx
              do k=1,kx
                do i=1,ix
                  UB0_io(i,k,j) = sav_0(i,      k,j)
                  VB0_io(i,k,j) = sav_0(i,kx   +k,j)
                  QB0_io(i,k,j) = sav_0(i,kx*2 +k,j)
                  TB0_io(i,k,j) = sav_0(i,kx*3 +k,j)
                end do
              end do
              do i=1,ix
                PS0_io(i,j) =    sav_0(i,kx*4+1 ,j)
                TS0_io(i,j) =    sav_0(i,kx*4+2 ,j)
              end do
            end do
          end if
          if(EHSO4) then
             do j=1,jendl
               do k=1,kx
                 do i=1,ix
                   sav0s(i,k,j) = SO0(i,k,j)
                 end do
               end do
             end do
             call MPI_GATHER(
     &                sav0s(1,1,1),ix*kx*jxp,MPI_DOUBLE_PRECISION,
     &               sav_0s(1,1,1),ix*kx*jxp,MPI_DOUBLE_PRECISION,
     &                                     0,MPI_COMM_WORLD,ierr)
             if(myid.eq.0) then
               do j=1,mjx
                 do k=1,kx
                   do i=1,ix
                     SO0_io(i,k,j) = sav_0s(i,k,j)
                   end do
                 end do
               end do
             end if
          end if
          do j=1,jendl
            do k=1,kx
              do i=1,ix
                sav0(i,      k,j) = ua(i,k,j)
                sav0(i,kx   +k,j) = ub(i,k,j)
                sav0(i,kx*2 +k,j) = va(i,k,j)
                sav0(i,kx*3 +k,j) = vb(i,k,j)
              end do
            end do
            do i=1,ix
              sav0(i,kx*4+1 ,j) = psa(i,j)
              sav0(i,kx*4+2 ,j) = psb(i,j)
            end do
          end do
          allrec = kx*4+2
          call MPI_GATHER(
     &             sav0(1,1,1),ix*allrec*jxp,MPI_DOUBLE_PRECISION,
     &            sav_0(1,1,1),ix*allrec*jxp,MPI_DOUBLE_PRECISION,
     &                                     0,MPI_COMM_WORLD,ierr)
          if(myid.eq.0) then
            do j=1,mjx
              do k=1,kx
                do i=1,ix
                  ua_io(i,k,j) = sav_0(i,      k,j)
                  ub_io(i,k,j) = sav_0(i,kx   +k,j)
                  va_io(i,k,j) = sav_0(i,kx*2 +k,j)
                  vb_io(i,k,j) = sav_0(i,kx*3 +k,j)
                end do
              end do
              do i=1,ix
                psa_io(i,j) =    sav_0(i,kx*4+1 ,j)
                psb_io(i,j) =    sav_0(i,kx*4+2 ,j)
              end do
            end do
          end if
          do j=1,jendl
            do k=1,kx
              do i=1,ix
                sav0(i,      k,j) = ta(i,k,j)
                sav0(i,kx   +k,j) = tb(i,k,j)
                sav0(i,kx*2 +k,j) =qva(i,k,j)
                sav0(i,kx*3 +k,j) =qvb(i,k,j)
              end do
            end do
            do i=1,ix
              sav0(i,kx*4+1 ,j) = tga(i,j)
              sav0(i,kx*4+2 ,j) = tgb(i,j)
            end do
          end do
          allrec = kx*4+2
          call MPI_GATHER(
     &             sav0(1,1,1),ix*allrec*jxp,MPI_DOUBLE_PRECISION,
     &            sav_0(1,1,1),ix*allrec*jxp,MPI_DOUBLE_PRECISION,
     &                                     0,MPI_COMM_WORLD,ierr)
          if(myid.eq.0) then
            do j=1,mjx
              do k=1,kx
                do i=1,ix
                  ta_io(i,k,j) = sav_0(i,      k,j)
                  tb_io(i,k,j) = sav_0(i,kx   +k,j)
                 qva_io(i,k,j) = sav_0(i,kx*2 +k,j)
                 qvb_io(i,k,j) = sav_0(i,kx*3 +k,j)
                end do
              end do
              do i=1,ix
                tga_io(i,j) =    sav_0(i,kx*4+1 ,j)
                tgb_io(i,j) =    sav_0(i,kx*4+2 ,j)
              end do
            end do
          end if
          do j=1,jendl
            do k=1,kx
              do i=1,ix
                sav0(i,      k,j) =qca(i,k,j)
                sav0(i,kx   +k,j) =qcb(i,k,j)
                sav0(i,kx*2 +k,j) =fcc(i,k,j)
              end do
            end do
            do i=1,ix
              sav0(i,kx*4+1 ,j) =  rainc(i,j)
              sav0(i,kx*4+2 ,j) = rainnc(i,j)
            end do
          end do
          do j=1,jendx
            do k=1,kx
              do i=1,ilx
                sav0(i,kx*3 +k,j) = heatrt(i,k,j)
              end do
            end do
          end do
          allrec = kx*4+2
          call MPI_GATHER(
     &             sav0(1,1,1),ix*allrec*jxp,MPI_DOUBLE_PRECISION,
     &            sav_0(1,1,1),ix*allrec*jxp,MPI_DOUBLE_PRECISION,
     &                                     0,MPI_COMM_WORLD,ierr)
          if(myid.eq.0) then
            do j=1,mjx
              do k=1,kx
                do i=1,ix
                 qca_io(i,k,j) = sav_0(i,      k,j)
                 qcb_io(i,k,j) = sav_0(i,kx   +k,j)
                 fcc_io(i,k,j) = sav_0(i,kx*2 +k,j)
                end do
              end do
              do i=1,ix
                rainc_io(i,j) =    sav_0(i,kx*4+1 ,j)
               rainnc_io(i,j) =    sav_0(i,kx*4+2 ,j)
              end do
            end do
            do j=1,mjx-1
              do k=1,kx
                do i=1,ilx
                  heatrt_io(i,k,j) = sav_0(i,kx*3 +k,j)
                end do
              end do
            end do
          end if
          do j=1,jendl
            do i=1,ix
              sav0a(i,1 ,j) =    hfx(i,j)
              sav0a(i,2 ,j) =    qfx(i,j)
              sav0a(i,3 ,j) = uvdrag(i,j)
              sav0a(i,4 ,j) =   tgbb(i,j)
            end do
            do n=1,NNSG
              do i=1,ix
                sav0a(i,4+n,j) = snowc(n,i,j)
              end do
            end do
          end do
          do j=1,jendx
            do k=1,kx+1
              do i=1,ilx
                sav0a(i,NNSG+4 +k,j) = o3prof(i,k,j)
              end do
            end do
          end do
          allrec = 5+NNSG+kx
          call MPI_GATHER(
     &             sav0a(1,1,1),ix*allrec*jxp,MPI_DOUBLE_PRECISION,
     &            sav_0a(1,1,1),ix*allrec*jxp,MPI_DOUBLE_PRECISION,
     &                                      0,MPI_COMM_WORLD,ierr)
          if(myid.eq.0) then
            do j=1,mjx
              do i=1,ix
                  hfx_io(i,j) =    sav_0a(i,1 ,j)
                  qfx_io(i,j) =    sav_0a(i,2 ,j)
               uvdrag_io(i,j) =    sav_0a(i,3 ,j)
                 tgbb_io(i,j) =    sav_0a(i,4 ,j)
              end do
              do n=1,NNSG
                do i=1,ix
                  snowc_io(n,i,j) = sav_0a(i,4+n,j)
                end do
              end do
            end do
            do j=1,mjx-1
              do k=1,kx+1
                do i=1,ilx
                  o3prof_io(i,k,j) = sav_0a(i,4+NNSG+k,j)
                end do
              end do
            end do
          end if
          if(iocnflx.eq.2) then
            call MPI_GATHER( zpbl(1,1),ix*jxp,MPI_DOUBLE_PRECISION,
     &                    zpbl_io(1,1),ix*jxp,MPI_DOUBLE_PRECISION,
     &                                        0,MPI_COMM_WORLD,ierr)
          end if
          if(icup.eq.1) then
            do j=1,jendl
              do k=1,kx
                do i=1,ix
                  sav0c(i,   k,j) = rsheat(i,k,j) 
                  sav0c(i,kx+k,j) =  rswat(i,k,j) 
                end do
              end do
            end do
            allrec = kx*2
            call MPI_GATHER(
     &               sav0c(1,1,1),ix*allrec*jxp,MPI_DOUBLE_PRECISION,
     &              sav_0c(1,1,1),ix*allrec*jxp,MPI_DOUBLE_PRECISION,
     &                                       0,MPI_COMM_WORLD,ierr)
            if(myid.eq.0) then
              do j=1,mjx
                do k=1,kx
                  do i=1,ix
                    rsheat_io(i,k,j) = sav_0c(i,   k,j)
                     rswat_io(i,k,j) = sav_0c(i,kx+k,j)
                  end do
                end do
              end do
            end if
          else if(icup.eq.3) then
            do j=1,jendl
              do k=1,kx
                do i=1,ix
                  sav0b(i,k,j) = TBASE(i,k,j)
                end do
              end do
              do i=1,ix
                sav0b(i,kx+1,j) = CLDEFI(i,j)
              end do
            end do
            allrec = kx+1
            call MPI_GATHER(
     &               sav0b(1,1,1),ix*allrec*jxp,MPI_DOUBLE_PRECISION,
     &              sav_0b(1,1,1),ix*allrec*jxp,MPI_DOUBLE_PRECISION,
     &                                       0,MPI_COMM_WORLD,ierr)
            if(myid.eq.0) then
              do j=1,mjx
                do k=1,kx
                  do i=1,ix
                    TBASE_IO(i,k,j) = sav_0b(i,   k,j)
                  end do
                end do
                do i=1,ix
                  CLDEFI_IO(i,j) = sav_0b(i,kx+1,j)
                end do
              end do
            end if
          else if(icup.eq.4) then
            call MPI_GATHER( cbmf2d(1,1),ix*jxp,MPI_DOUBLE_PRECISION,
     &                    cbmf2d_io(1,1),ix*jxp,MPI_DOUBLE_PRECISION,
     &                                        0,MPI_COMM_WORLD,ierr)
          end if
          do j=1,jendx
            do l=1,4
              do k=1,kx
                do i=1,ilx
                  sav1(i,(l-1)*kx+k,j)=absnxt(i,k,l,j)
                end do
              end do
            end do
          end do
          allrec = kx*4
          do j=1,jendx
            do l=1,kx+1
              do k=1,kx+1
                do i=1,ilx
                  sav1(i,allrec+(l-1)*(kx+1)+k,j)=abstot(i,k,l,j)
                end do
              end do
            end do
          end do
          allrec = allrec+(kx+1)*(kx+1)
          do j=1,jendx
            do k=1,kx+1
              do i=1,ilx
                sav1(i,allrec+k,j)=emstot(i,k,j)
              end do
            end do
          end do
          allrec = allrec+kx+1
          call MPI_GATHER(
     &             sav1(1,1,1),ilx*allrec*jxp,MPI_DOUBLE_PRECISION,
     &            sav_1(1,1,1),ilx*allrec*jxp,MPI_DOUBLE_PRECISION,
     &                                      0,MPI_COMM_WORLD,ierr)
          if(myid.eq.0) then
            do j=1,mjx-1
              do l=1,4
                do k=1,kx
                  do i=1,ilx
                    absnxt_io(i,k,l,j)=sav_1(i,(l-1)*kx+k,j)
                  end do
                end do
              end do
            end do
            allrec = kx*4
            do j=1,mjx-1
              do l=1,kx+1
                do k=1,kx+1
                  do i=1,ilx
                    abstot_io(i,k,l,j)=sav_1(i,allrec+(l-1)*(kx+1)+k,j)
                  end do
                end do
              end do
            end do
            allrec = allrec+(kx+1)*(kx+1)
            do j=1,mjx-1
              do k=1,kx+1
                do i=1,ilx
                  emstot_io(i,k,j)=sav_1(i,allrec+k,j)
                end do
              end do
            end do
          end if
          do j=1,jendx
            do n=1,NNSG
              do i=1,ilx
                sav2(i,       n,j) =  taf2d(n,i,j)
                sav2(i,NNSG  +n,j) = tlef2d(n,i,j)
                sav2(i,NNSG*2+n,j) =  ssw2d(n,i,j)
                sav2(i,NNSG*3+n,j) =  srw2d(n,i,j)
              end do
            end do
            do i=1,ilx
              sav2(i,NNSG*4+1 ,j) =   sol2d(i,j)
              sav2(i,NNSG*4+2 ,j) = solvd2d(i,j)
              sav2(i,NNSG*4+3 ,j) = solvs2d(i,j)
              sav2(i,NNSG*4+4 ,j) =   flw2d(i,j)
            end do
          end do
          allrec = NNSG*4+4
          call MPI_GATHER(
     &             sav2(1,1,1),ilx*allrec*jxp,MPI_DOUBLE_PRECISION,
     &            sav_2(1,1,1),ilx*allrec*jxp,MPI_DOUBLE_PRECISION,
     &                                      0,MPI_COMM_WORLD,ierr)
          if(myid.eq.0) then
            do j=1,mjx-1
              do n=1,NNSG
                do i=1,ilx
                  taf2d_io(n,i,j)  = sav_2(i,       n,j)
                  tlef2d_io(n,i,j) = sav_2(i,NNSG  +n,j)
                  ssw2d_io(n,i,j)  = sav_2(i,NNSG*2+n,j)
                  srw2d_io(n,i,j)  = sav_2(i,NNSG*3+n,j)
                end do
              end do
              do i=1,ilx
                sol2d_io(i,j)   = sav_2(i,NNSG*4+1 ,j)
                solvd2d_io(i,j) = sav_2(i,NNSG*4+2 ,j)
                solvs2d_io(i,j) = sav_2(i,NNSG*4+3 ,j)
                flw2d_io(i,j)   = sav_2(i,NNSG*4+4 ,j)
              end do
            end do
          end if
          do j=1,jendx
            do n=1,NNSG
              do i=1,ilx
                sav2(i,       n,j) =  tgb2d(n,i,j)
                sav2(i,NNSG  +n,j) =  swt2d(n,i,j)
                sav2(i,NNSG*2+n,j) =  scv2d(n,i,j)
                sav2(i,NNSG*3+n,j) = gwet2d(n,i,j)
              end do
            end do
            do i=1,ilx
              sav2(i,NNSG*4+1 ,j) = flwd2d(i,j)
              sav2(i,NNSG*4+2 ,j) =  fsw2d(i,j)
              sav2(i,NNSG*4+3 ,j) = sabv2d(i,j)
              sav2(i,NNSG*4+4 ,j) = sinc2d(i,j)
            end do
          end do
          allrec = NNSG*4+4
          call MPI_GATHER(
     &             sav2(1,1,1),ilx*allrec*jxp,MPI_DOUBLE_PRECISION,
     &            sav_2(1,1,1),ilx*allrec*jxp,MPI_DOUBLE_PRECISION,
     &                                      0,MPI_COMM_WORLD,ierr)
          if(myid.eq.0) then
            do j=1,mjx-1
              do n=1,NNSG
                do i=1,ilx
                  tgb2d_io(n,i,j)  = sav_2(i,       n,j)
                  swt2d_io(n,i,j)  = sav_2(i,NNSG  +n,j)
                  scv2d_io(n,i,j)  = sav_2(i,NNSG*2+n,j)
                  gwet2d_io(n,i,j) = sav_2(i,NNSG*3+n,j)
                end do
              end do
              do i=1,ilx
                flwd2d_io(i,j) = sav_2(i,NNSG*4+1 ,j)
                fsw2d_io(i,j)  = sav_2(i,NNSG*4+2 ,j)
                sabv2d_io(i,j) = sav_2(i,NNSG*4+3 ,j)
                sinc2d_io(i,j) = sav_2(i,NNSG*4+4 ,j)
              end do
            end do
          end if
          do j=1,jendx
            do n=1,NNSG
              do i=1,ilx
                sav2(i,       n,j) =  veg2d1(n,i,j)
                sav2(i,NNSG  +n,j) =   sag2d(n,i,j)
                sav2(i,NNSG*2+n,j) =  sice2d(n,i,j)
                sav2(i,NNSG*3+n,j) =   dew2d(n,i,j)
              end do
            end do
            do i=1,ilx
              sav2(i,NNSG*4+1 ,j) =   pptnc(i,j)
              sav2(i,NNSG*4+2 ,j) =    pptc(i,j)
              sav2(i,NNSG*4+3 ,j) =  prca2d(i,j)
              sav2(i,NNSG*4+4 ,j) = prnca2d(i,j)
            end do
          end do
          allrec = NNSG*4+4
          call MPI_GATHER(
     &             sav2(1,1,1),ilx*allrec*jxp,MPI_DOUBLE_PRECISION,
     &            sav_2(1,1,1),ilx*allrec*jxp,MPI_DOUBLE_PRECISION,
     &                                      0,MPI_COMM_WORLD,ierr)
          if(myid.eq.0) then
            do j=1,mjx-1
              do n=1,NNSG
                do i=1,ilx
                  veg2d1_io(n,i,j) = sav_2(i,       n,j)
                  sag2d_io(n,i,j)  = sav_2(i,NNSG  +n,j)
                  sice2d_io(n,i,j) = sav_2(i,NNSG*2+n,j)
                  dew2d_io(n,i,j)  = sav_2(i,NNSG*3+n,j)
                end do
              end do
              do i=1,ilx
                pptnc_io(i,j)   = sav_2(i,NNSG*4+1 ,j)
                pptc_io(i,j)    = sav_2(i,NNSG*4+2 ,j)
                prca2d_io(i,j)  = sav_2(i,NNSG*4+3 ,j)
                prnca2d_io(i,j) = sav_2(i,NNSG*4+4 ,j)
              end do
            end do
          end if
          do j=1,jendx
            do n=1,NNSG
              do i=1,ilx
                sav2a(i,       n,j) = ircp2d(n,i,j) 
                sav2a(i,NNSG  +n,j) = text2d(n,i,j)
                sav2a(i,NNSG*2+n,j) =  col2d(n,i,j)
                sav2a(i,NNSG*3+n,j) = ocld2d(n,i,j)
                sav2a(i,NNSG*4+n,j) =   tg2d(n,i,j)
              end do
            end do
            do i=1,ilx
              sav2a(i,NNSG*5+1 ,j) = veg2d(i,j)
            end do
          end do
          allrec = NNSG*5+1
          call MPI_GATHER(
     &             sav2a(1,1,1),ilx*allrec*jxp,MPI_DOUBLE_PRECISION,
     &            sav_2a(1,1,1),ilx*allrec*jxp,MPI_DOUBLE_PRECISION,
     &                                      0,MPI_COMM_WORLD,ierr)
          if(myid.eq.0) then
            do j=1,mjx-1
              do n=1,NNSG
                do i=1,ilx
                  ircp2d_io(n,i,j) = sav_2a(i,       n,j)
                  text2d_io(n,i,j) = sav_2a(i,NNSG  +n,j)
                  col2d_io(n,i,j)  = sav_2a(i,NNSG*2+n,j)
                  ocld2d_io(n,i,j) = sav_2a(i,NNSG*3+n,j)
                    tg2d_io(n,i,j) = sav_2a(i,NNSG*4+n,j)
                end do
              end do
              do i=1,ilx
                veg2d_io(i,j) = sav_2a(i,NNSG*5+1 ,j)
              end do
            end do
          end if
            
          if(ichem.eq.1) then
            do j=1,jendl
              do n=1,ntr
                do k=1,kx
                  do i=1,ix
                    sav4(i,         (n-1)*kx+k,j)=  chia(i,k,j,n)
                    sav4(i,ntr*kx+  (n-1)*kx+k,j)=  chib(i,k,j,n)
                    sav4(i,ntr*kx*2+(n-1)*kx+k,j)=remlsc(i,k,j,n)
                    sav4(i,ntr*kx*3+(n-1)*kx+k,j)=remcvc(i,k,j,n)
                  end do
                end do
              end do
            end do
            allrec = 4*ntr*kx
            do j=1,jendl
              do n=1,ntr
                do i=1,ix
                  sav4(i,allrec+n,j)=remdrd(i,j,n)
                end do
              end do
            end do
            allrec = allrec+ntr
            call MPI_GATHER(
     &               sav4(1,1,1),ix*allrec*jxp,MPI_DOUBLE_PRECISION,
     &              sav_4(1,1,1),ix*allrec*jxp,MPI_DOUBLE_PRECISION,
     &                                      0,MPI_COMM_WORLD,ierr)
            if(myid.eq.0) then
              do j=1,mjx
                do n=1,ntr
                do k=1,kx
                  do i=1,ix
                chia_io(i,k,j,n)  =sav_4(i,         (n-1)*kx+k,j)
                chib_io(i,k,j,n)  =sav_4(i,ntr*kx+  (n-1)*kx+k,j)
                remlsc_io(i,k,j,n)=sav_4(i,ntr*kx*2+(n-1)*kx+k,j)
                remcvc_io(i,k,j,n)=sav_4(i,ntr*kx*3+(n-1)*kx+k,j)
                  end do
                end do
                end do
              end do
              allrec = 4*ntr*kx
              do j=1,mjx
                do n=1,ntr
                  do i=1,ix
                    remdrd_io(i,j,n)=sav_4(i,allrec+n,j)
                  end do
                end do
              end do
            end if
            do j=1,jendx
              do i=1,ilx
                sav4a(i,1,j) = ssw2da(i,j)
                sav4a(i,2,j) = sdeltk2d(i,j)
                sav4a(i,3,j) = sdelqk2d(i,j)
                sav4a(i,4,j) = sfracv2d(i,j)
                sav4a(i,5,j) = sfracb2d(i,j)
                sav4a(i,6,j) = sfracs2d(i,j)
                sav4a(i,7,j) = svegfrac2d(i,j)
              end do
            end do
            call MPI_GATHER(
     &             sav4a(1,1,1),ilx*7*jxp,MPI_DOUBLE_PRECISION,
     &            sav_4a(1,1,1),ilx*7*jxp,MPI_DOUBLE_PRECISION,
     &                                  0,MPI_COMM_WORLD,ierr)
            if(myid.eq.0) then
              do j=1,mjx-1
                do i=1,ilx
                  ssw2da_io(i,j)     = sav_4a(i,1,j)
                  sdeltk2d_io(i,j)   = sav_4a(i,2,j)
                  sdelqk2d_io(i,j)   = sav_4a(i,3,j)
                  sfracv2d_io(i,j)   = sav_4a(i,4,j)
                  sfracb2d_io(i,j)   = sav_4a(i,5,j)
                  sfracs2d_io(i,j)   = sav_4a(i,6,j)
                  svegfrac2d_io(i,j) = sav_4a(i,7,j)
                end do
              end do
            end if
          end if
          do j=1,jendl
            do n=1,nsplit
              do i=1,ix
                sav0d(i,n       ,j) =  dstor(i,j,n)
                sav0d(i,n+nsplit,j) =  hstor(i,j,n)
              end do
            end do
          end do
          call MPI_GATHER(sav0d(1,1,1),ix*nsplit*2*jxp,
     &             MPI_DOUBLE_PRECISION,sav_0d(1,1,1),ix*nsplit*2*jxp,
     &             MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
          if(myid.eq.0) then
            do j=1,mjx
              do n=1,nsplit
                do i=1,ix
                  dstor_io(i,j,n) =  sav_0d(i,n       ,j)
                  hstor_io(i,j,n) =  sav_0d(i,n+nsplit,j)
                end do
              end do
            end do
          end if
          do j=1,jendl
            do k=1,kx
              sav6(k,1,j) = ui1(k,j)
              sav6(k,2,j) = ui2(k,j)
              sav6(k,3,j) =uilx(k,j)
              sav6(k,4,j) = uil(k,j)
              sav6(k,5,j) = vi1(k,j)
              sav6(k,6,j) = vi2(k,j)
              sav6(k,7,j) =vilx(k,j)
              sav6(k,8,j) = vil(k,j)
            end do
          end do
          call MPI_GATHER(sav6(1,1,1),kx*8*jxp,MPI_DOUBLE_PRECISION,
     &                   sav_6(1,1,1),kx*8*jxp,MPI_DOUBLE_PRECISION,
     &                                       0,MPI_COMM_WORLD,ierr)
          if(myid.eq.0) then
            do j=1,mjx
            do k=1,kx
              ui1_io(k,j) = sav_6(k,1,j)
              ui2_io(k,j) = sav_6(k,2,j)
              uilx_io(k,j)= sav_6(k,3,j)
              uil_io(k,j) = sav_6(k,4,j)
              vi1_io(k,j) = sav_6(k,5,j)
              vi2_io(k,j) = sav_6(k,6,j)
              vilx_io(k,j)= sav_6(k,7,j)
              vil_io(k,j) = sav_6(k,8,j)
            end do
            end do
          end if
          call MPI_BCAST(ujlx(1,1),ix*kx,MPI_DOUBLE_PRECISION,
     &                           nproc-1,MPI_COMM_WORLD,ierr)
          call MPI_BCAST(ujl(1,1),ix*kx,MPI_DOUBLE_PRECISION,
     &                          nproc-1,MPI_COMM_WORLD,ierr)
          call MPI_BCAST(vjlx(1,1),ix*kx,MPI_DOUBLE_PRECISION,
     &                           nproc-1,MPI_COMM_WORLD,ierr)
          call MPI_BCAST(vjl(1,1),ix*kx,MPI_DOUBLE_PRECISION,
     &                          nproc-1,MPI_COMM_WORLD,ierr)
          if(myid.eq.0) then
#endif
          close(iutsav)
          itype = 'SAV'
          write(newfil,199) itype,idatex
          inquire(file='output/'//newfil,exist=there)
          if(there) isystm=system('/bin/rm -f output/'//newfil)
          open(iutsav,file='output/'//newfil,status='unknown'
     &               ,form='unformatted')
          print*,'OPENING NEW SAV FILE: ',newfil
          call outsav(iutsav)
          print *,'restart written date = ',ldatez+xtime/1440.
          close(iutsav)
#ifdef MPP1
          end if
#endif
        else if (mod(ntime,nsavfrq).eq.0 .and.
     &          (.not.(jyear.eq.jyearr.and.ktau.eq.ktaur))) then
#ifdef MPP1
          if(myid.eq.0) then
#endif
          close(iutsav)
          itype = 'SAV'
 198      format(a3,'TMP.',i10)
          write(tmpfil,198) itype,idatex
          open(iutsav,file=tmpfil,status='unknown',form='unformatted')
          call outsav(iutsav)
          close(iutsav)
          print *,'SAVTMP RESTART WRITTEN: idatex=',idatex,'ktau=',ktau
          if(oldsav(1:3).eq.'SAV') then
            inquire(file=oldsav,exist=there)
            if (there) then
              write(rmsav,197) oldsav
 197          format('/bin/rm -f ',a50)
              isystm=system(rmsav)
              print*,rmsav
            end if
          end if
          oldsav = tmpfil
#ifdef MPP1
          end if
#endif
        end if
      end if
c
c-----printer output:
c
#ifdef MPP1
      if (myid.eq.0) then
#endif
      if (ifprt) then
        if ((jyear.eq.jyear0.and.ktau.eq.0) .or.
     &      mod(ntime,nprtfrq).eq.0) then
#ifdef MPP1
          write(*,*) 'outprt is not available for parallel code'
          idum = iexec
#else
          call outprt(iexec)
#endif
        end if
      end if

      if (lday.eq.1 .and. lhour.eq.0 .and. nint(xtime).eq.0
     &    .and. (.not.(jyear.eq.jyearr.and.ktau.eq.ktaur))
     &    .and. nnnnnn.ne.nnnend) then
        call mkfile
      end if
#ifdef MPP1
      end if
#endif

      return
      end subroutine output
