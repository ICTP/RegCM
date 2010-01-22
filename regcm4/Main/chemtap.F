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

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine chemtap
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                     c
c     this subroutine writes the model chem                           c
c                                                                     c
c     iutl : is the output unit number for large-domain variables.    c
c                                                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
 
#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
#     include "include/param1.cb"
#     include "include/param2.cb"
#     include "include/param3.cb"
#     include "include/iunits.cb"
#     include "include/date.cb"
#     include "include/cvaria.cb"
#     include "include/bats.cb2"
#     include "include/trachem.cb"

#ifdef MPP1
#     include "include/io.cb"
      real(kind=4)  fchem(mjx-2,ix-2)
#else
#     include "include/main.cb"
#     include "include/mainchem.cb"
      real(kind=4)  fchem(jx-2,ix-2)
#endif

      integer nrcchem   
      COMMON /CHEREC/ nrcchem
      COMMON /CHEOUT/ fchem
      integer i,j,k,itr
      real(kind=8)  xh, mpd

      xh=xtime/60.
      if (iotyp.eq.2) write (iutchem) idatex

c 3 d fields

      do itr =1,ntr 
       do k=kx,1,-1
         do i=1,ix-2 
#ifdef MPP1
         do j=1,mjx-2
            fchem(j,i) = chia_io(i+1,k,j+1,itr)/psa_io(i+1,j+1)
#else
         do j=1,jx-2
            fchem(j,i) = chia(i+1,k,j+1,itr)/psa(i+1,j+1)
#endif
         enddo
         enddo
         if (iotyp.eq.1) then
           nrcchem=nrcchem+1
           write(iutchem,rec=nrcchem) fchem
         elseif (iotyp.eq.2) then
           write(iutchem) fchem
         end if
       enddo
      enddo

c optical properties of the tracer mixing (independant of ntr)     

       do k=kx,1,-1
         do i=1,ix-2 
#ifdef MPP1
         do j=1,mjx-2
            fchem(j,i) = aerext_io(i+1,k,j+1)
#else
         do j=1,jx-2
            fchem(j,i) = aerext(i+1,k,j+1)
#endif
         enddo
         enddo
         if (iotyp.eq.1) then
           nrcchem=nrcchem+1
           write(iutchem,rec=nrcchem) fchem
         elseif (iotyp.eq.2) then
           write(iutchem) fchem
         end if
       enddo

       do k=kx,1,-1
         do i=1,ix-2 
#ifdef MPP1
         do j=1,mjx-2
            fchem(j,i) = aerssa_io(i+1,k,j+1)
#else
         do j=1,jx-2
            fchem(j,i) = aerssa(i+1,k,j+1)
#endif
         enddo
         enddo
         if (iotyp.eq.1) then
           nrcchem=nrcchem+1
           write(iutchem,rec=nrcchem) fchem
         elseif (iotyp.eq.2) then
           write(iutchem) fchem
         end if
       enddo

       do k=kx,1,-1
         do i=1,ix-2 
#ifdef MPP1
         do j=1,mjx-2
            fchem(j,i) = aerasp_io(i+1,k,j+1)
#else
         do j=1,jx-2
            fchem(j,i) = aerasp(i+1,k,j+1)
#endif
         enddo
         enddo
         if (iotyp.eq.1) then
           nrcchem=nrcchem+1
           write(iutchem,rec=nrcchem) fchem
         elseif (iotyp.eq.2) then
           write(iutchem) fchem
         end if
       enddo

c 2d fields
c instantaneous colum burden
      do itr=1,ntr
c     -------  
       do i=1,ix-2 
#ifdef MPP1
       do j=1,mjx-2
         fchem(j,i) = dtrace_io(i+1,j+1,itr)
#else
       do j=1,jx-2
         fchem(j,i) = dtrace(i+1,j+1,itr)
#endif
       enddo
       enddo
       if (iotyp.eq.1) then
        nrcchem=nrcchem+1
        write(iutchem,rec=nrcchem) fchem
       elseif (iotyp.eq.2) then
        write(iutchem) fchem
       end if
c     -------
c averaged deposition rates

c the variable are the cumul between two output time steps.
c in mg.m-2 ( cf tarcbud). We want to get the correspnding average
c rate between two time steps: so multiply by mpd to get a deposition
c /emissiom rate in mg.m-2.day ( consistant with rainfall output)
       mpd = 24. / chemfrq
c CARE here CUMUL
c        mpd = 1.

       do i=1,ix-2 
#ifdef MPP1
       do j=1,mjx-2
         fchem(j,i) = wdlsc_io(i+1,j+1,itr)*mpd
#else
       do j=1,jx-2
         fchem(j,i) = wdlsc(i+1,j+1,itr)*mpd
#endif
       enddo
       enddo
       if (iotyp.eq.1) then
        nrcchem=nrcchem+1
        write(iutchem,rec=nrcchem) fchem
       elseif (iotyp.eq.2) then
        write(iutchem) fchem
       end if
c    --------
       do i=1,ix-2 
#ifdef MPP1
       do j=1,mjx-2
         fchem(j,i) = wdcvc_io(i+1,j+1,itr)*mpd
#else
       do j=1,jx-2
         fchem(j,i) = wdcvc(i+1,j+1,itr)*mpd
#endif
       enddo
       enddo
       if (iotyp.eq.1) then
        nrcchem=nrcchem+1
        write(iutchem,rec=nrcchem) fchem
       elseif (iotyp.eq.2) then
        write(iutchem) fchem
       end if 
c   -------------
       do i=1,ix-2 
#ifdef MPP1
       do j=1,mjx-2
         fchem(j,i) = ddsfc_io(i+1,j+1,itr)*mpd
#else
       do j=1,jx-2
         fchem(j,i) = ddsfc(i+1,j+1,itr)*mpd
#endif
       enddo
       enddo
       if (iotyp.eq.1) then
        nrcchem=nrcchem+1
        write(iutchem,rec=nrcchem) fchem
       elseif (iotyp.eq.2) then
        write(iutchem) fchem
       end if
c----------------------
       do i=1,ix-2 
#ifdef MPP1
       do j=1,mjx-2
         fchem(j,i) = wxsg_io(i+1,j+1,itr)*mpd
#else
       do j=1,jx-2
         fchem(j,i) = wxsg(i+1,j+1,itr)*mpd
#endif
       enddo
       enddo
       if (iotyp.eq.1) then
        nrcchem=nrcchem+1
        write(iutchem,rec=nrcchem) fchem
       elseif (iotyp.eq.2) then
        write(iutchem) fchem
       end if 
c----------------------
       do i=1,ix-2 
#ifdef MPP1
       do j=1,mjx-2
         fchem(j,i) = wxaq_io(i+1,j+1,itr)*mpd
#else
       do j=1,jx-2
         fchem(j,i) = wxaq(i+1,j+1,itr)*mpd
#endif
       enddo
       enddo
       if (iotyp.eq.1) then
        nrcchem=nrcchem+1
        write(iutchem,rec=nrcchem) fchem
       elseif (iotyp.eq.2) then
        write(iutchem) fchem
       end if 
c----------------------
       do i=1,ix-2
#ifdef MPP1
       do j=1,mjx-2
         fchem(j,i) = cemtrac_io(i+1,j+1,itr)*mpd
#else
       do j=1,jx-2
         fchem(j,i) = cemtrac(i+1,j+1,itr)*mpd
#endif
       enddo
       enddo
       if (iotyp.eq.1) then
        nrcchem=nrcchem+1
        write(iutchem,rec=nrcchem) fchem
       elseif (iotyp.eq.2) then
        write(iutchem) fchem
       end if

c IF YOU WANT THE INSTANTANEOUS EMISSION IN THE SOURCE

c       do i=1,ix-2 
c       do j=1,jx-2
c         fchem(j,i) = chemsrc(i+1,j+1,lmonth,itr) * 1.e9
c        unit = microg/m2.s                   
c       enddo
c       enddo
c       if (iotyp.eq.1) then
c        nrcchem=nrcchem+1
c        write(iutchem,rec=nrcchem) fchem
c       elseif (iotyp.eq.2) then
c        write(iutchem) fchem
c       end if
c ------------------
c reinitialisation for avraged dposition rates ( a modifier)
#ifdef MPP1
       do i=1,ix-2 
       do j=1,mjx-2
         do k=1,kx
          remlsc_io(i,k,j,itr)=0.
          remcvc_io(i,k,j,itr)=0.
          rxsg_io(i,k,j,itr)=0.
          rxsaq1_io(i,k,j,itr)=0.
          rxsaq2_io(i,k,j,itr)=0.
         end do 
          cemtr_io(i,j,itr)=0.
          remdrd_io(i,j,itr)=0.

          wdlsc_io(i,j,itr)=0.
          wdcvc_io(i,j,itr)=0.
          ddsfc_io(i,j,itr)=0.
          wxsg_io(i,j,itr)=0.
          wxaq_io(i,j,itr)=0.
          cemtrac_io(i,j,itr) =0.
       end do
       end do
#else
       do i=1,ix
       do j=1,jx
         do k=1,kx
          remlsc(i,k,j,itr)=0.
          remcvc(i,k,j,itr)=0.
          rxsg(i,k,j,itr)=0.
          rxsaq1(i,k,j,itr)=0.
          rxsaq2(i,k,j,itr)=0.
         end do 
          cemtr(i,j,itr)=0.
          remdrd(i,j,itr)=0.

          wdlsc(i,j,itr)=0.
          wdcvc(i,j,itr)=0.
          ddsfc(i,j,itr)=0.
          wxsg(i,j,itr)=0.
          wxaq(i,j,itr)=0.
          cemtrac(i,j,itr) =0.
       end do
       end do
#endif

      enddo 

c----total aerosol TOA radiative forcing ( independant of the number of tracer) 
        
       do i=1,ix-2 
#ifdef MPP1
       do j=1,mjx-2
         fchem(j,i) = aertarf_io(i+1,j+1)
#else
       do j=1,jx-2
         fchem(j,i) = aertarf(i+1,j+1)                   
#endif
       enddo
       enddo
       if (iotyp.eq.1) then
        nrcchem=nrcchem+1
        write(iutchem,rec=nrcchem) fchem
       elseif (iotyp.eq.2) then
        write(iutchem) fchem
       end if


       do i=1,ix-2 
#ifdef MPP1
       do j=1,mjx-2
         fchem(j,i) = aersrrf_io(i+1,j+1)
#else
       do j=1,jx-2
         fchem(j,i) = aersrrf(i+1,j+1)                   
#endif
       enddo
       enddo
       if (iotyp.eq.1) then
        nrcchem=nrcchem+1
        write(iutchem,rec=nrcchem) fchem
       elseif (iotyp.eq.2) then
        write(iutchem) fchem
       end if

c    reset rad diag to 0 ( averaged between output time stepsin aerout)
#ifdef MPP1
       do j=1,mjx-1
       do i=1,ix-1
         aertarf_io(i,j)=0.
         aersrrf_io(i,j)=0.
       enddo
       enddo
#else
       do j=1,jx-1
       do i=1,ix-1
         aertarf(i,j)=0.
         aersrrf(i,j)=0.
       enddo
       enddo
#endif
      print*,'Chem variables written ',idatex

      return
      end
