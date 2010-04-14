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

      program sigmatop
      use mod_regcm_param , only : jx , iy , kz
      implicit none
!
! PARAMETER definitions
!
      integer , parameter :: npp = 11
!
! Local variables
!
      real(4) , dimension(jx,iy) :: a , fin , xlat , xlon
      real(4) , dimension(jx,iy,86) :: b
      character(6) :: cgtype
      real(4) :: clat , clon , delx , grdfac , plat , plon , xptop
      character(14) :: fn , fz
      real(4) , dimension(jx-2,iy-2) :: fout
      integer :: i , ibigend , idate , igrads , iyy , j , jxx , k ,     &
               & kzz , mrec , nday , ni , nj , nk , nrec , nslice , numb
      real(4) , dimension(kz+1) :: sigf

      real(4) , dimension(jx,iy,kz) :: h , pp , qv , t , u , v
      real(4) , dimension(jx,iy,npp) :: hp , qp , tp , up , vp
      real(4) , dimension(jx,iy) :: ht , ps , satbrt , slp1 , slp2 , tgb
      real(4) , dimension(npp) :: plev
      real(4) , dimension(kz) :: sig
      real(4) , dimension(jx,iy,kz+1) :: w
!
      data fn/'ICBC1993010100'/
      data fz/'IC_P1993010100'/
 
      data plev/1000. , 925. , 850. , 700. , 500. , 400. , 300. , 250. ,&
         & 200. , 150. , 100./
                            ! temperature of lower soil layer
!
      open (61,file='../DOMAIN.INFO',form='unformatted',access='direct',&
          & recl=iy*jx*4)
      read (61,rec=1) iyy , jxx , kzz , delx , clat , clon , plat ,     &
                    & plon , grdfac , cgtype , (sigf(k),k=1,kz+1) ,     &
                    & xptop , igrads , ibigend
      if ( iyy/=iy .or. jxx/=jx .or. kzz/=kz ) then
        write (*,*) 'There is inconsistence among parameters:'
        write (*,*) 'IY,JX,KZ,IYY,JXX,KZZ' , iy , jx , kz , iyy , jxx , &
                  & kzz
        stop
      end if
      do k = 1 , kz
        sig(k) = 0.5*(sigf(k)+sigf(k+1))
      end do
      read (61,rec=2) ht
      read (61,rec=5) xlat
      read (61,rec=6) xlon
      close (61)
      open (10,file=fz,form='unformatted',recl=(jx-2)*(iy-2)*4,         &
           &access='direct')
      nrec = 0
      open (64,file='../'//fn,form='unformatted',access='direct',       &
          & recl=iy*jx*4)
      mrec = 0
      do nday = 1 , 31
        nslice = 4
        if ( nday==1 ) nslice = 5
        do numb = 1 , nslice
          mrec = mrec + 1
          read (64,rec=mrec) idate , ni , nj , nk
          if ( ni/=jx .or. nj/=iy .or. nk/=kz ) then
            write (*,*) 'IDATE,NI,NJ,NK = ' , idate , ni , nj , nk
            stop
          else
            write (*,*) 'IDATE = ' , idate
          end if
          do k = nk , 1 , -1
            mrec = mrec + 1
            read (64,rec=mrec) ((u(i,j,k),i=1,ni),j=1,nj)
          end do
          do k = nk , 1 , -1
            mrec = mrec + 1
            read (64,rec=mrec) ((v(i,j,k),i=1,ni),j=1,nj)
          end do
          do k = nk , 1 , -1
            mrec = mrec + 1
            read (64,rec=mrec) ((t(i,j,k),i=1,ni),j=1,nj)
          end do
          do k = nk , 1 , -1
            mrec = mrec + 1
            read (64,rec=mrec) ((qv(i,j,k),i=1,ni),j=1,nj)
          end do
          mrec = mrec + 1
          read (64,rec=mrec) ps
          mrec = mrec + 1
          read (64,rec=mrec) tgb
 
          do j = 1 , iy
            do i = 1 , jx
              ps(i,j) = ps(i,j)*10.
              tgb(i,j) = t(i,j,kz)
            end do
          end do
!
!         to calculate Heights on sigma surfaces.
          call htsig(t,h,ps,ht,sig,jx,iy,kz)
!
!         to calculate Sea-Level Pressure using
!         1. ERRICO's solution described in height
!         2. a simple formulae
!         3. MM5 method
          call slpres(h,t,ps,ht,tgb,slp1,slp2,sig,jx,iy,kz)
 
!         to interpolate H,U,V,T,Q and QC
!         1. For Heights
          call height(hp,h,t,ps,ht,sig,jx,iy,kz,plev,npp)
!         2. For Zonal and Meridional Winds
          call intlin(up,u,ps,sig,jx,iy,kz,plev,npp)
          call intlin(vp,v,ps,sig,jx,iy,kz,plev,npp)
!         3. For Temperatures
          call intlog(tp,t,ps,sig,jx,iy,kz,plev,npp)
!         4. For Moisture qva & qca
          call humid1(t,qv,ps,sig,jx,iy,kz)
          call intlin(qp,qv,ps,sig,jx,iy,kz,plev,npp)
          call humid2(tp,qp,plev,jx,iy,npp)
          do k = 1 , npp
            do j = 1 , iy - 2
              do i = 1 , jx - 2
                fout(i,j) = hp(i+1,j+1,k)
              end do
            end do
            nrec = nrec + 1
            write (10,rec=nrec) fout
          end do
          do k = 1 , npp
            do j = 1 , iy - 2
              do i = 1 , jx - 2
                fout(i,j) = tp(i+1,j+1,k)
              end do
            end do
            nrec = nrec + 1
            write (10,rec=nrec) fout
          end do
          do k = 1 , npp
            do j = 1 , iy - 2
              do i = 1 , jx - 2
                fout(i,j) = up(i+1,j+1,k)
              end do
            end do
            nrec = nrec + 1
            write (10,rec=nrec) fout
          end do
          do k = 1 , npp
            do j = 1 , iy - 2
              do i = 1 , jx - 2
                fout(i,j) = vp(i+1,j+1,k)
              end do
            end do
            nrec = nrec + 1
            write (10,rec=nrec) fout
          end do
          do k = 1 , npp
            do j = 1 , iy - 2
              do i = 1 , jx - 2
                fout(i,j) = qp(i+1,j+1,k)
              end do
            end do
            nrec = nrec + 1
            write (10,rec=nrec) fout
          end do
          do j = 1 , iy - 2
            do i = 1 , jx - 2
              fout(i,j) = ps(i+1,j+1)
            end do
          end do
          nrec = nrec + 1
          write (10,rec=nrec) fout
          do j = 1 , iy - 2
            do i = 1 , jx - 2
              fout(i,j) = slp1(i+1,j+1)
            end do
          end do
          nrec = nrec + 1
          write (10,rec=nrec) fout
!         do j=1,iy-2
!           do i=1,jx-2
!             fout(i,j) = slp2(i+1,j+1)
!           end do
!         end do
!         nrec = nrec + 1
!         write (10,rec=nrec) fout
        end do
      end do
!
      end program sigmatop
