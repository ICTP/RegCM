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

      subroutine rdhead(clat,clon,ds,pt,sigf,sigh,sighrev,xplat,xplon,f,&
                      & xmap,dmap,xlat,xlon,zs,zssd,ls,mdate0,&
                      & iin,inhead,idirect)
 
      use mod_regcm_param , only : jx , iy , kz , ibyte
      use mod_postproc_param , only : dtout , dtbat , dtrad , dtche
      implicit none
!
! Dummy arguments
!
      real(4) :: clat , clon , ds , pt , xplat , xplon
      integer :: idirect , iin , mdate0
      character(70) :: inhead
      real(4) , dimension(jx,iy) :: dmap , f , ls , xlat , xlon , xmap ,&
                              &  zs , zssd
      real(4) , dimension(kz+1) :: sigf
      real(4) , dimension(kz) :: sigh , sighrev
      intent (in) iin , inhead
      intent (out) dmap , f , ls , sighrev , xlat , xlon , xmap , zs ,  &
                 & zssd
      intent (inout) clat , clon , ds , idirect , mdate0 , pt , sigf ,  &
                   & sigh , xplat , xplon
!
! Local variables
!
      real(4) :: dtb , dtc , dto , dtr
      integer :: ibltyp , iboudy , icup , ierr , imoist , k , kk , ni , &
               & nj , nk
      character(6) :: proj
      logical :: there
!
      inquire (file=inhead,exist=there)
      if ( .not.there ) then
        print * , 'OUT_HEAD FILE DOES NOT EXIST: ' , inhead
        stop 'SUBROUTINE RDHEAD'
      end if
      open (iin,file=inhead,status='old',form='unformatted',            &
          & recl=jx*iy*ibyte,access='direct')
      read (iin,rec=1,iostat=ierr) mdate0 , ibltyp , icup , imoist ,    &
                                 & iboudy , nj , ni , nk ,              &
                                 & (sigf(k),k=kz+1,1,-1) , ds , pt ,    &
                                 & clat , clon , xplat , xplon , proj , &
                                 & dto , dtb , dtr , dtc , idirect
      print * , 'mdate0,ibltyp,icup,imoist,iboudy,ni,nj,nk,ds='
      print * , mdate0 , ibltyp , icup , imoist , iboudy , ni , nj ,    &
          & nk , ds
      print * , 'sigf='
      print * , sigf
      print * , 'pt,clat,clon,xplat,xplon,proj,dto,dtb,dtr,dtc='
      print * , pt , clat , clon , xplat , xplon , proj , dto , dtb ,   &
          & dtr , dtc
      if ( ni/=iy .or. nj/=jx .or. kz/=nk ) then
        print * , 'Grid Dimensions DO NOT MATCH'
        print * , '  jx=' , jx , 'ix=' , iy , 'kx=' , kz
        print * , '  ni=' , ni , 'nj=' , nj , 'nk=' , nk
        print * , '  Also check ibyte in postproc.param: ibyte= ' ,     &
            & ibyte
        stop 'BAD DIMENSIONS (SUBROUTINE RDHEAD)'
      end if
      if ( dto/=dtout .or. dtb/=dtbat .or. dtr/=dtrad .or. dtc/=dtche ) &
         & then
        print * , 'OUTPUT INTERVALS ARE IMPROPERLY DEFINED'
        print * , 'dto=' , dto , 'dtout=' , dtout
        print * , 'dtb=' , dtb , 'dtbat=' , dtbat
        print * , 'dtr=' , dtr , 'dtrad=' , dtrad
        print * , 'dtc=' , dtc , 'dtche=' , dtche
        stop 'BAD TIME PARAMETERS (SUBROUTINE RDHEAD)'
      end if
      print * , 'Access type= ' , idirect
!     print*,'ZS'
      read (iin,rec=2,iostat=ierr) zs
!     print*,'ZSSD'
      read (iin,rec=3,iostat=ierr) zssd
!     print*,'LS'
      read (iin,rec=4,iostat=ierr) ls
!     print*,'SATBRT'
      read (iin,rec=5,iostat=ierr) ls
!     print*,'XLAT'
      read (iin,rec=6,iostat=ierr) xlat
!     print*,'XLON'
      read (iin,rec=7,iostat=ierr) xlon
!     print*,'XMAP'
      read (iin,rec=8,iostat=ierr) xmap
!     print*,'DMAP'
      read (iin,rec=9,iostat=ierr) dmap
!     print*,'F'
      read (iin,rec=10,iostat=ierr) f
 
      if ( ierr/=0 ) then
        print * , 'END OF FILE REACHED'
        print * , '  Check ibyte in postproc.param: ibyte= ' , ibyte
        stop 'EOF (SUBROUTINE RDHEAD)'
      end if
 
      do k = 1 , kz
        kk = kz - k + 1
        sigh(k) = (sigf(k)+sigf(k+1))/2.
        sighrev(kk) = sigh(k)
      end do
 
      end subroutine rdhead
