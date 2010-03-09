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

      subroutine cressmcr(b3,b2,alon,alat,glon,glat,jx,iy,i1ur,i1ul,    &
                        & i1dr,i1dl,j1ur,j1ul,j1dr,j1dl,d1xt,d1xa,d1xb, &
                        & d1xc,d1xd,nlon,nlat,nlev)
      use mod_mxncom
      implicit none
!
! Dummy arguments
!
      integer :: iy , jx , nlat , nlev , nlon
      real , dimension(jx,iy) :: alat , alon , d1xa , d1xb , d1xc ,     &
                               & d1xd , d1xt
      real , dimension(nlon,nlat,nlev,5) :: b2
      real , dimension(jx,iy,nlev,5) :: b3
      real , dimension(nlon,nlat) :: glat , glon
      integer , dimension(jx,iy) :: i1dl , i1dr , i1ul , i1ur , j1dl ,  &
                                  & j1dr , j1ul , j1ur
      intent (in) alat , alon , b2 , glat , glon , iy , jx , nlat ,     &
                & nlev , nlon
      intent (out) b3
      intent (inout) d1xa , d1xb , d1xc , d1xd , d1xt , i1dl , i1dr ,   &
                   & i1ul , i1ur , j1dl , j1dr , j1ul , j1ur
!
! Local variables
!
      real :: aaa , dist , dista , distb , distc , distd
      integer :: i , j , k , l , m , mdl , mdr , mul , mur , n , ndl ,  &
               & ndr , nul , nur
!
!     FIND THE FOUR CLOSEST POINTS TO THE GRID WE WANT TO HAVE VALUE,
!     THEN DO THE AVERAGE OF THOSE FOUR POINTS WEIGHTED BY THE DISTANCE.
!     THE CLOSEST ONE HAS BIGGEST WEIGHT.
!
!     B2(JX,IY,NLEV) IS THE INPUT FIELD ON PREVIOUS regular or
!     irregular GRID. B3(JX,IY,NLEV) IS THE OUTPUT FIELD ON new
!     (regular or irregular) GRID. GLON......LONGITUDE VALUES IN
!     DEGREES OF THE INTERMEDIATE GRID4. GLAT......LATITUDE VALUES IN
!     DEGREES OF THE INTERMEDIATE GRID4.
!
!
!     Find the Minimum and Maximum of GLON, GLAT, ALON and ALAT
!
      if ( imxmn==0 ) then
        pi = atan(1.)*4.
        glonmx = -361.
        glonmn = 361.
        do n = 1 , nlat
          do m = 1 , nlon
            if ( glonmx<glon(m,n) ) glonmx = glon(m,n)
            if ( glonmn>glon(m,n) ) glonmn = glon(m,n)
          end do
        end do
        alonmx = -361.
        alonmn = 361.
        do j = 1 , iy
          do i = 1 , jx
            if ( alonmx<alon(i,j) ) alonmx = alon(i,j)
            if ( alonmn>alon(i,j) ) alonmn = alon(i,j)
          end do
        end do
        glatmx = -91.
        glatmn = 91.
        do n = 1 , nlat
          do m = 1 , nlon
            if ( glatmx<glat(m,n) ) glatmx = glat(m,n)
            if ( glatmn>glat(m,n) ) glatmn = glat(m,n)
          end do
        end do
        alatmx = -91.
        alatmn = 91.
        do j = 1 , iy
          do i = 1 , jx
            if ( alatmx<alat(i,j) ) alatmx = alat(i,j)
            if ( alatmn>alat(i,j) ) alatmn = alat(i,j)
          end do
        end do
        write (*,*) 'GLONMN,ALONMN,ALONMX,GLONMX = '
        write (*,*) glonmn , alonmn , alonmx , glonmx
        write (*,*) 'GLATMN,ALATMN,ALATMX,GLATMX = '
        write (*,*) glatmn , alatmn , alatmx , glatmx
        imxmn = 1
      end if
      if ( lcross==0 ) then
        do j = 1 , iy
          do i = 1 , jx
 
            mur = 1000
            nur = 1000
            mul = 1000
            nul = 1000
            mdr = 1000
            ndr = 1000
            mdl = 1000
            ndl = 1000
 
            dista = 1.E8
            distb = 1.E8
            distc = 1.E8
            distd = 1.E8
            do n = 2 , nlat
              do m = 2 , nlon
                if ( (glon(m,n)>=alon(i,j) .and. glon(m,n)-alon(i,j)    &
                   & <10.) .and.                                        &
                   & (glat(m,n)>=alat(i,j) .and. glat(m,n)-alat(i,j)    &
                   & <10.) ) then
                  aaa = ((glon(m,n)-alon(i,j))                          &
                      & *cos((glat(m,n)+alat(i,j))/360.*pi))            &
                      & **2 + (glat(m,n)-alat(i,j))**2
                  if ( dista>aaa ) then
                    dista = aaa
                    mur = m
                    nur = n
                  end if
                end if
                if ( (glon(m,n)<alon(i,j) .and. alon(i,j)-glon(m,n)<10.)&
                   & .and.                                              &
                   & (glat(m,n)>=alat(i,j) .and. glat(m,n)-alat(i,j)    &
                   & <10.) ) then
                  aaa = ((glon(m,n)-alon(i,j))                          &
                      & *cos((glat(m,n)+alat(i,j))/360.*pi))            &
                      & **2 + (glat(m,n)-alat(i,j))**2
                  if ( distb>aaa ) then
                    distb = aaa
                    mul = m
                    nul = n
                  end if
                end if
                if ( (glon(m,n)>=alon(i,j) .and. glon(m,n)-alon(i,j)    &
                   & <10.) .and.                                        &
                   & (glat(m,n)<alat(i,j) .and. alat(i,j)-glat(m,n)<10.)&
                   & ) then
                  aaa = ((glon(m,n)-alon(i,j))                          &
                      & *cos((glat(m,n)+alat(i,j))/360.*pi))            &
                      & **2 + (glat(m,n)-alat(i,j))**2
                  if ( distc>aaa ) then
                    distc = aaa
                    mdr = m
                    ndr = n
                  end if
                end if
                if ( (glon(m,n)<alon(i,j) .and. alon(i,j)-glon(m,n)<10.)&
                   & .and.                                              &
                   & (glat(m,n)<alat(i,j) .and. alat(i,j)-glat(m,n)<10.)&
                   & ) then
                  aaa = ((glon(m,n)-alon(i,j))                          &
                      & *cos((glat(m,n)+alat(i,j))/360.*pi))            &
                      & **2 + (glat(m,n)-alat(i,j))**2
                  if ( distd>aaa ) then
                    distd = aaa
                    mdl = m
                    ndl = n
                  end if
                end if
              end do
            end do
            dist = amin1(dista,distb,distc,distd)
 
            i1ur(i,j) = mur
            j1ur(i,j) = nur
            i1ul(i,j) = mul
            j1ul(i,j) = nul
            i1dr(i,j) = mdr
            j1dr(i,j) = ndr
            i1dl(i,j) = mdl
            j1dl(i,j) = ndl
            d1xt(i,j) = dist
            d1xa(i,j) = dista
            d1xb(i,j) = distb
            d1xc(i,j) = distc
            d1xd(i,j) = distd
            do l = 1 , 4
              do k = 1 , nlev
                if ( dist>0.000001 ) then
                  b3(i,j,k,l) = (b2(mur,nur,k,l)/dista+b2(mul,nul,k,l)  &
                              & /distb+b2(mdr,ndr,k,l)                  &
                              & /distc+b2(mdl,ndl,k,l)/distd)           &
                              & /(1./dista+1./distb+1./distc+1./distd)
                else if ( dist==dista ) then
                  b3(i,j,k,l) = b2(mur,nur,k,l)
                else if ( dist==distb ) then
                  b3(i,j,k,l) = b2(mul,nul,k,l)
                else if ( dist==distc ) then
                  b3(i,j,k,l) = b2(mdr,ndr,k,l)
                else if ( dist==distd ) then
                  b3(i,j,k,l) = b2(mdl,ndl,k,l)
                else
                end if
              end do
            end do
          end do
        end do
        lcross = 1
      else
        do j = 1 , iy
          do i = 1 , jx
 
            mur = i1ur(i,j)
            nur = j1ur(i,j)
            mul = i1ul(i,j)
            nul = j1ul(i,j)
            mdr = i1dr(i,j)
            ndr = j1dr(i,j)
            mdl = i1dl(i,j)
            ndl = j1dl(i,j)
            dist = d1xt(i,j)
            dista = d1xa(i,j)
            distb = d1xb(i,j)
            distc = d1xc(i,j)
            distd = d1xd(i,j)
 
            do l = 1 , 4
              do k = 1 , nlev
                if ( dist>0.000001 ) then
                  b3(i,j,k,l) = (b2(mur,nur,k,l)/dista+b2(mul,nul,k,l)  &
                              & /distb+b2(mdr,ndr,k,l)                  &
                              & /distc+b2(mdl,ndl,k,l)/distd)           &
                              & /(1./dista+1./distb+1./distc+1./distd)
                else if ( dist==dista ) then
                  b3(i,j,k,l) = b2(mur,nur,k,l)
                else if ( dist==distb ) then
                  b3(i,j,k,l) = b2(mul,nul,k,l)
                else if ( dist==distc ) then
                  b3(i,j,k,l) = b2(mdr,ndr,k,l)
                else if ( dist==distd ) then
                  b3(i,j,k,l) = b2(mdl,ndl,k,l)
                else
                end if
              end do
            end do
          end do
        end do
      end if
 
      end subroutine cressmcr
