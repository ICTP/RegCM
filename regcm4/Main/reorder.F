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

      subroutine reorder(fdp,fsp,jx,ix,NZ)

      implicit none
      integer jx,ix,NZ,i,j,k,ii,jj
      real(kind=4)  fsp(jx*NZ,ix*NZ)
      real(kind=4)  fdp(NZ*NZ,jx,ix)

      do j=1,jx*NZ
      do i=1,ix*NZ
        jj=mod(j,NZ)
        if(jj.eq.0) jj=NZ
        ii=mod(i,NZ)
        if(ii.eq.0) ii=NZ
        k=(jj-1)*NZ+ii
        jj=(j+NZ-1)/NZ
        ii=(i+NZ-1)/NZ
        fsp(j,i) = fdp(k,jj,ii)
      end do
      end do
      
      return
      end
