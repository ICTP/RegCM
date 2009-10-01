
      subroutine reorder(fdp,fsp,jx,ix,NZ)

      implicit none
      integer jx,ix,NZ,i,j,k,ii,jj
      real*4  fsp(jx*NZ,ix*NZ)
      real*4  fdp(NZ*NZ,jx,ix)

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
