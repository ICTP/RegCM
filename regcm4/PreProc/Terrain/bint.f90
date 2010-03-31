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

      function bint(xx,yy,list,iii,jjj,flag)
      implicit none
!
! Dummy arguments
!
      logical :: flag
      integer :: iii , jjj
      real(8) :: xx , yy
      real(8) :: bint
      real(8) , dimension(iii,jjj) :: list
      intent (in) flag , iii , jjj , list , xx , yy
!
! Local variables
!
      real(8) :: a , b , c , d , e , f , g , h , x , y
      integer :: i , j , k , kk , knear , l , ll , lnear , n
      real(8) , dimension(4,4) :: stl
      real(8) , external :: oned
!
!-----bilinear interpolation among four grid values
!
      bint = 0.0
      n = 0
      i = int(xx+0.00001)
      j = int(yy+0.00001)
      x = xx - i
      y = yy - j
      if ( abs(x)>0.0001 .or. abs(y)>0.0001 ) then
        do k = 1 , 4
          kk = i + k - 2
          do l = 1 , 4
            stl(k,l) = 0.
            if ( .not.(flag .and. (l==1)) ) then
              if ( .not.(flag .and. (l==4)) ) then
                if ( .not.(flag .and. (k==1)) ) then
                  if ( .not.(flag .and. (k==4)) ) then
                    ll = j + l - 2
                    if ( kk>=1 .and. kk<=iii ) then
                      if ( ll<=jjj .and. ll>=1 ) then
                        stl(k,l) = list(kk,ll)
                        n = n + 1
                        if ( stl(k,l)<=0.0 ) stl(k,l) = -1.E-20
                      end if
                    end if
                  end if
                end if
              end if
            end if
          end do
        end do
!
!-----find index of closest point to xx,yy.
!
        knear = float(2) + x + 0.5
        lnear = float(2) + y + 0.5
        a = oned(x,stl(1,1),stl(2,1),stl(3,1),stl(4,1))
        b = oned(x,stl(1,2),stl(2,2),stl(3,2),stl(4,2))
        c = oned(x,stl(1,3),stl(2,3),stl(3,3),stl(4,3))
        d = oned(x,stl(1,4),stl(2,4),stl(3,4),stl(4,4))
        bint = oned(y,a,b,c,d)
!
!--------if closest point is ocean, automatically reset terrain to
!--------preserve coastline.
!
        if ( .not.flag .and. stl(knear,lnear)<=0.001 ) bint = -0.00001
        if ( n==16 ) return
        if ( flag .and. n==4 ) return
        e = oned(y,stl(1,1),stl(1,2),stl(1,3),stl(1,4))
        f = oned(y,stl(2,1),stl(2,2),stl(2,3),stl(2,4))
        g = oned(y,stl(3,1),stl(3,2),stl(3,3),stl(3,4))
        h = oned(y,stl(4,1),stl(4,2),stl(4,3),stl(4,4))
        bint = (bint+oned(x,e,f,g,h))/2.
        if ( .not.flag .and. stl(knear,lnear)<=0.001 ) bint = -0.00001
        go to 99999
      end if
      bint = list(i,j)
      return
99999 continue
      end function bint
