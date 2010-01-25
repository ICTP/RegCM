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

      subroutine initlk( veg2d, ix1, jx1 )

      implicit none
      integer ix1,jx1
#     include "parameter.inc"

c ******  lkpts = number of lake points desired

      real(kind=8)  t(400), veg2d(ix1,jx1)
      integer ilake, jlake, depth, freeze
      integer lcount,iin,iout,numpts
      common /ldata/ lcount, iin, iout, numpts

      real(kind=8)  hi,hice,hsnow,eta
      data hi, hice, hsnow, eta / .01, 0., 0., .5 /
      integer idata,lakeset,n,j


c ******  unit number containing lake pt locations, depths
      idata = 40

      lcount = 0
      iin = 41
      iout = 42
      numpts = lkpts


c read in vegetation type desired for lake points
c  (18 = mixed woodland (no lake model);  14 = inland water (lake model)

      read(idata,89)  lakeset
89    format(i2)
      print *, '*** lake points set to bats surface type ', lakeset



c initialize data for lake model
      do 20 n=1,lkpts
       read(idata,90)  ilake, jlake, depth
90     format(3i4)
       freeze = 1


c ******     Initially set eta to mean values
       if(depth.lt.50) then
        eta = .7
       else if (depth .gt. 100) then
        eta = .3
       else
          eta = .5
       end if

c ******     Initially set lake points isothermal at 6.0 C (June)
       do 10 j=1,depth
        t(j) = 6.0
10     continue

       veg2d(ilake,jlake) = dble(lakeset)
       write(iin) ilake, jlake, depth, freeze, hi, hice, hsnow
     &            , eta, (t(j),j=1,depth)
       print *, ilake, jlake, depth, freeze

20    continue

      rewind(iin)

      return
      end
