cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine chrsetc(ichar,len,isetch)
      implicit none
      integer len
c
c****     this routine sets to the first len characters in ichar
c****     to the character variable specified by isetch
c
c****     on input-
c****         len    - number of characters to be set in ichar (integer)
c****         isetch - character to which each character in ichar will
c****                  be set equal to if input is character. (char)
c
c****     on output-
c****         ichar  - len charcters in ichar is set to the character
c****                  specified by either isetch.  (character)
c
      character*(*) ichar
      character*1 isetch
      integer i
c
      do 20 i=1,len
       ichar(i:i) = isetch
20    continue

      return
      end
