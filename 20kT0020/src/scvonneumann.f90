      subroutine scvonneumann(isize, neighlist)
      implicit none
!     implicit real*8(a-h,o-z)
      integer      :: isize
      integer      :: neighlist(1000000, 30)
!
      integer      :: i1,i2,i3
      integer      :: j,j2
      do i1 =2, isize-1
        do i2 =2, isize-1
          do i3 =2, isize-1
              j = (i1-1)*isize
              j = (j+i2-1)*isize+i3
              neighlist(j, 1) =j-1
              neighlist(j, 2) =j+1

              j2 = (i1-1)*isize
              j2 = (j2+i2-2)*isize+i3
              neighlist(j, 3) =j2

              j2 = (i1-1)*isize
              j2 = (j2+i2)*isize+i3  
              neighlist(j, 4) =j2
    
              j2 = (i1-2)*isize
              j2 = (j2+i2-1)*isize+i3 
              neighlist(j, 5) =j2
 
              j2 = i1*isize
              j2 = (j2+i2-1)*isize+i3
              neighlist(j, 6) =j2

          enddo
        enddo
      enddo      

      return
      end
  
