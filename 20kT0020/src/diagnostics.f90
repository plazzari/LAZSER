module diagnomem
use latticemem
use statemem
contains

subroutine dump_diagnostics(CA_dom,CA_state,CA_rule)
use latticemem
use rulemem
implicit none
type(domain)   :: CA_dom
type(state )   :: CA_state
type(rule  )   :: CA_rule
!local variables
integer        :: i, j, k, j1
logical        :: isOccupied
logical        :: file_exist
INTEGER, allocatable, dimension(:)  :: distrib
REAL(8), allocatable, dimension(:)  :: average          ! indexes for seed
REAL(8), allocatable, dimension(:,:)  :: sigma 
character(len=100)  :: filename
character(len=1024)  :: cube_side
character(len=1024)  :: rule_name
character(len=10)    :: st_suffix ! state file suffix x,xy,xyz,xyzt,nD
allocate(distrib(0:CA_dom%nneigh))
allocate(average(1:CA_dom%D))
allocate(sigma(1:CA_dom%D,1:CA_dom%D))
distrib = 0
average =0.d0
sigma =0.d0

SELECT CASE (CA_dom%D)
   CASE(1)
       st_suffix='x'
   CASE(2)
       st_suffix='xy'
   CASE(3)
       st_suffix='xyz'
   CASE(4)
       st_suffix='xyzt'
   CASE DEFAULT
       st_suffix='nD'
END SELECT

write (cube_side, "(I6.6)") CA_dom%isize
write (rule_name, "(I15.15)") CA_rule%icrule

filename =  '../OUTPUT/'//TRIM(CA_dom%latticetype)//'_'//TRIM(cube_side)//'_r'//TRIM(rule_name)//'_diagno.'//TRIM(st_suffix)

inquire(FILE=TRIM(filename),exist=file_exist)
if(file_exist) then
    OPEN(UNIT=333,FILE=TRIM(filename),FORM="FORMATTED",STATUS="OLD",POSITION="APPEND",ACTION="WRITE")
else
    OPEN(UNIT=333,FILE=TRIM(filename),FORM="FORMATTED",STATUS="REPLACE",ACTION="WRITE")
endif

    write(unit=333,FMT=*) sum(CA_state%ipopulation(:))
    write(unit=333,FMT=*) 
    if( sum(CA_state%ipopulation(:)).eq.0) STOP 'empty lattice!'

    do i=1,CA_dom%S
!!!  For every cell, I calculate the number of occupied neighbours
       j=0
       do k=1,CA_dom%nneigh
          j = j + CA_state%ipopulation(CA_dom%neighlist(i,k))
       enddo
    distrib(j) = distrib(j) + CA_state%ipopulation(i) 
    enddo
    write(unit=333,FMT=*) distrib


SELECT CASE (CA_dom%latticetype)

     CASE('hc')
         if ( allocated(examap))        deallocate(examap)
         if ( allocated(v_aux))         deallocate(v_aux)
         allocate(examap(CA_dom%S,CA_dom%D))
         allocate(v_aux(CA_dom%S))
         examap = 0.
         v_aux  = .FALSE.
         call create_exa(CA_dom,1,0.,0.)
         do i=1,CA_dom%S
            do j =3, CA_dom%D
               examap(i,j) = CA_dom%v1DtoND(i,j)          
            enddo
!             isOccupied = (CA_state%ipopulation(i) .EQ. 1)
!             if (isOccupied) write(unit=333,FMT=*) 'Au', examap(i,1),examap(i,2), CA_dom%v1DtoND(i,3:CA_dom%D)
            do j=1, CA_dom%D
              average(j) = average(j) +CA_state%ipopulation(i)*examap(i,j)
              do j1=1,CA_dom%D
                 sigma(j,j1) = sigma(j,j1) + &
         CA_state%ipopulation(i)*examap(i,j)*CA_state%ipopulation(i)*examap(i,j1)  
              enddo 
            enddo

         enddo
  
     CASE DEFAULT
         do i=1,CA_dom%S
!             isOccupied = (CA_state%ipopulation(i) .EQ. 1)
!             if (isOccupied) write(unit=333,FMT=*) 'Au', CA_dom%v1DtoND(i,:)
             do j=1, CA_dom%D
              average(j) = average(j) +CA_state%ipopulation(i)*CA_dom%v1DtoND(i,j)
              do j1=1,CA_dom%D
                 sigma(j,j1) = sigma(j,j1) + &
       CA_state%ipopulation(i)*CA_dom%v1DtoND(i,j)*CA_state%ipopulation(i)*CA_dom%v1DtoND(i,j1)
              enddo
            enddo
         enddo

END SELECT

         average = average/sum(CA_state%ipopulation(:))
         sigma = sigma/sum(CA_state%ipopulation(:))

         do j=1, CA_dom%D
           do j1=1,CA_dom%D
              sigma(j,j1) = sigma(j,j1) - average(j)*average(j1)
           enddo
         enddo
         
         
         write(unit=333,FMT=*) 'average'
         write(unit=333,FMT=*) average
         write(unit=333,FMT=*)
         write(unit=333,FMT=*) 'sigma(x,:)'
         do j=1, CA_dom%D
            write(unit=333,FMT=*) sigma(j,:)
         enddo

CLOSE(UNIT=333)

deallocate(distrib)
deallocate(average)
deallocate(sigma)  

end subroutine
end module
