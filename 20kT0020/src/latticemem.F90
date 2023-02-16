module latticemem

type domain
INTEGER                               :: isize    ! length of N-dimensional cube edge
#include <define_array.h>
INTEGER                               :: S        ! Size of the vectorized domain
LOGICAL                               :: create_lattice! create lattice proporties or read from precomputed file
CHARACTER (LEN=10)                    :: latticetype
CHARACTER (LEN=30)                    :: neighbourhood
LOGICAL                               :: checkseed ! check flag for consistency with seed
LOGICAL                               :: shiftby1 = .FALSE.! Shift by 1 the lattice FCC
INTEGER                               :: nneigh       ! number of near nieghbours
INTEGER, allocatable, dimension(:,:)  :: v1DtoND  ! map of 1D to ND
INTEGER, allocatable, dimension(:,:)  :: neighlist      ! matrix of vectorized nearneighbours
end type domain

type(domain),save                     :: CA_dom1
REAL   , parameter                    :: PI=acos(-1.)
REAL(4), allocatable, dimension(:,:)  :: examap
LOGICAL, allocatable, dimension(:)    :: v_aux

contains

subroutine allocate_dom(CA_dom)

type(domain)                       :: CA_dom
! local variable
integer   :: i,j
integer   :: L,D
INTEGER, allocatable, dimension(:)    :: v        ! auxiliary index vector


D = CA_dom%D
L = CA_dom%isize

#include <allocate_array.h>


SELECT CASE (TRIM(CA_dom%latticetype))

   CASE('sc','hc')
 
       CA_dom%S=CA_dom%isize**CA_dom%D

       allocate(v(CA_dom%S))

       FORALL (i=1:CA_dom%S) v(i) = i

       CA_dom%m=UNPACK(v,CA_dom%m==CA_dom%m,CA_dom%m)! D-dimensional matrix of indexes

   CASE('fcc')


       CA_dom%n = 0

       do i=1,D
          CA_dom%m = 0
          do j=1,(CA_dom%isize/2)*2,2
             CA_dom%m=EOSHIFT(CA_dom%m,SHIFT= 1,BOUNDARY=1,DIM=i)
             CA_dom%m=CSHIFT(CA_dom%m,SHIFT= 1,DIM=i)
          enddo
          CA_dom%n = CA_dom%m + CA_dom%n 
       enddo

       CA_dom%n = mod(CA_dom%n,2)

       if (CA_dom%shiftby1)  CA_dom%n = -CA_dom%n + 1 ! Zero cells became 1 and viceversa

       CA_dom%S=SUM(CA_dom%n)

       allocate(v(CA_dom%S))

       FORALL (i=1:CA_dom%S) v(i) = i

       CA_dom%m = UNPACK(v,CA_dom%n == 1,CA_dom%n) ! D-dimensional matrix of indexes


   CASE DEFAULT

       STOP     'Lattice not understood, select sc, hc or fcc'

END SELECT

! construct map between 1D world and ND world

allocate(CA_dom%v1DtoND(CA_dom%S,D))

do i=1,D
   CA_dom%n = 0
     do j=1,CA_dom%isize
        CA_dom%n=EOSHIFT(CA_dom%n,SHIFT= 1,BOUNDARY=j,DIM=i)
     enddo
   CA_dom%v1DtoND(:,i) =PACK(CA_dom%n,CA_dom%m>0)
enddo


deallocate(v)

end subroutine

subroutine check_lattice_vs_seed(CA_dom,seedfile)
   implicit none
   type(domain)  :: CA_dom
   logical       :: check_seed
   character(30) :: seedfile
! local variables
   integer       :: nseed
   CHARACTER(30) :: origin_type
   integer       :: origin(CA_dom%D)
   real(8)       :: idxr(CA_dom%D)
   character(30) :: aux,aux1
   logical       :: match
   integer       :: k,shifted(CA_dom%D)

   SELECT CASE (TRIM(CA_dom%latticetype))
      CASE('fcc')
          open(12,file=TRIM(seedfile))

          read(12,*) nseed
          read(12,*) aux, aux1
          read(12,*) origin_type, origin(:)

          read(12,*) idxr(:)

          close(12)

          match = .FALSE.
   
          SELECT CASE (origin_type)

                CASE('manual')
                  shifted(:) = int(idxr(:)) + origin(:)

                CASE('center')
                   shifted = int(idxr(:)) + CA_dom%isize/2

          END SELECT

          do k =1, CA_dom%S

              IF(ALL(CA_dom%v1DtoND(k,:).EQ.shifted)) match = .TRUE.

          enddo

          if ( .NOT. match ) then
               write(*,*) 'shift fcc'
              call deallocate_dom(CA_dom)
              CA_dom%shiftby1=.TRUE.
              call allocate_dom(CA_dom)  
          endif

   END SELECT 

end subroutine

subroutine dump_lattice(CA_dom)
implicit none
!local variables
type(domain)   :: CA_dom
integer        :: i,j
character(len=1024)  :: filename
character(len=1024)  :: cube_side
character(len=10)    :: st_suffix ! state file suffix x,xy,xyz,xyzt,nD

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
filename =  TRIM(CA_dom%latticetype)//TRIM(cube_side)//'.'//TRIM(st_suffix)

OPEN(UNIT=333,FILE=TRIM(filename),FORM="FORMATTED",STATUS="REPLACE",ACTION="WRITE")

    write(unit=333,FMT=*) CA_dom%S
    write(unit=333,FMT=*) 

SELECT CASE (CA_dom%latticetype)

     CASE('hc')
         if ( allocated(examap))        deallocate(examap)
         if ( allocated(v_aux))         deallocate(v_aux)
         allocate(examap(CA_dom%S,2))
         allocate(v_aux(CA_dom%S))
         examap = 0.
         v_aux  = .FALSE.
         call create_exa(CA_dom,1,0.,0.)
         do i=1,CA_dom%S
             write(unit=333,FMT=*) 'Au', examap(i,1),examap(i,2), CA_dom%v1DtoND(i,3:CA_dom%D)
         enddo

     CASE DEFAULT
         do i=1,CA_dom%S
             write(unit=333,FMT=*) 'Au', CA_dom%v1DtoND(i,:)
         enddo

END SELECT

CLOSE(UNIT=333)

end subroutine

subroutine deallocate_dom(CA_dom)
implicit none
type(domain)   :: CA_dom

if ( allocated(CA_dom%m))         deallocate(CA_dom%m)
if ( allocated(CA_dom%n))         deallocate(CA_dom%n)
if ( allocated(CA_dom%v1DtoND))   deallocate(CA_dom%v1DtoND)
if ( allocated(CA_dom%neighlist)) deallocate(CA_dom%neighlist)

end subroutine

recursive subroutine create_exa(CA_dom,point0,x0,y0)
implicit none
type(domain)   :: CA_dom
integer    :: i,point0,point
real(4)    :: x0,y0
real(4)    :: x,y

if (point0 .eq. 0) then
   return
elseif (v_aux(point0)) then
   return
else
   v_aux(point0)    = .TRUE.
   examap(point0,1) = x0
   examap(point0,2) = y0
   do i =1,6
        point = CA_dom%neighlist(point0,i)
        x     = x0 + cos(-2*PI/6.*real(i-1) + 5./6. * PI )
        y     = y0 + sin(-2*PI/6.*real(i-1) + 5./6. * PI )
       call create_exa(CA_dom,point,x,y)
   enddo
   do i =7,CA_dom%nneigh
        point = CA_dom%neighlist(point0,i)
        x     = x0 
        y     = y0 
       call create_exa(CA_dom,point,x,y)
   enddo
endif

end subroutine

subroutine save_lattice(CA_dom)
implicit none
!local variables
type(domain)   :: CA_dom
integer        :: i,j
character(len=1024)  :: filename
character(len=1024)  :: cube_side
character(len=10)    :: st_suffix ! state file suffix x,xy,xyz,xyzt,nD

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
filename =  TRIM(CA_dom%latticetype)//TRIM(cube_side)//'.save.'//TRIM(st_suffix)

OPEN(UNIT=333,FILE=TRIM(filename),FORM="FORMATTED",STATUS="REPLACE",ACTION="WRITE")

    write(unit=333,FMT=*) CA_dom%isize
    write(unit=333,FMT=*) CA_dom%S
    write(unit=333,FMT=*) CA_dom%latticetype
    write(unit=333,FMT=*) CA_dom%neighbourhood
    write(unit=333,FMT=*) CA_dom%checkseed
    write(unit=333,FMT=*) CA_dom%shiftby1
    write(unit=333,FMT=*) CA_dom%nneigh
    do i=1,CA_dom%S
        write(unit=333,FMT=*) CA_dom%v1DtoND(i,:)
    enddo
    do i=1,CA_dom%S
        write(unit=333,FMT=*) CA_dom%neighlist(i,:)
    enddo

CLOSE(UNIT=333)

end subroutine

subroutine read_lattice(CA_dom)
implicit none
!local variables
type(domain)   :: CA_dom
integer        :: i,j
character(len=1024)  :: filename
character(len=1024)  :: cube_side
character(len=10)    :: st_suffix ! state file suffix x,xy,xyz,xyzt,nD

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
filename =  TRIM(CA_dom%latticetype)//TRIM(cube_side)//'.save.'//TRIM(st_suffix)

OPEN(UNIT=333,FILE=TRIM(filename),FORM="FORMATTED",STATUS="OLD",ACTION="READ")

    read(unit=333,FMT=*) CA_dom%isize
    read(unit=333,FMT=*) CA_dom%S
    read(unit=333,FMT=*) CA_dom%latticetype
    read(unit=333,FMT=*) CA_dom%neighbourhood
    read(unit=333,FMT=*) CA_dom%checkseed
    read(unit=333,FMT=*) CA_dom%shiftby1
    read(unit=333,FMT=*) CA_dom%nneigh

    allocate(CA_dom%v1DtoND(CA_dom%S, CA_dom%D))

    do i=1,CA_dom%S
        read(unit=333,FMT=*) CA_dom%v1DtoND(i,:)
    enddo
    
    allocate(CA_dom%neighlist(CA_dom%S,CA_dom%nneigh))

    do i=1,CA_dom%S
        read(unit=333,FMT=*) CA_dom%neighlist(i,:)
    enddo

CLOSE(UNIT=333)

end subroutine

end module
