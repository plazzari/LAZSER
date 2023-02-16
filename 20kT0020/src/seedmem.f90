module seedmem
use latticemem
use statemem
type seed
CHARACTER (LEN=30)                    :: seedfile
CHARACTER (LEN=10)                    :: latticetype
CHARACTER (LEN=30)                    :: neighbourhood
CHARACTER (LEN=30)                    :: stype         ! indices or coord.... 
CHARACTER (LEN=30)                    :: origin_type   ! manual or centered
INTEGER                               :: nseed         ! number of cells for seed
INTEGER, allocatable, dimension(:)    :: origin        ! indexes for seed
REAL(8), allocatable, dimension(:,:)  :: idxr          ! indexes for seed
INTEGER, allocatable, dimension(:,:)  :: idx           ! indexes for seed
end type seed

type(seed),save                       :: CA_seed1

contains

subroutine init_seed(CA_dom,CA_state,CA_seed)
   implicit none
   type(domain) :: CA_dom
   type(state ) :: CA_state
   type(seed  ) :: CA_seed
! local variables
   integer      :: i
      open(12,file=CA_seed%seedfile)

      read(12,*) CA_seed%nseed

      allocate(CA_seed%idxr   (CA_seed%nseed,CA_dom%D))
      allocate(CA_seed%idx    (CA_seed%nseed,CA_dom%D))
      allocate(CA_seed%origin (CA_dom%D))

      read(12,*) CA_seed%latticetype, CA_seed%stype
      read(12,*) CA_seed%origin_type, CA_seed%origin(:)

      write(*,*) 'Seed file: ',CA_seed%seedfile
      write(*,*) 'Number of cells in the seed: ',CA_seed%nseed
      write(*,*) 'From seed file: ',CA_seed%latticetype,'  ', CA_seed%stype
      write(*,*) 'Seed origin type: ',CA_seed%origin_type,' Seed origin position ', CA_seed%origin(:)

      if(CA_dom%latticetype.ne.CA_seed%latticetype)  STOP 'Seed lattice different from input lattice.'

      if(CA_seed%stype.eq.'indices') then
          if (CA_dom%D .ne. 3) STOP 'indices seed format valid only with 3D lattice, try with coord.'
      endif

      if(CA_seed%stype.eq.'indices') then
          if (CA_seed%latticetype .ne. 'sc') STOP 'indices seed format valid only with sc lattice, try with coord.'
      endif

      do i =1, CA_seed%nseed
          read(12,*) CA_seed%idxr(i,:)
      enddo

      close(12)

!   Now use this to initialize the array...
      SELECT CASE (CA_seed%stype)
          CASE ('indices')
            CA_seed%idx = int(CA_seed%idxr)
            call seed_in_state_ind(CA_dom,CA_state,CA_seed)

          CASE('coord') 
            CA_seed%idx = int(CA_seed%idxr)
            call seed_in_state_coord(CA_dom,CA_state,CA_seed)

          CASE DEFAULT
               STOP 'Wrong seedtype'

    END SELECT
!   End of seed initialization
end subroutine

subroutine seed_in_state_ind(CA_dom,CA_state,CA_seed)
   implicit none
   type(domain) :: CA_dom
   type(state ) :: CA_state
   type(seed  ) :: CA_seed
! local variables
   integer                          :: i,j,k,n,n1,m
   integer,dimension(CA_seed%nseed) :: checker

   checker(:) = 0 

   do n =1, CA_seed%nseed
      i = CA_seed%idx(n,1) + CA_seed%origin(1)
      j = CA_seed%idx(n,2) + CA_seed%origin(2)
      k = CA_seed%idx(n,3) + CA_seed%origin(3)
      if (match_ind(CA_dom,CA_seed,i,j,k)) then
          m = (k-1)*CA_dom%isize
          m = (m+j-1)*CA_dom%isize+i
          CA_state%ipopulation(m)=CA_state%ipopulation(m)+1
        CA_state%ipopulation0(m)=CA_state%ipopulation0(m)+1
          checker(n) = 1
      endif
   enddo

   if (minval(checker) .EQ. 0) STOP 'seed: point in seed not positioned in lattice'
   if (maxval(CA_state%ipopulation) .GT. 1) STOP 'seed: point redundant definition'

end subroutine

subroutine seed_in_state_coord(CA_dom,CA_state,CA_seed)
   implicit none
   type(domain) :: CA_dom
   type(state ) :: CA_state
   type(seed  ) :: CA_seed
! local variables
   integer      :: i,j
   integer      :: checker(CA_seed%nseed)

   checker(:) = 0 

   do i=1,CA_seed%nseed
      do j=1, CA_dom%s

        if (match_coord(CA_dom,j,CA_seed,i)) then
            CA_state%ipopulation(j) = 1
            CA_state%ipopulation0(j) = 1
            checker(i)            = checker(i) + 1
        endif

      enddo
   enddo

   if (minval(checker) .EQ. 0) STOP 'seed: point in seed not positioned in lattice'
   if (maxval(checker) .GT. 1) STOP 'seed: point redundant definition'
end subroutine

function  match_ind(CA_dom,CA_seed,i,j,k)
   implicit none
   type(domain) :: CA_dom
   type(seed)   :: CA_seed
   integer      :: i,j,k
   logical      :: match_ind
! local variables
   integer      :: shifted
   integer      :: i_s,j_s,k_s

   match_ind = .TRUE.
!  SELECT CASE (CA_seed%origin_type)
!        CASE('manual')
!            i_s =  i + CA_seed%origin(1)
!            j_s =  j + CA_seed%origin(2)
!            k_s =  k + CA_seed%origin(3)
!        CASE('center')
!            i_s = i + CA_dom%isize/2
!            j_s = j + CA_dom%isize/2
!            k_s = k + CA_dom%isize/2
!  END SELECT

   if ( (i.LT.1) .OR. (i.GT.CA_dom%isize)) match_ind = .FALSE.
   if ( (j.LT.1) .OR. (j.GT.CA_dom%isize)) match_ind = .FALSE.
   if ( (k.LT.1) .OR. (k.GT.CA_dom%isize)) match_ind = .FALSE.

end function

function  match_coord(CA_dom,j,CA_seed,i)
   implicit none
   type(domain) :: CA_dom
   type(seed  ) :: CA_seed
   integer      :: i,j
   logical      :: match_coord
! local variables
   integer      :: k,shifted
   
   match_coord = .TRUE.
   do k =1, CA_dom%D
      SELECT CASE (CA_seed%origin_type)
         CASE('manual')
             shifted = CA_seed%idx(i,k) + CA_seed%origin(k)
         CASE('center')
             shifted = CA_seed%idx(i,k) + CA_dom%isize/2
      END SELECT
      if (CA_dom%v1DtoND(j,k) .NE. shifted) match_coord = .FALSE.
   enddo
end function
    
end module

