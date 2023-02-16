subroutine calc_state_surface(CA_dom,CA_state_in,CA_rule)
use latticemem
use statemem
use rulemem
implicit none
type(domain), intent(in) :: CA_dom
type(state ), intent(in) :: CA_state_in
type(rule  ), intent(in) :: CA_rule
! local
!type(state ),save             :: CA_state_surf
integer  :: i

call init_state(CA_dom,CA_state_surf)

call init_state(CA_dom,CA_state_surf0)

call copy_state(CA_state_in,CA_state_surf)

call copy_state(CA_state_surf,CA_state_surf0)

!call fill_state(CA_dom,CA_state_surf)

call set_surface_BC(CA_dom,CA_state_surf,CA_rule)

do  i=1,2*CA_dom%s

    call copy_state(CA_state_surf,CA_state_surf0)

    call apply_surface_CA ( CA_dom,CA_state_surf,CA_rule )

    if ( compare_state(CA_dom,CA_state_surf,CA_state_surf0) )  exit

enddo

!call dump_state_surface(CA_dom,CA_state_surf,CA_rule)
call dump_state_surface_bin_2D(CA_dom,CA_state_surf,CA_rule)

end subroutine

logical function MOORE_FN(CA_dom_D,move)
implicit none
integer  :: CA_dom_D
integer  :: move(CA_dom_D)

!local variables
integer      :: i,counter
counter = 0
do i=1,CA_dom_D
  if  ( move(i) .NE. 0 ) counter =counter +1 
enddo
if (counter .NE. 1) then
   MOORE_FN = .FALSE.
else
   MOORE_FN = .TRUE.
endif

return
   
end function

subroutine fromAny_2_FN(CA_dom,near_n_idx)
use latticemem
use nearneighbourmem

implicit none
type(domain) :: CA_dom
integer      :: near_n_idx(2*CA_dom%D)
!local variables
integer      :: i,k,tot_filled,counter
INTEGER,DIMENSION(CA_dom%D)      :: move 
integer      :: a, res
logical      :: MOORE_FN

   SELECT CASE (TRIM(CA_dom%latticetype))

      CASE('sc')

           SELECT CASE (TRIM(CA_dom%neighbourhood))

              CASE('vonneumann')

                  do i = 1, 2*CA_dom%D 

                     near_n_idx(i) = i

                  enddo

              CASE('moore')

                  a = 0

                  counter =1

                  DO i=1,CA_dom%nneigh

                     res=CodeBase(REAL(i-1+a,4),3,move,CA_dom%D)

                     if (res .EQ. 0) then
                        a = 1
                        res=CodeBase(REAL(i-1+a,4),3,move,CA_dom%D)
                     endif

!                    write(*,*) 'move-->', move, 'TRUE', MOORE_FN(CA_dom%D,move) 

                     if (MOORE_FN(CA_dom%D,move)) then 

                         near_n_idx(counter) = i

                         counter = counter + 1

                     endif

                  ENDDO

              CASE DEFAULT
                  write(*,*) 'Scheme: '//TRIM(CA_dom%latticetype)//'-'//TRIM(CA_dom%latticetype)//' not implemented'
                  STOP 

           END SELECT

      CASE('hc')

           SELECT CASE (TRIM(CA_dom%neighbourhood))

              CASE('vonneumann')
                  write(*,*) 'Scheme: '//TRIM(CA_dom%latticetype)//'-'//TRIM(CA_dom%latticetype)//' not implemented'
                  STOP 


              CASE DEFAULT
                  write(*,*) 'Scheme: '//TRIM(CA_dom%latticetype)//'-'//TRIM(CA_dom%latticetype)//' not implemented'
                  STOP 

           END SELECT
          
      CASE('fcc')

           SELECT CASE (TRIM(CA_dom%neighbourhood))

              CASE('vonneumann')
                  write(*,*) 'Scheme: '//TRIM(CA_dom%latticetype)//'-'//TRIM(CA_dom%latticetype)//' not implemented'
                  STOP 

              CASE DEFAULT
                  write(*,*) 'Scheme: '//TRIM(CA_dom%latticetype)//'-'//TRIM(CA_dom%latticetype)//' not implemented'
                  STOP 

           END SELECT

      CASE DEFAULT

           write(*,*) 'Lattice: '//TRIM(CA_dom%latticetype)//' not implemented' 
           STOP 
   END SELECT



end subroutine

subroutine apply_surface_CA ( CA_dom,CA_state, CA_rule  )
use latticemem
use statemem
use rulemem
implicit none
type(domain)            :: CA_dom
type(state )            :: CA_state
type(rule  )            :: CA_rule

!local variables

integer  :: i,k,j,tot_filled,counter
integer  :: near_n_idx(2*CA_dom%D)
integer  :: collect_near_n_fn(2*CA_dom%D+1)
integer  :: collect_near_n(CA_dom%nneigh)
integer  :: rule_surface

  SELECT CASE (TRIM(CA_dom%latticetype))

      CASE('sc')
          call fromAny_2_FN(CA_dom,near_n_idx)

          do i=1, CA_dom%s

             collect_near_n_fn(1) = CA_state%ipopulation(i)

             do j=1,2*CA_dom%D!CA_dom%nneigh
   
                k = near_n_idx(j)

                collect_near_n_fn(j+1) = CA_state%ipopulation(CA_dom%neighlist(i,k))
     
             enddo

          CA_state%ipop(i) = rule_surface(CA_dom, collect_near_n_fn)

          enddo

      CASE('hc','fcc')
          do i=1, CA_dom%s

             collect_near_n(1) = CA_state%ipopulation(i)

             do j=1,CA_dom%nneigh

                collect_near_n(j+1) = CA_state%ipopulation(CA_dom%neighlist(i,j))

             enddo

          CA_state%ipop(i) = rule_surface(CA_dom, collect_near_n)

          enddo

      CASE DEFAULT

           write(*,*) 'Lattice: '//TRIM(CA_dom%latticetype)//' not implemented' 
           STOP 

   END SELECT

          
CA_state%ipopulation = CA_state%ipop

end subroutine

integer function rule_surface(CA_dom, collect_near_n)
use latticemem
implicit none
type(domain)            :: CA_dom
!local variables
integer  :: b
logical  :: c,d
integer  :: collect_near_n(2*CA_dom%D+1)
logical  :: IS_value

b = collect_near_n(1)

c = IS_value(CA_dom, collect_near_n,-1) 

d = IS_value(CA_dom, collect_near_n,1) 

   if (b  .EQ. -1)                                    rule_surface = -1

   if ( (b  .EQ.  0) .AND. (.NOT.c) .AND. (.NOT.d) )  rule_surface =  0
   if ( (b  .EQ.  0) .AND.       c  .AND. (.NOT.d) )  rule_surface = -1
   if ( (b  .EQ.  0) .AND. (.NOT.c) .AND.       d  )  rule_surface =  0 
   if ( (b  .EQ.  0) .AND.       c  .AND.       d  )  rule_surface = -1

   if ( (b  .EQ.  1) .AND. (.NOT.c) .AND. (.NOT.d) )  rule_surface =  1
   if ( (b  .EQ.  1) .AND.       c  .AND. (.NOT.d) )  rule_surface =  2
   if ( (b  .EQ.  1) .AND. (.NOT.c) .AND.       d  )  rule_surface =  1
   if ( (b  .EQ.  1) .AND.       c  .AND.       d  )  rule_surface =  2

   if   (b  .EQ.  2)                                  rule_surface =  2

return
end

logical function IS_value(CA_dom, collect_near_n,v)
use latticemem
implicit none
type(domain)            :: CA_dom
!local variables
integer  :: i,v
integer  :: collect_near_n(2*CA_dom%D+1)

do i=2, 2*CA_dom%D+1

   if (collect_near_n(i) .EQ. v) then

      IS_value = .TRUE. 
      return

   endif

enddo

IS_value = .FALSE.

return
end

subroutine set_surface_BC(CA_dom,CA_state,CA_rule)
use latticemem
use statemem
use rulemem
implicit none
type(domain)            :: CA_dom
type(state )            :: CA_state
type(rule  )            :: CA_rule
!local variables
integer                   :: i,j,tot_filled
integer, allocatable, dimension(:)   :: min_coord ,max_coord
integer, allocatable, dimension(:,:) :: collect_coord

allocate(min_coord(CA_dom%D))
allocate(max_coord(CA_dom%D))

tot_filled = 0

do i=1, CA_dom%s

  if ( CA_state%ipopulation(i) .EQ. 1 ) tot_filled = tot_filled + 1 

enddo

allocate(collect_coord( tot_filled,CA_dom%D))

j=1

do i=1, CA_dom%s

  if ( CA_state%ipopulation(i) .EQ. 1 ) then 
     collect_coord(j,:) = CA_dom%v1DtoND(i,:)   
     j = j + 1
  endif

enddo

do i=1,CA_dom%D
   min_coord(i) = MINVAL(collect_coord(:,i))
   max_coord(i) = MAXVAL(collect_coord(:,i))
enddo


do i=1, CA_dom%s

   do j=1,CA_dom%D
  
      if (  ( CA_dom%v1DtoND(i,j) .LT. min_coord(j) ) .OR. &
            ( CA_dom%v1DtoND(i,j) .GT. max_coord(j) ) ) then

          CA_state%ipopulation(i) =-1
          exit

      endif

   enddo

enddo

deallocate(min_coord,max_coord)
deallocate(collect_coord)

end subroutine

subroutine dump_state_surface(CA_dom,CA_state,CA_rule)
use latticemem
use statemem
use rulemem
implicit none
type(domain)   :: CA_dom
type(state )   :: CA_state
type(rule  )   :: CA_rule
!local variables
integer        :: i,counter
logical        :: isOccupied
logical        :: file_exist
character(len=100)  :: filename
character(len=1024)  :: cube_side
character(len=1024)  :: rule_name
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
write (rule_name, "(I15.15)") CA_rule%icrule

filename =  '../OUTPUT/'//TRIM(CA_dom%latticetype)//'_'//TRIM(cube_side)//'_r'//TRIM(rule_name)//'_state_surface.'//TRIM(st_suffix)

inquire(FILE=TRIM(filename),exist=file_exist)
if(file_exist) then
    OPEN(UNIT=333,FILE=TRIM(filename),FORM="FORMATTED",STATUS="OLD",POSITION="APPEND",ACTION="WRITE")
else
    OPEN(UNIT=333,FILE=TRIM(filename),FORM="FORMATTED",STATUS="REPLACE",ACTION="WRITE")
endif

    counter = 0
    do i=1,CA_dom%S

       if  (CA_state%ipopulation(i) .EQ. 2) counter = counter +1
       
    enddo

    write(unit=333,FMT=*) counter
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
             isOccupied = (CA_state%ipopulation(i) .EQ. 2)
             if (isOccupied) write(unit=333,FMT=*) 'Au', examap(i,1),examap(i,2), CA_dom%v1DtoND(i,3:CA_dom%D)
         enddo

     CASE DEFAULT
         do i=1,CA_dom%S
             isOccupied = (CA_state%ipopulation(i) .EQ. 2)
             if (isOccupied) write(unit=333,FMT=*) 'Au', CA_dom%v1DtoND(i,:)
         enddo

END SELECT

CLOSE(UNIT=333)

end subroutine

subroutine dump_state_surface_bin_2D(CA_dom,CA_state,CA_rule)
use latticemem
use statemem
use rulemem
implicit none
type(domain)   :: CA_dom
type(state )   :: CA_state
type(rule  )   :: CA_rule
!local variables
integer        :: i,counter
integer        :: Noccupied
integer,allocatable  :: occupied(:,:)
logical        :: isOccupied
logical        :: file_exist
character(len=100)  :: filename,filename00,filename01
character(len=1024)  :: cube_side
character(len=1024)  :: rule_name
character(len=1024)  :: frame
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
write (frame, "(I6.6)") frame_surf_count
write (rule_name, "(I15.15)") CA_rule%icrule

filename00 = '../OUTPUT/'//TRIM(CA_dom%latticetype)//'_'//TRIM(cube_side)
filename01 = '_r'//TRIM(rule_name)//'_f'//TRIM(frame)//'_state_surface.'//TRIM(st_suffix)
filename   = TRIM(filename00)//TRIM(filename01)

OPEN(UNIT=333,FILE=TRIM(filename),FORM="UNFORMATTED",STATUS="REPLACE",ACTION="WRITE")

    Noccupied = 0
    do i=1,CA_dom%S

       if  (CA_state%ipopulation(i) .EQ. 2) Noccupied = Noccupied +1
       
    enddo

allocate(occupied(Noccupied,2))

counter = 1

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
             isOccupied = (CA_state%ipopulation(i) .EQ. 2)
             if (isOccupied) then
                 occupied(counter,1)= int ( examap(i,1) )
                 occupied(counter,2)= int ( examap(i,2) )
                 counter = counter +1
             endif
         enddo

     CASE DEFAULT

         do i=1,CA_dom%S
             isOccupied = (CA_state%ipopulation(i) .EQ. 2)
             if (isOccupied) then
                 occupied(counter,1)= int( CA_dom%v1DtoND(i,1) )
                 occupied(counter,2)= int( CA_dom%v1DtoND(i,2) )
                 counter = counter +1
             endif
         enddo

END SELECT

write(unit=333) occupied

CLOSE(UNIT=333)

deallocate(occupied)

end subroutine
