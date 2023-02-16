!    This code performs cellular automata calculations
!   By Paolo Lazzari and Nicola Seriani
!   Originally developed in May 2016

      subroutine update(CA_dom, CA_state, CA_rule)
      use latticemem
      use statemem
      use rulemem
      use globalmem
      implicit none
      type(domain)         :: CA_dom
      type(state )         :: CA_state
      type(rule  )         :: CA_rule
!     Local variables
      integer i,j,k
!!      type(state )         :: CA_state1

       SELECT CASE (CA_rule%ruletype)
          CASE ('outer', 'outer2')

            do i=1,CA_dom%S
!!!  For every cell, I calculate the number of occupied neighbours
             j=0
             do k=1,CA_dom%nneigh
               j = j + CA_state%ipopulation(CA_dom%neighlist(i,k))
             enddo
!!  Apply update rule to site i, store new population in ipop
              CA_state%ipop(i)=ibits(CA_rule%ictilde, 2*j+CA_state%ipopulation(i),1)
!!              write(*,*) j, CA_state%ipopulation(i), CA_state%ipop(i), CA_rule%ictilde     
            enddo
!! Now update real population ipopulation
            do i=1,CA_dom%S
              CA_state%ipopulation(i) = CA_state%ipop(i)
            enddo            

      CASE ('totalistic')

            do i=1,CA_dom%S
!! This is a totalistic rule, so the central site has to be taken into account
             j=CA_state%ipopulation(i)
!!!  For every cell, I add the number of occupied neighbours
             do k=1,CA_dom%nneigh
               j = j + CA_state%ipopulation(CA_dom%neighlist(i,k))
             enddo
!!  Apply update rule to site i, store new population in ipop
              CA_state%ipop(i)=ibits(CA_rule%icrule, j,1)
!!              write(*,*) j, CA_state%ipopulation(i), CA_state%ipop(i), CA_rule%icrule
            enddo
!! Now update real population ipopulation
            do i=1,CA_dom%S
              CA_state%ipopulation(i) = CA_state%ipop(i)
            enddo

      CASE ('generic')
            STOP 'Update routine: Generic rule not yet implemented!'
      CASE DEFAULT
               STOP 'Update routine: Wrong ruletype'
      END SELECT

      end   
