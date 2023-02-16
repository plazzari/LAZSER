!    This code performs cellular automata calculations
!   By Paolo Lazzari and Nicola Seriani
!   Originally developed in May 2016

      subroutine cellular(CA_dom, CA_state, CA_rule)
      use latticemem
      use statemem
      use rulemem
      use globalmem
      use stepmem 
      use diagnomem
      implicit none
      type(domain)         :: CA_dom
      type(state )         :: CA_state
      type(rule  )         :: CA_rule
!     Local variabiles
      integer i, i2, j


      if(CA_rule%ruletype.eq.'alltotalistic') then

         CA_rule%ruletype = 'totalistic'
         do i2 =0, 2**(CA_dom%nneigh+2)-1
           write(*,*) 'rule: C = ', i2
           CA_rule%icrule = i2
           CA_state%ipopulation = CA_state%ipopulation0
           
           do i =1, CA_step%nstep
             
             call update(CA_dom,CA_state, CA_rule)
!            if(mod(i,CA_step%nprint).eq.0)  call dump_state(CA_dom,CA_state,CA_rule)
             if(mod(i,CA_step%nprint).eq.0)  call dump_state_bin_2D(CA_dom,CA_state,CA_rule)
             if(mod(i,CA_step%ndiagno).eq.0) call dump_diagnostics(CA_dom,CA_state,CA_rule)
             if(mod(i,CA_step%ndiagno).eq.0) call calc_state_surface(CA_dom,CA_state,CA_rule)

             call updtate_counters()
           enddo 
            
           call reset_counters()
         enddo

      elseif(CA_rule%ruletype.eq.'allgrowth') then
!         write(*,*) 'WARNING: All growth not tested!' 

         CA_rule%ruletype = 'outer'

         do i2 =0, 2**(CA_dom%nneigh+1)-1
!!!           write(*,*) 'rule: Ctilde = ', i2

           CA_rule%ictilde = 0   

!! In the growth rules, if the central site is full, it remains full:
!! This is why in the following the odd powers of 2 are always 1 in ictilde
!! (see the definition of C tilde from Packard and Wolfram, J Stat Phys 38,
!   901 (1985))
!! This definition must be kept compatible with ictilde definition in rulemem
!! and with CA dynamics as described in update.f90 for outer totalistic rules

! This condition is designed to avoid appearance of full cells from vacuum, i.e. without full neighbours
           if(ibits(i2, 0 ,1).eq.0) then


             do j=0, CA_dom%nneigh 
            
              CA_rule%ictilde=CA_rule%ictilde + ibits(i2, j ,1)*(2**(2*j))
              CA_rule%ictilde=CA_rule%ictilde + (2**(2*j+1))

             enddo
              write(*,*) 'rule: Ctilde = ', CA_rule%ictilde
              CA_rule%icrule = CA_rule%ictilde

             CA_state%ipopulation = CA_state%ipopulation0

             do i =1, CA_step%nstep

               call update(CA_dom,CA_state, CA_rule)
!              if(mod(i,CA_step%nprint).eq.0)  call dump_state(CA_dom,CA_state,CA_rule)
               if(mod(i,CA_step%nprint).eq.0)  call dump_state_bin_2D(CA_dom,CA_state,CA_rule)
               if(mod(i,CA_step%ndiagno).eq.0)  call dump_diagnostics(CA_dom,CA_state,CA_rule)
               if(mod(i,CA_step%ndiagno).eq.0) call calc_state_surface(CA_dom,CA_state,CA_rule)

             call updtate_counters()
             enddo
           endif
             call reset_counters()
           enddo
      else

        do i =1, CA_step%nstep 

         call update(CA_dom,CA_state, CA_rule)
!        if(mod(i,CA_step%nprint).eq.0)  call dump_state(CA_dom,CA_state,CA_rule)
         if(mod(i,CA_step%nprint).eq.0)  call dump_state_bin_2D(CA_dom,CA_state,CA_rule)
         if(mod(i,CA_step%ndiagno).eq.0)  call dump_diagnostics(CA_dom,CA_state,CA_rule)
         if(mod(i,CA_step%ndiagno).eq.0) call calc_state_surface(CA_dom,CA_state,CA_rule)

         call updtate_counters()

        enddo
      endif

      end
