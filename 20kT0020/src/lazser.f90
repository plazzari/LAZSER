!    This code performs cellular automata calculations
!   By Paolo Lazzari and Nicola Seriani
!   Originally developed in May 2016

      program lazser
      use latticemem             
      use statemem             
      use seedmem             
      use rulemem             
      implicit none
! local variables
      character(100)        :: filein

      filein = 'input.dat' 

      call init(filein, CA_dom1, CA_state1, CA_seed1, CA_rule1)

      call cellular(CA_dom1, CA_state1, CA_rule1)

!     call finalize(CA_dom1, CA_state1, CA_seed1, CA_rule1)

      write(*,*) 'End of calculation. Good bye.'
      end
