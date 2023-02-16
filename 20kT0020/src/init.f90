!    This reads input files and initializes variables
!   By Paolo Lazzari and Nicola Seriani
!   Originally developed in May 2016

      subroutine init(filein, CA_dom, CA_state, CA_seed, CA_rule) 
      use latticemem
      use nearneighbourmem
      use statemem
      use seedmem
      use stepmem
      use rulemem
      use globalmem
      implicit none
      character(100)        :: filein
      type(domain)         :: CA_dom
      type(seed  )         :: CA_seed
      type(state )         :: CA_state
      type(rule  )         :: CA_rule

!    Read main input, with information on lattice, rule, number of iterations, output options 
   
    write(*,*) 'Reading input file'
    flush(6)
   
    call init_global(filein,CA_dom, CA_state, CA_seed, CA_rule)

!   Initailize the mapping 1D <--> ND and the corresponding nearneihgbour vector

    write(*,*) 'Initialize lattice'
    flush(6)

    if (CA_dom%create_lattice) then

       call allocate_dom(CA_dom)

       call check_lattice_vs_seed(CA_dom,CA_seed%seedfile)

       call compute_nearneighbour(CA_dom)

       call save_lattice(CA_dom)

    else
       
       call read_lattice(CA_dom)

    endif
!   Initailize state

    write(*,*) 'Initialize state'
    flush(6)

    call init_state(CA_dom,CA_state)  

!   Read and initialize the seed

    if (CA_step%restart) then

        call read_restart(CA_dom,CA_state,CA_rule,CA_step%restart_frame)

    else

        call init_seed(CA_dom,CA_state,CA_seed)

        call write_restart(CA_dom,CA_state,CA_rule)

    endif

!   call dump_state(CA_dom,CA_state)

!   Read the rule

    write(*,*) 'Load rule'
    flush(6)

    call init_rule(CA_dom,CA_rule)

      return
      end
