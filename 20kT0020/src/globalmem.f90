!    This reads input files and initializes variables
!   By Paolo Lazzari and Nicola Seriani
!   Originally developed in May 2016

module globalmem

       integer :: nstep, nprint, ndiagno

contains

subroutine init_global(filein,CA_dom, CA_state, CA_seed, CA_rule)
      use latticemem
      use nearneighbourmem
      use statemem
      use seedmem
      use rulemem
      use stepmem
      implicit none
      character(100)       :: filein,seedfile, rulefile, neighbourhood 
      logical              :: create_lattice,restart
      character(3)         :: latticetype,latticetyper
      integer              :: idimension,isize
      integer              :: restart_frame, nstep, nprint, ndiagno
      type(domain)         :: CA_dom
      type(seed  )         :: CA_seed
      type(state )         :: CA_state
      type(rule  )         :: CA_rule

      namelist /files/ seedfile, rulefile
      namelist /lattice/ create_lattice, idimension, latticetype, neighbourhood, isize
      namelist /steps/ restart, restart_frame, nstep, nprint, ndiagno

!   Set some default values: 3 dimensions, simple cubic lattice, Moore's neighbourhood 
      CA_dom%latticetype   ='sc'
      CA_dom%neighbourhood = 'vonneumann'
      CA_dom%isize         = 100 

      nstep = 10
      nprint = 10
      ndiagno = 10

!    Read main input, with information on lattice, rule, number of iterations, output options 
      open(11,file=TRIM(filein))

      read(11,nml=files)

        CA_seed%seedfile      = seedfile 
        CA_rule%rulefile      = rulefile 

      read(11,nml=lattice)
        CA_dom%create_lattice = create_lattice
        CA_dom%latticetype    = latticetype
        CA_dom%neighbourhood  = neighbourhood
        CA_dom%isize          = isize

        write(*,*) 'Type of lattice:       ', CA_dom%latticetype
        write(*,*) 'Type of neighbourhood: ', CA_dom%neighbourhood

        if(idimension.ne.CA_dom%D) STOP 'Wrong dimension - check with value at compilation time'

        if ( CA_dom%latticetype .EQ. 'fcc') then
           if (CA_dom%D .ne. 3) STOP 'wrong dimension - fcc valid only for 3D lattice'
        endif 

      read(11,nml=steps)

       CA_step%restart       = restart
       CA_step%restart_frame = restart_frame
       CA_step%nstep         = nstep
       CA_step%nprint        = nprint
       CA_step%ndiagno       = ndiagno

      close(11)

      end subroutine
end module
