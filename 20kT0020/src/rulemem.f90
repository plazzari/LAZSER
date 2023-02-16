module rulemem
use latticemem
type rule
CHARACTER (LEN=100)                   :: rulefile
CHARACTER (LEN=10)                    :: latticetype
CHARACTER (LEN=30)                    :: ruletype
INTEGER                               :: ictilde         ! number of cells for seed
INTEGER                               :: icrule         ! number of cells for seed
INTEGER, allocatable, dimension(:,:)  :: irule           ! indexes for seed
end type rule

type(rule),save                       :: CA_rule1

contains

subroutine init_rule(CA_dom,CA_rule)
   implicit none
   type(domain) :: CA_dom
   type(rule )  :: CA_rule
! local variables
   integer      :: i
      open(12,file=CA_rule%rulefile)

      read(12,*) CA_rule%latticetype

      read(12,*) CA_rule%ruletype

      write(*,*) 'rule file: ', CA_rule%rulefile
      write(*,*) 'From rule file: ',CA_rule%latticetype,'  ', CA_rule%ruletype

      if(CA_rule%latticetype.ne.CA_dom%latticetype) STOP 'Lattice from rulefile different from input lattice'

!   Now use this to initialize the array...
      SELECT CASE (CA_rule%ruletype)
          CASE ('outer')
               read(12,*) CA_rule%ictilde
               CA_rule%icrule = CA_rule%ictilde
          CASE('outer2') 
               allocate(CA_rule%irule(2,0:CA_dom%nneigh))
!  irule(1,:): when the central site is empty
               read(12,*) CA_rule%irule(1,0:CA_dom%nneigh)
!  irule(2,:): when the central site is full
               read(12,*) CA_rule%irule(2,0:CA_dom%nneigh)
               CA_rule%ictilde = 0
!   This is the definition of C tilde from Packard and Wolfram, J Stat Phys 38,
!   901 (1985))
               do i=0,CA_dom%nneigh
                   CA_rule%ictilde = CA_rule%ictilde + CA_rule%irule(1,i)*(2**(2*i))
                   CA_rule%ictilde = CA_rule%ictilde + CA_rule%irule(2,i)*(2**(2*i+1))
               enddo
               write(*,*) ' Ctilde  ', CA_rule%ictilde
               CA_rule%icrule = CA_rule%ictilde

          CASE ('totalistic')
               read(12,*) CA_rule%icrule
!               STOP 'Non-outer totalistic rule not yet implemented!'
          CASE ('alltotalistic')
                write(*,*) ' Doing all totalistic rules  '
          CASE ('allgrowth')
! All outer totalistic rules that correspond to particle growth without dissolution 
               write(*,*) ' Doing all outer totalistic rules leading to growth  '
          CASE ('generic')
               STOP 'Generic rule not yet implemented!'

          CASE DEFAULT
               STOP 'Wrong ruletype'

    END SELECT
!   End of rule initialization
    close(12)
end subroutine

end module
