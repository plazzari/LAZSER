#! /usr/bin/python3

#     INTEGER, parameter        :: D=1        ! Spatial dimensions matrix RANK
#     INTEGER, allocatable,dimension(:)   :: m,n

import sys

cDD = int(sys.argv[1])

varLIST=['m','n']

S1 = 'INTEGER                   :: D=' + str(cDD) + '! Spatial dimensions matrix RANK' 
S2 = 'INTEGER, allocatable,dimension('
for i in range(cDD-1):
    S2 = S2+':,'
S2 = S2 + ':)'
S2 = S2 + '   :: ' + varLIST[0]
for var in varLIST[1:]:
    S2 = S2 + ', ' + var
print(S1)
print(S2)
