#! /usr/bin/python3

#     allocate(m(L,L,L,L,L,L)) 
#     allocate(n(L,L,L,L,L,L)) 

import sys

cDD = int(sys.argv[1])

varLIST=['CA_dom%m','CA_dom%n']

size = 'L'

S0='      allocate('
for var in varLIST:
    S1 = S0 + var + '(' + size 
    for i in range(cDD-1):
       S1 = S1 + ',' + size
    S1 = S1 + '))'
    print(S1)
