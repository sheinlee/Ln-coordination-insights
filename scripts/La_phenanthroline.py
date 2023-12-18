#!/usr/bin/env python
# This script is used for Ce/phenanthroline screening 

from ccdc import io 
import os
import numpy as np
from ccdc.molecule import Bond
from ccdc.search import QuerySubstructure

element = 'La'

single_bond = Bond.BondType('Single')
double_bond = Bond.BondType('Double')
aromatic_bond = Bond.BondType('Aromatic')
substructure_query = QuerySubstructure()
query_atom1 = substructure_query.add_atom(element)
query_atom2 = substructure_query.add_atom('N')
query_atom3 = substructure_query.add_atom('N')
query_atom4 = substructure_query.add_atom('C')
query_atom5 = substructure_query.add_atom('C')

query_bond1 = substructure_query.add_bond(single_bond, query_atom1, query_atom2)
query_bond2 = substructure_query.add_bond(single_bond, query_atom1, query_atom3)
query_bond3 = substructure_query.add_bond(aromatic_bond, query_atom2, query_atom4)
query_bond4 = substructure_query.add_bond(aromatic_bond, query_atom3, query_atom5)
query_bond5 = substructure_query.add_bond(aromatic_bond, query_atom4, query_atom5)


from ccdc.search import MoleculeSubstructure, SubstructureSearch

substructure_search = SubstructureSearch()

_ = substructure_search.add_substructure(substructure_query)
hits = substructure_search.search()
print(len(hits))

with open(element+'_phenathroline.txt','a') as f:
    f.close()

for h in hits:
    print(h.identifier)
    with open(element+'_phenanthroline.txt','a') as f:
        f.write(h.identifier)
        f.write("\n")

# 去除重复项
fi = open(element+'_phenanthroline.txt','r')
txt = fi.readlines()
with open(element+'_phenanthroline_OK.txt','a') as f:
    f.close()
for w in txt:
    fi2 = open(element+'_phenanthroline_OK.txt','r')
    txt2 = fi2.readlines()
    with open(element+'_phenanthroline_OK.txt','a') as f:
        if w not in txt2:
            f.write(w)
        else:
            print("去除已重复--> "+w)
        f.close()
fi.close()

