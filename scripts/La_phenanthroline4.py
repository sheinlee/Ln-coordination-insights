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

query_atom6 = substructure_query.add_atom('N')
query_atom7 = substructure_query.add_atom('N')
query_atom8 = substructure_query.add_atom('C')
query_atom9 = substructure_query.add_atom('C')

query_atom10 =substructure_query.add_atom('N')
query_atom11= substructure_query.add_atom('N')
query_atom12= substructure_query.add_atom('C')
query_atom13= substructure_query.add_atom('C')

query_atom14 =substructure_query.add_atom('N')
query_atom15= substructure_query.add_atom('N')
query_atom16= substructure_query.add_atom('C')
query_atom17= substructure_query.add_atom('C')

query_bond1 = substructure_query.add_bond(single_bond, query_atom1, query_atom2)
query_bond2 = substructure_query.add_bond(single_bond, query_atom1, query_atom3)
query_bond3 = substructure_query.add_bond(aromatic_bond, query_atom2, query_atom4)
query_bond4 = substructure_query.add_bond(aromatic_bond, query_atom3, query_atom5)
query_bond5 = substructure_query.add_bond(aromatic_bond, query_atom4, query_atom5)

query_bond6 = substructure_query.add_bond(single_bond, query_atom1, query_atom6)
query_bond7 = substructure_query.add_bond(single_bond, query_atom1, query_atom7)
query_bond8 = substructure_query.add_bond(aromatic_bond, query_atom6, query_atom8)
query_bond9 = substructure_query.add_bond(aromatic_bond, query_atom7, query_atom9)
query_bond10 = substructure_query.add_bond(aromatic_bond, query_atom8, query_atom9)

query_bond11 = substructure_query.add_bond(single_bond, query_atom1, query_atom10)
query_bond12 = substructure_query.add_bond(single_bond, query_atom1, query_atom11)
query_bond13 = substructure_query.add_bond(aromatic_bond, query_atom10, query_atom12)
query_bond14 = substructure_query.add_bond(aromatic_bond, query_atom11, query_atom13)
query_bond15 = substructure_query.add_bond(aromatic_bond, query_atom12, query_atom13)

query_bond16 = substructure_query.add_bond(single_bond, query_atom1, query_atom14)
query_bond17 = substructure_query.add_bond(single_bond, query_atom1, query_atom15)
query_bond18 = substructure_query.add_bond(aromatic_bond, query_atom14, query_atom16)
query_bond19 = substructure_query.add_bond(aromatic_bond, query_atom15, query_atom17)
query_bond20 = substructure_query.add_bond(aromatic_bond, query_atom16, query_atom17)

from ccdc.search import MoleculeSubstructure, SubstructureSearch

substructure_search = SubstructureSearch()

_ = substructure_search.add_substructure(substructure_query)
hits = substructure_search.search()
print(len(hits))

with open(element+'_phenanthroline4.txt','a') as f:
    f.close()

for h in hits:
    print(h.identifier)
    with open(element+'_phenanthroline4.txt','a') as f:
        f.write(h.identifier)
        f.write("\n")

# 去除重复项
fi = open(element+'_phenanthroline4.txt','r')
txt = fi.readlines()
with open(element+'_phenanthroline4_OK.txt','a') as f:
    f.close()
for w in txt:
    fi2 = open(element+'_phenanthroline4_OK.txt','r')
    txt2 = fi2.readlines()
    with open(element+'_phenanthroline4_OK.txt','a') as f:
        if w not in txt2:
            f.write(w)
        else:
            print("去除已重复--> "+w)
        f.close()
    fi2.close()
fi.close()

