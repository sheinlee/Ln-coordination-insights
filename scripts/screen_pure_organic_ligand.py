from ccdc import io 
import os
from ccdc.io import MoleculeReader
import numpy as np
import pandas as pd
from itertools import compress
from time import process_time

start=process_time()
LIST_OF_ELEMENT = ['La'
# ,'Ce','Pr','Nd','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu'
]

# 加入Tb 后会导致with open 报错，尚不知原因
# LIST_OF_ELEMENT = ['Tb']

writerA = pd.ExcelWriter('organic_ligand_percentage.xlsx')

for ElEMENT in LIST_OF_ELEMENT:
        A = np.empty(shape=(7000,30),dtype=object)
        COUNT = 0
        countno3d = 0

        filepath='refcode_'+ElEMENT+'.txt'
        csd_reader=io.EntryReader('csd')

        mol_reader = io.MoleculeReader(filepath, format='identifiers')
        for mol in mol_reader:
                COUNT += 1
                A[COUNT-1,0]=COUNT
                
                print(COUNT)
                molecule_name=mol.identifier
                #print(molecule_name)
                mol_name=csd_reader.molecule(molecule_name)
                A[COUNT-1,1]=mol_name.identifier
                
                #判断分子是否包含3D信息
                print(mol_name.is_3d)
                if mol_name.is_3d == False: 
                        countno3d += 1 
                        # print(countno3d)
                        with open(ElEMENT+"no3d.txt","a") as f:
                                f.write(molecule_name)
                                f.write("\n")
                        continue
                #判断ln是否为自由离子
        

                print(mol_name.smiles)
                print(mol_name.atoms)
                # count_for_all=0
                # bond_length_for_all=0

                # A[0,0]=mol_name.identifier
                # len_of_neighbours=len(mol_name.atom('Ce1').neighbours)

                # mol_name.assign_bond_types()
                # mol_name.add_hydrogens()

                len_of_rings=len(mol_name.rings)
                rings = mol_name.rings

                number_of_ln=0
                len_of_neighbours=0
                for atom in mol_name.atoms:
                        if atom.atomic_symbol == ElEMENT:
                                number_of_ln+=1
                                if len(atom.neighbours)>len_of_neighbours:
                                        len_of_neighbours=len(atom.neighbours)
                                # print(len_of_neighbours)
                        
                max_len_of_neighbours=len_of_neighbours
                # position_matrix=np.empty((len_of_neighbours+1,5,number_of_ln),dtype=object)
                connection_matrix=np.zeros((max_len_of_neighbours,max_len_of_neighbours,number_of_ln),dtype=object)
                # print(connection_matrix.shape)
                tag_matrix=np.zeros((max_len_of_neighbours+1,1,number_of_ln),dtype=object)
                
                x=0
                y=0
                z=0
                organic_percent=0
                for atom in mol_name.atoms:
                        if atom.atomic_symbol == ElEMENT:
                                # print('atom.neighbours=',atom.neighbours)
                                print('len(atom.neighbours)=',len(atom.neighbours))
                                #判断ln是否为自由离子
                                if len(atom.neighbours)==0:
                                        continue
                                else:
                                        # faster if narrow down the rings
                                        fil = [r.__contains__(atom) for r in rings]
                                        ring_has_ce = list(compress(rings, fil))
                                        # position_matrix[0,0,z]=mol_name.identifier
                                        # position_matrix[max_len_of_neighbours,1,z]=atom.atomic_number
                                        # position_matrix[max_len_of_neighbours,2,z]=atom.coordinates[0]
                                        # position_matrix[max_len_of_neighbours,3,z]=atom.coordinates[1]
                                        # position_matrix[max_len_of_neighbours,4,z]=atom.coordinates[2]
                                        x=0
                                        for atom1 in atom.neighbours:
                                                # position_matrix[x,1,z]=atom1.atomic_number
                                                # position_matrix[x,2,z]=atom1.coordinates[0]
                                                # position_matrix[x,3,z]=atom1.coordinates[1]
                                                # position_matrix[x,4,z]=atom1.coordinates[2]
                                                y=0
                                                for atom2 in atom.neighbours:
                                                        r=0
                                                        for ring in ring_has_ce:
                                                                
                                                                if atom1 in ring and atom2 in ring:
                                                                        
                                                                        connection_matrix[x,y,z]+=1
                                                        y+=1
                                                x+=1  
                                        tag=2
                                        tag_matrix[max_len_of_neighbours,0,z]=1
                                        tag_matrix[0,0,z]=2
                                        for i in range(0,max_len_of_neighbours):
                                                
                                                for j in range(0,max_len_of_neighbours):
                                                        if connection_matrix[i,j,z]>0:
                                                                if tag_matrix[i,0,z]==0 and tag_matrix[j,0,z]==0:
                                                                        tag+=1
                                                                        tag_matrix[i,0,z]=tag
                                                                        # print(tag)
                                                                        # print(tag_matrix[i,0,z])
                                                                tag_matrix[j,0,z]=tag_matrix[i,0,z]
                                        len_of_neighbours=len(atom.neighbours)
                                        
                                        for i in range(0,len_of_neighbours):
                                                if tag_matrix[i,0,z]==0:
                                                        tag+=1
                                                        tag_matrix[i,0,z]=tag
                                        # print('tag_matrix=\n',tag_matrix)
                                        number_of_organic_ligand=0
                                        number_of_ligand=0
                                        
                                        for i in range(2,10):
                                                ligand_is_organic=False
                                                i_is_tag=False
                                                for j in range(0,len_of_neighbours):
                                                        # print(i,j)
                                                        if tag_matrix[j,0,z]==i:
                                                                i_is_tag=True
                                                                atom3 = atom.neighbours[j]
                                                                # print(atom.neighbours)
                                                                # print('atom3=',atom3)
                                                                if atom3.atomic_symbol=='C':
                                                                        ligand_is_organic=True
                                                                        # print(ligand_is_organic)
                                                                else:
                                                                        for atom4 in atom3.neighbours:
                                                                                # print('atom4=',atom4)
                                                                                if atom4.atomic_symbol=='C':
                                                                                        # print(atom4.atomic_symbol=='C'or atom4.atomic_symbol=='P')
                                                                                        ligand_is_organic=True
                                                                                if atom4.atomic_symbol=='P':
                                                                                        ligand_is_organic=True
                                                                                if atom4.atomic_symbol=='Si':
                                                                                        ligand_is_organic=True
                                                                                        
                                                                if atom3.atomic_symbol=='O' or atom3.atomic_symbol=='N':
                                                                        if len(atom3.neighbours)==1:
                                                                                ligand_is_organic=True
                                                                if atom3.atomic_symbol=='N':
                                                                        for atom4 in atom3.neighbours:
                                                                                if atom4.atomic_symbol=='Si':
                                                                                        ligand_is_organic=True
                                                                if atom3.atomic_symbol=='O':
                                                                        for atom4 in atom3.neighbours:
                                                                                #not co3
                                                                                if atom4.atomic_symbol=='C':
                                                                                        count_for_O=0
                                                                                        for atom5 in atom4.neighbours:
                                                                                                if atom5.atomic_symbol=='O':
                                                                                                        count_for_O+=1
                                                                                        if count_for_O<3:
                                                                                                ligand_is_organic=True
                                                                                        if count_for_O==3:
                                                                                                ligand_is_organic=False        
                                                                                #not So4
                                                                                if atom4.atomic_symbol=='S':
                                                                                        count_for_O=0
                                                                                        for atom5 in atom4.neighbours:
                                                                                                if atom5.atomic_symbol=='O':
                                                                                                        count_for_O+=1
                                                                                        if count_for_O<4:
                                                                                                ligand_is_organic=True
                                                                                #not NO3
                                                                                if atom4.atomic_symbol=='N':
                                                                                        count_for_O=0
                                                                                        for atom5 in atom4.neighbours:
                                                                                                if atom5.atomic_symbol=='O':
                                                                                                        count_for_O+=1
                                                                                        if count_for_O<3:
                                                                                                ligand_is_organic=True
                                                                                if atom4.atomic_symbol=='Si':
                                                                                        ligand_is_organic=True
                                                                                if atom4.atomic_symbol=='B':
                                                                                        ligand_is_organic=True
                                                                                if atom4.atomic_symbol=='Pd':
                                                                                        ligand_is_organic=False
                                                                                if atom4.atomic_symbol=='Ca':
                                                                                        ligand_is_organic=False
                                                                                if atom4.atomic_symbol=='Cl':
                                                                                        ligand_is_organic=False
                                                                                if atom4.atomic_symbol=='Be':
                                                                                        ligand_is_organic=False
                                                                                
                                                                                
                                                                                
                                                                # # S-C     dulplicate
                                                                # if atom3.atomic_symbol=='S':
                                                                #         if len(atom3.neighbours)>1:
                                                                #                 for atom4 in atom3.neighbours:
                                                                #                         if atom4.atomic_symbol=='C':
                                                                #                                 ligand_is_organic=True

                                                                

                                                                if atom3.atomic_symbol=='Mo':
                                                                        ligand_is_organic=False
                                                                if atom3.atomic_symbol=='As':
                                                                        ligand_is_organic=False
                                                                if atom3.atomic_symbol=='Se':
                                                                        for atom4 in atom3.neighbours:
                                                                                if atom4.atomic_symbol=='P':
                                                                                        ligand_is_organic=True
                                                                                else: 
                                                                                        ligand_is_organic=False

                                                                if atom3.atomic_symbol=='H':
                                                                        for atom4 in atom3.neighbours:
                                                                                if atom4.atomic_symbol=='B':
                                                                                        ligand_is_inorganic = False
                                                                                else:
                                                                                        ligand_is_organic=True
                                                                print('ligand_is_organic=',ligand_is_organic)
                                                                
                                                                
                                                if i_is_tag==True:
                                                        number_of_ligand+=1
                                                        
                                                if ligand_is_organic==True:
                                                        number_of_organic_ligand+=1
                                                # print(ligand_is_organic)
                                                
                                                # print('number_of_organic_ligand:',number_of_organic_ligand)
                                                # print('number_of_ligand:',number_of_ligand)
                                        organic_percent=number_of_organic_ligand/number_of_ligand
                                        A[COUNT-1,z+2]=organic_percent          
                                        z+=1     
                        if z>20:
                                break
                
        # writerA = pd.ExcelWriter('in_organic_recognition.xlsx')
        pd.DataFrame(A).to_excel(writerA,sheet_name=ElEMENT,float_format='%.5f')

writerA.save()

#记录程序运行的时间
end=process_time()
print('Running time: %s Seconds'%(end-start))