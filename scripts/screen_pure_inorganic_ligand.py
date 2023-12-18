#this script is used to screen inorganic ligand in lanthanide complexes subset.  

from ccdc import io 
import numpy as np
import pandas as pd
from itertools import compress
from time import process_time

# def O_donor_is_organic()

start=process_time()
LIST_OF_ELEMENT = [
                    'La',
                    # 'Ce',
                    # 'Pr',
                    # 'Nd',
                    # 'Sm',
                    # 'Eu',
                    # 'Gd',
                    # 'Tb',
                    # 'Dy',
                    # 'Ho',
                    # 'Er',
                    # 'Tm',
                    # 'Yb',
                    # 'Lu'
                    ]


writerA = pd.ExcelWriter('inorganic_ligand_percentage.xlsx')

for ElEMENT in LIST_OF_ELEMENT:
    # percentage_of_organic_ligand=0
    # number_of_organic_ligand=0
    # number_of_inorganic_ligand=0
    A = np.empty(shape=(7000,30),dtype=object)
    COUNT = 0
    countno3d = 0

    filepath='refcode_'+ElEMENT+'.txt'
    csd_reader=io.EntryReader('csd')

    mol_reader = io.MoleculeReader(filepath, format='identifiers')
    for mol in mol_reader:
        # percentage_of_organic_ligand=0
        # number_of_organic_ligand=0
        # number_of_inorganic_ligand=0
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
        
        print(mol_name.smiles)
        print(mol_name.atoms)

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

                        number_of_inorganic_ligand=0
                        # number_of_organic_ligand=0
                        number_of_ligand=0

                        for i in range(2,10):
                            i_is_tag=False
                            ligand_is_inorganic=False
                            ligand_is_organic=False
                            
                            for j in range(0,len_of_neighbours):
                                
                                # print(i,j)
                                if tag_matrix[j,0,z]==i:
                                    i_is_tag=True
                                    atom3 = atom.neighbours[j]  #atom3 is donor atom in first shell
                                    if atom3.atomic_symbol=='H':
                                        for atom4 in atom3.neighbours:
                                            if atom4.atomic_symbol=='B':
                                                ligand_is_inorganic = True
                                            else:
                                                ligand_is_organic=True

                                    if atom3.atomic_symbol=='C':
                                        ligand_is_organic = True

                                    if atom3.atomic_symbol=='N':
                                        if len(atom3.neighbours)==1:
                                            ligand_is_organic=True
                                        if len(atom3.neighbours)>1:
                                            for atom4 in atom3.neighbours:   #atom4 is the atom in second shell
                                                if atom4.atomic_symbol=='C':
                                                    ligand_is_organic=True
                                                if atom4.atomic_symbol=='Si':
                                                    ligand_is_organic=True
                                                if atom4.atomic_symbol=='P':
                                                    ligand_is_organic=True
                                        if ligand_is_organic==False:
                                            ligand_is_inorganic=True

                                    if atom3.atomic_symbol=='O':
                                        if len(atom3.neighbours)==1:
                                            ligand_is_organic=True
                                        if len(atom3.neighbours)>1:
                                            #h2o
                                            number_of_H=0
                                            for atom4 in atom3.neighbours:
                                                
                                                if atom4.atomic_symbol=='H':
                                                    number_of_H+=1
                                            if number_of_H==2:
                                                ligand_is_inorganic=True
                                            
                                            #co3
                                            for atom4 in atom3.neighbours:
                                                if atom4.atomic_symbol=='C':
                                                    number_of_O=0
                                                    for atom5 in atom4.neighbours:
                                                        
                                                        if atom5.atomic_symbol=='O':
                                                            number_of_O+=1
                                                    if number_of_O==3:
                                                        ligand_is_inorganic=True
                                                    if number_of_O<3:
                                                        ligand_is_organic=True
                                                    
                                            #no3 or no4
                                            for atom4 in atom3.neighbours:
                                                if atom4.atomic_symbol=='N':
                                                    number_of_O=0
                                                    for atom5 in atom4.neighbours:
                                                        
                                                        if atom5.atomic_symbol=='O':
                                                            number_of_O+=1
                                                    if number_of_O>=3:
                                                        ligand_is_inorganic=True
                                                    if number_of_O<3:
                                                        ligand_is_organic=True
                                                    
                                            #so4
                                            for atom4 in atom3.neighbours:
                                                if atom4.atomic_symbol=='S':
                                                    number_of_O=0
                                                    for atom5 in atom4.neighbours:
                                                        
                                                        if atom5.atomic_symbol=='O':
                                                            number_of_O+=1
                                                    if number_of_O==4:
                                                        ligand_is_inorganic=True
                                                    if number_of_O<4:
                                                        ligand_is_organic=True
                                            #O2            
                                            for atom4 in atom3.neighbours:
                                                if atom4.atomic_symbol=='O':
                                                    ligand_is_inorganic=True

                                            #W,Se,Mo,V,Pb Ti acid
                                            for atom4 in atom3.neighbours:
                                                if atom4.atomic_symbol=='W':
                                                    ligand_is_inorganic=True
                                                if atom4.atomic_symbol=='Se':
                                                    ligand_is_inorganic=True
                                                if atom4.atomic_symbol=='Mo':
                                                    ligand_is_inorganic=True
                                                if atom4.atomic_symbol=='V':
                                                    ligand_is_inorganic=True
                                                if atom4.atomic_symbol=='Pb':
                                                    ligand_is_inorganic=True
                                                if atom4.atomic_symbol=='Ti':
                                                    ligand_is_inorganic=True
                                                if atom4.atomic_symbol=='Mn':
                                                    ligand_is_inorganic=True
                                                if atom4.atomic_symbol=='As':
                                                    ligand_is_inorganic=True
                                                if atom4.atomic_symbol=='Ni':
                                                    ligand_is_inorganic=True
                                                if atom4.atomic_symbol=='Cr':
                                                    ligand_is_inorganic=True
                                                if atom4.atomic_symbol=='Pd':
                                                    ligand_is_inorganic=True
                                                if atom4.atomic_symbol=='Be':
                                                    ligand_is_inorganic=True
                                                if atom4.atomic_symbol=='Cl':
                                                    ligand_is_inorganic=True
                                                if atom4.atomic_symbol=='B':
                                                    for atom5 in atom4.neighbours:
                                                        if atom5.atomic_symbol=='C':
                                                            ligand_is_organic=True
                                                        else:
                                                            ligand_is_inorganic=True


                                        

                                            # if ligand_is_inorganic==False:
                                            #     ligand_is_organic=True        
                                            
                                    
                                    if atom3.atomic_symbol=='P':
                                        ligand_is_organic = True
                                    
                                    if atom3.atomic_symbol=='S':
                                        # if len(atom3.neighbours)==1:
                                        #     ligand_is_inorganic=True
                                        if len(atom3.neighbours)>1: 
                                            #scn
                                            for atom4 in atom3.neighbours:
                                                if atom4.atomic_symbol=='C':
                                                    number_of_N=0
                                                    number_of_S=0
                                                    for atom5 in atom4.neighbours:
                                                        if atom5.atomic_symbol=='N':
                                                            number_of_N+=1
                                                        if atom5.atomic_symbol=='S':
                                                            number_of_S+=1
                                                    if number_of_S ==1 and number_of_N == 1:
                                                        ligand_is_inorganic=True
                                            
                                            #s-Sb,S,Cu,As,Sn,Ge,
                                            for atom4 in atom3.neighbours:
                                                if atom4.atomic_symbol=='Sb':
                                                    ligand_is_inorganic=True
                                                if atom4.atomic_symbol=='S':
                                                    ligand_is_inorganic=True
                                                if atom4.atomic_symbol=='Cu':
                                                    ligand_is_inorganic=True
                                                if atom4.atomic_symbol=='As':
                                                    ligand_is_inorganic=True
                                                if atom4.atomic_symbol=='Sn':
                                                    ligand_is_inorganic=True
                                                if atom4.atomic_symbol=='Ge':
                                                    ligand_is_inorganic=True
                                    
                                    if atom3.atomic_symbol=='F':
                                        if len(atom3.neighbours)==1:
                                            ligand_is_inorganic=True
                                        if len(atom3.neighbours)>1:
                                            for atom4 in atom3.neighbours:
                                                if atom4.atomic_symbol=='C':
                                                    ligand_is_organic=True
                                        if ligand_is_organic==False:
                                            ligand_is_inorganic=True
                                    
                                    if atom3.atomic_symbol=='Cl':
                                        ligand_is_inorganic=True

                                    if atom3.atomic_symbol=='Br':
                                        ligand_is_inorganic=True

                                    if atom3.atomic_symbol=='I':
                                        ligand_is_inorganic=True
                                    
                                    if atom3.atomic_symbol=='Sb':
                                        ligand_is_inorganic=True
                                    
                                    if atom3.atomic_symbol=='Sn':
                                        ligand_is_inorganic=True

                                    if atom3.atomic_symbol=='Bi':
                                        ligand_is_inorganic=True

                                    if atom3.atomic_symbol=='Pb':
                                        ligand_is_inorganic=True

                                    if atom3.atomic_symbol=='As':
                                        ligand_is_inorganic=True

                                    if atom3.atomic_symbol=='B':
                                        ligand_is_inorganic=True

                                    if atom3.atomic_symbol=='Te':
                                        ligand_is_inorganic=True
                                    
                                    if atom3.atomic_symbol=='Se':
                                        for atom4 in atom3.neighbours:
                                            if atom4.atomic_symbol=='P':
                                                ligand_is_organic=True
                                            else:
                                                ligand_is_inorganic=True
                                        


                            if i_is_tag==True:
                                number_of_ligand+=1              
                            if ligand_is_inorganic==True and ligand_is_organic==False:
                                number_of_inorganic_ligand+=1
                        inorganic_percent=number_of_inorganic_ligand/number_of_ligand
                        A[COUNT-1,z+2]=inorganic_percent          
                        z+=1  
                if z>20:
                        break
    # writerA = pd.ExcelWriter('inorganic_ligand_percentage.xlsx')
    pd.DataFrame(A).to_excel(writerA,sheet_name=ElEMENT,float_format='%.5f')

writerA.save()

#记录程序运行的时间
end=process_time()
print('Running time: %s Seconds'%(end-start))