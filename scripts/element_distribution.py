from ccdc import io 
import os
from ccdc.io import MoleculeReader
import numpy as np
import pandas as pd
from itertools import compress
from time import process_time

start=process_time()
LIST_OF_ELEMENT = ['La','Ce','Pr','Nd','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu']

# 加入Tb 后会导致with open 报错，尚不知原因
# LIST_OF_ELEMENT = ['Tb']

writerA = pd.ExcelWriter('element_distribution.xlsx')

for ElEMENT in LIST_OF_ELEMENT:
    A = np.empty(shape=(7000,21),dtype=object)
    COUNT = 0
    countno3d = 0

    filepath='mono_'+ElEMENT+'.txt'
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

        print(mol_name.identifier)
        print(mol_name.smiles)
        print(mol_name.atoms)

        rings = mol_name.rings

        counth = 0
        countc = 0
        countn = 0
        counto = 0
        counts = 0
        countp = 0
        countf = 0
        countcl = 0
        countbr = 0
        counti = 0

        # countc_is_cp=0
        # countc_is_carbene=0
        # countc_is_carbyne=0

        # countn_sp=0
        # countn_sp2=0
        # countn_sp3=0
        try:
            mol_name.assign_bond_types()
            mol_name.add_hydrogens()
        except:
            print('something wrong')
        else:
            counto_is_amide=0
            counto_is_carboxylate_acid=0
            counto_is_ketone=0
            counto_is_ether=0
            counto_is_ester=0
            counto_is_alkoxide=0
            counto_is_hydroxyl=0
            counto_is_acid_and_water=0  
            for atom in mol_name.atoms:
                if atom.atomic_symbol == ElEMENT:

                    # faster if narrow down the rings
                    fil = [r.__contains__(atom) for r in rings]
                    ring_has_ln = list(compress(rings, fil))

                    # print('len(atom.neighbours)=',len(atom.neighbours))
                    if len(atom.neighbours)==0:
                        continue
                    else:
                        range_all = range(len(atom.neighbours))
                        
                        A[COUNT-1,2]=len(atom.neighbours)
                        # for i in range_all:
                        #     #print(i)
                        #     if atom.neighbours[i].atomic_symbol == 'H': counth +=1
                        #     if atom.neighbours[i].atomic_symbol == 'C': countc +=1
                        #     if atom.neighbours[i].atomic_symbol == 'N': countn +=1
                        #     if atom.neighbours[i].atomic_symbol == 'O': counto +=1
                        #     if atom.neighbours[i].atomic_symbol == 'S': counts +=1
                        #     if atom.neighbours[i].atomic_symbol == 'P': countp +=1
                        #     if atom.neighbours[i].atomic_symbol == 'F': countf +=1
                        #     if atom.neighbours[i].atomic_symbol == 'Cl': countcl +=1
                        #     if atom.neighbours[i].atomic_symbol == 'Br': countbr +=1
                        #     if atom.neighbours[i].atomic_symbol == 'I': counti +=1
                            
                        for atom1 in atom.neighbours:
                            #print(i)
                            if atom1.atomic_symbol == 'H': counth +=1
                            if atom1.atomic_symbol == 'C': countc +=1
                            if atom1.atomic_symbol == 'N': countn +=1  
                            if atom1.atomic_symbol == 'O': 
                                counto +=1
                                print('len(atom1.neighbours)',len(atom1.neighbours))
                                print(atom1.neighbours)
                                if len(atom1.neighbours)==1:
                                    counto_is_acid_and_water+=1
                                if len(atom1.neighbours)>1:
                                    countnoc=0
                                    for atom2 in atom1.neighbours:
                                        if atom2.atomic_symbol!='C'and atom2.atomic_symbol!=ElEMENT:
                                            countnoc+=1
                                    if countnoc>1:
                                        counto_is_acid_and_water+=1
                                #search amide
                                for atom2 in atom1.neighbours:
                                    if atom2.atomic_symbol=='C':
                                        if  len(atom2.neighbours)==3:
                                            for atom3 in atom2.neighbours:
                                                if atom3.atomic_symbol=='N':
                                                    if len(atom3.neighbours)==3:
                                                        counto_is_amide+=1
                                #search carboxylate acid
                                for atom2 in atom1.neighbours:
                                    if atom2.atomic_symbol=='C':
                                        if len(atom2.neighbours)==3:
                                            counto_aroundc=0
                                            countc_aroundc=0
                                            #verify RCO2
                                            for atom3 in atom2.neighbours:
                                                
                                                if atom3.atomic_symbol=='O':
                                                    counto_aroundc+=1
                                                if atom3.atomic_symbol=='C':
                                                    countc_aroundc+=1
                                                    
                                            o_is_ionic=0
                                            if counto_aroundc==2 and countc_aroundc==1:
                                                #VERIFY RCO2-M/H
                                                
                                                for atom3 in atom2.neighbours:                                               
                                                    if atom3.atomic_symbol=='O':
                                                        for atom4 in atom3.neighbours:
                                                            # if atom4.atomic_symbol!=ElEMENT:
                                                            if atom4.atomic_symbol=='H' or atom4.is_metal==True:
                                                                o_is_ionic+=1
                                            if o_is_ionic>=2:
                                                counto_is_carboxylate_acid+=1
                                #search ketone
                                for atom2 in atom1.neighbours:
                                    if atom2.atomic_symbol=='C':
                                        for bond in mol_name.bonds:
                                            if (bond.atoms[0]==atom1 and bond.atoms[1]==atom2) or (bond.atoms[0]==atom2 and bond.atoms[1]==atom1):
                                                print(bond.bond_type)
                                                if bond.bond_type!=1:
                                                    if len(atom2.neighbours)==3:
                                                        counto_aroundc=0
                                                        countc_aroundc=0
                                                        for atom3 in atom2.neighbours:
                                                            if atom3.atomic_symbol=='O':
                                                                counto_aroundc+=1
                                                            if atom3.atomic_symbol=='C':
                                                                countc_aroundc+=1
                                                        if counto_aroundc==1 and countc_aroundc==2:
                                                            counto_is_ketone+=1
                                #search ester
                                for atom2 in atom1.neighbours:
                                    if atom2.atomic_symbol=='C':
                                        if len(atom2.neighbours)==3:
                                            counto_aroundc=0
                                            countc_aroundc=0
                                            for atom3 in atom2.neighbours:
                                                if atom3.atomic_symbol=='O':
                                                    counto_aroundc+=1
                                                if atom3.atomic_symbol=='C':
                                                    countc_aroundc+=1
                                            if counto_aroundc==2 and countc_aroundc==1:
                                                #判定O周围有几个非金属原子
                                                flag1=False
                                                flag2=False

                                                for atom3 in atom2.neighbours:
                                                    count_aroundo=0
                                                    if atom3.atomic_symbol=='O':
                                                        for atom4 in atom3.neighbours:
                                                            if atom4.is_metal==False and atom4.atomic_symbol!='H':
                                                                count_aroundo+=1
                                                    if count_aroundo==1:  flag1=True
                                                    if count_aroundo==2:  flag2=True
                                                if flag1 and flag2: 
                                                    counto_is_ester+=1   
                                if len(atom1.neighbours)==2:

                                    #search alkoxide
                                    for atom2 in atom1.neighbours:

                                        if atom2.atomic_symbol=='C':

                                            for bond in mol_name.bonds:
                                                if (bond.atoms[0]==atom1 and bond.atoms[1]==atom2) or (bond.atoms[0]==atom2 and bond.atoms[1]==atom1):
                                                    print(bond.bond_type)
                                                    if bond.bond_type==1:
                                                        counto_aroundc=0
                                                        for atom3 in atom2.neighbours:
                                                            if atom3.atomic_symbol=='O':
                                                                counto_aroundc+=1
                                                        if counto_aroundc==1:
                                                            counto_is_alkoxide+=1

                                if len(atom1.neighbours)==3:
                                    #search ether
                                    countc_aroundo=0
                                    for atom2 in atom1.neighbours:
                                        if atom2.atomic_symbol=='C':
                                            countc_aroundo+=1
                                    if countc_aroundo==2:
                                        counto_is_ether+=1
                                        
                                    #search hydroxyl
                                    for atom2 in atom1.neighbours:

                                        if atom2.atomic_symbol=='C':

                                            for bond in mol_name.bonds:
                                                if (bond.atoms[0]==atom1 and bond.atoms[1]==atom2) or (bond.atoms[0]==atom2 and bond.atoms[1]==atom1):
                                                    print(bond.bond_type)
                                                    if bond.bond_type==1:
                                                        counto_aroundc=0
                                                        counth_aroundc=0
                                                        for atom3 in atom2.neighbours:
                                                            if atom3.atomic_symbol=='O':
                                                                counto_aroundc+=1
                                                            if atom3.atomic_symbol=='H':
                                                                counth_aroundc+=1
                                                        if counto_aroundc==1 and counth_aroundc==1:
                                                            counto_is_hydroxyl+=1
                                # if len(atom1.neighbours)==2:
                                #     #search amide
                                #     for atom2 in atom1.neighbours:
                                #         if atom2.atomic_symbol=='C':
                                #             if  len(atom2.neighbours)==3:
                                #                 for atom3 in atom2.neighbours:
                                #                     if atom3.atomic_symbol=='N':
                                #                         if len(atom3.neighbours)==3:
                                #                             counto_is_amide+=1
                                    
                                #     #search carboxylate acid
                                #     for atom2 in atom1.neighbours:
                                #         if atom2.atomic_symbol=='C':
                                #             if len(atom2.neighbours)==3:
                                #                 counto_aroundc=0
                                #                 countc_aroundc=0
                                #                 #verify RCO2
                                #                 for atom3 in atom2.neighbours:
                                                    
                                #                     if atom3.atomic_symbol=='O':
                                #                         counto_aroundc+=1
                                #                     if atom3.atomic_symbol=='C':
                                #                         countc_aroundc+=1
                                                        
                                #                 o_is_ionic=0
                                #                 if counto_aroundc==2 and countc_aroundc==1:
                                #                     #VERIFY RCO2-M/H
                                                    
                                #                     for atom3 in atom2.neighbours:                                               
                                #                         if atom3.atomic_symbol=='O':
                                #                             for atom4 in atom3.neighbours:
                                #                                 # if atom4.atomic_symbol!=ElEMENT:
                                #                                 if atom4.atomic_symbol=='H' or atom4.is_metal==True:
                                #                                     o_is_ionic+=1
                                #                 if o_is_ionic==2:
                                #                     counto_is_carboxylate_acid+=1
                                #     #search ketone
                                #     for atom2 in atom1.neighbours:
                                #         if atom2.atomic_symbol=='C':
                                #             for bond in mol_name.bonds:
                                #                 if bond.atoms[0]==atom1 and bond.atoms[1]==atom2:
                                #                     print(bond.bond_type)
                                #                     if bond.bond_type!=1:
                                #                         if len(atom2.neighbours)==3:
                                #                             counto_aroundc=0
                                #                             countc_aroundc=0
                                #                             for atom3 in atom2.neighbours:
                                #                                 if atom3.atomic_symbol=='O':
                                #                                     counto_aroundc+=1
                                #                                 if atom3.atomic_symbol=='C':
                                #                                     countc_aroundc+=1
                                #                             if counto_aroundc==1 and countc_aroundc==2:
                                #                                 counto_is_ketone+=1
                                #     #search ester
                                #     for atom2 in atom1.neighbours:
                                #         if atom2.atomic_symbol=='C':
                                #             if len(atom2.neighbours)==3:
                                #                 counto_aroundc=0
                                #                 countc_aroundc=0
                                #                 for atom3 in atom2.neighbours:
                                #                     if atom3.atomic_symbol=='O':
                                #                         counto_aroundc+=1
                                #                     if atom3.atomic_symbol=='C':
                                #                         countc_aroundc+=1
                                #                 if counto_aroundc==2 and countc_aroundc==1:
                                #                     #判定O周围有几个非金属原子
                                #                     flag1=False
                                #                     flag2=False

                                #                     for atom3 in atom2.neighbours:
                                #                         count_aroundo=0
                                #                         if atom3.atomic_symbol=='O':
                                #                             for atom4 in atom3.neighbours:
                                #                                 if atom4.is_metal==False and atom4.atomic_symbol!='H':
                                #                                     count_aroundo+=1
                                #                         if count_aroundo==1:  flag1=True
                                #                         if count_aroundo==2:  flag2=True
                                #                     if flag1 and flag2: 
                                #                         counto_is_ester+=1                    
                                    
                                    
                                #     #search alkoxide
                                #     for atom2 in atom1.neighbours:

                                #         if atom2.atomic_symbol=='C':

                                #             for bond in mol_name.bonds:
                                #                 if bond.atoms[0]==atom1 and bond.atoms[1]==atom2:
                                #                     print(bond.bond_type)
                                #                     if bond.bond_type==1:
                                #                         counto_aroundc=0
                                #                         for atom3 in atom2.neighbours:
                                #                             if atom3.atomic_symbol=='O':
                                #                                 counto_aroundc+=1
                                #                         if counto_aroundc==1:
                                #                             counto_is_alkoxide+=1
                                #             # print(mol_name.bond(atom1,atom2))
                                #             # try:
                                #             #     mol_name.bond(atom1,atom2).BondType==1
                                #             # except: 
                                #             #     print("Something went wrong")
                                #             # else:  
                                #             #     if mol_name.bond(atom1,atom2).BondType==1:
                                #             #         counto_aroundc=0
                                #             #         for atom3 in atom2.neighbours:
                                #             #             if atom3.atomic_symbol=='O':
                                #             #                 counto_aroundc+=1
                                #             #         if counto_aroundc==1:
                                #             #             counto_is_hydroxyl+=1
                                # if len(atom1.neighbours)==3:
                                #     #search ether
                                #     countc_aroundo=0
                                #     for atom2 in atom1.neighbours:
                                #         if atom2.atomic_symbol=='C':
                                #             countc_aroundo+=1
                                #     if countc_aroundo==2:
                                #         counto_is_ether+=1
                                        
                                #     #search hydroxyl
                                #     for atom2 in atom1.neighbours:

                                #         if atom2.atomic_symbol=='C':

                                #             for bond in mol_name.bonds:
                                #                 if bond.atoms[0]==atom1 and bond.atoms[1]==atom2:
                                #                     print(bond.bond_type)
                                #                     if bond.bond_type==1:
                                #                         counto_aroundc=0
                                #                         counth_aroundc=0
                                #                         for atom3 in atom2.neighbours:
                                #                             if atom3.atomic_symbol=='O':
                                #                                 counto_aroundc+=1
                                #                             if atom3.atomic_symbol=='H':
                                #                                 counth_aroundc+=1
                                #                         if counto_aroundc==1 and counth_aroundc==1:
                                #                             counto_is_hydroxyl+=1



                            if atom1.atomic_symbol == 'S': counts +=1
                            if atom1.atomic_symbol == 'P': countp +=1
                            if atom1.atomic_symbol == 'F': countf +=1
                            if atom1.atomic_symbol == 'Cl': countcl +=1
                            if atom1.atomic_symbol == 'Br': countbr +=1
                            if atom1.atomic_symbol == 'I': counti +=1
                        A[COUNT-1,3]=counth
                        A[COUNT-1,4]=countc
                        A[COUNT-1,5]=countn
                        A[COUNT-1,6]=counto
                        A[COUNT-1,7]=counts
                        A[COUNT-1,8]=countp
                        A[COUNT-1,9]=countf
                        A[COUNT-1,10]=countcl
                        A[COUNT-1,11]=countbr
                        A[COUNT-1,12]=counti

                        A[COUNT-1,13]=counto_is_amide
                        A[COUNT-1,14]=counto_is_carboxylate_acid
                        A[COUNT-1,15]=counto_is_ketone
                        A[COUNT-1,16]=counto_is_ester
                        A[COUNT-1,17]=counto_is_ether
                        A[COUNT-1,18]=counto_is_alkoxide
                        A[COUNT-1,19]=counto_is_hydroxyl
                        A[COUNT-1,20]=counto_is_acid_and_water
                    # break   

                
    # writerA = pd.ExcelWriter('in_organic_recognition.xlsx')
    pd.DataFrame(A).to_excel(writerA,sheet_name=ElEMENT,float_format='%.5f')

writerA.save()

#记录程序运行的时间
end=process_time()
print('Running time: %s Seconds'%(end-start))