from ccdc import io 
import os
from ccdc.io import MoleculeReader
import numpy as np
import pandas as pd
from itertools import compress
from time import process_time

start=process_time()
LIST_OF_ELEMENT = ['La',
                   'Ce',
                   'Pr',
                   'Nd',
                   'Sm',
                   'Eu',
                   'Gd',
                   'Tb',
                   'Dy',
                   'Ho',
                   'Er',
                   'Tm',
                   'Yb',
                   'Lu'
                   ]

writerA = pd.ExcelWriter('element_distribution.xlsx')
A = np.zeros(shape=(16,41))
count_element=-1
for ElEMENT in LIST_OF_ELEMENT:
    count_element+=1
    filepath='corrected_asymmetric_mono_'+ElEMENT+'.txt'
    csd_reader=io.EntryReader('csd')
    mol_reader = io.MoleculeReader(filepath, format='identifiers')
    count_first_shell_atoms=0
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
    #org
    countb_o=0
    countc_o=0
    countn_o=0
    counts_o=0
    countp_o=0
    countsi_o=0
    #inorg
    counth_o=0
    countcl_o=0
    countbr_o=0
    counti_o=0
    countw_o=0
    countmo_o=0
    countpd_o=0
    countmn_o=0
    countv_o=0
    #O-C-R type of ligands 
    counto_is_amide=0
    counto_is_carboxylate_acid=0
    counto_is_ketone=0
    counto_is_ether=0
    counto_is_ester=0
    counto_is_alkoxide=0
    counto_is_hydroxyl=0

    counto_is_acid_and_water=0  
    counto_is_single_O=0
    counto_is_water=0
    counto_is_no3=0
    counto_is_so4=0
    for mol in mol_reader:
        molecule_name=mol.identifier
        mol_name=csd_reader.molecule(molecule_name)
        try:
            mol_name.assign_bond_types()
            mol_name.add_hydrogens()
        except RuntimeError as e:
            print(f"Error adding hydrogens to molecule {molecule_name}: {e}")
            continue
        find_element=False
        for atom in mol_name.atoms:
            if atom.atomic_symbol == ElEMENT:
                find_element=True

                if len(atom.neighbours)==0:
                    continue
                else:
                        
                    for atom1 in atom.neighbours:
                        count_first_shell_atoms+=1
                        if atom1.atomic_symbol == 'H': counth +=1
                        if atom1.atomic_symbol == 'C': countc +=1
                        if atom1.atomic_symbol == 'N': countn +=1  
                        if atom1.atomic_symbol == 'O': 
                            counto +=1
                            print('len(atom1.neighbours)',len(atom1.neighbours))
                            print(atom1.neighbours)
                            if len(atom1.neighbours)==1:
                                counto_is_single_O+=1
                            
                            if len(atom1.neighbours)>1:
                                #org
                                for atom2 in atom1.neighbours:
                                    if atom2.atomic_symbol=='B': 
                                        countb_o+=1 
                                        break
                                for atom2 in atom1.neighbours:
                                    if atom2.atomic_symbol=='C': 
                                        countc_o+=1
                                        break
                                for atom2 in atom1.neighbours:
                                    if atom2.atomic_symbol=='N': 
                                        countn_o+=1
                                        break
                                for atom2 in atom1.neighbours:
                                    if atom2.atomic_symbol=='S': 
                                        counts_o+=1
                                        break
                                for atom2 in atom1.neighbours:
                                    if atom2.atomic_symbol=='P': 
                                        countp_o+=1
                                        break
                                for atom2 in atom1.neighbours:
                                    if atom2.atomic_symbol=='Si': 
                                        countsi_o+=1
                                        break
                                    
                                #inorg
                                for atom2 in atom1.neighbours:
                                    if atom2.atomic_symbol=='H': 
                                        counth_o+=1
                                        break
                                for atom2 in atom1.neighbours:
                                    if atom2.atomic_symbol=='Cl': 
                                        countcl_o+=1
                                        break
                                for atom2 in atom1.neighbours:
                                    if atom2.atomic_symbol=='Br': 
                                        countbr_o+=1
                                        break
                                for atom2 in atom1.neighbours:
                                    if atom2.atomic_symbol=='I': 
                                        counti_o+=1
                                        break
                                for atom2 in atom1.neighbours:
                                    if atom2.atomic_symbol=='W': 
                                        countw_o+=1
                                        break
                                for atom2 in atom1.neighbours:
                                    if atom2.atomic_symbol=='Mo': 
                                        countmo_o+=1
                                        break
                                for atom2 in atom1.neighbours:
                                    if atom2.atomic_symbol=='Pd': 
                                        countpd_o+=1
                                        break
                                for atom2 in atom1.neighbours:
                                    if atom2.atomic_symbol=='Mn': 
                                        countmn_o+=1
                                        break
                                for atom2 in atom1.neighbours:
                                    if atom2.atomic_symbol=='V': 
                                        countv_o+=1
                                        break
                                #search h2o
                                counth_ard_o=0
                                for atom2 in atom1.neighbours:
                                    if atom2.atomic_symbol=='H':
                                        counth_ard_o+=1
                                if counth_ard_o==2:
                                    counto_is_water+=1
                                #search no3
                                counto_ard_n=0
                                for atom2 in atom1.neighbours:
                                    if atom2.atomic_symbol=='N':
                                        for atom3 in atom2.neighbours:
                                            if atom3.atomic_symbol=='O':
                                                
                                                counto_ard_n+=1
                                if counto_ard_n==3:
                                    counto_is_no3+=1
                                #search so4
                                counto_ard_n=0
                                for atom2 in atom1.neighbours:
                                    if atom2.atomic_symbol=='S':
                                        for atom3 in atom2.neighbours:
                                            if atom3.atomic_symbol=='O':
                                                
                                                counto_ard_n+=1
                                if counto_ard_n==4:
                                    counto_is_so4+=1
                                
                            
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
                                        countn_ard_c=0
                                        for atom3 in atom2.neighbours:
                                            if atom3.atomic_symbol=='N':
                                                if len(atom3.neighbours)==3:
                                                    countn_ard_c+=1
                                        if countn_ard_c>0:
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
                                                        if atom4.atomic_symbol=='H' or atom4.is_metal==True or len(atom3.neighbours)==1:
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
                            #search alkoxide
                            for atom2 in atom1.neighbours:

                                if atom2.atomic_symbol=='C':

                                    for bond in mol_name.bonds:
                                        if (bond.atoms[0]==atom1 and bond.atoms[1]==atom2) or (bond.atoms[0]==atom2 and bond.atoms[1]==atom1):
                                            print(bond.bond_type)
                                            if bond.bond_type==1:
                                                counto_aroundc=0
                                                flago_is_ionic=False
                                                for atom3 in atom2.neighbours:
                                                    if atom3.atomic_symbol=='O':
                                                        counto_aroundc+=1
                                                        
                                                        
                                                        
                                                counthormetal_aroundo=0
                                                counth_aroundo=0
                                                for atom4 in atom1.neighbours:
                                                    if atom4.atomic_symbol=='H'or atom4.is_metal==True:
                                                        counthormetal_aroundo+=1
                                                        counth_aroundo+=1
                                                if len(atom1.neighbours)-counthormetal_aroundo==1:
                                                    flago_is_ionic=True
                                                    print('counthormetal_aroundo',counthormetal_aroundo)
                                                if counto_aroundc==1 and flago_is_ionic:
                                                    counto_is_alkoxide+=1
                                                # if counto_aroundc==1 and flago_is_ionic==False and counth_aroundo==1:
                                                #     counto_is_hydroxyl+=1  


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
                        if atom1.atomic_symbol == 'S': counts +=1
                        if atom1.atomic_symbol == 'P': countp +=1
                        if atom1.atomic_symbol == 'F': countf +=1
                        if atom1.atomic_symbol == 'Cl': countcl +=1
                        if atom1.atomic_symbol == 'Br': countbr +=1
                        if atom1.atomic_symbol == 'I': counti +=1
            if find_element: break
    A[count_element,2]=count_first_shell_atoms                    
    A[count_element,3]=counth
    A[count_element,4]=countc
    A[count_element,5]=countn
    A[count_element,6]=counto
    A[count_element,7]=counts
    A[count_element,8]=countp
    A[count_element,9]=countf
    A[count_element,10]=countcl
    A[count_element,11]=countbr
    A[count_element,12]=counti

    A[count_element,13]=countb_o
    A[count_element,14]=countc_o
    A[count_element,15]=countn_o
    A[count_element,16]=counts_o
    A[count_element,17]=countp_o
    A[count_element,18]=countsi_o
    
    A[count_element,19]=counth_o
    A[count_element,20]=countcl_o
    A[count_element,21]=countbr_o
    A[count_element,22]=counti_o
    A[count_element,23]=countw_o
    A[count_element,24]=countmo_o
    A[count_element,25]=countpd_o
    A[count_element,26]=countmn_o
    A[count_element,27]=countv_o


    A[count_element,28]=counto_is_amide
    A[count_element,29]=counto_is_carboxylate_acid
    A[count_element,30]=counto_is_ketone
    A[count_element,31]=counto_is_ester
    A[count_element,32]=counto_is_ether
    A[count_element,33]=counto_is_alkoxide
    A[count_element,34]=counto_is_hydroxyl
    
    A[count_element,35]=counto_is_single_O
    A[count_element,36]=counto_is_water
    A[count_element,37]=counto_is_no3
    A[count_element,38]=counto_is_so4
    A[count_element,39]=counto_is_acid_and_water
    print(A)
    print(A[count_element,:])
    print(count_element)
    mol_reader.close()
pd.DataFrame(A).to_excel(writerA,sheet_name='element_distribution',float_format='%.5f')
writerA.save()
#record process time
end=process_time()
print('Running time: %s Seconds'%(end-start))