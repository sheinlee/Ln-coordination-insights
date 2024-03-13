from ccdc import io 
from ccdc.io import MoleculeReader, CrystalReader,EntryReader
import numpy as np
import pandas as pd
from itertools import compress
import re
from time import process_time
start=process_time()

def Count_Metal(molecule,element):
    countmetal=0
    for atom in molecule.atoms:
        if atom.atomic_symbol == element:
            countmetal+=1
    return countmetal

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

writerA = pd.ExcelWriter('monon_molecule.xlsx')
writerB = pd.ExcelWriter('mono_crystal.xlsx')
csd_reader=io.EntryReader('csd')
for ELEMENT in LIST_OF_ELEMENT:
    
    A = np.empty(shape=(12000,5),dtype=object)
    B = np.empty(shape=(12000,5),dtype=object)
    element_pattern = re.compile(rf'{ELEMENT}(\d+)')
    COUNT = 0
    counta=0
    countb=0
    counta_is_mono=False
    countb_is_mono=False
    filepath='refcode_'+ELEMENT+'.txt'
    # corrected_asymmetric_mono_filename = f'corrected_asymmetric_mono_{ELEMENT}.txt'
    
    mol_reader = MoleculeReader(filepath, format='identifiers')
    for mol in mol_reader:
        COUNT += 1
        A[COUNT-1,0]=COUNT
        print(COUNT)
        molecule_name=mol.identifier
        mol_name=csd_reader.molecule(molecule_name)
        A[COUNT-1,1]=mol_name.identifier
        
        #remove structures without 3D
        print(mol_name.is_3d)
        if mol_name.is_3d == False: 
            continue
        #add_hydrogens()
        # mol_name.assign_bond_types()
        # mol_name.add_hydrogens()

        if len(mol_name.components)==1:
            if Count_Metal(mol_name,ELEMENT)==1:
                counta_is_mono=True
            else:
                continue
        else:
            for com in mol_name.components:
                if Count_Metal(com,ELEMENT)==1:
                    counta_is_mono=True
                else:
                    continue
        if counta_is_mono: 
            counta+=1
            A[COUNT-1,2]=counta
            # with open('corrected_mono_'+ELEMENT+'.txt','a') as f:
            #     f.write(molecule_name + '\n')
    mol_reader.close()
    pd.DataFrame(A).to_excel(writerA,sheet_name=ELEMENT,float_format='%.5f')
    COUNT = 0

    crystal_reader = CrystalReader('CSD')
    cry_reader = CrystalReader(filepath, format='identifiers')
    entry_reader = EntryReader('CSD')
    # en_name = entry_reader.entry(name)
    for cry in cry_reader:
        COUNT += 1
        B[COUNT-1,0]=COUNT
        print(COUNT)
        crystal_name=cry.identifier
        cry_name=cry_reader.crystal(crystal_name)
        B[COUNT-1,1]=cry_name.identifier
        entry_name=csd_reader.entry(crystal_name)
        #remove structures without 3D
        print(entry_name.has_3d_structure)
        if entry_name.has_3d_structure == False: 

            continue
        #add_hydrogens()
        # mol_name.assign_bond_types()
        # mol_name.add_hydrogens()
        asy_cry_name=cry_name.asymmetric_unit_molecule
        if len(asy_cry_name.components)==1:
            if Count_Metal(asy_cry_name,ELEMENT)==1:
                countb_is_mono=True
            else:
                continue
        else:
            for com1 in asy_cry_name.components:
                if Count_Metal(com1,ELEMENT)==1:
                    countb_is_mono=True
                else:
                    continue
        if countb_is_mono: 
            formula = entry_name.formula
            # 使用正则表达式搜索当前元素及其数量
            match = element_pattern.search(formula)
            if match:
                element_count = int(match.group(1))  # 将匹配到的数字转换为整数
                if element_count == 1:
                    countb+=1
                    B[COUNT-1,2]=countb
                    with open('corrected_asymmetric_mono_'+ELEMENT+'.txt','a') as f:
                        f.write(cry_name.identifier + '\n')
                
    cry_reader.close()
    pd.DataFrame(B).to_excel(writerB,sheet_name=ELEMENT,float_format='%.5f')

writerA.save()
writerB.save()
print('counta',counta)
print('countb',countb)
#记录程序运行的时间
end=process_time()
print('Running time: %s Seconds'%(end-start))