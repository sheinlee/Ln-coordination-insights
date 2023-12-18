from ccdc import io 
from ccdc.io import MoleculeReader
import numpy as np
import pandas as pd
from itertools import compress
from time import process_time
start=process_time()

def Count_Metal(molecule,element):
    countmetal=0
    for atom in molecule.atoms:
        if atom.atomic_symbol == element:
            countmetal+=1
    return countmetal


#count_metal==1 component/molecule
def Component_CN(component,element):
    component_cn=0
    for atom in component.atoms:
        if atom.atomic_symbol==element:
            component_cn=len(atom.neighbours)
            # print(component_cn)
            break
    return component_cn

#count_metal==1
def Average_first_shell_distance(molecule,element):
    for atom in molecule.atoms:
        if atom.atomic_symbol==element:
            fil = [bond.atoms.__contains__(atom) for bond in molecule.bonds]
            bonds_has_Ln = list(compress(molecule.bonds, fil))
            bond_length=0
            count_for_bond=0
            for bond in bonds_has_Ln:
                try:
                    bond.length
                except:
                    print("ignore H")
                else:
                    count_for_bond+=1
                    bond_length+=bond.length
            break
    if count_for_bond==0:
        return None
    else:
        return bond_length/count_for_bond
    


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

#count,refcode, component index, component CN, Average first shell distance.          metal number==1


writerA = pd.ExcelWriter('total_dataset.xlsx')

for ELEMENT in LIST_OF_ELEMENT:
    
    A = np.empty(shape=(7000,5),dtype=object)
    COUNT = 0
    

    filepath='refcode_'+ELEMENT+'.txt'
    csd_reader=io.EntryReader('csd')
    mol_reader = MoleculeReader(filepath, format='identifiers')
    for mol in mol_reader:
        COUNT += 1
        A[COUNT-1,0]=COUNT
        print(COUNT)
        molecule_name=mol.identifier
        # print(molecule_name)
        mol_name=csd_reader.molecule(molecule_name)
        A[COUNT-1,1]=mol_name.identifier
        
        #remove structures without 3D
        print(mol_name.is_3d)
        if mol_name.is_3d == False: 
            # with open(ELEMENT+"no3d.txt","a") as f:
            #     f.write(molecule_name)
            #     f.write("\n")
            continue

        #add_hydrogens()
        # mol_name.assign_bond_types()
        # mol_name.add_hydrogens()

        print(mol_name.identifier)
        print('ELEMENT',ELEMENT)

        
        A[COUNT-1,2]=Count_Metal(mol_name,ELEMENT)

        for atom in mol_name.atoms:
            if atom.atomic_symbol==ELEMENT:
                A[COUNT-1,3]=len(atom.neighbours)
                fil = [bond.atoms.__contains__(atom) for bond in mol_name.bonds]
                bonds_has_Ln = list(compress(mol_name.bonds, fil))
                bond_length=0
                count_for_bond=0
                for bond in bonds_has_Ln:
                    try:
                        bond.length
                    except:
                        print("ignore H")
                    else:
                        count_for_bond+=1
                        bond_length+=bond.length
                break
        if count_for_bond!=0:
            A[COUNT-1,4]=bond_length/count_for_bond
              


        # A[COUNT-1,4]=Component_CN(mol_name,ELEMENT)
        # A[COUNT-1,5]=Average_first_shell_distance(mol_name,ELEMENT)

        # if len(mol_name.components)==1:
        #     if Count_Metal(mol_name,ELEMENT)==1:
        #         A[COUNT-1,2]=mol_name.components[0].identifier
        #         A[COUNT-1,3]=Component_CN(mol_name.components[0],ELEMENT)
        #         A[COUNT-1,4]=Average_first_shell_distance(mol_name,ELEMENT)
        #         COUNT+=1
        #     else:
        #         continue
        # else:
        #     for com in mol_name.components:
        #         if Count_Metal(com,ELEMENT)==1:
        #             A[COUNT-1,2]=com.identifier
        #             A[COUNT-1,3]=Component_CN(com,ELEMENT)
        #             A[COUNT-1,4]=Average_first_shell_distance(com,ELEMENT)
        #             COUNT+=1
        #             # break
        #         else:
        #             continue

    mol_reader.close()


    print(ELEMENT)
    pd.DataFrame(A).to_excel(writerA,sheet_name=ELEMENT,float_format='%.5f')

writerA.save()

#记录程序运行的时间
end=process_time()
print('Running time: %s Seconds'%(end-start))