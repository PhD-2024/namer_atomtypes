from rdkit import Chem
import numpy as np
import re

'''
This script will be split into a library file and main later
The functions will therefore be moved. If we want to later reuse stuff, we can also make a class based on struct of atoms and connections, which will get functions as methods
'''
def smiles_to_atoms_and_connections(smiles):
    '''
    takes a smile and returns a list of atoms and the connection matrix
    '''
    mol = Chem.MolFromSmiles(smiles)
    mol = Chem.AddHs(mol)
    atoms = [atom.GetSymbol() for atom in mol.GetAtoms()]
    natoms = len(atoms)
    
    connections = [[0] * natoms for _ in range(natoms)]
    for bond in mol.GetBonds():
        i = bond.GetBeginAtomIdx()
        j = bond.GetEndAtomIdx()
        connections[i][j] = 1
        connections[j][i] = 1
    
    return atoms, connections




def identify_from_connection_matrix_depth_2(atoms, connectivity):
    '''creates a dictionary (each) up to depth 2
    These are indexed by the atom indices as given in the atom list, 
    This by itself is insufficient, as they are not yet ordered'''
    atom_dict_0={}
    atom_dict_1={}
    atom_dict_2={}
    for i in range(len(atoms)):
        atom_dict_1[i]=""
        atom_dict_2[i]=""

    #connectivity=np.linalg.matrix_power(connectivity, 2)
    for i in range(len(atoms)):
        current_atom=atoms[i]
        #tmp=""
        #tmp+=current_atom+ " ("
        zero_string=current_atom
        atom_dict_0[i]=zero_string
        atom_dict_1[i]+=atoms[i]+"["
        atom_dict_2[i]+=atoms[i]+"["

        for j in range(len(atoms)):
            
            if connectivity[j][i]==1 and i!=j:
                
                atom_dict_1[i]+=atoms[j]
                atom_dict_2[i]+=atoms[j]+"("
                
                for k in range(len(atoms)):
                    if connectivity[k][j]==1 and j!=k:
                        atom_dict_2[i]+=atoms[k]
                        
                atom_dict_2[i]+=")"
                
        atom_dict_1[i]+="]"
        atom_dict_2[i]+="]"    

    return atom_dict_0, atom_dict_1, atom_dict_2

def sort_within_bracket(string_name):
    '''This function will convert the simple case of the alphabetical sort
    The ide is: string to list, then sort list and return as string'''
    unsorted=list(string_name)
    #check that the string is neither empty and does not contain chars, that are not letters
    if len(unsorted)==0:
        raise ValueError("String is empty")
    if not all(char.isalpha() for char in unsorted):
        print(unsorted)
        raise ValueError("String contains non-alphabetical characters")

    return "".join(sorted(unsorted))




def sort_by_name_and_length(collected):
    '''This is used as a self-defined sort key, to make the depth 2 unambigous.
    It sorts the collected list by the atom name and the length of the inner string'''
    return sorted(collected, key=lambda x: (x[0], x[1]))

def find_parentheses_substrings(string):
    '''finds all substrings in parentheses'''
    pattern = re.compile(r'\(.*?\)')
    return pattern.findall(string)

def sort_in_brackets(string_name):
    '''Sorts the string into its bracket contributions
    case 1 is dictionary after one iteration, case 2 
    is dictionary after two iterations'''
    #case 1
    if "[" in string_name: 
        if "(" not in string_name and ")" not in string_name:
            if "]"  not in string_name:
                raise ValueError("Brackets are not balanced")
            else:
                if string_name.count("[")!=string_name.count("]"):
                    raise ValueError("Brackets are not balanced")
                else:
                    print("assuming interation 1")
                    tmp=string_name.replace("]","").split(sep="[")[1:] 
                    return sort_within_bracket(tmp)
        #case2        
        elif "(" in string_name:
            if ")" not in string_name:
                raise ValueError("Brackets are not balanced")
            else:
                if string_name.count("(")!=string_name.count(")"):
                    raise ValueError("Brackets are not balanced")
                else:
                #first sort the inner brackets "( )"
                # then they should be equal and we can sort the outer brackets "[ ]"
                    #find all substrings in parentheses and sort them individually
                    all_substrings=find_parentheses_substrings(string_name)


                    for i in range(len(all_substrings)):
                        original_substring=all_substrings[i]
                        #print("to sort", original_substring)
                        sorted_substring="("+sort_within_bracket(all_substrings[i].replace("(","").replace(")",""))+")"
                        string_name=string_name.replace(original_substring,sorted_substring)
                        
                    #now sort outer brackets
                    #outer bracket will be sorted alphabetically and if the atoms are the same, the length of the inner string will be compared
                    
                    B,working_on=string_name.replace("]","").split(sep="[")
                    
                    all_substrings=working_on.split(sep=")")[:-1] # last is empty
                    
                    collected=[]
                    for item in all_substrings:
                        #print("sub",item)
                        atom_char, inner=item.split(sep="(")
                        length=len(inner)
                        collected.append((atom_char,length, inner))
                    sorted=sort_by_name_and_length(collected)
                    stringified=[a[0]+"("+a[2]+")" for a in sorted]
                    #print(sorted)
                    #print(stringified)

                    modified_outer="".join(str(x) for x in stringified)
                                           
                    final_string=B+"["+modified_outer+"]"
                    return final_string

def read_smiles_from_file(file:str="smiles.txt", header:bool=True):
    '''reads smiles from file, returns list of smiles'''
    all_smiles=[]
    with open(file, "r") as f:
        if header:
            f.readline()
        for line in f:
            if not line.startswith(";"):
                #if it is not an empty line
                if line.strip()!="":
                    all_smiles.append(line.strip())
        return all_smiles

    
def main():
    
    all_smiles=read_smiles_from_file("smiles.txt")
    all_types=[]

    #test with depth 2
    for smiles in all_smiles:
        print(smiles)

        atoms, connections = smiles_to_atoms_and_connections(smiles)
        a0,a1,a2=identify_from_connection_matrix_depth_2(atoms, connections)
        #das hier sind die sortierten atomtypen für den graph
        a2_0_sorted=   [ sort_in_brackets(a2[i]) for i in range(len(a2))]
        
        #print(a2_0_sorted)
        #all types sind NICHT sortiert
        all_types.extend(a2_0_sorted)
    print(set(all_types))
    print("Number of different atomtypes identified:", len(set(all_types)))
    print("from a total of", len(all_smiles), "molecules")  
    return set(all_types)




if __name__ == "__main__":
        main()