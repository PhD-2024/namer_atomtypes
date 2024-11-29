from rdkit import Chem
import numpy as np
import re

def smiles_to_atoms_and_connections(smiles):
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
    #creates a dictionary with the atom types up to depth 2
    #the dataframe is indexed by the atoms, the value is given to define the connectivity strings
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
    '''string to list, then sort list and return as string'''
    unsorted=list(string_name)
    #check that the string is neither empty and does not contain chars, that are not letters
    if len(unsorted)==0:
        raise ValueError("String is empty")
    if not all(char.isalpha() for char in unsorted):
        print(unsorted)
        raise ValueError("String contains non-alphabetical characters")

    return "".join(sorted(unsorted))




def sort_by_name_and_length(collected):
    '''sorts the collected list by the atom name and the length of the inner string'''
    return sorted(collected, key=lambda x: (x[0], x[1]))

def find_parentheses_substrings(string):
    pattern = re.compile(r'\(.*?\)')
    return pattern.findall(string)

def sort_in_brackets(string_name):
    '''splits string into its bracket contributions
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


#test with depth 2
smiles = "CCO"
atoms, connections = smiles_to_atoms_and_connections(smiles)
print("Atoms:", atoms)
print("Connections:", connections)



a0,a1,a2=identify_from_connection_matrix_depth_2(atoms, connections)
a2_0_sorted=   [ sort_in_brackets(a2[i]) for i in range(len(a2))]
print(a2_0_sorted)
exit()
def identify_from_connection_matrix_depth_1(atoms, connectivity):
    atom_types_0={}
    for i in range(len(atoms)):
        current_atom=atoms[i]
        tmp=""
        tmp+=current_atom+ " ("
        atom_types_0[i]=atoms[i]
        for j in range(len(atoms)):
            if connectivity[j][i]==1 and i!=j:
                tmp+=atoms[j]
        tmp+=")"
        atom_types_0[i]=tmp
    return atom_types_0    


def get_connections_via_n(connections,n):
    '''returns the connectivity matrix along n steps (n>=1) - connection to itself is written as one
    this trick works because A**n,  paths of length n between two nodes'''
    for i in range(len(connections)):
        for j in range(len(connections)):
            if j==i:
                connections[i][j]=1
    
    M= np.linalg.matrix_power(connections, n)
    #for i in range(len(M)):
    #    for j in range(len(M)):
    #        if M[i][j]>0:
    #            M[i][j]=1
    return M


print(get_connections_via_n(connections,2
                            ))


def identify_from_connection_matrix_up_to_depth_N(atoms, connectivity, N):
    '''returns a dictionary with the atom types up to depth N, without previously defined connections'''
    types_order={}
    assert isinstance(N, int) and N >= 0, "N must be a positive integer"
    for i in range(N):
        types_order[i]={}

    #manually zeroth_order_types    
    for i in range(len(atoms)):
        types_order[0][i]=atoms[i]
    if N==0:
        return types_order
    #connections_first order already defined
    types_order[1]=identify_from_connection_matrix_depth_1(atoms, connectivity)
    if N==1:
        return types_order
    

    
    #stattdessen lieber connections des vorherigen atoms suchen

    #connections second order and higher
    #for i in range(2,N+1):
    #    connections=get_connections_via_n(connectivity,i)
    #    #connections_one_less=get_connections_via_n(connectivity,i-1)#to subtract so we do not have stuff multiple times
    #    #connections=connections-connections_one_less
    #    for j in range(len(atoms)):
    #        print("Connections for atom",atoms[j],"at depth",i,":")
    #        for k in range(len(atoms)):
    #            if connections[j][k]==1:
    #                print()
    #                print(atoms[k],end="")

    return types_order

print(identify_from_connection_matrix_up_to_depth_N(atoms, connections, 2))

