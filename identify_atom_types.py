from rdkit import Chem
import numpy as np

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


smiles = "CC"
atoms, connections = smiles_to_atoms_and_connections(smiles)
print("Atoms:", atoms)
print("Connections:", connections)



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
        #print(current_atom)
        #atom_types_0[i]=atoms[i]
        for j in range(len(atoms)):
            
            if connectivity[j][i]==1 and i!=j:
                
                atom_dict_1[i]+=atoms[j]
                atom_dict_2[i]+=atoms[j]+"("
                #print (atoms[j], "(")
                for k in range(len(atoms)):
                    if connectivity[k][j]==1 and j!=k:
                        atom_dict_2[i]+=atoms[k]
                        #print("-------",atoms[k])
                atom_dict_2[i]+=")"
                #print(")")
        atom_dict_1[i]+="]"
        atom_dict_2[i]+="]"    
        #print("--------------------end")
            
                        #tmp+=atoms[k]
                #tmp+=atoms[j]
        #tmp+=")"
    #print(atom_dict_0)
    #print(atom_dict_1)
    #print(atom_dict_2)
    return atom_dict_0, atom_dict_1, atom_dict_2

#test with depth 2
print(identify_from_connection_matrix_depth_2(atoms, connections))

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

