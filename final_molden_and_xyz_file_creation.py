import numpy as np
import re

file = open("HOF_RHF.out", "r")
all_lines = file.readlines()

pat1 = re.compile(r"MO coeff matrix elements:..")
pat = re.compile(r"...The SCF basis contains...")
basis_line = []
for line in all_lines:
    for match in re.finditer(pat, line):
        basis_line.append(line)

nbas = int(basis_line[0].split()[5])
note_lines = []
for line in all_lines:
    for match in re.finditer(pat1, line):
        note_lines.append(line)

l = len(note_lines)
mo_coeff = []
for i in range(l):
    splitem = note_lines[i].strip().split()
    #print(splitem)
    temp = splitem[-1]
    mo_coeff.append(temp)
    
# since the mo coefficients are repeating therefore let's half it
l = len(mo_coeff)
n = int(l/2)
mo_coeff = mo_coeff[0:n]

            ### Converting mo_coeff list into the standard mo_ao matrix
#    The columns of this matrix would give us the required molecular orbitals
#    The design of matrix is such that a particular i,j index element would repesent, the contribution of jth ao in ith mo.
#    mo1 = linear combination of all aos
#    for 6 primitive basis, we will have six aos and therefore 6 mos
#    every column of the matrix would represent a particular mo.

def comp_coeff_parse(line):
    #c = line.strip().split()[-1]
    # let's find the comma inside the string to distinguish the real and imaginary part of the coeff
    s = line.rfind(',')
    # the real part will be from 2nd char to just before comma and imag part starts just after comma
    # and ends on the last second character of the string
    real_c = float(line[1:s])
    im_c = float(line[s+1:-1])
    comp_c = complex(real_c, im_c)
    
    return comp_c

mo_array = np.zeros(nbas*nbas, dtype = complex)
for i in range(len(mo_coeff)):
    mo_array[i] = comp_coeff_parse(mo_coeff[i])

mo_matrix = mo_array.reshape(nbas, nbas)

###  Parsing the molecular orbital energies and occupancy

pat_mo = re.compile(r"Orbital analysis...")
pat_end = re.compile(r"Properties:")

orbital_lines = []
for idx, line in enumerate(all_lines):
    if re.match(pat_mo, line.strip()):
        break
idx = idx+3
line = all_lines[idx]
while not(re.match(pat_end, line.strip())):
    orbital_lines.append(line.strip())
    idx += 1
    line = all_lines[idx]

# there have been some extra empty line addition from the output file , therefore delete the last 
# two lines
orbital_lines = orbital_lines[:-2]

orb_ene = []
occ_detail = []
for orb_line in orbital_lines:
    orb_ene.append(orb_line.strip().split()[0])
    occ_detail.append(orb_line.strip().split()[2])

#### So far  we have the mo_ao_matrix, orb_ene_list, occ_detail_list
#### Let's extract the charge info and co-ordinate infos

pat2 = re.compile(r"charge =")
with open("HOF_RHF.inp") as f:
    inp_lines = f.readlines()
charge_info = []
for line in inp_lines:
    for match in re.finditer(pat2, line):
        charge_info.append(line)

# Collecting charge information
molecular_charge = charge_info[0].strip()
atomic_charge = charge_info[1:]
charges = []
for line in atomic_charge:
    charges.append(line.strip().split()[2])

pat_st = re.compile(r"# geometry")
pat_en = re.compile(r"}")

for idx, line in enumerate(inp_lines):
    if re.match(pat_st, line.strip()):
        break       
        
# since coordinate info starts after two lines after of above mentioned starting pattern
idx = idx +2
line = inp_lines[idx]
atom_info = []
while not(re.match(pat_en, line.strip())):
    atom_info.append(line.strip())
    # we need alternate lines as in between coordinate info we have charge info as well
    # so we keep increment to "+2"
    idx = idx + 2
    line = inp_lines[idx]

labels = [] 
coord_info = []
for atom in atom_info:
    labels.append(atom.strip().split()[0])    
    coord_info.append(atom.strip()[1:])

#### Now we also have the labels_list, coordinate_list, charge_list separately

'''parsing_file = open('mo_info_extraction.molden', 'w')
parsing_file.write(f"[Molden Format]\nmade by MagChem_lab [0.0.0]\n")
parsing_file.write("[Atoms] (AU)\n")

for i in range(len(labels)):
    parsing_file.write(f"{labels[i]}   {i+1}   {charges[i]}    {coord_info[i]}\n")

parsing_file.write(f"[MO]\n")

for i in range(len(orb_ene)):
    parsing_file.write(f" Sym= A\n Ene= {orb_ene[i]}\n") 
    parsing_file.write(f" Spin= Alpha\n Occp= {occ_detail[i]}\n")    
    for j in range(nbas):
        parsing_file.write(f"   {j+1}     {mo_matrix[j,i].real:g}+{mo_matrix[j,i].imag:g}j\n")

parsing_file.close()'''

# writing the file the geom.xyz file
geom_xyz = open('geom.xyz', 'w')
geom_xyz.write(f"8\n\n")

for i in range(len(coord_info)):
    geom_xyz.write(f"{labels[i]}          {coord_info[i]}\n")
    
geom_xyz.close()



    ### Let's parse the basis info according to the format required

patst = re.compile(r"Basis file..")
patend = re.compile(r"User input successfully read..")
pat1 = re.compile("\$\s+..TYPE FUNCTIONS")   
# \$ : to read the '$',
# \s+ : to detect the just after space
# .. : to match presence of any word or letter
pat2 = re.compile(r"    charge =")
pat = re.compile(r"a \d")

# bulding the keys that can be used later for dictionary
atomic_numbers = []
shell_types = []
for line in all_lines:
    for match in re.finditer(pat1, line.strip()):
        shell_types.append(line.strip().split()[1][0])
        
    for match in re.finditer(pat2, line):
        number = line.strip().split()[-1]
        atomic_numbers.append(number)

def sent_prim_rows(prim_rows):
    n_columns = len(prim_rows[0].strip().split())
    n_rows = len(prim_rows)
    mat = np.zeros((n_rows, n_columns))
    for i in range(n_rows):
        for j in range(n_columns):
            mat[i,j] = prim_rows[i].strip().split()[j]
            
    return mat

# this line belongs to the S shell of first atom
def cal_inner_dict(idx, line, patst, patend, pat1):
    
    last_type = line.strip().split()[1]

    fun_dict = {}
    prim_rows = []
    while not re.match(patst, line.strip()):
        if re.match(pat1, line.strip()):
            type_ = line.strip().split()[1]
            idx +=1
            line = all_lines[idx]

            if prim_rows != []:
                mat = sent_prim_rows(prim_rows)
                fun_dict[f"{last_type[0]}"] = mat
                prim_rows = []
            last_type = type_
        
        elif re.match(patend, line.strip()):
            mat = sent_prim_rows(prim_rows)
            fun_dict[f"{last_type[0]}"] = mat
            prim_rows = []
            return fun_dict
            break
            
        else:
            #print(line.strip())
            prim_rows.append(line)

            idx+=1
            line = all_lines[idx]

    mat = sent_prim_rows(prim_rows)
    fun_dict[f"{last_type[0]}"] = mat
    
    return fun_dict

    # parsing for 6-31G
    # 6 exponents present in 1st column
    # and 6 contraction coeffs present in 2nd column

   # If we can read the "a .." line representing atomic number, then we can define a function which will send the this matching line into function, which can parse the shell info like exponent and coeffs etc. In the above code, we have performed parsing for "O" atom.

# Constructing the outer dictionary according to the atomic number
def outer_dict_atomic_num(pat,patst, patend, pat1, all_lines):
    # finding the atomic number line
    atomic_line_idx = []
    for idx, line in enumerate(all_lines):
        if re.match(pat, line.strip()):
            atomic_line_idx.append(idx)

    atomic_basisfun_nested_dict = {}
    for i in range(len(atomic_line_idx)):
        # we need to send the index and line of first basis function of atoms
        idx = atomic_line_idx[i] + 2
        line = all_lines[idx]

        atomic_number = atomic_numbers[i]

        atomic_basisfun_nested_dict[f"{atomic_number}"] = cal_inner_dict(idx, line, patst, patend, pat1)
        #print(cal_basisfun_dict(idx, line))
        
    return atomic_basisfun_nested_dict

def shell_gaussian_splitting(fun_mat,key):
    l = fun_mat.shape[1]
    coeffs = {}
    expos = {}
    t = 0
    for i in range(1,l):
        array = fun_mat[:,i]
        array = array[array != 0]
        coeffs[f"{key}{i}_coeff"] = array
        
        z = len(array)
        expos[f"{key}{i}_expo"] = fun_mat[t:t+z,0]
        
        t = t+z
        
    return coeffs, expos

def checkKey(dic, key):
    if key in dic.keys():
        return dic[key]
    else:
        return np.zeros(1)

def write_coeffs_expos(fun_mat, coeffs, expos, key, atom_ind, file):
    # l = number of columns in the shell,
    # where first column is expo and others are for number of slater functions
    # eg. if there are 3 types of contractions in "S" shell of "O" atom
    # 6 for core 1s and 3,1 for 2S valence shell , then l = 4 here
    
    l = fun_mat.shape[1]
    
    lowercase_key = key.lower()
    
    file.write(f"   {atom_ind} 0\n")
    for i in range(1,l):
        #print(f"{key}{i} coefficients:\n",coeffs[f"{key}{i}_coeff"])
        #print(f"{key}{i} expoenents:\n",expos[f"{key}{i}_expo"])
        
        coeff = coeffs[f"{key}{i}_coeff"]
        expo = expos[f"{key}{i}_expo"]
        u = len(coeff)
        
        file.write(f" {lowercase_key}    {u}  1.00\n")
        
        for j in range(u):
            file.write(f"         {expo[j]}             {coeff[j]}\n")



keyset1 = atomic_numbers
keyset2 = list(set(shell_types))
if keyset2[0] != "S":
    keyset2 = keyset2[::-1]
outer_dict = outer_dict_atomic_num(pat,patst, patend, pat1, all_lines)
#print(outer_dict)

'''file = open("basis_formatting.txt", "w")

for i in range(len(keyset1)):
    for j in range(len(keyset2)):
        # outer dict contains matrices according to atomic number ket
        # by inserting one key, we are narrowing the range to one atomic number and all shells
        outer_dict_keyi = outer_dict[keyset1[i]]
     
        # one inner dict will contain one shell and all the gaussian functions in that shell
        # but we must also check, whether some shells are not present in the some atom
        
        inner_dict_keyj = checkKey(outer_dict_keyi, keyset2[j])
        #print(inner_dict_keyj)
        # since we have zoomed in to one shell, lets split according to diff contractions
        if inner_dict_keyj.any() != 0:
            coeffs, expos = shell_gaussian_splitting(inner_dict_keyj,keyset2[j])
            write_coeffs_expos(inner_dict_keyj, coeffs, expos, keyset2[j], i+1, file)
            
file.close()'''

parsing_file = open('geom.molden', 'w')
parsing_file.write("[Molden Format]\nmade by MagChem_lab [0.0.0]\n")
parsing_file.write("[Atoms] (AU)\n")

for i in range(len(labels)):
    parsing_file.write(f"{labels[i]}   {i+1}   {charges[i]}    {coord_info[i]}\n")

parsing_file.write("[GTO]\n")

for i in range(len(keyset1)):
    for j in range(len(keyset2)):
        # outer dict contains matrices according to atomic number ket
        # by inserting one key, we are narrowing the range to one atomic number and all shells
        outer_dict_keyi = outer_dict[keyset1[i]]
     
        # one inner dict will contain one shell and all the gaussian functions in that shell
        # but we must also check, whether some shells are not present in the some atom
        
        inner_dict_keyj = checkKey(outer_dict_keyi, keyset2[j])
        #print(inner_dict_keyj)
        # since we have zoomed in to one shell, lets split according to diff contractions
        if inner_dict_keyj.any() != 0:
            coeffs, expos = shell_gaussian_splitting(inner_dict_keyj,keyset2[j])
            write_coeffs_expos(inner_dict_keyj, coeffs, expos, keyset2[j], i+1, parsing_file)
    

parsing_file.write("\n[MO]\n")

for i in range(len(orb_ene)):
    parsing_file.write(f" Sym= A\n Ene= {orb_ene[i]}\n") 
    parsing_file.write(f" Spin= Alpha\n Occp= {float(occ_detail[i])}\n")    
    for j in range(nbas):
        parsing_file.write(f"   {j+1}     {mo_matrix[j,i].real:g}+{mo_matrix[j,i].imag:g}j\n")

parsing_file.close()




