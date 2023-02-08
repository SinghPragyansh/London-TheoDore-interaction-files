import numpy as np
import re
with open("H2_RHF.out") as f:
    all_lines = f.readlines()
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
mo_list = []
for i in range(l):
    splitem = note_lines[i].strip().split()
    #print(splitem)
    temp = splitem[-1]
    mo_list.append(temp)
# Coverting mo coeff into proper complex numbered mo coeff
mo_coeff = []
for i in range(len(mo_list)):
    temp = mo_list[i]
    temp = re.sub(',', '+', temp)
    temp = temp[1:-1]
    temp += 'j'
    mo_coeff.append(temp)
    #print(mo_coeff)
# since the mo coefficients are repeating therefore let's half it
l = len(mo_coeff)
n = int(l/2)
mo_coeff = mo_coeff[0:n]
mo_array = np.zeros(nbas*nbas, dtype = complex)
for i in range(len(mo_coeff)):
    mo_array[i] = complex(mo_coeff[i])
    
mo_matrix = mo_array.reshape(nbas, nbas)
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
# there have been some extra empty line addition from the output file , therefore delte the last 
# two lines
orbital_lines = orbital_lines[:-2]
orb_ene = []
for orb_line in orbital_lines:
    orb_ene.append(orb_line.strip().split()[0])
occ_detail = []
for orb_line in orbital_lines:
    occ_detail.append(orb_line.strip().split()[2])
pat2 = re.compile(r"charge =")
with open("H2_RHF.inp") as f:
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
parsing_file = open('mo_info_extraction.molden', 'w')
parsing_file.write(f"[Molden Format]\nmade by MagChem_lab [0.0.0]\n")
parsing_file.write("[Atoms] (AU)\n")

for i in range(len(labels)):
    parsing_file.write(f"{labels[i]}   {i+1}   {charges[i]}    {coord_info[i]}\n")

parsing_file.write(f"[MO]\n")

for i in range(nbas):
    parsing_file.write(f" Sym= A\n Ene= {orb_ene[i]}\n") 
    parsing_file.write(f" Spin= Alpha\n Occp= {occ_detail[i]}\n")    
    for j in range(nbas):
        parsing_file.write(f"   {j+1}     {mo_matrix[j,i].real:g}+{mo_matrix[j,i].imag:g}j\n")

parsing_file.close()










