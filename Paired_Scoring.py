#!/usr/bin/env python
# coding: utf-8

# In[1]:


import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
import math
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdPartialCharges
from collections import defaultdict
from itertools import combinations
from itertools import product


# In[2]:


# empty list to read list from a file
smiles = []

# open file and read the content in a list
with open(r'smiles.txt', 'r') as fp:
    for line in fp:
        # remove linebreak from a current name
        # linebreak is the last character of each line
        x = line[:-1]

        # add current item to the list
        smiles.append(x)


# In[3]:





# In[4]:


#Function to calculate Gasteiger charges of each atom in each molecules
#output is a dictionary that pairs all atoms with their charges
def calculate_gasteiger_charges(smiles):
    molecule = Chem.MolFromSmiles(smiles)
    rdPartialCharges.ComputeGasteigerCharges(molecule)
    gasteiger_charges = defaultdict(list)
    for atom in molecule.GetAtoms():
        atom_symbol = atom.GetSymbol()
        atom_charge = atom.GetProp("_GasteigerCharge")
        gasteiger_charges[atom_symbol].append(float(atom_charge))
    return gasteiger_charges


# In[5]:


#Makes a list of all charges in all molecules
all_charges = defaultdict(list)
for smile in smiles:
    charges = calculate_gasteiger_charges(smile)
    for atom, charge in charges.items():
        all_charges[atom].extend(charge)


# In[6]:


#Removes all inf and NA
atom_dict = {}
for atom, charges in all_charges.items():
    charges = [x for x in charges if not math.isnan(x)]
    charges = [x for x in charges if not math.isinf(x)]
    atom_dict[atom] = charges


# In[7]:


elements = atom_dict.keys()
atom_pairs = list(product(elements, repeat=2))
charge_diffs = {}


# In[8]:


#itterates through all the dictionaries and calcualtes the differences
#in the respective atomic groups
for atom1, atom2 in atom_pairs:
    if atom1 in all_charges and atom2 in all_charges:
        combined_charges = all_charges[atom1] + all_charges[atom2]
        charge_diffs[(atom1, atom2)] = [abs(j-i) for i, j in combinations(combined_charges, 2)]


# In[9]:


pair_dif = {}
for pair, diffs in charge_diffs.items():
    pair_dif[pair] = diffs


# In[10]:


#Creates the scoring dictionary of dictionaries
score_pair = {}
for pair in pair_dif.keys():
    gast = pair_dif[pair]
    total_count = len(gast)
    score = {}
    
    for i in range(0,30):
        test = round(i*10**-1, 1)
        count = len([x for x in gast if x >= test])
        prob = count/total_count
        score[test] = np.log2(prob)
        
    score_pair[pair] = score
    


# In[12]:


#Function that changes inf to the largest penalty
def replace_neg_inf_with_max(dictionary):
    most_neg_value = min(value for inner_dict in dictionary.values() for value in inner_dict.values() if value != float('-inf'))
    for inner_dict in dictionary.values():
        for key, value in inner_dict.items():
            if value == float('-inf'):
                inner_dict[key] = most_neg_value

    return dictionary


# In[13]:


result = replace_neg_inf_with_max(score_pair)


# In[20]:


with open('score_pair_doubles.pkl', 'wb') as f:
    pickle.dump(result, f)


# In[43]:


#Refined dictionary for validation
refined_dict = {('C','C'): result['C','C'],('C','O'): result['C','O'],('O','O'): result['O','O']}


# In[40]:


with open('score_pair_doubles_refined.pkl', 'wb') as f:
    pickle.dump(refined_dict, f)

