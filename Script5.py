#!/usr/bin/env python
# coding: utf-8

# In[1]:


from Bio.KEGG import REST
import re
import pubchempy as pcp
from indigo import *


# In[2]:


#function that retreives the name of the compounds given a compound id
def get_compound_name(cpd_id):
    entry = REST.kegg_get(cpd_id).read()  # Fetch the compound entry from KEGG
    for line in entry.rstrip().split("\n"):  # Parse the entry line by line
        if line.startswith("NAME"):  # If the line contains the compound name
            name = line.split(" ", 1)[1]  # Extract the name
            return name.strip().rstrip(';')  # Remove leading/trailing whitespace and trailing semicolon


# In[3]:


#Retreives all the canonical SMILES in a given pathway id, and returns them as list
def pathway(kegg_map):
    result = REST.kegg_link('compound',kegg_map).read()
    
    lines = result.strip().split('\n')  # strip() removes leading/trailing whitespace (including newlines)
    # For each line, split by '\t' and get the second item
    cpd_codes = [line.split('\t')[1] for line in lines]
    
    names = []
    for codes in cpd_codes:
        names.append(get_compound_name(codes))
        
    smiles = []

    # Fetch the SMILES strings
    for compound in names:
        results = pcp.get_compounds(compound, 'name')
        if results:
            smiles.append(results[0].canonical_smiles)
            
    return(smiles)


# In[4]:


#pathway for citrate cycle
pathway('M00009')


# In[5]:


#SMILES list with co-enzymes removed
smiles = ('C(C(=O)[O-])C(CC(=O)[O-])(C(=O)[O-])O',
 'C(C(C(C(=O)O)O)C(=O)O)C(=O)O',
 'C(CC(=O)O)C(=O)C(=O)O',
 'C(CC(=O)[O-])C(=O)[O-]',
 'C(=CC(=O)[O-])C(=O)[O-]',
 'C(C(C(=O)O)O)C(=O)O','C(C(=O)C(=O)O)C(=O)O')


# In[6]:


#Re-canonization using indigo
indigo = Indigo()
for smile in smiles:
    mol1 = indigo.loadMolecule(smile)
    mol1.aromatize()
    print(mol1.canonicalSmiles())


# In[7]:


#Further SMILES processing for final canonized list of SMILES
canonized = ('[O-]C(=O)C(=O)CC([O-])=O',
'[O-]C(=O)C(O)(CC([O-])=O)CC([O-])=O',
'[O-]C(=O)C(CC([O-])=O)C(O)C([O-])=O',
'[O-]C(=O)CCC(=O)C([O-])=O',
'[O-]C(=O)CCC([O-])=O',
'[O-]C(=O)C=CC([O-])=O',
'[O-]C(=O)C(O)CC([O-])=O')


# In[8]:


pathway('M00004')


# In[197]:


paths = ['C(C1C(C(C(C(O1)O)O)O)O)OP(=O)(O)O',
 'C(C1C(C(C(C(=O)O1)O)O)O)OP(=O)(O)O',
 'C(C(C(C(C(C(=O)O)O)O)O)O)OP(=O)(O)O',
 'C(C(C(C(=O)CO)O)O)OP(=O)(O)O',
 'C(C(C(C(C=O)O)O)O)OP(=O)(O)O',         
 'C(C(C(C(C(C(=O)CO)O)O)O)O)OP(=O)(O)O',
 'C(C(C(C=O)O)O)OP(=O)(O)O',        
 'C(C1C(C(C(O1)(CO)O)O)O)OP(=O)(O)O',
 'C(C1C(C(C(C(O1)O)O)O)O)OP(=O)(O)O']


# In[200]:


#Re-canonization using indigo
ppp = []
indigo = Indigo()
for smile in paths:
    mol1 = indigo.loadMolecule(smile)
    mol1.aromatize()
    ppp.append(mol1.canonicalSmiles())


# In[201]:


ppp


# In[4]:


#pathway for Glycolysis (linear)
pathway('M00001')


# In[5]:


glyc = ('C(C1C(C(C(C(O1)O)O)O)O)O',
 'C(C1C(C(C(C(O1)O)O)O)O)OP(=O)(O)O',
 'C(C(C(C(C(=O)CO)O)O)O)OP(=O)(O)O',
 'C(C(=O)COP(=O)(O)O)O',
 'C(C(C=O)O)OP(=O)(O)O',
 'C(C(C(=O)OP(=O)(O)O)O)OP(=O)(O)O',
 'C(C(C(=O)O)O)OP(=O)(O)O',
 'C(C(C(=O)O)OP(=O)(O)O)O',
 'C=C(C(=O)O)OP(=O)(O)O',
 'CC(=O)C(=O)[O-]')


# In[9]:


#Re-canonization using indigo
glycolysis = []
indigo = Indigo()
for smile in glyc:
    mol1 = indigo.loadMolecule(smile)
    mol1.aromatize()
    glycolysis.append(mol1.canonicalSmiles())


# In[10]:


glycolysis


# In[ ]:




