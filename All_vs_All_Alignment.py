#!/usr/bin/env python
# coding: utf-8

# In[54]:


import numpy as np
import pubchempy
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
import math
import pickle


# In[55]:


#variable neccessary for list creation
CHAR = set(['(',')','{','}','[',']','=','@','.','#','+','-','/','\\''0','1','2','3','4','5','6','7','8','9'])
DOUB_ELEM = set(['Na','Mg','Al','Si','Cl','Ca','Cr','Mn', 'Fe','Co','Cu','Zn','Se','Mo','Cd','Sn','Br','As','Pb','Li'])
DOUBLE_TROUB = {'n':['M','Z','S'], 'o':['C','M']}
CHARs = set(['(',')','{','}','[',']','=','@','.','#','+','-','/','\\','0','1','2','3','4','5','6','7','8','9'])
elements = [
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", 
    "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", 
    "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", 
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", 
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn", 
    "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", 
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", 
    "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", 
    "Tl", "Pb", "Bi", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", 
    "Bk", "Cf", "Es", "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", 
    "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", 
    "Ts", "Og"
]


# In[56]:


#Opens scoring dictionary of all vs all
with open('score_all_v_all.pkl', 'rb') as f:
    score_dict = pickle.load(f)




# In[58]:

#adding threshold to dictionary
score_dict = {key: value + 1.524353523624236 for key, value in score_dict.items()}




# In[60]:


#function that changes a string of SMILES to a list
def parse_smiles(smile, reverse=False):
    smiles = []
    for i in range(len(smile) - 1):
        elem = smile[i]
        doub = smile[i:i+2]
        
        if elem.isdigit() or elem in CHAR:
            smiles.append(elem)

        elif doub in DOUB_ELEM:
            smiles.append(doub)

        elif elem in ('a','l','g','i','r','e','u','d','s','b'):
            continue

        elif elem in ('n','o',) and smile[i-1] in DOUBLE_TROUB[elem]:
            continue
        
        else:
            smiles.append(elem)

    # If the last character is not processed, process it
    if smile[-1] not in DOUB_ELEM:
        smiles.append(smile[-1])
        
    smiles = [i for i in smiles if i not in CHARs]

    if reverse:
        return "".join(smiles[::-1])
    else:
        return smiles


# In[61]:


#returns a list of gastieger charges in a smile, list includes the index of the SMILE
#where atoms are given a value, and characters are give '-'
def gast_smiles(smile):
    mol = Chem.MolFromSmiles(smile)
    AllChem.ComputeGasteigerCharges(mol)
    charges = [round(atom.GetDoubleProp('_GasteigerCharge'), 1) for atom in mol.GetAtoms()]
    
    sm1 = parse_smiles(smile)
    sm1_red = ['-' if i in CHARs else i for i in sm1]

    gast_smile = []
    count = 0
    for i in sm1_red:
        if i == '-':
            gast_smile.append(i)
        else:
            gast_smile.append(charges[count])
            count += 1

    return(gast_smile)


# In[62]:


#scoring function that returns score. characters are either matched or mismatched,
#and atoms are given a score based of the scoring matrix. alpha and beta are the 
#atoms or characters, and gast1 and gast2 are gasteiger charges (or '-' if characters)
def score(alpha, beta , gast1, gast2):
    alpha, beta = alpha.capitalize() , beta.capitalize()
    
    if alpha in elements and beta in elements:
        dif = round(abs(gast1-gast2), 1)
        return(score_dict[dif])
    
    elif alpha.isdigit() and beta.isdigit():
        return(match_score)
    
    elif alpha == beta:
        return(match_score)
    
    else: 
        return(mis_match_score)
    


# In[63]:

#Reintroduces non-atomic characters to aligned fingerprints
from collections import deque

def reconstruct_smiles(modified_smiles, original_smiles):
    modified_queue = deque(modified_smiles)
    original_queue = deque(original_smiles)

    result = []

    while modified_queue:
        char = modified_queue.popleft()

        # Directly append dashes
        if char == '-':
            result.append(char)
            continue

        # For atoms: match with original and append all characters till the next atom
        while original_queue:
            orig_char = original_queue.popleft()
            result.append(orig_char)
            if orig_char == char:
                break

    # Append any remaining characters from the original queue
    result.extend(original_queue)

    return ''.join(result)


# In[64]:


#alignment function based on Needelman-Wunsch
def align(seq1, seq2):
    
    seq1_orig = seq1
    seq2_orig = seq2
    
    gast1 , gast2 = gast_smiles(seq1) , gast_smiles(seq2)
    seq1, seq2 = parse_smiles(seq1) , parse_smiles(seq2)
    len1 , len2 = len(seq1) , len(seq2)
    
    M = np.zeros((len1+1, len2+1))
    X = np.zeros((len1+1, len2+1))
    Y = np.zeros((len1+1, len2+1))
    
    M[0, :] = -np.inf
    M[:, 0] = -np.inf
    M[0, 0] = 0

    X[:, 0] = gap_open_penalty + gap_extend_penalty * np.arange(len1+1)
    X[0, :] = -np.inf

    Y[0, :] = gap_open_penalty + gap_extend_penalty * np.arange(len2+1)
    Y[:, 0] = -np.inf
    
    #works
    for i in range(1, len1+1):
            for j in range(1, len2+1):
                 # Calculate match/mismatch score using the scoring matrix
                match = score(seq1[i-1], seq2[j-1], gast1[i-1] , gast2[j-1])
                M[i, j] = max(M[i-1, j-1] + match, X[i-1, j-1] + match, Y[i-1, j-1] + match)
                X[i, j] = max(M[i-1, j] + gap_open_penalty, X[i-1, j] + gap_extend_penalty, gap_open_penalty + gap_extend_penalty + Y[i-1,j])
                Y[i, j] = max(M[i, j-1] + gap_open_penalty, Y[i, j-1] + gap_extend_penalty)
                
                
    seq1_align, seq2_align = '', ''
    i, j = len1, len2
    while i > 0 or j > 0:
        if i > 0 and j > 0 and max(M[i, j], X[i, j], Y[i, j]) == M[i, j]:
            seq1_align += str(seq1[i-1])
            seq2_align += str(seq2[j-1])
            i -= 1
            j -= 1
        elif i > 0 and max(M[i, j], X[i, j], Y[i, j]) == X[i, j]:
            seq1_align += str(seq1[i-1])
            seq2_align += '-'
            i -= 1
        elif j > 0 and max(M[i, j], X[i, j], Y[i, j]) == Y[i, j]:
            seq1_align += '-'
            seq2_align += str(seq2[j-1])
            j -= 1
        else:
            break
    
    
    #return(seq1_align[::-1] , seq2_align[::-1])
    smiles_out1 = reconstruct_smiles(seq1_align[::-1], seq1_orig)
    smiles_out2 = reconstruct_smiles(seq2_align[::-1], seq2_orig)
    
    #Returns seq1 fingerprint alignemnt, seq2 fingerprint alignemnt, alignment score, seq1 reconstructed alignment, seq2 reconstructed alignemnt
    return(seq1_align[::-1] , seq2_align[::-1], M[len1,len2],smiles_out1,smiles_out2)
    


# In[65]:


#Levenshtein similairty for comparing the similarity of the truth and observed
def levenshtein_similarity(s1, s2):
    if len(s1) < len(s2):
        return levenshtein_similarity(s2, s1)

    # len(s1) >= len(s2)
    if len(s2) == 0:
        return len(s1)

    previous_row = range(len(s2) + 1)
    for i, c1 in enumerate(s1):
        current_row = [i + 1]
        for j, c2 in enumerate(s2):
            insertions = previous_row[j + 1] + 1
            deletions = current_row[j] + 1
            substitutions = previous_row[j] + (c1 != c2)
            current_row.append(min(insertions, deletions, substitutions))
        previous_row = current_row

    distance = previous_row[-1]
    max_len = max(len(s1), len(s2))
    return 1 - distance / max_len


# In[66]:

#exact string similiarity 
def string_similarity(truth, experimental):
    """
    Calculate the percentage similarity between two strings based on exact match. Extra characters in the 
    experimental string count against similarity.

    Parameters:
    - truth: The ground truth string.
    - experimental: The experimental string to compare against the truth.

    Returns:
    - Percentage similarity between the two strings.
    """
    max_len = max(len(truth), len(experimental))
    matches = sum(t == e for t, e in zip(truth, experimental))
    
    similarity_percentage = (matches / max_len)
    return similarity_percentage


# In[14]:

#Validation of Parameters


#canonized SMILES
canonized = ('[O-]C(=O)C(=O)CC([O-])=O',
'[O-]C(=O)C(O)(CC([O-])=O)CC([O-])=O',
'[O-]C(=O)C(CC([O-])=O)C(O)C([O-])=O',
'[O-]C(=O)CCC(=O)C([O-])=O',
'[O-]C(=O)CCC([O-])=O',
'[O-]C(=O)C=CC([O-])=O',
'[O-]C(=O)C(O)CC([O-])=O')


# In[15]:

#Hand aligned truth set
truth_set = ((('OCOCOCCOO', 'OCOCOCCOO'),
  ('OCOCO----CCOO', 'OCOCOCCOOCCOO'),
  ('OCO-CO---C-COO', 'OCOCC-COOCOCOO'),
  ('OCOCO-C-COO', 'OCOC-CCOCOO'),
  ('OCOCOCCOO', 'OCOC-CCOO'),
  ('OCOCOCCOO', 'OCOC-CCOO'),
  ('OCOCOCCOO', 'OCOCOCCOO')),
 [('OCOCOCCOOCCOO', 'OCOCO----CCOO'),
  ('OCOCOCCOOCCOO', 'OCOCOCCOOCCOO'),
  ('OCOCOCCOOC-COO', 'OCOC-CCOOCOCOO'),
  ('OCOCOCCOOC-COO', 'OCOC-C---COCOO'),
  ('OCOCOCCOOCCOO', 'OCOC----CCOO'),
  ('OCOCOCCOOCCOO', 'OCOC-----CCOO'),
  ('OCOCOCCOOCCOO', 'OCOCO----CCOO')],
 [('OCOCC-COOCOCOO', 'OCO-CO---C-COO'),
  ('OCOC-CCOOCOCOO', 'OCOCOCCOOC-COO'),
  ('OCOCCCOOCOCOO', 'OCOCCCOOCOCOO'),
  ('OCOCCCOOCOCOO', 'OCOCC---COCOO'),
  ('OCOCCCOOCOCOO', 'OCO-C---C-COO'),
  ('OCOCCCOOCOCOO', 'OCO-C---C-COO'),
  ('OCOC-CCOOCOCOO', 'OCOCOCCOO-----')],
 (('OCOC-CCOCOO', 'OCOCO-C-COO'),
  ('OCOC-C---COCOO', 'OCOCOCCOOC-COO'),
  ('OCOCC---COCOO', 'OCOCCCOOCOCOO'),
  ('OCOCCCOCOO', 'OCOCCCOCOO'),
  ('OCOCCCOCOO', 'OCOCC--COO'),
  ('OCOCCCOCOO', 'OCOCC--COO'),
  ('OCOC-CCOCOO', 'OCOCOC--COO')),
 (('OCOC-CCOO', 'OCOCOCCOO'),
  ('OCO---C--CCOO', 'OCOCOCCOOCCOO'),
  ('OCO-C---C-COO', 'OCOCCCOOCOCOO'),
  ('OCOCC--COO', 'OCOCCCOCOO'),
  ('OCOCCCOO', 'OCOCCCOO'),
  ('OCOCCCOO', 'OCOCCCOO'),
  ('OCOC-CCOO', 'OCOCOCCOO')),
 (('OCOC-CCOO', 'OCOCOCCOO'),
  ('OCO---C--CCOO', 'OCOCOCCOOCCOO'),
  ('OCO-C---C-COO', 'OCOCCCOOCOCOO'),
  ('OCOCC--COO', 'OCOCCCOCOO'),
  ('OCOCCCOO', 'OCOCCCOO'),
  ('OCOCCCOO', 'OCOCCCOO'),
  ('OCOC-CCOO', 'OCOCOCCOO')),
 (('OCOCOCCOO', 'OCOCOCCOO'),
  ('OCOCO----CCOO', 'OCOCOCCOOCCOO'),
  ('OCOCOCCOO-----', 'OCOC-CCOOCOCOO'),
  ('OCOCOC--COO', 'OCOC-CCOCOO'),
  ('OCOCOCCOO', 'OCOC-CCOO'),
  ('OCOCOCCOO', 'OCOC-CCOO'),
  ('OCOCOCCOO', 'OCOCOCCOO')))


# In[68]:


import pandas as pd
import time

start = time.time()
# Initialize an empty list to store the results
results = []
run = 1

# Loop over all combinations of parameters
for gap_open_penalty in range(-5, 6):
    for gap_extend_penalty in range(-5, 6):
        if gap_open_penalty <= gap_extend_penalty:
            similarities = []
            exact = []
            for j in range(0, len(canonized)):
                for i in range(0, len(canonized)):
                    # aligns
                    a = align(canonized[j], canonized[i])
                    # compares truth with levenshtein
                    similarity1 = levenshtein_similarity(a[0], truth_set[j][i][0])
                    similarity2 = levenshtein_similarity(a[1], truth_set[j][i][1])
                    similarities.extend((similarity1, similarity2))

                    # compares with exact
                    similarity_a = string_similarity(a[0], truth_set[j][i][0])
                    similarity_b = string_similarity(a[1], truth_set[j][i][1])
                    exact.extend((similarity_a, similarity_b))

            levenshtein_average = sum(similarities) / len(similarities)
            exact_average = sum(exact) / len(exact)
            results.append([gap_open_penalty, gap_extend_penalty, levenshtein_average, exact_average])

            
            
df = pd.DataFrame(results, columns=[ 'Gap Open Penalty', 'Gap Extend Penalty', 'Levenshtein Average', 'Exact Average'])

stop = time.time()

print(stop-start)

print(df)


# In[69]:


df_sorted = df.sort_values(by='Levenshtein Average', ascending=False)


# In[70]:


df_sorted.head(10)


# In[71]:


df_sorted = df.sort_values(by='Exact Average', ascending=False)


# In[72]:


df_sorted.head(10)


# In[73]:


import seaborn as sns
import matplotlib.pyplot as plt


# In[19]:




# In[26]:


# List of columns except 'Levenshtein Average'
columns = [col for col in df.columns if col != 'Levenshtein Average']

# Iterate over the columns and create a scatter plot with a trendline for each
for column in columns:
    plt.figure(figsize=(12, 6))
    sns.regplot(x=df[column], y=df['Levenshtein Average'], line_kws={"color": "red"})
    plt.title(f'Levenshtein Average vs {column} in All vs All Scoring with Trendline')
    plt.xlabel(column)
    plt.ylabel('Levenshtein Average')
    plt.show()


# In[27]:


# List of columns except 'Except Average'
columns = [col for col in df.columns if col != 'Exact Average']

# Iterate over the columns and create a scatter plot with a trendline for each
for column in columns:
    plt.figure(figsize=(12, 6))
    sns.regplot(x=df[column], y=df['Exact Average'], line_kws={"color": "red"})
    plt.title(f'Exact Average vs {column} in All vs All Scoring with Trendline')
    plt.xlabel(column)
    plt.ylabel('Exact Average')
    plt.show()


# Application of Pentose Phosphate Pathway

# In[28]:


ppp = ['OC1OC(COP(O)(O)=O)C(O)C(O)C1O',
 'OC1C(COP(O)(O)=O)OC(=O)C(O)C1O',
 'OC(COP(O)(O)=O)C(O)C(O)C(O)C(O)=O',
 'OCC(=O)C(O)C(O)COP(O)(O)=O',
 'OC(C(O)COP(O)(O)=O)C(O)C=O',
 'OCC(=O)C(O)C(O)C(O)C(O)COP(O)(O)=O',
 'OC(C=O)C(O)COP(O)(O)=O',
 'OCC1(O)OC(COP(O)(O)=O)C(O)C1O',
 'OC1OC(COP(O)(O)=O)C(O)C(O)C1O']


# In[29]:


gap_open = df_sorted['Gap Open Penalty'].head(10).tolist()
gap_extend = df_sorted['Gap Extend Penalty'].head(10).tolist()


# In[30]:


gap_open[0],gap_extend[0]


# In[31]:


result = []
n = len(ppp)

for i in range(n):
    # Creating a new list starting from the current item and cycling through the list
    cycled_list = ppp[i:] + ppp[:i]
    
    score_lst = []
    for j in range(0, len(cycled_list)):
        a = align(cycled_list[0], cycled_list[j])
        score_lst.append(a[2])
    
    result.append(score_lst)


# In[32]:


min_length = min(len(lst) for lst in result)

# Initialize a list to hold the sums
sums = [0] * min_length

# Sum up values at each index
for lst in result:
    for i in range(min_length):
        sums[i] += lst[i]

# Calculate the average for each index
averages = [sum_val / len(result) for sum_val in sums]


# In[52]:


data = []

for x in range(0,10):
    gap_open_penalty = gap_open[x]
    gap_extend_penalty = gap_extend[x]
    
    result = []
    
    n = len(ppp)

    for i in range(n):
        # Creating a new list starting from the current item and cycling through the list
        cycled_list = ppp[i:] + ppp[:i]

        score_lst = []
        for j in range(0, len(cycled_list)):
            a = align(cycled_list[0], cycled_list[j])
            score_lst.append(a[2])

        result.append(score_lst)
        
    min_length = min(len(lst) for lst in result)

    # Initialize a list to hold the sums
    sums = [0] * min_length

    # Sum up values at each index
    for lst in result:
        for i in range(min_length):
            sums[i] += lst[i]

    # Calculate the average for each index
    averages = [sum_val / len(result) for sum_val in sums]
    
    data.append(averages)


# In[53]:


for i, sublist in enumerate(data):
    plt.plot(sublist, label=f'List {i+1}')

# Adding titles and labels
plt.title('Top Paramters: Alignment of the Pentose Phosphate Pathway using All vs All Scoring')
plt.xlabel('Distance from Original Position')
plt.ylabel('Average Alignment Score')
plt.savefig('valid_all_v_all.png',bbox_inches='tight')
plt.show()


# In[38]:


plt.savefig('valid_all_v_all.png')


# In[37]:


import matplotlib.pyplot as plt
indexes = range(len(data[9]))

# Creating the plot
plt.figure(figsize=(8, 6))
plt.plot(indexes, data[0], marker='o')

# Adding title and labels
plt.title('Plot of Values by Index')
plt.xlabel('Index')
plt.ylabel('Value')

# Show the plot
plt.show()


# In[49]:


data = []

for x in range(0,10):
    gap_open_penalty = gap_open[x]
    gap_extend_penalty = gap_extend[x]
    
    result = []
    
    n = len(canonized)

    for i in range(n):
        # Creating a new list starting from the current item and cycling through the list
        cycled_list = canonized[i:] + canonized[:i]

        score_lst = []
        for j in range(0, len(cycled_list)):
            a = align(cycled_list[0], cycled_list[j])
            score_lst.append(a[2])

        result.append(score_lst)
        
    min_length = min(len(lst) for lst in result)

    # Initialize a list to hold the sums
    sums = [0] * min_length

    # Sum up values at each index
    for lst in result:
        for i in range(min_length):
            sums[i] += lst[i]

    # Calculate the average for each index
    averages = [sum_val / len(result) for sum_val in sums]
    
    data.append(averages)


# In[51]:


for i, sublist in enumerate(data):
    plt.plot(sublist, label=f'List {i+1}')

# Adding titles and labels
plt.title('Top 10 Paramters: Alignment of the Krebs Cycle using All vs All')
plt.xlabel('Distance from Original Position')
plt.ylabel('Average Alignment Score')
plt.savefig('valid_all_v_all_tca.png',bbox_inches='tight')
plt.show()


# In[ ]:





# In[67]:


gly = ('OCC1OC(O)C(O)C(O)C1O',
 'OC1OC(COP(O)(O)=O)C(O)C(O)C1O',
 'OCC(=O)C(O)C(O)C(O)COP(O)(O)=O',
 'OCC(=O)COP(O)(O)=O',
 'OC(COP(O)(O)=O)C=O',
 'OP(O)(=O)OCC(O)C(=O)OP(O)(O)=O',
 'OC(COP(O)(O)=O)C(O)=O',
 'OCC(OP(O)(O)=O)C(O)=O',
 'C=C(OP(O)(O)=O)C(O)=O',
 'CC(=O)C([O-])=O')


# In[79]:


data = []

for x in range(0,10):
    gap_open_penalty = gap_open[x]
    gap_extend_penalty = gap_extend[x]
    
    result = []
    
    for i in range(0, len(gly)):
        a = align(gly[0],gly[i])
        result.append(a[2])
        
    data.append(result)


# In[84]:


for i, sublist in enumerate(data):
    plt.plot(sublist, label=f'List {i+1}')

# Adding titles and labels
plt.title('Top 10 Paramters: Alignment of the Glycolsis (Linear Pathway) using All vs All')
plt.xlabel('Distance from Original Position')
plt.ylabel(' Alignment Score')
plt.savefig('valid_all_v_all_glyc.png',bbox_inches='tight')
plt.show()


# In[86]:






