'''
Created on Jun 8, 2018

@author: Tycho

Get the theoratical maximum residue’s solvent accessibility (ASA)- (Maximum Allowed Solvent
 
Accessibilites of Residues in proteins
Matthew Z. Tien1, Austin G. Meyer2,3, Dariya K. Sydykova3,4,5, Stephanie J. Spielman3,4,5, Claus O. Wilke3,4,5*)
'''

ASA_max = {'CYS': '167.0', 'ASP': '193.0', 'SER': '155.0', 'GLN': '225.0', 'LYS': '236.0', 'ILE': '197.0', 'PRO': '159.0', 'THR': '172.0', 'PHE': '240.0', 'ASN': '195.0', 'GLY': '104.0', 'HIS': '224.0', 'LEU': '201.0', 'ARG': '274.0', 'TRP': '285.0', 'ALA': '129.0', 'VAL': '174.0', 'GLU': '223.0', 'TYR': '263.0', 'MET': '224.0'}

def three_to_ASA(s):
    '''standard amino acid three letter code to thoeretical ASA-max.
    
    >>> whole_to_three('GLY')
    '104.0'
    >>> whole_to_three('TRP')
    '285.0'
    '''
    return float(ASA_max[s])