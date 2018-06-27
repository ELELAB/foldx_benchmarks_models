'''
Created on Jun 8, 2018

@author: Tycho

Retrieve the 3-letter amino-acid identifier from the whole name.
'''
from Bio.PDB import Polypeptide

wholeToThree = {'asparagine': 'ASN', 'alanine': 'ALA', 'glycine': 'GLY', 'phenylalanine': 'PHE', 'aspartate': 'ASP', 'isoleucine': 'ILE', 'glutamate': 'GLU', 'lysine': 'LYS', 'leucine': 'LEU', 'methionine': 'MET', 'glutamine': 'GLN', 'proline': 'PRO', 'serine': 'SER', 'tyrosine': 'TYR', 'arginine': 'ARG', 'glutamic acid': 'GLU', 'valine': 'VAL', 'tryptophan': 'TRP', 'threonine': 'THR', 'histidine': 'HIS', 'cysteine': 'CYS', 'aspartic acid': 'ASP'}
threeToWhole = {'CYS': 'cysteine', 'HIS': 'histidine', 'SER': 'serine', 'GLN': 'glutamine', 'LYS': 'lysine', 'TRP': 'tryptophan', 'PRO': 'proline', 'THR': 'threonine', 'PHE': 'phenylalanine', 'ALA': 'alanine', 'ASN': 'asparagine', 'GLY': 'glycine', 'ILE': 'isoleucine', 'LEU': 'leucine', 'ARG': 'arginine', 'ASP': 'aspartate', 'VAL': 'valine', 'GLU': 'glutamate', 'TYR': 'tyrosine', 'MET': 'methionine'}

def whole_to_three(s):
    '''Amino-acid name to three letter code.
    
    >>> whole_to_three('Aspartic acid')
    'ASP'
    >>> whole_to_three('Aspartate')
    'ASP'

    For non-standard amino acids, you get a KeyError:

    >>> whole_to_three('dehydralanine')
    Traceback (most recent call last):
       ...
    KeyError: 'dehydralanine'
    '''
    return wholeToThree[s.lower()]

def three_to_whole(s):
    '''Amino-acid three letter code to name.
    
    >>> three_to_whole('ASP')
    'aspartate'
    >>> three_to_whole('SER')
    'serine'

    For non-standard amino acids, you get a KeyError:

    >>> three_to_one('MSE')
    Traceback (most recent call last):
       ...
    KeyError: 'MSE'
    '''
    return threeToWhole[s.upper()]

def whole_to_one(s):
    '''Amino-acid name to one letter code.
    
    >>> whole_to_one('Alanine')
    'A'
    >>> whole_to_one('Tryptophan')
    'W'

    For non-standard amino acids, you get a KeyError:

    >>> whole_to_one('o-phosphotyrosine')
    Traceback (most recent call last):
       ...
    KeyError: 'o-phosphotyrosine'
    '''
    return Polypeptide.three_to_one(wholeToThree[s.lower()])

def one_to_whole(s):
    '''Amino-acid one letter code to name.
    
    >>> one_to_whole('A')
    'Alanine'
    >>> one_to_whole('W')
    'Tryptophan'

    For non-standard amino acids, you get a KeyError:

    >>> one_to_whole('y')
    Traceback (most recent call last):
       ...
    KeyError: 'y'
    '''
    return threeToWhole[Polypeptide.one_to_three(s)]

def get_whole_to_three():
    return wholeToThree