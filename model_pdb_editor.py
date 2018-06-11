'''
Created on Jun 1, 2018

@author: Tycho
'''

import os
import subprocess
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from Bio import PDB

def main():
    pdb_input_dir = "mutatex_project//P62942//model_pdbs"
    rmsd_infile = "rms_ddgs.csv"    
    groups = []
    for file_name in os.listdir(pdb_input_dir):
        groups.append([file_name.split("_")[1], int(file_name.split("_")[1].split("-")[0])])
    
#     print load_mutation_data(input_values)
    dssps = load_dssp_data(pdb_input_dir)
    print ""
    
    ddg_rmsds = load_mutation_data(rmsd_infile)
    print ""




def box_thing_thing():
    print "DO SOMETHING?"
    print "can't really remember what"



def scatter_rmsd_dssp(ddg_rmsds, dssps, groups):
    #Scatter ddg-rmsd values against ASA values
    fig, axs = plt.subplots(2, 3, figsize=(20,10))
    ax = axs.flat
    i = 0
    for group in sorted(groups, key=lambda x: x[1], reverse=True):
        group = group[0]
        print group
        x = []
        y = []
        for aa, asa in dssps[group].iteritems():
            rmsd = ddg_rmsds[group][aa]
            x.append(asa)
            y.append(rmsd)
        ax[i].set_title(group)
        ax[i].scatter(x,y, label = group, s = 8)
#         ax[i].axis('equal')
        ax[i].set_xlim(0,1)
        ax[i].set_ylim(0,20)
#         diffs = (np.diff(ax[i].get_ylim())[0], np.diff(ax[i].get_xlim())[0])
#         ax[i].set_aspect(np.diff(ax[i].get_xlim())[0] / np.diff(ax[i].get_ylim())[0])
        ax[i].set_xlabel("ASA-value")
        ax[i].set_ylabel("ddg-rmsd")
        i+=1        
    fig.delaxes(ax[-1])
    fig.suptitle("scatter ASA--ddgRmsd", fontsize=22, fontweight='bold')
    fig.savefig("ASA--ddgRmsd_scatter.pdf")
    fig.clf()
            
    
def pdb_bfactor_replace(pdb_input_dir, values, out_name):
    '''
    Replace the input-pdb-file(s)'s bfactor column values with position specific ddg-based
    discriptive values for the corresponding amino-acid position. One value per amino-acid position.
    These can be the mean ddg per amino-acid or root mean square deviation with control ddg values.
    
    Input is a dictionary with keys: group/pdb name
    values: another dict containing amino-acid position and the corresponding value. eg: {"GA1" : 2.5, "WA2" : 3.5}
    '''
    pdb_parser = PDB.PDBParser()
    for file_name in os.listdir(pdb_input_dir):
        group = file_name.split("_")[1]
        print group

        structure = pdb_parser.get_structure('s', pdb_input_dir + "//" + file_name)
        structure = structure.get_chains().next()
        for aa, rmsd in values[group].iteritems():
            for atom in structure[int(aa[2:])]:
                atom.set_bfactor(float(rmsd))
        pdb_io = PDB.PDBIO()
        pdb_io.set_structure(structure)
        pdb_io.save("modified_" + out_name + "_" + file_name)

def load_dssp_data(pdb_input_dir):
    '''
    Load position specific protein solvent accessibility values from DSSP 2.2.18 output.
    Returns a bfactor replacement dictionary.
    Requires DSSP to be installed on the system.
    '''
    ASA_max = {'CYS': 167.0, 'ASP': 193.0, 'SER': 155.0, 'GLN': 225.0, 'LYS': 236.0, 'PRO': 159.0, 'THR': 172.0, 'PHE': 240.0, 'ALA': 129.0, 'HIS': 224.0, 'GLY': 104.0, 'ILE': 197.0, 'GLU': 223.0, 'LEU': 201.0, 'ARG': 274.0, 'TRP': 285.0, 'VAL': 174.0, 'ASN': 195.0, 'TYR': 263.0, 'MET': 224.0}
    pdb_parser = PDB.PDBParser()
    dssp = {}
    for file_name in os.listdir(pdb_input_dir):
        group = file_name.split("_")[1]
        dssp[group] = {}
        print group
        try:
            dssp_output = subprocess.check_output(["dssp-2.2.8", "-i", "mutatex_project//P62942//model_pdbs//" + file_name]).split("\n")
        except OSError:
            print "Error: Could not load dssp-2.2.8. Make sure dssp is installed."
        except Exception as e:
            print e.output
            
        dssp_output.pop(-1)
        
        bool = False
        for line in dssp_output:
            if bool:
                values = [line[:5], line[5:10], line[10:12], line[12:14], line[14:17], line[17:25], line[25:29], line[29:33], line[33:34], line[34:38]]
                values = [value.strip() for value in values]
                position = values[1]
                aa = values[3]
                ASA = float(values[9])
                ASA = ASA / ASA_max[PDB.Polypeptide.one_to_three(aa)]
                dssp[group][aa + "A" + str(position)] = ASA
            if line.startswith("  #"):
                bool = True
    return dssp
    
def load_mutation_data(in_file):
    '''
    Load mutation-ddg based values from a comma separated csv file.

    The files-layout should contain a header: position,group1,group2,group3...
    Following lines should contain the corresponding comma-separated values per column.
    '''
    
    ddgs = {}
    with open(in_file, 'r') as f:
        groups = f.next().split(",")
        groups.pop(0)
        groups.pop(-1)
        groups = [group.split("_")[0] for group in groups]
        for group in groups:
            print group
            ddgs[group] = {}
        for line in f:
            line = line.split("\n")[0].split(",")
            line.pop(1)
            for i in range(len(line))[1:]:
                ddgs[groups[i]][line[0]] = float(line[i])
    return ddgs
        






if __name__ == '__main__':
    main()