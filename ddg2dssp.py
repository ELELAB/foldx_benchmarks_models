'''
Created on Jun 15, 2018

@author: Tycho
'''

import os
import subprocess
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from Bio import PDB
from scipy import stats

def main():
    pdb_input_dir = "mutatex_project/P62942/model_pdbs"
    rmsd_infile = "ddg_rmsds.csv"    
    groups = []
    for file_name in os.listdir(pdb_input_dir):
        groups.append([file_name.split("_")[1], int(file_name.split("_")[1].split("-")[0])])
    groups = [group[0] for group in sorted(groups, key=lambda x: x[1], reverse=True)]
    
    ASAs, sss = load_dssp_data(pdb_input_dir)
    ddg_rmsds = load_rmsd_data(rmsd_infile)
    
    os.listdir(pdb_input_dir)
    scatter_rmsd_dssp(ddg_rmsds, ASAs, groups)
    
    write_dssp_data(ASAs, groups)

def load_rmsd_data(in_file):
    '''Load mutation-ddg based values from a comma separated csv file.

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
            ddgs[group] = {}
        for line in f:
            line = line.split("\n")[0].split(",")
            for i in range(len(line))[1:]:
                ddgs[groups[i-1]][line[0]] = float(line[i])
    return ddgs

def load_dssp_data(pdb_input_dir):
    '''Load position specific protein solvent accessibility values from DSSP 2.2.18 output.
    Solvent accessibility values (ASA) are normelized by their amino-acid specific theoretical maximum ASA values.
    Returns a bfactor replacement dictionary.
    Requires DSSP to be installed on the system.
    '''
    ASA_max = {'CYS': 167.0, 'ASP': 193.0, 'SER': 155.0, 'GLN': 225.0, 'LYS': 236.0, 'PRO': 159.0, 'THR': 172.0, 'PHE': 240.0, 'ALA': 129.0, 'HIS': 224.0, 'GLY': 104.0, 'ILE': 197.0, 'GLU': 223.0, 'LEU': 201.0, 'ARG': 274.0, 'TRP': 285.0, 'VAL': 174.0, 'ASN': 195.0, 'TYR': 263.0, 'MET': 224.0}
    ASAs = {}   #solvent accessibility values
    sss = {}    #Secondary StructureS
    for file_name in os.listdir(pdb_input_dir):
        group = file_name.split("_")[1]
        ASAs[group] = {}
        sss[group] = {}
        try:
            dssp_output = subprocess.check_output(["dssp-2.2.8", "-i", "mutatex_project/P62942/model_pdbs/" + file_name]).split("\n")
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
                if values[4] == "G" or values[4] == "I": 
                    structure = "H"
                elif values[4] == "E":
                    structure = "B"
                elif values[4] == "S":
                    structure = "T"    
                elif values[4] == "":
                    structure = "--"
                else:
                    structure = values[4]
                    
                ASAs[group][aa + "A" + str(position)] = ASA
                sss[group][aa + "A" + str(position)] = structure
            if line.startswith("  #"):
                bool = True
    return ASAs, sss

def box_secondary_structure(ddg_rmsds, ss, groups):
    '''Create a box plot of the ddg-rmsd data corresponding to dssp predicted secondary structures.
    x-axis: Protein secondary structures
    y-axis: ddg-rmsd values
    '''

    structures = []
    for group in groups:
        for key, val in ss[group].iteritems():
            structures.append(val)
    structures = list(set(structures))
    print "ss" + "\t" + " " + str(structures)
    plt.figure(figsize=(20,10))
    for j, group in enumerate(sorted(groups, key=lambda x: x[1], reverse=True)):
        group = group[0]
        structure_rmsds = []
        counts = np.zeros(len(structures))
        for i, structure in enumerate(structures):
            rmsds = []
            for key, val in ss[group].iteritems():
                if val == structure:
                    rmsds.append(ddg_rmsds[group][key])
                    counts[i] += 1
            structure_rmsds.append(rmsds)

        print group + "\t" + " " + str(counts)
#         print structures
        plt.subplot(2, 3, j+1)
        plt.title(group)
        plt.boxplot(structure_rmsds)
        xt = structures
        xticks = np.arange(1,len(structure_rmsds)+1)
        rotation_xt=0
        plt.xticks(xticks,xt,rotation=rotation_xt)
        plt.xlabel("secondary structure")
        plt.ylabel("ddg-rmsd")
        plt.ylim(0,12)
        plt.suptitle("ddg-rmsd boxplot per dssp secondary structure", fontsize = 22)
    plt.savefig("secondary_structure_rmsd_boxplot.pdf")
                           
def scatter_rmsd_dssp(ddg_rmsds, ASAs, groups):
    '''Create a scatterplot of the ddg-rmsd data against the ASAs solvent accessebility data.
    x-axis: asa
    y-axis: rmsd
    '''
    fig, axs = plt.subplots(2, 3, figsize=(20,10))
    ax = axs.flat
    i = 0
    for group in groups:
        x = []
        y = []
        for aa, asa in ASAs[group].iteritems():
            rmsd = ddg_rmsds[group][aa]
            x.append(asa)
            y.append(rmsd)
        x = np.asarray(x, float)
        y = np.asarray(y, float)
        
        print x
        return
        xmin = -0.05
        xmax = 1
        ymin = 0
        ymax = 12
        X, Y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
        positions = np.vstack([X.ravel(), Y.ravel()])
        values = np.vstack([x, y])
        kernel = stats.gaussian_kde(values)
        Z = np.reshape(kernel(positions).T, X.shape)
        
        ax[i].imshow(np.rot90(Z), cmap=plt.cm.get_cmap("gist_earth_r"), extent=[xmin, xmax, ymin, ymax])
        ax[i].scatter(x,y, label = group, s = 12)
        ax[i].set_xlim(-0.05,1)
        ax[i].set_ylim(0,12)
        
        ax[i].set_title(group)
        ax[i].set_aspect(np.diff(ax[i].get_xlim())[0] / np.diff(ax[i].get_ylim())[0])
        ax[i].set_xlabel("ASA-normalized")
        ax[i].set_ylabel("ddg-rmsd")
        i+=1
                
    fig.delaxes(ax[-1])
    fig.suptitle("scatter ASA--ddgRmsd", fontsize=22, fontweight='bold')
    fig.savefig("ASA_ddgRmsd_scatter.pdf")
    fig.clf()
    
def write_dssp_data(dssp, groups):
    with open("model_ASAs.csv", 'w') as out_file:
        out_file.write(",")
        for group in groups:
            out_file.write(group + ",")
        out_file.write("\n")
        AAs = dssp[groups[0]].keys()
        for AA in AAs:
            outline = AA + ","
            for group in groups:
                outline += str(dssp[group][AA]) + ","
            out_file.write(outline + "\n")
        


if __name__ == '__main__':
    main()