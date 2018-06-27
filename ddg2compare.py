'''
Created on May 18, 2018

@author: Tycho
'''
import numpy as np
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis.rms import RMSD
from MDAnalysis.analysis import contacts
from MDAnalysis.tests.datafiles import PSF,DCD
from Bio import PDB
import subprocess

def main():
    np.set_printoptions(threshold=np.nan)
    mutation_list = "GAVLIMFWPSTCYNQDEKRH"
    dirs = ["mutatex_project/P62942/mutatex_data/100_experimental",
            "mutatex_project/P62942/mutatex_data/100_model",
            "mutatex_project/P62942/mutatex_data/90-70_model",
            "mutatex_project/P62942/mutatex_data/70-50_model",
            "mutatex_project/P62942/mutatex_data/50-30_model",
            "mutatex_project/P62942/mutatex_data/30-20_model",
            ]
    
    
    checked_pdbs = []
    repaired_pdbs = []
    for dir in dirs:
        dir += "/repair"
        dir = dir + "/" + os.listdir(dir)[0]
        checked_pdbs.append(dir + "/" + [file for file in os.listdir(dir) if file.endswith("checked_sorted.pdb")][0])
        repaired_pdbs.append(dir + "/" + [file for file in os.listdir(dir) if file.endswith("Repair_sorted.pdb")][0])
    
#     per_mutation_ddgs = get_per_mutation_ddgs(dirs, mutation_list)
    per_position_ddgs = get_per_position_ddg(dirs, mutation_list)
    per_position_rmsd = get_per_position_ddg_rmsd(dirs, mutation_list, per_position_ddgs)
#     per_position_fonc = get_per_positoin_fonc(repaired_pdbs, checked_pdbs[0])
#     per_position_self = get_self_mutation()
#     per_position_mean_ddg = get_per_position_mean_ddg(dirs)
#     write_per_position_rmsd(dirs, per_position_rmsd, "ddg_rmsds")
    
    per_position_ASAs, secondary_structures = get_per_position_dssp_data(checked_pdbs)
    
    
############## Self mutations ################
#     self_mutations_new = get_self_mutation_new(dirs)    
#     plot_per_position_rmsd_hist(dirs, mutation_list, self_mutations_new,
#                                 "per_group_self_mutations", x_limit = (-0.2,0.2), y_limit = (0, 60), 
#                                 bins = np.arange(-0.2, 0.2, 0.0108), x_label = "self-mut_ddg")
#     write_per_position_mean_ddgs(dirs, self_mutations_new, "self_mutation_ddgs")
    

############## Single Side chain root mean square deviations ################
#     per_position_sideChain_rmsd_checked = get_per_position_sideChain_rmsd(checked_pdbs, checked_pdbs[0])
#     per_position_sideChain_rmsd_repaired = get_per_position_sideChain_rmsd(repaired_pdbs, checked_pdbs[0])
#        
#     plot_per_position_scatterplot(dirs, per_position_sideChain_rmsd_checked, per_position_rmsd,
#                                         mutation_list, x_limit = (-0.05,10), y_limit = (-0.05,10),
#                                         out_file = "exp-model-sidechain-rmsd_ddg-rmsd_scatterplot", 
#                                         x_label = "side-chain-rmsd", y_label = "ddg_rmsd", color_data = per_position_ASAs)
#     plot_per_position_scatterplot(dirs, per_position_sideChain_rmsd_repaired, per_position_rmsd,
#                                         mutation_list, x_limit = (-0.05,10), y_limit = (-0.05,10),
#                                         out_file = "exp-repaired-sidechain-rmsd_ddg-rmsd_scatterplot", 
#                                         x_label = "side-chain-rmsd", y_label = "ddg_rmsd")
#     plot_per_position_scatterplot(dirs, per_position_sideChain_rmsd_repaired, per_position_rmsd,
#                                         mutation_list, x_limit = (-0.05,10), y_limit = (-0.05,10),
#                                         out_file = "exp-repaired-sidechain-rmsd_exp-model-sidechain-rmsd_ddg-rmsd_scatterplot", 
#                                         x_label = "side-chain-rmsd", y_label = "ddg_rmsd", add_data = per_group_data_checked)
#
#     plot_per_position_rmsd_hist(dirs, mutation_list, per_position_sideChain_rmsd_checked ,
#                                 "per_group_sidechain_rmsds_checked_hist", x_limit = (-0.05,7), y_limit = (0, 40), 
#                                 bins = np.arange(0, 7, 0.2), x_label = "side_chain_rmsd")
#     plot_per_position_rmsd_hist(dirs, mutation_list, per_position_sideChain_rmsd_repaired, 
#                                 "per_group_sidechain_rmsds_repaired_hist", x_limit = (-0.05,7), y_limit = (0, 40), 
#                                 bins = np.arange(0, 7, 0.2), x_label = "side_chain_rmsd")
#
#     write_per_position_mean_ddgs(dirs, per_position_sideChain_rmsd_checked, "side_chain_rmsd_checked")
#     write_per_position_mean_ddgs(dirs, per_position_sideChain_rmsd_repaired, "side_chain_rmsd_repaired")


###############Fraction of Native Contacts################
#     per_position_foncs_checked = get_per_positoin_fonc(checked_pdbs, checked_pdbs[0])
#     per_position_foncs_repaired = get_per_positoin_fonc(repaired_pdbs, checked_pdbs[0])
#     
#     
#     plot_per_position_scatterplot(dirs, per_position_foncs_checked, per_position_rmsd,
#                                   mutation_list, x_limit = (-0.05,1.05), y_limit = (-0.05,7),
#                                   out_file = "exp-model-fonc_ddg-rmsd_scatterplot", 
#                                   x_label = "fonc", y_label = "ddg_rmsd")
#     plot_per_position_scatterplot(dirs, per_position_foncs_repaired, per_position_rmsd,
#                                   mutation_list, x_limit = (-0.05,1.05), y_limit = (-0.05,7),
#                                   out_file = "exp-repaired-fonc_ddg-rmsd_scatterplot", 
#                                   x_label = "fonc", y_label = "ddg_rmsd")
#     plot_per_position_scatterplot(dirs, per_position_foncs_repaired, per_position_rmsd,
#                                   mutation_list, x_limit = (-0.05,1.05), y_limit = (-0.05,7),
#                                   out_file = "exp-repaired-fonc_exp-model-fonc_ddg-rmsd_scatterplot", 
#                                   x_label = "fonc", y_label = "ddg_rmsd", x2 = per_position_foncs_checked)
#         
#     plot_per_position_rmsd_hist(dirs, mutation_list, per_position_foncs_checked,
#                                 "per_group_fonc_checked_hist", x_limit = (-0.05,1.1), y_limit = (0, 40), 
#                                 bins = np.arange(0, 1.05, 0.03888), x_label = "side_chain_rmsd")      
#     plot_per_position_rmsd_hist(dirs, mutation_list, per_position_foncs_repaired, 
#                                 "per_group_fonc_repaired_hist", x_limit = (-0.05,1.1), y_limit = (0, 40), 
#                                 bins = np.arange(0, 1.05, 0.03888), x_label = "side_chain_rmsd")
#
#     write_per_position_mean_ddgs(dirs, per_position_foncs_checked, "fonc_checked")
#     write_per_position_mean_ddgs(dirs, per_position_foncs_repaired, "fonc_repaired")


def get_per_mutation_ddgs(dirs, mutation_list):
#     ddg = np.empty([len(mutation_list)], dtype=float)
    results = []
    for dir in dirs:
        dir += "/results/mutation_ddgs/final_averages/"
        ddgs = {}
        print dir.split("/")[-5]
        for file_name in os.listdir(dir):
            with open(dir + file_name, 'r') as f:
                f.next()
                for i, ddg in enumerate(f):
                    ddg = ddg.split(" ")[0]
                    if mutation_list[i] in ddgs:
                        ddgs[mutation_list[i]].append(ddg)
                    else:
                        ddgs[mutation_list[i]] = [ddg]
        results.append(ddgs)
    return results

def get_per_position_ddg(dirs, mutation_list):
    out = []
    for dir in dirs:
        dir += "/results/mutation_ddgs/final_averages/"
        group = {}
        print dir.split("/")[-5]
        for file_name in os.listdir(dir):
            with open(dir + file_name, 'r') as f:
                f.next()
                ddgs = np.empty([len(mutation_list)], dtype=float)
                for i, line in enumerate(f):
                    ddgs[i] = float(line.split(" ")[0])
                group[file_name] = ddgs
        out.append(group)
    return out

def get_per_position_mean_ddg(dirs):
    out = {}
    ddg = np.empty([20], dtype=float)
    for dir in dirs:
        dir += "/results/mutation_ddgs/final_averages/"
        print dir.split("/")[-5]
        for file_name in os.listdir(dir):
            with open(dir + file_name, 'r') as f:
                f.next()
                j = 0
                for line in f:
                    ddg[j] = float(line.split(" ")[0])
                    j += 1
                if file_name not in out:
                    out[file_name] = [np.mean(ddg)]
                else:
                    out[file_name].append(np.mean(ddg))
    return out

def get_per_position_ddg_rmsd(dirs, mutation_list, per_position_ddgs, index_reference_values = 0):
    rmsds = {}
    for i in range(len(per_position_ddgs)):
        for position, mutations in per_position_ddgs[i].iteritems():
            if position not in rmsds:
                rmsds[position] = [np.sqrt(np.mean(abs(mutations - per_position_ddgs[index_reference_values][position])**2))]
            else:
                rmsds[position].append(np.sqrt(np.mean(abs(mutations - per_position_ddgs[index_reference_values][position])**2))) 
    return rmsds

def get_per_position_sideChain_rmsd(pdbs, experimental):
    rmsds = {}
    for pdb in pdbs:
        structure = PDB.PDBParser().get_structure("structurw", pdb)
        length = (len(structure[0]['A']))
        resnums = range(1,length+1)
        
        print pdb
        structure1 = mda.Universe(experimental)
        structure2 = mda.Universe(pdb)
        
        for s in resnums:
            fit_g1  = "resnum %d and (name C or name CA or name N)" %s
            calc_g1 = "resnum %d and not backbone and not name OXT and prop mass >2" %s
            
            
            aa_name = PDB.Polypeptide.three_to_one(structure1.select_atoms(fit_g1).residues[0].resname) + "A" + str(s)
            
            RMSD(structure1, structure2, select=fit_g1, groupselections=[calc_g1])
            rmsd_obj = RMSD(structure1, structure2, select=fit_g1, groupselections=[calc_g1])
            out = rmsd_obj.run()
            
            
            if aa_name not in rmsds:
                rmsds[aa_name] = [out.rmsd[0,3]]
            else:
                rmsds[aa_name].append(out.rmsd[0,3])
    return rmsds

def get_per_positoin_fonc(pdbs, experimental):
    ''' Calculate the fraction of native contacts (FoNC) between a list of
    pdb structures and a reference experimental structure'''
    
    def radius_cut_q2(r, r0):
            return contacts.radius_cut_q(r, r0, 4.5)

    foncs = {}
    for pdb in pdbs:
        structure = PDB.PDBParser().get_structure("structurw", pdb)
        length = (len(structure[0]['A']))
        resnums = range(1,length+1)
        
        print pdb
        ur = mda.Universe(experimental)
        u  = mda.Universe(pdb)
        for r in resnums:
            residue_str ='resnum %d and not backbone and prop mass >2' % r
            prot_str ='not (resnum %d and not backbone and prop mass >2)' % r
            sel_ur_res = ur.select_atoms(residue_str)
            sel_ur_prot = ur.select_atoms(prot_str)
            
            aa_name = PDB.Polypeptide.three_to_one(u.residues[r-1].resname) + "A" + str(r)
            
            ca=contacts.Contacts(u, selection=(residue_str, prot_str), refgroup=(sel_ur_res, sel_ur_prot), radius=4.5, method=radius_cut_q2)
            ca.run()
            
            if aa_name not in foncs:
                foncs[aa_name] = [ca.timeseries[0,1]]
            else:
                foncs[aa_name].append(ca.timeseries[0,1])
    return foncs

def get_per_position_dssp_data(pdbs):
    '''Load position specific protein solvent accessibility values from DSSP 2.2.18 output.
    Solvent accessibility values (ASA) are normelized by their amino-acid specific theoretical maximum ASA values.
    Returns a bfactor replacement dictionary.
    Requires DSSP to be installed on the system.
    '''
    ASA_max = {'CYS': 167.0, 'ASP': 193.0, 'SER': 155.0, 'GLN': 225.0, 'LYS': 236.0, 'PRO': 159.0, 'THR': 172.0, 'PHE': 240.0, 'ALA': 129.0, 'HIS': 224.0, 'GLY': 104.0, 'ILE': 197.0, 'GLU': 223.0, 'LEU': 201.0, 'ARG': 274.0, 'TRP': 285.0, 'VAL': 174.0, 'ASN': 195.0, 'TYR': 263.0, 'MET': 224.0}
    ASAs = {}   #solvent accessibility values
    sss = {}    #Secondary StructureS
    for file_name in pdbs:
        try:
            dssp_output = subprocess.check_output(["dssp-2.2.8", "-i", file_name]).split("\n")
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
                aa_name = aa + "A" + str(position)
                
                if aa_name not in ASAs:
                    ASAs[aa_name] = [ASA]
                else:
                    ASAs[aa_name].append(ASA)
                if aa_name not in sss:
                    sss[aa_name] = [structure]
                else:
                    sss[aa_name].append(structure)
            if line.startswith("  #"):
                bool = True
    return ASAs, sss

def get_self_mutation_native(dirs, per_position_ddgs, mutation_list):
    groups = [dir.split("/")[-1] for dir in dirs]
    out = {}
    for i, group in enumerate(groups):
        self_muts = []
        for aa, ddgs in per_position_ddgs[i].iteritems():
            if aa not in out:
                out[aa] = [ddgs[mutation_list.index(aa[0])]]
            else:
                out[aa].append(ddgs[mutation_list.index(aa[0])])
    return out
        
def get_self_mutation_new(dirs):
    dirs = [dir.replace("mutatex_data", "mutatex_self_mutations") for dir in dirs]
    out = {}
    ddg = np.empty([20], dtype=float)
    for dir in dirs:
        dir += "/results/mutation_ddgs"
        data_folder = [i for i in os.listdir(dir) if i.endswith("Repair")][0]      #Get the folder that contains selfmutation.dat
        dir += "/" + data_folder
        print dir
        quit()
        
        with open(dir + "/selfmutation_energies.dat", 'r') as f:
            f.next()
            for line in f:
                line = line.split()
                
                aa = line[0]
                ddg = line[1]
                
                if aa not in out:
                    out[aa] = [float(ddg)]
                else:
                    out[aa].append(float(ddg))
    return out
        
  
def plot_per_position_rmsd_hist(dirs, mutation_list, per_position_rmsd,
                                out_file, x_limit = (-5,10), y_limit = (0,40),  
                                bins = np.arange(0, 20 + 0.5, 0.5), x_label = "ddg-rmsd"):
    
    groups = [dir.split("/")[-1] for dir in dirs]
    per_group_data = [[] for group in groups]
#     for group in groups:
#         per_group_data.append((group,[]))
        
    for aa, rmsds in per_position_rmsd.iteritems():
        for i in range(len(rmsds)):
            if not np.isnan(rmsds[i]):
                per_group_data[i].append(rmsds[i])
    
    plt.figure(figsize=(20,10))
    for i, group in enumerate(per_group_data):
        plt.subplot(2, 3, i+1)
        plt.title(dirs[i].split("/")[-1])
        plt.hist(group, bins=bins, rwidth = 0.8)
        plt.xlim(x_limit)
        plt.ylim(y_limit)
        plt.xlabel(x_label)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95], pad=2)
    plt.subplots_adjust(wspace = 0.1)
    plt.suptitle("%s\nbin-size: %.4f" %(out_file, np.diff(bins)[0]), fontsize=16, fontweight='bold')
    plt.savefig(out_file + ".pdf")
    plt.clf()

def plot_per_position_scatterplot(dirs, x, y, mutation_list,
                                  out_file = "rmsd-rmsd_scatterplot", x_limit = (0,10), y_limit = (0,10), 
                                  x_label = "x-values", y_label = "y-values",
                                  x2 = None, x2_name = "x2",  color_data = None):
    '''Create a scatterplot of two sets of per_position amino-acid data. 
    (eg. get_per_position_mean_ddg, get_per_position_ddg_rmsd, get_per_position_sideChain_rmsd or get_per_positoin_fonc)
    '''
    
    groups = [dir.split("/")[-1] for dir in dirs]
    per_group_x = [[] for group in groups]
    per_group_x2 = [[] for group in groups]
    per_group_y = [[] for group in groups]
    per_group_color_data = [[] for group in groups]
    
    for aa in x.keys():
        for i in range(len(x[aa])):
            if not np.isnan(x[aa][i]) and not np.isnan(y[aa][i]):
                per_group_x[i].append(x[aa][i])
                per_group_y[i].append(y[aa][i])
                if(x2 is not None and not np.isnan(x2[aa][i])):
                    per_group_x2[i].append(x2[aa][i])
                if(color_data is not None and not np.isnan(color_data[aa][i])):
                    per_group_color_data[i].append(color_data[aa][i])
    
    print color_data
    return

    fig, axs = plt.subplots(2, 3, figsize=(20,10))
    ax = axs.flat
    for i, group in enumerate(groups):
        per_group_x = np.asarray(per_group_x, float)
        per_group_y = np.asarray(per_group_y, float)
        ax[i].scatter(per_group_x[i], per_group_y[i], s = 12, label = out_file.split("_")[0], c = color_data[i], cmap = "autumn_r")
        
        if(x2 is not None):
            ax[i].scatter(per_group_x2[i], per_group_y[i], s = 12, color = "red", label = x2_name)
            ax[i].legend(loc='upper right');
        ax[i].set_xlim(x_limit)
        ax[i].set_ylim(y_limit)
        
        ax[i].set_title(group)
        ax[i].set_xlabel(x_label)
        ax[i].set_ylabel(y_label)
                
    fig.suptitle(out_file, fontsize=22, fontweight='bold')
    fig.savefig(out_file + ".pdf")

def plot_per_mutation_ddg_hist(dirs, mutation_list, ddg_data):
    for aa in mutation_list:
        plt.figure(figsize=(20,10))
        for i in range(len(ddg_data)):
            arr = [float(j) for j in ddg_data[i][aa]]
            plt.subplot(2, 3, i+1)
            plt.title(dirs[i].split("/")[-1])
            plt.hist(arr, bins=np.arange(min(arr), max(arr) + 0.5, 0.5), rwidth = 0.8)
            plt.xlim(-5,5)
            plt.ylim(0,30)
            plt.xlabel("ddg")
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.subplots_adjust(wspace = 0.1)
        plt.suptitle(aa, fontsize=22, fontweight='bold')
        plt.savefig(aa + "_ddg_per_mutation_hist.pdf")
        plt.clf()

def plot_per_mutation_ddg_scatter(dirs, mutation_list, per_residue_ddgs):
    for aa in mutation_list:
        fig, axs = plt.subplots(2, 3, figsize=(10,10))
        ax = axs.flat
        for i in range(len(per_residue_ddgs)):
            x = [float(j) for j in per_residue_ddgs[0][aa]]
            y = [float(j) for j in per_residue_ddgs[i][aa]]
            ax[i].set_title(dirs[i].split("/")[-1])
            ax[i].scatter(x,y, label = dirs[i].split("/")[-1], s = 5)
            ax[i].plot((-5,20),(-5,20), 'r--')
#             ax[i].legend(loc='upper right')
            ax[i].axis('equal')
            ax[i].set_aspect("equal", 'box')
            ax[i].set_xlim(-5,20)
            ax[i].set_ylim(-5,20)
            ax[i].set_xlabel(dirs[0].split("/")[-5])
            ax[i].set_ylabel(dirs[i].split("/")[-1])    
                             
        fig.tight_layout(rect=[0, 0.1, 1, 0.95])
#         fig.subplots_adjust(wspace = 0.01, hspace = 1)
        fig.suptitle(aa, fontsize=22, fontweight='bold')
        fig.savefig(aa + "_ddg_per_mutation_scatter.pdf")
        fig.clf()

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
    for i, group in enumerate(sorted(groups, key=lambda x: x[1], reverse=True)):
        group = group[0]
        structure_rmsds = []
        counts = np.zeros(len(structures))
        for j, structure in enumerate(structures):
            rmsds = []
            for key, val in ss[group].iteritems():
                if val == structure:
                    rmsds.append(ddg_rmsds[group][key])
                    counts[j] += 1
            structure_rmsds.append(rmsds)

        print group + "\t" + " " + str(counts)
#         print structures
        plt.subplot(2, 3, i+1)
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

def write_per_mutation_ddgs(dirs, mutation_list, per_mutation_ddgs):
    with open("per_mutation_ddgs.csv", 'w') as out_file:
        for aa in mutation_list:
            out_file.write(aa + "\n")
            for i in range(len(per_mutation_ddgs)):
                out_file.write(dirs[i].split("/")[-1] + "," + ",".join(per_mutation_ddgs[i][aa]) + "\n")

def write_per_mutation_corrcoef(dirs, mutation_list, per_mutation_ddgs):
    with open("per_mutation_correlation.csv", 'w') as out_file:
        out_file.write(",")
        for dir in dirs:
            out_file.write(dir.split("/")[-1]+ ",")
        out_file.write("\n")
        for aa in mutation_list:
            out_file.write(aa + ",")
            coor_matrix = []
            for i in range(len(per_mutation_ddgs)):
                coor_matrix.append([float(j) for j in per_mutation_ddgs[i][aa]])
            coors = np.corrcoef(coor_matrix)[0]
            coors = [str(i) for i in coors]
            out_file.write(",".join(coors) + "\n")

def write_per_position_mean_ddgs(dirs, mean_ddgs, out_file):
    with open(out_file + ".csv", 'w') as out_file:
        out_file.write(",")
        for dir in dirs:
            out_file.write(dir.split("/")[-1].split("_")[0]+ ",")
        out_file.write("\n")
        for key, val in mean_ddgs.iteritems():
            val = [str(i) for i in val]
            out_file.write(key + "," + ",".join(val) + '\n')

def write_per_position_rmsd(dirs, rmsds, out_file):
    with open(out_file + ".csv", 'w') as out_file:
        out_file.write(",")
        for dir in dirs[1:]:
            out_file.write(dir.split("/")[-1].split("_")[0]+ ",")
        out_file.write("\n")
        for key, val in rmsds.iteritems():
            val = [str(i) for i in val]
            val.pop(0)
            out_file.write(key + "," + ",".join(val) + '\n')

              
if __name__ == '__main__':
    main()