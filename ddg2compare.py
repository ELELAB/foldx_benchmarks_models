'''
Created on May 18, 2018

@author: Tycho
'''
import numpy as np
import os
import Bio
import matplotlib.pyplot as plt

def main():
    np.set_printoptions(threshold=np.nan)
    mutation_list = "GAVLIMFWPSTCYNQDEKRH"
    dirs = ["Z:\\mutatex_project\\P62942\\mutatex_data\\100_experimental\\results\\mutation_ddgs\\final_averages\\",
            "Z:\\mutatex_project\\P62942\\mutatex_data\\100_model\\results\\mutation_ddgs\\final_averages\\",
            "Z:\\mutatex_project\\P62942\\mutatex_data\\90-70_model\\results\\mutation_ddgs\\final_averages\\",
            "Z:\\mutatex_project\\P62942\\mutatex_data\\70-50_model\\results\\mutation_ddgs\\final_averages\\",
            "Z:\\mutatex_project\\P62942\\mutatex_data\\50-30_model\\results\\mutation_ddgs\\final_averages\\",
            "Z:\\mutatex_project\\P62942\\mutatex_data\\30-20_model\\results\\mutation_ddgs\\final_averages\\",
            ]

#     per_mutation_ddgs = get_per_mutation_ddgs(dirs, mutation_list)
    per_position_ddgs = get_per_position_ddg(dirs, mutation_list)
#     per_position_rmsd = get_per_position_rmsd(dirs, mutation_list, per_position_ddgs)
#     per_position_self = get_self_mutation()
#     per_position_mean_ddg = get_per_position_mean_ddg(dirs)
#     write_per_position_rmsd(dirs, per_position_rmsd, mutation_list)
    self_mutations_native = get_self_mutation_native(dirs, per_position_ddgs, mutation_list)
    self_mutations_new = get_self_mutation_new(dirs)
    
    amino_acids = [(k, float(k[2:])) for k in self_mutations_native.keys()]
    amino_acids = [k[0] for k in sorted(amino_acids, key=lambda x: x[1], reverse=False)]
    
    plot_per_position_rmsd_hist(dirs, mutation_list, self_mutations_new, 
                                "per_group_self_mutations",
                                 x_limit = (-0.2,0.2),
                                 y_limit = (0, 60), 
                                bins = np.arange(-0.2, 0.2, 0.0108), 
                                 x_label = "self-mut_ddg")
#     write_per_position_mean_ddgs(dirs, self_mutations_new, "self_mutation_ddgs")
    


def get_per_mutation_ddgs(dirs, mutation_list):
#     ddg = np.empty([len(mutation_list)], dtype=float)
    results = []
    for data_dir in dirs:
        ddgs = {}
        print data_dir.split("\\")[-5]
        for file_name in os.listdir(data_dir):
            with open(data_dir + file_name, 'r') as f:
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
    for data_dir in dirs:
        group = {}
        print data_dir.split("\\")[-5]
        for file_name in os.listdir(data_dir):
            with open(data_dir + file_name, 'r') as f:
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
    for data_dir in dirs:
        print data_dir.split("\\")[-5]
        for file_name in os.listdir(data_dir):
            with open(data_dir + file_name, 'r') as f:
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

def get_per_position_rmsd(dirs, mutation_list, per_position_ddgs, index_reference_values = 0):
    rmsds = {}
    for i in range(len(per_position_ddgs)):
        for position, mutations in per_position_ddgs[i].iteritems():
            if position not in rmsds:
                rmsds[position] = [np.sqrt(np.mean(abs(mutations - per_position_ddgs[index_reference_values][position])**2))]
            else:
                rmsds[position].append(np.sqrt(np.mean(abs(mutations - per_position_ddgs[index_reference_values][position])**2))) 
    return rmsds

def get_self_mutation_native(dirs, per_position_ddgs, mutation_list):
    groups = [dir.split("\\")[-5] for dir in dirs]
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
    #Over complicated way of replacing the sub proejct folder in every dir path.
    dirs = ["\\".join(dir.split("\\")[:3] + ["mutatex_self_mutations"] + dir.split("\\")[3+1:-2]) for dir in dirs]

    out = {}
    ddg = np.empty([20], dtype=float)
    for dir in dirs:
        data_folder = [i for i in os.listdir(dir) if not i.startswith("final")][0]      #Get the folder that contains selfmutation.dat
        dir += "\\" + data_folder
        
        with open(dir + "\\selfmutation_energies.dat", 'r') as f:
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
        
        
def plot_per_position_rmsd_hist(dirs, mutation_list, per_position_rmsd, out_file, 
                                x_limit = (-5,10), y_limit = (0,40), bins = np.arange(0, 20 + 0.5, 0.5), x_label = "ddg-rmsd"):
    groups = [dir.split("\\")[-5] for dir in dirs]
    per_group_data = []
    for group in groups:
        per_group_data.append((group,[]))
        
    for aa, rmsds in per_position_rmsd.iteritems():
        for i in range(len(rmsds)):
            per_group_data[i][1].append(rmsds[i])
            
    plt.figure(figsize=(20,10))
    for i, group in enumerate(per_group_data):
        plt.subplot(2, 3, i+1)
        plt.title(dirs[i].split("\\")[-5])
        plt.hist(group[1], bins=bins, rwidth = 0.8)
        plt.xlim(x_limit)
        plt.ylim(y_limit)
        plt.xlabel(x_label)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95], pad=2)
    plt.subplots_adjust(wspace = 0.1)
    plt.suptitle("%s\nbin-size: %.4f" %(out_file, np.diff(bins)[0]), fontsize=16, fontweight='bold')
    plt.savefig(out_file + ".pdf")
    plt.clf()
                
def plot_per_mutation_ddg_hist(dirs, mutation_list, ddg_data):
    for aa in mutation_list:
        plt.figure(figsize=(20,10))
        for i in range(len(ddg_data)):
            arr = [float(j) for j in ddg_data[i][aa]]
            plt.subplot(2, 3, i+1)
            plt.title(dirs[i].split("\\")[-5])
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
            ax[i].set_title(dirs[i].split("\\")[-5])
            ax[i].scatter(x,y, label = dirs[i].split("\\")[-5], s = 5)
            ax[i].plot((-5,20),(-5,20), 'r--')
#             ax[i].legend(loc='upper right')
            ax[i].axis('equal')
            ax[i].set_aspect("equal", 'box')
            ax[i].set_xlim(-5,20)
            ax[i].set_ylim(-5,20)
            ax[i].set_xlabel(dirs[0].split("\\")[-5])
            ax[i].set_ylabel(dirs[i].split("\\")[-5])    
                             
        fig.tight_layout(rect=[0, 0.1, 1, 0.95])
#         fig.subplots_adjust(wspace = 0.01, hspace = 1)
        fig.suptitle(aa, fontsize=22, fontweight='bold')
        fig.savefig(aa + "_ddg_per_mutation_scatter.pdf")
        fig.clf()  


def write_per_mutation_ddgs(dirs, mutation_list, per_mutation_ddgs):
    with open("per_mutation_ddgs.csv", 'w') as out_file:
        for aa in mutation_list:
            out_file.write(aa + "\n")
            for i in range(len(per_mutation_ddgs)):
                out_file.write(dirs[i].split("\\")[-5] + "," + ",".join(per_mutation_ddgs[i][aa]) + "\n")      

def write_per_mutation_corrcoef(dirs, mutation_list, per_mutation_ddgs):
    with open("per_mutation_correlation.csv", 'w') as out_file:
        out_file.write(",")
        for dir in dirs:
            out_file.write(dir.split("\\")[-5]+ ",")
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
            out_file.write(dir.split("\\")[-5].split("_")[0]+ ",")
        out_file.write("\n")
        for key, val in mean_ddgs.iteritems():
            val = [str(i) for i in val]
            out_file.write(key + "," + ",".join(val) + '\n')

def write_per_position_rmsd(dirs, rmsds, out_file):
    with open(out_file + ".csv", 'w') as out_file:
        out_file.write(",")
        for dir in dirs[1:]:
            out_file.write(dir.split("\\")[-5].split("_")[0]+ ",")
        out_file.write("\n")
        for key, val in rmsds.iteritems():
            val = [str(i) for i in val]
            val.pop(0)
            out_file.write(key + "," + ",".join(val) + '\n')


              
if __name__ == '__main__':
    main()