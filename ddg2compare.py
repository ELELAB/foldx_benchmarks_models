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
    dirs = ["Z:\\mutatex_project\\P62942\\mutatex_data\\hundred_experimental\\results\\mutation_ddgs\\final_averages\\"]
    dirs.append("Z:\\mutatex_project\\P62942\\mutatex_data\\hundred_model\\results\\mutation_ddgs\\final_averages\\")
    dirs.append("Z:\\mutatex_project\\P62942\\mutatex_data\\ninety_seventy_model\\results\\mutation_ddgs\\final_averages\\")
    dirs.append("Z:\\mutatex_project\\P62942\\mutatex_data\\seventy_fifty_model\\results\\mutation_ddgs\\final_averages\\")
    dirs.append("Z:\\mutatex_project\\P62942\\mutatex_data\\fifty_thirty_model\\results\\mutation_ddgs\\final_averages\\")
    dirs.append("Z:\\mutatex_project\\P62942\\mutatex_data\\thirty_twenty_model\\results\\mutation_ddgs\\final_averages\\")

#     ddgs_per_mutation = get_per_mutation_ddgs(dirs, mutation_list)
#     plot_per_mutation_hist(dirs, mutation_list, ddgs_per_mutation)
#     plot_per_mutation_scatter(dirs, mutation_list, ddgs_per_mutation)

    rms_ddgs = get_rmsd_ddgs(dirs)
    write_rms_ddgs(dirs, rms_ddgs)
                
def plot_per_mutation_hist(dirs, mutation_list, ddg_data):
    for aa in mutation_list:
        plt.figure(figsize=(20,10))
        for i in range(len(ddg_data)):
            arr = [float(j) for j in ddg_data[i][aa]]
            plt.subplot(2, 3, i+1)
            plt.title(dirs[i].split("\\")[-5])
            plt.hist(arr, bins=np.arange(min(arr), max(arr) + 0.5, 0.5), rwidth = 0.8)
            plt.xlim(-5,5)
            plt.ylim(0,30)
        plt.tight_layout(rect=[0, 0.03, 1, 0.95])
        plt.subplots_adjust(wspace = 0.1)
        plt.suptitle(aa, fontsize=22, fontweight='bold')
        plt.savefig(aa + "_ddg_per_mutation_hist.pdf")
        plt.clf()



def plot_per_mutation_scatter(dirs, mutation_list, ddg_data):
    for aa in mutation_list:
        plt.figure(figsize=(20,10))
        for i in range(len(ddg_data)):
            x = [float(j) for j in ddg_data[0][aa]]
            y = [float(j) for j in ddg_data[i][aa]]
#             plt.subplot(2, 3, i+1)
#             plt.title(dirs[i].split("\\")[-5])
#             plt.hist(arr, bins=np.arange(min(arr), max(arr) + 0.5, 0.5), rwidth = 0.8)

            plt.scatter(x,y, label = dirs[i].split("\\")[-5])
            plt.legend(loc='upper right')
            
            plt.xlim(-2,10)
            plt.ylim(-5,15)
        plt.tight_layout(rect=[0, 0.1, 1, 0.95])
        plt.subplots_adjust(wspace = 0.1)
        plt.suptitle(aa, fontsize=22, fontweight='bold')
        plt.savefig(aa + "_ddg_per_mutation_scatter.pdf")
        plt.clf()  

def get_rmsd_ddgs(dirs):
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
                rms = np.sqrt(np.mean(ddg**2))
                if file_name not in out:
                    out[file_name] = [rms]
                else:
                    out[file_name].append(rms)
    return out

def get_per_mutation_ddgs(dirs, mutation_list):
    ddg = np.empty([20], dtype=float)
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

      
    
def get_mean_ddgs(dirs):
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


def write_mean_ddg_corrcoef(dirs, mutation_list, mean_ddgs):
    #Calculate per AA mean ddg corr coeficient
    with open("coorelations.txt", 'w') as coors:
        for aa in mutation_list:
            ddgs = [[] for x in range(len(dirs))]
            for key, val in mean_ddgs.iteritems():
                if key.startswith(aa):
                    for i in range(len(ddgs)):
                        ddgs[i].append(val[i])
            print aa  
            coor = np.corrcoef(ddgs)[0]
            coor = [str(i) for i in coor]
            coors.write("\t".join(coor) + "\n")

def write_mean_ddgs(dirs, mean_ddgs):
    with open("mean_ddgs.csv", 'w') as out_file:
        for key, val in mean_ddgs.iteritems():
            val = [str(i) for i in val]
            out_file.write(key + "," + ",".join(val) + '\n')

def write_rms_ddgs(dirs, rms_ddgs):
    with open("rms_ddgs.csv", 'w') as out_file:
        out_file.write(",")
        for dir in dirs:
            out_file.write(dir.split("\\")[-5]+ ",")
        out_file.write("\n")
        for key, val in rms_ddgs.iteritems():
            val = [str(i) for i in val]
            out_file.write(key + "," + ",".join(val) + '\n')


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
                
                
if __name__ == '__main__':
    main()