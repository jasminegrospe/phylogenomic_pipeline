#import needed modules
import os
import sys
import glob
import subprocess
import matplotlib.pyplot as plt
from Bio import SeqIO
from Bio.Seq import Seq
from Bio import Phylo


#STEP 1: Set input/output paths

in_dir = "/shared/forsythe/BB485/Week06/Brass_CDS_seqs/"
out_dir = "/scratch/grospeja/BB485/Week06/out_phy_pipeline/"


#STEP 2: Multiple sequence alignment (MAFFT)

#get a list of all .fasta files from the input folder
fasta_files = glob.glob(in_dir+"*fasta")

#run MAFFT alignment 
for file in fasta_files:
    new_file_path = file.replace(in_dir, out_dir)

    #run MAFFT command
    aln_cmd = 'mafft --auto --quiet '+file+' > '+new_file_path
    os.system(aln_cmd)


#STEP 3: Run IQ-TREE for tree inference

#get a list of aligned .fasta files from the output folder
aln_files = glob.glob(out_dir+"*fasta")

#run IQtree on each aligned file
for aln in aln_files:
    tree_command = f"iqtree -s {aln} -m TEST -nt 16"

    #run IQtree command
    os.system(tree_command)


#STEP 4: Topology testing

#make list to collect all topology types
topo_list = []

#read and analyze each treefile from the output folder
tree_files = glob.glob(out_dir + "*.treefile")

for tree in tree_files:
    #read in the tree and store as Phylo object
    temp_tree = Phylo.read(tree, "newick")

    #loop through the tips to find Escherichia outgroup
    for tip in temp_tree.get_terminals():
        if "Es_" in tip.name:
            es_tip = tip
            #stop once Escherichia is found
            break
        
    #root the tree using the Escherichia outgroup
    temp_tree.root_with_outgroup(es_tip)
        
    #get all terminal (aka tips) branches 
    all_terminal_branches = temp_tree.get_terminals()
        
    #assign tree tips to their corresponding species
    for t in all_terminal_branches:
        if "Bs_" in t.name: 
            Bs_temp=t #Brassica
        elif "Cr_" in t.name:
            Cr_temp=t #Capsella
        elif "At_" in t.name:
            At_temp=t #Arabidopsis
        else:
            out_temp=t #outgroup tip
            
    #define possible monophyletic pairs
    P1_and_P2=[Bs_temp, Cr_temp]
    P1_and_P3=[Bs_temp, At_temp]
    P2_and_P3=[Cr_temp, At_temp]
        
    #determine which pair is monophyletic
    if bool(temp_tree.is_monophyletic(P1_and_P2)):
        topo_str = "12top"
    elif bool(temp_tree.is_monophyletic(P1_and_P3)):
        topo_str = "13top"
    elif bool(temp_tree.is_monophyletic(P2_and_P3)):
        topo_str = "23top"
    else:
        topo_str = "Unknown"

    #store topology results
    topo_list.append(topo_str)


#STEP 5: Print summary of topologies

print("\n--- Topology Summary ---")
print(f"A. tha and C. rub sister: {topo_list.count('23top')}")
print(f"B. str and C. rub sister: {topo_list.count('12top')}")
print(f"A. tha and B. str sister: {topo_list.count('13top')}")
print(f"Unknown topologies: {topo_list.count('Unknown')}")
print("\nPhylogenomic pipeline completed successfully!")


#STEP 6: Make a Bar Chart

#count topologies
labels = ["A. tha and C. rub sister", "B. str and C. rub sister", "A. tha and B. str sister", "Unknown"]
counts = [
    topo_list.count("23top"),
    topo_list.count("12top"),
    topo_list.count("13top"),
    topo_list.count("Unknown")]

#set custom colors
colors = ["#552583", "#FDB927", "#3B9AE1", "#4A4A4A"]

#create the bar chart
plt.figure(figsize=(8, 5))
plt.bar(labels, counts, color=colors)
plt.ylabel("Number of Trees")
plt.xlabel("Topology Type")
plt.title("Topology Count Summary")
plt.tight_layout()

#save figure
plt.savefig("topology_bar_chart.png")

#exit script
sys.exit()