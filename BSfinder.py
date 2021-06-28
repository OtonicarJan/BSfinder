import math
import argparse
import numpy as np
from Bio.Seq import Seq

infile = argparse.ArgumentParser(description='Set input files and destination file') 
#set the input genome file
infile.add_argument('genome', metavar='genome', help='Set the iniput genome', type=argparse.FileType('r'))
#set file with list of binding sites
infile.add_argument('binding_sites', metavar='list_of_binding_sites', help='file with aligned confirmed binding sites', type=argparse.FileType('r'))
#set the output file
infile.add_argument('dest_file', metavar='output_file', help='Set the name of your output file', type=argparse.FileType('w+'))
args = infile.parse_args()
Genome = args.genome
motifs = args.binding_sites.readlines()
output = args.dest_file

Motifs = []

for line in motifs:
    Motifs.append(line.strip())

#search window by window thorugh genome 
def Search(Motifs, Genome):
    tem = Genome.readlines()[1:]  
    genome = ""
    for line in tem: 
        genome = genome + line.strip()
    #genome = genome.replace("\n", "")
    #cutoff = percentile(Motifs, Genome)
    bindingsites = {}
    lenght = len(genome)
    k = len(Motifs[0])
    for i in range(lenght - k + 1):
        pattern = genome[i:i+k] 
        reverse = Seq(pattern).reverse_complement()
        score = BvH(pattern, Motifs) 
        scoreReverse = BvH(reverse, Motifs)
        if score < percentile(Motifs, Genome):
            bindingsites[i+1, i+k+1] = score 
        elif scoreReverse < percentile(Motifs, Genome): 
            bindingsites[i+1, i+k+1] = scoreReverse
    return bindingsites

#calculate the top 20% of hits based on scoring
def percentile(Motifs, Genome):
    profil = profile(Motifs)
    #Rseq = RSequence(Motifs, Genome)
    N = len(Motifs)
    k = len(Motifs[0])
    listofscores = []
    for pattern in Motifs:
        BvH = 0
        for bp in range(k):
            Pcons = 0
            Pobs = profil[pattern[bp]][bp]
            for base in "ACGT":
                if profil[base][bp] > Pcons:
                    Pcons = profil[base][bp]
            BvH += math.log(((Pcons + (1/N))/((Pobs + (1/N)))), math.e)
        listofscores.append(BvH)
    osem = np.percentile(listofscores, 20)
    return osem

#scoring function based on the calculated profile
def BvH(pattern, Motifs):
    profil = profile(Motifs)
    #Rseq = RSequence(Motifs, Genome)
    N = len(Motifs)
    k = len(pattern)
    #listofscores = []
    BvH = 0
    for bp in range(k):
        Pcons = 0
        Pobs = profil[pattern[bp]][bp]
        for base in "ACGT":
            if profil[base][bp] > Pcons:
                Pcons = profil[base][bp]
        BvH += math.log(((Pcons + (1/N))/((Pobs + (1/N)))), math.e)
        #listofscores.append(BvH)
    #osem = np.percentile(listofscores, 20)
    return BvH

#calculate the profile from aligned binding sites
def profile(Motifs):
    profile = Count(Motifs)
    t = len(Motifs)
    for letter in "ACGT":
        for i in range(14):
            profile[letter][i] = profile[letter][i]/t

    return profile

#count number of appearances of each letter at each position
def Count(Motifs):
    #for line in motifs:
        #alignedmotifs.append(line.strip())
    k = len(Motifs[0])
    count = {}
    for letter in "ACGT":
        count[letter] = []
        for j in range(k):
            count[letter].append(0)
    for string in Motifs:
        for i in range(k):
            if string[i] == "A":
                count["A"][i] += 1
            elif string[i] == "C":
                count["C"][i] += 1
            elif string[i] == "G":
                count["G"][i] += 1
            else:
                count["T"][i] += 1
    return count




print(Search(Motifs, Genome), file=output)
