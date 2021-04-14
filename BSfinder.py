import math
import argparse
import numpy as np
from Bio.Seq import Seq

infile = argparse.ArgumentParser(description='Set input files and destination file')
infile.add_argument('genome', metavar='genome', help='Set input csv depth file', type=argparse.FileType('r'))
#set the width of the window
infile.add_argument('binding_sites', metavar='binding sites', help='Determine the width of the window', type=argparse.FileType('r'))
#set the output file
infile.add_argument('dest_file', metavar='output file', help='Set the name of your output file', type=argparse.FileType('w+'))
args = infile.parse_args()
Genome = args.genome
motifs = args.binding_sites.readlines()
output = args.dest_file

Motifs = []

for line in motifs:
    Motifs.append(line.strip())


def Search(Motifs, Genome):
    genome = Genome.read().replace("\n", "")
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

def profile(Motifs):
    profile = Count(Motifs)
    t = len(Motifs)
    for letter in "ACGT":
        for i in range(14):
            profile[letter][i] = profile[letter][i]/t

    return profile

def Consensus(Motifs):
    counting = Count(Motifs)
    #t = len(alignedmotifs)
    k = len(Motifs[0])
    consensus = ""
    for i in range(k):
        m = 0
        frequentSymbol = ""
        for letter in "ACGT":
            if counting[letter][i] > m:
                m = counting[letter][i]
                frequentSymbol = letter
        consensus += frequentSymbol
    return consensus

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
