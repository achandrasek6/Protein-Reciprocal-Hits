"""
Best‐Reciprocal Homolog Finder

Calculates memoized global alignment scores (using BLOSUM62 and a constant gap penalty) between human and chicken protein sequences,
then identifies and reports best‐reciprocal‐hit gene pairs as candidate one‐to‐one orthologs.
"""

import sys
sys.setrecursionlimit(100000)
from humanChickenProteins import *
from blosum62 import *

# Memoized version of alignment score function to speed up computation
def memoAlignScore(S1, S2, gap, substitutionMatrix, memo):
    '''Memoized version of alignment score function to speed up computation'''
    if S1 == '': return len(S2) * gap # if S1 is empty will need to insert gap for every character in S2
    elif S2 == '': return len(S1) * gap # inverse of condition above
    else:
        if (S1, S2) in memo: return memo[(S1, S2)] # check cache
        else:
            option1 = substitutionMatrix[(S1[0], S2[0])] + memoAlignScore(S1[1:], S2[1:], gap, substitutionMatrix, memo)
            # Checking substitution matrix value at first position of each sequence; recursively moves to next index of both
            option2 = gap + memoAlignScore(S1[1:], S2, gap, substitutionMatrix, memo) # adds gap to S2 by only moving S1 to next index
            option3 = gap + memoAlignScore(S1, S2[1:], gap, substitutionMatrix, memo) # inverse of above
            solution = max(option1, option2, option3) # highest score is the winner
            memo[(S1, S2)] = solution # caching solution
            return solution

#print(memoAlignScore('','',-9,blosum62,{}))
# 0
#print(memoAlignScore('','FGTSK',-9,blosum62,{}))
# -45
#print(memoAlignScore('IVEKGYY','AVEYY',-9,blosum62,{}))
# 4
#print(memoAlignScore('CIEAFGTSKQKRALNSRRMNAVGNDIVSTAVTKAAADVIDAKGVTALIQDVAQD','RDLPIWTSVDWKSLPATEIFNKAFSQGSDEAMYDYMAVYKKSCPQTRR',-9,blosum62,{}))
# -48

# Getting best reciprocal hits

# The data
chromosome, startPosition, endPosition, proteinSequence=geneD['c6']
#print(chromosome)
# chr15
#print(startPosition)
# 791273
#print(endPosition)
# 791857
#print(proteinSequence[:25])
# MNSGILFLSLLGFLPSVIPTCPLPC

# Obtain the alignment score between all proteins in two species
def allScores(geneList1, geneList2):
    '''Function to obtain the alignment score between all proteins in two species'''
    results = {} # initializing results dictionary

    # comparing each gene in first list to each gene in second list (nested loop)
    for gene1 in geneList1:
        for gene2 in geneList2:
            score = memoAlignScore(geneD[gene1][3], geneD[gene2][3], -9, blosum62, {}) # alignment score
            results[(gene1, gene2)] = score # storing scores in dict
    return results

allScoresD=allScores(sampleHumanGeneList,sampleChickenGeneList) # test arguements
#print(len(allScoresD.keys()))
# 20
#print(allScoresD[('h4', 'c8')])
# -134

# Finding closest match
def closestMatch(geneName,allScoresD):
    '''Returns the gene from the other species which is most similar has the highest alignment score.'''
    highest_score = float('-inf')  # Initialize with negative infinity
    closest_gene = None

    # Iterate through the keys of the dictionary
    for key in allScoresD.keys():
        if geneName in key:  # Check if geneName is part of the key tuple
            gene1, gene2 = key
            score = allScoresD[key]  # Get the score for this pair

            # Determine the partner gene
            partner_gene = gene1 if gene2 == geneName else gene2

            if score > highest_score:
                highest_score = score
                closest_gene = partner_gene

    return closest_gene

#print(closestMatch('c19',allScoresD))
# h4
#print(closestMatch('h17',allScoresD))
# c8

# Finding closest reciprocal hit
def printBRH(geneName,allScoresD):
    '''Finds and prints the best reciprocal hit for geneName'''
    # Call closestMatch to get the best matching gene in the other species
    best_match = closestMatch(geneName, allScoresD)

    # Call closestMatch on the best match from the other species
    reciprocal_match = closestMatch(best_match, allScoresD)

    if reciprocal_match is geneName:
        # If the output matches geneName, print the best reciprocal hit
        gene1_info = geneD[geneName]
        gene2_info = geneD[best_match]

        # Concat to provide desired info
        print(f"{gene1_info[0]} {gene1_info[1]} {geneName} --- {gene2_info[0]} {gene2_info[1]} {best_match}")

#printBRH('c8',allScoresD)
# chr2 43123243 c8 --- chr3 45016733 h17
#printBRH('h7',allScoresD)
# None
#printBRH('h17',allScoresD)
# chr3 45016733 h17 --- chr2 43123243 c8

# Wrapper function on sample set
def runBRHSample():
    '''Print best reciprocal hits for sample data. First in human
    chromosome order, then in chicken chromosome order.'''
    allScoresD = allScores(sampleHumanGeneList, sampleChickenGeneList)
    print('human --- chicken')
    for geneName in sampleHumanGeneList:
        printBRH(geneName, allScoresD)
    print()
    print('chicken --- human')
    for geneName in sampleChickenGeneList:
        printBRH(geneName, allScoresD)

runBRHSample()
# human --- chicken
# chr11 118415243 h4 --- chr24 5629899 c19
# chr11 133938820 h6 --- chr24 2542440 c17
# chr14 72399156 h9 --- chr5 28862733 c22
# chr3 45016733 h17 --- chr2 43123243 c8
#
# chicken --- human
# chr2 43123243 c8 --- chr3 45016733 h17
# chr24 2542440 c17 --- chr11 133938820 h6
# chr24 5629899 c19 --- chr11 118415243 h4
# chr5 28862733 c22 --- chr14 72399156 h9

# Wrapper function on full dataset
def runBRH():
    '''Print best reciprocal hits for full data. First in human
    chromosome order, then in chicken chromosome order.'''
    allScoresD=allScores(humanGeneList,chickenGeneList)
    print ('human --- chicken')
    for geneName in humanGeneList:
        printBRH(geneName,allScoresD)
    print()
    print ('chicken --- human')
    for geneName in chickenGeneList:
        printBRH(geneName,allScoresD)

#runBRH()
# Based on the output, it does not appear that human chromosome X is homologous to chicken chromosome Z.
# This is because neither chromosome shares an ortholog with the other.

# human --- chicken
# chr1 8921061 h1 --- chr21 3197152 c15
# chr1 19967048 h2 --- chr21 4801383 c16
# chr1 165370159 h3 --- chr8 5341177 c26
# chr11 118415243 h4 --- chr24 5629899 c19
# chr11 122943035 h5 --- chr24 3068439 c18
# chr11 133938820 h6 --- chr24 2542440 c17
# chr14 59951161 h8 --- chr5 57373453 c25
# chr14 72399156 h9 --- chr5 28862733 c22
# chr14 77940740 h10 --- chr5 41664432 c23
# chr14 104552016 h11 --- chr5 53149708 c24
# chr22 19710468 h13 --- chr15 791273 c6
# chr22 35936915 h14 --- chr1 54068502 c1
# chr22 43088127 h15 --- chr1 70507590 c2
# chr3 45016733 h17 --- chr2 43123243 c8
# chr5 36248536 h18 --- chrZ 10378546 c27
# chr5 41925356 h19 --- chrZ 12697502 c28
# chr5 64813593 h20 --- chrZ 19960479 c29
# chr5 112196885 h22 --- chrZ 45168069 c31
# chr6 3118608 h23 --- chr2 67495019 c10
# chr6 12290596 h24 --- chr2 63816357 c9
# chr8 72753784 h25 --- chr2 121960191 c11
# chr8 109455830 h26 --- chr2 137539936 c12
# chr8 117654369 h27 --- chr2 140871172 c13
# chr9 37485932 h28 --- chrZ 74199054 c33
# chr9 80037995 h29 --- chrZ 37281239 c30
# chrX 13671225 h30 --- chr1 126408830 c5
# chrX 13789150 h31 --- chr1 126252299 c4
# chrX 15843929 h32 --- chr1 125276190 c3
# chrX 118533023 h34 --- chr4 16697566 c21
#
# chicken --- human
# chr1 54068502 c1 --- chr22 35936915 h14
# chr1 70507590 c2 --- chr22 43088127 h15
# chr1 125276190 c3 --- chrX 15843929 h32
# chr1 126252299 c4 --- chrX 13789150 h31
# chr1 126408830 c5 --- chrX 13671225 h30
# chr15 791273 c6 --- chr22 19710468 h13
# chr2 43123243 c8 --- chr3 45016733 h17
# chr2 63816357 c9 --- chr6 12290596 h24
# chr2 67495019 c10 --- chr6 3118608 h23
# chr2 121960191 c11 --- chr8 72753784 h25
# chr2 137539936 c12 --- chr8 109455830 h26
# chr2 140871172 c13 --- chr8 117654369 h27
# chr21 3197152 c15 --- chr1 8921061 h1
# chr21 4801383 c16 --- chr1 19967048 h2
# chr24 2542440 c17 --- chr11 133938820 h6
# chr24 3068439 c18 --- chr11 122943035 h5
# chr24 5629899 c19 --- chr11 118415243 h4
# chr4 16697566 c21 --- chrX 118533023 h34
# chr5 28862733 c22 --- chr14 72399156 h9
# chr5 41664432 c23 --- chr14 77940740 h10
# chr5 53149708 c24 --- chr14 104552016 h11
# chr5 57373453 c25 --- chr14 59951161 h8
# chr8 5341177 c26 --- chr1 165370159 h3
# chrZ 10378546 c27 --- chr5 36248536 h18
# chrZ 12697502 c28 --- chr5 41925356 h19
# chrZ 19960479 c29 --- chr5 64813593 h20
# chrZ 37281239 c30 --- chr9 80037995 h29
# chrZ 45168069 c31 --- chr5 112196885 h22
# chrZ 74199054 c33 --- chr9 37485932 h28