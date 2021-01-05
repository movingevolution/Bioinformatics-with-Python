import numpy as np
from collections import defaultdict
import itertools

def ReverseComplement(Pattern):
    comp_list = []
    for i in range(len(Pattern) - 1, -1, -1):
        if Pattern[i] == 'A':
            comp_list.append('T')
        elif Pattern[i] == 'T':
            comp_list.append('A')
        elif Pattern[i] == 'C':
            comp_list.append('G')
        elif Pattern[i] == 'G':
            comp_list.append('C')
    comp_str = ''.join(comp_list)
    return comp_str

def PatterntoNumber(kmer):
    total = 0
    nucleotide_to_number = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    for i in range(len(kmer)):
        total = (total * 4) + nucleotide_to_number[kmer[i]]
    return total

def NumberToPattern(numb,l):

    seq = []
    number_to_nucleotide = {0:'A', 1:'C', 2:'G', 3:'T'}
    for i in range(l):
        r = int(numb % 4)
        numb = int(numb / 4)
        seq.append(number_to_nucleotide[r])
    while len(seq) < l:
        seq.append(number_to_nucleotide[0])
    seq.reverse()
    string = ''.join(seq)
    return string


def MinimumSkew(Genome):

    nucleotide_to_number = {'A': 0, 'C': -1, 'G': 1, 'T': 0, '\n':0}
    holder = defaultdict(list)
    holder[0] = [0]
    skew = 0
    for i in range(len(Genome)):
        skew += nucleotide_to_number[Genome[i]]
        holder[skew].append(i+1)

    return holder[min(holder)]

def HammingDistance(p,q):
    count = 0
    for i in range(len(p)):
        if p[i] != q[i]:
            count += 1

    return count

def ApproximatePatternMatching(Text, Pattern, d):
    positions = []
    for i in range(len(Text)-len(Pattern) + 1):
        if HammingDistance(Text[i:i+len(Pattern)],Pattern) <= d:
            positions.append(str(i))
    pos_str = ' '.join(positions)
    return pos_str

def AproximatePatternCount(Text,Pattern,d):
    count = 0
    for i in range(len(Text) - len(Pattern) + 1):
        P_test = Text[i:i+len(Pattern)]
        if HammingDistance(Pattern,P_test) <= d:
            count += 1

    return count


def ImmediateNeighbors(pattern):
    neighbor = [pattern]
    for i in range(len(pattern)):
        symbol = pattern[i]
        for x in ['A','C','G','T']:
            if x != symbol:
                neighbor.append(pattern[:i]+x+pattern[i+1:])
    n_str = ' '.join(neighbor)
    return n_str


def neighbors(pattern,d):
    if d == 0:
        return pattern
    if len(pattern) == 1:
        return ['A','C','G','T']

    neighborhood = []
    suffixNeighbors = neighbors(pattern[1:],d)
    for text in suffixNeighbors:
        if HammingDistance(pattern[1:],text) < d:
            for x in ['A','C','G','T']:
                neighborhood.append(x+text)
        else:
            neighborhood.append(pattern[:1]+text)
    neighborhood = list(set(neighborhood))
    return neighborhood


def IterativeNeighbors(Pattern, d):
    Neighborhood = [Pattern]
    holder = 0
    for i in range(d):
        for k in range(len(Neighborhood)):
            holder = Neighborhood[k]
            for j in range(len(holder)):
                symbol = holder[j]
                for x in ['A', 'C', 'G', 'T']:
                    if x != symbol:
                        Neighborhood.append(holder[:j] + x + holder[j + 1:])
    Neighborhood = list(set(Neighborhood))

    return Neighborhood

def FrequentWordsWithMismatches(Text,k,d):
    FrequentPatterns =[]
    NB_hood = []
    for i in range(len(Text)-k+1):
        NB_hood.append(IterativeNeighbors(Text[i:i+k], d))

    NeighborhoodArray = list(itertools.chain(*NB_hood)) #converting listoflist to list

    Index = [0] * (len(NeighborhoodArray))
    Count = [0] * (len(NeighborhoodArray))
    for i in range(len(NeighborhoodArray)):
        Pattern = NeighborhoodArray[i]
        Index[i] = PatterntoNumber(Pattern)
        Count[i] = 1

    SortedIndex = sorted(Index)
    for i in range(len(NeighborhoodArray) - 1):
        if SortedIndex[i] == SortedIndex[i+1]:
            Count[i+1] = Count[i] + 1

    MaxCount = max(Count)

    for i in range(len(NeighborhoodArray)):
        if Count[i] == MaxCount:
            n_string = NumberToPattern(SortedIndex[i],k)
            FrequentPatterns.append(n_string)

    return FrequentPatterns






'''
f= open("submit.txt","w+")
f.write(r)
f.close()'''



