import numpy as np
import math
import random


def NumberToPattern(numb, l):
    seq = []
    number_to_nucleotide = {0: 'A', 1: 'C', 2: 'G', 3: 'T'}
    for i in range(l):
        r = int(numb % 4)
        numb = int(numb / 4)
        seq.append(number_to_nucleotide[r])
    while len(seq) < l:
        seq.append(number_to_nucleotide[0])
    seq.reverse()
    string = ''.join(seq)
    return string


def HammingDistance(p, q):
    count = 0
    for i in range(len(p)):
        if p[i] != q[i]:
            count += 1

    return count


def neighbors(pattern, d):
    if d == 0:
        return pattern
    if len(pattern) == 1:
        return ['A', 'C', 'G', 'T']

    neighborhood = []
    suffixNeighbors = neighbors(pattern[1:], d)
    for text in suffixNeighbors:
        if HammingDistance(pattern[1:], text) < d:
            for x in ['A', 'C', 'G', 'T']:
                neighborhood.append(x + text)
        else:
            neighborhood.append(pattern[:1] + text)
    neighborhood = list(set(neighborhood))
    return neighborhood


def DistanceBetweenPatternAndStrings(Pattern, Dna):
    Dna_list = Dna.split()
    k = len(Pattern)
    Distance = 0
    for Text in Dna_list:
        Hamming_Distance = float('inf')  # HammingDistance initialized to infinity
        for i in range(len(Text) - k + 1):
            if Hamming_Distance > HammingDistance(Pattern, Text[i:i + k]):
                Hamming_Distance = HammingDistance(Pattern, Text[i:i + k])
        Distance = Distance + Hamming_Distance

    return Distance


def MedianString(dna, k):
    distance = float('inf')
    for i in range((4 ** k) - 1):
        pattern = NumberToPattern(i, k)
        current_distance = DistanceBetweenPatternAndStrings(pattern, dna)
        if distance > current_distance:
            distance = current_distance
            median = pattern

    return median



def ProfileMostProbableKmer(text, k, profile):
    entropy_dict = {}

    for i in range(len(text) - k + 1):
        pattern = text[i:i + k]
        entropy = 1
        for j in range(len(pattern)):
            base = profile[pattern[j]][j]
            entropy = entropy * base
        if pattern not in entropy_dict:
            entropy_dict[pattern] = entropy

    max_list = [key for key, value in entropy_dict.items() if value == max(entropy_dict.values())]
    string = ''.join(max_list)
    return string[0:k]


def Count(Motifs):
    k = len(Motifs[0])
    t = len(Motifs)

    count = {new_list: np.zeros(k) for new_list in 'ACGT'}

    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1

    return count


def Profile(Motifs):
    k = len(Motifs[0])
    t = len(Motifs)

    profile = {new_list: np.zeros(k) for new_list in 'ACGT'}

    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            profile[symbol][j] += 1 / t

    return profile


def Consensus(Motifs):
    t = len(Motifs)
    k = len(Motifs[0])
    count = Count(Motifs)
    consensus = ""

    for j in range(k):
        m = 0
        most_frequent = ""
        for base in "ACGT":
            if count[base][j] > m:
                m = count[base][j]
                most_frequent = base
        consensus += most_frequent
    return consensus


def Score(Motifs):
    consensus = Consensus(Motifs)
    score = 0
    k = len(consensus)
    t = len(Motifs)
    for i in range(t):
        s = HammingDistance(consensus, Motifs[i])
        score += s

    return score


def GreedyMotifSearch(Dna, k, t):
    # type your GreedyMotifSearch code here.
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])
    n = len(Dna[0])
    for i in range(n - k + 1):
        Motifs = []
        Motifs.append(Dna[0][i:i + k])
        for j in range(1, t):
            P = Profile(Motifs[0:j])
            Motifs.append(ProfileMostProbableKmer(Dna[j], k, P))
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs

def ProfilePseudo(Motifs):
    k = len(Motifs[0])
    t = len(Motifs)

    profile = {new_list: np.ones(k) for new_list in 'ACGT'}

    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            profile[symbol][j] += 1 / (t+4)

    return profile


def GreedyMotifSearchPseudo(Dna, k, t):
    # type your GreedyMotifSearch code here.
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])
    n = len(Dna[0])
    for i in range(n - k + 1):
        Motifs = []
        Motifs.append(Dna[0][i:i + k])
        for j in range(1, t):
            P = ProfilePseudo(Motifs[0:j])
            Motifs.append(ProfileMostProbableKmer(Dna[j], k, P))
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs

def Motifs(profile,dna):
    k = len(profile['A'])
    k_mers = []

    for string in dna:
        new_kmer = ProfileMostProbableKmer(string,k,profile)
        k_mers.append(new_kmer)

    return k_mers

def RandomKmer(dna,k,t):
    M = len(dna[0]) - k
    random_kmers = []
    for string in dna:
        i = random.randint(0, M)
        random_kmer = string[i:i + k]
        random_kmers.append(random_kmer)
    return random_kmers

def RandomMotifSearch(dna,k,t):

    random_kmers = RandomKmer(dna,k,t)
    BestMotifs = random_kmers

    while True:
        P = ProfilePseudo(random_kmers)
        new_motifs = Motifs(P,dna)

        if Score(new_motifs) < Score(BestMotifs):
            BestMotifs = new_motifs
        else:
            return BestMotifs






