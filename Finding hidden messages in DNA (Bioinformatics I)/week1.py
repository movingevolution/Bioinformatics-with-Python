import numpy as np
from collections import defaultdict
np.set_printoptions(threshold=np.inf)


def PatternCount(Text, Pattern):
    count = 0
    for i in range(len(Text) + 1):
        if Text[i - len(Pattern):i] == Pattern:
            count = count + 1
    return count


def FrequentWords(Text, k):
    first_key = Text[0:k]
    k_dict = {first_key: 1}
    for i in range(1, len(Text) - k + 1):
        if Text[i:k + i] in k_dict:
            k_dict[Text[i:k + i]] += 1
        else:
            k_dict[Text[i:k + i]] = 1
    max_list = [key for key, value in k_dict.items() if value == max(k_dict.values())]
    return sorted(max_list)


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


def PatternMatching(Pattern, Genome):
    pos_list = []
    for i in range(len(Genome) + 1):
        if Genome[i - len(Pattern):i] == Pattern:
            pos_list.append(str(i - len(Pattern)))
    pos_str = ' '.join(pos_list)
    return pos_str



def FindingClump(genome, k, L, t):
    holder = defaultdict(list)
    result = set()

    for i in range(len(genome) - k + 1):
        seg = genome[i:i + k]

        while holder[seg] and i + k - holder[seg][0] > L:
            holder[seg].pop(0)

        holder[seg].append(i)
        if len(holder[seg]) == t:
            result.add(seg)

    return result

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

def ComputingFrequencies(Text,k):

    FrequencyArray = np.zeros([4**k-1], dtype=int)

    for i in range(len(Text)-k+1):
        Pattern = Text[i:i+k]
        j = PatterntoNumber(Pattern) #Using prior function
        FrequencyArray[j] = FrequencyArray[j]+1

    return FrequencyArray




'''li = ComputingFrequencies('',8)
st = str(li)
f= open("submit.txt","w+")
f.write(st)
f.close()'''

