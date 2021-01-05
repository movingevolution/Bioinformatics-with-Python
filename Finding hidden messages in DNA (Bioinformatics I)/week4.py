import random


def RandomMotifs(Dna, k, t):
    M = len(Dna[0]) - k
    randomKmers = []
    for item in Dna:
        i = random.randint(0, M)
        randomKmer = item[i:i + k]
        randomKmers.append(randomKmer)
    return randomKmers


def ProfileWithPseudocounts(Motifs):
    pseudoCount = CountWithPseudocounts(Motifs)
    numberStrings = len(Motifs[0])
    for key in pseudoCount:
        for j in range(numberStrings):
            pseudoCount[key][j] = pseudoCount[key][j] / (numberStrings + 4)
    return pseudoCount


def CountWithPseudocounts(Motifs):
    count = {}
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
            count[symbol].append(1)

    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    return count


def Motifs(Profile, Dna):
    k = len(Profile['A'])
    kmers = []
    for string in Dna:
        findKmer = ProfileMostProbablePattern(string, k, Profile)
        kmers.append(findKmer)
    return kmers


def ProfileMostProbablePattern(Text, k, Profile):
    mostProbable = ''
    biggestP = -1
    for i in range(len(Text) - k + 1):
        p = Pr(Text[i:i + k], Profile)
        if p > biggestP:
            biggestP = p
            mostProbable = Text[i:i + k]
    return mostProbable


def Pr(Text, Profile):
    p = 1
    for i in range(len(Text)):
        p = p * Profile[Text[i]][i]
    return p


def Score(Motifs):
    k = len(Motifs[0])
    t = len(Motifs)
    consensus = Consensus(Motifs)
    score = 0
    for i in range(t):
        for j in range(k):
            if consensus[j] != Motifs[i][j]:
                score += 1
    return score


def Consensus(Motifs):
    k = len(Motifs[0])
    count = Count(Motifs)

    consensus = ""
    for j in range(k):
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol
    return consensus


def Count(Motifs):
    count = {}
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
            count[symbol].append(0)

    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    return count


def RandomMotifSearch(Dna, k, t):
    M = RandomMotifs(Dna, k, t)
    BestMotifs = M
    while True:
        Profile = ProfileWithPseudocounts(M)
        M = Motifs(Profile, Dna)
        if Score(M) < Score(BestMotifs):
            BestMotifs = M
        else:
            return BestMotifs



'''with open('stepic.txt') as input_data:
    k, t = map(int, input_data.readline().split())
    dna = [line.strip() for line in input_data]
    score = math.inf
    min_motifs = []
    for i in range(1000):
        generated_motifs = RandomMotifSearch(dna, k, t)
        generated_score = Score(generated_motifs)
        if generated_score < score:
            score = generated_score
            min_motifs = generated_motifs

    for j in min_motifs:
        print(j)'''



