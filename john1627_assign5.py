'''
Assignment 5
CSE 256
Jenna Johnson
'''
# import relevant packages
import sys
import re
# define function that reads a fasta file
def readFASTA(file):
    names = []
    desc = []
    sequences = []

    with open(file, 'r') as fasta:
        SeqName = ""
        SeqDesc = ""
        Seq = ""

        for line in fasta:
            line = line.strip()
            if line.startswith(">"):
                if SeqName:
                    names.append(SeqName)
                    desc.append(SeqDesc)
                    sequences.append(Seq)
                defline = line[1:].split(" ", 1)
                SeqName = defline[0].replace("<", "")
                SeqDesc = defline[1] if len(defline) > 1 else ""
                Seq = ""
            else:
                Seq = Seq + line

        if SeqName:
            names.append(SeqName)
            desc.append(SeqDesc)
            sequences.append(Seq)

    return names, desc, sequences
# define function that translates a given sequence into protein
def translation(seq):
    # defining the genetic code dict for translation
    genetic_code = {'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
        'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
        'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
        'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
        'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
        'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
        'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
        'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
        'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
        'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
        'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
        'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
        'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
        'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
        'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
    }
    proteinSeq = ''
    codon = ''
    # iterate through the nucleotides in seq and assign translated var
    for nucleotides in seq:
        codon += nucleotides
        if len(codon) == 3:
            aminoAcid = genetic_code.get(codon, 'x')  # 'x' for unknown codon
            proteinSeq += aminoAcid
            codon = ''
    # return the translated protein out of the function
    return proteinSeq
# define function which calculates gc content in seq
def gcContent(seq):
    # calculate gc percentage
    gcPercent = (seq.count("G") + seq.count("C")) / len(seq) * 100
    # return the percentage out of the function
    return gcPercent
# define function which get num purines, pyrimidines
def get_Purine_Pyrimidine(seq):
    purines = 0
    pyrimidines = 0
    for nucleotides in seq:
        # assign purines and pyrimidines based on content
        if nucleotides in 'AG':
            purines += 1
        elif nucleotides in 'CT':
            pyrimidines += 1
    # return the purines, pyrimidines out of the function
    return purines, pyrimidines
# define function which finds maximum value in a list
def getMax(SeqNames, ListInt):
    value = ListInt[0]
    maxSeqName = SeqNames[0]
    for i in range(1, len(ListInt)):
        # iterate through the list
        if ListInt[i] > value:
            value = ListInt[i]
            maxSeqName = SeqNames[i]
    return (maxSeqName, value)

# define function which gets the minimum value in a list
def getMin(SeqNames, ListInt):
    val = ListInt[0]
    minSeqName = SeqNames[0]
    # iterate through the list
    for i in range(1, len(ListInt)):
        if ListInt[i] < val:
            val = ListInt[i]
            minSeqName = SeqNames[i]
    return (minSeqName, val)
# define a function which gets the average value in a list of integers
def getAverage(ListInt):
    #calculating average
    Total = sum(ListInt)
    average = Total / len(ListInt)
    return average
# a function which returns the counts of all nucleotides in each sequence
def nucleotide_counter_new(seq):
    d = {}  # declare empty dictionary
    # initialize all counts to 0
    seqA = 0
    seqT = 0
    seqG = 0
    seqC = 0
    for i in seq:  # loop to iterate over all elements, counting nucleotides
        if i == 'A':
            seqA = seqA + 1
        if i == 'T':
            seqT = seqT + 1
        if i == 'G':
            seqG = seqG + 1
        if i == 'C':
            seqC = seqC + 1
    # append dict to include counts
    d['A'] = seqA
    d['T'] = seqT
    d['G'] = seqG
    d['C'] = seqC
    return d  # return the dict out of the function
# function returns codons and occurrences of each
def get_codon_new(Seq):
    codon_list = []  # initialize empty list for codons
    for pos in range(0, len(Seq) - 2, 1):
        codon_list.append(Seq[pos: pos + 3])  # append to the list
    number_of_coden = len(codon_list)
    D = dict((codon, codon_list.count(codon)) for codon in set(codon_list))
    return D, number_of_coden
# function which finds common seq in files by set methods
def findCommon_1(sequences1, sequences2):
    seq1Set = set(sequences1)
    seq2Set = set(sequences2)
    if (seq1Set & seq2Set):
        common = list(seq1Set & seq2Set)
    else:
        print("No common elements")
    return common
# function which finds common seqs by for loop
def findCommon_2(sequences1, sequences2):
    commonSeq = []
    for i in sequences1:
        for z in sequences2:
            if i == z:
                commonSeq.append(i)
    return commonSeq
# function which finds unique seq in files by set methods
def findUnique_1(sequences1, sequences2):
    setSeq2 = set(sequences2)
    uniqueInSeq1 = [x for x in sequences1 if x not in setSeq2]
    seqSeq1 = set(sequences1)
    uniqueInSeq2 = [z for z in sequences2 if z not in seqSeq1]
    return uniqueInSeq1, uniqueInSeq2
# find unique names using for loops
def findUnique_2(sequences1, sequences2):
    uniqueS1 = [seq for seq in sequences1 if seq not in sequences2]
    uniqueS2 = [seq for seq in sequences2 if seq not in sequences1]
    return uniqueS1, uniqueS2
# returns the occurrences of a given restriction enzyme site and locations
def enzymeFinder(sequence, restrictionEnzyme):
    pattern = re.compile(restrictionEnzyme)
    matches = pattern.finditer(sequence)
    results = {}
    count = 1
    for match in matches:
        start = match.start()
        end = match.end()
        results[str(count)] = f"{start + 1}-{end}"
        count +=1
    return results

# Initialize lists to store data from each file separately
namesSeq1, descSeq1, sequencesSeq1 = readFASTA(sys.argv[1])
namesSeq2, descSeq2, sequencesSeq2 = readFASTA(sys.argv[2])

print("Welcome! I am Jenna Johnson.\n\nThe 1st File Name is ({0})\n\nThe 2nd File Name is ({1})"
      .format(sys.argv[1], sys.argv[2]))

print("\n[Part 1] Summary Information\n")
print("(1) The 1st File:")
TotalSeq1 = len(sequencesSeq1)
print("(A) There are a total of ({0}) sequences in the file".format(TotalSeq1))
# B
seqLengths1 = []
for i, seq in enumerate(sequencesSeq1):
    seqLen = len(seq)
    seqLengths1.append(seqLen)
avgSeqLen1 = getAverage(seqLengths1)
print("(B) The average length of all sequences in the file is ({0})". format(avgSeqLen1))

# C
seqName1, sequence1 = getMax(namesSeq1, sequencesSeq1)
print("(C) The longest sequence is ({0}), length = [{1}]".format(seqName1, len(sequence1)))

#D
seq_name1, seq1 = getMin(namesSeq1, sequencesSeq1)
print("(D) The shortest sequence is ({0}), length = [{1}]".format(seq_name1, len(seq1)))

#E, F, G, H
ACounts1 = []
TCounts1 = []
GCounts1 = []
CCounts1 = []
for seq in sequencesSeq1:
    nucleotideCounts = nucleotide_counter_new(seq)
    ACounts1.append(nucleotideCounts['A'])
    TCounts1.append(nucleotideCounts['T'])
    GCounts1.append(nucleotideCounts['G'])
    CCounts1.append(nucleotideCounts['C'])

Aseq1, Acount1 = getMax(namesSeq1, ACounts1)
print("(E) The sequence contains the most Adenine is ({0}), count = [{1}]".format(Aseq1, Acount1))

Tseq1, Tcount1 = getMax(namesSeq1, TCounts1)
print("(F) The sequence contains the most Thymine is ({0}), count = [{1}]".format(Tseq1, Tcount1))

Gseq1, Gcount1 = getMax(namesSeq1, GCounts1)
print("(G) The sequence contains the most Guanine is ({0}), count = [{1}]".format(Gseq1, Gcount1))

Cseq1, Ccount1 = getMax(namesSeq1, CCounts1)
print("(H) The sequence contains the most Cytosine is ({0}), count = [{1}]".format(Cseq1, Ccount1))

# I, J
PurineCountsS1 = []
PurineSeqsS1 = []
PyrimCountsS1 = []
PyrimSeqsS1 = []
for i, sequence in enumerate(sequencesSeq1):
    purines, pyrimidines = get_Purine_Pyrimidine(sequence)
    PurineCountsS1.append(purines)
    PyrimCountsS1.append(pyrimidines)
    PurineSeqsS1.append((namesSeq1[i], purines))
    PyrimSeqsS1.append((namesSeq1[i], pyrimidines))
PurineSeq, PurineCnt = getMax(namesSeq1, PurineCountsS1)
PySeq, PyCnt = getMax(namesSeq1, PyrimCountsS1)

print("(I) The sequence contains most Purines is ({0}) count = [{1}]".format(PurineSeq, PurineCnt))

print("(J) The sequence contains the most Pyrimidines is ({0}), count = [{1}]".format(PySeq, PyCnt))

# K
gcContentsS1 = []
gcSeqs1 = []
for i, sequence in enumerate(sequencesSeq1):
    gcCont = gcContent(sequence)
    gcContentsS1.append(gcCont)
    gcSeqs1.append((namesSeq1[i], gcCont))
MaxgcSeq1, MaxgcPercentS1 = getMax(namesSeq1, gcContentsS1)
print("(K) The sequence contains highest GC is ({0}) GC content = [{1}%]".format(MaxgcSeq1, MaxgcPercentS1))

# L
MingcSeq1, MingcPercentS1 = getMin(namesSeq1, gcContentsS1)
print("(L) The sequence contains lowest GC is ({0}) GC content = [{1}%]".format(MingcSeq1, MingcPercentS1))

# M
avgGC1 = getAverage(gcContentsS1)
print("(M) The average GC content for all sequences is ({0}%)".format(avgGC1))

# The second file summary info
print("\n(2) The 2nd File:")
TotalSeq2 = len(sequencesSeq2)
print("(A) There are a total of ({0}) sequences in the file".format(TotalSeq2))
# B
seqLengths2 = []
for i, seq in enumerate(sequencesSeq2):
    seqLen = len(seq)
    seqLengths2.append(seqLen)
avgSeqLen2 = getAverage(seqLengths2)
print("(B) The average length of all sequences in the file is ({0})". format(avgSeqLen2))

# C
seqName2, sequence2 = getMax(namesSeq2, sequencesSeq2)
print("(C) The longest sequence is ({0}), length = [{1}]".format(seqName2, len(sequence2)))

#D
seq_name2, seq2 = getMin(namesSeq2, sequencesSeq2)
print("(D) The shortest sequence is ({0}), length = [{1}]".format(seq_name2, len(seq2)))

#E, F, G, H
ACounts2 = []
TCounts2 = []
GCounts2 = []
CCounts2 = []
for seq in sequencesSeq2:
    nucleotideCounts = nucleotide_counter_new(seq)
    ACounts2.append(nucleotideCounts['A'])
    TCounts2.append(nucleotideCounts['T'])
    GCounts2.append(nucleotideCounts['G'])
    CCounts2.append(nucleotideCounts['C'])

Aseq2, Acount2 = getMax(namesSeq2, ACounts2)
print("(E) The sequence contains the most Adenine is ({0}), count = [{1}]".format(Aseq2, Acount2))

Tseq2, Tcount2 = getMax(namesSeq2, TCounts2)
print("(F) The sequence contains the most Thymine is ({0}), count = [{1}]".format(Tseq2, Tcount2))

Gseq2, Gcount2 = getMax(namesSeq2, GCounts2)
print("(G) The sequence contains the most Guanine is ({0}), count = [{1}]".format(Gseq2, Gcount2))

Cseq2, Ccount2 = getMax(namesSeq2, CCounts2)
print("(H) The sequence contains the most Cytosine is ({0}), count = [{1}]".format(Cseq2, Ccount2))

# I, J
PurineCountsS2 = []
PurineSeqsS2 = []
PyrimCountsS2 = []
PyrimSeqsS2 = []
for i, sequence in enumerate(sequencesSeq1):
    purines, pyrimidines = get_Purine_Pyrimidine(sequence)
    PurineCountsS2.append(purines)
    PyrimCountsS2.append(pyrimidines)
    PurineSeqsS2.append((namesSeq2[i], purines))
    PyrimSeqsS2.append((namesSeq2[i], pyrimidines))
PurineSeq2, PurineCnt2 = getMax(namesSeq2, PurineCountsS2)
PySeq2, PyCnt2 = getMax(namesSeq2, PyrimCountsS2)

print("(I) The sequence contains most Purines is ({0}) count = [{1}]".format(PurineSeq2, PurineCnt2))

print("(J) The sequence contains the most Pyrimidines is ({0}), count = [{1}]".format(PySeq2, PyCnt2))

# K
gcContentsS2 = []
gcSeqs2 = []
for i, sequence in enumerate(sequencesSeq2):
    gcCont2 = gcContent(sequence)
    gcContentsS2.append(gcCont2)
    gcSeqs2.append((namesSeq2[i], gcCont2))
MaxgcSeq2, MaxgcPercentS2 = getMax(namesSeq2, gcContentsS2)
print("(K) The sequence contains highest GC is ({0}) GC content = [{1}%]".format(MaxgcSeq2, MaxgcPercentS2))

# L
MingcSeq2, MingcPercentS2 = getMin(namesSeq2, gcContentsS2)
print("(L) The sequence contains lowest GC is ({0}) GC content = [{1}%]".format(MingcSeq2, MingcPercentS2))

# M
avgGC2 = getAverage(gcContentsS2)
print("(M) The average GC content for all sequences is ({0}%)".format(avgGC2))

# remove the brackets from the printout for common and unique

common = findCommon_1(namesSeq1, namesSeq2)
common_toStr = ', '.join(str(item) for item in common)
# common2 = findCommon_2(namesSeq1, namesSeq2) - This uses findCommon_2, gets same result
print("\n(3) The sequences common in both files:\n", common_toStr)

# uniqueInSeq1, uniqueInSeq2 = findUnique_1(namesSeq1, namesSeq2) - uses findUnique_1 - same result as below
uniq1, uniq2 = findUnique_2(namesSeq1, namesSeq2)
uniq1_toStr = ', '.join(str(item) for item in uniq1)
uniq2_toStr = ', '.join(str(item) for item in uniq2)
print("\n(4) The sequences unique in the 1st file:\n", uniq1_toStr)
print("\n(5) The sequences unique in the 2nd file: ", uniq2_toStr)

print("\n[Part 2] Examine individual sequence")
while True:
    seq = input("\nWhich sequence do you want to see ('q' for exit) ? ")
    # check for seq in names list first file, then assign to proper sequence contents
    if seq in namesSeq1:
        desiredSeqIndex = namesSeq1.index(seq)
        desiredSeq = sequencesSeq1[desiredSeqIndex]
        translated = translation(desiredSeq)
        gcCont = gcContent(desiredSeq)
        Dict, num_coden = get_codon_new(desiredSeq)
        nucleotideCounts = nucleotide_counter_new(desiredSeq)
        # print sequence data
        print("\n1. SeqName [{0}]\n\tDescription [{1}]\n\tSequence [{2}]\n\tTranslated Protein [{3}]\n\tNucleotide "
          "Counts: A = [{4}] T = [{5}] G = [{6}] C = [{7}]\n\tGC Content: ({8}%)\n\tCodon Profile:\n\t\tCodon\t\tCount"
          .format(seq, descSeq1[desiredSeqIndex], sequencesSeq1[desiredSeqIndex], translated, nucleotideCounts['A'],
                  nucleotideCounts['T'], nucleotideCounts['G'], nucleotideCounts['C'], gcCont, ))
        for key, values in Dict.items():
            print("\t\t", key, "\t\t", values)
    # check if seq in names list for second file, assign to proper sequence contents
    elif seq in namesSeq2:
        desiredSeqIndex = namesSeq2.index(seq)
        desiredSeq = sequencesSeq2[desiredSeqIndex]
        translated = translation(desiredSeq)
        gcCont = gcContent(desiredSeq)
        Dict, num_coden = get_codon_new(desiredSeq)
        nucleotideCounts = nucleotide_counter_new(desiredSeq)

        print("\n1. SeqName [{0}]\n\tDescription [{1}]\n\tSequence [{2}]\n\tTranslated Protein [{3}]\n\tNucleotide "
          "Counts: A = [{4}] T = [{5}] G = [{6}] C = [{7}]\n\tGC Content: ({8}%)\n\tCodon Profile:\n\t\tCodon\t\tCount"
          .format(seq, descSeq2[desiredSeqIndex], sequencesSeq2[desiredSeqIndex], translated, nucleotideCounts['A'],
                  nucleotideCounts['T'], nucleotideCounts['G'], nucleotideCounts['C'], gcCont, ))
        for key, values in Dict.items():
            print("\t\t", key, "\t\t", values)
    # quit if enter q
    elif seq == "q":
        exit()
    else:
        print("Sequence not found")
        continue

# while loop to keep asking to enter enzyme until enter q to quit, will loop back to statement above for false
    # ask for another sequence to check
    while True:
        restrictionEnzyme = input("\nQuestion 2: Which Restriction Enzyme to search? (e.q. GAAATC) ('q' for exit) ")
        if restrictionEnzyme == 'q':
            break
        detectedResults = enzymeFinder(desiredSeq, restrictionEnzyme)
        # if the site is not found
        if len(detectedResults) == 0:
            print("Not found")
        else:
            print("\nThere are {0} {1} restriction enzyme sites detected:".format(len(detectedResults), restrictionEnzyme))
            for site, position in detectedResults.items():
                start, end = position.split('-')
                print("({0}) {1} Start Position [{2}] End Position [{3}]".format(site, restrictionEnzyme, start, end))
