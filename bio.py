from Bio import SeqIO
from Bio.Seq import Seq

def sequenceFinder(found_start_codon_positions,found_stop_codon_positions,seq_record):
    found_sequence = []
    k = 0
    l = 0
    maxi = 0
    for k in range(len(found_start_codon_positions)):
        if(len(found_start_codon_positions) == 1 or len(found_stop_codon_positions) == 1):
            break
        for l in range(len(found_stop_codon_positions)):
            if(len(found_start_codon_positions) == 1 or len(found_stop_codon_positions) == 1):
                break
            if(k+1 > len(found_start_codon_positions or l > len(found_stop_codon_positions))):
                break
            if (found_start_codon_positions[k] < (found_stop_codon_positions[l] - 3) and (found_stop_codon_positions[l] + 3 - found_start_codon_positions[k])%3 == 0) and ((maxi - found_start_codon_positions[k])%3 != 0 or maxi < found_start_codon_positions[k]):
                maxi = found_stop_codon_positions[l]
                found_sequence.append(seq_record.seq[found_start_codon_positions[k]:found_stop_codon_positions[l]+3])
                found_stop_codon_positions[l] = 0
                found_start_codon_positions.pop(k)

    return found_sequence

def sequence2Finder(found_start_codon_positions,found_stop_codon_positions,seq_record):
    found_sequence = []
    k = 0
    l = 0
    maxl = l
    while l < len(found_stop_codon_positions) and k < len(found_start_codon_positions):
        if(found_start_codon_positions[k] + 3 - found_stop_codon_positions[l])%3 == 0:
            maxCodon = found_stop_codon_positions[l]
            maxl = l
        if (found_stop_codon_positions[l] > found_start_codon_positions[k] ) and maxl != 0:
            found_sequence.append(seq_record.seq[found_stop_codon_positions[maxl]:found_start_codon_positions[k]+3])
            found_stop_codon_positions.pop(maxl)
            l = 0
            k += 1
        else:
            l += 1

    return found_sequence

def toFindStartCodon(seq_record):
    start_codons = ['ATG']
    found_start_codon_positions = []
    n = len(seq_record)
    k = 0
    while k < n-2:
        possible_codon = str(seq_record.seq[k:k+3])
        if possible_codon in start_codons:
            found_start_codon_positions.append(k)
        
        k += 1
    return found_start_codon_positions

def toFindStopCodon(seq_record):
    stop_codons = ['TAG', 'TAA', 'TGA']
    found_stop_codon_positions = []
    n = len(seq_record)
    k = 0
    while k < n-2:
        possible_codon = str(seq_record.seq[k:k+3])
        if possible_codon in stop_codons:
            found_stop_codon_positions.append(k)
        k += 1

    return found_stop_codon_positions

def filterByHundread(seq_record):
    i=0
    while i < len(seq_record):
        if(len(seq_record[i]) < 99):
            seq_record.pop(i)
            i += -1
        i += 1
    return seq_record

def toFindCodon(seq_record):
    codons = []
    letters = ["A", "T", "C", "G"]
    newCodons = []
    for i in range(len(letters)):
        for j in range(len(letters)):
            for z in range(len(letters)):
                newCodons.append(letters[i]+letters[j]+ letters[z])
            
    frequncy = {}
    for str in newCodons:
        frequncy[str] = 0

    n = len(seq_record)
    k = 0
    for i in range(len(seq_record)):
        for j in range(len(seq_record[i])-2):
            k+=1
            frequncy[(seq_record[i][j:j+3])] += 1
    for str in newCodons:
        frequncy[str] = round( 100 * frequncy[str] / k, 3)
    return frequncy

def toFindDicodon(seq_record):
    letters = ["A", "T", "C", "G"]
    newDicodons = []
            
    for i in range(len(letters)):
        for j in range(len(letters)):
            for z in range(len(letters)):
                for a in range(len(letters)):
                    for b in range(len(letters)):
                        for c in range(len(letters)):
                            newDicodons.append(letters[i] + letters[j] + letters[z] + letters[a] + letters[b] + letters[c])
            
    frequncy = {}
    for str in newDicodons:
        frequncy[str] = 0

    k = 0
    for i in range(len(seq_record)):
        for j in range(len(seq_record[i])-5):
            k+=1
            frequncy[(seq_record[i][j:j+6])] += 1
    for str in newDicodons:
        frequncy[str] = round(100 * frequncy[str] / k, 7)
    return frequncy

def frequenciesComparasion(frequency):
    print(len(frequency))
    for i in frequency:
        sum = 0
        print(i, end = ' ')
        for z in frequency:
            for j in frequency[i]:
                sum += abs(frequency[i][j] - frequency[z][j])
                
            print(round(sum,4), end = ' ')
            sum = 0
        print("\n")

    return frequency

data = ["bacterial1.fasta", "bacterial2.fasta", "bacterial3.fasta", "bacterial4.fasta", "mamalian1.fasta", "mamalian2.fasta", "mamalian3.fasta", "mamalian4.fasta"]

start_codons = ['ATG']
stop_codons = ['TAG', 'TAA', 'TGA']

frequencyCod = {}
frequencyDicod = {}

for i in data:
    seq_record = SeqIO.read("data/" + i, "fasta")
    reversed_seq_record = seq_record.reverse_complement()
    ###################################################1############################################
    findSequence = sequenceFinder(toFindStartCodon(seq_record),toFindStopCodon(seq_record),seq_record)
    findReversedSequence = sequenceFinder(toFindStartCodon(reversed_seq_record),toFindStopCodon(reversed_seq_record),reversed_seq_record)

    ###################################################2############################################
    stopToStart = sequence2Finder(toFindStartCodon(seq_record),toFindStopCodon(seq_record),seq_record)
    allSequence = findSequence + findReversedSequence + stopToStart
    ###################################################3############################################
    allSequence = filterByHundread(allSequence)
    ###################################################4############################################
    frequencyCod[seq_record.id] = toFindCodon(allSequence)
    frequencyDicod[seq_record.id] = toFindDicodon(allSequence)

frqCodon = frequenciesComparasion(frequencyCod)
frqDicodon = frequenciesComparasion(frequencyDicod)                        