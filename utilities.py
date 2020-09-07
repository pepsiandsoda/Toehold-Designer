import os
import os.path
from subprocess import Popen, PIPE
import time


def split_sequence(sequence, window):
    '''
    Sliding window that splits the sequence into smaller sequences
    
    Args: 
        sequence(str): target RNA strand
        window(int): length of sliding window
        
    Returns(lst of strs): list of smaller sequence windows derived 
    from the original sequence
    '''
    
    sequences = []
    last_window_start_pos = len(sequence) - window + 1   #+1 as final index value = value - 1 
    
    for i in range(0, last_window_start_pos):
        sequences.append(sequence[i:window + i])

    return sequences


def reversed_complement(sequence):
    '''
    Creates the reverse complement sequence of an RNA molecule
    
    Sequence(str): input RNA strand
    
    Returns (str): reverse complementary strand of input RNA
    '''
    
    nt_pairing = {'A': 'U', 'G': 'C', 'U': 'A', 'C': 'G'}
    sequence_upper = sequence.upper()  #sequence has same format as strings in nt_pairing

    complement = ''
    for nt in sequence_upper:
        complement += nt_pairing[nt]  #add complementary nt from nt_pairing to complement sequence 

    return complement[::-1]  # reverse the sequence


def no_stop(sequence):
     '''
    Check for stop codons
    
    sequence(str): input RNA strand
    
    returns(Bool): whether the input strand contains a stop codon
    '''
        
    stop = ['UAA', 'UAG', 'UGA']

    for i in range(0, len(sequence), 3):  #only checks one reading frame - fix
        if sequence[i:i + 3] in stop:
            return False

    return True


def possible_toehold_B(reg, rev):
    stem_5 = 'UGCAUCCUCCUCCUCCU'
    stem_3 = 'AGGAGGAGAAAAAUGCA'
    linker = 'AACCUGGCGGCAGCGCAAAAG'
    
    toeholds = []

#for n in ['A', 'G', 'U', 'C']:  adding the nt if there is not a stop codon, otherwise, if there is a stop codon then it adds in a different nt there. This is probably to make the number of nt between the AUG and reporter a multiple of three - makes sure things are inline 
    
    if no_stop(stem_3[-2:] + linker):
        toehold.append(stem_5 + loop + stem_3 + linker)
        
        #toeholds.append(rev + loop + reg[0:11] + n + linker)  #if want to add in extra nt then include this line

    return toeholds


def single_streadness(sequence, result_path, wait=1):
    file = open('{}pipo.in'.format(result_path), 'w')
    file.write("{}\n".format(sequence))
    file.close()

    Popen(["pairs", "{}pipo".format(result_path)], stdout=PIPE)
    time.sleep(wait)
    with open("{}pipo.ppairs".format(result_path)) as res:
        parsed_res = parse_pairs_result(res, len(sequence))

    os.remove("{}pipo.ppairs".format(result_path,))
    os.remove("{}pipo.in".format(result_path))


    return parsed_res


def parse_pairs_result(res, length):
    final = []
    for r in res:
        r = r.strip('\n')
        if not r.startswith('%'):
            r = r.split('\t')
            if len(r) == 3:
                if r[1] == str(length+1):
                    final.append(float(r[2]))

    return final

def nupack_analysis(sequence, window, result_path):
    list_for_table = []

    processed_sequence = sequence.upper().replace('T', 'U').replace(' ', '')
    loop_sequences = split_sequence(processed_sequence, window)
    
    target_toehold_map = possible_toehold_B(loop_sequences)
    
    sequence = sequence.upper().replace('T', 'U')  #this seems like repetition
    single_streadness_sequence = single_streadness(sequence, result_path, wait=6)
    for target in target_toehold_map.items():
        id = sequence.index(target)

        target_defect = sum(single_streadness_sequence[id:id+15])/15

        score = 5*(1-target_defect)

        list_for_table.append(tuple([target[0:15], 1-target_defect]))

    return list_for_table
