from flask import Flask, render_template, request
import ast
from bs4 import BeautifulSoup
import requests as rq
import os
import os.path
from subprocess import Popen, PIPE
import time
import pandas as pd

app = Flask(__name__)

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
        sequences.append(sequence[i:(window+i)]) 

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


def possible_toehold_B(loop):
    #fill in

def nupack_analysis(sequence, window, result_path):
    list_for_table = []

    processed_sequence = sequence.upper().replace('T', 'U').replace(' ', '')
    d_1 = {'Triggers': split_sequence(processed_sequence, window)}  
    result = pd.DataFrame(data=d_1)
    
    for i in range(len(result.iloc[:,0])): 
        list_for_table.append((possible_toehold_B(reversed_complement(result.iloc[i,0])))[0])

    return list_for_table


@app.route('/', methods=['GET'])
def index():
    return render_template('index.html')


@app.route('/', methods=['POST'])
def toeholds():
    sequence=request.form.get('sequence')
    list = []
    if sequence != None:
        window = 18  
        result_path = ''
        list = sorted(nupack_analysis(sequence, window, result_path), key=lambda x: x[5])

    return render_template('toeholds.html', list=list)


@app.route('/home', methods=['POST', 'GET'])
def home():
    l = request.form.get('list')
    list=[]
    if l != None:
        list = ast.literal_eval(l)
    return render_template('toeholds.html', list=list)


@app.route('/details', methods=['POST'])
def details():
    list = ast.literal_eval(request.form.get('list'))
    index = int(request.form.get('index')) % len(list)
    element = list[index]

    target = element[0]
    toehold = element[1]
    target_ss = element[2]
    toehold_ss = element[3]
    defect = element[4]
    score = element[5]

    return render_template('details.html', details=(target, toehold, target_ss,toehold_ss, defect, score, list, index))


@app.route('/structure', methods=['POST'])
def structure():
    seq = request.form.get('sequence')
    details = ast.literal_eval(request.form.get('details'))
    d = {'commit': 'Analyze',
         'partition_job[dangle_level]': '1',
         'partition_job[dna_parameter_file]': 'dna1998',
         'partition_job[dotplot_target]': '0',
         'partition_job[email_address]': '',
         'partition_job[filter_max_number]': '',
         'partition_job[filter_min_fraction_of_max]': '',
         'partition_job[is_melt]': '0',
         'partition_job[max_complex_size]': '1',
         'partition_job[max_melt_temperature]': '',
         'partition_job[melt_temperature_increment]': '',
         'partition_job[mg_salt]': '0.0',
         'partition_job[min_melt_temperature]': '',
         'partition_job[na_salt]': '1.0',
         'partition_job[nucleic_acid_type]': 'RNA',
         'partition_job[num_sequences]': '1',
         'partition_job[predefined_complexes]': '',
         'partition_job[pseudoknots]': '0',
         'partition_job[rna_parameter_file]': 'rna1995',
         'partition_job[temperature]': '37.0',
         'partition_sequence[0][concentration]': '',
         'partition_sequence[0][contents]': '',
         'partition_sequence[0][name]': 'strand1',
         'partition_sequence[0][scale]': '-6'}

    d['partition_sequence[0][contents]'] = seq

    analyze = "/partition/new"
    base_url = 'http://www.nupack.org'

    page = rq.post(base_url + analyze, d)
    soup = BeautifulSoup(page.text, "html.parser")
    result = soup.find('a', {'title': 'Click to see results'})['href']

    img = None
    while img is None:
        page = rq.get(base_url + result)
        soup = BeautifulSoup(page.text, 'html.parser')
        img = soup.find('div', id='fullsize')

    img = img.find('img')['src']

    return render_template('details.html', img=base_url+img, details=details, title=seq)

if __name__ == '__main__':
    app.run()
