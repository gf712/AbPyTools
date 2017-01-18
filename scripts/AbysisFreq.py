#!usr/bin/env python

import re
import urllib2
import numpy as np
import pandas as pd
import StringIO
from tqdm import tqdm

def download(url, user_agent='wswp', num_retries=2, verbose=True):
    
    """
    Function to download contents from a given url
    
    Input:
            url: str
            string with the url to download from
            
            user_agent: str
            Default 'wswp'
            
            num_retries: int
            Number of times to retry downloading
            if there is an error
            
            verbose: bool
            Print out url and errors
            
    Output:
            returns: str
            string with contents of given url
    """
    
    error=False
    if verbose:
        print 'Downloading:', url
    headers = {'User-agent': user_agent}
    request = urllib2.Request(url, headers=headers)
    try:
        html = urllib2.urlopen(request).read()
    except urllib2.URLError as e:
        if verbose:
            print 'Download error:', e.reason
        html = None
        if num_retries > 0:
            if hasattr(e, 'code') and 500 <= e.code < 600:
                # retry 5XX HTTP errors
                return download(url, user_agent, num_retries-1)[0]
            #elif hasattr(e, 'code') and e.code == 404:
            else:
                error=True
    
    return html, error

def abysis_to_pandas(numbers):
    
    """
    Function to extract the amino acid frequency
    for a number of given positions by querying
    the Abysis server
    
    Input:
            numbers: list
            list of strings with the name of each
            position to obtain amino acid frequency
            from. Currently only added support for the
            Clothia numbering scheme.
            
    Output:
            returns: pandas.Dataframe
            table with a column for each position in
            numbers containing information about the
            count and frequency of each amino acid
            
    """
    
    # instantiate table_list as a list to store
    # a pandas.Dataframe for each Ab position
    table_list = list()
    
    # iterate through each position in the list numbers
    for position in tqdm(numbers):
        # first download the page
        query_result = download('http://www.bioinf.org.uk/abysis2.7/ws/resfreq.cgi?residue={}'.format(position),
                                verbose=False)[0]

        # for ease just read it in with pandas -> amino acid, count, %
        table = pd.read_csv(StringIO.StringIO(query_result),sep='\s+',skiprows=1)
        
        # add the table with the amino acid information for position
        # to table_list
        table_list.append(table)
    
    # crete pandas.MultiIndex object to name all the columns (position_i: ['AA','Count','%'])
    cols = pd.MultiIndex.from_product([numbers,['AA','Count', '%']])
    
    # return the pandas.Dataframe with the information for each position 
    return pd.DataFrame(np.column_stack(table_list),columns=cols)

def find_common_sequences_improved(seqs, aa_table,thresholds,tolerance=3):
    
    """
    Function to calculate the number of sequences 
    with common amino acids (i.e. frequency above
    elements in list thresholds). Additionally,
    the user can select a tolerance which corresponds
    to the maximum number of rare amino acids 
    tolerated in each sequence
    
    Input:
            seqs: list
            list of strings with the name of the
            sequence and the protein sequence itself
            separated by a single space ' '.
            
            aa_table: pandas.Dataframe
            table containing the amino acid frequency 
            at each position
            
            thresholds: list
            list of numbers (float or int) with the
            thresholds to be tested
            
            tolerance: int
            maximum number of rare amino acids per 
            sequence
    
    Output:
            threshold_dict: dict
            a dictionary with keys = thresholds as 
            strings and values that correspond to the 
            number of sequences for each threshold
    
    """
    
    # initiate variables
    common_sequences = list()
    threshold_dict = dict()

    for seq in tqdm(seqs):
        
        # parse sequence of format: "sequence_name sequence"
        # remove all the "-" as they cause an error with abnum
        seq_name = re.findall('\S+',seq)[0]
        sequence = re.findall('\S+',seq)[1].replace('-','')
        
        # generate the string to query abnum
        query = 'http://www.bioinf.org.uk/cgi-bin/abnum/abnum.pl?plain=1&aaseq={}&scheme={}'.format(sequence, '-c')
        # use the download function to obtain the text file with the sequence with the clothia numbering
        result, error = download(query, verbose = False)
        
        # result_list contains a list of strings of format: "position AminoAcid"
        # result_list has the same length as the sequence
        result_list = re.findall('[\S| ]+', result)
        
        # handle errors
        if result_list == ['Warning: Unable to number sequence']:
            continue
            
        # try each treshold and store result in threshold_dict
        # threshold_dict keys: threshold values as strings, e.g. "0.25"
        # threshold_dict values: number of sequences above threshold
        for threshold in thresholds:
            # initiate the tolerance value
            # the tolerance is the maximum value (inclusive) of unusual/rare amino acids in a sequence
            tolerance_i = 0
            
            # keep track if there is an unusual/rare amino acid
            found_rare_amino_acid = False
            
            # iterate through each amino acid contained in result_list
            for result_i in result_list:
                # extract the name of the position
                position_i = result_i[:-2]
                # find out the corresponding amino acid residue
                amino_acid_i = result_i[-1].upper()

                # skip CDRs
                #
                # aa_table is a pandas.Dataframe with the amino 
                # acid frequency at each non-CDR position
                #
                # therefore if the current position is not in aa_table
                # it will be part of a CDR region
                if position_i in aa_table.columns:
                    table_i = aa_table[position_i]
                else:
                    continue
                
                # find the row corresponding to the current amino acid
                # in table_i (which just contains information about the current position)
                aa_in_table = table_i['AA'] == amino_acid_i
                
                # find out if the frequency (in %) of the current amino acid 
                # is above the threshold
                #
                # TODO: clean up -> no need to sum
                # it will be either True or False, so 1 or 0
                # just put in if statement
                aa_above_threshold = (table_i['%'][aa_in_table] > threshold).sum()
                
                # if the amino acid freq is above threshold continue iterating 
                # over each amino acid in the sequence
                if aa_above_threshold > 0:
                    continue
                # if the previous statement is false it means that we move to the
                # following conditionals
                #
                # so if the number of rare amino acids we have found so far is below
                # the tolerance we continue iterating and add 1 to the number of 
                # rare amino acids found (which is tolerance_i)
                elif tolerance_i < tolerance:
                    tolerance_i += 1
                # otherwise we break from the loop and set found_rare_amino_acid to True
                else:
                    found_rare_amino_acid = True
                    break
            
            # if found_rare_amino_acid is False the current sequence is below the tolerance
            # and we keep track of the number of sequence below the threshold in treshold_dict
            if found_rare_amino_acid == False:
                if str(threshold) in threshold_dict:
                    threshold_dict[str(threshold)].append(seq)
                else:
                    threshold_dict[str(threshold)] = list()
    
    # after iterating through each sequence and threshold the function returns the threshold_dict
    # containing all the number of sequences for each threshold
    return threshold_dict