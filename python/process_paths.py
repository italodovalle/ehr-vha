#! /usr/bin/env python

import pandas as pd
import random
import sys
import argparse
from collections import defaultdict

def create_edge_list(infile, bigram_only = True):
    """Create edge lists from paths
    
    Create bigram and trigram edge lists.
    Bigram: each node is a disease and the links is one sequential step
    (i.e. disease_i -> disease j
    Trigram: each node is a pair of diseases and each link represent
    two sequential steps (i.e. disease_i -> disease_j -> disease_k).
    
    
    Args:
        infile (str): path to input file
    
    Returns:
        bigram (pd.DataFrame): bigram edge list
        trigram (optional, pd.DataFrame): trigram edge list 
    
    """
    
    data = open(infile, 'r').readlines()
    
    bigram = defaultdict(int)
    trigram = defaultdict(int)
    line = random.sample(data,1)[0]
    for line in data:
        seq = line.rstrip().split('-')
        edges = []
        for i in range(1,len(seq)):
            edge = '%s-%s'%(seq[i-1],seq[i])
            edges.append(edge)
            bigram[edge] += 1
        if not bigram_only:
            for i in range(1, len(edges)):
                st1 = edges[i-1]
                st1 = st1.replace('-',';')
                st2 = edges[i]
                st2 = st2.replace('-',';')
                state = '%s-%s'%(st1,st2)
                trigram[state] += 1
            
    bigram = pd.DataFrame.from_dict(bigram,orient='index')
    if not bigram_only:
        trigram = pd.DataFrame.from_dict(trigram,orient='index')
        return (bigram, trigram)
    else:
        return (bigram)
    
def format_edge_list(df):
    """Format output from create_edge_list function"""
    
    ## formating edge list
    nodes = list(df.index)
    df = df.reset_index()
    df.columns = ['index', 'count']
    df['source'] = [i.split('-')[0] for i in nodes]
    df['target'] = [i.split('-')[1] for i in nodes]
    
    return (df)


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Get edge lists')
    parser.add_argument('-infile', dest = 'infile',action='store',
                        help='infile', required=True)
    parser.add_argument('-b', dest='bigram', action='store',
                        help='output bigram', required=True)
    parser.add_argument('-t', dest='trigram', action='store',
                        help='output trigram (optional)',required=False)
    
    
    
    if len(sys.argv) <= 1:
        print (parser.print_help())
        sys.exit(1)
    else:
        args = parser.parse_args()
    
    
    if args.trigram:
        bigram_only = False
        bigram, trigram = create_edge_list(args.infile,
                                           bigram_only)
        outputs = [bigram , trigram]
        outfiles = [args.bigram, args.trigram]
    else:
        bigram = create_edge_list(args.infile)
        outputs = [bigram]
        outfiles = [args.bigram]
    
    for df,outfile in zip(outputs, outfiles):
        df = format_edge_list(df)
        df.to_csv(outfile, index=None)
            
        
    
    