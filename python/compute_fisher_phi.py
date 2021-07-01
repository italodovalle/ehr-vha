#! /usr/bin/env python

import pandas as pd
import numpy as np
import argparse
from collections import defaultdict
from progressbar import ProgressBar
from multiprocessing import Pool
import time
import networkx as nx
import os
import sys
import scipy.stats as stats
from statsmodels.stats.multitest import multipletests


class ComputeFisher():

    """Calculate phi coefficient and Fisher's exact test p-value
    for all the weighted edges in a graph"""
    
    def __init__(self, G, tmp_dir, ncpus = 10):
        
        """
        Args:
            G (nx.Graph): weighted graph
        """
        
        self.G = G
        self.S = self.G.size(weight='weight')
        pairs = list(self.G.edges())
        self.chunks = np.array_split(pairs, 100)
        self.tmp_dir = tmp_dir
        self.ncpus = 10
        
        
    def compute_fisher_directional(self, source, target):
        
        """
        Get phi coefficient and Fisher's Exact Test p-value
        Args:
            source (int): source node id
            target (int): target node id
        Returns:
            fe (tuple): Fisher's p-value and odds ratio
            phi_ij (float): phi coefficient
        """
        
        G = self.G
        S = self.S
        dij = G[source][target]['weight']
        ## i>any other
        di_ = (G.out_degree(source, weight='weight')) - dij 
        ## any other>j
        d_j = (G.in_degree(target, weight='weight')) - dij 
        ## any other, any other
        d_ = S - (G.out_degree(source, weight='weight')) - (G.in_degree(target, weight='weight')) + dij

        cont_table = [[dij,di_],[d_j,d_]]

        sj_in = G.in_degree(target, weight='weight')
        si_out = G.out_degree(source, weight='weight')
        phi_ij = ((dij * S) - (sj_in * si_out))/np.sqrt(si_out * sj_in  * (S - si_out) * (S - sj_in))

        fe = stats.fisher_exact(cont_table)
        return(fe,phi_ij)
    
    def run_parallel(self, pair):
        """Run compute_fisher_directional in parallel
        Args:
            pair (list): source and target node ids
        Returns:
            table (pd.DataFrame): fisher results
        """
        
        source, target = pair
        (odds,pvalue),phi_ij = self.compute_fisher_directional(source,target)
        d = (source,target,odds,pvalue,phi_ij)
        return(d)

    def prepare_table(self, res):
        """Organize results from run_parallel
        Args:
            res (list):
        Returns:
            table (pd.DataFrame): formatted results
        """
        table = defaultdict(dict)
        for i, elem in enumerate(res):
            table[i]['source'], table[i]['target'], table[i]['odds_ratio'], table[i]['pvalue'], table[i]['phi_ij'] = elem
        table = pd.DataFrame.from_dict(table, orient='index')
        return(table)

    
    def run(self):

        """Run Calculation"""
        
        pbar = ProgressBar()
        c = 0
        for chunk in pbar(self.chunks):
            s = time.time()
            p = Pool(self.ncpus)
            res = p.map(self.run_parallel, chunk)
            p.close()
            dres = self.prepare_table(res)
            dres.to_csv(self.tmp_dir + 'chunk_%d.csv'%c, index=None)
            c = c + 1
            e = time.time() - s

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Compute phi and Fisher')
    parser.add_argument('-infile', dest = 'infile',action='store',
                        help='infile', required=True)
    parser.add_argument('-tmp_dir', dest='tmp_dir', action='store',
                        help='tmp dir', required=True)
    parser.add_argument('-outfile', dest='outfile', action='store',
                        help='outfile', default='out.csv')
    parser.add_argument('-ncpus', dest='n cpus', action='store',default = 10,
                        help='number of cpus')
    
    if len(sys.argv) <= 1:
        print (parser.print_help())
        sys.exit(1)
    else:
        args = parser.parse_args()

    dt = pd.read_csv(args.infile, index_col = 0)
    wedges = [(i,j,k) for i,j,k in zip(dt.source, dt.target, dt.counts)]
    G = nx.DiGraph()
    G.add_weighted_edges_from(wedges)
    print ('## Load Network ##')
    print ('## {0} nodes and {1} edges'.format(len(G.nodes()), len(G.edges())))

    
    compute = ComputeFisher(G, args.tmp_dir)
    compute.run()
    
    
    res = []
    files = os.listdir(os.path.abspath(args.tmp_dir))
    for fi in files:
        if fi.startswith('chunk'):
            dx = pd.read_csv(os.path.abspath(args.tmp_dir) + '/' + fi)
            res.append(dx)
    
    res = pd.concat(res)
    res['fdr_bh_pvalue_adj'] = multipletests(res['pvalue'],
                                             alpha=0.05, 
                                             method='fdr_bh')[1]
    res.to_csv(args.outfile, index=None)
    
    
    