#!/data1/linhua/software/anaconda2/envs/hic/bin/python
# -*- coding: utf-8 -*-

## Pls Use Python3 with cooler installed
## Linhua Sun Modified the code at Thu Oct 24 20:57:07 CST 2019
## https://github.com/KaplanLab/Spermatogenesis
## Change: avgloglike = a*np.sum(Dnormed*log_dists) - sp.misc.logsumexp(a*log_dists) ==>> avgloglike = a*np.sum(Dnormed*log_dists) - sp.special.logsumexp(a*log_dists)
## It is useful to calculate scalings coee
## Example: python M_calc_scale_coeff.py -i Control_10000.cool -o test -show -ex_chr Chr1 -slope_range 20000 5000000

 
import argparse
import numpy as np

parser=argparse.ArgumentParser(description='Estimate power-law scaling coefficient of Hi-C matrix using maximum likelihood. Also calculates scaling plot.',formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('-i',help='input matrix (cool format)',dest='in_mat',type=str,required=True)
parser.add_argument('-o',help='output prefix (if supplied, saves scaling plot and estimated coefficient)',dest='out_prefix',type=str)
parser.add_argument('-head',help='Run headless (no terminal)',dest='headless',action='store_true')
parser.add_argument('-show',help='Show interactive plots',dest='show',action='store_true')
parser.add_argument('-ex_chr',help='exclude these chromosomes from analysis',dest='excluded_chrs',type=str,nargs="+",default=['chrM'])
parser.add_argument('-slope_range',help='calculate slope only for interactions within this distance range (give min and max, otherwise all distances>0). If 0 0 is given, the program will stop to prompt the user for the range after showing the scaling plot. ',dest='slope_range',type=float,nargs=2,default=[1,np.inf])

args=parser.parse_args()

in_mat = args.in_mat
out_prefix = args.out_prefix
show = args.show
headless=args.headless
excluded_chrs=args.excluded_chrs
slope_range=args.slope_range

import matplotlib
if headless:
    matplotlib.use('Agg')
import scipy as sp
import scipy.optimize
import scipy.misc
import cooler
import matplotlib.pyplot as plt
import pandas as pd

def main():
    global slope_range

    in_cool = cooler.Cooler(in_mat)
    chroms = in_cool.chromnames
    for i in excluded_chrs:
        chroms.remove(i)

    res=in_cool.binsize
   
    cis_vals=[]
    cis_dists=[]
    for c in chroms:
        chrom_ext = in_cool.extent(c)
        D = in_cool.matrix()[chrom_ext[0]:chrom_ext[1],chrom_ext[0]:chrom_ext[1]]

        n = D.shape[0]

        for i in range(n):

            x = np.diag(D,i)
            x = x[~np.isnan(x)]
            cis_vals += [x]
            cis_dists += [np.ones(x.shape[0])*i*res]

    cis_vals = np.concatenate(cis_vals)
    cis_dists = np.concatenate(cis_dists)

    n = cis_vals.shape[0]

    # scaling plot

    p = pd.DataFrame({'dist':cis_dists,'val':cis_vals})

    pmeans=p.groupby('dist',as_index=False).median()

    plt.loglog(pmeans['dist'],pmeans['val'],"o-")


    if out_prefix:
        np.save(out_prefix+"_scaling_plot.npy",pmeans)
        plt.savefig(out_prefix+"_scaling_plot.png")

    if show:
        if headless:
            figshow()
        else:
            plt.show()



    # estimate slope

    if slope_range == [0,0]:
        slope_range = list(map(float,input('please enter [min max] slope range: ').split(' ')))
    valid = (cis_dists>=slope_range[0]) & (cis_dists<=slope_range[1])
    Dnormed = cis_vals[valid]/np.sum(cis_vals[valid])
    log_dists = np.log(cis_dists[valid])

    res=sp.optimize.minimize(f,x0=-1.0,args=(Dnormed,log_dists),method="L-BFGS-B",bounds=[(None,0)])

    if out_prefix:
        with open(out_prefix+'_scaling_coeff.tab','w') as fh:
            print ("slope_range_min_dist",slope_range[0],"\nslope_range_max_dist",slope_range[1],"\nalpha",res['x'][0],file=fh)
    else:
        print("alpha = ",res['x'][0])



def f(a,Dnormed,log_dists):

    avgloglike = a*np.sum(Dnormed*log_dists) - sp.special.logsumexp(a*log_dists)
    return -avgloglike



main()

