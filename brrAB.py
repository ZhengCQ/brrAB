#!/usr/bin/env  python
# -*- coding:UTF-8 -*-
# @Author: Zheng ChenQing
# @Date: 2023.03.19
# @E-mail: zhengchenqing@qq.com


import argparse
import pandas as pd
import os
import time
import multiprocessing as mp
from bin.handle_data import chunk_handle
from bin.callrisk import *
import gc


ARGS = argparse.ArgumentParser(description="物种相对突变负荷AB法")
ARGS.add_argument(
	'-f', '--infile', dest='infile', required=True, help='snpEff注释后的vcf 或vcf解析后的gt_info_freq文件')
ARGS.add_argument(
	'-w', '--work_dir', dest='work_dir', default='.', help='工作目录，默认：当前目录。')
ARGS.add_argument(
	'-A', dest='A_population', required=True, help='A population的名称')
ARGS.add_argument(
	'-B', dest='B_population', required=True, help='B population的名称')
ARGS.add_argument(
	'-C', dest='C_population', required=False, help='C population的名称,外群')
ARGS.add_argument(
	'-G', dest='group_info', required=True, help='groups\tsample\n')
ARGS.add_argument(
	'--freq', dest='freq', type=float, help='population allele frequency,default is no filter ')
ARGS.add_argument(
	'--fix_sites', dest='fix_sites', type=int, help='fix sites for jackknifes,default is 1/5 of total sites')
ARGS.add_argument(
    '--n_cores', dest='n_cores', type=int,default=4,
    help='多进程数目, 默认为4')

 
def vcf2gtfreq(vcf,n_cores,sample_info,work_dir):
    header_lines = os.popen(f'tabix -H {vcf}').readlines()
    
    vcf_dtypes = {'#CHROM':'category',
    'POS':'int32',
    'REF':'category',
    'ALT':'category',
    'FORMAT':'category',
    'FILTER':'category'}

    reader = pd.read_csv(vcf, sep="\t",
                            compression='gzip',
                            skiprows=(len(header_lines) - 1),
                            iterator=True,
                            dtype=vcf_dtypes)
    loop = True
    chunkSize = 10000
    is_apd = False
    #is_vcfout=False
    pool = mp.Pool(n_cores) #启动多进程池
    while loop:
        try:
            chunk = reader.get_chunk(chunkSize)
            if is_apd:
                #chunk_handle(chunk,sample_info,is_apd,work_dir)
                try:
                    pool.apply_async(chunk_handle,(chunk,sample_info,is_apd,work_dir))
                except:
                    print('error')
                
            else:
                #chunk_handle(chunk,sample_info,is_apd,work_dir)
                try:
                    pool.apply_async(chunk_handle,(chunk,sample_info,is_apd,work_dir))
                except:
                    print('error')
                is_apd=True
            
            del chunk
            gc.collect()
        except StopIteration:
            loop = False
            print("Iteration is stopped.")
    print('Waiting for all subprocesses done...')
    pool.close()
    pool.join()
    print('All subprocesses done.')
    pool.terminate()


def read_gtfreq(gt_freq_info,sample_info):
    ## cols
    freq_dtypes = {'#CHROM':'category',
                'POS':'int32',
                'functional':'category',
                'func_cate':'category',
                'gene':'category',
                }
    target_cols = ['#CHROM','POS', 'gene','functional','func_cate','hgv.p']
    for i in sample_info.keys():
        freq_dtypes.update({f'{i}.hom_alt.freq':'float16'})
        freq_dtypes.update({f'{i}.het_alt.freq':'float16'})
        target_cols.append(f'{i}.hom_alt.freq')
        target_cols.append(f'{i}.het_alt.freq')

    df_info = pd.read_csv(gt_freq_info, sep='\t',usecols= target_cols,low_memory=False)
    ## 多个功能注释，只取前面一个，默认为功能危害最高
    df_info['functional'] = df_info['functional'].apply(lambda x:x.split('&')[0])
    return df_info
    


def call_risk(df_info,workdir,popA,popB,outgrp,freq=None,fix_sites=None):
    ## derived allele
    df_info_di = derived_allele(df_info,outgrp)
    del(df_info)
    ## add Gscores to missense 
    df_info_di =  add_Gscores(df_info_di)
    ## cal AB,BA
    df_info_di = pop_hom_het_freq(df_info_di,popA,popB)

    if freq:
        df_info_di = df_info_di[((df_info_di[f'{popA}.hom_alt.freq'] + df_info_di[f'{popA}.het_alt.freq'])<freq) & \
        ((df_info_di[f'{popB}.hom_alt.freq'] + df_info_di[f'{popB}.het_alt.freq'])<freq)]

    funcfi1 = CallfunctionalFi(df_info_di,fix_sites=fix_sites,boostrap_n=100)
    df_risk = CallBurdenRisk(funcfi1.AB,
                            funcfi1.BA,
                            norm_item='intergenic_region').df_risk

    if not os.path.exists(f'{workdir}/riskAB'):
        os.makedirs(f'{workdir}/riskAB')
    if freq:
        df_risk.to_csv(f'{workdir}/riskAB/results_freq_{freq}.tsv',sep='\t',index=False)
        plot_burden_risk(df_risk)
        plt.savefig(f'{workdir}/riskAB/Burden_risk_freq_{freq}.pdf',bbox_inches='tight')
    else:
        df_risk.to_csv(f'{workdir}/riskAB/results.tsv',sep='\t',index=False)
        plot_burden_risk(df_risk)
        plt.savefig(f'{workdir}/riskAB/Burden_risk.pdf',bbox_inches='tight')


def main():
    args = ARGS.parse_args()
    print (args)
    if not args.infile:
        print('Use --help for command line help')
        return
    try:
        os.makedirs(args.work_dir)
    except:
        print ('%s exists' %(args.work_dir))

    sample_info_r = pd.read_csv(args.group_info,sep='\t')
    sample_info = {}
    for idx,val in sample_info_r.groupby('Group'):
        for i in val['Sample']:
            sample_info.setdefault(idx,[]).append(i)
    if args.infile.endswith(('vcf.gz','vcf')):
        print(f'Your input {args.infile} is vcf file, starting vcf 2 genotype frequency ......')
        vcf2gtfreq(args.infile,args.n_cores,sample_info,args.work_dir)
        print(f'Reading genotype frequency ......')
        df_info = read_gtfreq(f'{args.work_dir}/gt_freq_info.tsv', sample_info)
    elif args.infile.endswith(('gt_freq_info.tsv.gz','gt_freq_info.tsv')):
        print(f'Your input {args.infile} is genotype freqency file')
        print(f'Reading genotype frequency ......')
        df_info = read_gtfreq(args.infile,sample_info)
    else:
        print('Please input vcf or gt_freq_info.tsv file!!')
        print('Use --help for command line help')
        return
    print(f'Starting call risk ......')
    call_risk(df_info, args.work_dir, 
              args.A_population, 
              args.B_population, 
              args.C_population,
              fix_sites = args.fix_sites,
              freq = args.freq)

if __name__ == '__main__':
	start = time.time()
	main()
	end = time.time()
	print ("时间总计%s"%(end - start))
