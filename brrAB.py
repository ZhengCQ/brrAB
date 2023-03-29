#!/usr/bin/env  python
# -*- coding:UTF-8 -*-
# @Author: Zheng ChenQing
# @Date: 2023.03.19
# @E-mail: zhengchenqing@qq.com


import argparse
import pandas as pd
import os
import re
import time
from collections import Counter
import warnings
from typing import Optional,Tuple
import multiprocessing as mp
from bin.handle_data import AnnoPopFrq


ARGS = argparse.ArgumentParser(description="物种相对突变负荷AB法")

ARGS.add_argument(
	'-v', '--vcf', dest='vcf', required=True, help='snpEff注释后的vcf')
ARGS.add_argument(
	'-w', '--work_dir', dest='work_dir', default='.', help='工作目录，默认：当前目录。')

ARGS.add_argument(
	'-A', dest='A_population', required=True, help='A population的名称')
ARGS.add_argument(
	'-B', dest='B_population', required=True, help='B population的名称')
ARGS.add_argument(
	'-C', dest='B_population', required=False, help='C population的名称,外群')
ARGS.add_argument(
	'-G', dest='group_info', required=True, help='groups\tsample\n')
ARGS.add_argument(
    '--n_cores', dest='n_cores', type=int,default=4,
    help='多进程数目, 默认为4')

def chunk_handle(chunk,sample_info,apd,vcfout,outdir):     
    chunk = AnnoPopFrq(chunk, sample_info).vcf_chunk
    if apd: 
        if vcfout:
            chunk.iloc[:,0:20].to_csv(f'{outdir}/{outvcf_fi}',sep='\t',index=False,header=None,mode='a') 
        pd.concat([chunk[['#CHROM','POS','REF','ALT']],chunk.iloc[:,20:]],axis=1).to_csv(f'{outdir}/gt_freq_info.tsv',
                         mode='a',
                         index=False,                                                           
                         sep='\t',
                         header=None)
    else:
        if vcfout:
            outvcf_fi = fi.replace('.vcf.gz','.recode.vcf')
            outf = open(f'{outdir}/{outvcf_fi}','w')
            for i in header_lines:
                outf.write(i)
            outf.close()
            chunk.iloc[:,0:20].to_csv(f'{outdir}/{outvcf_fi}',sep='\t',index=False,header=None,mode='a')
        pd.concat([chunk[['#CHROM','POS','REF','ALT']],chunk.iloc[:,20:]],axis=1).to_csv(
                     f'{outdir}/gt_freq_info.tsv',
                     index=False,
                     sep='\t',
                    )

def main():
    args = ARGS.parse_args()
    if not args.vcf:
        print('Use --help for command line help')
        return
    try:
        os.makedirs(args.work_dir)
    except:
        print ('%s exists' %(args.work_dir))

    header_lines = os.popen(f'tabix -H {args.vcf}').readlines()
    
    vcf_dtypes = {'#CHROM':'category',
    'POS':'int32',
    'REF':'category',
    'ALT':'category',
    'FORMAT':'category',
    'FILTER':'category'}

    sample_info_r = pd.read_csv(args.group_info,sep='\t')
    sample_info = {}
    for idx,val in sample_info_r.groupby('Group'):
        for i in val['Sample']:
            sample_info.setdefault(idx,[]).append(i)
        

    reader = pd.read_csv(args.vcf, sep="\t",
                            compression='gzip',
                            skiprows=(len(header_lines) - 1),
                            iterator=True,
                            dtype=vcf_dtypes)
    loop = True
    chunkSize = 10000
    chunks = []
    is_apd = False
    is_vcfout=False
    pool = mp.Pool(args.n_cores) #启动多进程池
    while loop:
        try:
            chunk = reader.get_chunk(chunkSize)
            if is_apd:
                #chunk_handle(chunk,sample_info,is_apd,is_vcfout,args.work_dir)
                try:
                    pool.apply_async(chunk_handle,(chunk,sample_info,is_apd,is_vcfout,args.work_dir))
                except:
                    print('error')
                
            else:
                #chunk_handle(chunk,sample_info,is_apd,is_vcfout,args.work_dir)
                try:
                    pool.apply_async(chunk_handle,(chunk,sample_info,is_apd,is_vcfout,args.work_dir))
                except:
                    print('error')
                #pool.apply_async(chunk_handle,(chunk,sample_info,is_apd,is_vcfout))
                #pool.map(chunk, (chunk,sample_info,is_apd,is_vcfout))
                #pool.apply_async(chunk_handle, args=(chunk,sample_info,apd))
                is_apd=True
        except StopIteration:
            loop = False
            print("Iteration is stopped.")
    print('Waiting for all subprocesses done...')
    pool.close()
    pool.join()
    print('All subprocesses done.')
    pool.terminate()

if __name__ == '__main__':
	start = time.time()
	main()
	end = time.time()
	print ("时间总计%s"%(end - start))
