import pandas as pd
import os
import re
from collections import Counter
import warnings
from typing import Optional,Tuple
import multiprocessing as mp
import gc 

warnings.filterwarnings('ignore')

class AnnoPopFrq():
    def __init__(self,
                 vcf_chunk,
                 sample_grp
                ): 
        self.vcf_chunk = vcf_chunk
        self.sample_grp = sample_grp
        self.main()
    
    
    def fetch_ann(self, info):
        try:
            for ann in info.split('ANN=')[1].split(','):
                infos = ann.split('|')
                if len(infos)<6:
                    continue
                return (infos[3],infos[6],infos[1],infos[2],infos[9],infos[10])
        except:
            print(info)
    
    def filtered_and_anno(self):
        ## filter ALT == *
        self.vcf_chunk = self.vcf_chunk[~(self.vcf_chunk['ALT']== '*')]
        ## filter multi allele
        self.vcf_chunk = self.vcf_chunk[~(self.vcf_chunk['ALT'].str.contains(','))]

        ## get annotation info
        self.vcf_chunk = self.vcf_chunk.assign(**{'gene':'.', 'transcript':'.', 'functional':'.',
                                                  'hgv.c':'.','hgv.p':'.'})

        self.vcf_chunk.loc[:,['gene','transcript','functional','func_cate','hgv.c','hgv.p']] = \
        [i for i in map(self.fetch_ann, self.vcf_chunk['INFO'].values)]

    
    def call_pop_freq(self,gt_info):
        info_dict ={'1/1':0,'0/0':0,'./.':0,'1/0':0,'0/1':0}
        gt_info = [i.split(':')[0].replace('|','/') for i in gt_info]

        info_dict.update(Counter(gt_info))
        n = len(gt_info)*2
        hom_alt_freq = info_dict['1/1']*2/n
        het_alt_freq = (info_dict['0/1'] + info_dict['1/0'])/n
        #hom_ref_freq = info_dict['0/0']*2/n
        miss_freq = info_dict['./.']/n

        return (hom_alt_freq, het_alt_freq, miss_freq, ','.join(gt_info))
    
    def get_pop_freq(self,pop):
        samples = self.sample_grp[pop]
        samples = [str(i) for i in samples]
        samples_name = ','.join(samples)
        self.vcf_chunk = self.vcf_chunk.assign(**{f'{pop}.hom_alt.freq':0, 
                                        f'{pop}.het_alt.freq':0, 
                                        f'{pop}.miss.freq':0, 
                                        f'gts:{samples_name}':'.',
                                       })

        self.vcf_chunk.loc[:,[f'{pop}.hom_alt.freq',
                         f'{pop}.het_alt.freq',
                         f'{pop}.miss.freq', 
                         f'gts:{samples_name}']] = [i for i in map(self.call_pop_freq, 
                                                          self.vcf_chunk[samples].values)]
    
    
    def pop_hom_het_freq(self,popA,popB):
        ## filter poplation A and B with hom freq == 1    
        self.vcf_chunk = self.vcf_chunk[~((self.vcf_chunk[f'{popA}.hom_alt.freq'] == 1) | 
                                          (self.vcf_chunk[f'{popB}.hom_alt.freq'] == 1))]
        
        self.vcf_chunk['HET_AB'] = (self.vcf_chunk[f'{popA}.het_alt.freq']*(1 - self.vcf_chunk[f'{popB}.het_alt.freq']))
        self.vcf_chunk['HOM_AB'] = (self.vcf_chunk[f'{popA}.hom_alt.freq']*(1 - self.vcf_chunk[f'{popB}.hom_alt.freq']))

        self.vcf_chunk['HET_BA'] = (self.vcf_chunk[f'{popB}.het_alt.freq']*(1 - self.vcf_chunk[f'{popA}.het_alt.freq']))
        self.vcf_chunk['HOM_BA'] = (self.vcf_chunk[f'{popB}.hom_alt.freq']*(1 - self.vcf_chunk[f'{popA}.hom_alt.freq']))

    
    def main(self):
        ## filter and anno
        self.filtered_and_anno()
        popList = list(self.sample_grp.keys())
        
        ## get each population hom and het freq
        for pop in popList:
            self.get_pop_freq(pop)
            
        #if len(popList)>=2:
        #    self.pop_hom_het_freq(popList[0],
        #                          popList[1])



def chunk_handle(chunk,sample_info,apd,outdir):     
    chunk = AnnoPopFrq(chunk, sample_info).vcf_chunk
    gene_idx = list(chunk.columns).index('gene')
    if apd: 
        tmp_chunk = pd.concat([chunk[['#CHROM','POS','REF','ALT']],chunk.iloc[:,gene_idx:]],axis=1)
        tmp_chunk.to_csv(f'{outdir}/gt_freq_info.tsv',
                         mode='a',
                         index=False,                                                           
                         sep='\t',
                         header=None)
        
    else:
        tmp_chunk = pd.concat([chunk[['#CHROM','POS','REF','ALT']],chunk.iloc[:,gene_idx:]],axis=1)
        tmp_chunk.to_csv(
                     f'{outdir}/gt_freq_info.tsv',
                     index=False,
                     sep='\t',
                    )
    del chunk
    gc.collect()
