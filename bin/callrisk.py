import pandas as pd
import re
import json
from typing import Optional,Tuple
import numpy as np
import random
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
import os
warnings.filterwarnings('ignore')

class CallfunctionalFi():
    def __init__(self,
                 df,
                 fix_sites: Optional[int] = None,
                 boostrap_n: Optional[int] = 100
        ) -> Optional[Tuple[np.ndarray, np.ndarray]]: 
        """\
        CallfunctionalFi(df)
        """
        
        if fix_sites:
            self.fix_sites = fix_sites
        else:
            self.fix_sites = int(df.shape[0]/5)
        self.boostrap_n = boostrap_n
        self.main(df)

        
    def trans_data(self,df):
        ## filtered frequency
        #df_di = pd.read_csv(file,sep='\t')
       # df[f'{self.popA}_pop_freq'] = df[f'{self.popA}.het_alt.freq'] + df[f'{self.popA}.het_alt.freq']
       # df[f'{self.popB}_pop_freq'] = df[f'{self.popB}.het_alt.freq'] + df[f'{self.popB}.het_alt.freq']
        
        #df_di = df[(df[f'{self.popA}_pop_freq'] <= self.max_freq) & (df[f'{self.popA}_pop_freq'] <= self.max_freq)]
        
        ## category delerious and ben 
        df['Gscores'] = df['Gscores'].fillna(0)
        df['Gscores'] = df['Gscores'].apply(lambda x:int(x))
        df.loc[((df['functional']=='missense_variant') & (df['Gscores']>=150 )),'functional'] = 'missense_del'
        df.loc[((df['functional']=='missense_variant') & (df['Gscores']<150 )),'functional'] = 'missense_ben'
        return df
    
    def call_spec_pop_freq(self, df):
        ## pop specific frequency 

        df_fi_AB = pd.DataFrame()
        df_fi_AB['All'] = (df['HET_AB'] + df['HOM_AB']).groupby(df['functional']).sum()
        df_fi_AB['HET'] = df['HET_AB'].groupby(df['functional']).sum()
        df_fi_AB['HOM'] = df['HOM_AB'].groupby(df['functional']).sum()
        
        df_fi_BA = pd.DataFrame()
        df_fi_BA['All'] = (df['HET_BA'] + df['HOM_BA']).groupby(df['functional']).sum()
        df_fi_BA['HET'] = df['HET_BA'].groupby(df['functional']).sum()
        df_fi_BA['HOM'] = df['HOM_BA'].groupby(df['functional']).sum()
        
        return df_fi_AB,df_fi_BA
    
    
    def _jackknifes(self,fix_sites, all_sites):
        """
        @param fix_sites: 固定位点数 
        @param all_sites: 总的数目
        """
        start = random.randint(0, all_sites) # get start from 0 and all_sites
        end = start + fix_sites
        if end > all_sites:
            start, end = self._jackknifes(fix_sites, all_sites) #重新得到的start, end
        return start, end
    
    #def get_regs(self,df):
        #regs = []
        #for idx,val in df[['#CHROM','POS']].groupby('#CHROM'):
         #   regs.append(f'{idx}:{min(val["POS"])}-{max(val["POS"])}')
        #return ','.join(regs)
    
    def _boostrap_run(self,df):
        n = 0
        fi_AB_infos = []
        fi_BA_infos = []
        #reg_n = []
        while n<self.boostrap_n:
            n = n + 1
            start,end = self._jackknifes(self.fix_sites, df.shape[0])
            df_di_n = df.iloc[start:end,]
            #reg_n.append([n,self.get_regs(df_di_n)])
            df_fi_AB,df_fi_BA = self.call_spec_pop_freq(df_di_n)
            df_fi_AB.columns = [i+ f'_{n}' for i in df_fi_AB.columns]
            df_fi_BA.columns = [i+ f'_{n}' for i in df_fi_BA.columns]
    
            fi_AB_infos.append(df_fi_AB)
            fi_BA_infos.append(df_fi_BA)
            
        df_fi_AB_bs  = pd.concat(fi_AB_infos,axis=1).fillna(0)
        df_fi_BA_bs  = pd.concat(fi_BA_infos,axis=1).fillna(0)
        #df_reg_bs =  pd.DataFrame(reg_n,columns=['N','Regions'])
        
        return df_fi_AB_bs,df_fi_BA_bs
    
    def main(self,df):
        self.df_di = self.trans_data(df)
        if self.boostrap_n > 0:
            self.AB, self.BA  = self._boostrap_run(self.df_di) 
        else:
            self.AB, self.BA = self.call_spec_pop_freq(self.df_di)


class CallBurdenRisk():
    def __init__(self,
             fi_AB,
             fi_BA,
             norm_item: Optional[str] = 'intergenic_region',
    ) -> Optional[Tuple[np.ndarray, np.ndarray]]:
        self.fi_AB = fi_AB
        self.fi_BA = fi_BA
        self.norm_item = norm_item
        self.main()
        
    def call_norm_risk(self, df_fi_AB, df_fi_BA, in_item, norm_item):
        if not isinstance(in_item,list):
            in_item = [in_item]
        if not isinstance(norm_item,list):
            norm_item = [norm_item]

        AB = (df_fi_AB.loc[list(set(in_item) & set(df_fi_AB.index))].sum()/df_fi_AB.loc[norm_item].sum())
        BA = (df_fi_BA.loc[list(set(in_item) & set(df_fi_BA.index))].sum()/df_fi_BA.loc[norm_item].sum())
        return AB/BA
    
    def norm_risk(self,
                  df_fi_AB,
                  df_fi_BA):
        lof_items=['start_lost','splice_acceptor_variant','stop_gained','stop_lost','splice_donor_variant']
        df_risk = pd.DataFrame()
        df_risk['missense_del'] = self.call_norm_risk(df_fi_AB,df_fi_BA,'missense_del',self.norm_item)
        df_risk['missense_ben'] = self.call_norm_risk(df_fi_AB,df_fi_BA,'missense_ben',self.norm_item)
        df_risk['LOF'] = self.call_norm_risk(df_fi_AB,df_fi_BA,lof_items,self.norm_item)
        return df_risk
    
    def main(self):
        self.df_risk = self.norm_risk(self.fi_AB,
                                      self.fi_BA).reset_index()
        try:
            self.df_risk[['Group','Bs_N']] = self.df_risk['index'].apply(lambda x:x.split('_')).apply(pd.Series)
        except:
            self.df_risk['Group'] = self.df_risk['index'].apply(lambda x:x)


def get_gscores(x,Gscores):
    aa1,aa2 = re.search(r'p\.(\w{3})\d+(\w{3})',x).groups()[:]
    try:
        missense_score = Gscores[aa1][aa2]
    except:
        print(x)
        missense_score = Gscores[aa2][aa1]
    return missense_score


def plot_burden_risk(df,kind='boxplot'):
    df_mt = df[['missense_del','missense_ben','LOF','Group']].melt(id_vars=['Group'])
    plt.rcParams['xtick.labelsize']=20
    plt.rcParams['ytick.labelsize']=20
    plt.rcParams['legend.fontsize']=16
    plt.rcParams['legend.title_fontsize']=16
    plt.rcParams['font.family'] =  'Helvetica'

    fig,ax = plt.subplots(1,1,figsize=(6,6))
    if kind == 'boxplot':
        g = sns.boxplot(x='Group',y='value',hue='variable'
                        ,data=df_mt,
                        hue_order=['LOF','missense_del','missense_ben'],
                        ax=ax,palette=['red','yellow','lightgreen'])
    elif kind == 'barplot':
        g = sns.barplot(x='Group',y='value',hue='variable'
                        ,data=df_mt,
                        hue_order=['LOF','missense_del','missense_ben'],
                        ax=ax,palette=['red','yellow','lightgreen'])

    handles, labels = ax.get_legend_handles_labels()
    plt.legend(handles=handles, labels = ['LOF','Missense deleterious','Missense benign'],
                          bbox_to_anchor=(1,0.65),frameon=False)
    g.set_ylabel('Relative Burden Risk$_{A/B}$',fontsize=20,labelpad=10)
    g.set_xlabel('')
    
    plt.axhline(y=1,ls='--',color='r')
    sns.despine()

def derived_allele(df, outgrp):
    if isinstance(outgrp,str):
        outgrp = [outgrp]
    flag_info = [(df[f'{grp}.hom_alt.freq'] == 0)  & (df[f'{grp}.het_alt.freq'] == 0) for grp in outgrp]

    is_di = flag_info[0]
    if len(flag_info) > 1:
        for each in flag_info[1:]:
            flag = np.logical_or(is_di,each)
    df_info_di = df[is_di]
    return df_info_di

def add_Gscores(df):
    binpath = os.path.split(os.path.realpath(__file__))[0]

    Gscores = json.load(open(f'{binpath}/Grantham_Scores.json','r'))
    df.loc[df['functional'] == 'missense_variant','Gscores'] = \
    df.loc[df['functional'] == 'missense_variant','hgv.p'].apply(lambda x:get_gscores(x,Gscores))
    return df


def pop_hom_het_freq(df,popA,popB):
    ## filter poplation A and B with hom freq == 1    
    df = df[~((df[f'{popA}.hom_alt.freq'] == 1) | 
                                      (df[f'{popB}.hom_alt.freq'] == 1))]

    df['HET_AB'] = df[f'{popA}.het_alt.freq']*(1 - df[f'{popB}.het_alt.freq'])
    df['HOM_AB'] = df[f'{popA}.hom_alt.freq']*(1 - df[f'{popB}.hom_alt.freq'])

    df['HET_BA'] = df[f'{popB}.het_alt.freq']*(1 - df[f'{popA}.het_alt.freq'])
    df['HOM_BA'] = df[f'{popB}.hom_alt.freq']*(1 - df[f'{popA}.hom_alt.freq']) 
    return df