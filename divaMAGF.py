import pandas as pd
import os
import matplotlib 
from matplotlib import pyplot as plt
import numpy as np
import random 

import random

def random_rgb():
    
    col = ["#"+''.join([random.choice('ABCDEF0123456789') for i in range(6)])]
    
    return col[0]


def load_experiment(DATA_DIR,experiment):
    
    print()
    print(3*" <",f' Loading: {experiment} experiment ','\n')

    B_DIR = os.path.join(DATA_DIR,experiment)
    
    for f in os.listdir(B_DIR):
        
        sheet = os.path.join(B_DIR,f)
        
        if os.path.isfile( sheet ):
            
            if 'scaffold' in f:

                scaffold=pd.read_csv(sheet,sep='\t')
                print(" >", 'scaffolds: ', f, ': loaded.')

            elif 'SNV' in f:

                snv=pd.read_csv(sheet,sep='\t')
                print(" >",'SNVs: ', f, ': loaded.')
            
    return scaffold,snv

def MAG2SNV(scaffold2bin,snv,scaffold,OUT_DIR,window_size,step_size,pc_threshold,which_MAGs=None):
    
    print()
    print(3*" <",f' SNV frequency calculation ','\n')

    w,d=window_size,step_size
    # coordinate relative to the entire mag, coord. rel. to the scaffold, freq of SNVs, scaffoldIds
    atlas_cols = ['genomic_coordinate','scaffold_coordinate','SNV_density','scaffold']
    
    scaffold_size = scaffold[['scaffold','length']].set_index('scaffold')['length']

    

    snv = snv[ snv['position_coverage']>pc_threshold ]

    #sf = snv['scaffold'].copy()
    #snv['MAG']=snv['scaffold'].map( scaffold2bin['MAG'].to_dict() )
    
    s2m = (snv['scaffold'].map(scaffold2bin['MAG'].to_dict()))
    snv = snv.join( pd.DataFrame(index = s2m.index, data=s2m.values, columns=['MAG']),on=None )    

    if which_MAGs==None:
        MAGs = list( scaffold2bin['MAG'].unique() )
    else:
        MAGs=which_MAGs
    
    for m in MAGs:
    
        name = m.replace('.fa','')
        print(" >",f' Processing MAG: {name}')
        # get scaffolld-mag mapping
        s2b_MAG = scaffold2bin[scaffold2bin['MAG']==m]

        # get SNV infos about the current MAG
        MAG_snv = snv[snv['MAG']==m]

        # get the list of unique scaffolds belonging to a MAG. Sort them as to have an "order"
        mag_scaffolds = list( MAG_snv['scaffold'].unique() )
        
        # get a scaffold-2-lenght mapping
        scf_size = scaffold_size[mag_scaffolds]
        mag_size = scf_size.sum()
        
        #mag_scaffolds = list( snv_MAG['scaffold'].unique() )
        mag_scaffolds = scaffold_size.loc[ mag_scaffolds ].sort_values(ascending=False).index
        mag_scaffolds = list( mag_scaffolds )

        # genomic scaffold shift: every time we  introduce a new scaffold the genomic coordinate will be shifted
        gsf = 0

        MAG_SNV_density = pd.DataFrame(columns=atlas_cols)

        # for each scaffold in the current mag, do the histogram of SNVs

        for ms in mag_scaffolds:

            # scaffold size
            l=scf_size.loc[ms]

            # discard scaffolds shorther than the window
            if l > w:

                SCA2MAG_snv = pd.DataFrame(columns=atlas_cols,index=range(gsf,gsf+int( (l-w)/d )+2))

                SCA_snv = MAG_snv[ MAG_snv['scaffold']==ms]

                snvsca = SCA_snv['position'].copy()
                snvsca.index = snvsca.values

                var = np.zeros( l )
                var[ SCA_snv['position'].values ] = 1

                snv_density=pd.Series(dtype=float)

                snv_density.loc[0.5*w]=var[0:w].sum() 

                for i in range(int( (l-w)/d )):
                    snv_density.loc[ w+(i+0.5)*d ]=var[ w+i*d : w+(i+1)*d ].sum()

                snv_density.loc[ 0.5*(w+(i+1)*d+l) ]=var[ w+(i+1)*d: l ].sum() 

                SCA2MAG_snv['scaffold_coordinate']= snv_density.index.values 
                SCA2MAG_snv['genomic_coordinate']=snv_density.index.values+gsf
                SCA2MAG_snv['SNV_density']=snv_density.values
                SCA2MAG_snv['scaffold']=ms

                MAG_SNV_density=MAG_SNV_density.append(SCA2MAG_snv)

                gsf+=l

        MAG_DIR = os.path.join(OUT_DIR,m.replace('.fa',''))
        
        if not os.path.isdir(MAG_DIR): os.system(f"mkdir {MAG_DIR}")
                
        MAG_SNV_density.to_csv(os.path.join(MAG_DIR,m.replace('.fa','.csv')),sep='\t',index=False)

def SNVIsual(OUT_DIR,which_MAGs):
    
    print()
    print(3*" <",f' MAG-SNV visualization ','\n')

    for mag in which_MAGs:

        name = mag.replace('.fa','')

        print(" >",f' Plotting MAG: {name}')
             
        try:

            MAG = pd.read_csv( os.path.join(OUT_DIR,name,name+'.csv' ),sep='\t')   
            
            fig,ax=plt.subplots(figsize=(40,5))
            for axis in ['top','bottom','left','right']:  ax.spines[axis].set_linewidth(3)
            ax.tick_params(axis='both', which='major', labelsize=20,length=12.5,width=3,direction='in')
                
            for s in list( MAG['scaffold'].unique() ):


                m = MAG[MAG['scaffold']==s]

                ax.plot(m['genomic_coordinate'],m['SNV_density'],color=random_rgb(),linewidth=2)

                #ax.text(x=m['genomic_coordinate'].iloc[0]+100,y=0.0020,s=s,rotation=30)

            ax.set_ylabel('SNV Frequency',fontsize=20)
            ax.set_xlabel('Genomic Coordinate',fontsize=20)
            ax.scatter([],[],label=mag,alpha=0)
            ax.legend(loc='best',fontsize=30)
            ax.set_facecolor('#F2F2F2')
            fig.savefig(os.path.join(OUT_DIR,name,name+'_snv'), dpi=150,bbox_inches='tight')

        except NotADirectoryError: 

            pass

    return

def compareMAGs(MAGs,experiments,height=0,collapse=True):
    
    print()
    print(3*" <",f' Compare MAGs SNVs between experiments ','\n')

    if not os.path.isdir('./compareMAGs'):
        
        os.system(f"mkdir ./compareMAGs")
    
    
    for m in MAGs:
        
        name = m.replace('.fa','')

        print(" >",f' Compare MAG: {name}')
        
        multiMAGs={}
        
        for e in experiments:
            
            MAG_DIR = os.path.join('DIVA-MAGs',e,'MAGsSNVs',m.replace('.fa',''))

            for f in os.listdir(MAG_DIR):

                if '.csv' in f: multiMAGs[e]=pd.read_csv(os.path.join(MAG_DIR,f),sep='\t')
    
        sharedScf = set(multiMAGs[ experiments[0] ]['scaffold'])

        for e in experiments[1:]: sharedScf = sharedScf.intersection( multiMAGs[ e ]['scaffold'] )

        E,S = len(experiments), len(sharedScf)

        fig,ax=plt.subplots(figsize=(40,10))
        for axis in ['top','bottom','left','right']:  ax.spines[axis].set_linewidth(3)
        ax.tick_params(axis='both', which='major', labelsize=20,length=12.5,width=3,direction='in')
        
        ax.xaxis.set_major_locator(plt.MaxNLocator(40))
        ax.yaxis.set_major_locator(plt.MaxNLocator(10))

        ax.set_facecolor('#F2F2F2')
        ax.set_facecolor('#F2F2F2')
        ax.set_ylabel('SNV frequency',fontsize=30)
        ax.set_xlabel('Genome Position (bp)',fontsize=30)

        fig.tight_layout()  
        
        e_color = {e:random_rgb() for e in experiments}
        
        for e,j in zip(experiments,range(E)):

            C = multiMAGs[e]
            
            for scf,i in zip(sharedScf,range(S)):
            
                V = C[C['scaffold']==scf]
                ax.plot(V['genomic_coordinate'],V['SNV_density']+height*j,color=e_color[e],linewidth=2,zorder=E+1-j) 
         
        vy = C['SNV_density'].max()
        
        for scf,i in zip(sharedScf,range(S)):    
            
            V = C[C['scaffold']==scf]
            v0 = V['genomic_coordinate'].iloc[0]
            vy = C['SNV_density'].max()+height*j
            ax.plot(v0*np.ones(50),np.linspace(-5,2*vy+j*height),color='black',linewidth=1,alpha=1,zorder=0) 
        
        ax.set_ylim(-5,1.2*vy+j*height)
        ax.tick_params(axis='x', labelsize=25,rotation=20)
        ax.tick_params(axis='y', labelsize=25)
        
        ax.scatter([],[],label=m,alpha=0)
        ax.legend(loc='best',fontsize=40)
        
        
        fig.savefig(os.path.join('./compareMAGs',f'compare_{name}'), dpi=150,bbox_inches='tight')

"""
def compareMAGs(MAGs,experiments,collapse=True):
    
    if not os.path.isdir('./compareMAGs'):
        
        os.system(f"mkdir ./compareMAGs")
    
    
    for m in MAGs:
        
        multiMAGs={}
        
        for e in experiments:
            
            MAG_DIR = os.path.join('DIVA-MAGs',e,'MAGsSNVs',m.replace('.fa',''))

            for f in os.listdir(MAG_DIR):

                if '.csv' in f: multiMAGs[e]=pd.read_csv(os.path.join(MAG_DIR,f),sep='\t')
    
        sharedScf = set(multiMAGs[ experiments[0] ]['scaffold'])

        for e in experiments[1:]: sharedScf = sharedScf.intersection( multiMAGs[ e ]['scaffold'] )

        E,S = len(experiments), len(sharedScf)

        print(S)
        
        
        if collapse == True:
            
            label='collapsed'

            fig,ax=plt.subplots(ncols=1,nrows=S,figsize=(7*E,2*S))

            for scf,i in zip(sharedScf,range(S)):

                color = random_rgb()

                for e,j in zip(experiments,range(E)):

                    C = multiMAGs[e]
                    V = C[C['scaffold']==scf]
                    ax[i].plot(V['scaffold_coordinate'],V['SNV_density'],color=random_rgb(),alpha=0.7,linewidth=1) 
                    ax[i].set_facecolor('#F2F2F2')

                ax[i].set_ylabel(scf+'\n'+' SNV frequency',fontsize=7.5)

        else:
                
            label='extended'

            fig,ax=plt.subplots(ncols=E,nrows=S,figsize=(7*E,5*S))

            for scf,i in zip(sharedScf,range(S)):

                color = random_rgb()

                for e,j in zip(experiments,range(E)):

                    C = multiMAGs[e]
                    V = C[C['scaffold']==scf]
                    ax[i][j].plot(V['scaffold_coordinate'],V['SNV_density'],color=color,linewidth=1) 
                    ax[i][j].set_facecolor('#F2F2F2')

                ax[i][0].set_ylabel(scf+'\n'+' SNV frequency',fontsize=7.5)
        
        fig.tight_layout()            
        fig.savefig(os.path.join('./compareMAGs'+f'compare_{m}_{label}'), dpi=150,bbox_inches='tight')
        
        break
"""