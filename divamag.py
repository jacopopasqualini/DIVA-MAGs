import divaMAGF as dm
import os
import pandas as pd
import json

print()
print(14*" >",13*" <")
print(5*" >"+35*" "+5*" <",'')
print(5*" >"+" DivaMag: a microdiversity explorer"+5*" <")
print(5*" >"+35*" "+5*" <",'')
print(14*" >",13*" <",'\n')


with open('config.json') as json_file:
    config = json.load(json_file)
    config['compare_experiments']=bool(config['compare_experiments'])
    
# setting environment
DATA_DIR = config['DATA_DIR']

# user defined variables
experiments = config['experiments']#['BP2','BP1']
pc_threshold = config['pc_threshold'] # position coverage
window_size = config['window_size'] # width of windows in which you are looking at snvs density
step_size= config['step_size']
compare_experiments= config['compare_experiments']

#general purpose data
s2b_file = pd.read_csv( os.path.join(DATA_DIR,'scaffold_to_bin.stb'),sep='\t',header=None)
s2b_file = s2b_file.rename(columns={0:'scaffold',1:'MAG'})
s2b_file = s2b_file.set_index('scaffold')['MAG']
s2b_rule = s2b_file.to_dict()


# load data
MultIScaffold = {}
MultISnv = {}

for BP in experiments:
    
    # load the single experiment
    # IS stands for InStrain
    ISscaffold, ISsnvs = dm.load_experiment(experiment=BP,DATA_DIR=DATA_DIR)
    
    if compare_experiments == True:

        MultIScaffold[BP]=ISscaffold
        MultISnv[BP]=ISscaffold

    # create an output directory inside each one relative to an experiment
    
    RESULTS_DIR = os.path.join(config['DATA_DIR'],BP,'MAGsSNVs')
    
    if not os.path.isdir(RESULTS_DIR): 
        os.system(f"mkdir {RESULTS_DIR}")

    # associate mag to scaffolds and clean snvs according to pc_threshold
    scaffold2bin = ISscaffold[['scaffold','length']]
    scaffold2bin = scaffold2bin.set_index('scaffold')
    scaffold2bin['MAG']=pd.Series(s2b_rule)
    
    mags=dm.MAG2SNV(which_MAGs=config['MAGs'],
                    scaffold2bin=scaffold2bin,
                    scaffold=ISscaffold,
                    snv=ISsnvs,
                    window_size=window_size,
                    step_size=step_size,
                    pc_threshold=pc_threshold,
                    OUT_DIR=RESULTS_DIR)

    dm.SNVIsual(which_MAGs=config['MAGs'],
                OUT_DIR=RESULTS_DIR)

if config['compare_experiments']==True:    
    
    dm.compareMAGs(MAGs=config['MAGs'],experiments=config['experiments'],height=config['comparison_height'])
