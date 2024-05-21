import pandas as pd 
import sys
import os 

# Main, handle passed arguments
if __name__ =='__main__':
    
    # TODO: argument handling, paths are hardcoded right now  
    pile_up = sys.argv[1]
    columns = ['Chromosome', 'Position', 'Reference_Base', 'Coverage', 'Bases', 'Qualities']
    df = pd.read_csv(pile_up, sep='\t', header=None, names=columns)
    print(df)
