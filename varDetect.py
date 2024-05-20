import pandas as pd 
import sys
import os 

# Main, handle passed arguments
if __name__ =='__main__':
    
    # TODO: argument handling, paths are hardcoded right now  
    pile_up = "/Users/Anthony/Desktop"
    
    df = pd.read_csv(pile_up, delimiter='\t')
    print(df.head(10))
    
