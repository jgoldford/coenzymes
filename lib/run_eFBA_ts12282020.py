import cobra
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-m", "--model", type=str,
                    help="path to cobra sbml file")
parser.add_argument("-e", "--media", type=str,
                    help="path to media text file")
parser.add_argument("-o", "--output", type=str,
                    help="path to output results")


args = parser.parse_args()
file_name = args.model;
media_file = args.media;
outFile = args.output

model = cobra.io.read_sbml_model(file_name)
media = pd.read_csv(media_file,sep='\t')

# get the base medium
medium = model.medium
medium.pop('EX_glc__D_e')
medium.pop('EX_nh4_e')
medium.pop('EX_o2_e')
medium_base = medium.copy()


obj_values = []
for idx,row in media.iterrows():
    m = medium_base.copy();
    m[row.carbon] = 10;
    m[row.nitrogen] = 1000;
    if row.electron != 'None':
        m[row.electron] = 1000;
    
    model.medium = m;
    f = model.slim_optimize(error_value=0)
    obj_values.append(f)
    df_tmp = pd.DataFrame({'growth_rate':obj_values})
    df_tmp.to_csv(outFile + '.simple.temp.txt',sep='\t')
    
media['growth_rate'] = obj_values;
media.to_csv(outFile,sep='\t')