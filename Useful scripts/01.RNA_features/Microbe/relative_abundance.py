import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import numpy as np
import pandas as pd
import argparse
import sys
parser = argparse.ArgumentParser(description='Calculate relative abundacnce')
parser.add_argument('--input', '-i', type=str, required=True, help='input count matrix')
parser.add_argument('--output','-o',type=str,required=True,help='output TPM matrix')
args = parser.parse_args()
print("Load data ...")
df = pd.read_csv(args.input,index_col=0,sep="\t")
print(df.columns)
print("Done .")
print("Calculate norm ...")
#length = df.index.map(lambda x:x.split("|")[-1]).astype('int')
lengthScaledDf = pd.DataFrame((df.values),index=df.index,columns=df.columns)
(100*lengthScaledDf.div(lengthScaledDf.sum(axis=0))).round(4).to_csv(args.output,sep="\t")
print("Done .")
