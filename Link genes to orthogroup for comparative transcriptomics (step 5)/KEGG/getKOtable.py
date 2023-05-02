import os 
import pandas as pd

output=pd.DataFrame()
output2=pd.DataFrame()

for f in os.listdir('.'):
	if '.emapper.annotations' in f:
		df=pd.read_csv(f,sep='\t',comment='#',header=None)
		df=df.groupby(11)[0].apply(list).to_dict()
		for ko in df:
			if ko!='-':
				output.loc[ko,f.split('_')[0]]=','.join(df[ko])
				output2.loc[ko,f.split('_')[0]]=len(df[ko])


output.to_csv('KO.csv',sep='\t')
output2.to_csv('KO.GeneCount.csv',sep='\t')

kolists=[kolist for kolist in output.index if ',' in kolist]
#print(kolists)
for ko in kolists:
	for bac in output.T[ko].dropna().index:
		for gene in output.loc[ko,bac].split(','):
			for k in ko.split(','):
				if k not in output.index: 
					output.loc[k,bac]=gene
					output2.loc[k,bac]=1
				elif str(output.loc[k,bac])=='nan':
					output.loc[k,bac]=gene
					output2.loc[k,bac]=1
				elif '|' in str(output.loc[k,bac]):
					output.loc[k,bac]=output.loc[k,bac]+','+gene
					output2.loc[k,bac]=output2.loc[k,bac]+1
				else:
					print('ERROR', output.loc[k,bac])


output=output.drop(index=kolists)
output2=output2.drop(index=kolists)
output.to_csv('KO_ungrouped.csv',sep='\t')
output2=output2.fillna(0)
output2.to_csv('KO_ungrouped.GeneCount.csv',sep='\t')

