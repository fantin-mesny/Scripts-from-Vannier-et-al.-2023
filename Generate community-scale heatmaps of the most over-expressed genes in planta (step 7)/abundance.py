import os
import sys
import argparse
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt


def get_params(argv):
	parser = argparse.ArgumentParser(description='Draw summary figure')
	parser.add_argument('-i', '--i', help="input file (strain.genes.qf)", required=True)
	parser.add_argument('-o', '--o', help="output file", required=True)
	a = parser.parse_args()
	return a


if __name__ == '__main__':
	a = get_params(sys.argv[1:])
	DF=pd.read_csv(a.i, sep='\t').set_index('Name')
	conditions={'soil':['A','B','C'],'roots':['D','E','Bpre']}
	#conditions={}
	#for c in DF.columns:
	#	if c not in ['Strain','OG']:
	#		try:
	#			if c.split('_')[1].split('.')[0][:-1] not in conditions:
	#				conditions[c.split('_')[1].split('.')[0][:-1]]=[]
	#			conditions[c.split('_')[1].split('.')[0][:-1]].append(c)
	#		except:
	#			print(c)

	for cond in conditions:
		df=DF[['Strain']+['NumReads_'+s+'.sf' for s in conditions[cond]]]
		samples={c:sum(df[c]) for c in df.columns if 'NumReads' in c}
		ab={}
		for i in set(df['Strain']):
			sub=df[df['Strain']==i]
			ab[i]={c.split('_')[1]:sum(sub[c])/samples[c] for c in samples}
		ab=pd.DataFrame(ab).T
		means=ab.mean(axis=1)
		stds=ab.std(axis=1)
		ab['RNA relative abundance']=means
		ab['Standard deviation']=stds
		ab=ab.reset_index(drop=False).rename(index=str, columns={'index':'Strain ID'}).sort_values(by=['RNA relative abundance'], ascending=False)
		ab[['Strain ID','RNA relative abundance','Standard deviation']].to_csv(a.o+'_'+cond+'.csv', sep='\t', index=False)
		plt.figure(figsize=(8,12))
		sns.barplot('RNA relative abundance','Strain ID',data=ab, xerr=ab['Standard deviation'])
		plt.savefig(a.o+'_'+cond+'.png', bbox_inches = 'tight')
		#print(ab)

