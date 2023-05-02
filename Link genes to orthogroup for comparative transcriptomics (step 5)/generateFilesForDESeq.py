#python ./scripts/generateFilesForDESeq.py -i /biodata/dep_psl/grp_hacquard/Nathan/Bioinfo_Fantin/filtered_quant_intrastrainDupKept -o /biodata/dep_psl/grp_hacquard/Nathan/Bioinfo_Fantin/testdir -g2k /biodata/dep_psl/grp_hacquard/Nathan/Bioinfo_Fantin/all_anno_strains_KO.txt

import argparse
import sys
import pandas as pd
import os
import re

def get_params(argv):
	parser = argparse.ArgumentParser(description='To generate files from filtered Salmon outputs')
	parser.add_argument('-i', '--i', help="Folder containing all Salmon outputs", required=True)
	parser.add_argument('-ext','--ext', help="Extension of files to parse in input folder", default='.sf')
	parser.add_argument('-og2g', '--og2g', help="OG to genes file", required=True)
	parser.add_argument('-o', '--o', help="Outputdir", required=True)
	parser.add_argument('-OG', '--OG', help="Generate OG-aggregated files (default='yes')", default='yes')
	parser.add_argument('-OGtype','--OGtype', help='OG table from OrthoFinder ("of") or OG table generated from KOs ("ko")',default='of')
	a = parser.parse_args()
	return a


if __name__ == '__main__':
	a = get_params(sys.argv[1:])
	#genetoKO=pd.read_csv(a.g2k,sep='\t')[['gene_ID','KO_name']]
	#genetoKO=genetoKO.rename(index=str, columns={"gene_ID": "Name", "KO_name": "KO"})
	#OGtoGene=pd.read_csv(a.og2g, sep='\t', dtype=str)
	DF=pd.DataFrame()
	for f in [f for f in os.listdir(a.i) if f.endswith(a.ext)]:
		print('Opening '+f)
		df=pd.read_csv(a.i+'/'+f, sep='\t')[['Name','NumReads']]
		df=df.rename(index=str, columns={"NumReads":"NumReads_"+f.replace('quant_','').replace(a.ext,'')})
		if len(DF)==0:
			#print('len DF==0 when ',f)
			DF=df[[c for c in df.columns]]
		else:
			DF=DF.merge(df, on='Name')
	DF['Strain']=DF['Name'].str.split('|').str[0]
	

	#g2og=[]
	#with open(a.og2g,'r') as f:
	#	n=0
	#	for line in f:
	#		for gene in re.split(', |\t',line.replace('\r\n','').replace('\n',''))[1:]:
	#			if gene!='' and n>0:
	#				g2og.append({'Name':gene,'OG':line.split('\t')[0]})
	#		n+=1
	#print(len(DF[DF['Strain'].str.contains('Morel')]))
	#genetoOG=pd.DataFrame(g2og)
	if a.OGtype=='of':
		g2og=pd.read_csv(a.og2g, sep='\t').set_index('Unnamed: 0').T
		l=[]
		for c in g2og.columns:
			new=pd.DataFrame()
			new['Name']=', '.join(list(g2og[c].dropna())).split(', ')
			new['OG']=c
			l.append(new)
		new=pd.read_csv('/'.join(a.og2g.split('/')[:-1])+'/Orthogroups_UnassignedGenes.csv',sep='\t',dtype=str).set_index('Unnamed: 0')
		new['Name']= new[list(new.columns)].apply(lambda row: ''.join(row.values.astype(str)), axis=1).str.replace('nan','')
		new=new.reset_index(drop=False).rename(index=str, columns={'Unnamed: 0':'OG'})[['Name','OG']]
		l.append(new)
		genetoOG=pd.concat(l)
		DF=genetoOG.merge(DF, on='Name', how='left')
		print(DF)
	elif a.OGtype=='ko':
		g2og=pd.read_csv(a.og2g, sep='\t').set_index('Unnamed: 0').T
		l=[]
		for c in g2og.columns:
			new=pd.DataFrame()
			new['Name']=','.join(list(g2og[c].dropna())).split(',')
			new['OG']=c
			l.append(new)
		genetoOG=pd.concat(l)
		DF=genetoOG.merge(DF, on='Name', how='left')
		print(DF)

	else:
		print('ERROR in --OGtype value')
	
	DF=DF.dropna()
	DF[['Name','Strain','OG']+sorted([c for c in DF.columns if 'NumReads' in c])].to_csv(a.o+'/strain.genes.qf', sep='\t', index=False)
	print('strains.genes.qf OK')
	if a.OG=='yes':
		allOGs=list(set(DF.OG))
		df=[]
		dfst=[]
		dfog={}
		with open(a.o+'/strainshavingOG.txt','w+') as txt:
			for og in allOGs:
				DFog=DF[DF['OG']==og]
				dic=dict(DFog.sum(axis=0))
				dic['OG']=og
				df.append(dic)
				strains=sorted(list(set(DFog['Strain'])))
				txt.write(str(og)+'\t'+','.join(strains)+'\n')
				dfog[og]={}
				for st in strains:
					dftmp=DFog[DFog['Strain']==st]
					dicst=dict(dftmp.sum(axis=0))
					dicst['Strain']=st
					dicst['OG']=og
					dfst.append(dicst)
					dfog[og][st]=len(dftmp)
		print('strainshavingOG.txt OK')
		DFogcom=pd.DataFrame(df)
		DFogcom[['OG']+sorted([c for c in DF.columns if 'NumReads' in c])].sort_values(by=['OG']).to_csv(a.o+'/community.OG.qf', sep='\t', index=False)
		print('community.og.qf OK')
		DFogst=pd.DataFrame(dfst)
		DFogst[['Strain','OG']+sorted([c for c in DF.columns if 'NumReads' in c])].sort_values(by=['Strain','OG']).to_csv(a.o+'/strain.OG.qf', sep='\t', index=False)
		print('strain.OG.qf OK')
		pd.DataFrame(dfog)[allOGs].fillna(0).T.to_csv(a.o+'/strain.OG.counts.csv', sep='\t')
		print('strain.OGcounts.csv OK')
