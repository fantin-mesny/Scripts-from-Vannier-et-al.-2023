import os
import sys
import argparse
import seaborn as sns
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import matplotlib.colors as colors
from matplotlib import gridspec
from Bio import Phylo

def get_params(argv):
	parser = argparse.ArgumentParser(description='Draw summary figure')
	parser.add_argument('-ps', '--ps', help="Directory containing all per-strain DESeq outputs", required=True)
	parser.add_argument('-phy', '--phy', help="Phylogenetic tree (.nwk)", required=True)
	parser.add_argument('-c', '--c', help="DESeq output for community", required=True)
	parser.add_argument('-ab', '--ab', help="Average abundance CSV file", required=True)
	parser.add_argument('-o', '--o', help="Output directory/basename for figures", required=True)
	parser.add_argument('-ann','--ann',help="Annotation table with OG numbers in column 'og' and names in column 'short_name'", required=False, default=None)
	a = parser.parse_args()
	return a

def get_y_positions(tree): 
	maxheight = tree.count_terminals() 
	heights = dict((tip, maxheight - i) for i, tip in enumerate(reversed(tree.get_terminals()))) 
	def calc_row(clade): 
		for subclade in clade: 
			if subclade not in heights: 
				calc_row(subclade)
		heights[clade] = (heights[clade.clades[0]] + heights[clade.clades[-1]]) / 2.0 
	if tree.root.clades: 
		calc_row(tree.root)
	nametopos={str(i):heights[i] for i in heights if str(i)!='Clade'}
	return nametopos

def get_x_positions(tree): 
	depths = tree.depths() 
	if not max(depths.values()): 
		depths = tree.depths(unit_branch_lengths=True)
	nametopos={str(i):depths[i] for i in depths if str(i)!='Clade'}
	return nametopos 


if __name__ == '__main__':
	a = get_params(sys.argv[1:])
	if a.ann:
		annot=pd.read_excel(a.ann, encoding='utf7')[['og','short_name']].dropna()
		annot['og']=annot['og'].str.replace('K','ko:K')
		annot=annot[annot['short_name']!='-']
		annot['og2']=annot['og'].str.replace('ko:','')+' ('+annot['short_name']+')'
		annot=annot.set_index('og')['og2'].to_dict()
		print(annot)
	com0=pd.read_csv(a.c, sep=',', decimal=".", dtype={'Unnamed: 0':str, 'log2FoldChange':float}).sort_values(by=['log2FoldChange'], ascending=False)
	com=com0[com0['padj']<0.05][['Unnamed: 0','log2FoldChange']]
	comP=com0[com0['padj']<0.05][['Unnamed: 0','log2FoldChange']]
	# Load per-OG per-strain DESeq outputs
	detectedStrains=[]
	for st in [s for s in os.listdir(a.ps) if 'deseq' in s and '.pdf' not in s]:
		dfst=pd.read_csv(a.ps+'/'+st)[['Unnamed: 0','log2FoldChange']].rename(index=str, columns={"log2FoldChange": st.split('OG')[0]})
		com=com.merge(dfst, on='Unnamed: 0', how='left')
		dfst=pd.read_csv(a.ps+'/'+st)[['Unnamed: 0','padj']].rename(index=str, columns={"padj": st.split('OG')[0]})
		comP=comP.merge(dfst, on='Unnamed: 0', how='left')
		detectedStrains.append(st.split('OG')[0])

	# Load phylogeny
	tree = Phylo.read(a.phy, 'newick')
	y=get_y_positions(tree) #might be a zero-problem here i.e. strain '0209' is called '209' in df 'com'
	x=get_x_positions(tree)


	# Load abundance data
	abundance=pd.read_csv(a.ab,sep='\t')
	abundance['Strain']=abundance['Strain ID'].map(y)
	abundance=abundance.sort_values(by=['Strain'], ascending=True)


	# Draw figure
	com=com.set_index('Unnamed: 0')
	comP=comP.set_index('Unnamed: 0')
	com.rename(index=str, columns={'Unnamed: 0':'OG'}).to_csv(a.o+'_dataframe.csv')
	for st in y:
		if st not in com.columns:
			com[st]=np.nan
			comP[st]=np.nan


	df=com.rename(index=str, columns=y)
	comP=comP.rename(index=str, columns=y)

	# Set pre-figure (necessary to get labels on top of figure)
	fig=plt.figure(figsize=(27,23))
	gs1 = gridspec.GridSpec(4, 3,width_ratios=[0.40,1,0.1], height_ratios=[1,20,0.5,2])
	gs1.update(wspace=0, hspace=0.05) # set the spacing between axes. 
	axbar = plt.subplot(gs1[10]) ##Colorbar
	ax0 = plt.subplot(gs1[1])
	ax1 = plt.subplot(gs1[3])
	ax2 = plt.subplot(gs1[4], sharey=ax1)
	ax3 = plt.subplot(gs1[5])
	ax4 = plt.subplot(gs1[7],sharex=ax1)
	axnamecom=plt.subplot(gs1[0], sharex=ax1)

	allLogs=df[['log2FoldChange']+sorted([c for c in df.columns if c in y.values()])]
	comP=comP[['log2FoldChange']+sorted([c for c in df.columns if c in y.values()])].drop(['log2FoldChange'], axis=1)
	counts=comP[comP<0.05]
	order=pd.DataFrame(counts.count(axis=1)).sort_values(by=[0], ascending=False).rename(index=str, columns={0:'Rank'})
	allLogs=allLogs.merge(order,how='left', left_index=True, right_index=True).sort_values(by=['Rank'], ascending=False)

	for r in [0]:
		logs=allLogs.drop(['Rank'], axis=1)[r:r+200]
		if a.ann:		
			logs=logs.rename(index=annot)
			print(logs)

		## Set figure
		fig=plt.figure(figsize=(54,19))
		gs1 = gridspec.GridSpec(4, 3,width_ratios=[0.40,2,0.1], height_ratios=[1,20,2,0.5])
		gs1.update(wspace=0, hspace=0.05) # set the spacing between axes. 

		## Set axes
		axbar = plt.subplot(gs1[10]) ##Colorbar
		ax0 = plt.subplot(gs1[1],sharex=ax2)
		ax1 = plt.subplot(gs1[3])
		ax2 = plt.subplot(gs1[4], sharey=ax1)
		ax3 = plt.subplot(gs1[5])
		ax4 = plt.subplot(gs1[7])
		axnamecom=plt.subplot(gs1[0], sharex=ax1)
		axnamecom2=plt.subplot(gs1[6], sharex=ax1)
		axnamecom3=plt.subplot(gs1[9], sharex=ax1)


		##Heatmap community
		ax0.xaxis.set_label_position('top')
		ax0.xaxis.tick_top()
		sns.heatmap(logs[['log2FoldChange']][r:r+200].T, cbar=False, vmin=-8, vmax=8, cmap='bwr', center=0, ax=ax0,xticklabels=True, yticklabels=False)
		ax0.set_xlabel('Top 100 OGs')

		logs=logs.drop(['log2FoldChange'], axis=1).T

		##Phylogenetic tree
		ax1.axis('off')
		Phylo.draw(tree, do_show=False,axes=ax1,show_confidence=False,label_colors={X:'white' for X in x})
		n=0
		maxX=ax1.get_xlim()[1]
		for X in x:
			n+=1
			ax1.axhline(n, xmin=x[X]/maxX+0.03, color='black',linestyle=':',linewidth=0.75, zorder=10000)
			ax1.text(maxX-0.005, n, X, horizontalalignment='right',verticalalignment='bottom')
		

		axnamecom.text(x=maxX-0.01,y=0.45,ha='right',s='Community log2FoldChange')
		axnamecom.axis('off')

		##Heatmap
		#print(df[sorted([c for c in df.columns if c!='log2FoldChange' and c in y.values()])][r:r+100].T)
		g=sns.heatmap(logs,cbar_kws={"orientation": "horizontal"}, cmap='bwr',cbar_ax=axbar, vmin=-8, vmax=8, center=0, ax=ax2,xticklabels=False, yticklabels=False)
		g.set_facecolor('#cccccc')
		ax2.set_xlabel('')
		ax2.tick_params(
			axis='x',          # changes apply to the x-axis
			which='both',      # both major and minor ticks are affected
			bottom=False,      # ticks along the bottom edge are off
			top=False,         # ticks along the top edge are off
			labelbottom=False) # labels along the bottom edge are off
		axnamecom3.text(x=maxX-0.01,y=0.45, ha='right',s='log2FoldChanges')
		axnamecom3.axis('off')

		## Abundance barplot
		abundance=abundance[abundance['Strain ID'].isin(y)]
		sns.barplot(y='Strain',x='RNA relative abundance', data=abundance, orient='h', ax=ax3, color='#66cc92', xerr=abundance['Standard deviation'])
		ax3.axis('off')
		ax3.set_xlim([-0.01, max(abundance['RNA relative abundance'])+max(abundance['Standard deviation'])])

		## Counts of over-expressing strains
		allLogs['Rank']=-allLogs['Rank']
		allLogs=allLogs.reset_index(drop=False)[r:r+200]
		sns.barplot(x='Unnamed: 0',y='Rank', data=allLogs, orient='v', ax=ax4, color='#CD5C5C')
		ax4.axis('off')
		axnamecom2.text(x=maxX-0.01,y=0.45, ha='right', s='Number of strains overexpressing OG')
		axnamecom2.axis('off')
		
		plt.savefig(a.o+'_'+str(r+200)+'_nbSpeciesOverExpressing.pdf', dpi=200)
	

	allLogs=df[['log2FoldChange']+sorted([c for c in df.columns if c in y.values()])]
	order=pd.DataFrame(-allLogs.drop(['log2FoldChange'], axis=1).sum(axis=1)).sort_values(by=[0], ascending=True).rename(index=str, columns={0:'CumulativeFC'})
	allLogs=allLogs.merge(order,how='left', left_index=True, right_index=True).sort_values(by=['CumulativeFC'], ascending=True)
	allLogs.to_csv(a.o+'_dataframe_cumulativeFCsorted.csv')

	for r in [0]:
		logs=allLogs.drop(['CumulativeFC'], axis=1)[r:r+200]
		if a.ann:
			logs=logs.rename(index=annot)

		## Set figure
		fig=plt.figure(figsize=(54,19))
		gs1 = gridspec.GridSpec(4, 3,width_ratios=[0.40,2,0.1], height_ratios=[1,20,2,0.5])
		gs1.update(wspace=0, hspace=0.05) # set the spacing between axes. 

		## Set axes
		axbar = plt.subplot(gs1[10]) ##Colorbar
		ax0 = plt.subplot(gs1[1],sharex=ax2)
		ax1 = plt.subplot(gs1[3])
		ax2 = plt.subplot(gs1[4], sharey=ax1)
		ax3 = plt.subplot(gs1[5])
		ax4 = plt.subplot(gs1[7])
		axnamecom=plt.subplot(gs1[0], sharex=ax1)
		axnamecom2=plt.subplot(gs1[6], sharex=ax1)
		axnamecom3=plt.subplot(gs1[9], sharex=ax1)


		##Heatmap community
		ax0.xaxis.set_label_position('top')
		ax0.xaxis.tick_top()
		sns.heatmap(logs[['log2FoldChange']][r:r+200].T, cbar=False, vmin=-8, vmax=8, cmap='bwr', center=0, ax=ax0,xticklabels=True, yticklabels=False)
		ax0.set_xlabel('Top 100 OGs')

		logs=logs.drop(['log2FoldChange'], axis=1).T

		##Phylogenetic tree
		ax1.axis('off')
		Phylo.draw(tree, do_show=False,axes=ax1,show_confidence=False,label_colors={X:'white' for X in x})
		n=0
		maxX=ax1.get_xlim()[1]
		for X in x:
			n+=1
			ax1.axhline(n, xmin=x[X]/maxX+0.03, color='black',linestyle=':',linewidth=0.75, zorder=10000)
			ax1.text(maxX-0.005, n, X, horizontalalignment='right',verticalalignment='bottom')
		axnamecom.text(x=maxX-0.01,y=0.45,ha='right',s='Community log2FoldChange')
		axnamecom.axis('off')

		##Heatmap
		#print(df[sorted([c for c in df.columns if c!='log2FoldChange' and c in y.values()])][r:r+100].T)
		g=sns.heatmap(logs,cbar_kws={"orientation": "horizontal"}, cmap='bwr',cbar_ax=axbar, vmin=-8, vmax=8, center=0, ax=ax2,xticklabels=False, yticklabels=False)
		g.set_facecolor('#cccccc')
		ax2.set_xlabel('')
		ax2.tick_params(
			axis='x',          # changes apply to the x-axis
			which='both',      # both major and minor ticks are affected
			bottom=False,      # ticks along the bottom edge are off
			top=False,         # ticks along the top edge are off
			labelbottom=False) # labels along the bottom edge are off
		axnamecom3.text(x=maxX-0.01,y=0.45, ha='right',s='log2FoldChanges')
		axnamecom3.axis('off')

		## Abundance barplot
		abundance=abundance[abundance['Strain ID'].isin(y)]
		sns.barplot(y='Strain',x='RNA relative abundance', data=abundance, orient='h', ax=ax3, color='#66cc92', xerr=abundance['Standard deviation'])
		ax3.axis('off')
		ax3.set_xlim([-0.01, max(abundance['RNA relative abundance'])+max(abundance['Standard deviation'])])

		## Cumulative LFC
		counts=pd.DataFrame(-logs.sum(axis='rows')).reset_index(drop=False)
		sns.barplot(x='Unnamed: 0',y=0, data=counts, orient='v', ax=ax4, color='#CD5C5C')
		ax4.axis('off')
		axnamecom2.text(x=maxX-0.01,y=0.45, ha='right', s='Cumulated Log2FC')
		axnamecom2.axis('off')
		
		plt.savefig(a.o+'_'+str(r+200)+'cumulativeFC_topOE.pdf', dpi=200)

	allLogs=df[['log2FoldChange']+sorted([c for c in df.columns if c in y.values()])]
	order=pd.DataFrame(allLogs.drop(['log2FoldChange'], axis=1).sum(axis=1)).sort_values(by=[0], ascending=True).rename(index=str, columns={0:'CumulativeFC'})
	allLogs=allLogs.merge(order,how='left', left_index=True, right_index=True).sort_values(by=['CumulativeFC'], ascending=True)

	for r in [0]:
		logs=allLogs.drop(['CumulativeFC'], axis=1)[r:r+200]
		print(logs.reset_index()[99:101])

		## Set figure
		fig=plt.figure(figsize=(54,19))
		gs1 = gridspec.GridSpec(4, 3,width_ratios=[0.40,2,0.1], height_ratios=[1,20,2,0.5])
		gs1.update(wspace=0, hspace=0.05) # set the spacing between axes. 

		## Set axes
		axbar = plt.subplot(gs1[10]) ##Colorbar
		ax0 = plt.subplot(gs1[1],sharex=ax2)
		ax1 = plt.subplot(gs1[3])
		ax2 = plt.subplot(gs1[4], sharey=ax1)
		ax3 = plt.subplot(gs1[5])
		ax4 = plt.subplot(gs1[7])
		axnamecom=plt.subplot(gs1[0], sharex=ax1)
		axnamecom2=plt.subplot(gs1[6], sharex=ax1)
		axnamecom3=plt.subplot(gs1[9], sharex=ax1)


		##Heatmap community
		ax0.xaxis.set_label_position('top')
		ax0.xaxis.tick_top()
		sns.heatmap(logs[['log2FoldChange']][r:r+200].T, cbar=False, vmin=-8, vmax=8, cmap='bwr', center=0, ax=ax0,xticklabels=True, yticklabels=False)
		ax0.set_xlabel('Top 200 OGs')

		logs=logs.drop(['log2FoldChange'], axis=1).T

		##Phylogenetic tree
		ax1.axis('off')
		Phylo.draw(tree, do_show=False,axes=ax1,show_confidence=False,label_colors={X:'white' for X in x})
		n=0
		maxX=ax1.get_xlim()[1]
		for X in x:
			n+=1
			ax1.axhline(n, xmin=x[X]/maxX+0.03, color='black',linestyle=':',linewidth=0.75, zorder=10000)
			ax1.text(maxX-0.005, n, X, horizontalalignment='right',verticalalignment='bottom')
		axnamecom.text(x=maxX-0.01,y=0.45,ha='right',s='Community log2FoldChange')
		axnamecom.axis('off')

		##Heatmap
		#print(df[sorted([c for c in df.columns if c!='log2FoldChange' and c in y.values()])][r:r+100].T)
		g=sns.heatmap(logs,cbar_kws={"orientation": "horizontal"}, cmap='bwr',cbar_ax=axbar, vmin=-8, vmax=8, center=0, ax=ax2,xticklabels=False, yticklabels=False)
		g.set_facecolor('#cccccc')
		ax2.set_xlabel('')
		ax2.tick_params(
			axis='x',          # changes apply to the x-axis
			which='both',      # both major and minor ticks are affected
			bottom=False,      # ticks along the bottom edge are off
			top=False,         # ticks along the top edge are off
			labelbottom=False) # labels along the bottom edge are off
		axnamecom3.text(x=maxX-0.01,y=0.45, ha='right',s='log2FoldChanges')
		axnamecom3.axis('off')

		## Abundance barplot
		abundance=abundance[abundance['Strain ID'].isin(y)]
		sns.barplot(y='Strain',x='RNA relative abundance', data=abundance, orient='h', ax=ax3, color='#66cc92', xerr=abundance['Standard deviation'])
		ax3.axis('off')
		ax3.set_xlim([-0.01, max(abundance['RNA relative abundance'])+max(abundance['Standard deviation'])])

		## Counts of over-expressing strains
		counts=pd.DataFrame(logs.sum(axis='rows')).reset_index(drop=False)
		
		sns.barplot(x='Unnamed: 0',y=0, data=counts, orient='v', ax=ax4, color='#0066cc')
		ax4.axis('off')
		axnamecom2.text(x=maxX-0.01,y=0.45, ha='right', s='Cumulated Log2FC')
		axnamecom2.axis('off')
		
		plt.savefig(a.o+'_'+str(r+200)+'cumulativeFC_topUE.pdf', dpi=200)



	allLogs=df[['log2FoldChange']+sorted([c for c in df.columns if c in y.values()])]
	order=pd.DataFrame(-allLogs.drop(['log2FoldChange'], axis=1).mean(axis=1)).sort_values(by=[0], ascending=True).rename(index=str, columns={0:'MeanFC'})
	allLogs=allLogs.merge(order,how='left', left_index=True, right_index=True).sort_values(by=['MeanFC'], ascending=True)

	for r in [0]:
		logs=allLogs.drop(['MeanFC'], axis=1)[r:r+200]

		## Set figure
		fig=plt.figure(figsize=(54,19))
		gs1 = gridspec.GridSpec(4, 3,width_ratios=[0.40,2,0.1], height_ratios=[1,20,2,0.5])
		gs1.update(wspace=0, hspace=0.05) # set the spacing between axes. 

		## Set axes
		axbar = plt.subplot(gs1[10]) ##Colorbar
		ax0 = plt.subplot(gs1[1],sharex=ax2)
		ax1 = plt.subplot(gs1[3])
		ax2 = plt.subplot(gs1[4], sharey=ax1)
		ax3 = plt.subplot(gs1[5])
		ax4 = plt.subplot(gs1[7])
		axnamecom=plt.subplot(gs1[0], sharex=ax1)
		axnamecom2=plt.subplot(gs1[6], sharex=ax1)
		axnamecom3=plt.subplot(gs1[9], sharex=ax1)


		##Heatmap community
		ax0.xaxis.set_label_position('top')
		ax0.xaxis.tick_top()
		sns.heatmap(logs[['log2FoldChange']][r:r+200].T, cbar=False, vmin=-8, vmax=8, cmap='bwr', center=0, ax=ax0,xticklabels=True, yticklabels=False)
		ax0.set_xlabel('Top 100 OGs')

		logs=logs.drop(['log2FoldChange'], axis=1).T

		##Phylogenetic tree
		ax1.axis('off')
		Phylo.draw(tree, do_show=False,axes=ax1,show_confidence=False,label_colors={X:'white' for X in x})
		n=0
		maxX=ax1.get_xlim()[1]
		for X in x:
			n+=1
			ax1.axhline(n, xmin=x[X]/maxX+0.03, color='black',linestyle=':',linewidth=0.75, zorder=10000)
			ax1.text(maxX-0.005, n, X, horizontalalignment='right',verticalalignment='bottom')
		axnamecom.text(x=maxX-0.01,y=0.45,ha='right',s='Community log2FoldChange')
		axnamecom.axis('off')

		##Heatmap
		#print(df[sorted([c for c in df.columns if c!='log2FoldChange' and c in y.values()])][r:r+100].T)
		g=sns.heatmap(logs,cbar_kws={"orientation": "horizontal"}, cmap='bwr',cbar_ax=axbar, vmin=-8, vmax=8, center=0, ax=ax2,xticklabels=False, yticklabels=False)
		g.set_facecolor('#cccccc')
		ax2.set_xlabel('')
		ax2.tick_params(
			axis='x',          # changes apply to the x-axis
			which='both',      # both major and minor ticks are affected
			bottom=False,      # ticks along the bottom edge are off
			top=False,         # ticks along the top edge are off
			labelbottom=False) # labels along the bottom edge are off
		axnamecom3.text(x=maxX-0.01,y=0.45, ha='right',s='log2FoldChanges')
		axnamecom3.axis('off')

		## Abundance barplot
		abundance=abundance[abundance['Strain ID'].isin(y)]
		sns.barplot(y='Strain',x='RNA relative abundance', data=abundance, orient='h', ax=ax3, color='#66cc92', xerr=abundance['Standard deviation'])
		ax3.axis('off')
		ax3.set_xlim([-0.01, max(abundance['RNA relative abundance'])+max(abundance['Standard deviation'])])

		## Counts of over-expressing strains
		counts=pd.DataFrame(-logs.mean(axis='rows')).reset_index(drop=False)
		
		sns.barplot(x='Unnamed: 0',y=0, data=counts, orient='v', ax=ax4, color='#CD5C5C')
		ax4.axis('off')
		axnamecom2.text(x=maxX-0.01,y=0.45, ha='right', s='Average Log2FC')
		axnamecom2.axis('off')
		
		plt.savefig(a.o+'_'+str(r+200)+'meanFC.pdf', dpi=200)



	allLogs=df[['log2FoldChange']+sorted([c for c in df.columns if c in y.values()])]
	order=pd.DataFrame(allLogs.drop(['log2FoldChange'], axis=1).mean(axis=1)).sort_values(by=[0], ascending=True).rename(index=str, columns={0:'MeanFC'})
	allLogs=allLogs.merge(order,how='left', left_index=True, right_index=True).sort_values(by=['MeanFC'], ascending=True)

	for r in [0]:
		logs=allLogs.drop(['MeanFC'], axis=1)[r:r+200]

		## Set figure
		fig=plt.figure(figsize=(54,19))
		gs1 = gridspec.GridSpec(4, 3,width_ratios=[0.40,2,0.1], height_ratios=[1,20,2,0.5])
		gs1.update(wspace=0, hspace=0.05) # set the spacing between axes. 

		## Set axes
		axbar = plt.subplot(gs1[10]) ##Colorbar
		ax0 = plt.subplot(gs1[1],sharex=ax2)
		ax1 = plt.subplot(gs1[3])
		ax2 = plt.subplot(gs1[4], sharey=ax1)
		ax3 = plt.subplot(gs1[5])
		ax4 = plt.subplot(gs1[7])
		axnamecom=plt.subplot(gs1[0], sharex=ax1)
		axnamecom2=plt.subplot(gs1[6], sharex=ax1)
		axnamecom3=plt.subplot(gs1[9], sharex=ax1)


		##Heatmap community
		ax0.xaxis.set_label_position('top')
		ax0.xaxis.tick_top()
		sns.heatmap(logs[['log2FoldChange']][r:r+200].T, cbar=False, vmin=-8, vmax=8, cmap='bwr', center=0, ax=ax0,xticklabels=True, yticklabels=False)
		ax0.set_xlabel('Top 100 OGs')

		logs=logs.drop(['log2FoldChange'], axis=1).T

		##Phylogenetic tree
		ax1.axis('off')
		Phylo.draw(tree, do_show=False,axes=ax1,show_confidence=False,label_colors={X:'white' for X in x})
		n=0
		maxX=ax1.get_xlim()[1]
		for X in x:
			n+=1
			ax1.axhline(n, xmin=x[X]/maxX+0.03, color='black',linestyle=':',linewidth=0.75, zorder=10000)
			ax1.text(maxX-0.005, n, X, horizontalalignment='right',verticalalignment='bottom')
		axnamecom.text(x=maxX-0.01,y=0.45,ha='right',s='Community log2FoldChange')
		axnamecom.axis('off')

		##Heatmap
		#print(df[sorted([c for c in df.columns if c!='log2FoldChange' and c in y.values()])][r:r+100].T)
		g=sns.heatmap(logs,cbar_kws={"orientation": "horizontal"}, cmap='bwr',cbar_ax=axbar, vmin=-8, vmax=8, center=0, ax=ax2,xticklabels=False, yticklabels=False)
		g.set_facecolor('#cccccc')
		ax2.set_xlabel('')
		ax2.tick_params(
			axis='x',          # changes apply to the x-axis
			which='both',      # both major and minor ticks are affected
			bottom=False,      # ticks along the bottom edge are off
			top=False,         # ticks along the top edge are off
			labelbottom=False) # labels along the bottom edge are off
		axnamecom3.text(x=maxX-0.01,y=0.45, ha='right',s='log2FoldChanges')
		axnamecom3.axis('off')

		## Abundance barplot
		abundance=abundance[abundance['Strain ID'].isin(y)]
		sns.barplot(y='Strain',x='RNA relative abundance', data=abundance, orient='h', ax=ax3, color='#66cc92', xerr=abundance['Standard deviation'])
		ax3.axis('off')
		ax3.set_xlim([-0.01, max(abundance['RNA relative abundance'])+max(abundance['Standard deviation'])])

		## Counts of over-expressing strains
		counts=pd.DataFrame(logs.mean(axis='rows')).reset_index(drop=False)
		
		sns.barplot(x='Unnamed: 0',y=0, data=counts, orient='v', ax=ax4, color='#0066cc')
		ax4.axis('off')
		axnamecom2.text(x=maxX-0.01,y=0.45, ha='right', s='Average Log2FC')
		axnamecom2.axis('off')
		
		plt.savefig(a.o+'_'+str(r+200)+'meanFC_UE.pdf', dpi=200)

