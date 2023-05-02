import argparse
import sys
import pandas as pd

def get_params(argv):
	parser = argparse.ArgumentParser(description='filter out duplicated genes')
	parser.add_argument('-dup', '--dup', help="Duplicates file", required=True)
	parser.add_argument('-i', '--i', help="input to filter", required=True)
	parser.add_argument('-o', '--o', help="output (filtering result)", required=True)
	parser.add_argument('-which', '--which', help="which duplicates to remove: inter-species, intra-species or both ? intra|inter|both", default='both')
	a = parser.parse_args()
	return a

if __name__ == '__main__':
	a = get_params(sys.argv[1:])
	## Parse duplicates list, create list of IDs to remove from file
	dup=pd.read_csv(a.dup, sep='\t')
	toDel=[]
	for g in set(dup['RetainedTxp']):
		l=[g]+list(dup[dup['RetainedTxp']==g]['DuplicateTxp'])
		l2=[ele.split('|')[0] for ele in l]
		if a.which=='both':
			toDel+=l
		elif a.which=='inter':
			if len(set(l2))>1:
				toDel+=l
			else:
				pass
		elif a.which=='intra':
			if len(set(l2))==1:
				toDel+=l
			else:
				pass
		else:
			print("ERROR: Wrong --which argument. Please select between 'inter','intra' and 'both'")

	print(str(len(toDel))+' lines will be deleted from input file')
	## Filter Salmon quantification files
	inp=pd.read_csv(a.i,sep='\t')
	inp[~inp['Name'].isin(toDel)].to_csv(a.o, index=False, sep='\t')
	print('DONE')
