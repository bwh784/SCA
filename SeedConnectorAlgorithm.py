#!/usr/bin/env python

import getopt, sys, urllib, time
import random as rd
import networkx as nx

def main():

	file0=open("PPI201307_proteins.txt","r")
	file1=open("GWAS_CAD_in_PPI201307.txt","r")
	file2=open("PPI201307.txt","r")
	file3=open("GWAS_CAD_in_PPI201307_SeedPred.txt","w")
	
 	G=nx.Graph()	
	subprot=[]
	for line in file1:
		tmp=line.split()
		for i in range(len(tmp)):
			subprot.append(tmp[i])
	file1.close()
	seedprot=subprot[:]
	
	prot=[]
	for line in file0:
		tmp=line.split()
		prot.append(tmp[0])
	file0.close()	

 	L1=0
	for line in file2:
		tmp=line.split()
		G.add_node(tmp[0])
		G.add_node(tmp[1])
		G.add_edge(tmp[0],tmp[1])	
	file2.close()
			
	
	CandGene=[]	
	for i in range(len(subprot)):
		N=G.neighbors(subprot[i])
		for j in range(len(N)):
			if(N[j] not in CandGene and N[j] not in subprot):
				CandGene.append(N[j])			

	while 1:
			
 		SubNet=G.subgraph(subprot)
 		H=sorted(nx.connected_component_subgraphs(SubNet),key = len, reverse=True)[0]		
 		LCC=len(set(H.nodes()).intersection(seedprot))	
 		RankCand=[0 for i in range(len(CandGene))]
 				
 		for k in range(len(CandGene)):
 			TmpGene=[]
 			TmpGene=subprot[:]
 			TmpGene.append(CandGene[k])
 			TmpNet=G.subgraph(TmpGene)  			

 			H1=sorted(nx.connected_component_subgraphs(TmpNet),key = len, reverse=True)[0]		
 			LCC1=len(set(H1.nodes()).intersection(seedprot))
  			RankCand[k]=LCC1-LCC
 			del TmpGene[:]
 		#print RankCand 
		idx=[s for s, x in enumerate(RankCand) if x == max(RankCand)]
		if(max(RankCand)==0):
			break
		
		s=-1
		r=-1
		for t in range(len(idx)):
			ratio=len(set(subprot).intersection(G.neighbors(CandGene[idx[t]])))/float(len(G.neighbors(CandGene[idx[t]])))
			if(r<ratio):
				r=ratio
				s=t			
		subprot.append(CandGene[idx[s]])
		print CandGene[idx[s]],RankCand[idx[s]],LCC

		
	
	for i in range(len(subprot)):
		file3.write("%s\n"%subprot[i])
	file3.close()	
	
	
if __name__ == "__main__":
    main()