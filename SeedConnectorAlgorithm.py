#!/usr/bin/env python

import networkx as nx
import requests
import json
import webbrowser

def open_drugstone_in_browser(nodes, edges, seed_num):
    node_objs = []
    edge_objs = []
    for node in nodes[:seed_num]:
        node_objs.append({'id': node, 'group': 'gene'})
    for node in nodes[seed_num:]:
        node_objs.append({'id': node, 'group': 'connector'})
    for (start, end) in edges:
        edge_objs.append({'from': start, 'to': end})
    payload = {
        'network': {'nodes': node_objs, 'edges': edge_objs}, 
        'groups': {
            'nodeGroups': {
                'connector': {
                  'groupName': "Connector",
                  'color': "#ADD8E6",
                  'shape': 'circle'
                    
                },
                'gene': {
                  'groupName': "Gene",
                  'color': "#AA4A44",
                  'shape': 'triangle',
		  		  'font':{'color': '#000000'}
                }
            }
        }
    }
    headers = {'Content-Type': 'application/json; charset=UTF-8'}
    response = requests.post('https://api.drugst.one/create_network', data=json.dumps(payload), headers=headers)
    
    url = f"https://drugst.one?id={response.json()}"
    webbrowser.open(url, new=0, autoraise=True)
    return


def main(show_drugstone=True):

	file1=open("GWAS_CAD_in_PPI201307.txt","r") #The seed gene list, the genes must be in the human interactome, i.e. covered by file2
	file2=open("PPI201307.txt","r")   #human protein-protein interactions
	file3=open("GWAS_CAD_in_PPI201307_SeedPred.txt","w")# seed genes and predicted seed connectors
	file4=open("GWAS_CAD_in_PPI201307_SeedPred_module.txt","w")
 	
	## reading seed genes
	subprot=[]
	for line in file1:
		tmp=line.split()
		for i in range(len(tmp)):
			subprot.append(tmp[i])
	file1.close()
	seedprot=subprot[:]
	seed_num=len(subprot)
	
	## reading human protein-protein interactions
	G=nx.Graph()	
	L1=0
	for line in file2:
		tmp=line.split()
		G.add_node(tmp[0])
		G.add_node(tmp[1])
		G.add_edge(tmp[0],tmp[1])	
	file2.close()
		
	
	##candidate seed connectors are the neighbors of the seed genes

	CandGene=[]	
	for i in range(len(subprot)):
		N=G.neighbors(subprot[i])
		NB=list(N)
		L=len(NB)
		for j in range(L):
			if(NB[j] not in CandGene and NB[j] not in subprot):
				CandGene.append(NB[j])			

	## ranking the candidates based on their roles in increasing the largest connected component induced by the seed genes
	while 1:
			
		SubNet=G.subgraph(subprot)
		Gcc = sorted(nx.connected_components(SubNet), key=len, reverse=True)
		H = G.subgraph(Gcc[0])	
		LCC=len(set(H.nodes()).intersection(seedprot))	
		RankCand=[0 for i in range(len(CandGene))]
 				
		for k in range(len(CandGene)):
			TmpGene=[]
			TmpGene=subprot[:]
			TmpGene.append(CandGene[k])
			TmpNet=G.subgraph(TmpGene)  			

			Gcc = sorted(nx.connected_components(TmpNet), key=len, reverse=True)
			H1 = G.subgraph(Gcc[0]) 
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
			L1=len(list(G.neighbors(CandGene[idx[t]])))
			ratio=len(set(subprot).intersection(G.neighbors(CandGene[idx[t]])))/float(L)
			if(r<ratio):
				r=ratio
				s=t
		subprot.append(CandGene[idx[s]])		#the candidate with highest rank is added into the seed gene list, and the process repeats
		print (CandGene[idx[s]])
	
	## output the seeds and the connectors
	for i in range(seed_num):
		file3.write("%s\t%s\n"%(subprot[i],'seed'))
	for i in range(seed_num,len(subprot)):
		file3.write("%s\t%s\n"%(subprot[i],'connector'))
	file3.close()	
	
	Module=G.subgraph(subprot)
	for i in range(len(list(Module.edges()))):
		print(list(Module.edges())[i])

	if show_drugstone:
		open_drugstone_in_browser(subprot, list(Module.edges()), seed_num)
	
if __name__ == "__main__":
    main()