#!/usr/bin/python3
import scipy.stats as ss
import numpy as np
import collections
import warnings
import sys,os
warnings.filterwarnings("ignore")
np.warnings.filterwarnings('ignore')

if not sys.warnoptions:
    warnings.simplefilter("ignore")

from argparse import ArgumentParser
from signal import signal, SIGPIPE, SIG_DFL
from statsmodels.sandbox.stats.multicomp import multipletests
from scipy.stats import spearmanr

def getGeneA_Del_GeneB_Amp(geneA, geneB, CNADict,geneExprDictTCGA, gTEXGeneExprDict):
	dCNAHigh=collections.OrderedDict()
	dCNALow=collections.OrderedDict()
	pValGeneA_Amp=1
	pvalGeneA_Del=1
	listOfexprHigh=list()
	listOfexprLow=list()
	listOfexprGTEX_A=list()
	listOfexprGTEX_B=list()
	listOfexprHighSample=list()
	listOfexprLowSample=list()
	
	for sampleList in CNADict[geneA]:
		for e in CNADict[geneA][sampleList]:
			if float(e) < -1:  
				dCNALow[sampleList]=float(e)

	for sampleList in CNADict[geneB]:
		if sampleList in dCNALow:
			for e in CNADict[geneB][sampleList]:
				if float(e) > 1:  
					dCNAHigh[sampleList]=float(e)

	for smpl, allExprValue in geneExprDictTCGA[geneB].items():
		if smpl in dCNAHigh:
			listOfexprHigh.append(float(list(allExprValue)[0]))
			listOfexprHighSample.append(smpl)

	for smpl, allExprValue in geneExprDictTCGA[geneA].items():
		if smpl in dCNALow:
			listOfexprLow.append(float(list(allExprValue)[0]))
			listOfexprLowSample.append(smpl)

	for s2, exprVal in gTEXGeneExprDict[geneA].items():
			listOfexprGTEX_A.append(float(list(exprVal)[0]))

	for s2, exprVal in gTEXGeneExprDict[geneB].items():
			listOfexprGTEX_B.append(float(list(exprVal)[0]))

	resultA_Del=ss.ranksums(listOfexprLow,listOfexprGTEX_A)
	pValGeneA_Del = resultA_Del.pvalue

	resultB_Amp=ss.ranksums(listOfexprHigh,listOfexprGTEX_B)
	pvalGeneB_Amp = resultB_Amp.pvalue
	
	if pValGeneA_Del <= 1 and pvalGeneB_Amp <= 1:
		tt=ss.ranksums(listOfexprHigh,listOfexprLow)
		return tt.pvalue 
	
	return 1

def getGeneA_Amp_GeneB_Del(geneA, geneB, CNADict,geneExprDictTCGA, gTEXGeneExprDict):
	dCNAHigh=collections.OrderedDict()
	dCNALow=collections.OrderedDict()
	pValGeneA_Amp=1
	pvalGeneA_Del=1
	listOfexprHigh=list()
	listOfexprLow=list()
	listOfexprGTEX_A=list()
	listOfexprGTEX_B=list()
	listOfexprHighSample=list()
	listOfexprLowSample=list()
	
	for sampleList in CNADict[geneA]:
		for e in CNADict[geneA][sampleList]:
			if float(e)> 1:  
				dCNAHigh[sampleList]=float(e)

	for sampleList in CNADict[geneB]:
		if sampleList in dCNAHigh:
			for e in CNADict[geneB][sampleList]:
				if float(e) < -1:  
					dCNALow[sampleList]=float(e)

	for smpl, allExprValue in geneExprDictTCGA[geneA].items():
		if smpl in dCNAHigh:
			listOfexprHigh.append(float(list(allExprValue)[0]))
			listOfexprHighSample.append(smpl)

	for smpl, allExprValue in geneExprDictTCGA[geneB].items():
		if smpl in dCNALow:
			listOfexprLow.append(float(list(allExprValue)[0]))
			listOfexprLowSample.append(smpl)

	for s2, exprVal in gTEXGeneExprDict[geneA].items():
			listOfexprGTEX_A.append(float(list(exprVal)[0]))

	for s2, exprVal in gTEXGeneExprDict[geneB].items():
			listOfexprGTEX_B.append(float(list(exprVal)[0]))

	resultA_Amp=ss.ranksums(listOfexprHigh,listOfexprGTEX_A)
	pValGeneA_Amp = resultA_Amp.pvalue

	resultB_Del=ss.ranksums(listOfexprLow,listOfexprGTEX_B)
	pvalGeneB_Del = resultB_Del.pvalue
	
	if pValGeneA_Amp <= 1  and pvalGeneB_Del <= 1:
		tt=ss.ranksums(listOfexprHigh,listOfexprLow)
		return tt.pvalue
	
	return 1

def getGeneB_Del(gene,CNASampleList, CNADict,geneExprDictTCGA, gTEXGeneExprDict):
	dCNALow=collections.OrderedDict()
	listOfexprLow=list()
	listOfexprHigh=list()
	listOfexprGTEX=list()
	pvalGeneB_Del=1

	for sampleList in CNASampleList:
		for e in CNADict[gene][sampleList]:
			if float(e) < -1:
				dCNALow[sampleList]=float(e)

	for smpl, allExprValue in geneExprDictTCGA[gene].items():
		if smpl in dCNALow:
			listOfexprLow.append(float(list(allExprValue)[0]))
			
		if smpl in CNASampleList:
			listOfexprHigh.append(float(list(allExprValue)[0]))

	for s2, exprVal in gTEXGeneExprDict[gene].items():
			listOfexprGTEX.append(float(list(exprVal)[0]))
	
	resultB_Del=ss.ranksums(listOfexprLow,listOfexprGTEX)
	pvalGeneB_Del = resultB_Del.pvalue

	return pvalGeneB_Del, len(dCNALow)
		
def getGeneB_Amp(gene,CNASampleList, CNADict,geneExprDictTCGA,gTEXGeneExprDict):
	dCNAHigh=collections.OrderedDict()
	listOfexprHigh=list()
	listOfexprLow=list()
	listOfexprGTEX=list()
	pvalGeneB_Amp=1

	for sampleList in CNASampleList:
		for e in CNADict[gene][sampleList]:
			if float(e) > 1:
				dCNAHigh[sampleList]=float(e)

	for smpl, allExprValue in geneExprDictTCGA[gene].items():
		if smpl in dCNAHigh:
			listOfexprHigh.append(float(list(allExprValue)[0]))
			
		if smpl in CNASampleList:
			listOfexprLow.append(float(list(allExprValue)[0]))

	for s2, exprVal in gTEXGeneExprDict[gene].items():
			listOfexprGTEX.append(float(list(exprVal)[0]))
	
	resultB_Amp=ss.ranksums(listOfexprHigh,listOfexprGTEX)
	pvalGeneB_Amp = resultB_Amp.pvalue

	return pvalGeneB_Amp,len(dCNAHigh)

def getGene_Amplification_Del(gene, CNADict, geneExprDictTCGA, gTEXGeneExprDict):
	dCNAHigh=collections.OrderedDict()
	dCNALow=collections.OrderedDict()
	pValGeneA_Amp=1
	pvalGeneA_Del=1
	listOfexprHigh=list()
	listOfexprLow=list()
	listOfexprGTEX=list()

	d = {}
	if gene not in d:
		d[gene] = {}
	
	#get sample in CNA high and low separtely
	for sampleList in CNADict[gene]:
		for e in CNADict[gene][sampleList]:
			if float(e) > 1:
				dCNAHigh[sampleList]=float(e)
			elif float(e) < -1:
				dCNALow[sampleList]=float(e)

	for smpl, allExprValue in geneExprDictTCGA[gene].items():
		if smpl in dCNAHigh:
			listOfexprHigh.append(float(list(allExprValue)[0]))
			
		if smpl in dCNALow:
			listOfexprLow.append(float(list(allExprValue)[0]))
			
	for s2, exprVal in gTEXGeneExprDict[gene].items():
			listOfexprGTEX.append(float(list(exprVal)[0]))

	resultA_Amp=ss.ranksums(listOfexprHigh,listOfexprGTEX)
	pValGeneA_Amp = resultA_Amp.pvalue

	resultA_Del=ss.ranksums(listOfexprLow,listOfexprGTEX)
	pvalGeneA_Del = resultA_Del.pvalue

	return pValGeneA_Amp, pvalGeneA_Del, dCNAHigh, dCNALow


def mkArgParser():
  parser = ArgumentParser()

  parser.add_argument("bc_all_gene_entrezid", help="input i.e bc_all_gene_entrezid file")
  parser.add_argument("linearCNA", help="input i.e data_CNA.txt normalized for certain cancer type file")
  parser.add_argument("geneExpression", help="input i.e data_RNA_Seq_v2_expression_median.txt for certain cancer type file")
  parser.add_argument("GTEXgeneExpression", help="input i.e GTEX breast_rnaseq_genexpr.txt for certain cancer type file")

  return parser

if __name__ == '__main__':
	signal(SIGPIPE,SIG_DFL)

	args = mkArgParser().parse_args()

	samples = collections.OrderedDict()
	linearCNADict = collections.OrderedDict()
	f=open(args.linearCNA, "r")
	line = f.readline().strip()
	line_ = line.split('\t')
	line_ = line_[2:]
	for smpl in line_:
		msmpl=smpl[0:15]
		samples[msmpl]=msmpl

	with open(args.linearCNA, "r") as f:
		next(f)
		for line in f:
			line_ = line.split('\t')
			geneName = line_[0].strip()
			entrezID = line_[1].strip()
			line_ = line_[2:]
			c=0
			if not entrezID in linearCNADict:
				linearCNADict[entrezID] = {}
				for smpl in samples:
					if not smpl in linearCNADict[entrezID]:
						linearCNADict[entrezID][smpl] = list()
					line_[c]=line_[c].strip()
					if line_[c]=="NA" or line_[c]=="None" :
						line_[c]=0

					linearCNADict[entrezID][smpl].append(float(line_[c]))
					c = c + 1
	
	samplesBC = collections.OrderedDict()
	geneExprDictBC =collections.OrderedDict()
	f=open(args.geneExpression, "r")
	line = f.readline().strip()
	line_ = line.split('\t')
	line_ = line_[2:]
	for smpl in line_:
		msmpl=smpl[0:15] # TCGA sample ID (the first 15 character)
		samplesBC[msmpl]=msmpl
	
	with open(args.geneExpression, "r") as f:
		next(f)
		for line in f:
			line_ = line.split('\t')
			geneName = line_[0].strip()
			entrezID = line_[1].strip()
			line_ = line_[2:]
			c=0

			if not entrezID in geneExprDictBC:
				geneExprDictBC[entrezID] = {}
			for smpl in samplesBC:
				if not smpl in geneExprDictBC[entrezID]:
					geneExprDictBC[entrezID][smpl] = list()
				line_[c]=line_[c].strip()
				if line_[c]=="NA" or line_[c]=="None" :
					line_[c]=0
				geneExprDictBC[entrezID][smpl].append(float(line_[c]))
				c = c + 1


	samplesGTEX = collections.OrderedDict()
	geneExprDictGTEX =collections.OrderedDict()
	f=open(args.GTEXgeneExpression, "r")
	line = f.readline().strip()
	line_ = line.split('\t')
	line_ = line_[2:]
	for smpl in line_:
		samplesGTEX[smpl]=smpl
	
	with open(args.GTEXgeneExpression, "r") as f:
		next(f)
		for line in f:
			line_ = line.split('\t')
			entrezID = line_[1].strip()
			line_ = line_[2:]
			c=0

			if not entrezID in geneExprDictGTEX:
				geneExprDictGTEX[entrezID] = {}
			for smpl in samplesGTEX:
				if not smpl in geneExprDictGTEX[entrezID]:
					geneExprDictGTEX[entrezID][smpl] = list()
				line_[c]=line_[c].strip()
				if line_[c]=="NA" or line_[c]=="None" :
					line_[c]=0
				geneExprDictGTEX[entrezID][smpl].append(float(line_[c]))
				c = c + 1
	
	AllBCGene = collections.OrderedDict()
	with open(args.bc_all_gene_entrezid, "r") as f:
		for line in f:
			line_ = line.split('\t')
			geneName1 = line_[0].strip()
			entrezID1 = line_[1].strip()
			geneName2 = line_[2].strip()			
			entrezID2 = line_[3].strip()

			if entrezID1 in linearCNADict and entrezID1 in geneExprDictBC and entrezID1 in geneExprDictGTEX:
				x,y,a,b=getGene_Amplification_Del(entrezID1, linearCNADict, geneExprDictBC, geneExprDictGTEX) #a return high sample list, b return low sample list
				print(str(entrezID1)+"\t"+str(geneName1), end='', flush=True)
				
				if entrezID2 in linearCNADict and entrezID2 in geneExprDictBC and entrezID2 in geneExprDictGTEX:
					b_del_pval,s_counts=getGeneB_Del (entrezID2,a,linearCNADict,geneExprDictBC, geneExprDictGTEX)
					b_amp_pval,s_counts_h=getGeneB_Amp (entrezID2,b,linearCNADict,geneExprDictBC, geneExprDictGTEX)
					geneAB_del_pval=getGeneA_Amp_GeneB_Del(entrezID1,entrezID2,linearCNADict,geneExprDictBC, geneExprDictGTEX)
					geneAB_amp_pval=getGeneA_Del_GeneB_Amp(entrezID1,entrezID2,linearCNADict,geneExprDictBC, geneExprDictGTEX)

					print("\t"+str(entrezID2)+"\t"+str(geneName2)+"\tAmplification\tDeletion\t"+str(x)+"\t"+str(b_del_pval)+"\t"+str(geneAB_del_pval)+"\t"+str(len(a))+"\t"+str(s_counts), flush=True)
					print(str(entrezID1)+"\t"+str(geneName1)+"\t"+str(entrezID2)+"\t"+str(geneName2)+"\tDeletion\tAmplification\t"+str(y)+"\t"+str(b_amp_pval)+"\t"+str(geneAB_amp_pval)+"\t"+str(len(b))+"\t"+str(s_counts_h))
				else:
					print("\t"+str(entrezID2)+"\t"+str(geneName2)+"\tAmplification\tDeletion")
					print(str(entrezID1)+"\t"+str(geneName1)+"\t"+str(entrezID2)+"\t"+str(geneName2)+"\tDeletion\tAmplification")
			else:
					print(str(entrezID1)+"\t"+str(geneName1)+"\t"+str(entrezID2)+"\t"+str(geneName2)+"\tAmplification\tDeletion")
					print(str(entrezID1)+"\t"+str(geneName1)+"\t"+str(entrezID2)+"\t"+str(geneName2)+"\tDeletion\tAmplification")

			if entrezID2 in linearCNADict and entrezID2 in geneExprDictBC and entrezID2 in geneExprDictGTEX:
				x,y,a,b=getGene_Amplification_Del(entrezID2, linearCNADict, geneExprDictBC, geneExprDictGTEX) #a return high sample list, b return low sample list
				print(str(entrezID2)+"\t"+str(geneName2), end='', flush=True)
				
				if entrezID1 in linearCNADict and entrezID1 in geneExprDictBC and entrezID1 in geneExprDictGTEX:
					b_del_pval,s_counts=getGeneB_Del (entrezID1,a,linearCNADict,geneExprDictBC, geneExprDictGTEX)
					b_amp_pval,s_counts_h=getGeneB_Amp (entrezID1,b,linearCNADict,geneExprDictBC, geneExprDictGTEX)
					geneAB_del_pval=getGeneA_Amp_GeneB_Del(entrezID2,entrezID1,linearCNADict,geneExprDictBC, geneExprDictGTEX)
					geneAB_amp_pval=getGeneA_Del_GeneB_Amp(entrezID2,entrezID1,linearCNADict,geneExprDictBC, geneExprDictGTEX)

					print("\t"+str(entrezID1)+"\t"+str(geneName1)+"\tAmplification\tDeletion\t"+str(x)+"\t"+str(b_del_pval)+"\t"+str(geneAB_del_pval)+"\t"+str(len(a))+"\t"+str(s_counts), flush=True)
					print(str(entrezID2)+"\t"+str(geneName2)+"\t"+str(entrezID1)+"\t"+str(geneName1)+"\tDeletion\tAmplification\t"+str(y)+"\t"+str(b_amp_pval)+"\t"+str(geneAB_amp_pval)+"\t"+str(len(b))+"\t"+str(s_counts_h))
				else:
					print("\t"+str(entrezID1)+"\t"+str(geneName1)+"\tAmplification\tDeletion")
					print(str(entrezID2)+"\t"+str(geneName2)+"\t"+str(entrezID1)+"\t"+str(geneName1)+"\tDeletion\tAmplification")
			else:
				print(str(entrezID2)+"\t"+str(geneName2)+"\t"+str(entrezID1)+"\t"+str(geneName1)+"\tAmplification\tDeletion")
				print(str(entrezID2)+"\t"+str(geneName2)+"\t"+str(entrezID1)+"\t"+str(geneName1)+"\tDeletion\tAmplification")
