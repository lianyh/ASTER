# ASTER: A Method to Predict Clinically ActionableSynthetic Lethal Interactions

--------------------------------
*Prerequisite*
--------------------------------
Please install:<br/>
python >= 3.6 version<br/>
python lib: scipy.stats <br/>
python lib: adaFDR <br/>

ASTER - Manual
------------------

Example Datasets
-------------------
Please refer to "data" folder for a list of files (or format) required for the input to ASTER program.

Command
-------------------
bash source/run.ASTER.sh -i data/genes_ID.txt -s data/data_CNA_TCGA.txt -g data/gene_expression_counts_TCGA.txt -n data/gene_expression_counts_GTEX.txt -c 0.05 

Input parameters:<br/>
i – A list of gene pairs with HUGO Gene Name and entrez ID <br/>
s – Normalized SCNA data in a range of -2 to 2 (TCGA samples)<br/>
g – RNA-Seq gene expression (rsem count) (TCGA) <br/>
n – GTEX gene expression (rsem count) (normal) <br/>
c - p-value cutoff, i.e. 0.05 <br/>


<br />
Please contact the author of this program Herty Liany, email: e0146315@u.nus.edu if you encountered any problems. <br /> Thank you.
<br />#(2020)

