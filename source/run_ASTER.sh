#!/bin/bash
usage() { echo "Usage: $0 [-i <genes_ID.txt file>] [-s <data_CNA_TCGA.txt file>] [-g <gene_expression_counts_TCGA.txt file>] [-n <gene_expression_counts_GTEX.txt file>] [-c <p-value cutoff>]" 1>&2; exit 1; }

Help()
{
   # Display Help
   echo "bash run_ASTER.sh -i genes_ID.txt -s data_CNA_TCGA.txt -g gene_expression_counts_TCGA.txt -n gene_expression_counts_GTEX.txt -c 0.05"
   echo
   echo "Syntax: scriptTemplate [-h|i|g|n|c]"
   echo "options:"
   echo "h     Print this help"
   echo "i     gene list and ID"
   echo "g     gene expression (rsem count) - cancer samples"
   echo "n     gene expression (rsem count) - normal samples (GTEX)"
   echo "c     p-value cutoff, i.e. 0.05"
   echo
}

while getopts ":h:i:s:g:n:c:" option; do
   case $option in
      h) # display Help
         Help
         exit;;
      i) 
         i=${OPTARG}
	 ;;
      s)
	 s=${OPTARG}
	 ;;
      g)
         g=${OPTARG}
	 ;;
      n)
	 n=${OPTARG}
	 ;;
      c)
         c=${OPTARG}
	 ;;
      \?) # incorrect option
         echo "Error: Invalid option"
         exit;;
            
   esac
done

if [ -z "${i}" ] || [ -z "${s}" ] || [ -z "${g}" ] || [ -z "${n}" ] || [ -z "${c}" ]; then
    usage
fi

echo "Running..."
python3 ASTER.py ${i} ${s} ${g} ${n} > temp_results


awk '$7 < '$c' && $8 < '$c' && $9 < '$c' && $7!=""' temp_results > my_results
rm temp_results
echo "Done"

