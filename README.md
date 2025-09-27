# epitope_analyzer
This tool is used to analyze the epitope prediction results from multiple immunoinformatics software programs, identify peptide regions predicted by at least three tools in common, and output them in a standardized format.
# Epitope Analyzer

A tool for analyzing epitope prediction results from multiple immunoinformatics software programs.  
It identifies peptide regions predicted by at least two tools, analyzes overlaps, and outputs results in a standardized format.

---

## Features

- Input epitope sequence sites predicted by multiple software tools  
- Analyze how many tools predict each amino acid position  
- Output continuous regions predicted by at least two software tools  
- Analyze overlaps among different software prediction intervals  
- Support providing a protein sequence to display peptide amino acid composition  
- Support reading predefined peptide data from a file  
- Create an example peptide data file  
- Save analysis results to a text file  

---

## Usage

This tool requires data to be provided either via parameters or interactive mode. Running without parameters is not supported.  

### Interactive Mode (Recommended)

```bash
python epitope_analyzer.py --interactive
