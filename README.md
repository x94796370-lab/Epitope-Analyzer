# Epitope Analyzer

A tool for analyzing epitope prediction results from multiple immunoinformatics software programs.  
It identifies peptide regions predicted by at least two tools, analyzes overlaps, and outputs results in a standardized format.

---

## Features

- Input epitope sequence sites predicted by multiple software tools  
- Analyze how many tools predict each amino acid position  
- Output continuous regions predicted by at least three software tools  
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
In interactive mode, the program will guide you through the process:

Ask whether to use predefined peptide data and its source

Ask whether to provide a protein sequence and how to input it

Allow input of prediction intervals from each software tool

Using a Predefined Peptide Data File
bash
python epitope_analyzer.py --peptides-file peptides.txt
This reads peptide data from the specified file and directly outputs the results, without calculating prediction intervals.

Creating an Example Peptide Data File
bash
python epitope_analyzer.py --create-example
This creates a file named peptides_example.txt, which contains a sample format for peptide data.
You can edit this file and then use it via the --peptides-file parameter.

Providing a Protein Sequence
To display peptide amino acid sequences, you can provide a protein sequence as follows:

bash
# Combine peptide data file with a protein sequence
python epitope_analyzer.py --peptides-file peptides.txt --sequence "MTTQAPTASDKGLVY..."

# Combine interactive mode with a sequence file
python epitope_analyzer.py --interactive --sequence-file protein.txt
Analyzing Overlaps
bash
python epitope_analyzer.py --interactive --analyze-overlap
Adding this parameter displays details of overlapping regions between each pair of software tools.

Peptide Data File Format
The peptide data file should contain one peptide per line, with four fields:
name, sequence, start position, end position

Supported formats:

Comma-separated:

nginx
Peptide 1,EQLEQKLRDVEETKAKAE,23,40
Space-separated:

nginx
Peptide 1 EQLEQKLRDVEETKAKAE 23 40
Lines starting with # are treated as comments.

Software Prediction Interval Input Format
Prediction intervals should be written as:

sql
start-end,start-end,...
Example:

10-40, 44-52, 57-63
The program supports multiple delimiters: commas, semicolons, and spaces.

Output
The program will output:

A list of linear epitopes (similar to “Table 4”), including sequence name, amino acid sequence (if provided), position, and length

(Optional) Details of overlapping regions among different software prediction intervals

All shared region results will be displayed on screen and saved to a file named:

linear_epitopes.txt
If no shared regions are found, the program will notify you and no result file will be generated.
python epitope_analyzer.py --interactive
