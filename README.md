# MISO-Utilities
A python package which provides a suite of functions to aid in the generation and interpretation of MISO data. Package provides both plotting features and functions to consolidate and extract psi values from summary data.

# misoWrapper
misoWrapper.py provides a minimal higher-level wrapper around MISO (Mixture of Isoforms) API. For ease of use, the script.py is the companion to the misoWrapper and allows for function calls from the command line. 

#Usage
##Basic Usage

After the python packages have been downloaded and MISO has been locally installed, we can envoke the run MISO command by simply typing:

```bash
1) Use MISO to quantitate psi values: 

python script.py misoWrapper.runMISO bams/ indexedMisoEvents/ overhanglen False outDir misoSettings.txt all numProcessors sge/torque

2) Use MISO to compare psi values: 

python script.py misoWrapper.compareMISO misoDir/ comparisonDir/ all True

3) Summarize MISO quantitations: 

python script.py misoWrapper.summarize misoDir/ summaryDir/ sge

4) Summarize MISO comparison: 

python script.py misoWrapper.summarizeVs comparisonDir/ mapFile summaryDir/

If you want to consolidate and run the "monotonicity" test to find consistently changing events across replicates, or use a time course, do this:

1) Consolidate the information

python script.py misoWrapper.consolidateSummaries summaryDir/ groups_f consolidated.txt 

groups_f needs to be a file that contains a line for each sample, with the group that each sample should be in. For example:
MEFWT_2C5 WT
MEFWT_2E4 WT
MEFWT_4G4 WT
MEFKO_1F11 KO
MEFKO_2E8 KO
MEFKO_4B10 KO

It should have 2 columns, were there is a space or tab between each column.

2) Run the monotonicity test
python cript.py misoWrapper.monotonic consolidated.txt groups_f minbf nshuffles monotonic.txt 
minbf can be 5
nshuffles can be 100
The Z-score gives you the metric for "monotonicity".

````


