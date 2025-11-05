
We provide instructions, codes and datasets for replicating the paper entitled "A Practical Guide to Endogeneity Correction Using Copulas." 

===== A guide on materials in the package and instructions to use the code to replicate results in the paper =====
=====
The package contains the following files:

1. Sim1-ModelwithIntercept.R
   This is the R code for replication of a simulation study to evaluate the proper construction of copula transformation. 

2. Sim2-TScopes.R
   This is the R code for a simulation study to demonstrate 2sCOPE/2sCOPE-np. 
  
4. Examples.R
   This is the R  Code for replication of empirical examples results in the main manuscript.  

4. Data Files for empirical examples: Example1.csv and Example2.csv

	Description of variables:
        logVol -  The weekly log sale volume;
        logPrice - The weekly log Sale price
        Fshare - The weekly Featuring intensity
        week - the week number
        Q2,Q3 and Q4 - indicator variables for quarters 2,3, and 4. 
Data are perturbed to protect original data values. They generate very close results but not exactly the same results as those presented in the paper.
Despite minor numerical differences, the main conclusions using these perturbed data remain the same as using the original data set.



To run the program, download all files and save them under the same folder in Windows operating system.
Then source Sim1-ModelwithIntercept.R, Sim2-TScopes.R, and Examples.R in R  to run the programs. 

===== A list of the version of R and computer specification used in the work =====
R version: 4.3.2 (2023-10-31 ucrt)
Processor: Intel(R) Core(TM) i9-9940X CPU @ 3.30GHZ
Installed RAM: 128.0 GB
To run the program, download all files and save them under the same folder. Then source the R files to run the programs. 
