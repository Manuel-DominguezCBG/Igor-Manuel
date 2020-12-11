######################################

 Igor-Manuel Programming Assessment

######################################


(THIS DOCUMENT IS A DRAFT, IT IS NOT A ULTIMATE README)


### Organization plan ###

This project involves different stages. 
Each stage has its own folder in this repository.
More specific details about each stage can be found in their respectively README file.
Short explanation of each stage can be read below.

Stages:

1) STAGE 1

I)Write a code that identify the 4 different ID Variants

II) This code needs to avoid repetition in case one CSV input may have more than one ID in the same row.

2) STAGE 2

I) Identify which column has Gene and GRCh38 coordinates.
To be able to correct Gene or coordinates, first we need a code that identify this values

3) STAGE 3

I) Check if Gene and coordinates are correct according some uptodate database
We will need to introduce error in this values.

II) Create a document that inform which values were identify wrongly 

III) Modify the incorretly values.

4) STAGE 4

I) Create a VCF with corrected data from the CSV.

5) STAGE 5 

I) Create a script that recognizes in CSV input and returns a VCF.

II) The script could ask which genome you want to use GRCh38 /GRCh37. (BONUS).

6) STAGE 6

I) Create a interface that performs the same execution that stage 5.



### IMPUT explanation ###

To achieve this project the next CSVs files have been created.
4 CSV files has been created, each one store the same 4 variants with different types of ID variant.
The columns of these files are
Gene, Coordinates, Variant_ID, Quality and Filter


Additionally, a fifth CSV has been created that contains the same variants with all 4 types of ID
Gene, Coordinates,Transcript,RefSeq Gene,Protein1, Protein2, Quality and Filter. 

The goal of creating these CVS is as follows.
Because in our labs, different clinicians store their variants with different IDs,
we are going to create a code that is capable of recognizing the most used types of ID.

If for some reason different types of ID are stored in the same row in the same CSV,
our code will be able to avoid repetitions. 
Therefore, we have created a fifth CSV with the 4 types of IDs used.
