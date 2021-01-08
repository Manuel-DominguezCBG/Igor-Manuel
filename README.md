######################################

 #Igor-Manuel Programming Assessment

######################################


(THIS DOCUMENT IS A DRAFT, IT IS NOT A ULTIMATE README)


### Organization plan ###

This project involves different stages. 
Each stage has its own folder in this repository.
More specific details about each stage can be found in their respectively README file.
Short explanation of each stage can be read below.

Stages:

1) STAGE 1

I) Write a code that identify 2 (if time maybe more) different ID Variants

II) This code needs to avoid repetition in case one CSV input may have more than one ID in the same row.

III) By using the variant ID, request Emsembl API to uptodate data.

IV) In case a variant ID is unknown in the API, a exception avoid this issue

2) STAGE 2

I) Compare gene symbol and location between input and API

II) Generater a txt file that informs where and what mismatchs has been found


3) STAGE 3

I) Create a VCF file with data of the input plus from the API

4) STAGE 4 

I) Put all together, create a script that taking a CSV file, generate the inform txt plus the VCF file

5) STAGE 5

I) Create a interface that performs the same execution that stage 4.



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
