# CSV2VCF
###### Version 1.0
![image](https://drive.google.com/uc?export=view&id=1qTdMNRkowLjhYSBPqZSy3Lp6gcv4V-A_)

 Welcome to the script CSV2VCF-converter (version 1)                        
 A collaborative assessment of the module Introduction to Programming      
 By Igor Malashchuk and Manuel Dominguez  


### Description

  - This script converts variants stored in a Comma Separated Values (CSV) file into Variant Call Format (VCF)

  - Also identify, inform and correct variant chararacteristic such as Gene_name and genome coordinates.
  
Following the agile methology, this first version has been splited in a few stages that can be read in the code.

### Motivation 
The motivation of this programme is to allow users without strong informatic background to convert variants stored in CSV files into VCF file. For this reason, our philosophy is to create a script as much flexible as possible trying to recognise the larges number of different variants ID. In this early version, 3 different variants ID are recognised (see details below), and the programme can be run by using the terminal of any Linux or Mac computer. Our aim is to keep working in this project (even if only at the beginning for a didactic objective) and be able of identify the vast majority of variants ID as well as recognise different input. This is important because we have seen how different clinical scientist store variant not only using different variant ID but also creating Excel files in different structure. Some don't create a column with gene symbol, column header can be called in many different ways, the order of the column can change, coordinates can be stores like chr3:123456 or like 3:123456. Our finall goal is to create an friendly interfae that recognise a great diversity of variants ID and different data frame. 

### Installation

 - CSV2VCF requires any Linux distributions or Mac OS and [Python 3](https://www.python.org/) to run.

Very likely that if you have a Mac or Linux computer you have already installed Python, anyway here we pass you a  [good tutorial](https://realpython.com/installing-python/).

 - Then, cloned of dowload this repository located the folder CSV2VCF in the directory you wish. Open the folder, you will see the you will see the script called CSVtoVCF_Run.py, an output and input folder and a README.txt
 
 - Leave the file you wish to convert in the input folder. At that point we remaind you to save the file in CSV formart. (open the excel file you want to convert, go to save as, and select .CSV). It is also important the order of the first 3 columns. Put gene symbol in the first column, coordinates in the second column and variant ID in the third column.
 Don't worry if your file doesn't have the gene symbol or gene name, put an empty column or any other column. The important thing is that. The rest of the data frame can have as many columns as you want and in the order you want.

 
 - Open the terminal and go to the directory where the script is located.
  If you don't know how to do this, have a look at  [this tutorial](https://www.youtube.com/watch?v=Vhcx4KJbtes&feature=emb_logo).

 - To ensure you are in the righ directory try this write ls in the terminal, you should see something like this.
```sh
(base) manolodominguez@Manolos-Mac-mini CSVtoVCF % ls
     CSVtoVCF_Run.py           input
     README.txt                output
```
If you can read the content of the folder CSVtoVCF, that means you are in the right place.

 - You can run the programme writing in the terminal **python CSVtoVCF** 
 
 - You will interact with the script, write a if you are working with the genome building GRCh38 - hg38 or be if you still work with the previous version.
 
 - We will create a CSV in the genome building you wish, so the second question we ask you is write a if you want the results in the genoma building GRCh38 - hg38 or b if you wish your coordinates results in Grch37. We will not change the variant ID, only the coordinates in the VCF.
 
 - Finally, go to the output folder where you will see the VCF file and a html file informing you about errors found in your input
