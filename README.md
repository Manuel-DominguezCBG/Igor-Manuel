# CSV2VCF

###### Version 1.0
![image](https://drive.google.com/uc?export=view&id=1qTdMNRkowLjhYSBPqZSy3Lp6gcv4V-A_)

 Welcome to the script CSV2VCF-converter (v1).                        
 A collaborative project of the module Introduction to Programming.      
 By Igor Malashchuk and Manuel Dominguez  


### Description

  - This script converts variants stored in a Comma Separated Values (CSV) file into Variant Call Format (VCF) files

  - Also this script identifies, informs and corrects variant characteristics such as Gene_name and genome coordinates.
  
Following the agile methodology, this first version has been split into a few stages that can be read in the code.

### Motivation 
The motivation of this programme is to allow users without strong informatics background to convert variants stored in CSV files into VCF file. For this reason, our philosophy is to create a script as much flexible as possible trying to recognise the larges number of different variants ID. In this early version, 3 different variants IDs are recognised (see details below), and the programme can be run by using the terminal of any Linux or Mac computer. Our aim is to keep working on this project (even if only at the beginning for a didactic objective) and be able to identify the vast majority of variants ID as well as recognise different input. This is important because we have seen how different clinical scientists store variants not only using different variant ID but also creating Excel files in different structures. Some don't create a column with gene symbol, column header can be called in many different ways, the order of the column can change, coordinates can be stores like chr3:123456 or like 3:123456. Our final goal is to create a friendly interface that recognise a great diversity of variants ID and different data frame. 

### Installation and how to used

 - CSV2VCF requires any Linux or Mac OS distributions and [Python 3](https://www.python.org/) to run.

Very likely that if you have a Mac or Linux computer you have already installed Python, anyway here we pass you a  [good tutorial](https://realpython.com/installing-python/) to install Python in case you need it.

 - Then, cloned or download this repository. Go to the folder called CSVtoVCF in which you will see an input and output folder and  the executable file called CSVtoVCF_Run.py
  
 - Leave the .csv file you wish to convert in the input folder. At that point we remind you a couple of things:

  1. Called the column with the gene_symbol/gene_name with the name "GENE" and called the column with the location "COORDINATES". Coordinates must be like the following example "1:123456" If your file doesn't store gene names, don't worry, call one column with that name. VCF file will be created anyway.

  2. Save the file as .csv to ensure successful conversion.
 
 - Open the terminal and go to the directory where the script is located.
  If you don't know how to do this, have a look at  [this tutorial](https://www.youtube.com/watch?v=Vhcx4KJbtes&feature=emb_logo).

 - To ensure you are in the right directory try this, write ls in the terminal, you should see something like this.

```sh
(base) manolodominguez@Manolos-Mac-mini CSVtoVCF % ls
     CSVtoVCF_Run.py           input
     output
```
If you can read the content of the folder CSVtoVCF, that means you are in the right place.

 - You can run the programme writing in the terminal **python CSVtoVCF_Run.py** 
 
 - You will interact with the script, write an "a" if you are working with the genome building GRCh37 - hg37 or write a "b" if you stored in your .csv file variant in the genome building  GRCh37.
 
 - We will create a CSV in the genome building you wish, so the second question we ask you is to write an "a" if you want the results in the genome building GRCh37 - hg37 or a "b" if you wish your coordinates results in Grch38. We will not change the variant ID, only the coordinates in the VCF.
 
 - Finally, go to the output folder where you will see the VCF file and an HTML file informing you about errors found in your input.

### Development

1. CSV2VCFv1 recognises 3 types of variant ID

 - dbsnp ID (i.e. rs12345)
 - HGVS notation coding such as AGT:c.803T>C or ESNT000000000304:c.1431_1433delTTC

 We are working on developing version 2 which will recognise genomic (:g.) HGVS notation.

2. Too much code in one script. For organization purposes and effectiveness, the next approach will be split the code in different scripts.



 ### Testing and issues.

A test is running at the end of the script. The costumer is informed when convertion is done so that their dont need to wait is test is not needed. 

The test used the output generated in the VCF and makes new request by using a different API. We have decided to use MyVariantINFO API. We compare some values between both API is they match we considerar that this variant is correct.

If you find any issue running this programme, please let us know in [Github issues](https://github.com/Manuel-DominguezCBG/Igor-Manuel/issues) and we will try to help you. The same, if you have any idea or suggestion or you would like to join our team, contact with us. 


### The output

This programme creates 2 files:

 - A VCF file version 4.2. WE have followed the recommendation that can read [here]((https://samtools.github.io/hts-specs/VCFv4.2.pdf). This document will looks like that:
 
 ```sh
 ##fileformat=VCFv4.2
##fileDate=20090805
#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT NA00001 NA00002 NA00003
20 14370 rs6054257 G A 29 PASS NS=3;DP=14;AF=0.5;DB;H2 GT:GQ:DP:HQ 0|0:48:1:51,51 1|0:48:8:51,51 1/1:43:5:.,.
20 17330. T A 3 q10 NS=3;DP=11;AF=0.017 GT:GQ:DP:HQ 0|0:49:3:58,50 0|1:3:5:65,3 0/0:41:3
20 1110696 rs6040355 A G,T 67 PASS NS=2;DP=10;AF=0.333,0.667;AA=T;DB GT:GQ:DP:HQ 1|2:21:6:23,27 2|1:2:0:18,2 2/2:35:4
 ```


 - An HTML file with the next information: I) What variants has not been recognised by Ensembl API, II) which variant is not recognised by our script, III) genes names that don't match with the API and which coordinates don't match with the coordinates found in the API. 

### Learning

Our first reflection is that this project is more complex than we originally thought. The different amount of variants ID made it impossible to recognize all of them in one script. To achieve this in future versions, we need to think first about the organisation of the programme. At this point of our compression in programming, we think that the best approach would be to split the code in many scripts resolving problem by problem and keep develop.

The second challenge we have faced during this project has been the diversity of input we have seen. Clinical scientists not only used different variant ID but also CSV with different structures. This forced us to restart writing the script trying to recognise as many as possible differences inputs.

In my lab, I took the opportunity of showing our progress to clinical scientists (professionals without a bioinformatics background). Another difficulty I have noticed is that people that don´t understand programming, sometimes ask things that look easy for them but then put in practice is really difficult to apply. For example indel variants. Originally sounds easy to add this type of variant in our repertory of recognizable variants, but when you think how this can be done, many problems emerge.

Another element that brought us many hours of work was the modules. Normally when you run this program on different computers, we noticed that on one computer some modules need to be installed in another computer with another characteristic need to install other modules. To solve this we write a code that checks if all modules need are installed if one of them is not installed, we install the module, and then we import all them. This is not an easy thing especially when we want to avoid that the user interacts with the terminal. We have been thinking about it and we will try to run a new env with the moduled we need first and then run the script. Still thinking about this issue.


### License

During the development of this first version, this repository has be privated. As soon as this project will be evaluated, this will be free code. See License below.

MIT License

Copyright (c) [2021] [Igor Malashchuk and Manuel Dominguez]

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

----
