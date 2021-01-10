#!/usr/bin/env python
# coding: utf-8 

''' 
Let's import all modules used in this script
Our goal is that the costumer uses this script with the minimun
intereaction with the terminal. For that reason, we have added 
the next lines that install modules if required

Note for next version: Instead of install modules, 
create a new env could be more pythonic. 
''' 
import os
import subprocess

try:
    exec("import Cython")
except:
    subprocess.call(['pip', 'install', 'Cython'])

modules = ['sys', 'requests', 'chardet', 'json', 'pandas', 'numpy', 'Bio', 're', 'liftover', 'jinja2', "warnings"]

for library in modules:
    try:
        exec("import {module}".format(module=library))
    except:
        subprocess.call(['pip', 'install', library])
import pandas as pd
import numpy as np
import inspect
from datetime import date
import jinja2
import requests
from liftover import get_lifter
from Bio import Entrez, SeqIO

'''
To make clear the interection user-cumputer 
User will be ask to download the folder called CSVtoVCF
that contains the  CSVtoVCF_Run.py and a input and output folder
The next lines take a csv file from the imput folder
'''
# Let´s delete one warning 
warnings.filterwarnings("ignore")

''' Read output '''
# This read the file located in the folder input
actual_path = os.path.dirname(os.path.abspath(inspect.stack()[0][1]))

folder_input = '/input/'
folder_output = '/output/'

input_file_name = os.listdir(actual_path+folder_input)
input_file_name= ''.join(input_file_name)
conv38to19 = get_lifter('hg38', 'hg19')

CSV_input = actual_path+folder_input+input_file_name


'''
Now, we convert this input in a data frame 
'''


# To avoid UnicodeDecodeError: 'utf-8'
with open(CSV_input, 'rb') as f:
    result = chardet.detect(f.read())

df = pd.read_csv(CSV_input, encoding=result['encoding'])

# We will need a extra copy to avoid repetitions
df = df.rename(columns=str.upper)
df_draft = df

''' 
The users will be asked to called GENE the columns that contains gene_symbol 
and COORDINATES the column with the coordinates
that way we can recognise this input
In the next lines also, we ask what genome building is used in the input
'''

if bool("GENE" in df.columns) == True and bool("COORDINATES" in df.columns) == True:
    is_stage2 = 'YES'
    # Select input genome build for the CSV file
    genome_build = input('Please select genome build for the CSV input file: (a) GRCh37 or (b) GRCh38 ')
    while not (genome_build == 'a' or genome_build == 'b'):
        print('Genome build not supported.')
        genome_build = input('Please select genome build for the CSV input file: (a) GRCh37 or (b) GRCh38 ')
    if genome_build == 'a':
        print('Genome build selected: GRCh37')
    if genome_build == 'b':
        print('Genome build selected: GRCh38')
else:
    is_stage2 = 'NO'


# Select which genome build to produce the  VCF file in
genome_vcf_out= input('Please select genome build for the VCF file: (a) GRCh37 or (b) GRCh38 ')
while not (genome_vcf_out == 'a' or genome_vcf_out == 'b'):
    print('Genome build not supported.')
    genome_vcf_out = input('Please select genome build for the VCF file: (a) GRCh37 or (b) GRCh38 ')
if genome_vcf_out == 'a':
    print('Genome build selected: GRCh37')
if genome_vcf_out == 'b':
    print('Genome build selected: GRCh38')

# Input email address to retrieve data from Entrez
Entrez_ID = input('Please provide an email address (this will be used to retrieve data from Entrez)')
while Entrez_ID.find('@')==-1:
    print('Email address not entered correctly')
    Entrez_ID = input('Please provide an email address (this will be used to retrieve data from Entrez):')

if Entrez_ID.find('@')>0:
    print('You entered: ' + Entrez_ID)

def get_html_gdrive(ID):
    URL = "https://drive.google.com/uc?export=download"
    session = requests.Session()
    response = session.get(URL, params = { 'id' : ID }, stream = True)
    token = get_confirm_token(response)
    if token:
        params = { 'id' : ID, 'confirm' : token }
        response = session.get(URL, params = params, stream = True)
    return re.sub('[\r\t]','', response.text)
def get_confirm_token(response):
    for key, value in response.cookies.items():
        if key.startswith('download_warning'):
            return value
    return None
html_template = get_html_gdrive('15djn2ZyaaYm9lIlupQ_9J_od9g4eC_-m')

# Reading CSV file values and looking for variants IDs

# Our programme recognise e.g AGT:c.803T>C, ENST00000003084:c.1431_1433delTTC and rs17289390)

def find_trans(data_frame):
    pattern = re.compile('(:c.)|(rs\d+)')
    return data_frame[data_frame.apply(lambda x: x.str.contains(pattern))].dropna(how='all').dropna(axis=1, how='all')

Transcript = find_trans(df)

# Now, we save the results found in a dict key=index and value=variand ID
if Transcript.empty == False:
    ind = Transcript.index.to_list()
    vals = list(Transcript.stack().values)
    row2Transcript = dict(zip(ind, vals))
    # We need to remove the row where rs has been found to avoid repetitions
    # In case in same row more than one kind of ID Variant is stored
    for index, Transcript in row2Transcript.items():
        # This will be done in df_draft
        df_draft = df_draft.drop(index)

#print('Is the DataFrame empty? ', df_draft.empty)
# df_draft empty means every row contained a rs or :c. ID
# print('Is the DataFrame empty? ', df_draft.empty)



def get_API(transcripts, data_frame):
    df2 = data_frame
    df2['exist'] = ''
    trans_list = df2['exist']
    trans_list = trans_list.to_frame()
    trans_list.columns = ['ID']
    decoded = dict()
    server1 = 'https://rest.ensembl.org/vep/human/hgvs/'
    server2 = "https://rest.ensembl.org/vep/human/id/"
    server_end = '/expand=1;content-type=application/json'
    interrogation = '?'
    for row, transcript in transcripts.items():
        trans1 = ':c.'
        trans2 = 'rs'
        if trans1 in transcript:
            #r = requests.get(server1 + transcript, headers={"Content-Type": "application/json"})
            r = requests.get(server1 + transcript + interrogation + server_end)
            if r.status_code == 200:
                decoded[row] = r.json()
                df2['exist'][row] = 'YES'
                trans_list['ID'][row] = transcript
            else:
                df2['exist'][row] = 'NO'
        elif trans2 in transcript:
            #r = requests.get(server2 + transcript + interrogation, headers={"Content-Type": "application/json"})
            r = requests.get(server2 + transcript + interrogation + server_end)
            if r.status_code == 200:
                decoded[row] = r.json()
                df2['exist'][row] = 'YES'
                trans_list['ID'][row] = transcript
            else:
                df2['exist'][row] = 'NO'
        else:
            print("unknown transcript")
    trans_df = trans_list[~trans_list.ID.str.contains('^$')]
    return decoded, df2, trans_df

decoded, df_exist, trans_df = get_API(row2Transcript, df)
#trans_df2 = trans_df[~trans_df.ID.str.contains('^$')]
##### GENERATE Table with API request errors ####

#### This function adds colour highlights to rows where variant has been succesfully retrievd from the API or not
def highlight_error(s, num_columns):
    if s.exist == 'YES':
        return ['background-color: #c0ffba']*num_columns
    elif s.exist == 'NO':
        return ['background-color: red']*num_columns
    elif s.exist == '':
        return ['background-color: orange']*num_columns
    else:
        return ['background-color: blue']*num_columns

style1 = df_exist.style.apply(highlight_error, num_columns=len(df_exist.columns), axis=1)
df_html = style1.render()

# This function adds the data frame to the html template
def generate_report(html_results_table, html_template):
    todays_date = date.today().strftime('%d/%m/%Y')
    html_text = jinja2.Template(html_template)
    #html_text = html_text.render(todays_date=todays_date)
    elements = html_results_table.split('</style><table')
    width='</style><table class=results'
    new_html_table=elements[0]+width+elements[1]
    html_text = html_text.render(todays_date=todays_date, content=new_html_table)
    #html_text += new_html_table
    #html_text += '</div>'
    return html_text
html_file = generate_report(df_html, html_template)

# Let's add gene symbol, chromosome, genome location, reference and alteration
#pattern2 = re.compile('(NO)|(^$)')
#df_exist2 = df_exist[~df_exist.exist.str.contains(pattern2)]
def  get_location(decoded, transcripts):
    column_names = ['CHROM', 'POS_START', 'POS_END', 'ID']
    df_location = pd.DataFrame(columns=column_names)
    for i in decoded:
        trans = decoded[i]
        df_location.loc[i, 'CHROM'] = trans[0]["seq_region_name"]
        start = trans[0]["start"]
        end = trans[0]["end"]
        if start > end:
            df_location.loc[i, 'POS_START'] = end
            df_location.loc[i, 'POS_END'] = start
        else:
            df_location.loc[i, 'POS_START'] = start
            df_location.loc[i, 'POS_END'] = end
        df_location.loc[i, 'CHROM'] = trans[0]["seq_region_name"]

        df_location.loc[i, 'ID'] = transcripts.loc[i, 'ID']
    return df_location

df_location = get_location(decoded, trans_df)

     # Always tell NCBI who you are


def call_Entrez(chrom_ID, start, end, Entrez_ID):
    '''
    Call Entrez to retrieve sequence in fasta format
    '''
    Entrez.email = Entrez_ID
    handle = Entrez.efetch(db="nucleotide", id=chrom_ID, rettype="fasta",strand=1,seq_start=start,seq_stop=end)
    record = SeqIO.read(handle, "fasta")
    handle.close()
    return str(record.seq)

def stage1(decoded, df_location, Entrez_ID):
    chromosomes = {"1": "CM000663.2", "2": "CM000664.2", "3": "CM000665.2", "4": "CM000666.2",
                   "5": "CM000667.2", "6": "CM000668.2", "7": "CM000669.2", "8": "CM000670.2",
                   "9": "CM000671.2", "10": "CM000672.2", "11": "CM000673.2", "12": "CM000674.2",
                   "13": "CM000675.2", "14": "CM000676.2", "15": "CM000677.2", "16": "CM000678.2",
                   "17": "CM000679.2", "18": "CM000680.2", "19": "CM000681.2", "20": "CM000682.2",
                   "21": "CM000683.2", "22": "CM000684.2", "X": "CM000685.2", "Y": "CM000686.2",
                   "M": "MF737176", "MT": "MF737176"}
    column_names = ["Gene_API", 'CHROM', 'POS_HG38', 'POS_HG19', 'ID', 'REF', 'ALT', 'STRAND']
    df_stage1 = pd.DataFrame(columns=column_names)
    # extract specific data from API and data frame into new data frame
    for i in decoded:
        trans = decoded[i]
        # For Gene_symbol
        try:
            df_stage1.loc[i, 'Gene_API'] = trans[0]["transcript_consequences"][0]['gene_symbol']
        except KeyError:
            df_stage1.loc[i, 'Gene_API'] = '.'
        # For CHROMOSOME
        df_stage1.loc[i,'CHROM'] = trans[0]["seq_region_name"]
        # For transcript
        df_stage1.loc[i, 'ID'] = df_location.loc[i, 'ID']
        # For POSITION in GRCh38
        start = trans[0]["start"]
        end = trans[0]["end"]
        if start > end:
            start = end
        # Strand
        df_stage1.loc[i, 'STRAND'] = trans[0]['strand']
        # REF and ALT on +tive stand
        if df_stage1.loc[i, 'STRAND'] == 1:
            # For REFERENCE SEQUENCE
            ref1 = trans[0]['allele_string'].split('/')[0]
            alt1 = trans[0]['allele_string'].split('/')[1]
            if ref1 == '-':
                ref1_seq = call_Entrez(chromosomes[df_location.loc[i, 'CHROM']], df_location.loc[i, 'POS_START'],
                                       df_location.loc[i, 'POS_START'], Entrez_ID)
                df_stage1.loc[i,'REF'] = ref1_seq
                df_stage1.loc[i,'ALT'] = ref1_seq + alt1
                df_stage1.loc[i, 'POS_HG38'] = start
            elif alt1 == '-':
                ref2_seq = call_Entrez(chromosomes[df_location.loc[i, 'CHROM']], (df_location.loc[i, 'POS_START']-1),
                                       df_location.loc[i, 'POS_END'], Entrez_ID)
                df_stage1.loc[i, 'REF'] = ref2_seq
                df_stage1.loc[i, 'ALT'] = ref2_seq[0]
                df_stage1.loc[i, 'POS_HG38'] = start-1
                #update_position = conv38to19[df_location.loc[i, 'CHROM']][(df_location.loc[i, 'POS_START'] - 1)]
                #df_stage1.loc[i, 'POS_HG19'] = update_position[0][1]
            else:
                df_stage1.loc[i, 'REF'] = ref1
                df_stage1.loc[i, 'ALT'] = alt1
                df_stage1.loc[i, 'POS_HG38'] = start
            df_stage1.loc[i, 'STRAND'] = '+'
        # REF and ALT on -tive strand
        if df_stage1.loc[i, 'STRAND'] == -1:
            ref2 = trans[0]['allele_string'].split('/')[0]
            alt2 = trans[0]['allele_string'].split('/')[1]
            complement = {'A': 'T', 'G': 'C', 'C': 'G', 'T': 'A', '-': '-'}
            # For REFERENCE SEQUENCE
            ref3 = "".join(complement[c] for c in ref2)[::-1]
            # For Variant Sequence
            alt3 = "".join(complement[d] for d in alt2)[::-1]
            if ref2 == '-':
                ref3_seq = call_Entrez(chromosomes[df_location.loc[i, 'CHROM']], df_location.loc[i, 'POS_START'],
                                       df_location.loc[i, 'POS_START'], Entrez_ID)
                df_stage1.loc[i, 'REF'] = ref3_seq
                df_stage1.loc[i, 'ALT'] = ref3_seq + alt3
                df_stage1.loc[i, 'POS_HG38'] = start
            elif alt2 == '-':
                ref4_seq = call_Entrez(chromosomes[df_location.loc[i, 'CHROM']], (df_location.loc[i, 'POS_START']-1),
                                       df_location.loc[i, 'POS_END'], Entrez_ID)
                df_stage1.loc[i, 'REF'] = ref4_seq
                df_stage1.loc[i, 'ALT'] = ref4_seq[0]
                df_stage1.loc[i, 'POS_HG38'] = start - 1
                
            else:
                df_stage1.loc[i, 'REF'] = ref3
                df_stage1.loc[i, 'ALT'] = alt3
                df_stage1.loc[i, 'POS_HG38'] = start
           
            df_stage1.loc[i, 'STRAND'] = '-'
       
    return df_stage1


df_stage1 = stage1(decoded, df_location, Entrez_ID)

def get_hg19(df):
    for index, row in df.iterrows():
        chrom = row['CHROM']
        if chrom == 'MT':
            chrom = 'M'
        pos = row['POS_HG38']
        converted = conv38to19[chrom][pos]
        row['POS_HG19'] = converted[0][1]
    return df

df_stage1 = get_hg19(df_stage1)


############ End of the STAGE 1 ######################


# Stage 2

# First, let's prepare where our information file is going to be create.

# Introce in path the directory where we wish to create the new CVF
# save_path = '/Users/manolodominguez/Desktop/git-repos/STP_mini_projects/Igor-Manuel/Outputs'
# save_path = '/Users/Igor/Google Drive/Clinical_Bioinformatics_STP/Igor-Manuel_Project/Outputs'
# Introduce in name, the name of your new file


txt_file = "CSVtoVCF_ErrorReport.html"

output_inform = actual_path + folder_output + txt_file

# This will  create directory + file name
# completeName = os.path.join(save_path, name_of_file + ".txt")

# To compare input with API result
# We need to create a new column from CHROM + ':' + POS

#df_stage1 = df_stage1.applymap(str)

def get_genome(data_frame, genome):
    df_stage1 = data_frame.applymap(str)
    if genome == 'a':
        df_stage1['Location'] = df_stage1['CHROM'].str.cat(df_stage1['POS_HG19'], sep=':')
    if genome == 'b':
        df_stage1['Location'] = df_stage1['CHROM'].str.cat(df_stage1['POS_HG38'], sep=':')
    return df_stage1

def get_stage2(df_stage1, df_filtered, genome):
    df_stage1_final = get_genome(df_stage1, genome)
#df_stage1['Location'] = df_stage1['CHROM'].str.cat(df_stage1['POS'], sep=':')

    df_comparation = pd.concat([df_filtered[['GENE', 'COORDINATES']], df_stage1_final[['ID', 'Gene_API', 'Location']]], axis=1)
# df_comparation = df_comparation.set_index("Transcript")

# New two columns as result of the comparation
    df_comparation['result1'] = np.where(df_comparation.iloc[:, 0] == df_comparation.iloc[:, 3], 'OK', 'ERROR')

    df_comparation['result2'] = np.where(df_comparation.iloc[:, 1] == df_comparation.iloc[:, 4], 'OK', 'ERROR')

# We save here the rows where errors in gene symbol have been found
    gene_error = df_comparation[df_comparation['result1'].str.match('ERROR')]
    gene_error = gene_error.drop(gene_error.columns[[1, 4, 5, 6]], axis=1)

# We save here the rows where errors in location have been found
    location_error = df_comparation[df_comparation['result2'].str.match('ERROR')]
    location_error = location_error.drop(location_error.columns[[0, 3, 5, 6]], axis=1)
    return gene_error, location_error

def gene_loc_error(html_file, location_error, gene_error):
    variant2locations = location_error.set_index('ID').T.to_dict('list')
    # we save this information in a dict with the variant where the error has been found as key
    # Plus in a list the gene simbol of the input and the gene simbol of the API as value.
    gene2gene_API = gene_error.set_index('ID').T.to_dict('list')

    # If there is not any mismatch found, one document will be created saying
    No_error_found = 'Congratulation, no errors have been found in your CSV with regard to gene symbol and location\n'
    # If there are mismathes found, one document will be created saying
    Errors_found = 'If you are reading this message, then the CSVtoVCF application detected some error(s) or mismatch(s) in the Gene Symbol and/or the location of your variants as compared to data retrieved from API Emsembl (VEP).\n Below the differences found are shown.\n\n'
    # Plus the mismatch

    # Number of locations mismatches
    how_many_location_error = 'Number of location mismatch(s) found:', str(len(location_error))
    how_many_location_error = ' '.join(how_many_location_error)


    # Now, show these mismatchs in this order
    header_location = '<th class="col_heading level0 col0" >Your location</th> <th class="col_heading level0 col1" >Variant</th>  <th class="col_heading level0 col2" >API location</th>'
    html = html_file

# Now with these dicts we can write the inform_error.txt
# If both dicts are empty
# then create the txt saying there is not error found
# If one or both dict contain data
# that means mismatch found in the input

# Now the same with gene_symbol

# Number of gene_symbol mismatchs
    how_many_GS_error = 'Number of Gene Symbol mismatch(s) found:', str(len(gene_error))
    how_many_GS_error = ' '.join(how_many_GS_error)

# Now, show these mismatchs in this order

    header_gene_symbol = '<th class="col_heading level0 col0" >Your Gene_symbol</th> <th class="col_heading level0 col1" >Variant</th>  <th class="col_heading level0 col2" >API Gene_symbol</th>'


    location_title = "############# LOCATION'S MISMATCH(S) ############## \n\n"
    gene_title = "\n\n############# GENE_SYMBOL'S MISMATCH(S) #############\n\n"

# Here we introduce a condition, if dicts are empty
# write No_error_found
# if dicts are not empty, write all this information


    if (bool(variant2locations) == False) and (bool(gene2gene_API) == False):
        html += '<div class=errors_summary><b>\n'+ No_error_found +'</b></div>\n'
        html += "</body>"
        html += "</html>"
    #    with open(output_inform, 'w') as out:
    #        out.writelines([No_error_found])
    else:
        print("errors found - refer to output error report")
    # Write location

    #with open(output_inform, 'w') as out:
    #    out.writelines([Errors_found, location_title, how_many_location_error, separator, header_location])
    #location_error.to_csv(output_inform,
    #                      header=None, index=None, sep='\t', mode='a')

    # Now, with Gene_symbol

    #with open(output_inform, 'a') as out:
    #    out.writelines([gene_title, how_many_GS_error, separator, header_gene_symbol])
    #gene_error.to_csv(output_inform,
    #                  header=None, index=None, sep='\t', mode='a')
        html += '<div class=errors_summary>\n'
        html += '<p>'+Errors_found+'</p>'
        html += '<div class=error_sum1>\n'
        html += '<p style="text-align:center">'+ gene_title +'</p>\n'
        html += '<p style="text-align:center">'+ how_many_GS_error +'</p>\n'
        html += '<table class=errors>\n'+ header_gene_symbol
        for j in range(len(gene_error)):
            html+= '<tr>\n'
            for col2 in gene_error.columns:
                value2 = gene_error.iloc[j][col2]
                html +='<td>\n'+str(value2)+'\n</td>\n'
            html+= '</tr>\n'
        html +='</table>\n'
        html += '</div>\n'
        html += '<div class=error_sum2>\n'
        html += '<p style="text-align:center">'+ location_title +'</p>\n'
        html += '<p style="text-align:center">'+ how_many_location_error +'</p>\n'
        html += '<table class=errors>\n'+ header_location
        for i in range(len(location_error)):
            html+= '<tr>\n'
            for col in location_error.columns:
                value = location_error.iloc[i][col]    
                html +='<td>\n'+str(value)+'\n</td>\n'
            html+= '</tr>\n'
        html +='</table>\n'
        html += '</div>\n'
        html += '</div>\n'
        html += '</body>\n'
        html += "</html>"
    return html

if is_stage2 == 'YES':
    pattern = re.compile('(NO)|(^$)')
    df_filtered = df_exist[~df_exist.exist.str.contains(pattern)]
    gene_error, location_error = get_stage2(df_stage1, df_filtered, genome_build)
    html_final = gene_loc_error(html_file, location_error, gene_error)
    with open (output_inform, 'w') as out:
        out.write(html_final)

#  End of stage 2

#################### Skip stage two #######################

# get error report as html file html file
def finish_html(html_file):
    text = "End of error report. CSVtoVCF can also check if your GENE symbol and chromosome " \
           "coordinates match with the GENE and location retrieved by API. " \
           "To enable this function make sure that your gene names are located under a column with a title 'GENE'" \
           "and chromosome coordinates matches example formal (e.g. 17:656767), the column title for this must be 'COORDINATES'."
    html_file += '<div class=errors_summary><b>\n'+ text +'</b></div>\n'
    html_file += "</body>"
    html_file += "</html>"
    return html_file
if is_stage2 == 'NO':
    finished_html = finish_html(html_file)
    with open (output_inform, 'w') as out:
        out.write(finished_html)




# # STAGE 3

# The VCF create for us follow the  structure explained in https://samtools.github.io/hts-specs/VCFv4.2.pdf
# We have created the minimum structure that a VCF document must have
# according to documentation (VCFv4.3 format)

# An example here

#      ##fileformat=VCFv4.3
#      #CHROM POS      ID         REF   ALT    QUAL  FILTER   INFO
#      20     14370    rs6054257  G     A      29    PASS    NS=3;DP=14;AF=0.5;DB;H2

# Here we explain how we are going to fill the columns of the  VCF file.

### Mandatory columns ###
# CHROM taken from df_stage1['CHROM']
# POS taken from df_stage1['POS']
# ID taken df_stage1['ID']
# REF taken from df_stage1['REF']
# ALT taken from df_stage1['ALT']
# QUAL taken from df_stage1['QUAL']
# FILTER taken from df_stage1['FILTER']
# INFO filled with (‘.’)

# First, let's create a empty document in the same directory

# Introce in path the directory where we wish to create the new CVF
# save_path = '/Users/manolodominguez/Desktop/git-repos/STP_mini_projects/Igor-Manuel/Outputs'

# Select which genome build to produce the  VCF file in





# Introduce in name, the name of your new file
name_of_VCF = "your_" + input_file_name + "_converted_into_VCF.vcf"

output_inform = actual_path + folder_output + name_of_VCF

# This will  create directory + file name
# completeName2 = os.path.join(save_path, name_of_file2 + ".vcf")

# Now, let's create the metadata information
file_format = '##fileformat=VCFv4.3\n'
today = date.today()
d1 = today.strftime("%Y-%m-%d")
file_Date = '##fileDate=' + d1 + '\n'

# And now the header of the columns
#header_line = '#CHROM\tPOS\tID\tREF\tALT\tSTRAND\tQUAL\tFILTER\tINFO\n'
header_line = '#CHROM\tPOS\tID\tREF\tALT\tSTRAND\tGENE\tINFO\n'

# Now, we write meta-information line and header line
with open(output_inform, 'w') as out:
    out.writelines([file_format, file_Date, header_line])

# Let's modify the df_stage1 to be added to the VCF file
df_stage2 = df_stage1
df_stage2['INFO'] = '.'

def convert(v):
    try:
        return int(v)
    except ValueError:
        return v


def get_vcf_genome(data_frame, genome_build):
    choice = {"a": "POS_HG19", "b": "POS_HG38"}
    data_frame2 = pd.DataFrame([convert(c) for c in l] for l in data_frame.values).sort_values([1, 2], ascending=(True, True))
    data_frame2.columns = data_frame.columns
    data_frame2 = data_frame2.sort_values(['CHROM', choice[genome_build]])
    #data_frame['CHROM'] = data_frame['CHROM'].astype(float)
    #data_frame[choice[genome_build]] = data_frame[choice[genome_build]].astype(float)
    #data_frame = data_frame.sort_values([ 'CHROM', choice[genome_build]], ascending=(True, True))
    #data_frame = data_frame[['CHROM', choice[genome_build], 'ID', 'REF', 'ALT', 'STRAND', 'QUAL', 'FILTER', 'INFO', 'Gene_API']]
    data_frame2 = data_frame2[
        ['CHROM', choice[genome_build], 'ID', 'REF', 'ALT', 'STRAND', 'Gene_API', 'INFO']]
    data_frame2 = data_frame2.replace(to_replace="nan",
                              value=".")
    return data_frame2

df_stage3 = get_vcf_genome(df_stage2, genome_vcf_out)

df_stage3.to_csv(output_inform,
                 header=None, index=None, sep='\t', mode='a')

print('Ta-daa!')
print('Have a look in the output folder')
print('Your files are ready!')



