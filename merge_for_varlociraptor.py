import os 
import fnmatch
import gzip 
from optparse import OptionParser


Parser = OptionParser()
Parser.add_option("--callers", action='store', dest='callers',help="A list of variant callers to be included.")
Parser.add_option("--output", action='store', dest='output_name', help = "Name of outfile")
Parser.add_option("--dir", action='store', dest='analysis_dirPath', help= "Tha PATH of the directory where vcf files locate.")
#Parser.add_option('-r', "--resultLUID", action='store', dest='resultLUID')

#(options, args) = Parser.parse_args()
options=  Parser.parse_args()[0]

callers= options.callers
callers = callers.split(',')
analysis_dirPath = options.analysis_dirPath
fout_name = options.output_name

print(analysis_dirPath, callers, fout_name)
all_vcfs =[]
for name in callers:
    pattern = '*' + name + '.vcf.gz'
    globals()[name] = fnmatch.filter(os.listdir(analysis_dirPath), pattern)
    
    files = fnmatch.filter(os.listdir(analysis_dirPath), pattern)
    all_vcfs = all_vcfs + files

#analysis_dirPath = '/fs1/results/myeloid/vcf/'
#freebayes = fnmatch.filter(os.listdir(analysis_dirPath), '*freebayes.vcf.gz')
#vardict = fnmatch.filter(os.listdir(analysis_dirPath), '*vardict.vcf.gz')
#tnscope = fnmatch.filter(os.listdir(analysis_dirPath), '*tnscope.vcf.gz')
#all_vcfs= freebayes +  vardict +  tnscope





header_lines=[]
lines= []
tags=[]

for  file in all_vcfs:
    file_path = analysis_dirPath + file
    print(file_path)
    with gzip.open(file_path, 'rt') as fin, open(fout_name , 'w') as merged_file:
        for line in fin:
                line = line.rstrip()
                #try to extract INFO/contig tag lines from headers.
                if  line.startswith('##'):
                    if line.startswith('##contig') or line.startswith('##INFO')  or  line.startswith('ALT'):
                        if line.split('<')[1].split(',')[0] not in tags:
                            tags.append(line.split('<')[1].split(',')[0])
                            header_lines.append(line)
                    
                elif  line.startswith('#CHROM'):
                        continue
                #try to get file contents
                else:
                    cols=line.split('\t')
                    new_line = '\t'.join(cols[0:5])
                    #check if  'SVLEN' is present
                    if  'SVLEN' in cols[7]:
                        new_line = new_line + '\t' + cols[7]
                    else:
                        new_line =  new_line + '\t' + '.'
                    if new_line not in lines:
                        lines.append(new_line)
        #collect all lines
        all_lines = header_lines + ['#CHROM\tPOS\tID\tREF\tALT\tINFO'] +  lines
        for line in all_lines:
            print(line, file= merged_file) #print lines into output file


print('#header lines' , len(header_lines))
print('#lines', len(lines))


'''
def main():
    
    This script merges the vcf files from different callers which are \n
    given as options. Please notic that the name of each vcf file has to satisfy this pattern.
    '*(nameofcaller).vcf.gz , e.g(*vardict.vcf.gz)'.
    The name of out file and the durectory where all vcf are located it are required to be passed as options.
    In the output file, All unique varinat lines and  contig,INFO heaer lines are saved. 
       
''' 
