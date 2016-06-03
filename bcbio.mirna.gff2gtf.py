#/usr/bin/python
import sys

inFile = open(sys.argv[1],'r')

for line in inFile:
#skip comment lines that start with the '#' character
  if line[0] != '#':
  #split line into columns by tab
    data = line.strip().split('\t')
            
    ID = ''
                
    #if the feature is a gene 
    if data[2] == "miRNA":
        #get the id
      ID = data[-1].split('ID=')[-1].split(';')[0]
                                   
      #modify the last column
      data[-1] = 'gene_id "' + ID + '"; transcript_id "' + ID + '"'
                                                                 
      #print out this new GTF line
      print '\t'.join(data)
