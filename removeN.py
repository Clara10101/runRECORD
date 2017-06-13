import argparse
import sys

parser = argparse.ArgumentParser()
parser.add_argument("-f", "--file", help="input file", type=argparse.FileType('r'), required=True)
args = parser.parse_args()

def removeExtension(file_name):
    return file_name.split('.')[0]

outFile = open(removeExtension(args.file.name) + '_n.fastq','w')

for i,line in enumerate(args.file):
    if i % 4 == 0:
        header = line
    elif i % 4 == 1:
        sequence = line
    elif i % 4 == 3:
        quality = line
        #write last four lines to output file
        #if they don't contain N
        if sequence.count('N') == 0:
            outFile.write(header + sequence + '+\n' + quality)

outFile.close()