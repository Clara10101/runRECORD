import os
import subprocess
import multiprocessing
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-f","--file",help="SAM file",type=str, required=True)
parser.add_argument("-r","--record",help="Path to RECORD",type=str, required=True)
parser.add_argument("-i","--imca",help="Path to IMCA",type=str, required=True)
parser.add_argument("-v","--velvet",help="Path to Velvet",type=str, required=True)
parser.add_argument("-m","--mummer",help="Path to MUMmer",type=str, required=True)
parser.add_argument("-o","--outFile",help="Output file name",type=str, required=True)
parser.add_argument("-r","--ref",help="Reference files directory",type=str, required=True)
args = parser.parse_args()

fileName, fileExtension = os.path.splitext(args.file)
if fileExtension and fileExtension not in ['.sam']:
    raise Exception('Error - not valid file type')

if not (os.path.isdir(args.ref) and os.path.isdir(args.imca) and os.path.isdir(args.record) and os.path.isdir(args.velvet) and os.path.isdir(args.mummer)):
    raise Exception('Error - not valid ref directory')

if args.ref.endswith('/'):
    args.ref = args.ref[:-1]

if args.imca.endswith('/'):
    args.imca = args.imca[:-1]

if args.record.endswith('/'):
    args.record = args.record[:-1]

if args.velvet.endswith('/'):
    args.velvet = args.velvet[:-1]

if args.mummer.endswith('/'):
    args.mummer = args.mummer[:-1]

#Plikiem wejsciowym dla programu jest wynik mapowania odczytow na genom referencyjny z zamaskowanymi obszarami repetytywnymi
chr = [str(i) for i in range(1,23)]
chr.extend(['X','Y','M'])

output_file = open(args.outFile,'w')
cwd = os.getcwd()
workspace = os.path.dirname(cwd + "/workspace/")
numProcessors = multiprocessing.cpu_count()


def runLn(file_from, file_to):
    command = 'ln ' + str(file_from) + ' ' + str(file_to)

    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
    process.wait()

def runProcess(command):
    process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
    process.wait()

    if process.returncode == 1:
        exit(1)

def removeExtension(file_name):
    return file_name.split('.')[0]

def runRECORD(c, file_name):
    #Przygotowanie pliku fastq tylko z obszarami nierepetytywnymi
    command = "samtools view -b -h " + file_name + "_sorted.bam chr" + c + " -o workspace/chr" + c + "/" + file_name + "_sorted_chr" + c + ".bam"
    runProcess(command)
    command = "bedtools bamtofastq -i workspace/chr" + c + "/" + file_name + "_sorted_chr" + c + ".bam -fq workspace/chr" + c + "/" + file_name + "_sorted_chr" + c + ".fastq"
    runProcess(command)

    #Polecenia wywolywane przez RECORDA oraz wywolanie IMCA i dodanie adnotacji dbSNP
    if os.path.isfile("workspace/chr" + c + "/ref/chr" + c + ".fa"):
        # Generowanie pseudoodczytow
        command = "java -Xmx3g " + args.record + "/SampleReference workspace/chr" + c + "/ok.txt -paired workspace/chr" + c + "/ref/chr" + c + ".fa 10 100 1000 30 workspace/chr" + c + "/results/pseudoreads1.fastq workspace/chr" + c + "/results/pseudoreads2-a.fastq"
        runProcess(command)

        command = "perl " + args.record + "/reverse_complement.pl workspace/chr" + c + "/ok.txt workspace/chr" + c + "/results/pseudoreads2-a.fastq workspace/chr" + c + "/results/pseudoreads2.fastq -r -c"
        runProcess(command)

        # usuniecie pseudoodczytow zawierajacych N
        command = "python removeN.py -f workspace/chr" + c + "/results/pseudoreads1.fastq"
        runProcess(command)
        command = "python removeN.py -f workspace/chr" + c + "/results/pseudoreads2.fastq"
        runProcess(command)

        # Velveth
        command = args.velvet + "/velveth workspace/chr" + c + "/results/velvet_assembly/ 27 -fastq -short workspace/chr" + c + "/" + file_name + "_sorted_chr" + c + ".fastq -shortPaired workspace/chr" + c + "/results/pseudoreads1_n.fastq workspace/chr" + c + "/results/pseudoreads2_n.fastq"
        runProcess(command)

        #Velvetg
        command = args.velvet + "/velvetg workspace/chr" + c + "/results/velvet_assembly/ -exp_cov auto -read_trkg yes"
        runProcess(command)

        # MUMmer
        command = args.mummer + "/nucmer -p workspace/chr" + c + "/results/alignment/alignment workspace/chr" + c + "/results/velvet_assembly/contigs.fa workspace/chr" + c + "/ref/chr" + c + ".fa"
        runProcess(command)
        command = args.mummer + "/nucmer -p workspace/chr" + c + "/results/alignment/alignment workspace/chr" + c + "/ref/chr" + c + ".fa workspace/chr" + c + "/results/velvet_assembly/contigs.fa"
        runProcess(command)

        #Show-coords
        command = args.mummer + "/show-coords workspace/chr" + c + "/results/alignment/alignment.delta > workspace/chr" + c + "/results/alignment/showcoords"
        runProcess(command)

        # IMCA
        command = args.imca + "find_variants.sh -w " + cwd + "/workspace/chr" + c + " -m " + args.mummer + " -o " + cwd + "/workspace/chr" + c + "/wyniki_mica_chr" + c + ".vcf"
        runProcess(command)

        print 'ok'


#Przygotowanie struktury katalogow workspace
#Dla kazdego chromosomu tworzony osobny katalog w ktorym zawarty jest katalog ref i results
for c in chr:
    directory = workspace + "/chr" + c
    if not os.path.exists(directory):
        os.makedirs(directory)
        os.makedirs(directory + "/ref")

        if os.path.exists(args.ref + "/chr" + str(c) + ".fa"):
            print 'ln'
            runLn(args.ref + "/chr" + str(c) + ".fa",workspace + "/chr" + str(c) + "/ref/chr" + str(c) + ".fa")

        os.makedirs(directory + "/results")
        os.makedirs(directory + "/results/velvet_assembly")
        os.makedirs(directory + "/results/alignment")


#Dla kazdego chromosomu z pliku SAM wybierane sa odczyty zmapowane tylko na ten chromosom i zapisywane do pliku fastq
#Przekonwertowaniu pliku SAM na BAM
file_name = removeExtension(args.file)
runProcess("samtools view -bS " + file_name + ".sam -o " + file_name + ".bam")
#Sortowanie pliku BAM
runProcess("samtools sort " + file_name + ".bam -o " + file_name + "_sorted.bam")
runProcess("samtools index " + file_name + "_sorted.bam")

#Wywolanie funkcji runRECORD dla kazdego z chromosomow
pool = multiprocessing.Pool(numProcessors)

tasks = []
for c in chr:
    tasks.append((c,file_name))

results = [pool.apply_async(runRECORD, t) for t in tasks]

for i, result in enumerate(results):
    print("Result for chr %d" % (i))

pool.close()
pool.join()

#Zebranie i zapisanie wynikow IMCA
output_file.write('chr\twszystkie warianty\tSNP\tindele\n')
for c in chr:
    if os.path.isfile("workspace/chr" + c + "/wyniki_mica_chr" + c + ".vcf"):
        command = "bcftools stats workspace/chr" + c + "/wyniki_mica_chr" + c + ".vcf | grep 'records\|SNPs\|indels' | head -3"
        output = subprocess.check_output(command, shell=True).split('\n')

        records = output[0].split('\t')[3]
        snps = output[1].split('\t')[3]
        indels = output[2].split('\t')[3]
        output_file.write('chr' + c + '\t' + str(records) + '\t' + str(snps) + '\t' + str(indels) + '\n')
