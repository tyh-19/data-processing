#!/usr/bin/env python
import subprocess
import argparse
import os
import io

def align(fastq,index,bam,unmapped,log,threads=1):
    mapping_cmd = ["bowtie2","-p",str(threads),"--norc","--sensitive","--no-unal","--un-gz",unmapped,"-x",index,fastq,"-S","-"]
    #print(" ".join(mapping_cmd))
    sort_cmd = ["samtools","sort","--output-fmt","BAM","-o",bam]
    ps = subprocess.Popen(mapping_cmd, stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    _ = subprocess.check_output(sort_cmd, stdin=ps.stdout)
    ps.wait()
    _ = subprocess.run(["samtools","index",bam])
    with open(log,"w") as f:
        for line in io.TextIOWrapper(ps.stderr, encoding="utf-8"):
            f.write(line)
    return ps.poll() 
 
def main():
    parser = argparse.ArgumentParser(description='Calculate Secondary Structure Entropy of Given Regions')
    parser.add_argument('--fastq','-f',required=True,help="Input single end cleaned fastq")
    parser.add_argument('--bam-dir','-bd',required=True,help="Dir for output bam files")
    parser.add_argument('--unmapped-dir','-fd',required=True,help="Dir for output unmapped fastq files")
    parser.add_argument('--index-dir','-id',default="index",help="Dir contains bowtie index")
    parser.add_argument('--log-dir','-ld',required=True,help="Dir for output log files")
    parser.add_argument('--priority','-p',default="spikein_small,univec,rRNA,mtRNA,miRNA,lncRNA,mRNA,piRNA,snoRNA,snRNA,srpRNA,tRNA,tucpRNA,Y_RNA,genome,circRNA",help="Priority of mapping")
    parser.add_argument('--threads','-t',type=int,default=1,help="Number of threads for mapping")
    args = parser.parse_args()
    fastq = args.fastq
    if not os.path.exists(args.bam_dir):
        os.makedirs(args.bam_dir)
    if not os.path.exists(args.unmapped_dir):
        os.makedirs(args.unmapped_dir)
    if not os.path.exists(args.log_dir):
        os.makedirs(args.log_dir)
    for sequence in args.priority.split(","):
        print("Start align {}".format(sequence))
        bam = os.path.join(args.bam_dir,sequence+".bam")
        unmapped = os.path.join(args.unmapped_dir,sequence+".fastq.gz")
        index = os.path.join(args.index_dir,sequence)
        log = os.path.join(args.log_dir,sequence+".log")
        if not os.path.exists(bam):
            assert align(fastq,index,bam,unmapped,log,threads=args.threads) == 0
            print("Done .")
        fastq = unmapped
        

if __name__ == "__main__":
    main() 
