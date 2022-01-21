#!/bin/bash
#SBATCH -J bigwig
#SBATCH -p CN_BIOT
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --output=%j.out
#SBATCH --error=%j.err
dataset=$1
input="/data/taoyuhuan/projects/exOmics_RNA/${dataset}/output/bam"
chrom_sizes="/data/taoyuhuan/projects/exOmics_RNA/bin/genome/genome/chrom_sizes/genome"
output="/data/taoyuhuan/projects/exOmics_RNA/level_4/${dataset}"
mkdir ${output}
mkdir ${output}/bam_sortbycoord
mkdir ${output}/bedgraph
mkdir ${output}/bigwig
temp="/data/taoyuhuan/projects/exOmics_RNA/level_4/${dataset}/temp"
mkdir ${temp}

echo "Start bam to bigwig convert for ${dataset} at `date`."
for sample in $(ls ${input})
do
echo "Start convert ${sample} at `date`."
#sort bam
echo "Sorting ${sample} bam..."
samtools sort -@ 16 -T ${temp}/${sample} -o ${output}/bam_sortbycoord/${sample}.genome_rmdup.sorted.coord.bam ${input}/${sample}/genome_rmdup.bam

#index bam
echo "Indexing ${sample} bam..."
samtools index -@ 16 ${output}/bam_sortbycoord/${sample}.genome_rmdup.sorted.coord.bam

#bam to bedgraph
echo "${sample} to bedgrapgh..."
bedtools genomecov -ibam ${output}/bam_sortbycoord/${sample}.genome_rmdup.sorted.coord.bam -bg -split | LC_COLLATE=C sort -T ${temp} -k1,1 -k2,2n > ${output}/bedgraph/${sample}.bedgraph
bedtools genomecov -ibam ${output}/bam_sortbycoord/${sample}.genome_rmdup.sorted.coord.bam -bg -split -strand + | LC_COLLATE=C sort -T ${temp} -k1,1 -k2,2n > ${output}/bedgraph/${sample}.+.bedgraph
bedtools genomecov -ibam ${output}/bam_sortbycoord/${sample}.genome_rmdup.sorted.coord.bam -bg -split -strand - | LC_COLLATE=C sort -T ${temp} -k1,1 -k2,2n > ${output}/bedgraph/${sample}.-.bedgraph

#bedgraph to bigwig
echo "${sample} to bigwig..."
bedGraphToBigWig ${output}/bedgraph/${sample}.bedgraph ${chrom_sizes} ${output}/bigwig/${sample}.bigwig
bedGraphToBigWig ${output}/bedgraph/${sample}.+.bedgraph ${chrom_sizes} ${output}/bigwig/${sample}.+.bigwig
bedGraphToBigWig ${output}/bedgraph/${sample}.-.bedgraph ${chrom_sizes} ${output}/bigwig/${sample}.-.bigwig

#normalize bigwig
echo "${sample} normalizing..."
read_depth=`bamtools count -in ${input}/${sample}/genome_rmdup.bam`
bigWigToBedGraph ${output}/bigwig/${sample}.bigwig stdout | awk -v d="$read_depth" 'BEGIN{{OFS="\t";FS="\t"}}{{print $1,$2,$3,1000000.0*2*$4/d}}' > ${output}/bedgraph/${sample}.normalized.bedgraph
bedGraphToBigWig ${output}/bedgraph/${sample}.normalized.bedgraph ${chrom_sizes} ${output}/bigwig/${sample}.normalized.bigwig

bigWigToBedGraph ${output}/bigwig/${sample}.+.bigwig stdout | awk -v d="$read_depth" 'BEGIN{{OFS="\t";FS="\t"}}{{print $1,$2,$3,1000000.0*2*$4/d}}' > ${output}/bedgraph/${sample}.+.normalized.bedgraph
bedGraphToBigWig ${output}/bedgraph/${sample}.+.normalized.bedgraph ${chrom_sizes} ${output}/bigwig/${sample}.+.normalized.bigwig

bigWigToBedGraph ${output}/bigwig/${sample}.-.bigwig stdout | awk -v d="$read_depth" 'BEGIN{{OFS="\t";FS="\t"}}{{print $1,$2,$3,1000000.0*2*$4/d}}' > ${output}/bedgraph/${sample}.-.normalized.bedgraph
bedGraphToBigWig ${output}/bedgraph/${sample}.-.normalized.bedgraph ${chrom_sizes} ${output}/bigwig/${sample}.-.normalized.bigwig

echo "End convert ${sample} at `date`."
done

echo "Finish convert for ${dataset} at `date`."
