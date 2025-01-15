# Estimation for gene expression

#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH -o ./tmp/exp.o.%j
#SBATCH --err ./tmp/exp.e.%j

ml star
ml rsem

samp=${i##*/}
samp=${samp%%.bam}
fq1=${samp}_1.fq
fq2=${samp}_2.fq

rsem-calculate-expression --no-bam-output --paired-end -p 12 --star 2.fastq/1.fq/${fq1} 2.fastq/1.fq/${fq2} \
  ./4.rsem/WithOutpatch/rsem.hg19 ./5.rsem_result/${samp}.result
