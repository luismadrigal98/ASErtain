#!/usr/bin/env bash
# ASErtain end-to-end driver for the kunthii x amphorellae anthocyanin/FLS study.
#
# Run from the data directory (For_Luis_from_Haylee/). Does: index BAMs, merge
# RNA flowers per genotyping unit, build a candidate-gene regions BED, joint
# variant-call over those regions, then run the ASErtain pipeline (report-mode
# bias) restricted to the anthocyanin loci.
#
# Usage:
#   cd /path/to/For_Luis_from_Haylee
#   bash run_pkunthii_ase.sh /path/to/ASErtain/examples/pkunthii_anthocyanin_ase.yaml
set -euo pipefail

CONFIG="${1:?pass the ASErtain YAML config as the first argument}"
THREADS="${THREADS:-8}"

REF=PGA_assembly_shortnames.fasta
GFF=Pkunthii_annotation_viaHelixer.gff
BAMS=BAM_files
LOCI=anthocyanin_loci.txt
WORK=work
mkdir -p "$WORK"

echo "############ 0. reference + flower-BAM indexes ############"
samtools faidx "$REF"
for b in "$BAMS"/*_Aligned.bam; do
    [ -f "$b.bai" ] || samtools index -@ "$THREADS" "$b"
done

echo "############ 1. merge RNA flowers per genotyping unit ############"
# One merged BAM per VCF column. Pooling flowers deepens coverage so the
# outbred parents (and the F1 genomes) are genotyped more accurately.
merge () { local out="$1"; shift; samtools merge -f -@ "$THREADS" "$WORK/$out.bam" "$@"; samtools index "$WORK/$out.bam"; }
merge k2          "$BAMS"/k2_cor_29_Aligned.bam "$BAMS"/k2_cor_29B_Aligned.bam "$BAMS"/k2_cor_31_Aligned.bam "$BAMS"/k2_cor_33_Aligned.bam
merge k3          "$BAMS"/k3_cor_23_Aligned.bam "$BAMS"/k3_cor_25_Aligned.bam "$BAMS"/k3_cor_27_Aligned.bam
merge amphorellae "$BAMS"/a_cor_1B_Aligned.bam "$BAMS"/a_cor_3_Aligned.bam "$BAMS"/a_cor_5_Aligned.bam "$BAMS"/a_cor_7_Aligned.bam "$BAMS"/a_cor_9_Aligned.bam "$BAMS"/a_cor_10_Aligned.bam
merge F1_k2a1     "$BAMS"/k2a_f1_cor_11_Aligned.bam "$BAMS"/k2a_f1_cor_13_Aligned.bam "$BAMS"/k2a_f1_cor_15B_Aligned.bam
merge F1_k3a1     "$BAMS"/k3a_f1_cor_17_Aligned.bam "$BAMS"/k3a_f1_cor_19B_Aligned.bam "$BAMS"/k3a_f1_cor_21_Aligned.bam
printf "%s\n" "$WORK"/k2.bam "$WORK"/k3.bam "$WORK"/amphorellae.bam "$WORK"/F1_k2a1.bam "$WORK"/F1_k3a1.bam > "$WORK/bamlist.txt"

echo "############ 2. candidate-gene regions BED ############"
# anthocyanin_loci.txt lines look like  <mRNA>_gene:<geneID>  -> take <geneID>,
# then look its coordinates up in the Helixer GFF gene features.
cut -d: -f2 "$LOCI" | tr -d ' \r' | sed '/^$/d' > "$WORK/cand_ids.txt"
awk -F'\t' 'NR==FNR{c[$1]=1; next}
            $3=="gene"{id=$9; sub(/^ID=/,"",id); sub(/;.*/,"",id);
                       if(id in c) printf "%s\t%d\t%d\t%s\n",$1,$4-1,$5,id}' \
    "$WORK/cand_ids.txt" "$GFF" | sort -k1,1 -k2,2n > "$WORK/anthocyanin_genes.bed"
echo "  located $(wc -l < "$WORK/anthocyanin_genes.bed") / $(wc -l < "$WORK/cand_ids.txt") candidate genes"

echo "############ 3. joint variant calling over candidate genes ############"
# RNA-aware: unique reads only (STAR uniques are MAPQ 255 >= 20); AD+DP kept for
# ASErtain's genotype calling; biallelic SNPs only.
bcftools mpileup -f "$REF" -R "$WORK/anthocyanin_genes.bed" -b "$WORK/bamlist.txt" \
        -a AD,DP --min-MQ 20 --min-BQ 20 --max-depth 5000 -Ou \
  | bcftools call -mv -Ou \
  | bcftools norm -f "$REF" -Ou \
  | bcftools view -v snps -m2 -M2 -Oz -o "$WORK/anthocyanin.vcf.gz"
bcftools index -t "$WORK/anthocyanin.vcf.gz"
echo "  VCF samples : $(bcftools query -l "$WORK/anthocyanin.vcf.gz" | tr '\n' ' ')"
echo "  biallelic SNPs: $(bcftools view -H "$WORK/anthocyanin.vcf.gz" | wc -l)"

echo "############ 4. ASErtain pipeline (report-mode bias) ############"
asertain run --config "$CONFIG" \
    --vcf "$WORK/anthocyanin.vcf.gz" --out "$WORK/ase" \
    --bias-mode report \
    --min-parent-depth 10 --maf-threshold 0.10 \
    --min-count-depth 10 --min-mapq 20 --min-baseq 20

echo
echo "Done. Key outputs in $WORK/:"
echo "  ase.informative_snps.tsv   phased informative SNPs per F1 genome"
echo "  ase.allele_counts.tsv      per-flower allele counts (+ variable_is_ref)"
echo "  ase.gene_ase.tsv           per-gene ASE calls (check fixed_allele_seen!)"
echo "  ase.report.html            summary"
echo
echo "Reference is kunthii, so first inspect the bias: in ase.allele_counts.tsv,"
echo "compare variable_count vs fixed_count where variable_is_ref=True vs False."
echo "If the kunthii(reference) allele is systematically inflated, re-run with the"
echo "N-masked reference:  asertain mask-reference --config $CONFIG \\"
echo "    --snps $WORK/ase.informative_snps.tsv --out-fasta $WORK/ref.Nmasked.fa \\"
echo "    --wasp-dir $WORK/wasp_snps   (then re-align F1 reads, count --bias-mode nmask)"
