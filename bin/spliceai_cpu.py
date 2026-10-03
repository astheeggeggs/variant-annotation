#!/usr/bin/env python
"""
SpliceAI on CPU, one record at a time. This is the non-batched loop of the pinned SpliceAI fork
(spliceai/__main__.py:run_spliceai, used when -B is 1), which at that revision crashes on a typo
(`args_output_data`). Same scoring code (spliceai.utils), header line and output as that loop.

  spliceai_cpu.py -I in.vcf -O out.vcf -R ref.fa -A annotation [-D 50] [-M 0]
"""

import argparse
import logging

import pysam
from spliceai.utils import Annotator, get_delta_scores

HEADER = ('##INFO=<ID=SpliceAI,Number=.,Type=String,Description="SpliceAIv1.3.1 variant '
          'annotation. These include delta scores (DS) and delta positions (DP) for '
          'acceptor gain (AG), acceptor loss (AL), donor gain (DG), and donor loss (DL). '
          'Format: ALLELE|SYMBOL|DS_AG|DS_AL|DS_DG|DS_DL|DP_AG|DP_AL|DP_DG|DP_DL">')


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-I", required=True, help="input VCF")
    p.add_argument("-O", required=True, help="output VCF")
    p.add_argument("-R", required=True, help="reference FASTA")
    p.add_argument("-A", required=True, help="gene annotation file")
    p.add_argument("-D", type=int, default=50, help="max distance to gained/lost splice site (default 50)")
    p.add_argument("-M", type=int, default=0, choices=[0, 1], help="mask scores (default 0)")
    args = p.parse_args()
    logging.basicConfig(format="%(asctime)s %(levelname)s %(name)s: - %(message)s",
                        datefmt="%Y-%m-%d %H:%M:%S", level=logging.INFO)

    ann = Annotator(args.R, args.A)
    vcf = pysam.VariantFile(args.I)
    vcf.header.add_line(HEADER)
    out = pysam.VariantFile(args.O, mode="w", header=vcf.header)
    for record in vcf:
        scores = get_delta_scores(record, ann, args.D, args.M)
        if len(scores) > 0:
            record.info["SpliceAI"] = scores
        out.write(record)
    vcf.close()
    out.close()


if __name__ == "__main__":
    main()
