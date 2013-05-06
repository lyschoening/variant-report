# -*- coding: utf-8 -*-

import argparse
from operator import attrgetter
import re
import subprocess
import math
import vcf
from annotation.refgene import RefGene, VariantDescription
from annotation.template import get_template
import os

__author__ = 'lyschoening'


def main():
    parser = argparse.ArgumentParser(description='Generate a PDF report of Variants in genes.')
    parser.add_argument('genes', metavar='gene', type=str, nargs='*', help='RefSeq gene name(s)')
    parser.add_argument('-v', '--variants', type=str, help='Variant file (VCF), indexed')
    parser.add_argument('-G', type=str, help='File listing gene accession numbers, one per line')
    parser.add_argument('-r', '--refgene', type=str, help='RefGene (genePredExt) table, e.g. -r hg19.refGene')
    parser.add_argument('-o', '--output', type=str, default='Report')
    args = parser.parse_args()

    refgene = RefGene(args.refgene)

    print args
    genes = set(args.genes)

    if args.G:
        for line in open(args.G, 'r'):
            genes.add(line.strip().split("\t")[0])

    genes = list(refgene.filter(genes))

    print "Genes:", ", ".join(map(str, genes))

    vcf_reader = vcf.Reader(open(args.variants, 'r'))

    first_record = vcf_reader.next()

    samples = map(attrgetter('sample'), first_record.samples)

    for sample, sample_name in enumerate(samples, start=0):

        tex_file_prefix = '%s_%s' % (args.output, sample_name)
        tex_file_name = '%s.tex' % tex_file_prefix

        template = get_template()

        with open(tex_file_name, 'w') as tex_file:

            for gene in genes:
                try:
                    print gene
                    for variant in vcf_reader.fetch(gene.chrom, gene.start, gene.end):
                        print variant, variant.samples
                    print
                except:
                    print "Error..", gene



            def objects(genes):
                for gene in genes:
                    try:
                        coverage_table_rows = int(math.ceil(len(gene.exons) / 15.0))
                        coverage_tuples = []

                        exon_names = gene.get_exon_names()

                        for row in range(coverage_table_rows):
                            row_exon_names = exon_names[15 * row:15 * (row + 1)]
                            row_exons = gene.exons[15 * row: 15 * (row + 1)]
                            #row_exons_coverage = None
                            # TODO coverage
                            coverage_tuples.append(zip(row_exon_names, [0] * len(row_exons)))

                        variants = vcf_reader.fetch(gene.chrom, gene.start, gene.end)

                        yield (gene, lambda: (VariantDescription(v, gene) for v in variants), coverage_tuples)
                    except:
                        print "Error..", gene


            tex_file.write(template.render(objects=objects(genes), sample=sample, sample_name=sample_name))

        for i in range(2): # call twice for proper table layout.
            subprocess.call(('pdflatex', '-output-directory=%s' % os.path.dirname(tex_file_name), tex_file_name))
        subprocess.call(('latexmk', '-output-directory=%s' % os.path.dirname(tex_file_name), '-c', tex_file_name))

        for extension in ('tex', 'aux', 'log', 'out'): # latexmk old versions
            os.unlink('.'.join((tex_file_prefix, extension)))

if __name__ == '__main__':
    main()