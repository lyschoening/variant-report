from __future__ import division
import argparse
import math
import subprocess
import numpy
import os.path
import pysam
from annotation.coverage_template import get_template
from annotation.refgene import RefGene
from os.path import basename

__author__ = 'lyschoening'



def main():

    parser = argparse.ArgumentParser(description='Generate a PDF report of Coverage in genes.')
    parser.add_argument('genes', metavar='gene', type=str, nargs='*', help='RefSeq gene name(s)')
    parser.add_argument('-c' '--coverage', type=str, help='Bam file (*.bam), indexed')
    parser.add_argument('-G', type=str, help='File listing gene accession numbers, one per line')
    parser.add_argument('-r', '--refgene', type=str, help='RefGene (genePredExt) table, e.g. -r hg19.refGene')
    parser.add_argument('-o', '--output', type=str, default='Report')
    # parser.add_argument('--tables', type=str, default='IMPORTANT_VARIANTS,ALL_VARIANTS')

    args = parser.parse_args()

    print args

    refgene = RefGene(args.refgene)

    print args
    genes = set(args.genes)

    if args.G:
        for line in open(args.G, 'r'):
            genes.add(line.strip().split("\t")[0])

    genes = list(refgene.filter(genes))

    print "Genes:", ", ".join(map(str, genes))

    sample = pysam.Samfile(args.c__coverage, 'rb')
    sample_name = basename(args.c__coverage).split('.')[0]


    tex_file_prefix = '%s_%s_coverage' % (args.output, sample_name)
    tex_file_name = os.path.abspath('%s.tex' % tex_file_prefix)

    template = get_template()

    with open(tex_file_name, 'w') as tex_file:


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


                    pileups = numpy.zeros(gene.end - gene.start)
                    coverage = numpy.zeros(len(gene.exons))

                   # short_chrom = gene.chrom[3:] if gene.chrom.startswith('chr') else gene.chrom

                    for pileup in sample.pileup(gene.chrom, gene.start, gene.end):
                        try:
                            pileups[pileup.pos - gene.start] = pileup.n
                        except IndexError:
                            pass


                    def exons():
                        for j, exon in enumerate(gene.exons):
                            exon_start, exon_end = exon
                            exon_pileup = pileups[max(0, exon_start - gene.start):min(gene.end - gene.start, exon_end - gene.start)]


                            if gene.is_reverse:
                                points = [(i, exon_pileup[len(exon_pileup) // (20 - 20 * i) if i != 20 else 0]) for i in range(21)]
                            else:
                                points = [(i, exon_pileup[len(exon_pileup) // 20 * i if i != 0 else 0]) for i in range(21)]

                            #coverage[i][j] = max(exon_pileup)
                            # TODO reverse points if on minus strand
                            # TODO exon.coding False if non-coding
                            # TODO define exon numbers along with genes.



                            yield {'points': points, 'name': exon_names[j], 'coding': exon_start > gene.coding_start and exon_end < gene.coding_end}





                    #variants = list(vcf_reader.fetch(gene.chrom, gene.start, gene.end))


                    yield (gene, list(exons()), max(20, pileups.max()))
                except Exception as e:
                    print e
                    print "Error..", gene


        tex_file.write(template.render(objects=objects(genes), sample=sample, sample_name=sample_name))


    print ('pdflatex', '-output-directory=%s' % os.path.dirname(tex_file_name), tex_file_name)

    for i in range(2): # call twice for proper table layout.
        subprocess.call(('pdflatex', '-output-directory=%s' % os.path.dirname(tex_file_name), tex_file_name))

    try:
        subprocess.call(('latexmk', '-output-directory=%s' % os.path.dirname(tex_file_name), '-c', tex_file_name))
    except OSError:
        print "Could not find latexmk executable, skipping."

    for extension in ('tex', 'aux', 'log', 'out'): # latexmk old versions
        os.unlink('.'.join((tex_file_prefix, extension)))

if __name__ == '__main__':
    main()