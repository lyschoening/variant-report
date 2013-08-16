from __future__ import division
import argparse
import math
from operator import attrgetter, itemgetter
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

    tex_file_prefix = args.output
    tex_file_name = os.path.abspath('%s.tex' % tex_file_prefix)

    template = get_template()

    with open(tex_file_name, 'w') as tex_file:
        def objects(genes):
            for gene in genes:
                try:
                    COLS_PER_ROW = 20
                    coverage_table_rows = int(math.ceil(len(gene.exons) / COLS_PER_ROW))
                    coverage_tuples = []

                    coverage_exon_rows = []

                    exon_names = gene.get_exon_names()

                    for row in range(coverage_table_rows):
                        row_exon_names = exon_names[COLS_PER_ROW * row:COLS_PER_ROW * (row + 1)]
                        row_exons = gene.exons[COLS_PER_ROW * row: COLS_PER_ROW * (row + 1)]
                        #row_exons_coverage = None
                        # TODO coverage
                        coverage_exon_rows.append(zip(range(COLS_PER_ROW * row, min(len(gene.exons), COLS_PER_ROW * (row + 1))), row_exon_names))
                        coverage_tuples.append(zip(row_exon_names, [0] * len(row_exons)))

                    pileups = numpy.zeros(gene.end - gene.start)
                    coverage = numpy.zeros(len(gene.exons))

                    # short_chrom = gene.chrom[3:] if gene.chrom.startswith('chr') else gene.chrom

                    print gene.chrom, gene.start, gene.end

                    for pileup in sample.pileup(gene.chrom, gene.start, gene.end):
                        try:
                            pileups[pileup.pos - gene.start] = pileup.n
                        except IndexError:
                            pass

                    print 'x', sum(pileups), sample


                    def exons():

                        print gene.name,

                        for j, exon in enumerate(gene.exons):
                            exon_start, exon_end = exon
                            exon_pileup = pileups[max(0, exon_start - gene.start):min(gene.end - gene.start,
                                                                                      exon_end - gene.start)]


                            all_otr = exon_start > gene.coding_start and exon_end < gene.coding_end

                            print exon_names[j], all_otr, exon_end <= gene.coding_start, exon_start >= gene.coding_end, sum(exon_pileup), sum(pileups)

                            if not all_otr:
                                if exon_end <= gene.coding_start or exon_start >= gene.coding_end:
                                    yield { 'name': '--', 'skip': True}
                                    continue



                                if exon_start < gene.coding_start:
                                    print 'start', exon_start, gene.coding_start, gene.coding_start - exon_start + 1
                                    exon_pileup = exon_pileup[gene.coding_start - exon_start + 1:]
                                if exon_end > gene.coding_end:
                                    print 'end', exon_end, gene.coding_end, -(gene.coding_end - exon_end)
                                    exon_pileup = exon_pileup[:gene.coding_end - exon_end]

                            if gene.is_reverse:
                                points = [(i, exon_pileup[((len(exon_pileup) - 1) // 20) * (20 - i) if i != 20 else 0]) for i in
                                          range(21)]
                            else:
                                points = [(i, exon_pileup[((len(exon_pileup) - 1) // 20) * i if i != 0 else 0]) for i in
                                          range(21)]

                            #coverage[i][j] = max(exon_pileup)
                            # TODO reverse points if on minus strand
                            # TODO exon.coding False if non-coding
                            # TODO define exon numbers along with genes.

                            yield {
                                'points': points,
                                'name': exon_names[j],
                                'min': exon_pileup.min(),
                                'max': exon_pileup.max(),
                                'mean': exon_pileup.mean(),
                                'skip': False,
                                'coding': exon_start > gene.coding_start and exon_end < gene.coding_end
                            }


                    #variants = list(vcf_reader.fetch(gene.chrom, gene.start, gene.end))


                    exons_lst = list(exons())
                    yield (gene, exons_lst, len(exons_lst) - sum(map(itemgetter('skip'), exons_lst)), max(20, pileups.max()), coverage_exon_rows)
                except Exception as e:
                    print e
                    print "Error..", gene

        gene_objs = objects(genes)

        tex_file.write(template.render(objects=gene_objs, sample=sample, sample_name=sample_name))


    print ('pdflatex', '-output-directory=%s' % os.path.dirname(tex_file_name), tex_file_name)


    for i in range(2): # call twice for proper table layout.
        subprocess.call(('pdflatex', '-output-directory=%s' % os.path.dirname(tex_file_name), tex_file_name))

    exit()

    try:
        subprocess.call(('latexmk', '-output-directory=%s' % os.path.dirname(tex_file_name), '-c', tex_file_name))
    except OSError:
        print "Could not find latexmk executable, skipping."

#    for extension in ('tex', 'aux', 'log', 'out'): # latexmk old versions
#        os.unlink('.'.join((tex_file_prefix, extension)))

if __name__ == '__main__':
    main()