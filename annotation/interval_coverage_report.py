from __future__ import division
import argparse
import subprocess
import numpy
import os.path
import pysam
from annotation.utils import BEDReader
from annotation.interval_coverage_template import get_template
from os.path import basename

__author__ = 'lyschoening'


def main():
    parser = argparse.ArgumentParser(description='Generate a PDF report of Coverage in genes.')
    parser.add_argument('genes', metavar='gene', type=str, nargs='*', help='RefSeq gene name(s)')
    parser.add_argument('-c' '--coverage', type=str, help='Bam file (*.bam), indexed')
    parser.add_argument('-i', '--intervals', type=str, help='Bed file (*.bed)')
    parser.add_argument('-o', '--output', type=str, default='Report')
    parser.add_argument('-minimum-coverage', type=int, default=30)
    parser.add_argument('--use-chr-prefix', type=bool, default=False)

    # parser.add_argument('--tables', type=str, default='IMPORTANT_VARIANTS,ALL_VARIANTS')

    args = parser.parse_args()

    print args

    sample = pysam.Samfile(args.c__coverage, 'rb')
    sample_name = basename(args.c__coverage).split('.')[0]

    tex_file_prefix = args.output
    tex_file_name = os.path.abspath('%s.tex' % tex_file_prefix)

    template = get_template()

    # TODO logic to switch between 'chr1' and '1' naming.



    with open(tex_file_name, 'w') as tex_file:

        def objects(intervals):
            for line in BEDReader(intervals):
                chrom, start, end = line['chrom'], int(line['chromStart']), int(line['chromEnd'])

                pileups = numpy.zeros(end - start + 1)

                if args.use_chr_prefix:
                    chrom = chrom[3:] if chrom.startswith('chr') else chrom
                else:
                    chrom = chrom if chrom.startswith('chr') else 'chr{}'.format(chrom)

                for pileup in sample.pileup(chrom, start, end):
                    try:
                        pileups[pileup.pos - start] = pileup.n
                    except IndexError:
                        pass

                quants = numpy.percentile(pileups[1:], [10, 30, 50, 70, 90])

                yield line['name'], chrom, start, end, numpy.min(pileups[1:]), numpy.max(pileups), quants

        tex_file.write(
            template.render(objects=objects(args.intervals), minimum_coverage=args.minimum_coverage, sample=sample,
                            sample_name=sample_name))

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