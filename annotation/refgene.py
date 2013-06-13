import csv
import string
import math
import re

"""

Nomenclature from http://www.hgvs.org/mutnomen/recs.html

"""

__author__ = 'lyschoening'


class RefGene(object):
    def __init__(self, filename):
        self.__filename = filename


    def __iter__(self):
        with open(self.__filename, 'rb') as file:
            for row in csv.reader(file, delimiter='\t', quotechar='"'):
                unknown_1, accession, chrom, strand, chrom_start, chrom_end, thick_start, thick_end, block_count, block_starts, block_ends, unknown_5, name, unknown_2, unknown_3, unknown_4 = row

                block_starts = map(int, block_starts[:-1].split(","))
                block_ends = map(int, block_ends[:-1].split(","))
                exons = zip(block_starts, block_ends)

                if strand == '-':
                    exons.reverse()

                yield Gene(name, chrom,
                           int(chrom_start),
                           int(chrom_end),
                           exons,
                           strand=strand,
                           coding_start=int(thick_start),
                           coding_end=int(thick_end), refseq_name=accession)

    def filter(self, gene_names):
        for gene in self:
            if gene.accession in gene_names:
                yield gene


EFF_REGEX_PATTERN = r'(.*)\(([^\|]*)\|([^\|]*)\|([^\|]*)\|([^\|]*)\|([^\|]*)\|([^\|]*)\|([^\|]*)\|([^\|]*)\|([^\|]*)\|(.*)\)'

EFF_TYPE_TEXTS = {
    'UTR_5_PRIME': "UTR 5'",
    'UTR_3_PRIME': "UTR 3'",
    'UTR_5_DELETED': "UTR 5' Deleted",
    'UTR_3_DELETED': "UTR 3' Deleted",
    'NON_SYNONYMOUS_CODING': "Nonsynonymous Coding",
}

AMINO_ACID_SYMBOLS = {
    'A': ('Ala', 'Alanine'),
    'R': ('Arg', 'Arginine'),
    'N': ('Asn', 'Asparagine'),
    'D': ('Asp', 'Aspartic acid'),
    'B': ('Asx', 'Asn or Asp'),
    'C': ('Cys', 'Cysteine'),
    'Q': ('Gln', 'Glutamine'),
    'E': ('Glu', 'Glutamic acid'),
    'Z': ('Glx', 'Gln or Glu'),
    'G': ('Gly', 'Glycine'),
    'H': ('His', 'Histidine'),
    'I': ('Ile', 'Isoleucine'),
    'L': ('Leu', 'Leucine'),
    'K': ('Lys', 'Lysine'),
    'M': ('Met', 'Methionine'),
    'F': ('Phe', 'Phenylalanine'),
    'P': ('Pro', 'Proline'),
    'S': ('Ser', 'Serine'),
    'T': ('Thr', 'Threonine'),
    'W': ('Trp', 'Tryptophan'),
    'Y': ('Tyr', 'Tyrosine'),
    'V': ('Val', 'Valine'),
}

def aa_symbol_to_name(sym):
    return AMINO_ACID_SYMBOLS[sym][2]


class Eff(object):
    _eff_matcher = re.compile(EFF_REGEX_PATTERN)

    def __init__(self, effect_str):
        result = Eff._eff_matcher.match(effect_str)

        print result.groups()

        self.effect_type, self.impact, self.functional_class, self.codon_change, self.aa_change, _unknown_number_, self.gene_name,\
        self.gene_biotype, self.coding, self.transcript_id, self.exon_id = result.groups()

        # ('INTRON', 'MODIFIER', '', '', '', '403', 'PTEN', '', 'CODING', 'NM_000314', '6')
        # UTR_5_PRIME(MODIFIER||||403|PTEN||CODING|NM_000314|1)

        if self.effect_type in EFF_TYPE_TEXTS:
            self.effect_type_text = EFF_TYPE_TEXTS[self.effect_type]
        else:
            self.effect_type_text = string.capwords(self.effect_type.replace("_", " "))

        self.aa_change_text = None

        if self.aa_change:
            self.aa_change = (self.aa_change[0], int(self.aa_change[1:-1]), self.aa_change[-1])

            pos = self.aa_change[1]
            ref, sub = tuple(
                AMINO_ACID_SYMBOLS[aa][0] if aa in AMINO_ACID_SYMBOLS else aa
                    for aa in (self.aa_change[0], self.aa_change[2])
            )

            self.aa_change_text = 'p.%s%s%s' % (ref, pos, sub)


dna_complement_trans = string.maketrans('TAGCtagc', 'ATCGATCG')

class VariantDescription(object):
    def __init__(self, variant, gene):
        self.variant = variant
        self.gene = gene

        if 'EFF' in variant.INFO:
            self.effects = gene.parse_effs(variant.INFO['EFF'])
        else:
            self.effects = []

        self.impacts = [effect.impact for effect in self.effects]
        self.aa_change_texts = set([effect.aa_change_text for effect in self.effects if effect.aa_change_text])

        self.exon, self.exon_offset = gene.get_position_details(variant.POS)
        self.pos = "%s:%s" % (variant.CHROM, variant.POS)

        self.type = variant.var_type
        self.REF = variant.REF # TODO implement get fwd.
        self.ID = variant.ID
        self.count = sum(call.is_variant == True for call in variant.samples)


    def __getattr__(self, item):
        return getattr(self.variant, item)

    def get_call(self, sample):
        return self.variant.samples[sample]

    def get_base_change(self, call):
        return self.gene.get_variant_alleles(call)

    def get_mut_ref(self, call):
        return self.gene.get_variant_description(call)

    def get_genotype_call_accuracy(self, call):
        quality = call.data.GQ
        nines = math.ceil(quality / 10.0)
        accuracy = (1. - 10. ** (-quality / 10.)) * 100
        if accuracy > 99.999:
            return '>99.999%'
        else:
            return ("{0:.%if}%%" % nines).format(accuracy)


class Gene(object):
    def __init__(self, name, chrom, start, end, exons, strand, coding_start=None, coding_end=None, refseq_name=None):
        self.name = name
        self.accession = refseq_name
        self.chrom = chrom
        self.start = start
        self.end = end
        self.exons = tuple(exons)
        exons.reverse()
        self.exons_reverse = tuple(exons)
        self.coding_start = coding_start
        self.coding_end = coding_end
        self.is_reverse = True if strand == '-' else False

        self.coding_start_offset = self.__get_base_offset(self.coding_start)
        self.coding_end_offset = self.__get_base_offset(self.coding_end)

    def get_exon_names(self):
        return range(1, len(self.exons) + 1)

    def get_named_exons(self):
        return zip(self.get_exon_names(), self.exons)


    def parse_effs(self, eff_str_list, transcript_only=True):
        effs = map(Eff, eff_str_list)

        if transcript_only:
            effs = filter(lambda eff: eff.transcript_id == self.accession, effs)
        return effs

    def __get_base_offset(self, pos):
        exon_offset, prev_exon_end, prev_exon_start = 1, 0, 0

        if self.is_reverse:
            for i, (exon_start, exon_end) in enumerate(self.exons, start=1):
                if exon_start <= pos <= exon_end: # position inside exon.
                    return i, exon_offset + (exon_end - pos + 1), None
                elif pos > exon_end:
                    if i == 1: # position before first exon.
                        return 1, 0, pos - exon_start - 1
                    else: # position before Nth exon; N > 1
                    #                        if pos == 41208692:
                    #                            print pos, exon_start, exon_end, pos, exon_start - pos, prev_exon_start - pos, exon_end - pos
                    #                            exit()
                        if prev_exon_start - pos < pos - exon_end:
                            return i - 1, exon_offset, prev_exon_start - pos + 1
                        else:
                            return i, exon_offset + 2, exon_end - pos - 1
                else:
                    prev_exon_start = exon_start
                    exon_offset += exon_end - exon_start
            return len(self.exons), exon_offset, pos - prev_exon_end
        else:
            for i, (exon_start, exon_end) in enumerate(self.exons, start=1):
                if exon_start <= pos <= exon_end:
                    # position inside exon.
                    return i, exon_offset + (pos - exon_start + 1), None
                elif pos < exon_start:
                    # position before exon.
                    if prev_exon_end != 0:
                        if pos - prev_exon_end < exon_start - pos:
                            return i - 1, exon_offset + 1, pos - prev_exon_end
                        else:
                            return i, exon_offset + 2, pos - exon_start - 1
                    return 1, 0, pos - exon_start - 1
                prev_exon_end = exon_end
                exon_offset += exon_end - exon_start
            return len(self.exons), exon_offset, pos - prev_exon_end

    def get_coding_regions(self):
        for i, (exon_start, exon_end) in enumerate(self.exons):
            if exon_start >= self.coding_start:
                if exon_end <= self.coding_end:
                    yield (exon_start, exon_end, i)
                elif exon_start <= self.coding_end:
                    yield (exon_start, self.coding_end, i)
                    # else: coding ends before exon starts
            else: # coding starts after exon starts
                if exon_end >= self.coding_start:
                    if exon_end <= self.coding_end:
                        yield (self.coding_start, exon_end, i)
                    else:
                        yield (self.coding_start, self.coding_end, i)
                        # else: coding starts after exon ends


    def get_position_details(self, position):
        exon, exon_offset, intron_offset = self.__get_base_offset(position)

        if intron_offset is None:
            return exon, str(exon_offset)
        if intron_offset < 0:
            return exon, str(intron_offset).replace('-', '--')
        else:
            return exon, '+' + str(intron_offset)

    def get_variant_alleles(self, call):
        site = call.site
        return list(set(bases for bases in (str(site.alleles[int(a)]) for a in call.gt_alleles) if bases != site.REF))

    def get_variant_description(self, call):
        site = call.site

        assert call.called
        # TODO FIX for reverse strand genes.

        def str_for_pos(position):
            exon, exon_offset, intron_offset = self.__get_base_offset(position)

            if self.is_reverse:
                coding_exon_offset = self.coding_end_offset[1] - 1
            else:
                coding_exon_offset = self.coding_start_offset[1]

            if exon_offset < coding_exon_offset:
                coding_offset = exon_offset - coding_exon_offset - 1
            else:
                coding_offset = exon_offset - coding_exon_offset

            if position > self.coding_end and not self.is_reverse:
                return '*%s' % (exon_offset - self.coding_end_offset[1])
            if position < self.coding_start and self.is_reverse:
                return '*%s' % (exon_offset - self.coding_end_offset[1])

            if exon_offset == 0:
                return '%s' % intron_offset
            if intron_offset is None:
                return coding_offset

            if intron_offset < 0:
                return '%s%s' % (coding_offset, intron_offset)
            else:
                return '%s+%s' % (coding_offset, intron_offset)

        def str_for_range(a, b, change):
            if self.is_reverse:
                a, b = b, a
            return 'c.%s_%s%s' % (str_for_pos(a), str_for_pos(b), change)

        reference_sequence, allele_sequences = site.REF, (str(site.alleles[int(a)]) for a in call.gt_alleles)
        allele_sequences = [sequence for sequence in allele_sequences if sequence != reference_sequence]

        if self.is_reverse:
            reference_sequence = reference_sequence.translate(dna_complement_trans)
            allele_sequences = [sequence.translate(dna_complement_trans) for sequence in allele_sequences]

        allele_sequences = list(set(allele_sequences))

        if site.is_snp:
            return 'c.%s%s>%s' % (str_for_pos(site.POS), reference_sequence, ', '.join(map(str, allele_sequences)))
        elif site.is_indel:
            if not allele_sequences:
                deletion_length = len(site.REF)
            else:
                deletion_length = len(site.REF) - len(allele_sequences[0])

            if deletion_length > 0:
                if deletion_length == 1:
                    return 'c.%sdel' % str_for_pos(site.POS + 1)
                else:
                    return str_for_range(site.POS + deletion_length - 1, site.POS + deletion_length * 2 - 2,
                                         'del%i' % deletion_length)
            else: # is insertion
                return str_for_range(site.POS, site.POS + 1, 'ins' + allele_sequences[0][len(reference_sequence):])

        # TODO insertions

        elif site.is_sv:
            if site.INFO['SVTYPE'] in ('DUP', 'INV'):
                variation_length = len(site.REF)
                if variation_length == 1:
                    return 'c.%sdup' % str_for_pos(site.POS)
                else:
                    return str_for_range(site.POS, site.POS + variation_length, site.INFO['SVTYPE'].lower())
            return 'Unknown: %s' % site.INFO['SVTYPE']

            # TODO delins

    def __repr__(self):
        return "<%s %s %s %s..%s, exons: %s>" % (
        self.name, self.accession, self.chrom, self.start, self.end, len(self.exons))