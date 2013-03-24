import vcf
import vcf.model

d400 = vcf.Reader(open("../400Danes_exome_common.vcf.gz", 'r'))
d400af = vcf.Writer(open("../400Danes_exom_common_af.vcf", 'w'), d400)

for variant in d400:

    db_snp = 'rs%s' % variant.INFO['DBSNP'][0] if 'DBSNP' in variant.INFO else None
    print variant

    if variant.samples:
        record = vcf.model._Record(variant.CHROM, variant.POS, db_snp, variant.REF, variant.ALT, variant.QUAL, variant.FILTER,
                         {'AF': variant.aaf}, variant.FORMAT, [], [])

        d400af.write_record(record)

d400af.close()