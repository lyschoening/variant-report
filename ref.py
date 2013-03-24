__author__ = 'lyschoening'


exons = ((41277287, 41277500), (41276033, 41276132), (41267742, 41267796), (41258472, 41258550),
         (41256884, 41256973), (41256138, 41256278), (41251791, 41251897), (41249260, 41249306),
         (41247862, 41247939), (41243451, 41246877), (41242960, 41243049), (41234420, 41234592),
         (41228504, 41228631), (41226347, 41226538), (41222944, 41223255), (41219624, 41219712),
         (41215890, 41215968), (41215349, 41215390), (41209068, 41209152), (41203079, 41203134),
         (41201137, 41201211), (41199659, 41199720), (41196311, 41197819))

gene_start = 41196311
gene_end = 41277500



thick_start = 41197694
thick_end = 41276113

pos = 41251931
pos =41197708
pos = 41208692

reverse_strand = True

def ref(pos):

    exon_offset = 0
    prev_exon_start = 0

    if pos < thick_start:
        return "before coding"
    if pos > thick_end:
        # todo get length of all coding
        return "after coding"


    for i, (exon_start, exon_end) in enumerate(exons, start=1):
        if exon_start <= pos <= exon_end:
            # position inside exon.
            return "in exon", i, exon_offset + (exon_end - pos - 1), None


        if pos < exon_start:
            print 'exon after', i
            prev_exon_start = exon_start
            exon_offset += exon_end - exon_start


        if pos > exon_end:
            if i == 1: # upstream first exon
                return 'before first exon'
            else:
                if prev_exon_start - pos < exon_end - pos:
                    # closer to previous exon:
                    return "after exon", i - 1
                else:
                    return "before exon", i, exon_offset + 2, exon_end - pos

def rref():

    coding_start = ref(thick_start)
    coding_end = ref(thick_end)

    print coding_start, coding_end
    coding_start = coding_start[2]
    coding_end = coding_end[2]

    description, exon, inner, outer = ref(pos)
    print description, exon, inner, outer
    print description, exon, inner - coding_end, outer

#
#
#        elif pos > exon_end:
#            # position before exon.
#            if prev_exon_start != 0:
#
#                if prev_exon_start - pos < pos - exon_end:
#                    return i - 1, exon_offset + 1, pos - prev_exon_start
#                else:
#                    return i, exon_offset + 2, exon_end - pos - 1
#            return 1, 0, pos - exon_start - 1
#        else: # position after exon
#            prev_exon_start = exon_start
#            exon_offset += exon_end - exon_start
#        return len(exons), exon_offset, pos - prev_exon_start
#
print rref()

total = 0

for i, exon in enumerate(exons, start=1):
    length = exon[1] - exon[0]
    print i, exon, length, '%s - %s' % (total, total + length)
    total += length