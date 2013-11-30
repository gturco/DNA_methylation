from flatfeature import Bed
from jgi_metyl import Methyl
from intersection import Intersecter, Feature





def loadintointersect(bed_file):
        """loads features of bedfile into a tree for finding their freq within genomic regions"""
        query_list = {}
        feature_list = Bed(bed_file)
        for feature in feature_list:
	    if float(feature['accn']) < .4 : continue
	    if feature['seqid'] not in list(query_list):
                query_list[feature['seqid']] = Intersecter()
            query_list[feature['seqid']].add_interval(Feature(int(feature['start']),
            int(feature['end']),name=['name']))
        print "done"
        return query_list




def get_pos(gene, i):
    three_prime = []
    three_ds = []
    five_prime = []
    five_ds = []
    exonic = []
    eds = []
    intronic = []
    ids = []
    exons = gene['locs']
    introns = [(e,exons[n+1][0]) for n,(s,e) in enumerate(exons) if n != len(exons)-1]
    introns = [(gene['start'],exons[0][0])] + introns + [(exons[-1][1],gene["end"])]
    for methyl in i:
	methyl_pos = methyl[0]
        if methyl_pos < exons[0][0]:
            three_d = exons[0][0] - methyl_pos
            three_prime.append(methyl)
            three_ds.append(three_d)
        elif methyl_pos > exons[-1][1]:
            five_d = methyl_pos  = exons[-1][1]
            five_prime.append(methyl)
            five_ds.append(five_d)
        elif methyl_pos > exons[0][0] and methyl_pos < exons[-1][1]:
            e = [(methyl,n) for n,ex in enumerate(exons) if methyl_pos >= ex[0] and methyl_pos <= ex[1]]
            if len(e) == 0:
		intronic.append(methyl_pos)
	    else:
            	exonic.append(e[0][0])
            	eds.append(e[0][1])
	    #if len(three_ds) + len(five_ds)+ len(eds) != len(i):
            #	intron = [(methyl,n) for n,intr in enumerate(introns) if methyl_pos >= intr[0] and methyl_pos <= intr[1]]
            # 	if len(intron) == 0 : 
	    #		print "not finding hit"
	    #		continue
	    #	intronic.append(intron[0][0])
            #	ids.append(intron[0][1])
    return three_prime,map(str,three_ds), five_prime,map(str,five_ds), exonic,map(str,eds), intronic, ids


def main(filename,gene_file,mtype,padding,outfile):
    out = open(outfile,"wb")
    genome = Bed(gene_file)
    methyl_points = loadintointersect(filename)
    for gene in genome:
	start = gene['start'] - padding
        end = gene['end'] + padding
        matches = methyl_points[gene['seqid']].find(start, end)
        i = [(m.start,m.stop) for m in matches]
        three_p,three_d,five_p,five_d,exon,eds, intronic, ids = get_pos(gene, i)
        if len(i) > 0:
		newline = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\n".format(gene['accn'],len(i),len(three_p),",".join(three_d),len(five_p),",".join(five_d),len(exon),",".join(eds),len(intronic),",".join(ids))
        	out.write(newline)
    out.close()

### test
##main("../NUXC_CHG_raw_head.bed","../NUXC_CHG_raw_head.bed","CpG",1000,"../NUXC_CHG_raw_head.genelist")
#main("../sorg.test.bed","../sorg.test.bed","CpG",1000,"../sorg_head.genelist")
main("../test_methyl.bed","../test.bed","CpG",10000,"../test.genelist")
#
#main("../NUXC_CpG_raw.bed","../sorg.bed","CpG",1000,"../NUXC_CpG_raw.genelist")
#main("../NUXC_CHG_raw.bed","../sorg.bed","CHG",1000,"../NUXC_CHG_raw.genelist")
### prefiltered out 0.4
#main("../NUXC_CHH_raw_filtered.bed","../sorg.bed","CHH",1000,"../NUXC_CHH_raw.genelist")





