



def main(methyl_file,out_file):
	out = open(out_file, wb)
	for line in open(methyl_file):
		### if bismark header skip
		if line[0] == "B": continue
		name, m, seqid, pos, mtype = line.strip("\n").split("\t")
		### if not methylated skip
		if m == "-" : continue
		#### if scaffold skip
		if "uper" in seqid : continue
		new_line = "{0}\t{1}\t{2}\t{3}\n".format(seqid.strip("chromosome_"), int(pos)-1,int(pos)+1, name)
		out.write(new_line)
	out.close()

def jgi_format(methyl_file,out_file):
	out = open(out_file, "wb")
	for line in open(methyl_file):
		### if bismark header skip
		seqid, start, score = line.strip("\n").split("\t")[:3]
		### if not methylated skip
		#### if scaffold skip
		if "uper" in seqid : continue
		if float(score) < .4: continue
		new_line = "{0}\t{1}\t{2}\t{3}\n".format(seqid.strip("chromosome_"),start,int(start)+1, score)
		out.write(new_line)
	out.close()


#jgi_format("../jgi_raw/NUXC_CHH_raw.txt","../jgi_raw/NUXC_CHH_raw_filtered.bed")
#jgi_format("../jgi_raw/NUXC_CHG_raw.txt","../jgi_raw/NUXC_CHG_raw_filtered.bed")
jgi_format("../jgi_raw/NUXC_CpG_raw.txt","../jgi_raw/NUXC_CpG_raw_filtered.bed")





