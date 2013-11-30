class Methyl(object):
    
    def __init__(self,filename,mtype):
        self.filename = filename
        self.mtype = mtype

    def get_metyl(self):
        d = []
        for  line in open(self.filename):
            metyline = MethylLine(line)
            d.append(metyline)
        return d

class MethylLine(object):
    
    def __init__(self,line):
        args = line.strip().split("\t")
        self.seqid = args[0].strip("chromosome_")
        self.pos = args[1]
        self.per = float(args[2])
        self.name = "{0}_{1}".format(self.seqid,self.pos)

    def __str__(self):
        return "\t".join(self.dups)



