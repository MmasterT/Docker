from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Restriction import *

rb=
rb = DpnII + HinfI + DdeI + MseI
dpnii_dig=[]
for record in SeqIO.parse("GCA_905146935.1_idScaPyra1.1_genomic.fna", "fasta"):
	my_seq = Seq(str(record.seq))
	dpnii_dig.append(DpnII.catalyse(my_seq))

test_1=list(sum(dpnii_dig, ()))


hinfi_dig=[]
for i in range(len(test_1)):
	newSeq= test_1[i]
	hinfi_dig.append(HinfI.catalyse(my_seq))
print(hinfi_dig)

test_2=list(sum(hinfi_dig, ()))

ddei_dig=[]
for i in range(len(test_1)):
	newSeq= test_1[i]
	hinfi_dig.append(HinfI.catalyse(my_seq))
