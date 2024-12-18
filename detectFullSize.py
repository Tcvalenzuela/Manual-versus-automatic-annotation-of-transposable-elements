from Bio import SeqIO
import argparse
parser = argparse.ArgumentParser(description= "Select Kimura under param and create bed")
parser.add_argument("--S", "-SeqConse",  help="Sequence consensus")
parser.add_argument("--T", "-tableInsertion",  help="Out of RepeatMasker table")
parser.add_argument("--Ot", "-OutputTable",  help="Output table with the detailed information for each insertion")
parser.add_argument("--St", "-StatsTable",  help="Output Stats")

arg = parser.parse_args()
ListSeq={}
Insertions={}
Divergence={}
ListNames=[]
TEorganizer = {}
counter = {}

consensus = SeqIO.parse(open(arg.S),"fasta")
for fasta in consensus:
	name = fasta.id
	seq = fasta.seq
	leng = len(seq)
	ListSeq[name] = str(leng)
	ListNames.append(name)

Reps = open(arg.T)
for insertion in Reps:
	if insertion.startswith("Scaff"):
		continue
	stripped=insertion.split("\t")
	Scaff=stripped[0]
	start=stripped[1]
	end=stripped[2]
	size=stripped[3]
	TEtype=stripped[4].split("#")[0]
	Kval=stripped[5]
	if float(ListSeq[TEtype]) * 0.9 <=  float(size) <= float(ListSeq[TEtype])+(float(ListSeq[TEtype]) * 0.1):
		if float (Kval) < 20:
			#print (str(ListSeq[TEtype])+"\t"+insertion.rstrip("\n"))
			TEorganizer[insertion] = TEtype
			if TEtype in counter:
				counter[TEtype]=float(counter[TEtype])+1
			else:
				counter[TEtype]= 0


oWriter = open(arg.Ot)

for name in ListNames:
	value={i for i in TEorganizer if TEorganizer[i]==name}
	for ele in value:
		owriter.write(ele)
oWriter.close()

oWriter = open(arg.St)
for i in counter:
	oWriter.write(i+"\t"+counter[i]+"\n")
