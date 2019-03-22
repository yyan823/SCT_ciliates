
import argparse, os, sys
from distutils import spawn

from Bio import SeqIO, AlignIO
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio.Data.CodonTable import CodonTable
from collections import Counter

#------------------------------ special codon table ------------------------------#

fake_table = CodonTable(forward_table={
	'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
	'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
	'TAT': 'Y', 'TAC': 'Y', 'TAA': '#', 'TAG': '+',           
	'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
	'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
	'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
	'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
	'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
	'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
	'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
	'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
	'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
	'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
	'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
	'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
	'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'},
	start_codons = [ 'ATG'],
	stop_codons = [])

#------------------------------ Colors For Print Statements ------------------------------#

class color:
   PURPLE = '\033[95m'
   CYAN = '\033[96m'
   DARKCYAN = '\033[36m'
   ORANGE = '\033[38;5;214m'
   BLUE = '\033[94m'
   GREEN = '\033[92m'
   YELLOW = '\033[93m'
   RED = '\033[91m'
   BOLD = '\033[1m'
   UNDERLINE = '\033[4m'
   END = '\033[0m'


def check_usearch_path():

	usearch_path = ''

	if usearch_path == '':
		usearch_path = spawn.find_executable("usearch")
	else:
		pass

	if usearch_path == None:
		print (color.BOLD + '\n\nPlease open this script and check that you have included'\
		+' the PATH to the'+color.BLUE+' "usearch" '+color.END+color.BOLD+'executable.\n\n'+color.END)
		print (color.BOLD+color.BLUE+'LOOK FOR:\n\n'+color.RED\
		+'#------------------------------ UPDATE USEARCH PATH BELOW! -------------------------------#'\
		+color.BLUE+'\n\nThis is somewhere around lines 50 - 80...\n\n'+color.END)

		sys.exit()
	else:
		pass
	
	return usearch_path


def translate(infile):
	infasta = [i for i in SeqIO.parse(infile,'fasta')]
	with open(infile.replace('.fasta','_universal_ORF.aa.fasta'),'w+') as w:
		for i in infasta:
			if 'TAG' in str(i.seq)[-3:] or 'TAA' in str(i.seq)[-3:] or 'TGA' in str(i.seq)[-3:]:
				new_seq = Seq(str(i.seq)[:-3],IUPAC.unambiguous_dna)
				w.write('>'+i.id+'\n'+str(new_seq.translate(table=fake_table)).upper()+'\n')
			else:
				new_seq = Seq(str(i.seq),IUPAC.unambiguous_dna)
				w.write('>'+i.id+'\n'+str(new_seq.translate(table=fake_table)).upper()+'\n')


def ublast(usearch_path,infile):
	os.system(usearch_path+' -ublast '+infile.replace('.fasta','_universal_ORF.aa.fasta') + ' -db ../db_StopFreq/RepEukProts.udb -evalue 1e-20 -target_cov 0.5 -maxhits 100 -blast6out '+infile.replace('.fasta','_universal_ORF.RepEukProts.tsv')+' -dbmatched '+infile.replace('.fasta','_dbmatched.fasta'))


def toAln(infile):
	in_tsv = [i for i in open(infile.replace('.fasta','_universal_ORF.RepEukProts.tsv'),'r').read().split('\n') if i != '']
	in_dbmatched = [i for i in SeqIO.parse(infile.replace('.fasta','_dbmatched.fasta'),'fasta')]
	in_fasta = [i for i in SeqIO.parse(infile.replace('.fasta','_universal_ORF.aa.fasta'),'fasta')]	
	grab = {}
	for line in in_tsv:
		grab[line.split('\t')[0]] = ''
	for line in in_tsv:
		grab[line.split('\t')[0]] = grab[line.split('\t')[0]] + '\t' + line.split('\t')[1]


	outlist = []
	for key, value in grab.items():
		outlist.append(key + ''.join(grab[key]))
	for i in outlist:
		toaln = []
		for j in in_fasta:
			if j.id in i.split('\t'):
				toaln.append('>'+j.id+'\n'+str(j.seq)+'\n')
			else:
				pass
		for h in in_dbmatched:
			if h.id in i.split('\t'):
				toaln.append('>'+h.id+'\n'+str(h.seq)+'\n')
			else:
				pass
		with open(i.split('\t')[0] + '.fasta', 'w+') as w:
			for x in toaln:
				w.write(x)
				
def Aln(infile):
	for file in os.listdir('.'):
		if file.startswith(infile.split('_')[0]) and '_OG5_' in file:
			os.system('mafft --anysymbol --inputorder '+file+' > '+file.split('_Len')[0]+'_mafft.fasta')
	
def GrabCol(infile,aalist):

### count TGA ###
	for file in os.listdir('.'):
		if file.endswith('_mafft.fasta') and file.startswith(infile.split('_')[0]):
			tableout = open(infile.split('_InTree')[0]+'_StopReassign_TGA.tsv','a')
			inmft = AlignIO.read(file,'fasta')
			ColCount = 0
			outCol = []
			outlist = []
			
			for SeqRecord in inmft:
				if SeqRecord.id.startswith(infile.split('_')[0]):
					position = [idx for idx,w in enumerate(str(SeqRecord.seq)) if w == '*']
					ColCount = ColCount + len(position)
#					print SeqRecord.id + '\t' + str(position)
					for column in position:
						k = inmft[:, column]
						if int(str(Counter(k).most_common(1)).split(',')[1].split(')')[0])>float(len(k)*0.5):
							for i in str(k):
								if i != str(Counter(k).most_common(1)).split(',')[0].split('[(')[1].replace("'","") and i != '*':
									k = k.replace(i,'-')
								else:
									pass
							outCol.append(k)
						else:
							pass
							
			inlist = [i for i in aalist if i != '']
			inlist = sorted(inlist)
			outCol2 = str(outCol).replace("[","").replace("]","").replace("', '","").replace("'","")
			codonCount = []
			for i in inlist:
				codonCount.append(outCol2.count(i))	
#			print str(sum(codonCount))					
#---- write out counting of each seq						
			if open(infile.split('_InTree')[0]+'_StopReassign_TGA.tsv','r').readline().split('\t')[0] == 'Contig':
				tableout.write(file.split('_mafft.fasta')[0]+'\t'+ str(codonCount).replace('[','').replace(']','').replace(',','\t') + '\t' + str(sum(codonCount)) + '\n')	
			else:
				tableout.write('Contig\t'+ str(inlist).replace('[','').replace(']','').replace('\n','').replace("'","").replace(',','\t')+'\tcodon counted\n')
				tableout.write(file.split('_mafft.fasta')[0]+'\t'+ str(codonCount).replace('[','').replace(']','').replace(',','\t') + '\t' + str(sum(codonCount)) + '\n')	

### count TAA ###
	for file in os.listdir('.'):
		if file.endswith('_mafft.fasta') and file.startswith(infile.split('_')[0]):
			tableout = open(infile.split('_InTree')[0]+'_StopReassign_TAA.tsv','a')
			inmft = AlignIO.read(file,'fasta')
			ColCount = 0
			outCol = []
			outlist = []
			
			for SeqRecord in inmft:
				if SeqRecord.id.startswith(infile.split('_')[0]):
					position = [idx for idx,w in enumerate(str(SeqRecord.seq)) if w == '#']
					ColCount = ColCount + len(position)
#					print SeqRecord.id + '\t' + str(position)
					for column in position:
						k = inmft[:, column]
						if int(str(Counter(k).most_common(1)).split(',')[1].split(')')[0])>float(len(k)*0.75):
							for i in str(k):
								if i != str(Counter(k).most_common(1)).split(',')[0].split('[(')[1].replace("'","") and i != '#':
									k = k.replace(i,'-')
								else:
									pass
							outCol.append(k)
						else:
							pass
							
			inlist = [i for i in aalist if i != '']
			inlist = sorted(inlist)
			outCol2 = str(outCol).replace("[","").replace("]","").replace("', '","").replace("'","")
			codonCount = []
			for i in inlist:
				codonCount.append(outCol2.count(i))						
#---- write out counting of each seq						
			if open(infile.split('_InTree')[0]+'_StopReassign_TAA.tsv','r').readline().split('\t')[0] == 'Contig':
				tableout.write(file.split('_mafft.fasta')[0]+'\t'+ str(codonCount).replace('[','').replace(']','').replace(',','\t') + '\t' + str(sum(codonCount)) + '\n')	
			else:
				tableout.write('Contig\t'+ str(inlist).replace('[','').replace(']','').replace('\n','').replace("'","").replace(',','\t')+'\tcodon counted\n')
				tableout.write(file.split('_mafft.fasta')[0]+'\t'+ str(codonCount).replace('[','').replace(']','').replace(',','\t') + '\t' + str(sum(codonCount)) + '\n')	

### count TAG ###
	for file in os.listdir('.'):
		if file.endswith('_mafft.fasta') and file.startswith(infile.split('_')[0]):
			tableout = open(infile.split('_InTree')[0]+'_StopReassign_TAG.tsv','a')
			inmft = AlignIO.read(file,'fasta')
			ColCount = 0
			outCol = []
			outlist = []
			
			for SeqRecord in inmft:
				if SeqRecord.id.startswith(infile.split('_')[0]):
					position = [idx for idx,w in enumerate(str(SeqRecord.seq)) if w == '+']
					ColCount = ColCount + len(position)
#					print SeqRecord.id + '\t' + str(position)
					for column in position:
						k = inmft[:, column]
						if int(str(Counter(k).most_common(1)).split(',')[1].split(')')[0])>float(len(k)*0.75):
							for i in str(k):
								if i != str(Counter(k).most_common(1)).split(',')[0].split('[(')[1].replace("'","") and i != '+':
									k = k.replace(i,'-')
								else:
									pass
							outCol.append(k)
						else:
							pass
							
			inlist = [i for i in aalist if i != '']
			inlist = sorted(inlist)
			outCol2 = str(outCol).replace("[","").replace("]","").replace("', '","").replace("'","")
			codonCount = []
			for i in inlist:
				codonCount.append(outCol2.count(i))						
### write out counting of each seq						
			if open(infile.split('_InTree')[0]+'_StopReassign_TAG.tsv','r').readline().split('\t')[0] == 'Contig':
				tableout.write(file.split('_mafft.fasta')[0]+'\t'+ str(codonCount).replace('[','').replace(']','').replace(',','\t') + '\t' + str(sum(codonCount)) + '\n')	
			else:
				tableout.write('Contig\t'+ str(inlist).replace('[','').replace(']','').replace('\n','').replace("'","").replace(',','\t')+'\tcodon counted\n')
				tableout.write(file.split('_mafft.fasta')[0]+'\t'+ str(codonCount).replace('[','').replace(']','').replace(',','\t') + '\t' + str(sum(codonCount)) + '\n')	
		
		
def cleanup(infile):
	taxon_folder = infile.split('_InTree')[0] + '/'
	orig_folder = infile.split('_InTree')[0] + '/Original_files'
	result_folder = infile.split('_InTree')[0] + '/Result'
	if os.path.isdir(taxon_folder) != True:
		os.system('mkdir '+taxon_folder)
	if os.path.isdir(orig_folder) != True:
		os.system('mkdir '+orig_folder)
	if os.path.isdir(result_folder) != True:
		os.system('mkdir '+result_folder)
		
	for file in os.listdir('.'):
		if file.startswith(infile.split('_InTree')[0]):
			if 'StopReassign' in file or '_ORF' in file or 'dbmatched' in file:
				os.system('mv '+file+' '+infile.split('_InTree')[0] + '/Result')			
			elif file.endswith('NTD.ORF.fasta'):
				os.system('mv '+file+' '+infile.split('_InTree')[0] + '/Original_files')
			elif '_OG5_' in file or 'mafft' in file:
				os.system('rm '+file)
			else:
				pass
	
def main():
	usearch_path = check_usearch_path()
	aalist = [i for i in open('aalist.txt','r').read().split('\n')]
	for infile in os.listdir('.'):
#		if infile.endswith('XContClean.fasta'):
		if infile.endswith('NTD.ORF.fasta'):
			translate(infile)
			ublast(usearch_path,infile)
			toAln(infile)
			Aln(infile)
			GrabCol(infile,aalist)
			cleanup(infile)
main()
