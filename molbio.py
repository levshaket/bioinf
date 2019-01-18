#! python2.7
#programs for molecular biology

def motifs(sequence,Framesize):
# look into regex options 
#produce list of repeating motifs in sequence
	matches = []
	for framesize in range(Framesize,0,-1):
		for i in range(len(sequence)-framesize+1):
			query = sequence[i:i+framesize]
			if query in sequence[i+1:]:
				matches.append(query)
				a = True
				for match in matches[:-1]:
					if query in match:
						a = False
						break
				if a:
					print('*Framesize: %d* *Repeats found for: %s*' %(framesize,query))

def blat(q,s):
	import re
	x=re.compile(q)
	mo=re.search(x,s)
	return(mo.span())

def blast(q,s=None,f=None):
	from math import log
	if f:
		f=open('/genomes/{}'.format(f)); s=f.read().lower(); f.close();i=0; e=20000
	r=[];u=[q];t=[];b=log(len(s),4)
	while u:
		q=u[0];del u[0]
		try:
			r.append(blat(q,s))
			if f:
				i+=1
				if i==1:
					c=max(r[0][0]-e,0)
					s=s[c:min(r[0][1]+e,len(s))]
					b=log(len(s),4); r[0]=(r[0][0]-c,r[0][1]-c)
		except AttributeError:
			m=len(q)/2
			if m>b:
				u.extend((q[0:m],q[m:]))
			else:
				t.append(q)
			continue
	return([c,r])
###pblastn: possibly regex-based protein sequence search where each amino acids is broken into list of possible codons,
#Need dict of amino acid: codon list to generate regex string ((?:ata|atc|att)(?:gaa|gac)(?:ttt|ttc|ttg|tta)...etc)

###################################primer design#############################################
def primers(template,temp,minlen=18,maxlen=30,minsize=500,maxsize=1000):
	import regex
	primer_dict={}
	for primer in regex.finditer(primer_regex_string(min_gc(temp,minlen,maxlen)),template,overlapped=True):
		len_primer = len(primer.group(0))
		if len_primer in primer_dict.keys():
			primer_dict[len_primer].append(primer.span())
		else:
			primer_dict[len_primer]=[primer.span()]
	primer_indices = []
	for item in sorted(primer_dict.items()):
		primer_indices.append(item[1])
	for i,j in scan_order(primer_indices):
		for k in range(len(primer_indices[j])):
			product_size = primer_indices[j][-k][1]-primer_indices[i][0][0]
			if (product_size<=maxsize and product_size>=minsize):
				a,b = primer_indices[i][0]; c,d= primer_indices[j][-k]
				return([template[a:b],template[c:d],len(template[a:d]),template[a:d]])
		if i!=j:
			for k in range(len(primer_indices[i])):
				product_size = primer_indices[i][-k][1]-primer_indices[j][0][0]
				if (product_size<=maxsize and product_size>=minsize):
					a,b = primer_indices[j][0]; c,d = primer_indices[i][-k]
					return([template[a:b],rc(template[c:d]),len(template[a:d]),template[a:d]])

def scan_order(list_):
	from num_theory import partition_generator
	scan_order = []; maxrow = len(list_)-1
	for i in range(2*maxrow+1):
		partitions=[]
		for partition in partition_generator(i,2,0):
			a,b = partition
			if (a<maxrow and b<maxrow):
				partitions.append(partition)
		for reversed_partition in sorted(partitions,reverse=True):
			scan_order.append(reversed_partition)
	return(scan_order)

def primer_regex_string(min_gc):
	regex_list=[]
	for item in sorted(min_gc.items(),reverse=True):
		regex_list.append('(?:[gc]{%i}){%i<=i<=%i}'%(item[0],min(item[1])-item[0],max(item[1])-item[0]))
	regex_string= '('+'|'.join(regex_list)+')'
	return(regex_string)

def min_gc(temp,minlen,maxlen):#minlen,maxlen set to defaults in primers()
	out={}
	for i in range(minlen,maxlen+1):
		for j in range(0,i+1):
			if (41*j-672.4)/i >= temp-64.9:
				if j in out.keys():
					out[j].append(i)
				else:
					out[j]=[i]
				break
	return(out)

def Tm(primer):
	primer = primer.lower()
	gc=0
	for base in primer:
		if (base=='c' or base=='g'):
			gc+=1
	at=len(primer)-gc
	Tm = 64.9 +41*(gc-16.4)/(gc+at)
	return(Tm)
#########################################################################################################

####################sequence_editing###################
def rc(sequence):
	sequence = sequence[::-1]
	dic={'a':'t','t':'a','g':'c','c':'g'}
	rc=[]
	for base in sequence:
		rc.append(dic[base])
	rc = ''.join(rc)
	return(rc)
