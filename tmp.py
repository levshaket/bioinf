#! python2.7
#temporary programs in development

	
def primers_manual(template,mintemp,minspacer):#minframe,maxspacer,maxframe #not necessary to finish but may do so to compare efficiency of program to the regex one to be pasted here presently
	framesize=18#minframe
	rctemplate = rc(template)
	while framesize<=40:#maxframe
		i=0
		while i<len(template)-minspacer:
			if Tm(template[i:i+framesize])>mintemp:
				j=1000#maxspacer
				while j>=minspacer:
					1000
				
				try:
					pos = i+j
					if template[pos-framesize:position]:
						j-=1
				except IndexError:
					j-=1
			i+=1
		framesize+=1

############################## regex_search_bestmatch based primers ##############################
def primers(template,temp,minlen=18,maxlen=30,minspacer=500,maxspacer=1000):
	import regex
	primer = regex.search(primer_regex_string(min_gc(temp,minlen,maxlen),minspacer,maxspacer),template,overlapped=True)
	return([primer.group(1),primer.group(2),len(primer.group(0)),primer.group(0)])

def primer_regex_string(min_gc,minspacer,maxspacer):#minspacer,maxspacer set to defaults in primers()
	regex_list=[]
	for item in sorted(min_gc.items(),reverse=True):
		regex_list.append('(?:[gc]{%i}){%i<=i<=%i}'%(item[0],min(item[1])-item[0],max(item[1])-item[0]))
	regex_string= '('+'|'.join(regex_list)+')'
	regex_string= '(?b)'+regex_string+('.{%i,%i}?'%(minspacer,maxspacer))+regex_string#(?b) is a best-match flag supported by python-regex demanding minimal mismatches
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
###################################################################################################################





def blat(query,sequence):
	import re
	rgx = re.compile(query)
	bld,red,end=('\033[1m','\033[91m','\033[0m')
	flnk=6
#	print('\nQuery: %s\nSequence: %s\n'%(query,sequence))
	for match in rgx.finditer(sequence):
		a,b = match.span()
		c,d=sequence[a-flnk:a],sequence[b:b+flnk]
		print('Positions %d..%d: %s%s%s%s%s%s'%(a+1,b,c,bld,red,query,end,d))

def fastblast(Query,Subject,minMatch=0.99,visualize=True):#to be improved by allowing query overhangs or block mapping
	import regex,os
	Query=Query.lower(); Subject=Subject.lower()
	Query_Len, Sub_Len = len(Query),len(Subject)
	#	score=[0,0,(0,0),(0,0)]#for stepping down minMatch values in new function
	for i in xrange(Query_Len):
		for j in xrange(Query_Len,0,-1):
			if j>i:
				query_len=j-i
				query=Query[i:j]
				errors = int((1-minMatch)*query_len)
				rgx = regex.compile('(?:%s){s<=%i}'%(query,errors))
				mo = regex.search(rgx,Subject)
				try:
					subject = mo.group()
				except AttributeError:
					continue
				n = 0
				for k in range(query_len):
					if query[k]==subject[k]:
						n+=1
				p=n/float(query_len);s=n*p
				score=[s,p,(i,j),mo.span()]
				break
		else:
			continue
		break
	s,p,(qi,qj),(si,sj)=score
	if visualize:
		query=Query[qi:qj]; subject=Subject[si:sj]
		matched_len=len(subject); pipes=[]
		for i in range(matched_len):
			if query[i] == subject[i]:
				pipes.append('|')
			else:
				pipes.append(' ')
		pipes=''.join(pipes)

		c=int(os.popen('stty size','r').read().split()[1])
		print('Query_Len={}\tSubject_Len={}\tMatched={}\tIdentity={:3.2f}\tScore={}'.format(Query_Len,Sub_Len,matched_len,p,s))
		if c-24<matched_len:
			for i in range(0,matched_len-c,c):
				print('QUERY\t{}\t{}\t{}\n\t\t{}\nSUBJECT\t{}\t{}\t{}\n\n'.format(qi+i,query[i:i+c],qi+i+c-1,pipes[i:i+c],si+i,subject[i:i+c],si+i+c-1))
			print('QUERY\t{}\t{}\t{}\n\t\t{}\nSUBJECT\t{}\t{}\t{}\n'.format(qi+i+c,query[i+c:],qj-1,pipes[i+c:],si+i+c,subject[i+c:],sj-1))
		else:
			print('QUERY\t{}\t{}\t{}\n\t\t{}\nSUBJECT\t{}\t{}\t{}\n'.format(qi,query,qj-1,pipes,si,subject,sj-1))
	return(score)

