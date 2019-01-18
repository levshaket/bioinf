from datetime import date
import pandas as pd#import openpyxl as xl
import numpy as np
import molbio as mb
import os, re, pyperclip

##preparing dataframe for transcripts
os.chdir('/home/neo/Desktop/nbscraper'); filenames = os.listdir('.');
descriptors = ['Locus','Description','Cds','Cds_Len','Sb_Id','Sb_Dna','Forward','Reverse','Product_Size','Product_Expected']
df = pd.DataFrame(index=filenames,columns=descriptors)
ew = pd.ExcelWriter('/home/neo/Desktop/primers.xlsx')

template_regex = re.compile('[atcgATCG]{100,}'); cds_regex = re.compile(r'name="(LOC\d*)?" directionality="(1|2)" translationMW="[0-9.]+" type="CDS".*?Segment range="(\d+)-(\d+).*?text="(.*?)"')

##filename loop
for filename in filenames:
	try:
		file_= open(filename); filecontents=file_.read(); file_.close()
		template = template_regex.search(filecontents).group()
		name,directionality,cds_start,cds_end,text = cds_regex.search(filecontents).groups()
		if directionality=='1':
			cds=template[int(cds_start)-1:int(cds_end)]
		elif directionality=='2':
			cds=mb.rc(template[int(cds_start)-1:int(cds_end)])
		cds_len=len(cds)
		df.loc[filename]['Locus':'Cds_Len']=[name,text,cds,cds_len]
		df.loc[filename]['Forward':'Product_Expected']=mb.primers(cds,60,maxlen=34)
#		pyperclip.copy('>%s_cds_%i\n%s\n'%(filename,cds_len,cds))#copies cds of modeled transcript
#		queryEntry.send_keys(Keys.CONTROL+'v')
	except AttributeError:
		continue
df.index.name='Name'; df=df.sort_values(['Cds_Len'],ascending=0)
df.to_excel(ew,sheet_name=str(date.today())); ew.save()
















