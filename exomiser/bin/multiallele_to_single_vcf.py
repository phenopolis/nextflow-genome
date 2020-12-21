#! /usr/bin/env python
from __future__ import print_function
import sys
import argparse

#these 9 column headers are standard to all VCF files
#the samples headers depend on the number of samples in the file
#which we find out once we read the #CHROM line
SAMPLE_HEADERS=[]

parser=argparse.ArgumentParser(description='Arguments to multiallele_to_single_gvcf.py')
parser.add_argument('--GQ',default=None,type=int)
parser.add_argument('--DP',default=None,type=int)
parser.add_argument('--rowname_file', default=None)
parser.add_argument('--colname_file', default=None)
parser.add_argument('--calls_file', default=None)
parser.add_argument('--depth_file', default=None)
parser.add_argument('--file_for_vep', default=None)
parser.add_argument('--headers_for_vep',default=','.join(['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT']),type=str)
parser.add_argument('--suppress', action="store_true", default=False)
args=parser.parse_args()

# these files are for R
if args.rowname_file:
    ROWNAMES_FILE=open(args.rowname_file, 'w')
else:
    ROWNAMES_FILE=None
if args.colname_file:
    COLNAMES_FILE=open(args.colname_file, 'w')
else:
    COLNAMES_FILE=None
if args.calls_file:
    CALLS_FILE=open(args.calls_file, 'w')
else:
    CALLS_FILE=None
if args.depth_file:
    DEPTH_FILE=open(args.depth_file,'w')
else:
    DEPTH_FILE=None
# this file is for VEP
#args.GQ
#args.DP
#STD_HEADERS=['CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO','FORMAT']
STD_HEADERS=args.headers_for_vep.strip().split(',')
if args.file_for_vep:
    VEP_FILE=open(args.file_for_vep,'w')


'''
Potentially this could also normalise indel, but needs reference to reference.
Also when 2 variants are 1 bp apart they should be merged into one (this would need another pass).
this function cleans messy variant format
    var = '1-94512001-GTT-GAT'
    print clean_variant(var)
prints
    1-94512002-T-A
'''
def clean_variant(s):
    chrom=s['CHROM']
    pos=s['POS']
    ref=s['REF']
    alt=s['ALT']
    pos = int(pos)
    if len(ref) < len(alt):
        ran = range(len(ref))
    else:
        ran = range(len(alt))
    # insert
    for e in ran:
        ref_e = len(ref) - e - 1
        alt_e = len(alt) - e - 1
        if ref[ref_e] != alt[alt_e]: break
    for b in ran:
        if ref[b] != alt[b] or len(ref[b:ref_e+1]) == 1 or len(alt[b:alt_e+1]) == 1:
            break
    #return '-'.join([chrom,str(pos+b),ref[b:ref_e+1],alt[b:alt_e+1]])
    #return [chrom,str(pos+b),ref[b:ref_e+1],alt[b:alt_e+1]]
    s['CHROM']=chrom
    s['POS']=str(pos+b)
    s['REF']=ref[b:ref_e+1]
    s['ALT']=alt[b:alt_e+1]
    return s

def recode_gt(gt):
    return {'./.':'NA', './0':'0','./1':'1','./2':'2','0/0':'0', '0/1':'1', '1/1':'2', '0/2':'1'}[gt]

def print_header(s):
    if s.startswith('#CHROM') and COLNAMES_FILE:
        for sh in SAMPLE_HEADERS: print(sh,file=COLNAMES_FILE)
        COLNAMES_FILE.close()
    if s.startswith('#CHROM') and VEP_FILE:
        print('##fileformat=VCFv4.2',file=VEP_FILE)
        print('#'+'\t'.join(STD_HEADERS),file=VEP_FILE)
    if s.startswith('#') and not COLNAMES_FILE: print(s)

def print_line(s):
    #s=clean_variant(s)
    #s['ID']='-'.join([s['CHROM'],s['POS'],s['REF'][:4],s['ALT'][:4]])
    if VEP_FILE:
        print('\t'.join([s[h] for h in STD_HEADERS]),file=VEP_FILE)
    if ROWNAMES_FILE:
        print('-'.join([s['CHROM'],s['POS'],s['REF'],s['ALT']]), file=ROWNAMES_FILE)
    if CALLS_FILE:
        print( *([recode_gt(s[h].split(':')[0]) for h in SAMPLE_HEADERS]), sep='\t', file=CALLS_FILE)
    if DEPTH_FILE:
        print( *([int(s[h].split(':')[2] if s[h].split(':')[2]!='.' else 0) for h in SAMPLE_HEADERS]), file=DEPTH_FILE)
    if not CALLS_FILE and not ROWNAMES_FILE and not args.suppress:
        print( *([s[h] for h in STD_HEADERS+SAMPLE_HEADERS]), sep='\t' )

for line in sys.stdin:
    line=line.strip()
    #remove beginning
    #lines which start with '###'
    #are not tab separated
    if line.startswith('##'):
        print_header(line)
        continue
    #this is tab separated line
    s=line.split("\t")
    #header line yay!: #CHROM ...
    if line.startswith('#'):
        headers=s
        headers[0]=headers[0].strip('#')
        #the first 9 names in the header are standard (see above)
        if (headers[0:len(STD_HEADERS)] != STD_HEADERS):
            print( headers[0:len(STD_HEADERS)],  STD_HEADERS )
            raise 'hell'
        #everything else in the header is a sample name
        SAMPLE_HEADERS=headers[len(STD_HEADERS):]
        print_header(line)
        continue
    #you can now access elements by their header name
    s=dict(zip(headers,s))
    #split alternate alleles
    alternative_alleles=s['ALT'].split(',')
    #I would expect each sample to be formatted according to the format column.
    #However I found in practise this not true.
    #Only the first two fields: GT (genotype) and AD (allele depth)
    #are compulsory, the remaining fields, including DP (total depth),
    #can be missing.
    #So I calculate DP from AD instead.
    #I ignore the other fields for now.
    #Actually, if the genotype is missing then none of the fields are compulsory.
    #In fact AD can be '.' even when geno is 0/0 if using HaplotypeCaller with multiple BAMs.
    for h in SAMPLE_HEADERS:
        d=dict(zip(s['FORMAT'].split(':'),s[h].split(':')))
        GT=d['GT'].replace('|','/')
        AD=d.get('AD',','.join(['0']*(len(alternative_alleles)+1)))
        if 'DP' in d:
            DP=(d['DP'])
        else:
            DP=str(sum([int(x) if x!='.' else 0 for x in AD.split(',')]))
        GQ=d.get('GQ','')
        # if GQ is '.' the must be missing
        if GQ == '.' and GT != './.': raise 'hell'
        # if fails QC set to missing
        if ((args.DP and int(DP) < args.DP) or (args.GQ and GQ!='.' and int(GQ) < args.GQ)): GT='./.'
        PL=d.get('PL','')
        s[h]=':'.join([GT,AD,DP,GQ,PL])
    #if single alternate allele
    #then just print out as normal
    if len(alternative_alleles)<=1:
        s['FORMAT']='GT:AD:DP:GQ:PL'
        print_line(s)
        continue
    #otherwise split over as many lines as there are alternative alleles
    #each genotype then takes a 1 if matches the alternative or a 0 otherwise
    alleles=[s['REF']]+alternative_alleles
    n2geno=dict(zip(['.']+[str(i) for i in range(0,len(alleles))],['.']+alleles))
    #sys.stderr.write(line+'\n')
    for idx, alt in enumerate(alternative_alleles):
        s1=s.copy()
        s['ALT']=alt
        # split info fields on ',' across lines
        def f(info):
            k,v,=info.split('=')
            if ',' in v:
                return '='.join([k,v.split(',')[idx]])
            else:
                return '='.join([k,v])
        s1['INFO']=';'.join(map(f, [x for x in s['INFO'].split(';') if '=' in x]))
        #recode GT
        for h in SAMPLE_HEADERS:
            d=dict(zip(s1['FORMAT'].split(':'),s1[h].split(":")))
            #length of allele depth is 2 where first is always REF allele depth
            #and second can be either ALT
            if 'AD' in d and d['AD']!='' and d['AD']!='.':
                AD=d['AD'].split(',')
                AD[1]=AD[idx+1]
                AD=','.join(AD[:2])
            else:
                AD='0'
            # matches REF -> 0
            # matches ALT -> 1
            # matches anything else -> .
            #sorted so that "." precedes number
            GT='/'.join(sorted([{s['REF']:'0',s['ALT']:'1'}.get(n2geno[g],'.') for g in d['GT'].split('/')]))
            if AD!='.':
                DP=str(sum([int(x) for x in AD.split(',')]))
            else:
                DP=0
            s1[h]=':'.join( [GT, AD, DP] )
            s1['FORMAT']='GT:AD:DP'
        s1['ALT']=alt
        print_line(s1)
        #
        continue
        if alt=='*': alt='-'
        pos=int(s1['POS'])
        ref=s1['REF']
        info=process_info(s1['INFO'],idx)
        print_line(s1['#CHROM'],pos,id,ref,alt,qual,filter,info,sep='\t')
        if len(alt)>len(ref) and ref in alt :
            #insertion
            ind=alt.index(ref)
            pos=pos+ind+len(ref)-1
            alt=alt.replace(ref,'')
            ref='-'
            print_line(s1['#CHROM'],pos,id,ref,alt,qual,filter,info,sep='\t')
        elif len(alt)<len(ref) and alt in ref :
            #deletion
            ind=ref.index(alt)
            pos=pos+ind+len(alt)
            ref=ref.replace(alt,'')
            alt='-'
            print_line(s1['#CHROM'],pos,id,ref,alt,qual,filter,info,sep='\t')



if ROWNAMES_FILE: ROWNAMES_FILE.close()
if CALLS_FILE: CALLS_FILE.close()
if DEPTH_FILE: DEPTH_FILE.close()
if VEP_FILE: VEP_FILE.close()



