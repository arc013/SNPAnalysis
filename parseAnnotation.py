import sys 
import numpy as np
import csv
from snpedia import getGenos, getPubs 
from bs4 import BeautifulSoup
import urllib
import json
import sys
import csv
import codecs
import operator 
#
# it gets all the positions from the VCF file
#

#
# to run python parseClinvar.py samples/HG00097.vcf clinvar_with_exac.tsv
#
##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: 'Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO' ">

#Annotation_Impacts = {}
#Gene_INFOs = {}
#Feature_Types = {}
#Transcript_BioTypes = {}
#Regions = {}

def populationAverage(line, multiple=False):
    infos = line.split(';')
    
    pops = {}
    if multiple == True:
        pops1 = {}
        pops2 = {}
    
    for info in infos:
        check = info.split('=')
        
        if len(check) == 1:
            continue 
        pop = check[0]
        #print check
        if pop in ['EAS_AF','AMR_AF','AFR_AF','EUR_AF','SAS_AF']:
            if multiple == False: 
                freqs = check[1].split(',')
                if len(freqs) == 1:
                    freq = float(check[1])
                else:
                    freq = [float(x) for x in freqs]
                pops[pop] = freq
            else: 
                freqs = check[1].split(',')
                pops1[pop] = float(freqs[0])
                pops2[pop] = float(freqs[1])
    if multiple == True:
        return (pops1, pops2)
    return pops 

def inputkeywords(filename, order_level):
    keywords = {}
    for line in open(filename):
        #print line 
        order = line.split('. ')
        keyword = order[1].strip()
        #print keyword 
        #print order
        level = int(order[0])
        if level <= order_level:
            keywords[keyword] = level
        
    return keywords 

def filterAnnotation(ann):
    #print len(ann.split('|'))
    info = filter(lambda x: x not in '', ann.split('|'))
    for x in info:
        check = x.split('&')
        if len(check) > 1:
            info.append(check[0])
            info.append(check[1])
    #print len(info)
    return info 
    """
    Allele = info[0]
    Annotation_Impact = info[1] 
    Gene_INFO = info[2]
    Gene_Name = info[3]
    Gene_ID = info[4]
    Feature_Type = info[5]
    Feature_ID = info[6]
    Transcript_BioType = info[7]
    if len(info) > 14:
        Rank = info[9]
        Region = info[11]
        Distance_IDs = info[14]
    
    if Annotation_Impact not in Annotation_Impacts:
        Annotation_Impacts[Annotation_Impact] = True
    if Gene_INFO not in Gene_INFOs:
        Gene_INFOs[Gene_INFO] = True
    if Feature_Type not in Feature_Types:
        Feature_Types[Feature_Type] = True
    if Transcript_BioType not in Transcript_BioTypes:
        Transcript_BioTypes[Transcript_BioType] = True
    if len(info) > 14:
        if Region not in Regions:
            Regions[Region] = True
    """
    #print info 

def snpedia_nothing(soup):
    div = soup.find("div", attrs= { "class" : "noarticletext mw-content-ltr"})
    return div

def getrsIDfromVCFFile(filename, file_keywords, order_level=46, check_keyword='', every=False, get_pop=False):
    keywords = inputkeywords(file_keywords, order_level)
    #print keywords 
    chrNums = {}
    #print every 
    for line in open(filename):
        if '#' not in line:
            info = line.split('\t')
            chrNum = info[0]
            rsID = info[2]
            #there are multiple rsID for same position and chromosome 
            multi_check = rsID.split(';')
            multiple = False
            #print multi_check
            if len(multi_check) > 1:
                multiple = True 
                multiple_rsID = multi_check

            ann = info[7]
            info = filterAnnotation(ann)
            #with open('test.txt', 'a') as text:
            #    text.write('\t'.join(info[1:-1]))
            #print info[1:10]
            keyword = ''
            order = 46
            if every == False:
                #print 'not here'
                if check_keyword == '':
                    for word in info:
                        if word in keywords:
                            keyword = word 
                            order = keywords[word]
                            #print order
                            continue
                else:
                    #print 'here'
                    #print check_keyword
                    if check_keyword in info:
                        #print 'yes'
                        keyword = check_keyword
                        order = keywords[check_keyword]
                    else: 
                        continue 
            else:
                #print 'here'
                for word in info:
                    if word in keywords:
                        keyword = word 
                        order = keywords[word]
                        #print order
                        continue
            
             
            if get_pop == True:   
                pops = populationAverage(info[0], multiple)
            
            if multiple == False:
                alist = [rsID, keyword]
                if get_pop == True:
                    alist.append(pops)
            else:
                alist1 = [multiple_rsID[0], keyword]
                alist2 = [multiple_rsID[1], keyword]
                if get_pop == True:
                    alist1.append(pops[0])
                    alist2.append(pops[1])
            
            if order <= order_level:
                if chrNum not in chrNums:
                    if multiple == False:
                        chrNums[chrNum] = [alist]
                    else:
                        chrNums[chrNum] = [alist1, alist2]
                else:
                    if multiple == False:
                        chrNums[chrNum].append(alist)
                    else:
                        chrNums[chrNum].append(alist1)
                        chrNums[chrNum].append(alist2)
    #print chrNums.keys()
    return chrNums
        
	
if __name__ == "__main__":
    vcf_file = sys.argv[1]
    keywords_file = sys.argv[2] 
    #chrNums = getrsIDfromVCFFile(vcf_file, keywords_file, 10)
    #chrNums = getrsIDfromVCFFile(vcf_file, keywords_file, order_level=10, get_pop=True)
    #chrNums = getrsIDfromVCFFile(vcf_file, keywords_file, every=True)
    chrNums = getrsIDfromVCFFile(vcf_file, keywords_file)
    if len(chrNums.keys()) != 0:
        print chrNums['22']
    """
    print Annotation_Impacts.keys()
    print Gene_INFOs.keys()
    print Feature_Types.keys()
    print Transcript_BioTypes.keys()
    print Regions.keys()
    """
    # sorting the genes based on the occurance
    #sorted_x = sorted(keywords.items(), key=operator.itemgetter(1), reverse=True)
    #print sorted_x[0:100]