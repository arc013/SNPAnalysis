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
#
# it gets all the positions from the VCF file
#

#
# to run python parseClinvar.py samples/HG00097.vcf clinvar_with_exac.tsv
#
def snpedia_nothing(soup):
    div = soup.find("div", attrs= { "class" : "noarticletext mw-content-ltr"})
    return div

def getrsIDfromVCFFile(filename):
    chrNums = {}
    for line in open(filename):
        if '#' not in line:
            info = line.split('\t')
            chrNum = info[0]
            rsID = info[2]
            if chrNum not in chrNums:
                chrNums[chrNum] = [rsID]
            else:
                chrNums[chrNum].append(rsID)
    #print chrNums.keys()
    return chrNums
        
	
if __name__ == "__main__":
    vcf_file = sys.argv[1]
    rsIDs_by_chr = getrsIDfromVCFFile(vcf_file)
    for geneId in rsIDs_by_chr['22']:
        print 'http://www.snpedia.com/index.php/'+geneId
        r = urllib.urlopen('http://www.snpedia.com/index.php/'+geneId).read()
        
        soup = BeautifulSoup(r, "html.parser")
        check = snpedia_nothing(soup)
        if check != None:
            continue 
        genos = getGenos(soup)
        #writeFile("genes/"+geneId+"_genos.csv", ["geno", "mag", "summary"], genos)
        pubs = getPubs(soup, True)
        #writeFile("genes/"+geneId+"_pubs.csv", "wb", "utf-8", ["id", "title"], pubs)
        print genos
        print genos
        #break 
    