#
# Author: Normand Overney
#
    
import sys
import sys
from numbers import Number
from collections import Set, Mapping, deque
import operator

#
# The function below was grabbed from stack overflow to verfiy the size
# of the dictionary in bytes that was to be stored in memory
#
try: # Python 2
    zero_depth_bases = (basestring, Number, xrange, bytearray)
    iteritems = 'iteritems'
except NameError: # Python 3
    zero_depth_bases = (str, bytes, Number, range, bytearray)
    iteritems = 'items'

def getsize(obj_0):
    """Recursively iterate to sum size of object & members."""
    def inner(obj, _seen_ids = set()):
        obj_id = id(obj)
        if obj_id in _seen_ids:
            return 0
        _seen_ids.add(obj_id)
        size = sys.getsizeof(obj)
        if isinstance(obj, zero_depth_bases):
            pass # bypass remaining control flow and return
        elif isinstance(obj, (tuple, list, Set, deque)):
            size += sum(inner(i) for i in obj)
        elif isinstance(obj, Mapping) or hasattr(obj, iteritems):
            size += sum(inner(k) + inner(v) for k, v in getattr(obj, iteritems)())
        # Now assume custom object instances
        elif hasattr(obj, '__slots__'): 
            size += sum(inner(getattr(obj, s)) for s in obj.__slots__ if hasattr(obj, s))
        else: 
            attr = getattr(obj, '__dict__', None)
            if attr is not None:
                size += inner(attr)
        return size
    return inner(obj_0)

#
# given the OMIM datafile it parses through the file with the chromosme number
# as the key and the values being the unique sorted lists of start and end with
# another dictionary with the start position and end position as a the keys and
# the values are the Gene symbols and the gene name
#
def parseOMIMFile(filename):
    dictOMIM = {}
    chromDict = {}
    prevChrom =""
    index = 0
    for line in open(filename):
        #if index == 10:
        #    break
        #index += 1
        # the infomation that is deemed important 
        if "#" not in line:
            info = line.split("\t")
            chromNum = info[0]
            startPos = int(info[1])
            endPos = int(info[2])
            geneSym = info[6]
            geneName = info[7]
            

            #laoding the dictioanry 
            if chromNum not in dictOMIM:
                chromDict = {}
                chromDict[(startPos,endPos)] = [(geneSym, geneName)]
                dictOMIM[chromNum] = [set([(startPos, endPos)]), chromDict] 
            else:
                if (startPos, endPos) not in chromDict:
                    chromDict[(startPos, endPos)] = [(geneSym, geneName)]
                else:
                    chromDict[(startPos, endPos)].append((geneSym, geneName))
                dictOMIM[chromNum][0].add((startPos, endPos))
     
                dictOMIM[chromNum][1] = chromDict
    #sort the keys 
    for key in dictOMIM:
        items = dictOMIM[key]
        sortedStartEnds = list(sorted(items[0]))
        startlist, endlist = zip(*sortedStartEnds)
        dictOMIM[key][0] = list(endlist)
        dictOMIM[key].insert(0,list(startlist))

    #print dictOMIM
    #print getsize(dictOMIM)
    return dictOMIM

#
# creates a gene ID to phenotype dictionary with the OMIM database as the
# bases for it 
#
def createGeneIDtoPheno(filename):
    GeneIDtoPheno = {}
    for line in open(filename):
        if "#" not in line:
            info = line.split("\t")
            geneSym = info[6]
            # makes sure the information is actually in the column
            if len(info) > 12 and info[12] != "":
                pheno = info[12]
                if geneSym not in GeneIDtoPheno:
                    GeneIDtoPheno[geneSym] = [pheno]
                else:
                    GeneIDtoPheno[geneSym].append(pheno)
    return GeneIDtoPheno

#
# using the dictioanry of gene IDs to PHenotype it goes through the sorted
# dictionary for the phenotype and the value for how many times you want to beat
# it
#
def findPhenotypeFromOMIM(dictGeneIDPheno, sortedGenes, value):
    relaventPheno = []
    for gene in sortedGenes:
        geneID = gene[0]
        if len(relaventPheno) >= value:
            break
        if geneID in dictGeneIDPheno:
            relaventPheno.append((dictGeneIDPheno[geneID], gene[1]))
    return relaventPheno
        
#
# a simple modified binary search algorithm that instead of just returning
# a macth to the data it will return the closest smaller value to the item
# when looking for it in the list
#
def modifiedBinarySearch(alist, item):
    first = 0
    last = len(alist) - 1
    found = False
    #print 'item:',item
    while first <= last and not found:
        midpoint = (first + last)//2
        if alist[midpoint] == item:
             found = True
        else:
            if item < alist[midpoint]:
                last = midpoint - 1
            else:
                first = midpoint+1
		#print alist[midpoint]
    #print alist[midpoint-1]
    #making sure the return references and index that cannot be references
    if midpoint == 0:
        return 0
    return midpoint-1

#
# the actual running of the algorithm which takes in the dict of OMIM 
# where it will count occurances for each range of positions which it will
# verify with a simple binary search. It then returns a list of sorted genes
# with the gene IDs as the key and the value being the number of occurances
#
def binarySearchOMIM(dictOMIM, data):
    chroms = {}
    for chrNum in data:
        #print chrNum
        #print dictOMIM.keys()
        listpos = data[chrNum]
        if 'chr' not in chrNum:
            chrNum = 'chr'+chrNum
        info = dictOMIM[chrNum]
        startPositions = info[0]
        endPositions = info[1]
        geneDict = info[2]
        #print geneDict 
        genestoPos = {}
        for pos in listpos:
            #print startPositions
            index = modifiedBinarySearch(startPositions, pos)
            checkEnd = endPositions[index]
            checkEnd2 = endPositions[index-1]
            checkStart =  startPositions[index]
            checkStart2 =  startPositions[index-1]
            while checkEnd >= pos and pos >= checkStart:
                if pos not in genestoPos:
                    genestoPos[pos] = [geneDict[(checkStart, checkEnd)]]
                else:
                    genestoPos[pos].append(geneDict[(checkStart,checkEnd)])
                    index += 1
                    if index < len(endPositions):
                        checkEnd = endPositions[index]
                        checkStart =  startPositions[index]
                    else:
                        break
            #break
        #print genestoPos
        geneids = {}
        #adding hte counts of the genes to the geneids
        for key in genestoPos:
            for items in genestoPos[key]:
                for item in items:
                    if item[0] not in geneids:
                        geneids[item[0]] = 1
                    else:
                        geneids[item[0]]+=1
        # sorting the genes based on the occurance
        sorted_x = sorted(geneids.items(), key=operator.itemgetter(1), reverse=True)
        #print sorted_x[0:20]
        #print genestoPos
        
        chroms[chrNum] = sorted_x
    return chroms

#
# it gets all the positions from the VCF file
#
def getPositionsfromVCFFile(filename):
    chrNums = {}
    for line in open(filename):
        if '#' not in line:
            info = line.split('\t')
            chrNum = info[0]
            position = int(info[1])
            if chrNum not in chrNums:
                chrNums[chrNum] = [position]
            else:
                chrNums[chrNum].append(position)
    #print chrNums.keys()
    return chrNums


# the main fuctio whihc takes in the filename and the samplename and returns
# the phenotypes it asks for in the top number 
#
if __name__ == "__main__":
    filename = sys.argv[1]
    chrfile = sys.argv[2]
    dictOMIM = parseOMIMFile(filename)
    #print dictOMIM['chr22'][0]
    
    #testlist = [100,200,300,300,300,305,306,400,500]
    #print modifiedBinarySearch(testlist, 11122151)
    #""" 
    data = getPositionsfromVCFFile(chrfile)
    sorted_genes = binarySearchOMIM(dictOMIM, data)
    #"""
    #adding this now 
    sorted_genes = sorted_genes['chr22']
    #print sorted_genes
    #"""
    dictGeneIDtoPheno = createGeneIDtoPheno(filename)
    phenotypes = findPhenotypeFromOMIM(dictGeneIDtoPheno, sorted_genes, 20)
    #print len(phenotypes)
    for pheno in phenotypes:
        print pheno
    #"""
