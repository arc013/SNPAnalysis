import sys
from numbers import Number
from collections import Set, Mapping, deque
import operator
import numpy as np
import time

try: # Python 2
    zero_depth_bases = (basestring, Number, xrange, bytearray)
    iteritems = 'iteritems'
except NameError: # Python 3
    zero_depth_bases = (str, bytes, Number, range, bytearray)
    iteritems = 'items'

individuals = {}

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
	
header = ""
pos_dict = {}

def parseChromosome(filename, folder):
    start = time.clock()
    
    samplenames = []
    comments = False

    for line in open(filename):
        info = line.split("\t")
        if len(info) > 10 and comments == False:
            header = header + line + '\n'
			samplenames = info[9:]
            comments = True
			
        if "#" not in line:
            chrPos = info[1]
			pos_dict[chrPos] = info[0:9]
            for i, sample in enumerate(info[9:]):
                data = info[i+9]
                if data != "0|0":
                    samplename = samplenames[i]
                    if samplename not in individuals:
                        individuals[samplename] = [(chrPos,data)]
                    else:
                        individuals[samplename].append((chrPos,data))
    end = time.clock()
    print "time:",end - start 
    for key in individuals:
        print "sample:", key
        with open(folder + key+".vcf", "w") as text_file:
			text_file.write(header)
            values = individuals[key]
            for position,data in values:
                text_file.write("\t".join(pos_dict[position])+'\t'+data+"\n")
    

    return
                                                



if __name__ == "__main__":
    filename = sys.argv[1]
    folder = sys.argv[2]
    parseChromosome(filename, folder)
    print getsize(individuals)
