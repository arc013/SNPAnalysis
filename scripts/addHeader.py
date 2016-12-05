import sys
import glob 
import os
import subprocess


def getVCFfiles(directory):
    filenames = []
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith(".vcf"):
                filenames.append(os.path.join(root,file))   
    return filenames

def addHeader(orig_file, directory):
    header_process = subprocess.Popen(["grep",'#',orig_file], stdout=subprocess.PIPE)
    output, error = header_process.communicate()
    filenames = getVCFfiles(directory)
    
    for file in filenames:
        
        ## print file
        f = open(file, 'r')
        temp = f.read()
        sample_name = file.split('/')[-1].rstrip('.vcf')
        # print sample_name
        f.close()
        ## print output
        keep_top = output.split('\n')[:-2]
        ## print keep_top
        new_last_line = '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t' + sample_name + '\n'
        ## print new_last_line
        keep_top.append(new_last_line)
        f = open(file, 'w')
        f.write('\n'.join(keep_top))
        f.write(temp)
        f.close()

        #break
if __name__ == "__main__":
    orig_file = sys.argv[1]
    directory = sys.argv[2]
    addHeader(orig_file, directory)














