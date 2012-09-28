import os
import logging

from fsl.melodic import *
from optparse import OptionParser

class Classification():
    class_name = ""
    ic_number = 0
    filter = bool(False)
    
    def __init__(self, class_name, ic_number, filter):
        self.class_name = class_name
        self.ic_number = ic_number
        self.filter = filter
        
    def __str__(self):
        return "{0}, {1}, {2}".format(self.ic_number, self.class_name, self.filter)
        
class ClassificationFile():
    
    def __init__(self, filename, **kwargs):
        self.path = os.path.join('.', filename)
        
    def write(self, class_list, dirpath):
        try:
            f = open(self.path, 'w')
    
            f.write(dirpath + '\n')
#            for c in [ c for c in class_list if c.class_name != 'Unknown' ]:
            for c in class_list:
                f.write( '{0}\n'.format(c) )    
            f.write( str([ int(ic.ic_number) for ic in class_list if ic.filter == True]) + '\n')
    
            f.close()

        except Exception, e:
            print e
            logging.error('Error writing file: {0}'.format(e))


    def read(self):    
        try:
            f = open(self.path, 'r')
    
            cl = []

            lines = f.readlines()
            icadirpath=lines[0].strip()
            mfname = "{0}/melodic_mix".format(icadirpath)
            mix = genfromtxt(mfname)
            npts, nics = mix.shape

            for ic in arange(0, nics):
                cl.append( Classification(ic_number=ic+1, class_name='Unknown', filter=False) )

            doFix = True
            for i,fv,c in [ (int(i),(fv.strip() == 'True'),c.strip()) 
                           for (i,c,fv) in [l.split(',') for l in lines[1:-1]] ] :
    #        print i, fv, c
                cl[i-1] = Classification(ic_number=i, class_name=c, filter=fv)
                doFix = False

            last_line=lines[-1]
            
            if doFix:
                import re
                if re.match('\[[0-9, ]*\]\n$', last_line):
                    
                    for ic in eval(last_line):
                        cl[ic-1].class_name='Unclassified Noise'
    
                else:
                    raise Exception("Bad input in last line of classifications file!")
                        
            return cl, icadirpath
    
        except Exception, e:
            print e
            logging.error('Error reading file: {0}'.format(e))

def main():
    parser = OptionParser()

    parser.add_option("-f", "--infile", dest="infilename",
                      help="read from FILE", metavar="FILE")
    parser.add_option("-o", "--outfile", dest="outfilename",
                      help="write to FILE", metavar="FILE")
    
    (options, args) = parser.parse_args()
    
    inf = ClassificationFile(options.infilename)
    class_list, dirpath = inf.read()
    outf = ClassificationFile(options.outfilename)
    outf.write(class_list, dirpath)
    
if __name__ == "__main__":
    main()

