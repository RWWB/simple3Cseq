#/bin/env python


import os
import os.path
import sys
import getopt
from libutil import BedWriter, BedReader



def readSeqSizes( fname ):
    '''
    Read the genome sizes from a file

    ''' 
    f  = open( fname, 'rU' )

    # get the size from the 
    rv = {}
    for l in f:
        if l[0] != '#':
            l = l.rstrip()
            a = l.split( '\t' )
            rv[a[0]] = int( a[1] )
    
    f.close()

    # return the sizes dict
    return rv


def imputeRegions(  bedin, bedout, sizes={} ):
    '''  
    Impute the in-between regions

    '''

    pr = None
    for reg in bedin:
         # for the first region
         if not pr:
             bedout.add( reg[0], 0, reg[1] )
         else: 

             # we have traversed to another chromosome
             if pr[0] != reg[0]:
                 # check whether the chromosome lengths are present
                 if sizes.has_key( pr[0] ): bedout.add( pr[0], pr[2], sizes[pr[0]] )

                 # add the first region again if the chromosome does not start with a restriction site
                 if reg[1] != 0 : 
                     bedout.add( reg[0], 0, reg[1] )
             else:
                 # add the in-between region to the output 
                 bedout.add( reg[0], pr[2], reg[1] )
         # assign the current read to the previous read
         pr = reg

    # return nothing
    return None



def usage( mess='', error=None ):
    '''
    Print the usage information

    '''
    print '''
%s

Usage:

python %s -i [input bed file] -o [output bed file] -s [chromosome sizes]

Options:
-i/--input     the input bedfile
-o/--out       the output bedfile 
-s/--sizes     the chromosome sizes
-h/--help      print this message

Description:

''' % (mess, sys.argv[0])

    #
    if error != None:
        sys.exit( error )




if __name__ == '__main__':
    pabedin  = None
    pabedout = None 
    pasizes  = None 

    # parse the options
    try:
        opts, args = getopt.getopt( sys.argv[1:], 'i:o:s:h',  ['input', 'output','sizes', 'help' ] )

        for o,a in opts:
            if o in ('-i','--input'):
                pabedin = os.path.abspath( a )
            elif o in ('-o','--output'):
                pabedout = os.path.abspath( a )
            elif o in ('-s','--sizes'):
                pasizes = os.path.abspath( a )
            elif o in ('-h', '--help'):
                usage( "Help was asked", error=0)

    except getopt.GetoptError, err:
        usage( str(err), error=2 )


    # check whether everything is defined
    if not pabedin or not os.path.exists( pabedin ) :
        usage( "Could not find the input BED file",  error=2)
    if not pabedout:
        pabedout = "-"
    if not pasizes or not os.path.exists( pasizes ) :
        usage( "Could not find the chromosome sizes file", error=2 )

    # open the  
    sizes  = readSeqSizes( pasizes )
    bedin  = BedReader( pabedin ) 
    bedout = BedWriter( pabedout )

    #
    imputeRegions( bedin, bedout, sizes )

    #
    bedin.close()
    bedout.close()
    



