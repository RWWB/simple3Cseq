#/bin/env python

__author__ = 'r.w.w.brouwer'

import os
import os.path
import sys
import getopt
from libutil import FastQ

def usage( mess='', error=None ):
    '''
    Print the usage information

    '''
    print '''
%s

Usage:

python %s -i [fastq file] -o [output tsv file]

Options:
-i/--input     the input FastQ file,default stdin
-o/--output    the output tab delimited file, default stdout
-h/--help       print this message

Description:

''' % (mess, sys.argv[0])

    if error != None:
        sys.exit( error )

def main( fin=sys.stdin, fout=sys.stdout ):
    '''

    :param fin:
    :param fout:
    :return:
    '''
    fastq = iter( FastQ( fin ) )
    for id, seq, qual in fastq:
        fout.write( "%s\t%s\t%s\n" % (id, seq, qual) )
    fout.close()

if __name__ == '__main__':

    # set our input field
    fin  = sys.stdin
    fout = sys.stdout

    # parse the options
    try:
        opts, args = getopt.getopt( sys.argv[1:], 'i:o:h',  ['input', 'output','help' ] )

        for o,a in opts:
            if o in ('-i','--input'):
                fin = open( os.path.abspath( a ), 'r')
            elif o in ('-o','--output'):
                fout = open( os.path.abspath( a ), 'w')
            elif o in ('-h', '--help'):
                usage( "Help was asked", error=0)

    except getopt.GetoptError, err:
        usage( str(err), error=2 )

    # convert the fastq format
    main( fin, fout )