#/bin/env python

__author__ = 'r.w.w.brouwer'

import os
import os.path
import sys
import getopt
import re
from libutil import FastQ
from findSequence import IUPAC_nucleotide_decode

def usage( mess='', error=None ):
    '''
    Print the usage information

    '''
    print '''
%s

Description:

Extracts reads with the recognition sequence from the input file and
 writes them to  the output file.

Usage:

python %s -i [input fastq file] -o [output fastq file]

Options:
-i/--input     the input FastQ file,default stdin
-o/--output    the output FastQ file, default stdout
-s/--sequence  the recognition sequence in IUPAC bases, default None
-h/--help       print this message

Description:

''' % (mess, sys.argv[0])

    if error != None:
        sys.exit( error )

def main( fin=sys.stdin, fout=sys.stdout, sequence=None ):
    '''

    :param fin:
    :param fout:
    :return:
    '''

    #
    assert sequence is not None

    regexp  = IUPAC_nucleotide_decode( sequence )
    matcher = re.compile( regexp )

    fastq = iter( FastQ( fin ) )
    for id, seq, qual in fastq:

        # only print reads that match
        if matcher.search(seq):
            fout.write( "%s\t%s\t%s\n" % (id, seq, qual) )
    fout.close()

if __name__ == '__main__':

    # set our input field
    fin  = sys.stdin
    fout = sys.stdout
    seq  = None

    # parse the options
    try:
        opts, args = getopt.getopt( sys.argv[1:], 'i:o:s:h',  ['input', 'output','sequence', 'help' ] )

        for o,a in opts:
            if o in ('-i','--input'):
                fin = open( os.path.abspath( a ), 'r')
            elif o in ('-o','--output'):
                fout = open( os.path.abspath( a ), 'w')
            elif o in ('-s','--sequence'):
                seq = a
            elif o in ('-h', '--help'):
                usage( "Help was asked", error=0)

    except getopt.GetoptError, err:
        usage( str(err), error=2 )

    if seq is None:
        usage( 'Recognition sequenc not provided', error=101 )

    # convert the fastq format
    main( fin, fout, sequence=seq )