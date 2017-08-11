#!/bin/env python

__author__ = 'r.w.w.brouwer'
__version__ = '0.1'

import os
import os.path
import sys
import getopt
from libutil import FastQ

def usage( mess='', error=None):
    sys.stdout.write("""
%s

called as: %s

Version: %s

Usage:
%s -i <infile> -o <outfile> -s <int> -e <int>

Options:
    -i / --input        the input FastQ file, default stdin
    -o / --output       the output FastQ file, default stdout
    -s / --start        the first base to take, default 0
    -e / --end          the last base to take, default 20
    -h / --help         prints this message
    --table             print the data in tabular format

Please note that reads with fewer bases than start and end are
not written to the output.
""" % (mess, ' '.join(sys.argv), __version__, sys.argv[0] ) )

    if error is not None:
        sys.exit( int(error) )


def main( fin=sys.stdin, fout=sys.stdout, first=0, last=20, table=False):
    '''
    The main program loop

    :param fin: the input stream
    :param fout: the output stream
    :param first: the first base to take
    :param last: the last base to take
    :param table: should the output be written as tabular output
    :return: nothing
    '''
      # process the reads in the FastQ file
    fastq = iter( FastQ( fin ) )
    for id, seq, qual in fastq:

        # initialize the sequence
        #  and quality to get
        subseq  = None
        subqual = None

        # If we can extract a sub-sequence
        #  from the read.
        if len(seq) > first and len(seq) > last:
            subseq  = seq[first:last]
            subqual = qual[first:last]

        # If the sequence could not be obtained make
        #  a virtual sequence consisting of low-quality
        #  N bases
        if subseq is None:
            subseq  = 'N' * (last - first)
            subqual = '!' * (last - first)

        # Write the data as in FastQ format or
        #  as a table if the parameter switch was specified
        if not table:
            fout.write( "%s\n%s\n+\n%s\n" % (id, subseq, subqual) )
        else:
            fout.write( "%s\t%s\t%s\n" % (id, subseq, subqual) )

    # close the output file
    fout.close()


if __name__ == '__main__':

    # prepare variables
    fin   = sys.stdin
    fout  = sys.stdout
    first = 0
    last  = 20
    table = False

    # option parsing
    try:
        opts, args = getopt.getopt( sys.argv[1:], 'i:o:s:e:h',  ['input=','output=', 'start=', 'end=', 'table', 'help' ] )
        for o,a in opts:
            if o in ('-i','--input'):
                fin = open( os.path.abspath( a ), 'r')
            elif o in ('-o','--output'):
                fout = open( os.path.abspath( a ), 'w' )
            elif o in ('-s','--start'):
                first = int(a)
            elif o in ('-e','--end'):
                last = int(a)
            elif o in ('-h', '--help'):
                usage( "Help was asked", error=0 )
            elif o == '--table':
                table = True
    except getopt.GetoptError as err:
        usage( str(err), error=2 )

    # check the coordinates
    if last < first:
        usage( "The last base to take is smaller than the first base (%d-%d)" % (last, first), error=101 )

    # run the program
    main( fin=fin, fout=fout, first=first, last=last, table=table )