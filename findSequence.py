#/bin/env python 


import os
import os.path
import sys
import getopt
import re
from libutil import BedWriter, FastA


def findSequence(  fasta, bedwriter, sequence='ACTG' ):
    '''
    find all occurences of sequence in a FastA file and 
    write their coordinates to a BED file

    Args:
       fasta     --   a fasta file with the sequences to search in 
       bedwriter --   a file 
       sequence  --   a regular expression with the sequence to search for

    Returns:
       the number of sequences found

    '''
    occur  = 0
    regexp = re.compile( sequence, flags=re.IGNORECASE )
    for chr,seq in fasta:

            # for each found site
            for grp in regexp.finditer( seq ):
                bedwriter.add( chr, grp.start(), grp.end() )
                occur += 1
    # return the number of occurences of sequence
    return occur

def IUPAC_nucleotide_decode( seq ):
    '''
    decodes a sequence which is IUPAC and convert this to a regular expression friendly sequence

    '''
    # convert IUPAC bases
    seq = seq.replace( 'R', '[A,G]' )
    seq = seq.replace( 'Y', '[C,T]' )
    seq = seq.replace( 'S', '[G,C]' )
    seq = seq.replace( 'W', '[A,T]' )
    seq = seq.replace( 'K', '[G,T]' )
    seq = seq.replace( 'M', '[A,C]' )
    seq = seq.replace( 'B', '[C,G,T]' )
    seq = seq.replace( 'D', '[A,G,T]' )
    seq = seq.replace( 'H', '[A,C,T]' )
    seq = seq.replace( 'V', '[A,C,G]' )
    seq = seq.replace( 'N', '[A,C,T,G]' )

    # returns the sequence
    return seq


def usage( mess='', error=None ):
    '''
    Print the usage information

    '''
    print '''
%s

Usage:

python %s -f [fasta file] -b [output bed file] -s [sequence] 

Options:
-f/--fasta     a fasta file with the sequences to search in.
-b/--bed       a bed file to which to write the found occurences 
-s/--sequence  the sequence to search for
-h/--help       print this message

Description:

''' % (mess, sys.argv[0])

    #
    if error != None:
        sys.exit( error )



if __name__ == '__main__':
    '''

    '''
 
    # set our input field
    pafasta  = None
    pabed    = None
    sequence = None

    # parse the options
    try:
        opts, args = getopt.getopt( sys.argv[1:], 'f:b:s:h',  ['fasta', 'bed','sequence', 'help' ] )

        for o,a in opts:
            if o in ('-f','--fasta'):
                pafasta = os.path.abspath( a )
            elif o in ('-b','--bed'):
                pabed = os.path.abspath( a )
            elif o in ('-s','--sequence'):
                sequence = a
            elif o in ('-h', '--help'):
                usage( "Help was asked", error=0)

    except getopt.GetoptError, err:
        usage( str(err), error=2 )


    # check whether everything is defined
    if not pafasta or not os.path.exists(pafasta) :
        usage( "Could not find the input FastA file",  error=2)
    if not pabed: 
        pabed = "-"
    if not sequence:
        usage( "A sequence to find has to be specified", error=2 )


    # convert the IUPAC sequence to a regular expression
    sequence = IUPAC_nucleotide_decode( sequence )


    # finally start to do some work
    fasta = FastA( pafasta )
    bedwr = BedWriter( pabed )
    occur = findSequence( fasta, bedwr, sequence)

    bedwr.close()
    fasta.close()

    sys.stderr.write( "[%d occurences of %s found in %s]\n" % (occur, sequence, pafasta) )        


