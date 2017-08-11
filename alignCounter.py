#/bin/env python

import os
import os.path
import sys
import pysam
import getopt
from libutil import BedReader

class ReadCounter:

    def __init__( self, chr=None, start=None, end=None, padding=110 ):

        if chr == None: raise ValueError( 'chromosome must be defined' )
        if start == None: raise ValueError( 'start position must be defined' )
        if end == None: raise ValueError( 'end position must be defined' )

        # initialize the counters
        self.counter = 0
        self.forward = 0
        self.reverse = 0

        # initialize position specific counters
        self.p5  = [0,0,0]
        self.p3  = [0,0,0]
        self.oth = [0,0,0]

        # set the chromosome of the current region
        self.chr     = chr
        self.start   = start
        self.end     = end

        self.padding = padding

    def __call__( self, aln ): 
        ''' 
        The counter function to determine the reads on the regions

        Args:
            self    --   
            aln     --
        
        Returns:
            nothing
        ''' 

        # increase the general read counter
        self.counter += 1
        

        # increase the strand of the read counters
        if aln.is_reverse:
            self.reverse += 1
        else:
            self.forward += 1
 
        # 
        loc = self.alignmentLocation( aln )
        if loc:            
            f = 1
            if aln.is_reverse: f = 2 
            if loc == 5:
                self.p5[0] += 1
                self.p5[f] += 1
            elif loc == 3:
                self.p3[0] += 1
                self.p3[f] += 1
            else:
                self.oth[0] += 1
                self.oth[f] += 1

 
    def __str__( self ):
        '''
        Represent the region as a string 

        '''
        #
        sreg = "%s\t%d\t%d" % (self.chr, self.start, self.end)
        sgen = "%d\t%d\t%d" % (self.counter, self.forward, self.reverse)
        spr5 = "%d\t%d\t%d" % (self.p5[0], self.p5[1], self.p5[2])
        spr3 = "%d\t%d\t%d" % (self.p3[0], self.p3[1], self.p3[2])
        soth = "%d\t%d\t%d" % (self.oth[0], self.oth[1], self.oth[2])

        #
        return "%s\t%s\t%s\t%s\t%s" % ( sreg, sgen, spr5, spr3, soth )


    def alignmentLocation( self, aln ):
        '''
        Determine whether the read is on the 5' or 3' end of the region

        Args:
            self   --  
            aln    --   the alignment

        Returns:
            5, 3, -1 or None

        '''

        if self.start and self.end:

            # determine the length of the region
            l = self.end - self.start

            # impossible to make a good call due to a too small region size
            if l < (2 * self.padding):
                return None

            # get the relative position of the read to the region 
            relpos = aln.pos - self.start

            # check for the 5' 
            if relpos < self.padding: return 5
        
            # check for the 3'
            if relpos > l-self.padding: return 3

            # if not 5' or 3' return none 
            return -1 

        # end 
        return None    





def usage( mess='', error=None ):
    '''
    Print the usage information

    '''
    print '''
%s

Usage:

python %s -b [input bam file] -r [input bed file] -o [output file] 

Options:
-b/--bam       the input bam file
-r/--region    the input bed file
-o/--output    the output file
-h/--help      print this message

Description:

''' % (mess, sys.argv[0])

    #
    if error != None:
        sys.exit( error )



if __name__ == '__main__':

    #
    #
    #
    bamfname = None #sys.argv[1] 
    bedfname = None #sys.argv[2] 
    outfname = "-"

    # parse the options
    try:
        opts, args = getopt.getopt( sys.argv[1:], 'b:r:o:h',  ['bam', 'region','output', 'help' ] )

        for o,a in opts:
            if o in ('-b','--bam'):
                bamfname = os.path.abspath( a )
            elif o in ('-r','--region'):
                bedfname = os.path.abspath( a )
            elif o in ('-o','--output'):
                outfname = os.path.abspath( a )
            elif o in ('-h', '--help'):
                usage( "Help was asked", error=0)

    except getopt.GetoptError, err:
        usage( str(err), error=2 )


    # check whether everything is defined
    if not bamfname or not os.path.exists( bamfname ) :
        usage( "Could not find the input BAM file",  error=2)
    if not bedfname or not os.path.exists( bedfname ) :
        usage( "Could not find the input BED file",  error=2)

    #
    if outfname == "-":
        out = sys.stdout
    else:
        out = open( outfname, 'w' )

    # open the bam file and the bed
    samfile = pysam.Samfile( bamfname ,'rb' )
    bedfile = BedReader( bedfname ) 
    padding = 40

    # write the header to the file 
    out.write( '#bamfile=%s\n' % bamfname )
    out.write( '#bedfile=%s\n' % bedfname )
    out.write( '#padding=%d\n' % padding ) 
    out.write( '#chromosome\tstart\tend\ttotal.counts\ttotal.forward\ttotal.reverse\t5prime.counts\t5prime.forward\t5prime.reverse\t3prime.counts\t3prime.forward\t3prime.reverse\tmid.counts\tmid.forward\tmid.reverse\n' )
    for (chr, start, end, fld) in bedfile:
        
        if chr in samfile.references:

            # initialize a readcounter object 
            rcnt = ReadCounter( chr=chr, start=start, end=end, padding=padding )

            # count the matching reads
            it = samfile.fetch( chr, start, end )
            for aln in it:
                rcnt( aln ) 

            # print the output
            out.write( "%s\n" % rcnt ) #"%s\t%d\t%d\t%d\t%d\t%d" % ( chr, start, end, rcnt.counter, rcnt.forward, rcnt.reverse )
    #
    if outfname != "-" : out.close()



















