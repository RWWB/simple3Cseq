

import sys

'''
A utility library for Python 

''' 


class BedWriter:
    '''
     A very simple bed writer

    '''
    def __init__( self, fname ) :

        if fname == '-':
            self._f = sys.stdout
        else:
            self._f = open( fname, 'w' )

    def add( self, chr, start, end, fields=[] ):
        self._f.write(  '%s\t%d\t%d' % (chr, start, end)  )
        if len( fields ) > 0: self._f.write( "\t%s" % '\t'.join( fields ) )
        self._f.write( "\n" )


    def close(self):
        self._f.close()



# this class has been copied over from the NARWHAL pipeline
class FastA:
    ''' A class to iterate over the fasta file  '''
    def __init__( self, fname ):
        self._fn = fname
        self._f  = open( fname, 'rU' )
        self._l = None

    def __iter__( self ):
        return self

    def next( self ):
        '''

        '''
        seqid = None
        seq   = None

        # determine the next identifier
        if self._l == None:
            self._l = self._f.next().rstrip()

        if self._l != '' and self._l[0] == '>':
            # assign the identifier
            seqid = self._l[1:]

            # add the sequence
            seq   = ''
            self._l = self._f.next().rstrip()
            while self._l != '' and self._l[0] != '>':
                seq += self._l
                self._l = self._f.next().rstrip()

        # return the sequence id and the sequence
        return( seqid, seq )

    def rewind( self ):
        '''

        '''
        self._f.close()
        self._f = open( self._fn, 'rU' )

    def close( self ):
        self._f.close()

class BedReader:
    '''
    A very simple bed reader

    '''

    def __init__( self, fname, sep="\t" ):
        '''
        
        '''
        self._f   = open( fname, 'rU' )
        self._sep = sep

    def __iter__( self ):
        '''
        Return the iterator for this object

        Args:
            self

        Returns:
            an iterator (self)
        '''
        return self

    def next( self ):
        '''
        Reads the next entry in the BED file
       
        It uses the StopIteration exception from the file.next function
 
        Args:
          self  --  

        Returns:
          a tuple with sequence name, start, end, additional fields
        '''
        l = self._f.next()
        l = l.rstrip()
        a = l.split( self._sep )

        # get the core region representation
        snm = a[0]
        sta = int( a[1] )
        end = int( a[2] )

        # get the additional fields        
        fld = []
        if len(a) > 3: fld = a[3:]
        
        # get the region tuple
        rval = ( snm, sta, end, fld )
        
        # return the bed representation
        return( rval )

    def close( self ):
        self._f.close()

class FastQ(object):
    '''


    '''
    def __init__(self, fin=sys.stdin):
        '''

        :param fin:
        :return:
        '''
        self._fin = fin

    def __del__(self):

        # self closing object
        if not self._fin.closed:
            self._fin.close()

    def __iter__(self):

        lineno   = 0
        recordno = 0
        while True:

            # check the id
            lineno += 1
            line    = self._fin.next().rstrip()
            if len(line) == 0 or line[0] != '@':
                raise ValueError( 'invalid header at line %d' % lineno)
            id = line[1:]

            # check the sequence
            lineno += 1
            line    = self._fin.next().rstrip()
            seq     = line

            # check the secondary header
            lineno += 1
            line    = self._fin.next().rstrip()
            if len(line) == 0 or line[0] != "+":
                raise ValueError( 'invalid secondary header at line %d' % lineno)

            # get the quality string
            lineno += 1
            line    = self._fin.next().rstrip()
            qual    = line

            recordno += 1
            if len(seq) != len(qual):
                raise ValueError( 'the sequence and quality strings are not of equal length for record %d' % recordno)
            # return the data
            yield( id, seq, qual)

if __name__ == '__main__':
    pass

