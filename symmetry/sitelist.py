import numpy as np

class SiteList( list ):

    @classmethod
    def read_from_file( cls, filename ):
        return cls( np.loadtxt( filename ).tolist() )

    def print_site( self, index ):
        print( ' '.join( [ str(e) for e in self[ index ] ] ) )
