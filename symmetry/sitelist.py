import numpy as np

class SiteList( list ):

    @classmethod
    def read_from_file( cls, filename ):
        """defines a SiteList instance by reading a set of coordinates x1, y1, z1, etc."""
        return cls( np.loadtxt( filename ).tolist() )

    def print_site( self, index ):
        """prints the list of sites as x1, y1, z1, etc."""
        print( ' '.join( [ str(e) for e in self[ index ] ] ) )
