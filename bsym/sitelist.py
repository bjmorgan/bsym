import numpy as np

class SiteList( list ):

    @classmethod
    def read_from_file( cls, filename ):
        """
        Generate a SiteList instance by reading a set of coordinates x1, y1, z1, etc.

        Args:
            filename (str): name of the file to be read.

        Returns:
            a new SiteList instance
        """
        return cls( np.loadtxt( filename ).tolist() )

    def print_site( self, index ):
        """
        Print a subset of sites as x1, y1, z1, etc.

        Args:
            index (list): list of indices for which site coordinates to print

        Returns:
            None
        """
        print( ' '.join( [ str(e) for e in self[ index ] ] ) )
