import numpy as np

class SiteList( list ):
    """
    A :any:`SiteList` object is a list of coordinates (these can be Cartesian or fractional)
    that describe some set of sites. This can be used with a particular :any:`Configuration`
    to generate the coordinates of a particular arrangement of atoms.
    """

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
        Print coordinates of a single site: x1, y1, z1.

        Args:
            index (Int): index for the site coordinates to print

        Returns:
            None
        """
        print( ' '.join( [ str(e) for e in self[ index ] ] ) )
