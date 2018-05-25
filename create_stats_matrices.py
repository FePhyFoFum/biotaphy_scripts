"""
@summary: This script generates matrix objects containing phylo stats for a PAM 
             and tree
"""
import argparse
import numpy as np

from lm_munge import munge_lm_tree
from matrix import Matrix
from pd import pd
from tree_reader import read_tree_string

# .............................................................................
def get_site_statistics(tree, pam, squid_dict):
    """
    @summary: Calculate site statistics
    @note: To add additional stats:
              1. Add an additional numpy array
              2. Add stat in "Site statistics"
              3. Add stat array and header in return
    """
    num_rows, num_cols = pam.data.shape
    squid_headers = pam.getColumnHeaders()
    
    # Stat arrays
    # ..............
    pd_stat = np.zeros((num_rows, 1), dtype=np.float)
    
    for i in xrange(num_rows):
        # Get the species names that are present at this site
        present_squids = [
           squid_headers[j] for j in xrange(num_cols) if int(pam.data[i,j])]
        # Get the species names for the present squids
        sp_tips = [squid_dict[squid] for squid in present_squids]
        
        # Site statistics
        # ..............
        pd_stat[i,0] = pd(tree, sp_tips)
        
        # TODO: Add any additional statistics
   
    # Return statistics
    # Site statistics np.arrays
    site_stats = [pd_stat]
    # Site statistic headers
    stat_headers = ['PD']
    return Matrix(np.concatenate(site_stats, axis=1),
                  headers={'0' : pam.getRowHeaders(),
                           '1' : stat_headers})

# .............................................................................
def get_species_statistics(tree, pam, squid_dict):
    """
    @todo: Fill in if we have any species based statistics
    """
    return None

# .............................................................................
if __name__ == '__main__':
   
    parser = argparse.ArgumentParser(
      description='Calculate phylogenetic statistics for a PAM and a tree')
    parser.add_argument('pam_matrix_filename', type=str, 
                       help='File path to PAM matrix (in LM format)')
    parser.add_argument('tree_filename', type=str, 
                       help='File path to tree (NEXUS)')
    parser.add_argument('-sites', '--sites_filename', type=str, 
                  help='If provided, write site stats to this location (CSV)')
    parser.add_argument('-sp', '--species_filename', type=str, 
               help='If provided, write species stats to this location (CSV)')
   
    args = parser.parse_args()
    
    # Read in data
    pam = Matrix.load(args.pam_matrix_filename)
    # We use a hack function to get the data out of the tree, in the LM code
    #    this is done with dendropy.  I just didn't want to add any dependencies
    tree_newick, squid_dict = munge_lm_tree(args.tree_filename)
    
    print tree_newick
    print len(squid_dict.keys())
    
    # Read the newick tree
    tree = read_tree_string(tree_newick)
    
    # Generate stats matrices
    site_stats = get_site_statistics(tree, pam, squid_dict)
    
    species_stats = get_species_statistics(tree, pam, squid_dict)
    
    
    # Write matrices if desired
    if args.sites_filename is not None:
       with open(args.sites_filename, 'w') as outF:
          site_stats.writeCSV(outF)
    
    if args.species_filename is not None:
       raise Exception, 'No species stats are currently implemented'
       with open(args.species_filename, 'w') as outF:
          species_stats.writeCSV(outF)
          