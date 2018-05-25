"""
@summary: Hack functions to get data out of files that Lifemapper uses without
             adding dependencies.  
"""

def munge_lm_tree(tree_filename):
   """
   Return the newick string and a dictionary of squid : species pairs
   """
   squid_dict = {}
   with open(tree_filename) as inF:
      line = inF.next()
      while line.strip().lower() != 'taxlabels':
         line = inF.next()
      
      # Create squid dictionary
      taxon = inF.next()
      while taxon.strip() != ';':
         annotations = inF.next()
         # Squid is after "squid=" and before next comma
         squid = annotations.split('squid=')[1].split(',')[0]
         
         # Note: May need to get rid of quotations too
         squid_dict[squid] = taxon.strip()
         taxon = inF.next()
   
      # Get to tree
      line = inF.next()
      while line.strip().lower() != 'begin trees;':
         line = inF.next()

      # Get tree
      line = inF.next()
      # Start with first open paren and assume the rest is tree until end
      newick = line[line.find('('):]
      line = inF.next()
      while line.strip().lower() != 'end;':
         newick += line.strip()
         line = inF.next()
   
   return newick.strip(), squid_dict

# .............................................................................
if __name__ == '__main__':
   fn = 'examples/Geometridae.nex'
   newick, squid_dict = munge_lm_tree(fn)
   print newick
   print len(squid_dict.keys())