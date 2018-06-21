import sys

import numpy as np

from mpd import phylo_dist
import tree_reader
import tree_utils

# sum of phylogenetic distances
# mpd * number of species pairs
def spd(tree,tips):
    dists = []
    ncompared = 0
    su = 0
    lt = list(tips)
    for i in range(len(lt)):
        for j in range(len(lt)):
            if i < j:
                m = phylo_dist(lt[i],lt[j],tree)
                dists.append(m)
                ncompared += 1
    measure = ncompared * np.mean(dists)
    return measure

"""
this main is not checking lots of error things because it is just a
demo of how to use the function
"""
if __name__ == "__main__":
    if len(sys.argv) != 3:
        print "python "+sys.argv[0]+" newick.tre comma_seperated_taxa"
        sys.exit(0)
    
    tree = tree_reader.read_tree_file_iter(sys.argv[1]).next()
    ndsdict = {} # key is name, value is node object
    for i in tree.leaves():
        ndsdict[i.label] = i

    taxa = sys.argv[2].split(",")
    tips = set()
    for i in taxa:
        try:
            tips.add(ndsdict[i])
        except:
            print i,"not in tree"
            print "exiting..."
            sys.exit(1)

    print "spd:",spd(tree,tips)