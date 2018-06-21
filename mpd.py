import sys

import numpy as np

import tree_reader
import tree_utils

# mean phylogenetic distance between all pairs of species
def mpd(tree,tips):
    dists = []
    lt = list(tips)
    for i in range(len(lt)):
        for j in range(len(lt)):
            if i < j:
                dists.append(phylo_dist(lt[i],lt[j],tree))
    measure = np.mean(dists)
    return measure

def phylo_dist(tip1,tip2,tree):
    mrca = tree_utils.get_mrca([tip1,tip2],tree)
    d1 = 0
    d2 = 0
    cn = tip1
    while cn != mrca:
        d1 += cn.length
        cn = cn.parent
    cn = tip2
    while cn != mrca:
        d2 += cn.length
        cn = cn.parent
    return d1+d2

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

    print "mpd:",mpd(tree,tips)