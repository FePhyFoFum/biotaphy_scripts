import sys
import dendropy as dp


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print "python "+sys.argv[0]+" newick.tre comma_seperated_taxa"
        sys.exit(0)

    tree = dp.Tree.get(path=sys.argv[1],schema="newick")
    pdm = tree.phylogenetic_distance_matrix()
    comtax = sys.argv[2].split(",")
    n = len(comtax)
    numcom = ((n*n)-n)/2.
    comfil = lambda taxon : taxon.label in comtax
    print "spd:",pdm.mean_pairwise_distance(filter_fn=comfil)*numcom

