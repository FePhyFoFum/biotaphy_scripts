import sys
import numpy as np
import math
import scipy.linalg as la
import dendropy as dp
import tree_reader
import aln_reader

"""
sq_change which produces the same results as ML without SE
"""
def calc_cont_anc_states(tree):
    df = 0
    nodenum = {}
    count = 0
    for i in tree.iternodes(order="postorder"):
        if i.istip:
            i.data['val'] = float(i.data['cont_values'][0])
            i.data['valse'] = float(i.data['cont_values'][0])
        else:
            nodenum[i] = count
            count += 1
            df += 1
            i.data['val'] = 0.
            i.data['valse'] = 0.
    df -= 1
    #compute the mlest of the root
    fullMcp = np.zeros((df+1,df+1))
    fullVcp = np.zeros(df+1)
    count = 0
    for i in tree.iternodes(order="postorder"):
        if i.istip == False:
            nni = nodenum[i]
            for j in i.children:
                tbl = 2./j.length
                fullMcp[nni][nni] += tbl;
                if j.istip:
                    fullVcp[nni] += (j.data['val'] * tbl)
                else:
                    nnj = nodenum[j]
                    fullMcp[nni][nnj] -= tbl;
                    fullMcp[nnj][nni] -= tbl;
                    fullMcp[nnj][nnj] += tbl;
            count += 1
    b = la.cho_factor(fullMcp)
    #these are the ML estimates for the ancestral states
    mle = la.cho_solve(b,fullVcp)
    sos = 0
    for i in tree.iternodes(order="postorder"):
        if i.istip == False:
            i.data['val'] = mle[nodenum[i]]
            #print i.data['val']
            i.label = str(mle[nodenum[i]])
            for j in i.children:
                temp = (i.data['val'] - j.data['val'])
                sos += temp*temp / j.length
    #print "Square Length: ",sos
    #calcSE
    """
    for i in tree.iternodes(order="postorder"):
        if i.istip == False:
            qpq = fullMcp[nodenum[i]][nodenum[i]]
            tm1 = np.delete(fullMcp,(nodenum[i]),axis=0)
            tm = np.delete(tm1,(nodenum[i]),axis=1)
            b = cho_factor(tm)
            sol = cho_solve(b,tm1[:,nodenum[i]])
            tempse = qpq - np.inner(tm1[:,nodenum[i]],sol)
            i.data['valse'] = math.sqrt(2*sos/(df*tempse))
    """
    return 0

def match_tips_and_cont_values(tree,seqs):
    lvs = tree.leaves()
    for i in lvs:
        test = False
        for j in seqs:
            if j.name == i.label:
                i.data['cont_values'] = j.cont_values
                test = True
                break
        if test == False:
            print "can't find "+i.label+" in cont_values"
            return False

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print "python "+sys.argv[0]+" newick.tre dataasphy"
        sys.exit(0)
    infile = open(sys.argv[1],"r")
    tree = tree_reader.read_tree_string(infile.readline())
    infile.close()    
    seqs = aln_reader.read_phylip_cont_file(sys.argv[2])
    match_tips_and_cont_values(tree,seqs)
    sqch = calc_cont_anc_states(tree)
    outfile = open("contanc.tre","w")
    outfile.write(tree.get_newick_repr(True)+";\n")
    outfile.close()
