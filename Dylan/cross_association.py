""" Implmentation of Cross-Association Clustering

Follows the paper 
http://www.cs.cmu.edu/~christos/PUBLICATIONS/kdd04-cross-assoc.pdf


Takes a matrix, M, with the following functions defined:

Given a 5 x 3 matrix with 3 x 2 clusters in the row. e.g.,

    M.rows_by_clusters() => [[5,2],[1,4],[3]]
    M.columns_by_clusters() => [[1],[2],[3]]
    M.clusters() => ([[5,2],[1,4],[3]], [[1],[2],[3]])
    M.members_in_row_cluster(2) => [1,4]
    M.members_in_column_cluster(2) => [2]
    M.move_to_new_row_cluster([1,2,4])
    M.move_to_new_colum_cluster([1,2,4])


"""
import random, math
import scipy

def argmax_random_tie(seq, fn):
    """Return an element with highest fn(seq[i]) score; break ties at random.
    Thus, for all s,f: argmin_random_tie(s, f) in argmin_list(s, f)
    
    from http://aima-python.googlecode.com/svn/trunk/utils.py """
    best_score = fn(seq[0]); n = 0
    for x in seq:
        x_score = fn(x)
        if x_score > best_score:
            best, best_score = x, x_score; n = 1
        elif x_score == best_score:
            n += 1
            if random.randrange(n) == 0:
                    best = x
    return best

def entropy(values):
    "Number of bits to represent the probability distribution in values."
    # If the values do not sum to 1, normalize them to make them a Prob. Dist.
    values = removeall(0, values)
    s = float(sum(values))
    if s != 1.0: values = [v/s for v in values]
    return sum([- v * log2(v) for v in values])

def reassignment_decreases_per_row_entropy(before_cluster, after_cluster):
    pass


def regroup(step=0, k, l):
    """ Find the cross associates (cluster assignments for rows and columns)
    using an alternating minimization algorithm. """
   
    # cross_associations
    row_clusters, col_clusters = M.clusters()



def cross_association_search(M, step=0, k=1, l=0):
    print "Cross Association Iteration %i" % (step)
    
    previous_cost = compute_cost(M)
    # at each iteration, find the row group with the maximum entropy 
    max_entropy_row_group = arg_max_randon_tie(M.rows_by_clusters(), entropy)
    # split that group, add some members to new group
    threshold_function = reassignment_decreases_per_row_entropy(before_cluster,after_cluster)
    to_be_reassigned = filter(threshold_function, M.members_in_row_cluster(max_entropy_row_group))
    if len(to_be_reassigned) > 0:
        M.move_to_new_row_cluster(to_be_reassigned)



