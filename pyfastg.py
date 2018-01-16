import sys
import os
import re

f = "./assembly_graph.fastg"

class FastgNode(object):
    def __init__(self, name=None, length=None, cov=None, reverse_complimented = None, neighbor_list=None):
        self.name = name
        self.length = length
        self.cov = cov
        self.neighbor_list = neighbor_list
        self.reverse_complimented = reverse_complimented # only used during creation to combine neighbors
    def __str__(self):
        return("Node: {0}\nNeighbors: {1}\nLength: {2}\nCoverage: {3}\nReverse_Complimented?: {4}".format(
            self.name,
            "None" if self.neighbor_list is None else ",".join([str(x.name) for x in self.neighbor_list]),
            "None" if self.length is None else self.length,
            "None" if self.cov is None else self.cov,
            str(self.reverse_complimented)
        )
    )
class FastgNodeNeighbor(object):
    def __init__(self, name=None,reverse_complimented=None):
        self.name = name
        self.reverse_complimented = reverse_complimented
    def __str__(self):
        return("NeighborNode: {0}\nReverse_Complimented?: {1}".format(
            self.name,
            "Yes" if self.reverse_complimented else "No"
        )
    )


def extract_node_len_cov_rc(node_name):
    rc = False
    if node_name.endswith("'"):
        rc = True
        node_name = node_name[0:-1]
    p = re.compile(r'EDGE_(?P<node>\d*?)_length_(?P<length>\d*?)_cov_(?P<cov>[\d|\.]*)')
    m = p.search(node_name)
    return (m.group("node"), int(m.group("length")), float(m.group("cov")), rc)


def make_Node(name):
    node_name, length, cov, rc = extract_node_len_cov_rc(name)
    new_node = FastgNode(
        name=node_name,
        length=length,
        cov=cov,
        reverse_complimented = rc
    )
    return(new_node)

def make_Neighbors(neighbor):
    node_name, length, cov, rc = extract_node_len_cov_rc(neighbor)
    new_neigh = FastgNodeNeighbor(
        name=node_name,
        reverse_complimented=rc
    )
    return(new_neigh)


def parse_fastg(f):
    node_neighs = []
    with open(f, "r") as inf:
        for line in inf:
            if line.startswith(">"):
                colons = sum([1 for x in line if x == ":" ])
                if colons > 1:
                    sys.stderr.write("multiple ':'s found in line, and can only " +
                                     "be used to separate nodes from neighbor " +
                                     "list\n")
                elif colons == 0:
                    # orphaned node
                    pass
                    # sys.stderr.write("Header does not contain a colon!\n")
                else:
                    node_neighs.append(line.strip().split(":"))

    node_list = []
    for node, neighs in node_neighs:
        new_node = make_Node(node)
        new_node.neighbor_list = [make_Neighbors(x) for x in neighs.split(",")]
        node_list.append(new_node)
    # # now, filter outthe duplicate nodes by combining the forward and rc hits
    # new_node_list = []
    # for node in node_list:
    #     # if name already in new list, add the neighbors from the reverse
    #     if node.name in [x.name for x in new_node_list]:

    return node_list

def pathfind(node_list, top_parent, parent, prev_path, prev_length, path_list, thresh=1000, found_exit=False):
# def pathfind(node_list, parent, prev_path, prev_length, path_list, reject_direction_node, thresh=1000):
    """Returns deeper path, or None if path is already in path list
    """
    print("Prev: " + prev_path)
    print("Prev_len: " + str(prev_length))
    # look for both forward and rc matches
    parent_opposite_strand = [x for x in node_list if x.name == parent.name and x.reverse_complimented != parent.reverse_complimented][0]
    poss_f = [x for x in parent.neighbor_list if  x.name not in prev_path.split(":")]
    poss_rc = [x for x in parent_opposite_strand.neighbor_list if  x.name not in prev_path.split(":")]
    # this ensures we only get single direction hits from the topmost parent
    if parent == top_parent:
        print("THALKJSNHVKAJNVAWE")
        possible_neighbors = poss_f
    else:
        possible_neighbors = poss_f + poss_rc
    print("Poss:" + " ".join([x.name for x in possible_neighbors]))

    for node in possible_neighbors:
        # here we assume only one hit
        this_node = [x for x in node_list if x.name == node.name and x.reverse_complimented == node.reverse_complimented][0]
        if this_node.length > 250:
            found_exit = True
        print("This: " + this_node.name + " RC: " + str(this_node.reverse_complimented))
        this_length = prev_length + this_node.length
        this_path = prev_path + ":" + this_node.name
        # are we outside the zone of flanking similarity? usually 1kb?
        # and have we hit a decent stretch of sequence? one node longer than 3x a tRNA, or about 270
        if this_length >= thresh and found_exit:
            if this_path in path_list:
                pass
            else:
                path_list.append(this_path)
                # return True
        else:
            # further in, further in!
            pathfind(
                node_list=node_list,
                top_parent=top_parent,
                parent=this_node,
                prev_path=this_path,
                prev_length=this_length,
                path_list=path_list,
                found_exit=found_exit,
                thresh=thresh)
        print("And on to the next one!")
    return path_list


if __name__ == "__main__" or True:
    nodes = parse_fastg(f=f)
    start = [x for x in nodes if x.name == "2" and  x.reverse_complimented ][0]
    print(start)
    print(start.neighbor_list[1])
    l = pathfind(
        node_list=nodes,
        top_parent=start,
        parent=start,
        prev_path="init",
        prev_length=0,
        thresh=1000,
        path_list=[],
        found_exit=False)
    print(l)
    # keep_pathfinding = True
    # while keep_pathfinding:
    #     for node in nodes:
    #         print(node)
