import sys
import os
import re


class FastgNode(object):
    def __init__(self, name=None, length=None, cov=None,
                 reverse_complimented=None, neighbor_list=None):
        self.name = name
        self.length = length
        self.cov = cov
        self.neighbor_list = neighbor_list
        self.reverse_complimented = reverse_complimented
    def __str__(self):
        return str("Node: {0}\nNeighbors: {1}\nLength: {2}\nCoverage: " +
                   "{3}\nReverse_Complimented?: {4}").format(
                       self.name,
                       "None" if self.neighbor_list is None else
                           ",".join([str(x.name) for x in self.neighbor_list]),
                       "None" if self.length is None else self.length,
                       "None" if self.cov is None else self.cov,
                       str(self.reverse_complimented)
                   )

class FastgNodeNeighbor(object):
    def __init__(self, name=None,reverse_complimented=None):
        self.name = name
        self.reverse_complimented = reverse_complimented
    def __str__(self):
        return "NeighborNode: {0}\nReverse_Complimented?: {1}".format(
            self.name,
            "Yes" if self.reverse_complimented else "No"
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
    return new_node

def make_Neighbors(neighbor):
    node_name, length, cov, rc = extract_node_len_cov_rc(neighbor)
    new_neigh = FastgNodeNeighbor(
        name=node_name,
        reverse_complimented=rc
    )
    return new_neigh


def parse_fastg(f):
    node_neighs = []
    with open(f, "r") as inf:
        for line in inf:
            if line.startswith(">"):
                colons = sum([1 for x in line if x == ":" ])
                if colons > 1:
                    raise ValueError(
                        "multiple ':'s found in line, and can only " +
                        "be used to separate nodes from neighbor " +
                        "list\n")
                elif colons == 0:
                    # orphaned node or terminal node
                    node_neighs.append([line.strip().split(";")[0], None])
                else:
                    node_neighs.append(line.strip().split(":"))

    node_list = []
    for node, neighs in node_neighs:
        new_node = make_Node(node)
        if neighs is None:
            new_node.neighbor_list = None
        else:
            new_node.neighbor_list = [make_Neighbors(x) for x in neighs.split(",")]
        node_list.append(new_node)

    return node_list


if __name__ == "__main__":
    f = "./assembly_graph.fastg"
    nodes = parse_fastg(f=f)
    for node in nodes:
        print(node)
