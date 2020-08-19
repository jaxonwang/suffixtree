import sys
from graphviz import Digraph

if __name__ == "__main__":

    dot = Digraph(comment='Suffix Tree')

    node_ids = set()
    edges = []
    for line in sys.stdin:
        words = line.strip().split(' ')
        if len(words)!=3:
            continue
        node_ids.add(words[0])
        node_ids.add(words[2])
        edges.append(words)

    for node_id in node_ids:
        dot.node(node_id)

    for edge in edges:
        dot.edge(edge[0], edge[2], label=edge[1])

    dot.render('/tmp/suffix_tree.png', view=True)
