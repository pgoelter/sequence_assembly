import graph
from graphviz import Digraph

if __name__ == "__main__":
    F = graph.read_fragments("data/frag.dat")
    overlap_graph = graph.OverlapGraph.build_from_fragments(F)

    # Get edges and vertices for plotting
    edges = overlap_graph.get_edges()
    vertices = overlap_graph.get_vertices()

    # Create directed graph (Requires a valid installation of graphviz)
    dot = Digraph(comment='Overlap Graph')

    # Add vertices to directed graph
    for v in vertices:
        dot.node(str(v.get_id()), v.get_value())

    # Add edges to directed graph
    for e in edges:
        dot.edge(str(e.source.get_id()), str(e.sink.get_id()), label=str(e.weight))

    # Render graph and show it in browser
    dot.render('overlap_graph.gv', view=True)