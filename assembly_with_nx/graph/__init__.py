import copy

import networkx as nx

from .utils import overlap, show_graph, read_fragments, hamilton


def build_overlap_graph(fragments):
    """Constructs an overlap graph based on a list of fragments (DNA reads).
    Determines the edges automatically by searching the largest overlapping suffix/prefix pairs for every
    combination of fragments. The weight of the edges is defined by the length of the overlapping string.
    Args:
        fragments: A list of fragments. Each fragment is a string based on the alphabet {A, T, C, G}.

    Returns:
        An OverlapGraph object with edges and vertices.
    """
    overlap_graph = nx.DiGraph()
    for i, fragment in enumerate(fragments):
        for j, _fragment in enumerate(fragments):
            if i != j:
                overlap_info = overlap(fragment, _fragment)
                overlap_exists = overlap_info["weight"] != 0
                if overlap_exists:
                    overlap_graph.add_node(i, read=fragment)
                    overlap_graph.add_node(j, read=_fragment)

                    overlap_graph.add_edge(i,
                                           j,
                                           weight=overlap_info["weight"],
                                           match=overlap_info["overlap"],
                                           pos_start=overlap_info["prefix_start"],
                                           pos_end=overlap_info["prefix_end"])

    return overlap_graph


def merge_fragments(overlap_graph: nx.DiGraph, source_id: int, sink_id: int, slice_positions):
    """Merge two fragments from corresponding vertices representing source and sink of an edge.
    Args:
        overlap_graph: The graph containing the nodes
        source_id: Origin node id of an edge.
        sink_id:  Target node of id an edge.
        slice_positions: Array indices for the starting and ending array indices of sink.value.

    Returns:
        A string representing the merge fragments of source and sink.
    """
    merged_reads = overlap_graph.nodes[source_id]["read"] + overlap_graph.nodes[sink_id]["read"][slice_positions[1]:]

    return merged_reads


def assembly_greedy(overlap_graph: nx.DiGraph, print_only_result=False, print_graph=False):
    """Merge vertices, starting with the vertex pairs connected by the edge with the highest weight. Then update all
    edges accordingly. Repeat this step until no more vertices can be merged.
    Args:
        overlap_graph: The graph in which to assemble the fragments to a single sequence if possible.
        print_only_result: If set to True prints the resulting graph and stores it as a .pdf file.
        print_graph: If set to True prints all graphs and stores them as pdf files.
    Returns:
        The resulting sequence or a list of sequences left if not everything could be assembled.
    """
    # Only for graph output
    _p = 0
    if print_only_result or print_graph:
        show_graph(edges=list(overlap_graph.edges(data=True)), vertices=list(overlap_graph.nodes(data=True)),
                   name=f"graph_{_p}")

    while overlap_graph.number_of_edges() != 0:
        max_edge = sorted(overlap_graph.edges(data=True), key=lambda e: e[2]["weight"])[-1]

        old_source = max_edge[0]
        old_sink = max_edge[1]

        merged_fragment = merge_fragments(overlap_graph, source_id=old_source, sink_id=old_sink,
                                          slice_positions=(max_edge[2]["pos_start"], max_edge[2]["pos_end"]))
        _id_merged = max(overlap_graph.nodes()) + 1
        merged_node = overlap_graph.add_node(
            _id_merged,
            read=merged_fragment
        )

        # Delete all outgoing edges from source
        outgoing_edges = list(overlap_graph.out_edges(nbunch=old_source))
        for source_node, sink_node in outgoing_edges:
            overlap_graph.remove_edge(u=source_node, v=sink_node)

        # If an edge from sink to source exists delete it
        if overlap_graph.has_edge(old_sink, old_source):
            overlap_graph.remove_edge(old_sink, old_source)

        # All outgoing edges from sink node remain.

        # Delete all incoming edges to sink node
        incoming_edges = list(overlap_graph.in_edges(nbunch=old_sink))
        for source_node, sink_node in incoming_edges:
            overlap_graph.remove_edge(u=source_node, v=sink_node)

        # Replace the ids of source node and sink node in all edges with the id of the new node and update
        # the weights if necessary.

        # Update the graph
        for edge in list(overlap_graph.edges(data=True)):
            if edge[1] == old_sink or edge[1] == old_source:
                ov = overlap(overlap_graph.nodes[edge[0]]["read"], overlap_graph.nodes[_id_merged]["read"])

                overlap_graph.remove_edge(edge[0], edge[1])
                overlap_graph.add_edge(edge[0],
                                       _id_merged,
                                       weight=ov["weight"],
                                       match=ov["overlap"],
                                       pos_start=ov["prefix_start"],
                                       pos_end=ov["prefix_end"])

            if edge[0] == old_sink or edge[0] == old_source:
                ov = overlap(overlap_graph.nodes[_id_merged]["read"], overlap_graph.nodes[edge[1]]["read"])

                overlap_graph.remove_edge(edge[0], edge[1])
                overlap_graph.add_edge(_id_merged,
                                       edge[1],
                                       weight=ov["weight"],
                                       match=ov["overlap"],
                                       pos_start=ov["prefix_start"],
                                       pos_end=ov["prefix_end"])

        # Delete sink and source from graph as a new node will be added to replace both
        overlap_graph.remove_node(old_source)
        overlap_graph.remove_node(old_sink)

        _p += 1
        if print_graph:
            show_graph(edges=list(overlap_graph.edges(data=True)), vertices=list(overlap_graph.nodes(data=True)),
                       name=f"graph_{_p}")

    if overlap_graph.number_of_nodes() == 1:
        return list(overlap_graph.nodes(data=True))[0][1]["read"]
    elif overlap_graph.number_of_nodes() > 1:

        return list(overlap_graph.nodes(data=True))
    else:
        raise ValueError("No nodes in the Graph left. That shouldn't happen!")


def assembly_hamilton(overlap_graph: nx.DiGraph, print_only_result=False, print_graph=False):
    """Merge nodes by first calculation the hamiltonian paths and pick one with the maximum summed up weight.
           Then merge all nodes together.
    Args:
        overlap_graph: The graph represented as an object of type networkx.DiGraph
        print_only_result: If True show resulting graph as output and save as .pdf file.
        print_graph: If true show all resulting graphs for each step where a node is merged and save as .pdf file.

    Note:
       1. Get search for hamiltonian paths done.
       2. Update graph accordingly
       3. Merge nodes of hamiltonian path to get the resulting DNA sequence.

    Returns:
        The result sequence constructed from all given DNA fragments (reads)
           """
    hamiltonian_paths = []

    max_node_count = overlap_graph.number_of_nodes()

    for node in list(overlap_graph.nodes()):
        none_or_hampath = hamilton(overlap_graph, max_node_count, node, path=[node], visited=[node])
        hamiltonian_paths.append(none_or_hampath)

    # Clean from None types
    hamiltonian_paths = [h for h in hamiltonian_paths if h is not None]

    hampaths_sorted_by_weight = []

    for ham_nr, ham in enumerate(hamiltonian_paths):
        tmp = {"summed_weight": 0,
               "hampath_index": ham_nr,
               "node_pairs": []}
        for i, h in enumerate(ham):
            if i < max_node_count - 1:
                tmp["summed_weight"] += overlap_graph[ham[i]][ham[i + 1]]["weight"]
                tmp["node_pairs"].append((ham[i], ham[i + 1]))
        hampaths_sorted_by_weight.append(tmp)

    hampaths_sorted_by_weight = sorted(hampaths_sorted_by_weight, key=lambda element: element["summed_weight"])

    # todo: Orientation

    max_path = hampaths_sorted_by_weight[-1]

    all_edges = copy.deepcopy(list(overlap_graph.edges()))

    for edge in all_edges:
        if edge not in max_path["node_pairs"]:
            overlap_graph.remove_edge(edge[0], edge[1])

    # FIXME: mergin wont work yet
    def merge_edges(e1, e2):
        data = overlap_graph.get_edge_data(e1[0], e1[1])

        _overlap = merge_fragments(overlap_graph, e1[0], e1[1],
                                   slice_positions=(data["pos_start"], data["pos_end"]))

        ov = overlap(_overlap, overlap_graph.nodes[e2[0]]["read"])
        _id_merged = overlap_graph.number_of_nodes() + 1
        overlap_graph.add_node(_id_merged, read=_overlap)

        overlap_graph.remove_edge(e1[0], e1[1])
        overlap_graph.remove_node(e1[0])
        overlap_graph.remove_node(e1[1])

        overlap_graph.add_edge(_id_merged,
                               e2[1],
                               weight=ov["weight"],
                               match=ov["overlap"],
                               pos_start=ov["prefix_start"],
                               pos_end=ov["prefix_end"])

    while len(graph.edges()) > 1:
        edges = list(overlap_graph.edges())
        new_edge = merge_edges(edges[0], edges[1])
        edges[0] = new_edge
        edges.remove(edges[1])


if __name__ == "__main__":
    reads = read_fragments("../data/frag.dat")
    graph = build_overlap_graph(reads)
    show_graph(edges=list(graph.edges(data=True)), vertices=list(graph.nodes(data=True)), name="dies")
    assembly_hamilton(graph)
    show_graph(edges=list(graph.edges(data=True)), vertices=list(graph.nodes(data=True)), name="das")

    print()
