from graphviz import Digraph


def read_fragments(filename: str):
    """Reads fragments from a given file.
    Args:
        filename: The name of the file containing the fragments, separated by linebreaks.
    Returns:
        fragments: A list of strings whereas every string is a fragment.
    """
    with open(filename, "r") as fd:
        # Read fragments and remove linebreaks from string
        fragments = [frag.strip() for frag in fd.readlines()]
    return fragments


def is_suffix(suffix: str, word: str):
    """Checks whether suffix is an actual suffix of a word.
    Args:
        suffix: String that potentially is a suffix.
        word: String which may ends with the suffix.

    Returns:
        True if suffix, False otherwise.
    """
    return word.endswith(suffix)


def is_prefix(prefix: str, word: str):
    """Checks whether prefix is an actual prefix of a word.
    Args:
        prefix: String that potentially is a prefix.
        word: String which may start with the prefix.

    Returns:
        True if prefix, False otherwise.
    """
    return word.startswith(prefix)


def overlap(string_one: str, string_two: str):
    """Returns the largest prefix of string two that overlaps with a respective suffix from string one.
    Args:
        string_one: String providing suffixes for comparison.
        string_two: String providing prefixes for comparison.

    Returns:
        Dictionary with information about the overlapping prefix.
    """
    len_s_one = len(string_one)
    len_s_two = len(string_two)

    # Content of the tuple => (start_index, stop_index, string from start_index to stop_index)
    largest_overlap = {"suffix_string": string_one,
                       "prefix_string": string_two,
                       "overlap": None,
                       "weight": 0}
    len_overlap = 0
    # Iterate over every suffix from string one from right to left
    for i in range(len_s_one, 0, -1):
        # Current suffix from string one
        current_suffix_one = string_one[i - 1:len_s_one]

        # Iterate over every prefix of string two from left to right
        for j in range(len_s_two + 1):
            current_prefix_two = string_two[:j]
            if current_suffix_one == current_prefix_two:
                tmp_len = len(current_prefix_two)
                if tmp_len > len_overlap:
                    len_overlap = tmp_len

                    largest_overlap["suffix_start"] = i - 1
                    largest_overlap["suffix_end"] = len_s_one

                    largest_overlap["prefix_start"] = 0
                    largest_overlap["prefix_end"] = j

                    largest_overlap["weight"] = tmp_len
                    largest_overlap["overlap"] = current_suffix_one
    return largest_overlap


def find_largest_overlaps(string_one, string_two):
    """Searches for overlaps between two strings.
    1. Compare all suffixes of string one with all prefixes of string two.
    2. Compare all suffixes of string two with all prefixes of string one.
    """
    two_overlaps_suffix_one = overlap(string_one, string_two)
    one_overlaps_suffix_two = overlap(string_two, string_one)

    return two_overlaps_suffix_one, one_overlaps_suffix_two


def search_hamilton_path(graph, vertex_count, start_vertex, path=None, edges_visited=None, path_weight=None):
    """Find a hamilton path in a given graph from a given starting vertex. Search by traversing the edges of the graph
    Args:
        graph: The graph to search in
        vertex_count: The total amount of vertices in the graph.
        start_vertex: The vertex to start the search from.
        path: Variable to hold the path of vertices which form a potential hamiltonian path.
        edges_visited: A list of edges already visited.
        path_weight: The summed up path weight of a path.

    Note: # TODO
        Currently Broken! Does not work as intended!
    Returns:
        A list of lists, whereas each element is a list of Vertices which form a hamilton path.
    """
    if edges_visited is None:
        edges_visited = []

    if path is None:
        path = []

    if start_vertex not in path:
        path.append(start_vertex)

    if path_weight is None:
        path_weight = 0

    if len(path) == vertex_count:
        path.append(path_weight)
        return path, edges_visited

    outgoing_edges = graph.find_outgoing_edges(start_vertex)

    for visited in edges_visited:
        if visited in outgoing_edges:
            outgoing_edges.remove(visited)

    if outgoing_edges:
        for edge in outgoing_edges:
            if edge not in edges_visited and edge.sink not in path:

                edges_visited.append(edge)

                path_weight += edge.weight

                candidate = search_hamilton_path(graph, vertex_count, edge.sink, path, edges_visited, path_weight)

                if candidate:
                    return candidate
    else:

        print("Dead end!")


def hamilton(graph, vertex_count, start_vertex, path=None, visited=None):
    """Find a hamilton path in a given graph from a given starting vertex. Search by traversing the nodes.
    Args:
        graph: The graph to search in
        vertex_count: The total amount of vertices in the graph.
        start_vertex: The vertex to start the search from.
        path: Variable to hold the path of vertices which form a potential hamiltonian path.
        visited: A list of vertices that were already visited.

    Note: # TODO
        Currently Broken! Does not work as intended!

    Returns:
        A list of lists, whereas each element is a list of Vertices which form a hamilton path.
    """
    if visited is None:
        visited = []

    if path is None:
        path = []

    if start_vertex not in path:
        path.append(start_vertex)

    if len(path) == vertex_count:
        return path

    sucessor_vertices = [edge.sink for edge in graph.find_outgoing_edges(start_vertex)]

    if sucessor_vertices:
        for successor in sucessor_vertices:
            if successor not in set(visited):
                # Then visit
                visited.append(successor)

                res_path = [element for element in path]  #

                candidate = hamilton(graph, vertex_count, successor, res_path, visited)

                if candidate:
                    return candidate
                print("Dead end!")


def show_graph(edges, vertices, name):
    """Prints the given graph out by utilizing the python graphviz interface.
    Args:
        edges: A list of objects of type Edge.
        vertices: A list of objects of type Vertex.
        name: The name of the file where the graph gets stored.

    Returns:
        None
    """
    dot = Digraph(comment=name)

    # Add vertices to directed graph
    for v in vertices:
        dot.node(str(v.get_id()), v.get_value())

    # Add edges to directed graph
    for e in edges:
        dot.edge(str(e.source.get_id()), str(e.sink.get_id()), label=str(e.weight))

    # Render graph and show it in browser
    dot.render(name, view=True)
