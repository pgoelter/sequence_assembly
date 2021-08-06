import copy

import networkx as nx
from graphviz import Digraph

COMPLEMENTS = {
    "A": "T",
    "G": "C",
    "C": "G",
    "T": "A"
}


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


def weight(fragment_one, fragment_two):
    """Returns the weight of edge from fragment_one to fragment_two or the length of the longest suffix from
    fragment_one, which is also prefix to fragment_two.
    Args:
        fragment_one: Fragment to read.
        fragment_two: Fragment to read.
    Returns:
        The weight of an edge defined by two fragments.
    """
    return overlap(fragment_one, fragment_two)["weight"]


def complement(fragment: str):
    """Returns the complement of a given fragment.
    Args:
        fragment: Fragment to read.

    Returns:
        Complement of fragment.
    """
    c = reversed([COMPLEMENTS[l] for l in list(fragment)])
    return "".join(c)


def same(fragment_one: str, fragment_two: str):
    """Largest possible edge weight if fragment_one and fragment_two belong to one orientation.
    Args:
        fragment_one: Fragment to read.
        fragment_two: Fragment to read.

    Returns:
        Largest possible edge weight.
    """
    return max(weight(fragment_one, fragment_two), weight(fragment_two, fragment_one))


def opp(fragment_one: str, fragment_two: str):
    """Largest possible edge weight if fragment_one and complement of fragment_two belong to one orientation.
    Args:
        fragment_one: Fragment to read.
        fragment_two: Fragment to read.

    Returns:
        Largest possible edge weight.
    """
    c_two = complement(fragment_two)

    return max(weight(fragment_one, c_two), weight(c_two, fragment_one))


def _get_orientations_input(fragments: list):
    combinations = []
    for f in fragments:
        _copy_fragments = copy.deepcopy(fragments)
        _copy_fragments.remove(f)
        _tuple = (f, _copy_fragments)
        combinations.append(_tuple)
    return combinations


def get_good_orientation(fragments: list):
    to_calculate = _get_orientations_input(fragments)

    orientations_with_weight = []

    for _input in to_calculate:
        orientation, _weight = calc_orientation(_input[0], _input[1])
        orientations_with_weight.append((orientation, _weight))

    max_weight_orientation, _w = sorted(orientations_with_weight, key=lambda e: e[1])[-1]

    return max_weight_orientation


def calc_orientation(start_fragment: str, fragments: list):
    """Calculate an orientation from a given start fragment and a list of fragments.
    Args:
        start_fragment: Fragment to start with (first element in set O)
        fragments: Fragments to check.
    Returns:
        Set of fragments which represents an Orientation.
    """
    O = [start_fragment]
    others = fragments

    while others:
        sum_same = 0
        sum_opp = 0
        fragment_to_test = others.pop(0)
        for fragment in O:
            sum_same += same(fragment, fragment_to_test)
            sum_same += same(fragment_to_test, fragment)

            sum_opp += opp(fragment, fragment_to_test)
            sum_opp += opp(fragment_to_test, fragment)
        if sum_same < sum_opp:
            O.append(complement(fragment_to_test))
        elif sum_opp <= sum_same:
            O.append(fragment_to_test)
    return O, sum_same + sum_opp


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


def hamilton(graph: nx.DiGraph, vertex_count, start_vertex, path=None, visited=None):
    """Find a hamilton path in a given graph from a given starting vertex. Search by traversing the nodes.
    Args:
        graph: The graph to search in
        vertex_count: The total amount of vertices in the graph.
        start_vertex: The vertex to start the search from.
        path: Variable to hold the path of vertices which form a potential hamiltonian path.
        visited: A list of vertices that were already visited.

    Returns:
        A list of lists, whereas each element is a list of Vertices which form a hamilton path.
    """
    if visited is None:
        visited = []

    if path is None:
        path = []

    # if start_vertex not in path:
    #     path.append(start_vertex)

    if len(path) == vertex_count:
        return path

    sucessor_vertices = graph.successors(start_vertex)

    for successor in sucessor_vertices:
        if successor not in set(visited):
            # Then visit
            visited.append(successor)
            path.append(successor)
            candidate = hamilton(graph, vertex_count, successor, path, visited)

            if candidate:
                return candidate

            visited.pop()
            path.pop()


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
        dot.node(str(v[0]), v[1]["read"])

    # Add edges to directed graph
    for i, e in enumerate(edges):
        dot.edge(str(e[0]), str(e[1]), label=f"{str(e[2]['weight'])}: {e[2]['match']}")

    # Render graph and show it in browser
    dot.render(name, view=True)
