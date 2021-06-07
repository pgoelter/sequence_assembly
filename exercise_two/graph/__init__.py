import random
from typing import List

from graphviz import Digraph

from .utils import read_fragments, find_largest_overlaps, overlap, search_hamilton_path, hamilton, show_graph


class Vertex:
    """Class representing a single vertex of a graph.
    """

    def __init__(self, id, value: str):
        """Construct a vertex by passing an id a fragment and a list of edges.
        Args:
            id: A unique identifier.
            value: The value the vertex holds. in this case e.g. a DNA  (e.g. GATCGTACTGACT).

        """
        self.id = id
        self.value = value

    def set_value(self, value: str):
        """Setter for the attribute Vertex.value.
        Args:
            value: A value as a string.

        Returns:
            None
        """
        self.value = value

    def get_value(self):
        """Return the value of the respective node.
        Returns:
            value
        """
        return self.value

    def get_id(self):
        """Getter returning the id of the vertex.
        Returns:
            Integer representing the id of a vertex.
        """
        return self.id

    def info(self):
        """Returns information about the vertex.
        Returns:
            A tuple containing the id and the value of the vertex.
        """
        return self.id, self.value


class Edge:
    """Class representing a single edge of a graph. An edge is defined by two adjacent vertices.
    """

    def __init__(self, id: int, source: Vertex, sink: Vertex, weight=None, match=None, pos_start=None, pos_end=None):
        """Construct an edge of a directed graph.
        Args:
            id: A unique identifier.
            source: Element of type Vertex. Represents the start node of an edge.
            sink: Element of type Vertex. Represents the end node of an edge.
            weight: A weight assigned to an edge.
            match: Largest overlapping prefix/suffix pair of the values in source and sink. Specific to DNA sequencing.
        """
        self.id = id

        self.source = source
        self.sink = sink

        self.weight = weight

        self.match = match

        # Array indices where the match can be found in the string of the sink (sink contains the matching prefix)
        self.pos_start = pos_start
        self.pos_end = pos_end

    def set_weight(self, weight):
        """Set the weight of an edge.
        Args:
            weight: Weight of the edge.
        Returns:
            None
        """
        self.weight = weight

    def set_match(self, match):
        """Set the attribute match of an edge. In this case this means the overlapping part of the two strings of a
        sink and a source vertex.
        Args:
            match: Largest overlapping prefix/suffix pair of the values in source and sink. Specific to DNA sequencing.

        Returns:
            None
        """
        self.match = match

    def set_match_position(self, match_position):
        self.pos_start = match_position[0]
        self.pos_end = match_position[1]

    def set_sink(self, sink):
        """Set the current sink to another vertex.
        Args:
            sink: A vertex type object.

        Returns:
            None
        """
        self.sink = sink

    def set_source(self, source):
        """Set the current source to another vertex.
        Args:
            source: A vertex type object.

        Returns:
            None
        """
        self.source = source

    def get_weight(self):
        """Getter returning the weight of the edge.
        Returns:
            Integer representing the weight of the edge.
        """
        return self.weight

    def get_sink(self):
        """Getter for returning the sink vertex.
        Returns:
            Vertex representing the sink of the edge.
        """
        return self.sink

    def get_source(self):
        """Getter for returning the source vertex.
        Returns:
            Vertex representing the source of the edge.
        """
        return self.source

    def get_match(self):
        """Getter for returning the overlapping part of two fragments.
        Returns:
            String representing the overlapping slice (prefix) of two fragments.
        """
        return self.match

    def get_match_position(self):
        return self.pos_start, self.pos_start


class OverlapGraph:
    """Class representing a graph with edges and vertices.
    """

    def __init__(self, edges: list = None, vertices: list = None, random_tiebreak: bool = True,
                 print_graph: bool = False, print_only_result: bool = False):
        """Creates an overlap graph object.
        Args:
            edges: A list of elements of type Edge, representing the edges of a graph.
            vertices: A list of elements of type Vertex, representing the vertices of a graph.
        """
        self.edges = edges if edges else []
        self.vertices = vertices if vertices else []
        self.random_tiebreak = random_tiebreak

        self.print_graph = print_graph
        self.print_only_result = print_only_result

    @staticmethod
    def build_from_fragments(fragments):
        """Constructs an overlap graph based on a list of fragments (DNA reads).
        Determines the edges automatically by searching the largest overlapping suffix/prefix pairs for every
        combination of fragments. The weight of the edges is defined by the length of the overlapping string.
        Args:
            fragments: A list of fragments. Each fragment is a string based on the alphabet {A, T, C, G}.

        Returns:
            An OverlapGraph object with edges and vertices.
        """
        graph = OverlapGraph()

        for i, fragment in enumerate(fragments):
            graph.add_vertex(value=fragment)

        graph.determine_edges()
        return graph

    def determine_edges(self):
        """Determines the edges automatically by searching the largest overlapping suffix/prefix pairs for every
        combination of fragments. The weight of the edges is defined by the length of the overlapping string.
        Note:
            String matching uses a naive approach by comparing each suffix and prefix of each string.
        Returns:
            None
        """
        for vertex in self.vertices:
            for _vertex in self.vertices:
                if vertex.id != _vertex.id:
                    ov = overlap(vertex.get_value(), _vertex.get_value())

                    if ov["weight"] != 0:
                        self.add_edge(source=vertex,
                                      sink=_vertex,
                                      weight=ov["weight"],
                                      match=ov["overlap"],
                                      pos_start=ov["prefix_start"],
                                      pos_end=ov["prefix_end"]
                                      )

    @staticmethod
    def merge_fragments(source, sink, slice_positions):
        """Merge two fragments from corresponding vertices representing source and sink of an edge.
        Args:
            source: Origin vertex of an edge.
            sink:  Target vertex of an edge.
            slice_positions: Array indices for the starting and ending array indices of sink.value.

        Returns:
            A string representing the merge fragments of source and sink.
        """
        merged_reads = source.get_value() + sink.get_value()[slice_positions[1]:]

        return merged_reads

    def merge_by_hamiltonian_path(self):
        """Merge nodes by first calculation the hamiltonian paths and pick one with the maximum summed up weight.
        Then merge all nodes together.

        Note:
            CURRENTLY BROKEN. Hamilton path searching does not work!
            # Todo:
            1. Get search for hamiltonian paths done.
            2. Update graph accordingly
            3. Merge nodes of hamiltonian path to get the resulting DNA sequence.

        Returns:
            The result sequence constructed from all given DNA fragments (reads)
        """
        hamiltonian_paths = self.get_all_hamilton_paths()
        for i, g in enumerate(hamiltonian_paths):
            r = self.find_edges_to_path(g)
            dot = Digraph(comment='Overlap Graph')

            # Add vertices to directed graph
            for v in g:
                dot.node(str(v.get_id()), v.get_value())

            # Add edges to directed graph
            for e in r:
                dot.edge(str(e.source.get_id()), str(e.sink.get_id()), label=str(e.weight))

            # Render graph and show it in browser
            dot.render(f'overlap_graph{i}.gv', view=True)
            print()
        # Find hamilton path with max weight. The last element of a path is the summed up weight.
        max_vertices, max_edges = sorted(hamiltonian_paths, key=lambda path: path[0][-1])[-1]
        max_weight = max_vertices[-1]

        all_max_weight_hampaths = [h for h in hamiltonian_paths if h[0][-1] == max_weight]

        def merge_edges(edge_one, edge_two):
            # Todo: Update Graph
            overlap = OverlapGraph.merge_fragments(edge_one.source, edge_one.sink,
                                                   slice_positions=(edge_one.pos_start, edge_one.pos_end))
            ov = utils.overlap(overlap, edge_two.sink.value)

            return Edge(id="", source=Vertex(id="", value=overlap), sink=edge_two.sink, weight=ov["weight"],
                        match=ov["overlap"], pos_start=ov["prefix_start"], pos_end=ov["prefix_end"])

        final_sequence = None

        for _, edges in all_max_weight_hampaths:

            while len(edges) > 1:
                new_edge = merge_edges(edges[0], edges[1])
                edges[0] = new_edge
                edges.remove(edges[1])
                if len(edges) == 2:
                    print()

    def find_edges_to_path(self, path: List[Vertex]):
        """Given a ordered list of vertices that represent a valid path in the graph, search for the according edges.
        Args:
            path: A valid list of vertices representing a sequence in the graph.

        Returns:
            A list of edges.
        """

        def find_edge(source, sink):
            return next((edge for edge in self.edges if edge.sink.id == sink.id and edge.source.id == source.id), None)

        all_edges = []

        for i, e in enumerate(path):
            if i < len(path) - 1:
                all_edges.append(find_edge(path[i], path[i + 1]))

        return all_edges

    def merge_by_arbitrary_tiebreaks(self):
        """Merge vertices, starting with the vertex pairs connected by the edge with the highest weight. Then update all
        edges accordingly. Repeat this step until no more vertices can be merged.

        Returns:
            None
        """
        # Only for graph output
        _p = 0

        if self.print_graph:
            print("Do I Print? ", self.print_graph)
            edges = self.get_edges()
            vertices = self.get_vertices()
            show_graph(edges=edges, vertices=vertices, name=f"graph_{_p}")

        while self.edges:
            max_edge = sorted(self.edges, key=lambda e: e.weight)[-1]

            if self.random_tiebreak:
                #########
                max_edges = [edge for edge in self.edges if edge.weight == max_edge.weight]
                max_edge = random.choice(max_edges)
                #########

            source = max_edge.get_source()
            sink = max_edge.get_sink()

            old_source_id = source.get_id()
            old_sink_id = sink.get_id()

            self.edges.remove(max_edge)

            merged_fragment = OverlapGraph.merge_fragments(source=source, sink=sink,
                                                           slice_positions=(max_edge.pos_start, max_edge.pos_end))

            merged_vertex = self.add_vertex(value=merged_fragment)

            # Delete all outgoing edges from source node
            source_outgoing = self.find_outgoing_edges(source)
            for edge in source_outgoing:
                self.edges.remove(edge)

            # All incoming edges to source node remain except edges from sink node
            edges_to_remove = [edge for edge in self.edges if edge.sink.id == source.id and edge.source.id == sink.id]
            for _edge in edges_to_remove:
                self.edges.remove(_edge)

            # All outgoing edges from sink node remain except for outgoing edge to source node.

            # Delete all incoming edges to sink node
            sink_incoming = self.find_incoming_edges(sink)
            for edge in sink_incoming:
                self.edges.remove(edge)

            self.vertices.remove(source)
            self.vertices.remove(sink)

            # Replace the ids of source node and sink node in all edges with the id of the new node and update
            # the weights if necessary.

            for i, edge in enumerate(self.edges):
                if edge.sink.get_id() == old_sink_id or edge.sink.get_id() == old_source_id:
                    self.edges[i].sink = merged_vertex

                    ov = overlap(self.edges[i].source.get_value(), merged_vertex.get_value())

                    self.edges[i].set_weight(ov["weight"])
                    self.edges[i].set_match(ov["overlap"])
                    self.edges[i].set_match_position((ov["prefix_start"], ov["prefix_end"]))

                if edge.source.get_id() == old_sink_id or edge.source.get_id() == old_source_id:
                    self.edges[i].source = merged_vertex

                    ov = overlap(merged_vertex.get_value(), self.edges[i].sink.get_value())

                    self.edges[i].set_weight(ov["weight"])
                    self.edges[i].set_match(ov["overlap"])
                    self.edges[i].set_match_position((ov["prefix_start"], ov["prefix_end"]))
            _p += 1
            if self.print_only_result or self.print_graph:
                edges = self.get_edges()
                vertices = self.get_vertices()
                show_graph(edges=edges, vertices=vertices, name=f"graph_{_p}")

        if len(self.vertices) == 1:
            return self.vertices[0].value
        elif len(self.vertices) > 1:
            return self.vertices
        else:
            print("Something went wrong!?")

    def find_outgoing_edges(self, vertex, only_id=False):
        """Find all edges starting from a vertex.
        Args:
            vertex: An object of type Vertex.
        Returns:
            A list of objects of type Edge.
        """
        vertex_id = vertex.get_id()
        if only_id:
            return [edge.id for edge in self.edges if edge.get_source().get_id() == vertex_id]
        return [edge for edge in self.edges if edge.get_source().get_id() == vertex_id]

    def find_incoming_edges(self, vertex):
        """Find all edges leading to a vertex.
        Args:
            vertex: An object of type Vertex.
        Returns:
            A list of objects of type Edge.
        """
        vertex_id = vertex.get_id()

        return [edge for edge in self.edges if edge.get_sink().get_id() == vertex_id]

    def add_edge(self, source, sink, weight, match, pos_start, pos_end):
        """Adds an edge to the graph.
        Args:
            source: Origin vertex of the edge.
            sink: Target vertex of the edge.
            weight: Weight that can be added to the edge.
            match: In the context of DNA assembly, the overlapping string from both fragments.
            pos_start: Array index prefix start in sink.value.
            pos_end: Array index prefix end in sink.value.
        Returns:
            The id of the edge.
        """
        edge_ids = self.get_edge_ids()
        _id = 1 if not edge_ids else max(edge_ids) + 1

        edge = Edge(id=_id, source=source, sink=sink, weight=weight, match=match, pos_start=pos_start, pos_end=pos_end)
        self.edges.append(edge)
        return edge

    def add_vertex(self, value):
        """Adds a vertex to the graph.
        Args:
            value: A value associated with the vertex.

        Returns:
            The id of the vertex.
        """
        vertex_ids = self.get_vertex_ids()
        _id = 1 if not vertex_ids else max(vertex_ids) + 1

        vertex = Vertex(id=_id, value=value)
        self.vertices.append(vertex)

        return vertex

    def remove_edge(self, id):
        """Removes an edge with the given id from the graph.
        Args:
            id: Id of the edge to remove.

        Returns:
            True when removing was successful, False otherwise.
        """
        try:
            self.edges.remove(self.get_edge(id))
            return True
        except ValueError:
            return False

    def remove_vertex(self, id):
        """Removes an vertex with the given id from the graph.
        Args:
            id: Id of the vertex to remove.

        Returns:
            True when removing was successful, False otherwise.
        """
        try:
            self.edges.remove(self.get_edge(id))
            return True
        except ValueError:
            return False

    def set_random(self, val: bool):
        """Setter. If this value is set to true the choice made in method self.merge_by_arbitrary_tiebreaks is random
        among those edges with the same max weight.
        Args:
            val: Boolean value

        Returns:
            None
        """
        self.random_tiebreak = val

    def set_print(self, val: bool):
        """Setter. If this value is set to true the graph output is printed with graphviz.
        Args:
            val: Boolean value

        Returns:
            None
        """
        self.print_graph = val

    def set_print_result_only(self, val: bool):
        """Setter. If this value is set to true only the resulting graph gets printed..
        Args:
            val: Boolean value

        Returns:
            None
        """
        self.print_only_result = val

    def get_edge(self, id):
        """Returns the edge with a given id.
        Args:
            id: Id of the edge.

        Returns:
            An object of type Edge. Returns None if id could not be found.
        """
        return next((edge for edge in self.edges if edge.id == id), None)

    def get_vertex(self, id):
        """Returns the vertex with a given id.
        Args:
            id: Id of the vertex.

        Returns:
            An object of type Vertex. Returns None if id could not be found.
        """
        return next((edge for edge in self.edges if edge.id == id), None)

    def get_vertex_ids(self):
        """Get a list of ids of all existing vertices.
        Returns:
            List of vertex ids.
        """
        if self.vertices:
            return [v.id for v in self.vertices]
        return None

    def get_edge_ids(self):
        """Get a list of ids of all existing edges.
        Returns:
            List of edge ids.
        """
        if self.edges:
            return [e.id for e in self.edges]
        return None

    def get_edges(self):
        """Get all existing edges.
        Returns:
            A list of elements of type Edge.
        """
        return self.edges

    def get_vertices(self):
        """Get all existing vertices.
        Returns:
            A list of elements of type Vertex.
        """
        return self.vertices

    def get_fragments(self):
        """Get all existing fragments.
        Returns:
            A list of strings whereas each string represents a read (fragment).
        """
        return [vertex.get_value() for vertex in self.vertices]

    def get_all_hamilton_paths(self):
        """Find all hamilton paths in the graph.
        Returns:
            A list of hamilton paths, whereas each element in the graph is a single vertex.
        """
        hamiltonian_paths = []
        vertex_count = len(self.vertices)
        for vertex in self.vertices:
            res = hamilton(self, vertex_count, vertex)

            hamiltonian_paths.append(res)
        # Remove None values
        return [h for h in hamiltonian_paths if h is not None]
