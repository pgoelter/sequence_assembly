from typing import Tuple
from .utils import read_fragments, find_largest_overlaps, overlap


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

    def __init__(self, id: int, source: Vertex, sink: Vertex, weight=None, match=None):
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


class OverlapGraph:
    """Class representing a graph with edges and vertices.
    """

    def __init__(self, edges: list = None, vertices: list = None):
        """Creates an overlap graph object.
        Args:
            edges: A list of elements of type Edge, representing the edges of a graph.
            vertices: A list of elements of type Vertex, representing the vertices of a graph.
        """
        self.edges = edges if edges else []
        self.vertices = vertices if vertices else []

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
                        self.add_edge(source=vertex, sink=_vertex, weight=ov["weight"], match=ov["overlap"])

    def merge(self):
        """Merge vertices, starting with the vertex pairs connected by the edge with the highest weight. Then update all
        edges accordingly. Repeat this step until no more vertices can be merged.

        Returns:
            None
        """
        # todo
        pass

    def add_edge(self, source, sink, weight, match):
        """Adds an edge to the graph.
        Args:
            source: Origin vertex of the edge.
            sink: Target vertex of the edge.
            weight: Weight that can be added to the edge.
            match: In the context of DNA assembly, the overlapping string from both fragments.

        Returns:
            The id of the edge.
        """
        edge_ids = self.get_edge_ids()
        _id = 1 if not edge_ids else max(edge_ids) + 1

        edge = Edge(id=_id, source=source, sink=sink, weight=weight, match=match)
        self.edges.append(edge)
        return _id

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

        return _id

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
