from typing import Tuple
from .utils import read_fragments, find_largest_overlaps, overlap


class Vertex:
    """Class representing a single vertex of a graph.
    """

    def __init__(self, id, value: str):
        """Construct a vertex by passing an id a fragment and a list of edges.
        Args:
            id: A unique identifier.
            value: The value the nod holds. in this casse e.g. a DNA  (e.g. GATCGTACTGACT).
            edges: A list of elements of type Edge.
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

        """
        return self.id, self.value


class Edge:
    """Class representing a single edge of a graph. An edge is defined by two adjacent vertices.
    """

    def __init__(self, id: int, source: Vertex, sink: Vertex, weight=None, match=None):
        """Construct an edge by passing an id and a tuple containing two vertices of type Vertex.
        Args:
            id: A unique identifier.
            adjacent_vertices: Two elements of type Vertex. These to vertices must be adjacent as they define the edge.
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
        self.match = match

    def get_weight(self):
        """Getter returning the weight of the edge.
        Returns:
            Integer representing the weight of the edge.
        """

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
        self.edges = edges if edges else []
        self.vertices = vertices if vertices else []

    @staticmethod
    def build_from_fragments(fragments):
        """
        Args:
            fragments:

        Returns:

        """
        graph = OverlapGraph()

        for i, fragment in enumerate(fragments):
            graph.add_vertex(value=fragment)

        graph.determine_edges()
        return graph

    def determine_edges(self):
        """Determine the edges with respective weight. Based on the fragments every vertex holds.
        For every fragment find the largest overlap to another node and add an edge if an overlap exists.
        Returns:

        """
        for vertex in self.vertices:
            for _vertex in self.vertices:
                if vertex.id != _vertex.id:
                    ov = overlap(vertex.get_value(), _vertex.get_value())

                    if ov["weight"] != 0:
                        self.add_edge(source=vertex, sink=_vertex, weight=ov["weight"], match=ov["overlap"])

    def merge(self):
        pass

    def add_edge(self, source, sink, weight, match):
        edge_ids = self.get_edge_ids()
        _id = 1 if not edge_ids else max(edge_ids) + 1

        edge = Edge(id=_id, source=source, sink=sink, weight=weight, match=match)
        self.edges.append(edge)

    def add_vertex(self, value):
        vertex_ids = self.get_vertex_ids()
        _id = 1 if not vertex_ids else max(vertex_ids) + 1

        vertex = Vertex(id=_id, value=value)
        self.vertices.append(vertex)

    def get_vertex_ids(self):
        if self.vertices:
            return [v.id for v in self.vertices]
        return None

    def get_edge_ids(self):
        if self.edges:
            return [e.id for e in self.edges]
        return None

    def get_edges(self):
        return self.edges

    def get_vertices(self):
        return self.vertices

    def get_fragments(self):
        return [vertex.get_value() for vertex in self.vertices]
