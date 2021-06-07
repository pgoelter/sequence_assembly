import graph

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description='Assemble DNA fragments from a *.frag file.')

    parser.add_argument('path', type=str,
                        help='The path containing the DNA fragments (reads).')

    parser.add_argument('--print_graphs', action='store_true', default=False,
                        help='Print all graphs in between each merge of two nodes.')

    parser.add_argument('--print_only_result',
                        action='store_true', default=False,
                        help='Prints only the resulting graph. Should be a single node if everything worked.')

    parser.add_argument('--assemble_greedy', action='store_true', default=True,
                        help='Assembles the fragments by building the overlap graph and merging the nodes afterward by '
                             'picking the edges with the biggest weight and breaking ties arbitrarily.')

    parser.add_argument('--assemble_hamilton', action='store_true', default=False,
                        help='NOTE: CURRENTLY NOT WORKING AS INTENDED! '
                             'Assembles the fragments by building the overlap graph, finding a hamilton path with max '
                             'summed up weight. Then merges all nodes of the path together.')

    parser.add_argument('--random', action='store_true', default=False,
                        help='When choosing option assemble_greedy this option can be turned on so each time for '
                             'choosing nodes to merge a random edge with maximum weight gets picked if there are more '
                             'than one edges with the same weight among those with the highest weight in the graph.')

    args = parser.parse_args()

    # Whether and how to print graphs if graphviz is installed
    print_only_result = args.print_only_result
    print_graphs = args.print_graphs

    # How to assemble the sequence
    assemble_hamilton = args.assemble_hamilton
    assemble_greedy = args.assemble_greedy

    # Whether to make a random choice which edges with maximum weight should get picked.
    random = args.random

    # The path to file containing the fragments
    path = args.path

    # - - - - - - Initialization - - - - - -

    # Load fragments from file
    fragments = graph.read_fragments(args.path)

    # Build overlap graph
    overlap_graph = graph.OverlapGraph.build_from_fragments(fragments=fragments)

    # Only works with a valid installation of graphviz (see: https://graphviz.org/download/)
    overlap_graph.set_print(print_graphs)
    overlap_graph.set_print_result_only(print_only_result)


    overlap_graph.set_random(random)

    # - - - - - - Assembly - - - - - -
    if assemble_greedy:
        assembled_sequence = overlap_graph.merge_by_arbitrary_tiebreaks()
        if isinstance(assembled_sequence, str):
            print(f"Resulting sequence: {assembled_sequence}")
        elif isinstance(assembled_sequence, list):
            print("Could not assemble all fragments!")
            print("Those are left:")
            for seq in assembled_sequence:
                print(f"Node {seq.id}: ", seq.value)

    if assemble_hamilton:
        # NOTE: Broken does not work!
        assembled_sequence = overlap_graph.merge_by_hamiltonian_path()
        if isinstance(assembled_sequence, str):
            print(f"Resulting sequence: {assembled_sequence}")
        elif isinstance(assembled_sequence, list):
            print("Could not assemble all fragments!")
            print("Those are left:")
            for seq in assembled_sequence:
                print(f"Node {seq.id}: ", seq.value)
