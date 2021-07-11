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
    parser.add_argument('--consider_orientation', action='store_true', default=True,
                        help='If the input fragments are not only from a single strand, activate this option to assign \
                        the fragments to an orientation, which then serves as input to the sequencer.')
    parser.add_argument('--assemble_hamilton', action='store_true', default=False,
                        help='NOTE: CURRENTLY NOT WORKING! Todos: Calculation orientation; Updating the graph after finding the hamilton path; '
                             'Assembles the fragments by building the overlap graph, finding a hamilton path with max '
                             'summed up weight. Then merges all nodes of the path together.')

    parser.add_argument('--assemble_greedy', action='store_true', default=True,
                        help='Assembles the fragments by building the overlap graph and merging the nodes afterward by '
                             'picking the edges with the biggest weight and breaking ties arbitrarily.')

    args = parser.parse_args()

    # Whether and how to print graphs if graphviz is installed
    print_only_result = args.print_only_result
    print_graphs = args.print_graphs

    consider_orientation = args.consider_orientation

    # How to assemble the sequence
    assemble_hamilton = args.assemble_hamilton
    assemble_greedy = args.assemble_greedy

    # The path to file containing the fragments
    path = args.path

    # - - - - - - Initialization - - - - - -

    # Load fragments from file
    fragments = graph.read_fragments(args.path)
    print("Without orientation: ", fragments)
    if consider_orientation:
        print("Calculate good orientation...")
        fragments = graph.get_good_orientation(fragments)
        print("Orientation calculated!")
        print("With orientation: ", fragments)

    # Build overlap graph
    overlap_graph = graph.build_overlap_graph(fragments=fragments)

    # - - - - - - Assembly - - - - - -
    if assemble_greedy:
        assembled_sequence = graph.assembly_greedy(overlap_graph, print_only_result=print_only_result,
                                                   print_graph=print_graphs)
        if isinstance(assembled_sequence, str):
            print(f"Resulting sequence: {assembled_sequence}")
        elif isinstance(assembled_sequence, list):
            print("Could not assemble all fragments!")
            print("Those are left:")
            for seq in assembled_sequence:
                print(f"Node {seq[0]}: ", seq[1]["read"])

    if assemble_hamilton:
        # NOTE: Broken does not work!
        raise NotImplementedError(
            "Work in progress! Todos: Calculation orientation; Updating the graph after finding the hamilton path; ")
        assembled_sequence = graph.assembly_hamilton(overlap_graph, print_only_result=print_only_result,
                                                     print_graph=print_graphs)
        if isinstance(assembled_sequence, str):
            print(f"Resulting sequence: {assembled_sequence}")
        elif isinstance(assembled_sequence, list):
            print("Could not assemble all fragments!")
            print("Those are left:")
            for seq in assembled_sequence:
                print(f"Node {seq.id}: ", seq.value)
