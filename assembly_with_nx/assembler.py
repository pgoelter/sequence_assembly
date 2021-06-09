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

    parser.add_argument('--assemble_hamilton', action='store_true', default=False,
                        help='NOTE: CURRENTLY NOT WORKING AS INTENDED! '
                             'Assembles the fragments by building the overlap graph, finding a hamilton path with max '
                             'summed up weight. Then merges all nodes of the path together.')

    parser.add_argument('--assemble_greedy', action='store_true', default=True,
                        help='Assembles the fragments by building the overlap graph and merging the nodes afterward by '
                             'picking the edges with the biggest weight and breaking ties arbitrarily.')

    args = parser.parse_args()

    # Whether and how to print graphs if graphviz is installed
    print_only_result = args.print_only_result
    print_graphs = args.print_graphs

    # How to assemble the sequence
    assemble_hamilton = args.assemble_hamilton
    assemble_greedy = args.assemble_greedy

    # The path to file containing the fragments
    path = args.path

    # - - - - - - Initialization - - - - - -

    # Load fragments from file
    fragments = graph.read_fragments(args.path)

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
        raise NotImplementedError("Work in progress! Todos: Calculation orientation; Updating the graph after finding the hamilton path; ")
        assembled_sequence = graph.assembly_hamilton(overlap_graph, print_only_result=print_only_result,
                                                     print_graph=print_graphs)
        if isinstance(assembled_sequence, str):
            print(f"Resulting sequence: {assembled_sequence}")
        elif isinstance(assembled_sequence, list):
            print("Could not assemble all fragments!")
            print("Those are left:")
            for seq in assembled_sequence:
                print(f"Node {seq.id}: ", seq.value)