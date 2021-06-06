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