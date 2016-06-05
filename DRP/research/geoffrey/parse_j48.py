#!/usr/bin/env python

import json
import re
import os
import sys

re_head = re.compile("J48 (un)?pruned tree")
re_divider_line = re.compile("^-*\n$")
re_blank_line = re.compile("^[ \t\n]*$")
re_splitter = re.compile("[ :]")
re_range = re.compile(
    r"^'\("
    r"(-inf|-?[0-9]+(\.[0-9]+)?)"
    r"-"
    r"(-?[0-9]+(\.[0-9]+)?\]|inf\))"
    r"'$")


def parse_value(token):
    """Returns an float if the token represents a number, a range if the token
    represents a range of numbers, otherwise return the token as is."""
    try:
        return float(token)
    except ValueError:
        # Look for ranges of the form '(start-end]', ' included
        if re_range.match(token):
            range_str = token[2:-2]

            # Careful not to use a minus sign as a dash.
            separator_dash = range_str.find("-", 1)
            return (parse_value(range_str[:separator_dash]),
                    parse_value(range_str[separator_dash + 1:]))
        else:
            # Not a number or range - so it must be nominal, leave it as it.
            return token


def parse_count(count_str):
    content = count_str[1:-1]
    vs = content.split('/')
    if len(vs) == 1:
        return float(vs[0]), 0.0
    else:
        return float(vs[0]), float(vs[1])


def parse_line(line):
    """Split the line into a tuple
    (depth, feature, comparator, value, classification/None)"""
    # Avoid empty strings from double whitespaces and the likes.
    split = [l for l in re_splitter.split(line) if len(l) > 0]
    depth = 0
    for part in split:
        if part == "|":
            depth += 1
        else:
            break
    has_classif = ':' in line
    if not has_classif:
        return (depth, split[depth], split[depth + 1],
                parse_value(split[depth + 2]),
                None,
                None)
    else:
        no_nonsense = line.replace('|   ', '')
        fcv = no_nonsense[:no_nonsense.find(':')]
        fcv = fcv.split()
        classif = no_nonsense[no_nonsense.find(':') + 1:no_nonsense.find('(')].strip()
        return (depth, fcv[0], fcv[1], " ".join(fcv[2:]), classif,
                parse_count(split[-1]))


def parse_tree(lines):
    """Parses input lines into a decision tree."""
    current_index = [0]  # need mutable container because of closure limitations

    def parse(current_depth):
        """Helper recursive function."""
        node_feature = None
        children = []
        while current_index[0] < len(lines):
            line = lines[current_index[0]]
            depth, feature, comparator, value, classif, count = parse_line(line)
            if depth < current_depth:
                # Finished parsing this node.
                break
            elif depth == current_depth:
                if node_feature is None:
                    node_feature = feature
                elif node_feature != feature:
                    raise Exception("Error : Feature mismatch - expected %s"
                                    "but got : \n%s"
                                    % (node_feature, line))

                # Another branch
                current_index[0] += 1
                if classif is None:
                    children.append((comparator, value,
                                     parse(current_depth + 1)))
                else:
                    children.append((comparator, value, classif, count))
            else:
                sys.stderr.write("Parse line output: %s\n" % ((depth, feature, comparator, value, classif, count),))
                sys.stderr.write("Depth: %s\n" % depth)
                sys.stderr.write("Current Depth: %s\n" % current_depth)
                raise Exception("Error : Input jumps two levels at once\n%s."
                                % line)

        return (node_feature, children)

    return parse(0)


def get_tree_lines(lines):
    """Return the lines of the input that correspond to the decision tree."""
    tree_lines = []
    for i in range(len(lines) - 2):
        if re_head.match(lines[i]):
            assert (re_divider_line.match(lines[i + 1]) and
                    re_blank_line.match(lines[i + 2])), \
                "Input not in expected format."
            for l in lines[i + 3:]:
                if re_blank_line.match(l):
                    return tree_lines
                else:
                    tree_lines.append(l[:-1])  # remove newline at the end

    raise Exception("Error : Failed to find tree in input.")


def main(argv):
    if argv:
        input_filename = argv[0]
        if os.path.isfile(input_filename):
            f = open(input_filename)
            lines = f.readlines()
            f.close()
        else:
            raise Exception("Error : File %s not found!" % input_filename)
    else:
        lines = sys.stdin.readlines()

    if not lines:
        raise Exception("Error : Empty input!")

    tree_lines = get_tree_lines(lines)
    tree = parse_tree(tree_lines)
    print json.dumps(tree)


main(sys.argv[1:])
