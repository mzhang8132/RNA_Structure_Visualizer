import copy
import matplotlib.pyplot as plt
from tqdm import trange
import networkx as nx
import plotly.graph_objects as go
import csv
import argparse

pairs = [("A", "U"), ("U", "A"), ("C", "G"), ("G", "C")]

"""
Implementation of the nussinov algorithm with the recurrence as depicted in class

Input:
    seq: RNA sequence

Returns:
    filled out Nussinov DP table
"""
def nussinov(seq):
    s = [[0 for j in range(len(seq))] for i in range(len(seq))]

    for j in range(len(seq)):
        for i in range(len(seq)-1, -1, -1):
            if i >= j:
                continue
            else:
                temp_scores = {"redundant": s[i+1][j-1], 
                            "i unpaired": s[i+1][j], 
                            "j unpaired": s[i][j-1]}

            if (seq[i], seq[j]) in pairs:
                temp_scores["pair"] = s[i+1][j-1] + 1

            bifurcations = sorted({k: s[i][k] + s[k+1][j] for k in range(i+1, j) if k > i and k < j}.items(), key = lambda x: x[1], reverse = True)
            if bifurcations != []:
                temp_scores["bifurcation"] = bifurcations[0][1]

            temp_scores = sorted(temp_scores.items(), key = lambda item: item[1], reverse = True)

            s[i][j] = temp_scores[0][1]

    return s

"""
Implementation of the altered Wuchty algorithm for the multiple traceback 
of the filled out Nussinov DP table without the use of backpointers

Input:
    seq: RNA sequence
    P: filled out Nussinov DP table 
    limit: max number of optimal structures/solutions to return (default is 5 million)

Returns:
    list of optimal structures/solutions as list of pairs of nucleotides
"""
def multi_trace(seq, P, limit = 5000000):
    traces = []
    n = len(seq)

    sigma = [(0, n-1)]
    p = []

    L = (sigma, p)
    R = []

    R.append(L)

    index = 0

    while R != []:
        numfound = 0
        sigma, p = R.pop()

        if index == limit:
            break

        if sigma == []:
            traces.append(p)
            index += 1
            continue

        a, b = sigma.pop()

        if a < b and P[a][b] == P[a][b-1]:
            new_sigma = sigma.copy() + [(a, b-1)]
            L_prime = (new_sigma, p)
            R.append(L_prime)
            numfound += 1

        for l in range(a,b):
            if (seq[l], seq[b]) in pairs:
                temp_val = 1 + P[l+1][b-1]
                if a < l-1:
                    temp_val += P[a][l-1]
                if P[a][b] == temp_val:
                    new_sigma = sigma.copy() + [(l+1, b-1)]
                    new_sigma += [(a, l-1)]
                    new_p = p.copy() + [(l, b)]
                    L_prime = (new_sigma, new_p)
                    R.append(L_prime)
                    numfound += 1

        if numfound == 0:
            L_prime = (sigma, p)
            R.append(L_prime)

    return traces

"""
Converts a list of pairs of nucleotide (stem pairings) into a full dot-parenthesis notation

Input:
    seq: RNA sequence
    trace: list of pairs of nucleotides (stem pairings) of a RNA secondary structure

Returns:
    dot-parenthesis notation of the secondary structure (same length as RNA sequence)
"""
def dot_parenthesis_notation(seq, trace):
    converted = ["." for i in range(len(seq))]

    for left, right in trace:
        converted[left] = "("
        converted[right] = ")"

    return "".join(converted)

"""
Converts a RNA secondary structure from dot-parenthesis notation to a list of pairs of nucleotide (stem pairings)

Input:
    notation:  dot-parenthesis notation of the secondary structure

Returns:
    list of pairs of nucleotides (stem pairings) of a RNA secondary structure
"""
def notation_to_pairs(notation):
    unpaired = []
    pairs = []

    for i in range(len(notation)):
        if notation[i] == ".":
            continue
        elif notation[i] == "(":
            unpaired.append(i)
        elif notation[i] == ")":
            pairs.append((unpaired.pop(), i))

    return pairs

"""
Computes the accuracy of RNA secondary structure with a reference structure based on the stem pairs

Input:
    trace: list of pairs of nucleotides (stem pairings) of a predicted RNA secondary structure
    expected: reference dot-parenthesis notation of the secondary structure

Returns:
    accuracy score of the predicted structure against the reference structure
"""
def accuracy(trace, expected):
    actual = notation_to_pairs(expected)

    tp = 0
    tn = 0
    fp = 0
    fn = 0

    for pair in trace:
        x, y = pair

        if (x, y) in actual:
            tp += 1
        elif (y, x) in actual:
            tp += 1
        else:
            fp += 1

    fn = len(actual) - tp

    return (tp + tn) / (tp + tn + fp + fn)

"""
Parse RNA sequence and dot-bracket notation into graph components.

Input:
    sequence: RNA sequence
    dot_bracket: dot-parenthesis notation of the secondary structure

Returns:
    a NetworkX graph representation of the RNA
"""
def parse_dot_bracket(sequence, dot_bracket):
    stack = []
    graph = nx.Graph()

    for i, nucleotide in enumerate(sequence):
        graph.add_node(i, label=nucleotide)

    for i, symbol in enumerate(dot_bracket):
        if symbol == "(":
            stack.append(i)
        elif symbol == ")":
            if stack:
                j = stack.pop()
                graph.add_edge(i, j, type="pair")
        if i > 0:
            graph.add_edge(i, i - 1, type="covalent")

    return graph

"""
Visualize RNA secondary structure in 3D.

Input:
    graph: a NetworkX graph representation of the RNA
"""
def plot_rna_3d(graph):
    pos = nx.spring_layout(graph, dim=3)

    x_nodes = [pos[i][0] for i in graph.nodes()]
    y_nodes = [pos[i][1] for i in graph.nodes()]
    z_nodes = [pos[i][2] for i in graph.nodes()]

    edge_x = []
    edge_y = []
    edge_z = []
    for edge in graph.edges():
        x0, y0, z0 = pos[edge[0]]
        x1, y1, z1 = pos[edge[1]]
        edge_x.extend([x0, x1, None])
        edge_y.extend([y0, y1, None])
        edge_z.extend([z0, z1, None])

    color_map = {
        "A": "red",
        "U": "blue",
        "G": "yellow",
        "C": "green"
    }
    node_colors = [color_map[graph.nodes[i]["label"]] for i in graph.nodes()]

    edge_trace = go.Scatter3d(
        x=edge_x, y=edge_y, z=edge_z,
        mode='lines',
        line=dict(color='black', width=2),
        hoverinfo='none'
    )

    node_trace = go.Scatter3d(
        x=x_nodes, y=y_nodes, z=z_nodes,
        mode='markers+text',
        marker=dict(
            size=8,
            color=node_colors,
            line=dict(width=1, color='black')
        ),
        text=[graph.nodes[i]["label"] for i in graph.nodes()],
        textposition="top center",
        hoverinfo='text'
    )

    fig = go.Figure(data=[edge_trace, node_trace])
    fig.update_layout(
        title="RNA Secondary Structure in 3D",
        showlegend=False,
        scene=dict(
            xaxis=dict(showbackground=False),
            yaxis=dict(showbackground=False),
            zaxis=dict(showbackground=False),
        ),
    )
    fig.show()

"""
Wrapper function to combine Nussinov, multitrace, and visualization of the best predicted structure

Input:
    sequence: RNA sequence
    structure: reference dot-parenthesis notation of the secondary structure
    limit_structures: max number of optimal structures/solutions to return (default is 5 million)
    vis: boolean on if you want to visualize the best predicted structure or not
"""
def multi_nussinov(sequence, structure, limit_structures = 5000000, vis = True):
    s = nussinov(sequence)

    traces = multi_trace(sequence, s, limit_structures)

    print("Generated", len(traces), "optimal structures.")

    sum_accuracy = 0
    max_accuracy = 0

    best_struct = None

    for i in trange(len(traces)):
        trace = traces[i]

        score = accuracy(trace, structure)
        converted = dot_parenthesis_notation(sequence, trace)

        if score > max_accuracy:
            max_accuracy = score
            best_struct = converted

        sum_accuracy += score


    print("Best predicted structure accuracy:", max_accuracy)
    print("Average predicted structure accuracy:", sum_accuracy / len(traces))

    if vis:
        print("Best structure:", best_struct)
        rna_graph = parse_dot_bracket(sequence, best_struct)
        plot_rna_3d(rna_graph)

if __name__ == "__main__":
    ids = ["" for _ in range(8)]
    sequences = ["" for _ in range(8)]
    structures = ["" for _ in range(8)]

    with open("test_data.csv", "r") as data:
        csvreader = csv.reader(data)
        iter = 0

        for row in csvreader:
            ids[iter] = row[0]
            sequences[iter] = row[1]
            structures[iter] = row[2]

            iter += 1

    parser = argparse.ArgumentParser(description='RNA Structure Predictor and Visualizer')
    parser.add_argument('--ref',dest="reference", type=int, default=-1,
                        help='Test/reference sequences from RNA Central and RFam database (0-7)')
    parser.add_argument('--seq',dest="sequence", type=str, default="AU",
                        help='RNA sequence')
    parser.add_argument('--struct',dest="structure", type=str, default = "()",
                        help='RNA structure in dot-parenthesis notation')
    parser.add_argument('--limit',dest="limit", type=int, default = 5000000,
                        help="Limit on the number of optimal solutions to return (smaller numbers run faster and won't crash due to memory usage)")
    parser.add_argument('--vis',dest="visualize", type=bool, default = True,
                        help='Visualize the best predicted structure')
    args = parser.parse_args()

    if args.reference > -1:
        print("ID:", ids[args.reference])
        multi_nussinov(sequences[args.reference], structures[args.reference], args.limit, args.visualize)
    else:
        multi_nussinov(args.sequence, args.structure, args.limit, args.visualize)

