from flask import Flask, request, jsonify, render_template
import numpy as np
import scipy as sp
from scipy.sparse import coo_matrix
import copy
import os
import io
import base64
import networkx as nx
import pandas as pd
import heapq
from itertools import islice
import matplotlib
#matplotlib.use('Agg')  # Use a non-interactive backend
import matplotlib.pyplot as plt
from collections import Counter

app=Flask(__name__)
application=app

adj_matrix = []
source_array = []
sink_array = []
#387,388,389,389,390,391,392
#328,329,334,338,378,348

@app.route('/')
def index():
    return render_template('index.html')

def process_dat_file(file):
    data = np.loadtxt(file)
    num_nodes = data.shape[0]
    rows = []
    cols = []
    correlations = []

    for i in range(num_nodes):
        for j in range(num_nodes):
            if i != j:  # Exclude self-loops
                mutual_info = data[i, j]
                if mutual_info > 0 and not np.isnan(mutual_info) and not np.isinf(mutual_info):  # Filter invalid weights:
                    # Add both (i, j) and (j, i) to ensure bidirectional edges
                    rows.extend([i, j])
                    cols.extend([j, i])
                    correlations.extend([mutual_info, mutual_info])

    return rows, cols, correlations


@app.route('/upload', methods=['POST'])
def upload_file():
    global adj_matrix
    global source_array
    global sink_array

    file = request.files['file']
    source_array = list(map(int, request.form['source'].split(',')))
    sink_array = list(map(int, request.form['sink'].split(',')))

    all = False #calculate average or all
    
    k = 51 #by default k is 10

    largest_betweenness = 0

    if request.form['k'] != "": #if k has input
        k = int(request.form['k'])
    if request.form['average'] == "No":
        all = True
    
    if file.filename.endswith('.csv'):
        data = np.loadtxt(file, delimiter=',')
        filtered_data = data[data[:, 2] != 0]
        rows, cols, correlations = filtered_data[:, 0].astype(int), filtered_data[:, 1].astype(int), filtered_data[:, 2]
    elif file.filename.endswith('.dat'):
        rows, cols, correlations = process_dat_file(file)
    else:
        return jsonify({'error': 'Unsupported file format'}), 400

    adj_matrix = coo_matrix((correlations, (rows, cols)))

     # Create graph
    G = nx.from_scipy_sparse_array(adj_matrix)

    # label nodes to match their index in the adjacency matrix
    mapping = {i: i for i in range(adj_matrix.shape[0])}
    G = nx.relabel_nodes(G, mapping)

    tempAdjDense = adj_matrix.todense()
    
    # Convert adjacency matrix to Laplacian matrix
    tempLapDense = -copy.deepcopy(tempAdjDense)
    for ii, irow in enumerate(tempLapDense):
        tempLapDense[ii, ii] = -np.sum(irow)
    
    # Compute pseudoinverse of Laplacian matrix
    tempLinv = np.linalg.pinv(tempLapDense)

    #variable to check if input is out of bounds
    incorrect_input = False

    #creating multiple graphs so you aren't just using average betweenness
    array_of_graphs = []
    if all:
        for so in source_array:
            for si in source_array:
                array_of_graphs.append(G.copy())
    print(len(array_of_graphs))

    #need the average graph because this is the graph we are drawing
    for u, v, data in G.edges(data=True):
        # Calculate based on the indices of the source (u) and target (v)
        betw = get_betw_value(u, v, tempLinv, tempAdjDense)
        if betw is None:
            incorrect_input = True
            response_data = {
                'incorrect_input': incorrect_input
            }
            return jsonify(response_data)
        edge_length = -np.log(betw) #edge length is equal to -ln(|betw|) 
        edge_length2 = -np.log(data['weight'])

        data['betw'] = betw 
        data['edge_length'] = edge_length
        data['edge_length2'] = edge_length2

        #also want largest betweenness value
        if largest_betweenness < betw:
            largest_betweenness = betw

    if not all:
        top_paths, top_paths2, top_paths_lengths, top_paths2_lengths, most_important_nodes, most_important_nodes2 = generateTopPaths(G, k, tempLinv, tempAdjDense)

    else:
        top_paths, top_paths2, top_paths_lengths, top_paths2_lengths, most_important_nodes, most_important_nodes2 = generateTopPaths2(array_of_graphs, k, tempLinv, tempAdjDense)

    #save histograms to different page
    img_data, img_data2 = histograms(top_paths_lengths, top_paths2_lengths)

    frequencyGraph(most_important_nodes, most_important_nodes2)

    # Create the top_paths_data using path_lengths_edge_weights and top_paths_2
    top_paths_data = [
        {'edge_length': top_paths_lengths[i], 'nodes': top_paths[i]}
        for i in range(len(top_paths))  
    ]
    top_paths_data2 = [
        {'edge_length': top_paths2_lengths[i], 'nodes': top_paths2[i]}
        for i in range(len(top_paths2))  
    ]

    ranked_nodes_data = [
        {'node': node, 'frequency': freq}
        for node, freq in most_important_nodes
    ]

    ranked_nodes_data2 = [
        {'node': node, 'frequency': freq}
        for node, freq in most_important_nodes2
    ]


    response_data = {
        'graph_data': nx.node_link_data(G),
        'top_paths': top_paths_data,
        'top_paths2': top_paths_data2,
        'histogram1': img_data,
        'histogram2': img_data2,
        'ranked_nodes_data': ranked_nodes_data,
        'ranked_nodes_data2': ranked_nodes_data2,
        'incorrect_input': incorrect_input,
        'largest_betweenness': largest_betweenness
    }
    
    return jsonify(response_data)

def get_betw_value(u, v, tempLinv, tempAdjDense):
    global source_array
    global sink_array

    total_betweenness_score = 0

    try:
        for s in source_array:
            for t in sink_array:
                # Compute flow betweenness for the given edge
                v_source_sink_resist1 = tempLinv[u, s] - tempLinv[u, t]
                v_source_sink_resist2 = tempLinv[v, s] - tempLinv[v, t]
                b_resist1_resist2 = tempAdjDense[u, v] * (v_source_sink_resist1 - v_source_sink_resist2)
                total_betweenness_score += b_resist1_resist2

        #divide by number of combinations
        num_of_combinations = len(source_array) * len(sink_array)
        total_betweenness_score /= num_of_combinations

        betweenness_score = total_betweenness_score.item() # Convert to a standard Python type
        if betweenness_score < 0:
            betweenness_score *= -1

        return betweenness_score
    except IndexError as e:
        return None

def get_betw_value2(u, v, tempLinv, tempAdjDense, source, sink):
    # Compute flow betweenness for the given edge
    v_source_sink_resist1 = tempLinv[u, source] - tempLinv[u, sink]
    v_source_sink_resist2 = tempLinv[v, source] - tempLinv[v, sink]
    b_resist1_resist2 = tempAdjDense[u, v] * (v_source_sink_resist1 - v_source_sink_resist2)

    betweenness_score = b_resist1_resist2.item() # Convert to a standard Python type
    if betweenness_score < 0:
        betweenness_score *= -1

    return betweenness_score

@app.route('/calculate', methods=['POST'])
def compute_flow_betweenness():
    global adj_matrix
    global source_array
    global sink_array

    resist_1 = int(request.form['resist1'])
    resist_2 = int(request.form['resist2'])
    
    tempAdjDense = adj_matrix.todense()
    
    # Convert adjacency matrix to Laplacian matrix
    tempLapDense = -copy.deepcopy(tempAdjDense)
    for ii, irow in enumerate(tempLapDense):
        tempLapDense[ii, ii] = -np.sum(irow)
    
    # Compute pseudoinverse of Laplacian matrix
    tempLinv = np.linalg.pinv(tempLapDense)
    
    total_betweenness_score = 0

    for s in source_array:
        for t in sink_array:
            # Compute flow betweenness for the given edge
            v_source_sink_resist1 = tempLinv[resist_1, s] - tempLinv[resist_1, t]
            v_source_sink_resist2 = tempLinv[resist_2, s] - tempLinv[resist_2, t]
            b_resist1_resist2 = tempAdjDense[resist_1, resist_2] * (v_source_sink_resist1 - v_source_sink_resist2)
            total_betweenness_score += b_resist1_resist2

    #divide by number of combinations
    num_of_combinations = len(source_array) * len(sink_array)
    total_betweenness_score /= num_of_combinations

    betweenness_score = total_betweenness_score.item() # Convert to a standard Python type
    if betweenness_score < 0:
        betweenness_score *= -1
    
    return jsonify({'betweenness_score': betweenness_score})

def histograms(path_lengths, path_lengths2):
    weights = np.ones_like(path_lengths) / len(path_lengths)
    bins = int(len(path_lengths) / 2 )
    if len(path_lengths) > 51:
        bins = int(np.sqrt(len(path_lengths)))

    # Generate histogram
    plt.figure()
    plt.hist(path_lengths, bins=bins, weights=weights, alpha=0.5, label='Path Length(bet)')
    plt.xlabel('Path Length')
    plt.ylabel('Probability')
    plt.legend(loc='upper right')
    plt.title('Path Length(flow-betweenness)')

    # Save to a bytes buffer
    buf = io.BytesIO()
    plt.savefig(buf, format='png')
    buf.seek(0)
    img_data = base64.b64encode(buf.getvalue()).decode('utf8')
    buf.close()
    plt.close()

    weights = np.ones_like(path_lengths2) / len(path_lengths2)
    # Generate histogram
    plt.figure()
    plt.hist(path_lengths2, bins=bins, weights=weights, alpha=0.5, label='Path Length(cor)')
    plt.xlabel('Path Length')
    plt.ylabel('Probability')
    plt.legend(loc='upper right')
    plt.title('Path Length(correlation)')

    # Save to a bytes buffer
    buf = io.BytesIO()
    plt.savefig(buf, format='png')
    buf.seek(0)
    img_data2 = base64.b64encode(buf.getvalue()).decode('utf8')
    buf.close()
    plt.close()

    return img_data, img_data2

def frequencyGraph(most_important_nodes, most_important_nodes2):
    # Extract nodes and frequencies
    nodes1 = [data[0] for data in most_important_nodes]
    frequencies1 = [data[1] for data in most_important_nodes]

    nodes2 = [data[0] for data in most_important_nodes2]
    frequencies2 = [data[1] for data in most_important_nodes2]

    #normalize between 0 and 1 for frequencies
    min_val = min(frequencies1)
    max_val = max(frequencies1)
    frequencies1_normalized = [(v - min_val) / (max_val - min_val) for v in frequencies1]
    min_val = min(frequencies2)
    max_val = max(frequencies2)
    frequencies2_normalized = [(v - min_val) / (max_val - min_val) for v in frequencies2]

    # Handle the mismatch in lengths by padding the shorter list with zeros
    max_length = max(len(nodes1), len(nodes2))
    
    nodes1.extend([None] * (max_length - len(nodes1)))
    frequencies1_normalized.extend([0] * (max_length - len(frequencies1_normalized)))

    nodes2.extend([None] * (max_length - len(nodes2)))
    frequencies2_normalized.extend([0] * (max_length - len(frequencies2_normalized)))

    # Plotting
    plt.figure(figsize=(14, 7))

    # Plot for most_important_nodes
    bar_width = 0.35
    index = range(max_length)

    plt.bar(index, frequencies1_normalized, width=bar_width, color='blue', label='Ranked Nodes 1')
    plt.bar([i + bar_width for i in index], frequencies2_normalized, width=bar_width, color='red', label='Ranked Nodes 2')

    # Adding labels and title
    plt.xlabel('Nodes')
    plt.ylabel('Normalized Frequency')
    plt.title('Frequency of Most Important Nodes')
    plt.xticks([i + bar_width / 2 for i in index], [str(n) for n in nodes1], rotation=90)
    plt.legend()
    plt.grid(True)

    # Show plot
    # plt.tight_layout()
    # plt.show()




def generateTopPaths(G, k, tempLinv, tempAdjDense):
    global largest_betweenness
    
    top_paths = []
    top_paths_lengths = []
    top_paths2 = []
    top_paths2_lengths = []

    for so in source_array:
        for si in sink_array:
            # Find the top k optimal paths from source to sink
            paths = list(islice(nx.shortest_simple_paths(G, so, si, weight="edge_length"), k))
            paths2 = list(islice(nx.shortest_simple_paths(G, so, si, weight="edge_length2"), k))
            
            # Calculate path lengths for the first set of paths
            lengths = [sum(G[u][v]["edge_length"] for u, v in zip(path[:-1], path[1:])) for path in paths]
            lengths2 = [sum(G[u][v]["edge_length2"] for u, v in zip(path[:-1], path[1:])) for path in paths2]

            # Store the top paths and their lengths
            top_paths.extend(paths)
            top_paths_lengths.extend(lengths)
            top_paths2.extend(paths2)
            top_paths2_lengths.extend(lengths2)

    # Remove duplicates by converting paths to a tuple and using a set
    unique_paths_with_lengths = list({tuple(path): length for length, path in zip(top_paths_lengths, top_paths)}.items())
    unique_paths2_with_lengths = list({tuple(path): length for length, path in zip(top_paths2_lengths, top_paths2)}.items())
 
    # Sort the unique paths and lengths by length
    sorted_paths_with_lengths = sorted(unique_paths_with_lengths, key=lambda x: x[1])[:k]
    sorted_paths2_with_lengths = sorted(unique_paths2_with_lengths, key=lambda x: x[1])[:k]

    # Unpack the sorted pairs back into the arrays
    top_paths, top_paths_lengths = zip(*sorted_paths_with_lengths) if sorted_paths_with_lengths else ([], [])
    top_paths2, top_paths2_lengths = zip(*sorted_paths2_with_lengths) if sorted_paths2_with_lengths else ([], [])

    # Convert the results back to lists if needed
    top_paths = list(top_paths)
    top_paths_lengths = list(top_paths_lengths)
    top_paths2 = list(top_paths2)
    top_paths2_lengths = list(top_paths2_lengths)

    # Create an array of the most important nodes based on frequency
    all_nodes_in_paths = [node for path in top_paths for node in path]  # Flatten the list of top paths
    all_nodes_in_paths2 = [node for path in top_paths2 for node in path]  # Include nodes from the second set of paths
    node_frequencies = Counter(all_nodes_in_paths)  # Count the frequency of each node
    node_frequencies2 = Counter(all_nodes_in_paths2)

    # Create a list of (node, frequency_value) tuples sorted by frequency
    most_important_nodes = [(node, freq) for node, freq in node_frequencies.most_common()]
    most_important_nodes2 = [(node, freq) for node, freq in node_frequencies2.most_common()]

    print(most_important_nodes)
    print(most_important_nodes2)

    return top_paths, top_paths2, top_paths_lengths, top_paths2_lengths, most_important_nodes, most_important_nodes2

def generateTopPaths2(array_of_graphs, k, tempLinv, tempAdjDense):
    array_index = 0
    for so in source_array:
        for si in sink_array:
            for u, v, data in array_of_graphs[array_index].edges(data=True):
                # Calculate based on the indices of the source (u) and target (v)
                epsilon = 1e-10 #prevent betw being so small that recorded as 0 which is bad for ln
                betw = get_betw_value2(u, v, tempLinv, tempAdjDense, so, si)
                edge_length = -np.log(max(betw, epsilon)) #edge length is equal to -ln(|betw|) 
                edge_length2 = -np.log(data['weight'])

                data['betw'] = betw 
                data['edge_length'] = edge_length
                data['edge_length2'] = edge_length2
            array_index += 1

    top_paths = []
    top_paths2 = []
    top_paths_lengths = []
    top_paths2_lengths = []

    array_index = 0
    for so in source_array:
        for si in sink_array:
            top_path, top_path2, length, length2 = miniGenerateTopPaths(array_of_graphs[array_index], k, so, si)
            top_paths.extend(top_path)
            top_paths_lengths.extend(length)
            top_paths2.extend(top_path2)
            top_paths2_lengths.extend(length2)
            array_index += 1
    #For betweenness
    unique_paths_with_lengths = list({tuple(path): length for length, path in zip(top_paths_lengths, top_paths)}.items())
    sorted_paths_with_lengths = sorted(unique_paths_with_lengths, key=lambda x: x[1])[:k]
    top_paths, top_paths_lengths = zip(*sorted_paths_with_lengths) if sorted_paths_with_lengths else ([], [])

    #for correlation
    unique_paths_with_lengths2 = list({tuple(path): length for length, path in zip(top_paths2_lengths, top_paths2)}.items())
    sorted_paths_with_lengths2 = sorted(unique_paths_with_lengths2, key=lambda x: x[1])[:k]
    top_paths2, top_paths2_lengths = zip(*sorted_paths_with_lengths2) if sorted_paths_with_lengths2 else ([], [])

    top_paths = list(top_paths)
    top_paths_lengths = list(top_paths_lengths)
    top_paths2 = list(top_paths2)
    top_paths2_lengths = list(top_paths2_lengths)

    # Create an array of the most important nodes based on frequency
    all_nodes_in_paths = [node for path in top_paths for node in path]  # Flatten the list of top paths
    all_nodes_in_paths2 = [node for path in top_paths2 for node in path]  # Include nodes from the second set of paths
    node_frequencies = Counter(all_nodes_in_paths)  # Count the frequency of each node
    node_frequencies2 = Counter(all_nodes_in_paths2)

    # Create a list of (node, frequency_value) tuples sorted by frequency
    most_important_nodes = [(node, freq) for node, freq in node_frequencies.most_common()]
    most_important_nodes2 = [(node, freq) for node, freq in node_frequencies2.most_common()]

    print(most_important_nodes)
    print(most_important_nodes2)

    return top_paths, top_paths2, top_paths_lengths, top_paths2_lengths, most_important_nodes, most_important_nodes2

def miniGenerateTopPaths(graph, k, so, si):
    top_paths = list(islice(nx.shortest_simple_paths(graph, so, si, weight="edge_length"), k))
    top_paths2 = list(islice(nx.shortest_simple_paths(graph, so, si, weight="edge_length2"), k))
    top_paths_lengths = [sum(graph[u][v]["edge_length"] for u, v in zip(path[:-1], path[1:])) for path in top_paths]
    top_paths2_lengths = [sum(graph[u][v]["edge_length2"] for u, v in zip(path[:-1], path[1:])) for path in top_paths2]

    return top_paths, top_paths2, top_paths_lengths, top_paths2_lengths


if __name__ == '__main__':
    app.run(debug=True)
    # port = int(os.environ.get("PORT", 5000))
    # app.run(host='0.0.0.0', port=port)
