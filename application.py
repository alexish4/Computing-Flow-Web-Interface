from flask import Flask, request, jsonify, render_template
import numpy as np
import scipy as sp
from scipy.sparse import coo_matrix
import copy
import os
import networkx as nx
import pandas as pd
import heapq
from itertools import islice

app=Flask(__name__)

adj_matrix = []
source_array = []
sink_array = []
betweenness_values = []
#387,388,389,389,390,391,392
#328,329,334,338,378,348

@app.route('/betweenness-histogram', methods=['POST'])
def betweenness_histogram():
    global betweenness_values

    # Prepare histogram data
    histogram, bin_edges = np.histogram(betweenness_values, bins='auto')

    # Create histogram data structure for JSON response
    histogram_data = {
        'histogram': histogram.tolist(),
        'bin_edges': bin_edges.tolist()
    }

    return jsonify(histogram_data)

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
    
    k = 10 #by default k is 10

    if request.form['k'] != "":
        k = int(request.form['k'])
    print(k, " is k")
    
    if file.filename.endswith('.csv'):
        data = np.loadtxt(file, delimiter=',')
        rows, cols, correlations = data[:, 0].astype(int), data[:, 1].astype(int), data[:, 2]
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

    # Each edge will also have betweenness score
    largest_betweenness = 0

    tempAdjDense = adj_matrix.todense()
    
    # Convert adjacency matrix to Laplacian matrix
    tempLapDense = -copy.deepcopy(tempAdjDense)
    for ii, irow in enumerate(tempLapDense):
        tempLapDense[ii, ii] = -np.sum(irow)
    
    # Compute pseudoinverse of Laplacian matrix
    tempLinv = np.linalg.pinv(tempLapDense)

    for u, v, data in G.edges(data=True):
        # Calculate based on the indices of the source (u) and target (v)
        betw = get_betw_value(u, v, tempLinv, tempAdjDense)
        edge_length = -np.log(betw) #edge length is equal to -ln(|betw|) 
        edge_length2 = -np.log(data['weight'])
        print(edge_length2)

        data['betw'] = betw 
        data['edge_length'] = edge_length
        data['edge_length2'] = edge_length2

        betweenness_values.append(betw)

        #also want largest betweenness value
        if largest_betweenness < betw:
            largest_betweenness = betw

    # Find the top k optimal paths from source to sink
    top_paths = list(islice(nx.shortest_simple_paths(G, source_array[0], sink_array[0], "edge_length"), k))
    top_paths2 = list(islice(nx.shortest_simple_paths(G, source_array[0], sink_array[0], "edge_length2"), k))
    print(top_paths2, " is top paths 2")
    
    #calculate path length
    path_lengths_edge_weights = []
    path_lengths_edge_weights2 = []
    for path in top_paths:
        path_length = 0
        path_length2 = 0
        for i in range(len(path) - 1):
            path_length += G[path[i]][path[i + 1]]["edge_length"]
            path_length2 += G[path[i]][path[i + 1]]["edge_length2"]
        path_lengths_edge_weights.append(path_length)
        path_lengths_edge_weights2.append(path_length2)
    print(top_paths, " is top paths")
    print(path_lengths_edge_weights, " is from betweenness")
    print(path_lengths_edge_weights2, " is from correlation")

    # Create the top_paths_data using path_lengths_edge_weights and top_paths_2
    top_paths_data = [
        {'edge_length': path_lengths_edge_weights[i], 'nodes': top_paths[i]}
        for i in range(len(top_paths))  # Limit to top 4 paths
    ]
    top_paths_data2 = [
        {'edge_length': path_lengths_edge_weights2[i], 'nodes': top_paths2[i]}
        for i in range(len(top_paths2))  # Limit to top 4 paths
    ]

    response_data = {
        'graph_data': nx.node_link_data(G),
        'top_paths': top_paths_data,
        'top_paths2': top_paths_data2,
        'largest_betweenness': largest_betweenness
    }
    
    return jsonify(response_data)

def get_betw_value(u, v, tempLinv, tempAdjDense):
    global source_array
    global sink_array

    total_betweenness_score = 0

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

if __name__ == '__main__':
    app.run(debug=True)
    # port = int(os.environ.get("PORT", 5000))
    # app.run(host='0.0.0.0', port=port)
