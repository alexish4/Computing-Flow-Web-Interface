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
source = 0
sink = 0

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
    global source
    global sink

    file = request.files['file']
    source = int(request.form['source']) - 1
    sink = int(request.form['sink']) - 1
    
    if file.filename.endswith('.csv'):
        data = np.loadtxt(file, delimiter=',')
        rows, cols, correlations = data[:, 0].astype(int), data[:, 1].astype(int), data[:, 2]
    elif file.filename.endswith('.dat'):
        rows, cols, correlations = process_dat_file(file)
    else:
        return jsonify({'error': 'Unsupported file format'}), 400

    adj_matrix = coo_matrix((correlations, (rows, cols)))
    print(adj_matrix)

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
        betw = get_betw_value(u, v, tempLinv, tempAdjDense, source, sink)
        edge_length = -np.log(betw) #edge length is equal to -ln(|betw|) 

        data['betw'] = betw 
        data['edge_length'] = edge_length

        #also want largest betweenness value
        if largest_betweenness < betw:
            largest_betweenness = betw
    print(largest_betweenness, " is largest betweenness")

    # Find the top 4 optimal paths from source to sink
    top_paths = list(islice(nx.shortest_simple_paths(G, source, sink, "edge_length"), 4))
    
    #calculate path length
    path_lengths_edge_weights = []
    for path in top_paths:
        path_length = 0
        for i in range(len(path) - 1):
            path_length += G[path[i]][path[i + 1]]["edge_length"]
        path_lengths_edge_weights.append(path_length)
    print(top_paths)

    # Create the top_paths_data using path_lengths_edge_weights and top_paths_2
    top_paths_data = [
        {'edge_length': path_lengths_edge_weights[i], 'nodes': top_paths[i]}
        for i in range(len(top_paths))  # Limit to top 4 paths
    ]

    response_data = {
        'graph_data': nx.node_link_data(G),
        'top_paths': top_paths_data,
        'largest_betweenness': largest_betweenness
    }
    
    return jsonify(response_data)

def get_betw_value(u, v, tempLinv, tempAdjDense, source, sink):
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
    global source
    global sink

    resist_1 = int(request.form['resist1'])
    resist_2 = int(request.form['resist2'])
    
    tempAdjDense = adj_matrix.todense()
    
    # Convert adjacency matrix to Laplacian matrix
    tempLapDense = -copy.deepcopy(tempAdjDense)
    for ii, irow in enumerate(tempLapDense):
        tempLapDense[ii, ii] = -np.sum(irow)
    
    # Compute pseudoinverse of Laplacian matrix
    tempLinv = np.linalg.pinv(tempLapDense)
    
    # Compute flow betweenness for the given edge
    v_source_sink_resist1 = tempLinv[resist_1, source] - tempLinv[resist_1, sink]
    v_source_sink_resist2 = tempLinv[resist_2, source] - tempLinv[resist_2, sink]
    b_resist1_resist2 = tempAdjDense[resist_1, resist_2] * (v_source_sink_resist1 - v_source_sink_resist2)

    betweenness_score = b_resist1_resist2.item() # Convert to a standard Python type
    if betweenness_score < 0:
        betweenness_score *= -1
    
    return jsonify({'betweenness_score': betweenness_score})

if __name__ == '__main__':
    app.run(debug=True)
    # port = int(os.environ.get("PORT", 5000))
    # app.run(host='0.0.0.0', port=port)
