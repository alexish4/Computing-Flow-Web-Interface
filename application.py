from flask import Flask, request, jsonify, render_template
import numpy as np
import scipy as sp
from scipy.sparse import coo_matrix
import copy
import os
import networkx as nx
import pandas as pd
import heapq

app=Flask(__name__)

adj_matrix = []
source = 0
sink = 0

@app.route('/')
def index():
    return render_template('index.html')

def find_top_k_paths(G, source, sink, k=10):
    # Priority queue to store paths with their total betweenness value
    queue = [(0, source, [])]  # (total_betweenness, current_node, path)
    paths = []

    while queue and len(paths) < k:
        total_betw, current_node, path = heapq.heappop(queue)
        path = path + [current_node]

        if current_node == sink:
            paths.append((total_betw, path))
            continue

        for neighbor, edge_data in G[current_node].items():
            if neighbor not in path:  # Avoid cycles
                edge_betw = edge_data['betw']
                heapq.heappush(queue, (total_betw + edge_betw, neighbor, path))

    return paths

@app.route('/upload', methods=['POST'])
def upload_file():
    global adj_matrix
    global source
    global sink

    file = request.files['file']
    source = int(request.form['source']) - 1
    sink = int(request.form['sink']) - 1
    
    data = np.loadtxt(file, delimiter=',')
    rows, cols, weights = data[:, 0].astype(int), data[:, 1].astype(int), data[:, 2]
    adj_matrix = coo_matrix((weights, (rows, cols)))

     # Create graph
    G = nx.from_scipy_sparse_array(adj_matrix)

    # label nodes to match their index in the adjacency matrix
    mapping = {i: i for i in range(adj_matrix.shape[0])}
    G = nx.relabel_nodes(G, mapping)

    # Each edge will also have betweenness score
    largest_betweenness = 0
    for u, v, data in G.edges(data=True):
        # Calculate based on the indices of the source (u) and target (v)
        betw = get_betw_value(u, v)
        data['betw'] = betw 

        #also want largest betweenness value
        if largest_betweenness < betw:
            largest_betweenness = betw

    # Find the top 4 optimal paths from source to sink
    top_paths = find_top_k_paths(G, source, sink)
    top_paths = top_paths[::-1] #reversing order
    print (len(top_paths), " is length of top paths")

    # Convert top_paths to a format that is easy to send to the front-end
    top_paths_data = [
        {'total_betw': path[0], 'nodes': path[1]}
        for path in top_paths[:4]  # Limit to top 4 paths
    ]

    response_data = {
        'graph_data': nx.node_link_data(G),
        'top_paths': top_paths_data,
        'largest_betweenness': largest_betweenness
    }
    
    return jsonify(response_data)

def get_betw_value(u, v):
    global adj_matrix
    global source
    global sink

    tempAdjDense = adj_matrix.todense()
    
    # Convert adjacency matrix to Laplacian matrix
    tempLapDense = -copy.deepcopy(tempAdjDense)
    for ii, irow in enumerate(tempLapDense):
        tempLapDense[ii, ii] = -np.sum(irow)
    
    # Compute pseudoinverse of Laplacian matrix
    tempLinv = np.linalg.pinv(tempLapDense)
    
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
