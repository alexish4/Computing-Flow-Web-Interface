from flask import Flask, request, jsonify, render_template
import numpy as np
import scipy as sp
from scipy.sparse import coo_matrix
import copy
import os
import networkx as nx
import pandas as pd

app=Flask(__name__)

adj_matrix = []
source = 0
sink = 0

@app.route('/')
def index():
    return render_template('index.html')

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
    mapping = {i: i + 1 for i in range(adj_matrix.shape[0])}
    G = nx.relabel_nodes(G, mapping)

    # Convert graph to D3.js compatible format
    data = nx.node_link_data(G)
    return jsonify(data)

@app.route('/calculate', methods=['POST'])
def compute_flow_betweenness():
    global adj_matrix
    global source
    global sink

    resist_1 = int(request.form['resist1']) - 1
    resist_2 = int(request.form['resist2']) - 1
    
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
    #app.run(debug=True)
    port = int(os.environ.get("PORT", 5000))
    app.run(host='0.0.0.0', port=port)
