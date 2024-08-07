from flask import Flask, request, jsonify, render_template
import numpy as np
import scipy as sp
from scipy.sparse import coo_matrix
import copy
import os

app=Flask(__name__)

@app.route('/')
def index():
    return render_template('index.html')

@app.route('/upload', methods=['POST'])
def upload_file():
    print("test")
    file = request.files['file']
    source = int(request.form['source'])
    sink = int(request.form['sink'])
    
    data = np.loadtxt(file, delimiter=',')
    rows, cols, weights = data[:, 0].astype(int), data[:, 1].astype(int), data[:, 2]
    adj_matrix = coo_matrix((weights, (rows, cols)))
    print(adj_matrix)
    
    betweenness_score = compute_flow_betweenness(adj_matrix, source, sink)
    
    return jsonify({'betweenness_score': betweenness_score})

def compute_flow_betweenness(adj_matrix, source, sink):
    tempAdjDense = adj_matrix.todense()
    print(tempAdjDense)
    
    # Convert adjacency matrix to Laplacian matrix
    tempLapDense = -copy.deepcopy(tempAdjDense)
    for ii, irow in enumerate(tempLapDense):
        tempLapDense[ii, ii] = -np.sum(irow)
    
    # Compute pseudoinverse of Laplacian matrix
    tempLinv = np.linalg.pinv(tempLapDense)
    
    # Compute flow betweenness for the given edge
    v_1_10_source = tempLinv[source, 0] - tempLinv[source, -1]
    v_1_10_sink = tempLinv[sink, 0] - tempLinv[sink, -1]
    b_source_sink = adj_matrix[source, sink] * (v_1_10_source - v_1_10_sink)
    
    return b_source_sink.item()  # Convert to a standard Python type

if __name__ == '__main__':
    app.run(debug=True)
    # port = int(os.environ.get("PORT", 5000))
    # app.run(host='0.0.0.0', port=port)
