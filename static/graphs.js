let resistor_1 = null;
let resistor_2 = null;

let source_id_array = [];
let sink_id_array = [];

let large_bet = 0;
let startingIndexValue = 1;

const svg = d3.select("#graph-svg");
const width = parseInt(svg.style("width"));
const height = parseInt(svg.style("height"));

// Group for zoomable content
const zoomGroup = svg.append("g");

// Zoom and pan functionality
const zoom = d3.zoom()
    .scaleExtent([0.1, 10])
    .on("zoom", event => {
        zoomGroup.attr("transform", event.transform);
    });

svg.call(zoom);

// Zoom control buttons
d3.select("#zoom-in-button").on("click", () => zoom.scaleBy(svg.transition().duration(500), 1.2));
d3.select("#zoom-out-button").on("click", () => zoom.scaleBy(svg.transition().duration(500), 0.8));
d3.select("#reset-zoom-button").on("click", () => svg.transition().duration(500).call(zoom.transform, d3.zoomIdentity));

// Function to display histogram images
function displayHistogram(imgData, imgData2) {
    const img1 = document.createElement('img');
    img1.src = 'data:image/png;base64,' + imgData;
    img1.alt = 'Histogram 1';
    
    const histogramContainer = document.getElementById('histogram-container');
    histogramContainer.innerHTML = ''; // Clear previous image
    histogramContainer.appendChild(img1); // Add new image

    const img2 = document.createElement('img');
    img2.src = 'data:image/png;base64,' + imgData2;
    img2.alt = 'Histogram 2';
    
    const histogramContainer2 = document.getElementById('histogram-container2');
    histogramContainer2.innerHTML = ''; // Clear previous image
    histogramContainer2.appendChild(img2); // Add new image
}

document.getElementById('upload-form').onsubmit = function(e) {
    e.preventDefault();
    
    const fileInput = document.getElementById('csv-file');
    const sourceResid = document.getElementById('source-resid').value;
    const sinkResid = document.getElementById('sink-resid').value;
    const k = document.getElementById('k-resid').value;

    //resetting resistor values
    resistor_1 = null
    resistor_2 = null

    // Get the selected radio button value
    startingIndexValue = document.querySelector('input[name="option"]:checked').value;
    if (startingIndexValue === "1") {
        console.log("Selected Starting Index:", startingIndexValue);
    }

    let sourceInput = document.getElementById('source-resid').value;
    let sinkInput = document.getElementById('sink-resid').value;

    // Process source input
    source_id_array = sourceInput
        .split(',')
        .map(node => node.trim())
        .filter(node => node !== '');

    if (source_id_array.length > 1) {
        document.getElementById("source-node-label").innerText = "Source Nodes: " + source_id_array.join(", ");
    } else{
        document.getElementById("source-node-label").innerText = "Source Node: " + source_id_array[0];
    }

    // Process sink input
    sink_id_array = sinkInput
        .split(',')
        .map(node => node.trim())
        .filter(node => node !== '');
        
    if (sink_id_array.length > 1) {
        document.getElementById("sink-node-label").innerText = "Sink Nodes: " + sink_id_array.join(", ");
    } else  {
        document.getElementById("sink-node-label").innerText = "Sink Node: " + sink_id_array[0];
    } 

    if (startingIndexValue === "1") {
        source_id_array = source_id_array.map(node => parseInt(node, 10) - 1); // Convert each to an integer
        sink_id_array = sink_id_array.map(node => parseInt(node, 10) - 1); // Convert each to an integer
    } else {
        source_id_array = source_id_array.map(node => parseInt(node, 10)); // Convert each to an integer
        sink_id_array = sink_id_array.map(node => parseInt(node, 10)); // Convert each to an integer
    }
    // Print arrays to console
    console.log("Source Nodes Array:", source_id_array);
    console.log("Sink Nodes Array:", sink_id_array);
    
    const formData = new FormData();
    formData.append('file', fileInput.files[0]);
    formData.append('source', source_id_array);
    formData.append('sink', sink_id_array);
    formData.append('k', k);
    fetch('/upload', {
        method: 'POST',
        body: formData
    })
    .then(response => response.json())
    .then(data => {
        if(data.incorrect_input === true) {
            alert("Input Out of Bounds")
        } else {
            large_bet = data.largest_betweenness;
            updateSourceNodeLabel();
            displayTopPaths(data.top_paths, data.top_paths2); // Display the top paths0                                
            drawGraph(data.graph_data);
            setupColorScaleAndEdges();
            drawColorScale();
            displayHistogram(data.histogram1, data.histogram2);
            drawCorrelationMatrix(data.graph_data);
        }
    }); 
};

function drawGraph(graph) {
    zoomGroup.selectAll("*").remove(); // Clear previous graph
    
    // Position the source and sink nodes
    const sourceNode = graph.nodes.find(d => d.id === source_id_array[0]);
    const sinkNode = graph.nodes.find(d => d.id === sink_id_array[0]);

    if (sourceNode) sourceNode.fx = 200; // Fixed x position for source node
    if (sinkNode) sinkNode.fx = width - 200; // Fixed x position for sink node

    const simulation = d3.forceSimulation(graph.nodes)
        .force("link", d3.forceLink(graph.links).id(d => d.id).distance(50).strength(1))
        .force("charge", d3.forceManyBody().strength(-1500))
        .force("center", d3.forceCenter(width / 2, height / 2));

    const link = zoomGroup.append("g")
        .attr("class", "links")
        .selectAll("line")
        .data(graph.links)
        .enter().append("line")
        .attr("stroke-width", 8)
        .attr("stroke", "#999");

    // Label edges with weights
    const linkLabels = zoomGroup.append("g")
        .attr("class", "link-labels")
        .selectAll("text")
        .data(graph.links)
        .enter().append("text")
        .attr("text-anchor", "middle")
        .attr("dy", -5)
        .attr("font-size", "12px")
        .attr("fill", "black")
        .text(d => `${d.weight.toFixed(3)} (${d.betw.toFixed(3)})`);

    const node = zoomGroup.append("g")
        .attr("class", "nodes")
        .selectAll("g")
        .data(graph.nodes)
        .enter().append("g");

    const circles = node.append("circle")
        .attr("r", 10)
        .attr("fill", function(d) {
            if (source_id_array.includes(d.id)) {
                return "green";
            } else if (sink_id_array.includes(d.id)) {
                return "red";
            } else {
                return "black";
            }
        })
        .on("click", function(event, d) {
            if (resistor_1 === d.id) {
                resistor_1 = null;
                d3.select(this).attr("fill", "black");
            } else if (resistor_2 === d.id) {
                resistor_2 = null;
                d3.select(this).attr("fill", "black");
            } else if (resistor_1 === null) {
                resistor_1 = d.id;
                d3.select(this).attr("fill", "yellow");
            } else if (resistor_2 === null) {
                resistor_2 = d.id;
                d3.select(this).attr("fill", "yellow");
            }

            updateCalculateButtonVisibility();

            console.log("Resistor 1 ID:", resistor_1, "Resistor 2 ID:", resistor_2);
        });

    const labels = node.append("text")
        .attr("x", 15)
        .attr("y", 5)
        .attr("fill", "black")
        .attr("font-size", "12px")
        .attr("text-anchor", "middle")
        .text(d => startingIndexValue === "1" ? d.id + 1 : d.id);

    simulation
        .nodes(graph.nodes)
        .on("tick", ticked);

    simulation.force("link")
        .links(graph.links);

    function ticked() {
        link
            .attr("x1", d => d.source.x)
            .attr("y1", d => d.source.y)
            .attr("x2", d => d.target.x)
            .attr("y2", d => d.target.y);

        linkLabels
            .attr("x", d => (d.source.x + d.target.x) / 2)
            .attr("y", d => (d.source.y + d.target.y) / 2);

        node
            .attr("transform", d => `translate(${d.x},${d.y})`);
    }

    function updateCalculateButtonVisibility() {
        const calculateButton = document.getElementById('calculate-button');
        if (resistor_1 !== null && resistor_2 !== null) {
            calculateButton.style.display = 'inline'; // Show button
        } else {
            calculateButton.style.display = 'none'; // Hide button
        }
    }

    // Handle Calculate button click
    document.getElementById('calculate-button').onclick = async function() {
        if (resistor_1 !== null && resistor_2 !== null) {
            const resist1Resid = resistor_1;
            const resist2Resid = resistor_2;

            const formData = new FormData();
            formData.append('resist1', resist1Resid);
            formData.append('resist2', resist2Resid);
            const response = await fetch('/calculate', {
                method: 'POST',
                body: formData
            });

            const result = await response.json();
            if (result.betweenness_score === 0) {
                alert("Please select 2 Resistors that have a direct connection");
            } else {
                document.getElementById('output').textContent = JSON.stringify(result, null, 2);
            }
        }
    };
    document.getElementById('refresh-button').onclick = async function() {
        resistor_1 = null;
        resistor_2 = null;
        drawGraph(graph);
        setupColorScaleAndEdges();
    };

}        

function displayTopPaths(paths, paths2) {
    const container = d3.select('#top-paths');
    container.selectAll("*").remove(); // Clear previous paths

    // Create a flex container to hold paths and paths2 side-by-side
    const flexContainer = container.append("div")
        .style("display", "flex"); // Use flexbox for horizontal alignment

    // Create a container for paths
    const pathsContainer = flexContainer.append("div")
        .html("<strong>Edge Length = -ln(betweenness)</strong>")
        .style("flex", "1");
        //.style("margin-right", "20px"); // Optional: Add space between paths and paths2


    paths.forEach((path, index) => {
        // Create a div to hold the path information and button
        const pathDiv = pathsContainer.append("div")
            .style("display", "flex") // Use flexbox for horizontal alignment
            .style("align-items", "center") // Center items vertically
            .style("margin-bottom", "10px"); // Optional: Add space between rows

        const formattedNodes = path.nodes.map(node => startingIndexValue === "1" ? node + 1 : node).join(" -> ");
        
        pathDiv.append("div")
            .text(`Path ${index + 1}: Total Path Length From Betweenness: ${path.edge_length}, Path: ${formattedNodes}`)
            .style("margin-right", "10px"); // Optional: Add space between text and button

        // Add a button to the div
        pathDiv.append("button")
            .text("Highlight") // Change this text to whatever you want the button to display
            .on("click", () => {
                // Define what should happen when the button is clicked
                //alert(`Button for Path ${index + 1} clicked!`);
                highlightPathEdges(path);
            });
    });

    const paths2Container = flexContainer.append("div")
        .html("<strong>Edge Length = -ln(correlation)</strong>")
        .style("flex", "1");

    paths2.forEach((path, index) => {
        // Create a div to hold the path information and button
        const pathDiv = paths2Container.append("div")
            .style("display", "flex") // Use flexbox for horizontal alignment
            .style("align-items", "center") // Center items vertically
            .style("margin-bottom", "10px"); // Optional: Add space between rows

        const formattedNodes = path.nodes.map(node => startingIndexValue === "1" ? node + 1 : node).join(" -> ");
        
        pathDiv.append("div")
            .text(`Path ${index + 1}: Total Path Length From Correlation: ${path.edge_length}, Path: ${formattedNodes}`)
            .style("margin-right", "10px"); // Optional: Add space between text and button

        // Add a button to the div
        pathDiv.append("button")
            .text("Highlight") // Change this text to whatever you want the button to display
            .on("click", () => {
                // Define what should happen when the button is clicked
                //alert(`Button for Path ${index + 1} clicked!`);
                highlightPathEdges(path);
            });
    });
}

function highlightPathEdges(path) {
    // Convert the path's nodes into a set of edges
    const edges = [];
    for (let i = 0; i < path.nodes.length - 1; i++) {
        edges.push({
            source: path.nodes[i],
            target: path.nodes[i + 1]
        });
    }

    // Update the styles of the links to highlight the edges in the path
    d3.selectAll(".links line")
        .attr("stroke", d => {
            // Check if the current link is part of the highlighted path
            const isHighlighted = edges.some(edge =>
                (d.source.id === edge.source && d.target.id === edge.target) ||
                (d.source.id === edge.target && d.target.id === edge.source)
            );
            return isHighlighted ? "orange" : "#999";
        })
        .attr("stroke-width", d => {
            // Increase the stroke width for highlighted edges
            const isHighlighted = edges.some(edge =>
                (d.source.id === edge.source && d.target.id === edge.target) ||
                (d.source.id === edge.target && d.target.id === edge.source)
            );
            return isHighlighted ? 16 : 8;
        });
}

// Function to update the source node label
function updateSourceNodeLabel() {
    const sourceNodeLabel = document.getElementById('source-node-label');
    const sinkNodeLabel = document.getElementById('sink-node-label');
    const dirLabel = document.getElementById('directions-label');
    const refreshButton = document.getElementById('refresh-button');
    
    // Show the label
    sourceNodeLabel.style.display = 'inline';
    sinkNodeLabel.style.display = 'inline'; 
    dirLabel.style.display = 'inline';
    refreshButton.style.display = 'inline';

}

function setupColorScaleAndEdges() {
    // Define the color scale
    const colorScale = d3.scaleLinear()
        .domain([0, large_bet]) // Using the maximum edge frequency
        .range(["lightgray", "red"]); // Low frequency -> black, High frequency -> red

    // Implement the colorEdges Function
    function colorEdges() {
        console.log("Color scale domain:", colorScale.domain());
        console.log("Color scale range:", colorScale.range());

        d3.selectAll(".links line") // Select all line elements in the links group
            .attr("stroke-width", 8)
            .attr("stroke", d => colorScale(d.betw)); // Set stroke color based on betweenness
    }
    colorEdges();
}

function drawColorScale() {
    const svg = d3.select("svg");
    
    // Remove any existing color scales
    svg.selectAll(".color-scale").remove();
    
    // Define dimensions and position for the color scale
    const width = 20;
    const height = 300;
    const x = svg.attr("width") - 60; // Positioning on the right side
    const y = 50;
    
    // Create a group for the color scale
    const colorScaleGroup = svg.append("g")
        .attr("class", "color-scale")
        .attr("transform", `translate(${x}, ${y})`);
    
    // Define the gradient
    const gradient = colorScaleGroup.append("defs")
        .append("linearGradient")
        .attr("id", "color-gradient")
        .attr("x1", "0%")
        .attr("y1", "100%")
        .attr("x2", "0%")
        .attr("y2", "0%");
    
    gradient.append("stop")
        .attr("offset", "0%")
        .attr("stop-color", "lightgray"); // Lightest
    
    gradient.append("stop")
        .attr("offset", "100%")
        .attr("stop-color", d3.interpolateReds(1)); // Darkest
    
    // Draw the rectangle filled with the gradient
    colorScaleGroup.append("rect")
        .attr("width", width)
        .attr("height", height)
        .style("fill", "url(#color-gradient)");
    
    // Define the scale for the color scale's axis
    const scale = d3.scaleLinear()
        .domain([0, large_bet])
        .range([height, 0]);
    
    // Define the axis for the color scale
    const axis = d3.axisRight(scale)
        .ticks(6); // Adjust the number of ticks as needed
    
    // Draw the axis
    colorScaleGroup.append("g")
        .attr("class", "axis")
        .attr("transform", `translate(${width}, 0)`)
        .call(axis);
}

function drawCorrelationMatrix(graph) {
    const correlationSvg = d3.select("#correlation-svg")
                             .attr("width", 500)
                             .attr("height", 500);
    const gridSize = 20;
    const nodes = graph.nodes;  
    const links = graph.links;

    const colorScale = d3.scaleSequential()
        .domain([0, 1])
        .interpolator(d3.interpolateRdBu);

    const cells = correlationSvg.selectAll("rect")
        .data(links)
        .enter().append("rect")
        .attr("x", d => nodes.findIndex(n => n.id === d.source.id) * gridSize)
        .attr("y", d => nodes.findIndex(n => n.id === d.target.id) * gridSize)
        .attr("width", gridSize)
        .attr("height", gridSize)
        .attr("fill", d => colorScale(d.weight))
        .on("mouseover", function(event, d) {
            d3.select(this).attr("stroke", "black").attr("stroke-width", 2);

            d3.select("#tooltip")
                .style("left", (event.pageX + 10) + "px")
                .style("top", (event.pageY - 10) + "px")
                .style("opacity", 1)
                .html(`Nodes: ${nodes.find(n => n.id === d.source.id).name} & ${nodes.find(n => n.id === d.target.id).name}<br>Correlation: ${d.weight.toFixed(3)}`);
        })
        .on("mouseout", function(d) {
            d3.select(this).attr("stroke", "none");

            d3.select("#tooltip")
                .style("opacity", 0);
        });

    correlationSvg.selectAll(".rowLabel")
        .data(nodes)
        .enter().append("text")
        .attr("x", 0)
        .attr("y", (d, i) => i * gridSize + gridSize / 2)
        .attr("dy", ".35em")
        .attr("text-anchor", "end")
        .text(d => d.name);

    correlationSvg.selectAll(".colLabel")
        .data(nodes)
        .enter().append("text")
        .attr("x", (d, i) => i * gridSize + gridSize / 2)
        .attr("y", 0)
        .attr("dy", ".35em")
        .attr("text-anchor", "middle")
        .text(d => d.name)
        .attr("transform", "rotate(-90)");

    d3.select("body").append("div")
        .attr("id", "tooltip")
        .attr("class", "tooltip")
        .style("opacity", 0);
}


// Function to open a tab
function openTab(evt, tabName) {
    let i, tabcontent, tablinks;
    tabcontent = document.getElementsByClassName("tabcontent");
    for (i = 0; i < tabcontent.length; i++) {
        tabcontent[i].style.display = "none";
    }
    tablinks = document.getElementsByClassName("tablinks");
    for (i = 0; i < tablinks.length; i++) {
        tablinks[i].className = tablinks[i].className.replace(" active", "");
    }
    document.getElementById(tabName).style.display = "block";
    evt.currentTarget.className += " active";
}

// Open the first tab by default
document.querySelector(".tablinks").click();