use std::fs;
use std::io;
use rand::distributions::Uniform;
use rand::distributions::Distribution;
use std::fs::File;
use std::io::{Write};
use std::collections::HashMap;
use plotters::prelude::*;
use std::error::Error;
use csv::Writer;
use petgraph::graph::{UnGraph};
use petgraph::dot::{Dot, Config};



//Generating K-mers
//writing a function to generate k-mers from a given DNA sequence
fn generate_kmers(dna_sequence: &str, k: usize) -> Vec<String> {
    let mut kmers = Vec::new();
    for i in 0..=dna_sequence.len() - k {
        kmers.push(dna_sequence[i..i + k].to_string());
    }
    kmers
}


//Counting K-mers Using Hashing
//use a HashMap to count the frequency of each k-mer

fn count_kmers(kmers: Vec<String>) -> HashMap<String, usize> {
    let mut kmer_counts = HashMap::new();
    for kmer in kmers {
        *kmer_counts.entry(kmer).or_insert(0) += 1;
    }
    kmer_counts
}

// Define the structure for a De Bruijn graph.
struct DeBruijnGraph {
    // Use a HashMap to represent the graph.
    // Keys are String representing a k-1-mer (a node),
    // Values are Vec<String> representing a list of adjacent k-1-mers (edges).
    edges: HashMap<String, Vec<String>>,
}

// Implement methods for the DeBruijnGraph structure.
impl DeBruijnGraph {
    // Constructor method to create a new De Bruijn graph from a list of k-mers.
    fn new(kmers: &[String]) -> Self {
        let mut edges = HashMap::new();

        // Iterate over each k-mer in the given list.
        for kmer in kmers {
            // Split the k-mer into two parts:
            // - The first part (node) is all but the last character.
            // - The second part (next) is all but the first character.
            // This split creates an overlap between the k-1-mers.
            let (node, next) = kmer.split_at(kmer.len() - 1);

            // Insert the node into the HashMap if it doesn't exist,
            // and append the adjacent node to the list of edges.
            edges.entry(node.to_string()).or_insert_with(Vec::new).push(next.to_string());
        }

        // Return a new De Bruijn graph with these edges.
        DeBruijnGraph { edges }
    }

    // Method to display the graph, useful for debugging and visualization.
    fn display(&self) {
        // Iterate over each node and its adjacent nodes in the graph.
        for (node, next_nodes) in &self.edges {
            // Print the current node and its connected nodes.
            println!("{} -> {:?}", node, next_nodes);
        }
    }
}


// Function to generate a random DNA sequence of a given length
fn generate_random_dna_sequence(length: usize) -> String {
    // Define the DNA bases
    let bases = ['A', 'C', 'G', 'T'];

    // Initialize a random number generator
    let mut rng = rand::thread_rng();

    // Create a uniform distribution to randomly select from the DNA bases
    let between = Uniform::from(0..bases.len());

    // Generate the DNA sequence by randomly selecting bases
    (0..length).map(|_| bases[between.sample(&mut rng)]).collect()
}

// Function to write the generated DNA sequence to a file
fn write_to_file(contents: &str, file_name: &str) -> io::Result<()> {
    // Create or open the file
    let mut file = File::create(file_name)?;

    // Write the DNA sequence to the file
    file.write_all(contents.as_bytes())?;

    // Return the result of the file operation
    Ok(())
}
// Function to read the DNA sequence from a file
fn read_dna_sequence_from_file(file_name: &str) -> io::Result<String> {
    fs::read_to_string(file_name)
}

// Function to write k-mer counts to a CSV file
fn write_kmer_counts_to_csv(kmer_counts: &HashMap<String, usize>, file_name: &str) -> Result<(), Box<dyn Error>> {
    let mut wtr = Writer::from_path(file_name)?;

    // Write header
    wtr.write_record(&["K-mer", "Count"])?;

    // Write k-mer data
    for (kmer, count) in kmer_counts {
        wtr.write_record(&[kmer, &count.to_string()])?;
    }

    wtr.flush()?;
    Ok(())
}

fn write_graph_to_csv(graph: &DeBruijnGraph, file_name: &str) -> Result<(), Box<dyn Error>> {
    let mut wtr = Writer::from_path(file_name)?;

    // Write header
    wtr.write_record(&["Node", "Connected Nodes"])?;

    // Write graph data
    for (node, edges) in &graph.edges {
        let connected_nodes = edges.join(", ");
        wtr.write_record(&[node, &connected_nodes])?;
    }

    wtr.flush()?;
    Ok(())
}


// Function to plot a histogram of k-mer counts
fn plot_kmer_histogram(kmer_counts: &HashMap<String, usize>, output_file: &str) -> Result<(), Box<dyn std::error::Error>> {
    // Create a drawing area for the plot, specifying the output file and dimensions
    let root_area = BitMapBackend::new(output_file, (640, 480)).into_drawing_area();
    // Fill the drawing area with a white background
    root_area.fill(&WHITE)?;

    // Determine the maximum count for the y-axis
    let max_count = *kmer_counts.values().max().unwrap();

    // Build the chart with specified caption, margins, and axis sizes
    let mut chart = ChartBuilder::on(&root_area)
        .caption("K-mer Frequency Histogram", ("sans-serif", 40))
        .margin(5)
        .x_label_area_size(30)
        .y_label_area_size(40)
        .build_cartesian_2d(0usize..kmer_counts.len(), 0usize..max_count)?;



    // Configure and draw the mesh (grid lines) for the chart
    chart.configure_mesh().draw()?;

    // Convert k-mer count data into a format suitable for the histogram
    let data: Vec<(usize, usize)> = kmer_counts.iter().enumerate().map(|(index, (_kmer, &count))| (index, count)).collect();

    // Draw the histogram series on the chart
    chart.draw_series(
        Histogram::vertical(&chart)
            .style(RED.filled())
            .data(data.iter().map(|&(index, count)| (index, count))),
    )?;


    // Present the final plot
    root_area.present()?;
    Ok(())
}

fn create_petgraph(de_bruijn_graph: &DeBruijnGraph) -> UnGraph<String, ()> {
    let mut graph = UnGraph::<String, ()>::new_undirected();

    let mut index_map = std::collections::HashMap::new();

    for (node, edges) in &de_bruijn_graph.edges {
        let node_index = *index_map.entry(node.clone()).or_insert_with(|| graph.add_node(node.clone()));
        for edge in edges {
            let edge_index = *index_map.entry(edge.clone()).or_insert_with(|| graph.add_node(edge.clone()));
            graph.add_edge(node_index, edge_index, ());
        }
    }

    graph
}

fn save_graph_dot(graph: &UnGraph<String, ()>, file_name: &str) -> Result<(), std::io::Error> {
    let dot = Dot::with_config(&graph, &[Config::EdgeNoLabel]);
    let mut file = File::create(file_name)?;
    writeln!(file, "{:?}", dot)
}


fn main() {
    // Prompt for and read the DNA sequence length from the user
    println!("Enter the length of the DNA sequence:");
    let mut dna_length_str = String::new();
    io::stdin().read_line(&mut dna_length_str).expect("Failed to read line");
    let dna_length: usize = dna_length_str.trim().parse().expect("Please type a number!");

    // Generate a random DNA sequence and save it to a file
    let dna_sequence = generate_random_dna_sequence(dna_length);
    match write_to_file(&dna_sequence, "random_dna_sequence.txt") {
        Ok(_) => println!("DNA sequence saved to random_dna_sequence.txt"),
        Err(e) => eprintln!("Failed to write DNA sequence to file: {}", e),
    }

    // Read the DNA sequence from the file
    let dna_sequence = match read_dna_sequence_from_file("random_dna_sequence.txt") {
        Ok(sequence) => sequence,
        Err(e) => {
            eprintln!("Failed to read DNA sequence from file: {}", e);
            return;
        }
    };

    // Prompt for and read the k-mer size from the user
    println!("Enter the size of k-mer:");
    let mut k_str = String::new();
    io::stdin().read_line(&mut k_str).expect("Failed to read line");
    let k: usize = k_str.trim().parse().expect("Please type a number!");

    // Generate k-mers from the DNA sequence
    let kmers = generate_kmers(&dna_sequence, k);

    // Count the frequency of each k-mer
    let kmer_counts = count_kmers(kmers.clone()); // Clone kmers for further use

    // Plot the k-mer histogram
    match plot_kmer_histogram(&kmer_counts, "kmer_histogram.png") {
        Ok(_) => println!("K-mer histogram plotted in kmer_histogram.png"),
        Err(e) => eprintln!("Failed to plot k-mer histogram: {}", e),
    }
    // Optionally, display k-mer counts
    for (kmer, count) in &kmer_counts {
        println!("{}: {}", kmer, count);
    }

    // Create a De Bruijn graph from the k-mers
    let dbg = DeBruijnGraph::new(&kmers);

    // Display the De Bruijn graph
    dbg.display();
    // Convert to petgraph graph
    let graph = create_petgraph(&dbg);

    // Save graph to DOT file
    if let Err(e) = save_graph_dot(&graph, "de_bruijn_graph.dot") {
        eprintln!("Failed to save graph to DOT file: {}", e);
    }
    // Write k-mer counts to an Excel file
    match write_kmer_counts_to_csv(&kmer_counts, "kmer_counts.csv") {
        Ok(_) => println!("K-mer counts saved to kmer_counts.csv"),
        Err(e) => eprintln!("Failed to write k-mer counts to CSV: {}", e),
    }
    match write_graph_to_csv(&dbg, "de_bruijn_graph.csv") {
        Ok(_) => println!("De Bruijn graph saved to de_bruijn_graph.csv"),
        Err(e) => eprintln!("Failed to write De Bruijn graph to CSV: {}", e),
    }
}
