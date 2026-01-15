use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use serde_json::Result;
use std::collections::HashMap;
use std::collections::HashSet;
use std::mem;
use std::io::Write;
use std::error::Error;
use std::fs::File;
use std::fs;
use std::io::{self, BufRead, BufReader};
use std::path::Path;
use std::time::{Duration, Instant};
use std::sync::{Arc, RwLock};
use rayon::ThreadPoolBuilder;
use clap::Parser;

extern crate kodama;
use kodama::{linkage, Method, Step};
// Define a structure for command line arguments


#[derive(Parser, Debug)]
#[command(version, about, long_about = "Hierarchical clustering of genomes based on Nucleotide Average Identity (ANI)")]
struct Args {
    /// path of Nucleotide Average Identity (ANI) JSON file, format is {"genome1": {"genome2": 0.99, "genome3": 0.98, ...}, "genome2": {"genome1": 0.99, "genome3": 0.97, ...}, ...}
    #[clap(short, long, help="path of Nucleotide Average Identity (ANI) JSON file")]
    ani_file: String,
    /// path of genome list file, each line is a genome name to cluster, format is "genome1\n genome2\n genome3\n ..."
    #[clap(short, long, help="path of genome list file, each line is a genome name to cluster.")]
    glist_file: String,
    /// path of output cluster dictionary JSON file, format is {"cluster1": ["genome1", "genome2"], "cluster2": ["genome3", "genome4"], ...}
    #[clap(short, long, help="path of output for cluster dictionary, Dendrogram and Newick file")]
    outdir: String,
    // distance cutoff 
    #[clap(short, long, default_value_t = 0.05, help="distance cutoff criteria default 0.05.")]
    cutoff: f32,
    /// number of threads to use for parallel processing to create distance matrix
    #[clap(short, long, default_value_t = 1, help="number of threads to use for parallel processing to create distance matrix")]
    threads: usize,
}

/// Tree node that can hold labels

#[derive(Debug)]
pub enum TreeNode {
    Leaf(String),
    Node {
        left: Option<Box<TreeNode>>, // MODIFIED
        right: Option<Box<TreeNode>>, // MODIFIED
        distance: f64,
    },
}

impl Drop for TreeNode {
    fn drop(&mut self) {
        // This implementation correctly works if TreeNode's children are Option<Box<TreeNode>>
        if let TreeNode::Node { ref mut left, ref mut right, .. } = *self {
            // Now `left` and `right` are `&mut Option<Box<TreeNode>>`.
            // So, `left.take()` and `right.take()` will call `Option::take()`.
            let mut stack = Vec::new();

            if let Some(left_child_box) = left.take() { // This will now work
                stack.push(left_child_box);
            }
            if let Some(right_child_box) = right.take() { // This will now work
                stack.push(right_child_box);
            }

            while let Some(mut boxed_node) = stack.pop() { // boxed_node is Box<TreeNode>
                // When boxed_node is processed, its children are added to the stack.
                // Then boxed_node (the Box) goes out of scope and is dropped.
                // Dropping the Box<TreeNode> calls TreeNode::drop on its content.
                if let TreeNode::Node { left: ref mut l_child, right: ref mut r_child, .. } = *boxed_node {
                    // l_child and r_child are &mut Option<Box<TreeNode>>
                    if let Some(child_box) = l_child.take() { // This will now work
                        stack.push(child_box);
                    }
                    if let Some(child_box) = r_child.take() { // This will now work
                        stack.push(child_box);
                    }
                }
            }
        }
        // If self is a Leaf, its String will be dropped automatically.
        // If self was a Node, its left and right are now None (after .take()). The distance f64 is dropped.
    }
}
// ...existing code...
// Define a struct that matches the JSON data structure
#[derive(Serialize, Deserialize, Debug)]
struct NestedDict {
    #[serde(flatten)]
    outer_dict: HashMap<String, HashMap<String, f32>>, // Outer dictionary with inner dictionaries of floats
}

pub struct Linkage<T> {
    pub dendrogram: Vec<Step<T>>,
}


fn welcome() {
    println!("Welcome to the hierarchical clustering of genomes based on Nucleotide Average Identity (ANI)");
}

fn read_txt(filename: impl AsRef<Path>) -> io::Result<Vec<String>> {
    BufReader::new(File::open(filename)?).lines().collect()
}

fn create_folder_if_missing(path: &str) -> std::io::Result<()> {
    let folder = Path::new(path);
    if !folder.exists() {
        fs::create_dir_all(folder)?; // create_dir_all creates parent folders if needed
        println!("Created directory: {}", path);
    } else {
        println!("Directory already exists: {}", path);
    }
    Ok(())
}

fn comb2(n: u32) -> u64 {
    n as u64 * (n as u64 - 1) / 2
}

fn get_ani(ani_info: &NestedDict, qry1: String, qry2: String) -> f32 {
    let mut ani: f32 = 0.0;
    if let Some(in_d) = ani_info.outer_dict.get(&qry1) {
        if let Some(mani) = in_d.get(&qry2) {
            ani = mani.clone();
        } else if let Some(in_d) = ani_info.outer_dict.get(&qry2) {
            if let Some(mani) = in_d.get(&qry1) {
                ani = mani.clone();
            }
        }
    } else if let Some(in_d) = ani_info.outer_dict.get(&qry2) {
        if let Some(mani) = in_d.get(&qry1) {
            ani = mani.clone();
        } else if let Some(in_d) = ani_info.outer_dict.get(&qry1) {
            if let Some(mani) = in_d.get(&qry2) {
                ani = mani.clone();
            }
        }
    }
    ani
}

fn get_darray_parallel(
    ani_info: &Arc<NestedDict>,
    glist_chunks: Vec<Vec<String>>,
    darray: Arc<RwLock<Vec<f32>>>,
    g2idx: &Arc<HashMap<String, u32>>,
) {
    let glist_ln = g2idx.len() as u32;
    glist_chunks.into_par_iter().for_each(|chunk| {
        for g1 in &chunk {
            for g2 in g2idx.keys() {
                if g1 == g2 {
                    continue;
                }
                let idx1 = g2idx.get(g1).unwrap();
                let idx2 = g2idx.get(g2).unwrap();
                let dist: f32 = 1.0 - get_ani(&ani_info, g1.clone(), g2.clone());
                if idx1 >= idx2 {
                    continue;
                }
                let comb2_ln1 = comb2(glist_ln);
                let comb2_ln2 = comb2(glist_ln - idx1);
                let comb2_ln3 = (idx2 - idx1 - 1) as u64;

                let da_i = comb2_ln1 - comb2_ln2 + comb2_ln3;

                let mut darray_lock = darray.write().unwrap();
                darray_lock[da_i as usize] = dist;
            }
        }
    });

    let is_present = darray.read().unwrap().contains(&-1.0);
    if is_present {
        println!("Some values are missing in the distance array");
    } else {
        println!("All values are present in the distance array");
    }
}

fn check_duplicates<K: Eq + std::hash::Hash, V>(
    map1: &HashMap<K, V>,
    map2: &HashMap<K, V>,
) -> Vec<K>
where
    K: Clone,
{
    map2.keys()
        .filter(|key| map1.contains_key(*key))
        .cloned()
        .collect()
}

fn get_cludict(
    dend: &kodama::Dendrogram<f32>,
    cutoff: f32,
    glist_ln: u32,
    g2idx: &HashMap<String, u32>,
) -> HashMap<u32, HashSet<String>> {
    let idx2g: HashMap<u32, String> = g2idx.iter().map(|(k, v)| (*v, k.clone())).collect();
    let mut clusters: HashMap<u32, HashSet<String>> = HashMap::new();
    let mut merged: HashMap<u32, HashSet<String>> = HashMap::new();
    let valid_index: u32 = (glist_ln - 1) as u32;
    let mut accu_index = glist_ln - 1;

    for step in dend.steps() {
        let cluster1 = step.cluster1 as u32;
        let cluster2 = step.cluster2 as u32;
        let dissimilarity = step.dissimilarity as f32;

        if dissimilarity <= cutoff {
            accu_index += 1;
            if cluster1 <= valid_index {
                if cluster2 <= valid_index {
                    let mut new_cluster = HashSet::new();
                    new_cluster.insert(idx2g[&cluster1].clone());
                    new_cluster.insert(idx2g[&cluster2].clone());
                    merged.insert(accu_index, new_cluster);
                } else {
                    let mut new_cluster = HashSet::new();
                    new_cluster.insert(idx2g[&cluster1].clone());
                    new_cluster = new_cluster
                        .union(merged.get(&cluster2).unwrap())
                        .cloned()
                        .collect();
                    merged.insert(accu_index, new_cluster);
                    merged.remove(&cluster2);
                }
            } else {
                if cluster2 <= valid_index {
                    let mut new_cluster = HashSet::new();
                    new_cluster.insert(idx2g[&cluster2].clone());
                    new_cluster = new_cluster
                        .union(merged.get(&cluster1).unwrap())
                        .cloned()
                        .collect();
                    merged.insert(accu_index, new_cluster);
                    merged.remove(&cluster1);
                } else {
                    let mut new_cluster = merged.get(&cluster1).unwrap().clone();
                    new_cluster = new_cluster
                        .union(merged.get(&cluster2).unwrap())
                        .cloned()
                        .collect();
                    merged.insert(accu_index, new_cluster);
                    merged.remove(&cluster1);
                    merged.remove(&cluster2);
                }
            }
        } else {
            if cluster1 <= valid_index {
                clusters.insert(cluster1, HashSet::from([idx2g[&cluster1].clone()]));
                if cluster2 <= valid_index {
                    clusters.insert(cluster2, HashSet::from([idx2g[&cluster2].clone()]));
                }
            }
        }
    }

    let duplicates = check_duplicates(&clusters, &merged);
    if duplicates.is_empty() {
        println!("No duplicates found in the clusters");
    } else {
        println!("Duplicates found in the clusters");
    }

    clusters.extend(merged.into_iter());
    clusters = clusters
        .into_iter()
        .map(|(k, v)| (k, v.into_iter().collect::<HashSet<String>>()))
        .collect::<HashMap<u32, HashSet<String>>>();
    clusters
}

fn build_labeled_tree(
    dend: &kodama::Dendrogram<f32>,
    g2idx: &HashMap<String, u32>,
    glist_ln: u32,
) -> TreeNode {
    let idx2g: HashMap<u32, String> = g2idx.iter().map(|(k, v)| (*v, k.clone())).collect();
    let mut nodes: HashMap<u32, TreeNode> = idx2g.iter().map(|(&i, label)| (i, TreeNode::Leaf(label.clone()))).collect();
    // ...existing code...
    let mut next_index = glist_ln as u32;

    for step in dend.steps() {
        let cluster1 = step.cluster1 as u32;
        let cluster2 = step.cluster2 as u32;
        let dissimilarity = step.dissimilarity as f32;

        let left_node = nodes.remove(&cluster1).unwrap();
        let right_node = nodes.remove(&cluster2).unwrap();

        nodes.insert(
            next_index,
            TreeNode::Node {
                left: Some(Box::new(left_node)), // MODIFIED: Wrap in Some()
                right: Some(Box::new(right_node)), // MODIFIED: Wrap in Some()
                distance: dissimilarity as f64,
            },
        );
        next_index += 1;
    }

    // The final tree node is the root
    nodes.into_iter().next().unwrap().1
}

fn to_newick(node: &TreeNode) -> String {
    let denom = 2.0 as f64;
    match node {
        TreeNode::Leaf(label) => label.clone(),
        TreeNode::Node { left, right, distance } => {
            // Need to handle Option here
            let l_str = left.as_ref().map_or_else(|| String::from(""), |l_box| to_newick(l_box));
            let r_str = right.as_ref().map_or_else(|| String::from(""), |r_box| to_newick(r_box));
            format!("({}:{},{}:{})", l_str, distance/denom, r_str, distance/denom)
        }
    }
}


fn main() -> Result<()> {
    // Open the file
    welcome();
    let Initialization = Instant::now();
    let args = Args::parse();
    let ani_file = args.ani_file;
    let glist_file = args.glist_file;
    let out_dir = args.outdir;
    let threads: usize = args.threads;
    let cut_off: f32 = args.cutoff;
    println!("NOTICE :: Starting the hierarchical clustering of genomes based on Nucleotide Average Identity (ANI)");
    let start_time = Instant::now();
    println!("NOTICE :: Loading all the necessary files");
    let ani_d: File = File::open(&ani_file).map_err(serde_json::Error::io)?;
    println!("NOTICE :: Trying to load ANI dictionary file");
    // Deserialize the JSON content into a ani_info struct
    let ani_info: NestedDict = serde_json::from_reader(ani_d)?;
    println!("NOTICE :: ANI dictionary file has been successfully loaded");
    // Create a buffered reader to read the file line by line
    let glist: Vec<String> = read_txt(&glist_file).expect("Could not load lines");
    let glist_ln = glist.len() as u32;
    // mapping numeric index to genome name
    let mut g2idx: HashMap<String, u32> = HashMap::new();
    for i in 0..glist.len() {
        g2idx.insert(glist[i].clone(), i as u32);
    }
    let duration_load = start_time.elapsed().as_secs() / 60;
    println!("NOTICE :: Genomes list file has been successfully loaded and Number of genoms to indexed: {} (Elapsed time {:?} min", g2idx.len(), duration_load);
    
    // Calculate the number of combinations
    let da_dim = comb2(glist_ln);
    let da_size_bytes: usize = da_dim as usize * std::mem::size_of::<f32>();
    // let da_size_bytes: usize = da_dim as usize * std::mem::size_of::<f8>();
    let start_to_alloc = Instant::now();
    println!(
        "NOTICE :: Checking and Attempting to allocate distance array for {} genomes, estimated size is (~{:.3} GiB)",
        glist.len(),
        da_size_bytes as f32 / 1024.0 / 1024.0 / 1024.0
    );

    // Try to create the vector for distance array
    // Try to allocate the vector (this will likely fail on most systems due to memory limits)

    let alloc_result = std::panic::catch_unwind(|| {
        let mut darray: Vec<f32> = Vec::with_capacity(da_dim as usize);

        // Fill the vector with -1 values
        let base_value: f32 = -1.0;
        for _ in 0..da_dim {
           darray.push(base_value);
        }

        println!(
            "NOTICE :: Successfully allocated memory for the vector with {} elements",
            darray.len()
        );
        // Start filling the distance array with distance values with multithreading
        //Constraint: The number of threads
        let pool = ThreadPoolBuilder::new().num_threads(threads).build().unwrap();
        // 1.split genome list in to chunks
        let glist_chunks: Vec<Vec<String>> =
            glist.chunks(threads).map(|chunk| chunk.to_vec()).collect();
        // 2. create a thread based on chunk size
        let ani_info = Arc::new(ani_info);
        let darray = Arc::new(RwLock::new(darray));
        let g2idx = Arc::new(g2idx);
        let fill_duration = start_to_alloc.elapsed().as_secs() / 60;
        pool.install(|| {
            get_darray_parallel(&ani_info, glist_chunks, darray.clone(), &g2idx);
            println!("NOTICE :: Successfully filled the distance array with distance value (Elapsed time {:?} min", fill_duration);
        });
        // Perform hierarchical clustering
        let mut darray_mut = darray.write().unwrap();
        let darray_vec: &mut Vec<f32> = &mut *darray_mut;
        // Create a dendrogram from the distance array
        let clus_time = Instant::now();
        println!("NOTICE :: Creating the dendrogram from the distance array");
        let dend = linkage(darray_vec, glist_ln as usize, Method::Average);
        println!("NOTICE :: Successfully created the dendrogram with {} elements", dend.len());
        
        assert_eq!(dend.len(), glist_ln as usize - 1);
        // println!("Dendrogram: {:?}", dend);

        println!("NOTICE :: Creating the cluster dictionary from the dendrogram");
        let clu_dict = get_cludict(&dend, cut_off, glist_ln, &g2idx);
        fs::create_dir_all(&out_dir).expect("Unable to create output directory");
        let out_clud = Path::new(&out_dir).join("clus_dict.json");
        let json_file = File::create(&out_clud).expect("Unable to create file");
        serde_json::to_writer(json_file, &clu_dict).expect("Unable to write data");
        println!("NOTICE :: Successfully created the JSON file for cluster dictionary in {}", &out_clud.display());
        let dend_file = Path::new(&out_dir).join("dendrogram.txt");
        let mut dfile = File::create(&dend_file).expect("Unable to create file");
        for step in dend.steps() {
            writeln!(dfile, "{:?}", step).expect("Unable to write data");
        }
        println!("NOTICE :: Successfully created the dendrogram file in {}", &dend_file.display());
        let clus_duration = clus_time.elapsed().as_secs() / 60;
        println!("NOTICE :: Successfully created the cluster dictionary with {} elements (Elapsed time {:?} min", clu_dict.len(), clus_duration);
        //println!("Cluster dictionary: {:?}", clu_dict);
        println!("NOTICE :: Creating the Phylogentic Newick format from the dendrogram");
        let tree = build_labeled_tree(&dend, &g2idx, glist_ln);
        let newick = to_newick(&tree);
        println!("NOTICE :: Successfully created the Phylogentic Newick format");
        create_folder_if_missing(&out_dir).expect("Unable to create directory");

        let tree_file = Path::new(&out_dir).join("tree.newick");

        println!("NOTICE :: Saving the cluster dictionary to a JSON file");

        let mut file = File::create(&tree_file).expect("Unable to create file");
        file.write_all(newick.as_bytes()).expect("Unable to write data");

        println!("NOTICE :: Successfully created the Newick file in {}", &tree_file.display());

    });
    let total_duration = Initialization.elapsed().as_secs() / 60;
    println!("NOTICE :: ALL modules have been successfully completed without any errors (Total elapsed time {:?} min", total_duration);

    // Check if the allocation was successful or if we encountered an error
    match alloc_result {
        Ok(_) => println!("Condensed distance array for {} genoems with {} genome pairs has been successfully created", glist.len(), da_dim),
        Err(_) => println!("Failed to allocate the vector. Memory limit exceeded or other issue."),
    }

    Ok(())
}
