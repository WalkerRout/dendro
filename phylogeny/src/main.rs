use std::collections::HashMap;
use std::fmt::Display;
use std::fs;
use std::path::Path;

use anyhow::Context;

use lib_kruskal::{Cluster, Kruskal, Species};

// we will be anaylzing COX3 (Cytochrome c oxidase subunit III)
fn load_species_from_file(path: impl AsRef<Path>) -> Result<Vec<Species>, anyhow::Error> {
  let json_str = fs::read_to_string(path)?;
  let data: HashMap<String, String> = serde_json::from_str(&json_str)?;
  let species = data
    .into_iter()
    .map(|(k, v)| Species::new(k, v.chars()))
    .collect();
  Ok(species)
}

fn emit_graphviz<T: Display>(cluster: Cluster<T>) -> String {
  // we need to create some nodes with some ids, so we better keep track of what
  // id we are on...
  fn traverse<T: Display>(cluster: &Cluster<T>, counter: &mut usize, output: &mut String) -> usize {
    let node_id = *counter;
    *counter += 1;
    match cluster {
      Cluster::Leaf(val) => {
        // our leaves have a box shape to distinguish them
        output.push_str(&format!(
          "  node{} [label=\"{}\", shape=box];\n",
          node_id, val
        ));
      }
      Cluster::Node {
        left,
        right,
        similarity,
      } => {
        // define an internal node for merge similarity
        output.push_str(&format!(
          "  node{} [label=\"sim: {}\"];\n",
          node_id, similarity
        ));
        // find kids
        let left_id = traverse(left, counter, output);
        let right_id = traverse(right, counter, output);
        // connect parent to kids
        output.push_str(&format!("  node{} -> node{};\n", node_id, left_id));
        output.push_str(&format!("  node{} -> node{};\n", node_id, right_id));
      }
    }
    node_id
  }

  // count unique node ids
  let mut counter = 0;

  let mut output = String::new();
  output.push_str("digraph ClusterTree {\n");
  output.push_str("  node [fontname=\"Helvetica\"];\n");

  // traverse from root
  traverse(&cluster, &mut counter, &mut output);

  output.push_str("}\n");
  output
}

// we dont really need this now but we might do some stuff before we save later...
fn write_graph_to_file(graph_dot: String) -> Result<(), anyhow::Error> {
  // save graph dot with .gv extension
  fs::write("phylogeny.gv", graph_dot)?;
  Ok(())
}

fn main() -> Result<(), anyhow::Error> {
  let species = load_species_from_file("pull-species/cox3_translations.json")?;
  // the root `Cluster` represents a dendrogram with all species...
  let dendrogram = species.cluster().context("no graph found")?;
  let graph = emit_graphviz(dendrogram);
  write_graph_to_file(graph)?;
  Ok(())
}
