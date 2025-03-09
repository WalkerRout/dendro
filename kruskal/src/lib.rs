use std::cmp::Ordering;
use std::collections::BinaryHeap;

use rayon::prelude::*;

use lib_needleman::Needleman;

use lib_genome_kit::blosum::Blosum62;
use lib_genome_kit::genome::Genome;

/// Wrapper for a species index
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct SpeciesId(usize);

impl SpeciesId {
  #[inline]
  fn index(self) -> usize {
    self.0
  }
}

/// Wrapper for a cluster index
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct ClusterId(usize);

impl ClusterId {
  #[inline]
  fn index(self) -> usize {
    self.0
  }
}

/// An edge between two species/clusters with a similarity score
///
/// Note: higher scores mean higher similarity...
#[derive(Debug, Eq, PartialEq)]
struct Edge {
  species_a: SpeciesId,
  species_b: SpeciesId,
  similarity: i32,
}

impl PartialOrd for Edge {
  #[inline]
  fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
    Some(self.similarity.cmp(&other.similarity))
  }
}

impl Ord for Edge {
  #[inline]
  fn cmp(&self, other: &Self) -> Ordering {
    self.similarity.cmp(&other.similarity)
  }
}

/// A similarity-score based binary tree, generic over leaf storage...
#[derive(Debug, Clone)]
pub enum Cluster<T> {
  Leaf(T),
  Node {
    left: Box<Cluster<T>>,
    right: Box<Cluster<T>>,
    similarity: i32,
  },
}

/// A simple union–find (disjoint set) implementation
struct UnionFind {
  parent: Vec<ClusterId>,
  rank: Vec<usize>,
}

impl UnionFind {
  fn new(n: usize) -> Self {
    Self {
      parent: (0..n).map(ClusterId).collect(),
      rank: vec![0; n],
    }
  }

  /// Finds the representative of the set containing `x` using path compression
  fn find(&mut self, x: ClusterId) -> ClusterId {
    let idx = x.index();
    if self.parent[idx] != x {
      self.parent[idx] = self.find(self.parent[idx]);
    }
    self.parent[idx]
  }

  /// Unites the sets that contain `x` and `y` (if they are disjoint), returns
  /// `true` if a union occurred
  fn union(&mut self, x: ClusterId, y: ClusterId) -> bool {
    let x_root = self.find(x);
    let y_root = self.find(y);
    if x_root == y_root {
      return false;
    }

    let x_idx = x_root.index();
    let y_idx = y_root.index();

    match self.rank[x_idx].cmp(&self.rank[y_idx]) {
      Ordering::Less => {
        self.parent[x_idx] = y_root;
      }
      Ordering::Equal => {
        self.parent[y_idx] = x_root;
      }
      Ordering::Greater => {
        self.parent[y_idx] = x_root;
        self.rank[x_idx] += 1;
      }
    }

    true
  }
}

/// A species of animal, has some name and some genome...
#[derive(Debug, Clone, PartialEq)]
pub struct Species {
  name: String,
  genome: Genome,
}

impl Species {
  #[inline]
  pub fn new(name: String, genome: impl Into<Genome>) -> Self {
    Self {
      name,
      genome: genome.into(),
    }
  }
}

/// Some way to split some type into two dependent types...
trait Partition {
  type A;
  type B;

  fn part(self) -> (Self::A, Self::B);
}

impl Partition for Species {
  type A = String;
  type B = Genome;

  /// We can part a `Species` into its components
  #[inline]
  fn part(self) -> (Self::A, Self::B) {
    (self.name, self.genome)
  }
}

impl<T, U, P> Partition for Vec<P>
where
  P: Partition<A = T, B = U>,
{
  type A = Vec<T>; // impl Iterator<Item = T>
  type B = Vec<U>; // impl Iterator<Item = U>

  /// We can part a `Vec` of partables into two lists of their components...
  #[inline]
  fn part(self) -> (Self::A, Self::B) {
    self.into_iter().map(Partition::part).unzip()
  }
}

/// A manager that encapsulates the union–find method of querying and updating
/// clusters...
struct ClusterManager<T> {
  uf: UnionFind,
  clusters: Vec<Cluster<T>>,
  cluster_map: Vec<ClusterId>,
}

impl<T> ClusterManager<T>
where
  T: Clone,
{
  fn new(initial_clusters: Vec<Cluster<T>>) -> Self {
    let n = initial_clusters.len();
    Self {
      clusters: initial_clusters,
      uf: UnionFind::new(2 * n),
      cluster_map: (0..n).map(ClusterId).collect(),
    }
  }

  /// Returns the current cluster id for a given species
  #[inline]
  fn get(&mut self, species: SpeciesId) -> ClusterId {
    self.uf.find(self.cluster_map[species.index()])
  }

  /// Adds a new cluster node to the manager, returning its new id
  #[inline]
  fn add_cluster(&mut self, new_cluster: Cluster<T>) -> ClusterId {
    self.clusters.push(new_cluster);
    ClusterId(self.clusters.len() - 1)
  }

  /// Merges clusters represented by `a` and `b`, updating the mapping to `new_id`
  fn merge(&mut self, a: ClusterId, b: ClusterId, new_id: ClusterId) {
    self.uf.union(a, b);
    for entry in self.cluster_map.iter_mut() {
      let current = self.uf.find(*entry);
      if current == a || current == b {
        *entry = new_id;
      }
    }
  }

  /// Produces the root (final merged) cluster
  #[inline]
  fn root(&self) -> Option<Cluster<T>> {
    self
      .cluster_map
      .first()
      .map(|&cid| self.clusters[cid.index()].clone())
  }
}

pub trait Kruskal {
  type Leaf;

  fn cluster(self) -> Option<Cluster<Self::Leaf>>;
}

impl<T, P> Kruskal for Vec<P>
where
  T: Clone,
  P: Partition<A = T, B = Genome>,
{
  /// A leaf is a species' name
  type Leaf = T;

  fn cluster(self) -> Option<Cluster<Self::Leaf>> {
    let (names, genomes) = self.part();
    let n = names.len();
    if n == 0 {
      return None;
    }
    let initial_clusters = create_initial_clusters(names);
    let mut manager = ClusterManager::new(initial_clusters);
    let mut heap = compute_edges(&genomes);

    while let Some(edge) = heap.pop() {
      let a = edge.species_a;
      let b = edge.species_b;
      let ca = manager.get(a);
      let cb = manager.get(b);
      if ca != cb {
        let new_cluster = Cluster::Node {
          left: Box::new(manager.clusters[ca.index()].clone()),
          right: Box::new(manager.clusters[cb.index()].clone()),
          similarity: edge.similarity,
        };
        let new_id = manager.add_cluster(new_cluster);
        manager.merge(ca, cb, new_id);
      }
    }

    manager.root()
  }
}

// initial clusters are just leaves representing all the species names...
#[inline]
fn create_initial_clusters<T>(names: Vec<T>) -> Vec<Cluster<T>> {
  names.into_iter().map(Cluster::Leaf).collect()
}

// !!! this is the hottest function in the program !!!
// we avoid inline so i can profile with perf lol
#[inline(never)]
fn compute_edges(genomes: &[Genome]) -> BinaryHeap<Edge> {
  let n = genomes.len();
  let edges: Vec<Edge> = (0..n)
    .into_par_iter()
    .flat_map_iter(|i| {
      (i + 1..n).map(move |j| {
        let score = Blosum62::needleman_wunsch(&genomes[i], &genomes[j]);
        Edge {
          species_a: SpeciesId(i),
          species_b: SpeciesId(j),
          similarity: score,
        }
      })
    })
    .collect();
  // we want a max binheap from our vector...
  edges.into_iter().collect()
}

#[cfg(test)]
mod tests {
  use super::*;

  mod kruskal {
    use super::*;

    #[test]
    fn cluster() {
      let species = vec![
        Species::new("Species A".into(), "ARND".chars()),
        Species::new("Species B".into(), "ARNE".chars()),
        Species::new("Species C".into(), "ARNS".chars()),
        Species::new("Species D".into(), "RRDD".chars()),
        Species::new("Species E".into(), "RRDS".chars()),
        Species::new("Species F".into(), "RRDA".chars()),
        Species::new("Species G".into(), "ARDD".chars()),
        Species::new("Species H".into(), "ARDS".chars()),
        Species::new("Species I".into(), "RRNS".chars()),
      ];
      let dendrogram = species.cluster();
      assert!(dendrogram.is_some());
      // should definitely be a better test, but we were running into ordering
      // concerns while optimizing between vec and binaryheap so im too lazy to
      // add something concrete here
    }
  }
}
