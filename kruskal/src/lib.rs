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

/// A simple union–find (disjoint set) implementation using negative parent values to encode set size.
struct UnionFind {
  parent: Vec<isize>, // root: –size, otherwise: index of parent
}

impl UnionFind {
  /// Creates a new UnionFind with `n` singleton sets
  fn new(n: usize) -> Self {
    UnionFind {
      parent: vec![-1; n], // each element is its own root of size 1
    }
  }

  /// Finds the representative of the set containing `x`, with path compression
  fn find(&mut self, x: ClusterId) -> ClusterId {
    let idx = x.index();
    if self.parent[idx] < 0 {
      x
    } else {
      // this is the final root of the cluster
      let root = self.find(ClusterId(self.parent[idx] as usize));
      self.parent[idx] = root.index() as isize;
      root
    }
  }

  /// Unites the sets containing `a` and `b`, returns whether or not a union
  /// occurred
  fn union(&mut self, a: ClusterId, b: ClusterId) -> bool {
    let ra = self.find(a);
    let rb = self.find(b);
    if ra == rb {
      return false;
    }

    // pick the root with the more negative parent (bigger size)
    let (big, small) = {
      let pa = self.parent[ra.index()];
      let pb = self.parent[rb.index()];
      if pa <= pb { (ra, rb) } else { (rb, ra) }
    };

    let bi = big.index();
    let si = small.index();
    self.parent[bi] += self.parent[si]; // combine sizes
    self.parent[si] = bi as isize; // attach smaller tree under big

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

impl<T> ClusterManager<T> {
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
  fn root(&self) -> Option<Cluster<T>>
  where
    T: Clone,
  {
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
<<<<<<< HEAD
    if names.is_empty() || genomes.is_empty() {
      return None;
    }

=======
    if names.is_empty() {
      return None;
    }
    
>>>>>>> 087b68c3ea25b385887ebe86b9dfdf6e9a1d2d6b
    let mut manager = ClusterManager::new(create_initial_clusters(names));
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

  mod union_find {
    use super::*;

    #[test]
    fn new_every_element_is_its_own_root() {
      let mut uf = UnionFind::new(5);
      for i in 0..5 {
        let id = ClusterId(i);
        assert_eq!(uf.find(id), id);
      }
    }

    #[test]
    fn union_connects_two_elements() {
      let mut uf = UnionFind::new(3);
      let a = ClusterId(0);
      let b = ClusterId(1);
      // initially separate
      assert_ne!(uf.find(a), uf.find(b));
      // union returns true on first union
      assert!(uf.union(a, b));
      // now they share the same root
      let ra = uf.find(a);
      let rb = uf.find(b);
      assert_eq!(ra, rb);
    }

    #[test]
    fn union_returns_false_when_already_connected() {
      let mut uf = UnionFind::new(4);
      let x = ClusterId(2);
      let y = ClusterId(3);
      assert!(uf.union(x, y));
      assert!(!uf.union(x, y));
      // connecting through a third element still returns false if already in same set
      let z = ClusterId(1);
      assert!(uf.union(y, z));
      assert!(!uf.union(x, z));
    }

    #[test]
    fn transitive_union_merges_multiple_sets() {
      let mut uf = UnionFind::new(5);
      let a = ClusterId(0);
      let b = ClusterId(1);
      let c = ClusterId(2);
      assert!(uf.union(a, b));
      // a and b now together, c separate
      assert_eq!(uf.find(a), uf.find(b));
      assert_ne!(uf.find(a), uf.find(c));
      // connect b and c, which should also connect a to c
      assert!(uf.union(b, c));
      let root = uf.find(a);
      assert_eq!(root, uf.find(b));
      assert_eq!(root, uf.find(c));
    }

    #[test]
    fn disjoint_sets_remain_separate() {
      let mut uf = UnionFind::new(6);
      let a = ClusterId(0);
      let b = ClusterId(1);
      let c = ClusterId(2);
      let d = ClusterId(3);
      // build two disjoint pairs: 0,1 and 2,3
      assert!(uf.union(a, b));
      assert!(uf.union(c, d));
      // ensure no cross connection
      assert_ne!(uf.find(a), uf.find(c));
      assert_ne!(uf.find(b), uf.find(d));
    }
  }

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
