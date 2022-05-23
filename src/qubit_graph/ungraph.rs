use itertools::Itertools;
use rand::{rngs::SmallRng, Rng, SeedableRng};
use std::cell::Cell;
use std::collections::HashMap;
use std::rc::Rc;

#[derive(Debug)]
struct ClassicalRegister {
    coord_to_index: HashMap<(i32, i32, i32), usize>,
    index_to_coord: HashMap<(usize, usize), (i32, i32, i32)>,
    register: Vec<Vec<Rc<Cell<u8>>>>,
}

impl ClassicalRegister {
    fn new(round: usize) -> Self {
        let coord_to_index = HashMap::new();
        let index_to_coord = HashMap::new();
        let register = vec![Vec::new(); round + 1];

        Self {
            coord_to_index,
            index_to_coord,
            register,
        }
    }

    /// set Rc<Cell<u8>>
    fn set_classical_register(&mut self, node: (i32, i32, i32), classical_register: Rc<Cell<u8>>) {
        self.index_to_coord
            .insert((self.register[node.2 as usize].len(), node.2 as usize), node);
        self.coord_to_index
            .insert(node, self.register[node.2 as usize].len());
        self.register[node.2 as usize].push(classical_register);
    }

    /// get register
    fn register_from_coord(&self, node: &(i32, i32, i32)) -> Option<&Rc<Cell<u8>>> {
        self.register
            .get(node.2 as usize)?
            .get(*self.coord_to_index.get(node)?)
    }

    /// get register vec mut
    fn register_mut(&mut self) -> &mut Vec<Vec<Rc<Cell<u8>>>> {
        &mut self.register
    }

    /// iter
    fn iter(&self) -> Iter {
        Iter {
            register: &self.register,
            hash_iter: self.coord_to_index.iter(),
        }
    }

    /// iterate only value
    fn iter_value(&self) -> core::slice::Iter<Vec<Rc<Cell<u8>>>> {
        self.register.iter()
    }

    /// return register
    fn register(&self) -> &Vec<Vec<Rc<Cell<u8>>>> {
        &self.register
    }

    /// return
    fn index_to_coord(&self) -> &HashMap<(usize, usize), (i32, i32, i32)> {
        &self.index_to_coord
    }
}

pub struct Iter<'a> {
    register: &'a Vec<Vec<Rc<Cell<u8>>>>,
    hash_iter: std::collections::hash_map::Iter<'a, (i32, i32, i32), usize>,
}

impl<'a> Iterator for Iter<'a> {
    type Item = (&'a (i32, i32, i32), &'a Rc<Cell<u8>>);

    #[cfg_attr(feature = "inline-more", inline)]
    fn next(&mut self) -> Option<Self::Item> {
        match self.hash_iter.next() {
            Some((coord, index)) => {
                let register = &self.register[coord.2 as usize][*index];
                Some((coord, register))
            }
            None => None,
        }
    }
}

#[derive(Debug)]
pub struct UnGraph {
    network: HashMap<(i32, i32, i32), Vec<(i32, i32, i32)>>,
    edge_weight: HashMap<((i32, i32, i32), (i32, i32, i32)), f64>,
    node_is_boundary: HashMap<(i32, i32, i32), bool>,
    classical_register: ClassicalRegister,
    rng: rand::rngs::SmallRng,
}

impl UnGraph {
    pub fn new(round: usize, seed: u64) -> Self {
        let network = HashMap::new();
        let edge_weight = HashMap::new();
        let node_is_boundary = HashMap::new();
        let classical_register = ClassicalRegister::new(round);
        let rng = SmallRng::seed_from_u64(seed);

        Self {
            network,
            edge_weight,
            node_is_boundary,
            classical_register,
            rng,
        }
    }

    /// make graph from edges
    pub fn from_edges(
        edges: &[((i32, i32, i32), (i32, i32, i32))],
        round: usize,
        seed: u64,
    ) -> Self {
        let mut network = HashMap::new();
        let classical_register = ClassicalRegister::new(round);

        for &(u, v) in edges {
            network.entry(u).or_insert_with(Vec::new).push(v);
            network.entry(v).or_insert_with(Vec::new).push(u);
        }

        let edge_weight = HashMap::new();
        let node_is_boundary = HashMap::new();
        let rng = SmallRng::seed_from_u64(seed);

        Self {
            network,
            edge_weight,
            node_is_boundary,
            classical_register,
            rng,
        }
    }

    /// add edges from vec
    pub fn add_edges_from(&mut self, edges: &[((i32, i32, i32), (i32, i32, i32))]) {
        for &(u, v) in edges {
            self.network.entry(u).or_insert_with(Vec::new).push(v);
            self.network.entry(v).or_insert_with(Vec::new).push(u);
        }
    }

    /// add edge
    pub fn add_edge_from(&mut self, edge: &((i32, i32, i32), (i32, i32, i32))) {
        self.network
            .entry(edge.0)
            .or_insert_with(Vec::new)
            .push(edge.1);
        self.network
            .entry(edge.1)
            .or_insert_with(Vec::new)
            .push(edge.0);
    }

    /// set all edge weight equal p
    pub fn set_all_edge_weight(&mut self, p: f64) {
        let Self {
            network,
            edge_weight,
            ..
        } = self;

        for (u, value) in network.iter() {
            for v in value.iter() {
                edge_weight.insert((*u, *v), p);
            }
        }
    }

    /// set edge weight
    pub fn set_edge_weight(&mut self, edge: &((i32, i32, i32), (i32, i32, i32)), p: f64) {
        self.edge_weight.insert(*edge, p);
        self.edge_weight.insert((edge.1, edge.0), p);
    }

    /// set edges weight equal p
    pub fn set_edges_weight(&mut self, edges: &[((i32, i32, i32), (i32, i32, i32))], p: f64) {
        for &(u, v) in edges.iter() {
            self.edge_weight.insert((u, v), p);
            self.edge_weight.insert((v, u), p);
        }
    }

    /// set is_boundary
    pub fn set_is_boundary(&mut self, node: (i32, i32, i32), is_boundary: bool) {
        self.node_is_boundary.insert(node, is_boundary);
    }

    /// set all is_boundary
    pub fn set_all_is_boundary(&mut self, nodes: &[(i32, i32, i32)], is_boundary: bool) {
        for &n in nodes {
            self.node_is_boundary.insert(n, is_boundary);
        }
    }

    /// return is_boundary
    pub fn is_boundary(&self, coord: &(i32, i32, i32)) -> Option<bool> {
        self.node_is_boundary.get(coord).cloned()
    }

    /// set classical register
    pub fn set_classical_register(&mut self, node: (i32, i32, i32), register: Rc<Cell<u8>>) {
        self.classical_register
            .set_classical_register(node, register);
    }

    /// flip classical register
    pub fn flip_classical_register(&self, coord: &(i32, i32, i32), value: u8) {
        self.classical_register
            .register_from_coord(coord)
            .unwrap()
            .set(value);
    }

    /// iterate classical register
    pub fn iter_classical_register(&self) -> Iter {
        self.classical_register.iter()
    }

    /// return classical register
    pub fn get_register(&self, node: &(i32, i32, i32)) -> Option<&Rc<Cell<u8>>> {
        self.classical_register.register_from_coord(node)
    }

    /// get classical register mut
    pub fn classical_register_mut(&mut self) -> &mut Vec<Vec<Rc<Cell<u8>>>> {
        self.classical_register.register_mut()
    }

    /// get classical register 
    pub fn classical_register(&self) -> &Vec<Vec<Rc<Cell<u8>>>> {
        self.classical_register.register()
    }

    /// tales one parameter, iterating the nodes connected to it by an outbound edge
    pub fn neighbors(&self, node: &(i32, i32, i32)) -> Option<&Vec<(i32, i32, i32)>> {
        self.network.get(node)
    }

    /// return edge weight
    pub fn edge_weight(&self, edge: &((i32, i32, i32), (i32, i32, i32))) -> Option<f64> {
        self.edge_weight.get(edge).cloned()
    }

    /// Iterates all nodes
    pub fn nodes(&self) -> std::collections::hash_map::Keys<(i32, i32, i32), Vec<(i32, i32, i32)>> {
        self.network.keys()
    }

    /// Returns the number of edges
    pub fn size(&self) -> usize {
        let mut size = 0;
        for value in self.network.values() {
            size += value.len();
        }
        size / 2
    }

    /// Returns the number of nodes
    pub fn order(&self) -> usize {
        self.network.keys().len()
    }

    /// 前の時間とxor
    pub fn xor_to_last_time(&mut self) {
        let mut temp_mat: Vec<Vec<u8>> = vec![
            vec![0; self.classical_register.register()[0].len()];
            self.classical_register.register().len() - 2
        ];
        for ((back, forth), temp_vec) in self
            .classical_register_mut()
            .iter()
            .tuple_windows()
            .zip(temp_mat.iter_mut())
        {
            forth
                .iter()
                .zip(back.iter())
                .zip(temp_vec.iter_mut())
                .for_each(|((f, b), temp)| *temp = f.get() ^ b.get());
        }

        // replace
        for (register, temp) in self
            .classical_register_mut()
            .iter()
            .skip(1)
            .flatten()
            .zip(temp_mat.into_iter().flatten())
        {
            register.set(temp);
        }
    }

    /// get edge and return qubit coord
    pub fn edge_to_qubit(edge: ((i32, i32, i32), (i32, i32, i32))) -> (i32, i32) {
        let u = edge.0;
        let v = edge.1;

        debug_assert!((u.0 + v.0) % 2 == 0);
        debug_assert!((u.1 + v.1) % 2 == 0);

        ((u.0 + v.0) / 2, (u.1 + v.1) / 2)
    }

    /// show all defect
    pub fn show_all_defect(&self) {
        print!("defect");
        for (coord, register) in self.classical_register.iter() {
            if register.get() == 1 {
                print!("{:?}, ", coord);
            }
        }
        println!();
    }

    /// reset all register
    pub fn reset_register(&self) {
        for value in self.classical_register.iter_value().flatten() {
            value.set(0);
        }
    }

    pub fn index_to_coord(&self) -> &HashMap<(usize, usize), (i32, i32, i32)> {
        self.classical_register.index_to_coord()
    }
}
