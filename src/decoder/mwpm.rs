use crate::qubit_graph::ungraph::UnGraph;

use itertools::Itertools;
use petgraph::graphmap::GraphMap;
use petgraph::graphmap::UnGraphMap;
use petgraph::Undirected;
use std::cmp::Ordering;
use std::collections::BinaryHeap;

use retworkx_core::max_weight_matching::max_weight_matching;
use retworkx_core::petgraph as rpet;
use retworkx_core::Result;

use hashbrown::{HashMap, HashSet};

#[derive(Clone, PartialEq, Debug)]
struct State {
    distance: f64,
    coord: (i32, i32, i32),
    predecessor: (i32, i32, i32),
}

impl Eq for State {}

impl Ord for State {
    fn cmp(&self, other: &Self) -> Ordering {
        other
            .distance
            .partial_cmp(&self.distance)
            .unwrap()
            .then_with(|| Ordering::Less)
    }
}

impl PartialOrd for State {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

pub fn local_dijkstra(
    graph: &UnGraph,
    m: usize,
    s: &(i32, i32, i32),
) -> Vec<(Vec<(i32, i32, i32)>, f64)> {
    let mut state_list = graph
        .nodes()
        .map(|n| {
            (
                *n,
                State {
                    distance: f64::MAX,
                    coord: *n,
                    predecessor: *n,
                },
            )
        })
        .collect::<HashMap<_, _>>();

    state_list.get_mut(&s).unwrap().distance = 0.0;
    let mut heap = BinaryHeap::new();
    heap.push(State {
        distance: 0.0,
        coord: *s,
        predecessor: *s,
    });

    let mut m_nearest_node = Vec::new();
    let mut visited = HashSet::new();

    while m_nearest_node.len() < m + 1 {
        if let Some(State {
            distance,
            coord,
            predecessor,
        }) = heap.pop()
        {
            if visited.contains(&coord) {
                continue;
            }
            if graph
                .get_register(&coord)
                .unwrap_or_else(|| panic!("{:?} isn't exist", coord))
                .get()
                == 1
            {
                m_nearest_node.push(State {
                    distance,
                    coord,
                    predecessor,
                });
            }
            for v in graph.neighbors(&coord).unwrap() {
                let new_distance = graph.edge_weight(&(coord, *v)).unwrap() + distance;
                if new_distance < state_list.get(v).unwrap().distance {
                    let state = state_list.get_mut(v).unwrap();
                    state.distance = new_distance;
                    state.predecessor = coord;
                    heap.push(State {
                        distance: new_distance,
                        coord: *v,
                        predecessor: coord,
                    });
                }
            }
            visited.insert(coord);
        } else {
            break;
        }
    }

    // make m path
    let mut m_nearest_path = Vec::new();
    for n in m_nearest_node.into_iter().skip(1) {
        let mut path = vec![n.coord];
        let mut coord = n.coord;
        loop {
            let next = state_list.get(&coord).unwrap();
            if next.coord == next.predecessor {
                break;
            }

            path.push(next.predecessor);
            coord = next.predecessor;
        }
        m_nearest_path.push((path, n.distance));
    }

    m_nearest_path
}

fn construct_syndrome_graph(
    graph: &UnGraph,
    m: usize,
) -> (
    GraphMap<(i32, i32, i32), f64, Undirected>,
    HashMap<((i32, i32, i32), (i32, i32, i32)), Vec<(i32, i32, i32)>>,
) {
    let mut paths = Vec::new();
    let mut path_detail = HashMap::new();

    for (&coord, defect) in graph.iter_classical_register() {
        if defect.get() == 1 {
            let path = local_dijkstra(graph, m, &coord);

            for (p, d) in path.into_iter() {
                paths.push((coord, p[0], -d));
                path_detail.insert((coord, p[0]), p.clone());
                path_detail.insert((p[0], coord), p);
            }
        }
    }

    // construct local matching graph
    let local_graph = UnGraphMap::<(i32, i32, i32), f64>::from_edges(&paths);

    // dbg!(&local_graph);

    (local_graph, path_detail)
}

fn minimum_weight_perfect_matching(
    local_graph: GraphMap<(i32, i32, i32), f64, Undirected>,
) -> Vec<((i32, i32, i32), (i32, i32, i32))> {
    let mut coord_to_index = HashMap::new();
    let mut index_to_coord = HashMap::new();

    for (i, coord) in local_graph.nodes().enumerate() {
        coord_to_index.insert(coord, i);
        index_to_coord.insert(i, coord);
    }
    let mut edges = Vec::new();

    for (u, v, &w) in local_graph.all_edges() {
        let &start = coord_to_index.get(&u).unwrap();
        let &end = coord_to_index.get(&v).unwrap();

        let weight = (w * 100000000.0) as i128;
        edges.push((start as u32, end as u32, weight));
    }

    let g = rpet::graph::UnGraph::<u32, i128>::from_edges(&edges);

    /*
    let mut f = File::create("example.dot").unwrap();
            let output = format!("{:?}", Dot::new(&g));
            f.write_all(&output.as_bytes())
                .expect("could not write file"); */

    let res: Result<HashSet<(usize, usize)>> =
        max_weight_matching(&g, true, |e| Ok(*e.weight()), true);

    let matching_index = res.unwrap();

    let mut matching = Vec::new();

    for (u, v) in matching_index {
        let &start = index_to_coord.get(&u).unwrap();
        let &end = index_to_coord.get(&v).unwrap();
        matching.push((start, end));
    }

    matching
}

fn decide_correction_qubit(
    graph: &UnGraph,
    matching: Vec<((i32, i32, i32), (i32, i32, i32))>,
    path_detail: HashMap<((i32, i32, i32), (i32, i32, i32)), Vec<(i32, i32, i32)>>,
) -> Vec<(i32, i32)> {
    let mut correction_qubit = Vec::new();

    // 空間方向にedgeが存在するものだけを抽出
    for (u, v) in matching.into_iter() {
        if cfg!(debug_assertions) {
            println!("edge {:?}, {:?}", u, v);
        }
        if (u.0 != v.0) || (u.1 != v.1) {
            let correction_path = path_detail
                .get(&(u, v))
                .unwrap_or_else(|| panic!("edge: {:?} is not exist", (u, v)));

            if cfg!(debug_assertions) {
                let weight: f64 = correction_path
                    .iter()
                    .tuple_windows()
                    .map(|(&u, &v)| graph.edge_weight(&(u, v)).unwrap())
                    .sum();
                println!("correction path: {:?}, weight: {}", correction_path, weight);
            }

            correction_path
                .into_iter()
                .tuple_windows()
                .filter(|(&u, &v)| (u.0 != v.0) || (u.1 != v.1))
                .filter(|(u, v)| !(graph.is_boundary(u).unwrap() && graph.is_boundary(v).unwrap()))
                .for_each(|(&u, &v)| correction_qubit.push(UnGraph::edge_to_qubit((u, v))));
        }
    }
    correction_qubit
}

/// decode
pub fn decode(graph: &UnGraph, m: usize) -> Vec<(i32, i32)> {
    let (local_graph, path_detail) = construct_syndrome_graph(graph, m);
    let matching = minimum_weight_perfect_matching(local_graph);
    decide_correction_qubit(graph, matching, path_detail)
}
