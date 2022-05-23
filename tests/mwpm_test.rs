use knot::decoder::mwpm;
use knot::qubit_graph::ungraph::UnGraph;
use std::cell::Cell;
use std::rc::Rc;

#[test]
fn test_local_dijkstra() {
    let v = vec![
        ((1, 1, 0), (2, 2, 0)),
        ((1, 1, 0), (3, 3, 0)),
        ((2, 2, 0), (4, 4, 0)),
        ((3, 3, 0), (5, 5, 0)),
        ((6, 6, 0), (7, 7, 0)),
        ((3, 3, 0), (8, 8, 0)),
        ((7, 7, 0), (4, 4, 0)),
        ((2, 2, 0), (9, 9, 0)),
        ((10, 10, 0), (1, 1, 0)),
        ((5, 5, 0), (4, 4, 0)),
    ];
    let mut graph = UnGraph::from_edges(&v, 1, 0);
    for (i, edge) in v.iter().enumerate() {
        graph.set_edge_weight(edge, i as f64);
    }
    graph.set_edge_weight(&((5, 5, 0), (4, 4, 0)), 1 as f64);
    for i in 1..11 {
        graph.set_classical_register((i, i, 0), Rc::new(Cell::new(1)));
    }
    let ans = mwpm::local_dijkstra(&graph, 5, &(1, 1, 0));

    assert_eq!(
        ans,
        vec![
            (vec![(2, 2, 0), (1, 1, 0)], 0.0),
            (vec![(3, 3, 0), (1, 1, 0)], 1.0),
            (vec![(4, 4, 0), (2, 2, 0), (1, 1, 0)], 2.0),
            (vec![(5, 5, 0), (4, 4, 0), (2, 2, 0), (1, 1, 0)], 3.0),
            (vec![(8, 8, 0), (3, 3, 0), (1, 1, 0)], 6.0)
        ]
    );
}