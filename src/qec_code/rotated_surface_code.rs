use itertools::Itertools;
use std::cell::Cell;
use std::collections::HashMap;
use std::rc::Rc;

use crate::decoder::mwpm;
use crate::noise::noise_model::NoiseType;
use crate::qec_code::stabilizer::Stabilizer;
use crate::qubit_graph::ungraph::UnGraph;
use crate::qubit_network::{ErrorDistribution, QubitNetwork};
use crate::simulator::{frame::PauliFrame, Type};

pub struct RotatedSurfaceCode {
    distance: usize,
    round: usize,
    network: QubitNetwork,
    z_stabilizers: Vec<Stabilizer>,
    x_stabilizers: Vec<Stabilizer>,
    data_qubit: Vec<(i32, i32)>,
    classical_register: Vec<Vec<Rc<Cell<u8>>>>,
    measurement_graph_z: UnGraph,
    measurement_graph_x: UnGraph,
    single_round_measurement_graph_z: UnGraph,
    pauli_frame: PauliFrame,
    error_rate_mean: f64,
    measurement_error_rate: f64,
}

impl RotatedSurfaceCode {
    pub fn new(
        distance: usize,
        round: usize,
        qubit_distribution: ErrorDistribution,
        p_m: f64,
        seed: u64,
    ) -> Self {
        if distance % 2 == 0 {
            panic!("distance must be odd number.");
        }

        // measurement bitの生成
        let (mut measurement_qubit_z, mut measurement_qubit_x) =
            Self::gen_measurement_qubit(distance);
        // 測定ビットを効率的に測定できるように並びかえ
        measurement_qubit_x.sort_by(|l, r| match l.1.cmp(&r.1) {
            std::cmp::Ordering::Equal => l.0.cmp(&r.0),
            other => other,
        });
        measurement_qubit_z.sort_by(|l, r| match l.0.cmp(&r.0) {
            std::cmp::Ordering::Equal => l.1.cmp(&r.1),
            other => other,
        });

        // data bitの生成
        let mut data_qubit = vec![];
        for x in (0..distance as i32 * 2).step_by(2) {
            for y in (0..distance as i32 * 2).step_by(2) {
                data_qubit.push((x, y));
            }
        }

        // qubit networkの作成
        let p_mean = qubit_distribution.mean();
        let network = QubitNetwork::new_rotated_planer_lattice_from_vec(
            data_qubit.clone(),
            measurement_qubit_z.clone(),
            measurement_qubit_x.clone(),
            qubit_distribution,
            p_mean,
            Type::CHPSimulator,
            seed,
        );

        // make syndrome graph
        let measurement_graph_z = Self::gen_measurement_graph(
            &network,
            &measurement_qubit_z,
            round,
            distance,
            'Z',
            p_mean,
            seed,
        );
        let measurement_graph_x = Self::gen_measurement_graph(
            &network,
            &measurement_qubit_x,
            round,
            distance,
            'X',
            p_mean,
            seed,
        );
        let single_round_measurement_graph_z = Self::gen_measurement_graph(
            &network,
            &measurement_qubit_z,
            1,
            distance,
            'Z',
            p_mean,
            seed,
        );

        // make pauli frame
        let pauli_frame = PauliFrame::new_rotated_surface_code(distance);

        // make stabilizers
        let z_stabilizers = Self::gen_stabilizer(&measurement_qubit_z, &data_qubit, 'Z');
        let x_stabilizers = Self::gen_stabilizer(&measurement_qubit_x, &data_qubit, 'X');

        // data qubit の測定結果を格納する行列
        let classical_register = (0..distance)
            .map(|_| (0..distance).map(|_| Rc::new(Cell::new(0))).collect())
            .collect();

        Self {
            distance,
            round,
            network,
            z_stabilizers,
            x_stabilizers,
            data_qubit,
            classical_register,
            measurement_graph_z,
            measurement_graph_x,
            single_round_measurement_graph_z,
            pauli_frame,
            error_rate_mean: p_mean,
            measurement_error_rate: p_m,
        }
    }

    /// generate stabilizer
    fn gen_stabilizer(
        measurement_qubit: &[(i32, i32)],
        data_qubit: &[(i32, i32)],
        mode: char,
    ) -> Vec<Stabilizer> {
        let order = match mode {
            'X' => [(1, 1), (1, -1), (-1, 1), (-1, -1)],
            'Z' => [(1, 1), (-1, 1), (1, -1), (-1, -1)],
            _ => panic!("mode must be X or Z"),
        };

        let mut stabilizers = Vec::new();

        // 各measurement bitについて対応するstabilizerを作成
        for ancilla in measurement_qubit.iter() {
            let mut pauli_product = Vec::new();
            for d in order.iter() {
                match (ancilla.0 + d.0, ancilla.1 + d.1) {
                    v if data_qubit.contains(&v) => {
                        pauli_product.push(Some(v));
                    }
                    _ => pauli_product.push(None),
                }
            }
            stabilizers.push(Stabilizer::new(*ancilla, pauli_product));
        }
        stabilizers
    }

    /// generate measurement qubits
    fn gen_measurement_qubit(distance: usize) -> (Vec<(i32, i32)>, Vec<(i32, i32)>) {
        let mut measurement_qubit_temp = Vec::new();

        for y in (1..distance as i32 * 2 - 1).step_by(4) {
            for x in (1..distance as i32 * 2).step_by(2) {
                measurement_qubit_temp.push((x, y));
            }
        }
        for y in (3..distance as i32 * 2 - 1).step_by(4) {
            for x in (-1..distance as i32 * 2 - 1).step_by(2) {
                measurement_qubit_temp.push((x, y));
            }
        }

        // 上下のqubitを追加
        for x in (1..distance as i32 * 2 - 2).step_by(4) {
            measurement_qubit_temp.push((x, -1));
        }
        for x in (3..distance as i32 * 2 - 2).step_by(4) {
            measurement_qubit_temp.push((x, distance as i32 * 2 - 1));
        }

        // Z, Xに振り分ける
        let mut measurement_qubit_z = Vec::new();
        let mut measurement_qubit_x = Vec::new();

        let mut qubit_is_z = true;
        for y in (-1..distance as i32 * 2).step_by(2) {
            for x in (-1..distance as i32 * 2).step_by(2) {
                if measurement_qubit_temp.contains(&(x, y)) {
                    match qubit_is_z {
                        true => measurement_qubit_z.push((x, y)),
                        false => measurement_qubit_x.push((x, y)),
                    }
                }
                // 交互に入れ替える
                qubit_is_z = !qubit_is_z;
            }
            qubit_is_z = !qubit_is_z;
        }

        (measurement_qubit_z, measurement_qubit_x)
    }

    /// generate measurement graph
    fn gen_measurement_graph(
        qubit_network: &QubitNetwork,
        measurement_qubit: &[(i32, i32)],
        round: usize,
        distance: usize,
        mode: char,
        p_mean: f64,
        seed: u64,
    ) -> UnGraph {
        let mut graph = if round == 1 {
            UnGraph::new(0, seed)
        } else {
            UnGraph::new(round, seed)
        };
        let direction = [(2, 2), (-2, 2)];
        let boundary_direction = [(2, 2), (-2, 2), (-2, -2), (2, -2)];

        let mut edges = Vec::new();

        // 通常のedgeを追加
        for &u in measurement_qubit.iter() {
            for d in direction.iter() {
                match (u.0 + d.0, u.1 + d.1) {
                    v if measurement_qubit.contains(&v) => {
                        edges.push((u, v));
                    }
                    _ => (),
                }
            }
        }

        // boundary nodeを追加
        let boundary_num = (distance / 2 + 1) as i32;
        let mut boundary_node = Vec::new();

        match mode {
            'Z' => {
                let x_start_list = [-1, 1];
                let height = (distance - 1) * 2;
                let y_list = [-1, height as i32 + 1];
                for (y, x_start) in y_list.iter().zip(x_start_list.iter()) {
                    for i in 0..boundary_num {
                        let x = x_start + 4 * i;
                        boundary_node.push((x, *y))
                    }
                }
            }
            'X' => {
                let y_start_list = [1, -1];
                let width = (distance - 1) * 2;
                let x_references = [-1, width as i32 + 1];
                for (x, y_start) in x_references.iter().zip(y_start_list.iter()) {
                    for i in 0..boundary_num {
                        let y = y_start + 4 * i;
                        boundary_node.push((*x, y))
                    }
                }
            }
            _ => panic!("Invalid mode"),
        }

        // boundary nodeと普通のnodeを結ぶ
        for &u in boundary_node.iter() {
            for d in boundary_direction.iter() {
                match (u.0 + d.0, u.1 + d.1) {
                    v if measurement_qubit.contains(&v) => {
                        edges.push((u, v));
                    }
                    _ => (),
                }
            }
        }

        for t in 0..round as i32 {
            // time boundary
            if round != 1 {
                // 次のroundの同座標に対するedgeを追加
                let mut time_edge = Vec::new();
                for &(x, y) in measurement_qubit.iter() {
                    time_edge.push(((x, y, t), (x, y, t + 1)));
                }
                // 最後のroundでは、次の時間はboundaryなので、boundary同士をweight0のedgeで繋ぐ
                if t == (round as i32 - 1) {
                    let time_boundary_edge: Vec<_> = measurement_qubit
                        .iter()
                        .tuple_windows()
                        .map(|(&(u_x, u_y), &(v_x, v_y))| ((u_x, u_y, t + 1), (v_x, v_y, t + 1)))
                        .collect();
                    graph.add_edges_from(&time_boundary_edge);
                    graph.set_edges_weight(&time_boundary_edge, 0.0);

                    // 時間のboundaryと空間のboundaryを繋ぐ
                    let t_boundary_to_boundary = (
                        (measurement_qubit[0].0, measurement_qubit[0].1, t + 1),
                        (boundary_node[0].0, boundary_node[0].1, t),
                    );
                    graph.add_edge_from(&t_boundary_to_boundary);
                    graph.set_edge_weight(&t_boundary_to_boundary, p_mean / 10.0);
                }

                graph.add_edges_from(&time_edge);
                graph.set_edges_weight(&time_edge, p_mean);
            }

            for &((u_x, u_y), (v_x, v_y)) in edges.iter() {
                graph.add_edge_from(&((u_x, u_y, t), (v_x, v_y, t)));
                let qubit_coord = ((u_x + v_x) / 2, (u_y + v_y) / 2);
                graph.set_edge_weight(
                    &((u_x, u_y, t), (v_x, v_y, t)),
                    qubit_network.qubit_error_rate(&qubit_coord),
                );
            }

            // 最後のround以外は次のboundary nodeとも繋ぐ(weightは0ではない)
            if t != (round as i32 - 1) {
                for &(x, y) in boundary_node.iter() {
                    let edge = ((x, y, t), (x, y, t + 1));
                    graph.add_edge_from(&edge);
                    graph.set_edge_weight(&edge, p_mean / 10.0);
                }
            }

            // boundary node 同士をweight0のedgeで繋ぐ
            let boundary_edge: Vec<_> = boundary_node
                .iter()
                .tuple_windows()
                .map(|(&u, &v)| ((u.0, u.1, t), (v.0, v.1, t)))
                .collect();
            graph.add_edges_from(&boundary_edge);
            graph.set_edges_weight(&boundary_edge, 0.0);

            // boundaryかどうかとregisterを設定
            for &(x, y) in measurement_qubit.iter() {
                graph.set_is_boundary((x, y, t), false);
                graph.set_classical_register((x, y, t), Rc::new(Cell::new(0))); // 順番が大事
                if (t == (round as i32 - 1)) && (round != 1) {
                    // 最後のroundでは、時間方向のboundaryを設定
                    graph.set_is_boundary((x, y, t + 1), true);
                    graph.set_classical_register((x, y, t + 1), Rc::new(Cell::new(0)));
                }
            }
            for &(x, y) in boundary_node.iter() {
                graph.set_is_boundary((x, y, t), true);
                graph.set_classical_register((x, y, t), Rc::new(Cell::new(0)));
            }
        }

        graph
    }

    /// syndrome measurement
    pub fn syndrome_measurement(&mut self) {
        let Self {
            round,
            network,
            z_stabilizers,
            x_stabilizers,
            measurement_graph_z,
            measurement_graph_x,
            data_qubit,
            ..
        } = self;

        for t in 0..*round as i32 {
            // 仮 現象論的ノイズ
            for c in data_qubit.iter() {
                let noise_type = NoiseType::Depolarizing(network.qubit_error_rate(c));
                network.insert_noise(*c, noise_type);
            }

            // XスタビライザーにHゲートを作用させる
            for Stabilizer { ancilla, .. } in x_stabilizers.iter() {
                network.h(*ancilla);
                // network.insert_noise(*ancilla, noise_type); // circuit noise
            }

            // CNOT
            for i in 0..4 {
                for (x_stab, z_stab) in x_stabilizers.iter().zip(z_stabilizers.iter()) {
                    // data bitが存在するときのみCNOT
                    match z_stab.pauli_product().get(i).unwrap() {
                        Some(data_coord) => {
                            network.cx(*data_coord, z_stab.ancilla);
                            // network.insert_noise(z_stab.ancilla, noise_type); // circuit noise
                            // network.insert_noise(*data_coord, noise_type);
                        }
                        None => (),
                    }
                    match x_stab.pauli_product().get(i).unwrap() {
                        Some(data_coord) => {
                            network.cx(x_stab.ancilla, *data_coord);
                            // network.insert_noise(*data_coord, noise_type); // circuit noise
                            // network.insert_noise(x_stab.ancilla, noise_type);
                        }
                        None => (),
                    }
                }
            }

            // XスタビライザーにHゲートを作用させる
            for Stabilizer { ancilla, .. } in x_stabilizers.iter() {
                network.h(*ancilla);
                // network.insert_noise(*ancilla, noise_type); // circuit noise
            }

            // measurement qubitの測定
            // Z
            for Stabilizer { ancilla, .. } in z_stabilizers.iter() {
                network.measurement_and_reset(
                    *ancilla,
                    Rc::clone(
                        measurement_graph_z
                            .get_register(&(ancilla.0, ancilla.1, t))
                            .unwrap(),
                    ),
                    self.measurement_error_rate,
                );
            }
            // X
            for Stabilizer { ancilla, .. } in x_stabilizers.iter() {
                network.measurement_and_reset(
                    *ancilla,
                    Rc::clone(
                        measurement_graph_x
                            .get_register(&(ancilla.0, ancilla.1, t))
                            .unwrap(),
                    ),
                    self.measurement_error_rate,
                );
            }
        }
    }

    /// encoding logical one
    pub fn initialize(&mut self) {
        let Self {
            network,
            x_stabilizers,
            ..
        } = self;

        // XスタビライザーにHゲートを作用させる
        for Stabilizer { ancilla, .. } in x_stabilizers.iter() {
            network.h(*ancilla);
        }
        // CNOT
        for x_stab in x_stabilizers.iter() {
            for x_data_coord in x_stab.pauli_product().iter() {
                // data bitが存在するときのみCNOT
                match x_data_coord {
                    Some(data_coord) => {
                        network.cx(x_stab.ancilla, *data_coord);
                    }
                    None => (),
                }
            }
        }

        // XスタビライザーにHゲートを作用させる
        for Stabilizer { ancilla, .. } in x_stabilizers.iter() {
            network.h(*ancilla);
        }

        // ancilla qubit の測定 (強制的に固有値+1に射影する)
        for Stabilizer { ancilla, .. } in x_stabilizers.iter() {
            network.measurement_to_zero(*ancilla);
        }
    }

    /// logical z measurement
    fn logical_measurement(&mut self) {
        let Self {
            network,
            data_qubit,
            classical_register,
            ..
        } = self;

        for &(x, y) in data_qubit.iter() {
            debug_assert!(x >= 0, "data coord must not be negative number");
            debug_assert!(y >= 0, "data coord must not be negative number");

            network.measurement_direct(
                (x, y),
                Rc::clone(&classical_register[(x / 2) as usize][(y / 2) as usize]),
                self.measurement_error_rate,
            );
        }
    }

    /// correct z error
    fn correct_z_error(&mut self) {
        let Self {
            classical_register,
            pauli_frame,
            ..
        } = self;

        pauli_frame
            .x_frame_mut()
            .iter()
            .zip(classical_register.iter().flatten())
            .for_each(|(&frame, register)| register.set(register.get() ^ frame));
    }

    /// decode logical value
    fn decode_logical_value(&mut self) {
        // single_round_measurement_graph_xにparityの情報を書き込む
        for z_stab in self.z_stabilizers.iter() {
            let parity =
                z_stab
                    .pauli_product()
                    .iter()
                    .filter_map(|n| *n)
                    .fold(0, |parity, (x, y)| {
                        parity ^ self.classical_register[x as usize / 2][y as usize / 2].get()
                    });

            let ancilla = z_stab.ancilla;
            self.single_round_measurement_graph_z
                .get_register(&(ancilla.0, ancilla.1, 0))
                .unwrap()
                .set(parity);
        }

        if cfg!(debug_assertions) {
            self.single_round_measurement_graph_z.show_all_defect();
        }

        Self::flip_defect(1, &mut self.single_round_measurement_graph_z);
        let correction_qubit_x = mwpm::decode(&self.single_round_measurement_graph_z, 10);
        if cfg!(debug_assertions) {
            println!("correction_qubit_x {:?}", correction_qubit_x);
        }
        // correction
        for (x, y) in correction_qubit_x.into_iter() {
            let register = &self.classical_register[x as usize / 2][y as usize / 2];
            register.set(register.get() ^ 1);
        }
    }

    /// return logical value
    pub fn logical_value(&mut self) -> u8 {
        self.logical_measurement();
        self.correct_z_error();
        if cfg!(debug_assertions) {
            println!("start logical decode");
        }
        self.decode_logical_value();

        let result = self.classical_register();

        // Rc<Cell<u8>>をほどいて縦方向に足す
        let result_vec = result
            .iter()
            .map(|row| row.iter().map(|value| value.get()).collect::<Vec<u8>>())
            .reduce(|row_a, row_b| {
                row_a
                    .iter()
                    .zip(row_b.iter())
                    .map(|(&a, &b)| a + b)
                    .collect()
            })
            .unwrap();

        let logical_value: usize = result_vec.into_iter().map(|v| (v % 2) as usize).sum();

        match logical_value {
            n if n == self.distance => 1,
            0 => 0,
            _ => u8::MAX,
        }
    }

    /// decode by mwpm
    pub fn decode_mwpm(&mut self, m: usize) {
        if cfg!(debug_assertions) {
            print!("before xor z: ");
            self.measurement_graph_z.show_all_defect();
            print!("before xor x: ");
            self.measurement_graph_x.show_all_defect();
        }

        self.measurement_graph_z.xor_to_last_time();
        self.measurement_graph_x.xor_to_last_time();

        if cfg!(debug_assertions) {
            print!("after xor z: ");
            self.measurement_graph_z.show_all_defect();
            print!("after xor x: ");
            self.measurement_graph_x.show_all_defect();
        }

        Self::flip_defect(self.round, &mut self.measurement_graph_z);
        Self::flip_defect(self.round, &mut self.measurement_graph_x);

        if cfg!(debug_assertions) {
            print!("after flip z: ");
            self.measurement_graph_z.show_all_defect();
            print!("after flip x: ");
            self.measurement_graph_x.show_all_defect();
        }

        let correction_qubit_z = mwpm::decode(&self.measurement_graph_x, m);
        let correction_qubit_x = mwpm::decode(&self.measurement_graph_z, m);

        if cfg!(debug_assertions) {
            println!("correction_qubit_x {:?}", correction_qubit_x);
            println!("correction_qubit_z {:?}", correction_qubit_z);
        }

        // pauli frameに設定
        let mut z_frame = self.pauli_frame.z_frame_mut();
        for (x, y) in correction_qubit_z.into_iter() {
            debug_assert!(
                x % 2 == 0,
                "data coord must not be odd number x: {:?}",
                (x, y)
            );
            debug_assert!(
                y % 2 == 0,
                "data coord must not be odd number y: {:?}",
                (x, y)
            );
            z_frame[((x / 2) as usize, (y / 2) as usize)] ^= 1;
        }
        let mut x_frame = self.pauli_frame.x_frame_mut();
        for (x, y) in correction_qubit_x.into_iter() {
            debug_assert!(
                x % 2 == 0,
                "data coord must not be odd number: {:?}",
                (x, y)
            );
            debug_assert!(
                y % 2 == 0,
                "data coord must not be odd number: {:?}",
                (x, y)
            );
            x_frame[((x / 2) as usize, (y / 2) as usize)] ^= 1;
        }
    }

    /// 必要に応じてboundaryを反転させる
    fn flip_defect(round: usize, measurement_graph: &mut UnGraph) {
        for t in 0..round {
            let defect_num = measurement_graph.classical_register()[t]
                .iter()
                .filter(|defect| defect.get() == 1)
                .count();
            let index_to_coord = measurement_graph.index_to_coord();
            if defect_num % 2 == 1 {
                for (coord, defect) in measurement_graph.classical_register()[t]
                    .iter()
                    .enumerate()
                    .map(|(index, defect)| (index_to_coord.get(&(index, t)).unwrap(), defect))
                    .filter(|&(coord, _)| measurement_graph.is_boundary(coord).unwrap())
                {
                    if defect.get() == 0 {
                        measurement_graph.flip_classical_register(coord, 1);
                        break;
                    }
                }
            }
        }
    }

    /// run circuit
    pub fn run(&mut self) {
        self.network.run();
    }

    /// reset code
    pub fn reset(&mut self) {
        self.measurement_graph_x.reset_register();
        self.measurement_graph_z.reset_register();
        self.single_round_measurement_graph_z.reset_register();
        self.pauli_frame.reset();
        self.network.reset();
    }

    pub fn classical_register(&self) -> &Vec<Vec<Rc<Cell<u8>>>> {
        &self.classical_register
    }

    pub fn index_to_sim(&self) -> &HashMap<(i32, i32), usize> {
        self.network.index_to_sim()
    }
}
