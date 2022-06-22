use rand::distributions::Distribution;
use rand::{rngs::SmallRng, SeedableRng, Rng};
use statrs::distribution::Normal;
use std::cell::Cell;
use std::collections::HashMap;
use std::rc::Rc;

use crate::noise::noise_model::NoiseType;
use crate::simulator::{self, SimulatorInterface, SimulatorWrapper, Type};

pub struct QubitNetwork {
    network: HashMap<(i32, i32), Vec<(i32, i32)>>,
    bit_error_map: HashMap<(i32, i32), f64>,
    connection_error_map: HashMap<((i32, i32), (i32, i32)), f64>,
    index_to_sim: HashMap<(i32, i32), usize>,
    sim: SimulatorWrapper,
    error_rate: f64,
}

impl QubitNetwork {
    /// 縦横の大きさからrotated surface codeに適したlatticeを作成する
    pub fn new_rotated_planer_lattice_from_rectangle(
        vertical: usize,
        horizontal: usize,
        p: f64,
        qubit_distribution: ErrorDistribution,
        sim_type: simulator::Type,
        seed: u64,
    ) -> Self {
        // node一覧を作成
        // data_qubitを追加
        let mut qubit_index: Vec<(i32, i32)> = vec![];
        for x in (0..horizontal as i32 * 2).step_by(2) {
            for y in (0..vertical as i32 * 2).step_by(2) {
                qubit_index.push((x, y));
            }
        }
        // measurement qubitを追加 (dual lattice)
        for x in (-1..horizontal as i32 * 2 + 1).step_by(2) {
            for y in (-1..vertical as i32 * 2 + 1).step_by(2) {
                qubit_index.push((x, y));
            }
        }

        Self::new_rotated_planer_lattice(qubit_index, qubit_distribution, p, sim_type, seed)
    }

    /// 受け取ったvecを元にrotated surface codeに適したlatticeを作成する
    pub fn new_rotated_planer_lattice_from_vec(
        data_qubit: Vec<(i32, i32)>,
        measurement_qubit_z: Vec<(i32, i32)>,
        measurement_qubit_x: Vec<(i32, i32)>,
        qubit_distribution: ErrorDistribution,
        p: f64,
        sim_type: simulator::Type,
        seed: u64,
    ) -> Self {
        let mut qubit_index = Vec::new();
        qubit_index.extend(data_qubit);
        qubit_index.extend(measurement_qubit_z);
        qubit_index.extend(measurement_qubit_x);

        Self::new_rotated_planer_lattice(qubit_index, qubit_distribution, p, sim_type, seed)
    }

    /// gen rotated_surface_lattice
    fn new_rotated_planer_lattice(
        qubit_index: Vec<(i32, i32)>,
        qubit_distribution: ErrorDistribution,
        p: f64,
        sim_type: simulator::Type,
        seed: u64,
    ) -> Self {
        let mut network = HashMap::new();

        // 斜めにedgeを追加
        let direction = [(1, 1), (-1, 1), (-1, -1), (1, -1)];
        for &u in qubit_index.iter() {
            for &d in direction.iter() {
                match (u.0 + d.0, u.1 + d.1) {
                    v if qubit_index.contains(&v) => {
                        network.entry(u).or_insert_with(Vec::new).push(v);
                    }
                    _ => (),
                }
            }
        }

        // qubitのerror rate dictを作成
        let mut bit_error_map = HashMap::new();
        let mut generator = qubit_distribution.generator();
        for &qubit in qubit_index.iter() {
            bit_error_map.insert(qubit, generator.gen());
        }
        // connectionのerror rate dictを作成
        let mut connection_error_map = HashMap::new();
        for &u in qubit_index.iter() {
            for &v in qubit_index.iter() {
                connection_error_map.insert((u, v), p);
            }
        }

        // simulatorを生成
        let rng_sim = SmallRng::seed_from_u64(seed + 1);

        let sim = match sim_type {
            Type::CHPSimulator => simulator::SimulatorWrapper::CHPSimulator(
                simulator::chp_simulator::CHPSimulator::new(qubit_index.len(), rng_sim),
            ),
        };
        // (i32, i32)のindexからsimulatorのindexのusizeへのmap
        let mut index_to_sim = HashMap::new();
        for (i, &coord) in qubit_index.iter().enumerate() {
            index_to_sim.insert(coord, i);
        }
        QubitNetwork {
            network,
            bit_error_map,
            connection_error_map,
            index_to_sim,
            sim,
            error_rate: p,
        }
    }

    /// 指定したqubitのerror rateを返す
    pub fn qubit_error_rate(&self, index: &(i32, i32)) -> f64 {
        let p = self
            .bit_error_map
            .get(index)
            .unwrap_or_else(|| panic!("index does not exist: {:?}", index));
        *p
    }

    /// 指定したconnectionのerror rateを返す
    pub fn connection_error_rate(&self, index: &((i32, i32), (i32, i32))) -> f64 {
        let p = self
            .connection_error_map
            .get(index)
            .expect("index does not exist");
        *p
    }

    /// ゲート操作を追加する
    /// CNOT gate
    pub fn cx(&mut self, a: (i32, i32), b: (i32, i32)) {
        debug_assert!(self.connection_error_map.contains_key(&(a, b)));
        self.sim.add_cx(
            *self.index_to_sim.get(&a).expect("index does not exist"),
            *self.index_to_sim.get(&b).expect("index does not exist"),
        );
    }

    /// H gate
    pub fn h(&mut self, a: (i32, i32)) {
        self.sim
            .add_h(*self.index_to_sim.get(&a).expect("index does not exist"));
    }

    /// S gate
    pub fn s(&mut self, a: (i32, i32)) {
        self.sim
            .add_s(*self.index_to_sim.get(&a).expect("index does not exist"));
    }

    /// x gate
    pub fn x(&mut self, a: (i32, i32)) {
        self.sim
            .add_x(*self.index_to_sim.get(&a).expect("index does not exist"));
    }

    /// z gate
    pub fn z(&mut self, a: (i32, i32)) {
        self.sim
            .add_z(*self.index_to_sim.get(&a).expect("index does not exist"));
    }

    /// measurement
    pub fn measurement(&mut self, a: (i32, i32), register: Rc<Cell<u8>>, error_rate: f64) {
        self.sim.add_measurement(
            *self.index_to_sim.get(&a).expect("index does not exist"),
            register,
            error_rate,
        );
    }

    ///measurement direct
    pub fn measurement_direct(&mut self, a: (i32, i32), register: Rc<Cell<u8>>, error_rate: f64) {
        self.sim.measurement(
            *self.index_to_sim.get(&a).expect("index does not exist"),
            register,
            error_rate,
        );
    }

    /// measurement to zero
    pub fn measurement_to_zero(&mut self, a: (i32, i32)) {
        self.sim
            .add_measurement_to_zero(*self.index_to_sim.get(&a).expect("index does not exist"));
    }

    /// measurement and reset
    pub fn measurement_and_reset(
        &mut self,
        a: (i32, i32),
        register: Rc<Cell<u8>>,
        error_rate: f64,
    ) {
        self.sim.add_measurement_and_reset(
            *self.index_to_sim.get(&a).expect("index does not exist"),
            register,
            error_rate,
        );
    }

    pub fn insert_noise(&mut self, a: (i32, i32), noise_type: NoiseType) {
        self.sim.add_noise(
            *self.index_to_sim.get(&a).expect("index does not exist"),
            noise_type,
        )
    }

    /// 指定された座標がネットワークに存在するかを判定する
    pub fn check_contains(&self, a: (i32, i32)) -> bool {
        self.network.contains_key(&a)
    }

    /// 回路を実行する
    pub fn run(&mut self) {
        self.sim.run();
    }

    /// get index_to_sim
    pub fn index_to_sim(&self) -> &HashMap<(i32, i32), usize> {
        &self.index_to_sim
    }

    /// return error rate
    pub fn error_rate(&self) -> f64 {
        self.error_rate
    }

    /// reset tableau
    pub fn reset(&mut self) {
        self.sim.reset();
    }
}

pub enum ErrorDistribution {
    Equal(f64),
    TruncNormal {
        mean: f64,
        std_dev: f64,
        seed: u64,
    },
    InsertingDefect {
        error_rate: f64,
        inserting_defect_rate: f64,
        defect_rate: f64,
        seed: u64,
    },
}

impl ErrorDistribution {
    pub fn mean(&self) -> f64 {
        match self {
            Self::Equal(p) => *p,
            Self::TruncNormal { mean, .. } => *mean,
            Self::InsertingDefect { error_rate, .. } => *error_rate,
        }
    }

    /// return rng generator
    pub fn generator(&self) -> Generator {
        match self {
            Self::Equal(p) => Generator::Equal(*p),
            Self::TruncNormal {
                mean,
                std_dev,
                seed,
            } => {
                let distribution = Normal::new(*mean, *std_dev).unwrap();
                let rng = SmallRng::seed_from_u64(*seed);
                Generator::TruncNormal { distribution, rng }
            }
            Self::InsertingDefect {
                error_rate,
                inserting_defect_rate,
                defect_rate,
                seed,
            } => {
                let rng = SmallRng::seed_from_u64(*seed);
                Generator::InsertingDefect {
                    error_rate: *error_rate,
                    inserting_defect_rate: *inserting_defect_rate,
                    defect_rate: *defect_rate,
                    rng,
                }
            }
        }
    }
}

pub enum Generator {
    Equal(f64),
    TruncNormal {
        distribution: Normal,
        rng: SmallRng,
    },
    InsertingDefect {
        error_rate: f64,
        inserting_defect_rate: f64,
        defect_rate: f64,
        rng: SmallRng,
    },
}

impl Generator {
    pub fn gen(&mut self) -> f64 {
        match self {
            Self::Equal(p) => *p,
            Self::TruncNormal { distribution, rng } => {
                let mut p = distribution.sample(rng);
                while !(0.0..1.0).contains(&p) {
                    // pが0以上1未満ではないといきは再生成
                    p = distribution.sample(rng);
                }
                p
            }
            Self::InsertingDefect {
                error_rate,
                inserting_defect_rate,
                defect_rate,
                rng,
            } => {
                if rng.gen::<f64>() < *inserting_defect_rate {
                    *defect_rate
                } else {
                    *error_rate
                }
            }
        }
    }
}
