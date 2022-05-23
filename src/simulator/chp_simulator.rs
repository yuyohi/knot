use super::{
    core::{Dispatcher, SimulatorCore},
    Operation, SimulatorInterface,
};
use crate::noise::noise_model::NoiseType;
use ndarray::*;
use rand::{rngs::SmallRng, Rng};
use std::cell::Cell;
use std::rc::Rc;

pub struct CHPSimulatorCore {
    qubit_num: usize,
    stabilizer_tableau: Array2<u8>,
    rng: SmallRng,
    classical_register: Vec<u8>,
}

pub struct CHPSimulator {
    core: CHPSimulatorCore,
    dispatcher: Dispatcher,
}

impl CHPSimulator {
    pub fn new(qubit_num: usize, rng: SmallRng) -> Self {
        let size = qubit_num * 2;
        let stabilizer_tableau: Array2<u8> =
            concatenate![Axis(1), Array::eye(size), Array::zeros((size, 1))];

        let operations = Vec::new();
        let classical_register = vec![0; qubit_num];
        let round = 1; // デフォルト値

        CHPSimulator {
            core: CHPSimulatorCore {
                qubit_num,
                stabilizer_tableau,
                rng,
                classical_register,
            },
            dispatcher: Dispatcher::new(operations, round),
        }
    }

    pub fn result(&self) -> &Vec<u8> {
        &self.core.classical_register
    }
}

impl CHPSimulatorCore {
    fn g(&self, x1: u8, z1: u8, x2: u8, z2: u8) -> i8 {
        match (x1, z1, x2, z2) {
            (0, 0, _, _) => 0,
            (1, 1, _, _) => (z2 as i8 - x2 as i8),
            (1, 0, _, _) => z2 as i8 * (2 * x2 as i8 - 1),
            (0, 1, _, _) => x2 as i8 * (1 - 2 * z2 as i8),
            _ => panic!("Tableau parameters must be 0 or 1"),
        }
    }

    fn row_sum(&mut self, h: usize, i: usize) {
        let mut g_sum = 0;

        for j in 0..self.qubit_num {
            g_sum += self.g(
                self.stabilizer_tableau[[i, j]],
                self.stabilizer_tableau[[i, self.qubit_num + j]],
                self.stabilizer_tableau[[h, j]],
                self.stabilizer_tableau[[h, self.qubit_num + j]],
            );
        }

        let checker = 2 * self.stabilizer_tableau[[h, self.qubit_num * 2]] as i8
            + 2 * self.stabilizer_tableau[[i, self.qubit_num * 2]] as i8
            + g_sum;

        match checker % 4 {
            2 | -2 => self.stabilizer_tableau[[h, self.qubit_num * 2]] = 1,
            0 => self.stabilizer_tableau[[h, self.qubit_num * 2]] = 0,
            _ => panic!("Error at row_sum value: {}", checker % 4),
        }

        let (mut row_h, row_i) = self
            .stabilizer_tableau
            .multi_slice_mut((s![h, ..self.qubit_num * 2], s![i, ..self.qubit_num * 2]));

        row_h ^= &row_i;
    }

    fn row_sum_temp(&self, i: usize, temp: &mut Array1<u8>) {
        let mut g_sum = 0;

        for j in 0..self.qubit_num {
            g_sum += self.g(
                self.stabilizer_tableau[[i, j]],
                self.stabilizer_tableau[[i, self.qubit_num + j]],
                temp[j],
                temp[self.qubit_num + j],
            );
        }

        let checker = 2 * temp[self.qubit_num * 2] as i8
            + 2 * self.stabilizer_tableau[[i, self.qubit_num * 2]] as i8
            + g_sum;

        match checker % 4 {
            2 | -2 => temp[self.qubit_num * 2] = 1,
            0 => temp[self.qubit_num * 2] = 0,
            _ => panic!("Error at row_sum"),
        }
    }
}

impl SimulatorCore for CHPSimulatorCore {
    /// CNOT gate
    fn cx(&mut self, a: usize, b: usize) {
        let (mut r, x_a, mut x_b, mut z_a, z_b) = self.stabilizer_tableau.multi_slice_mut((
            s![.., -1],
            s![.., a],
            s![.., b],
            s![.., self.qubit_num + a],
            s![.., self.qubit_num + b],
        ));

        // rを計算
        let buf = &x_a * &z_b * (&x_b ^ &z_a ^ 1);
        r ^= &buf;

        // xを計算
        x_b ^= &x_a;

        // zを計算
        z_a ^= &z_b;
    }

    /// Hadamard gate
    fn h(&mut self, a: usize) {
        let (mut r, x_a, z_a) = self.stabilizer_tableau.multi_slice_mut((
            s![.., -1],
            s![.., a],
            s![.., self.qubit_num + a],
        ));

        let buf = &x_a & &z_a;
        r ^= &buf;

        // XとZを入れ替える
        Zip::from(x_a).and(z_a).for_each(::std::mem::swap);
    }

    /// S gate (Phase gate)
    fn s(&mut self, a: usize) {
        let (mut r, x_a, mut z_a) = self.stabilizer_tableau.multi_slice_mut((
            s![.., -1],
            s![.., a],
            s![.., self.qubit_num + a],
        ));

        let buf = &x_a & &z_a;
        r ^= &buf;

        z_a ^= &x_a;
    }

    ///X gate
    fn x(&mut self, a: usize) {
        let (mut r, z_a) = self
            .stabilizer_tableau
            .multi_slice_mut((s![.., -1], s![.., self.qubit_num + a]));

        r ^= &z_a;
    }

    /// Z gate
    fn z(&mut self, a: usize) {
        let (mut r, x_a) = self
            .stabilizer_tableau
            .multi_slice_mut((s![.., -1], s![.., a]));
        r ^= &x_a;
    }

    /// measurement
    fn measurement(&mut self, a: usize, register: &Rc<Cell<u8>>, error_rate: f64) -> bool {
        let p = self
            .stabilizer_tableau
            .slice(s![self.qubit_num.., a])
            .iter()
            .enumerate()
            .filter(|(_, &x)| x == 1)
            .map(|(i, _)| i + self.qubit_num)
            .collect::<Vec<usize>>();

        // 一つでもXpa = 1のとき、結果はランダム
        if !p.is_empty() {
            let p_destabilizer = self
                .stabilizer_tableau
                .slice(s![..self.qubit_num, a])
                .iter()
                .enumerate()
                .filter(|(_, &x)| x == 1)
                .map(|(i, _)| i)
                .collect::<Vec<usize>>();

            p_destabilizer.iter().for_each(|&i| self.row_sum(i, p[0]));
            p.iter().skip(1).for_each(|&i| self.row_sum(i, p[0]));

            // (p[0] - qubit_num) 行目をp[0]行目に置換
            let (mut q_n, mut q) = self
                .stabilizer_tableau
                .multi_slice_mut((s![p[0] - self.qubit_num, ..], s![p[0], ..]));

            q_n.assign(&q);
            for i in q.iter_mut() {
                *i = 0;
            }
            q[self.qubit_num + a] = 1;

            // rpを1/2でセットし、これが観測結果となる
            if self.rng.gen::<f64>() < 0.5 {
                self.stabilizer_tableau[[p[0], self.qubit_num * 2]] = 1;
            } else {
                self.stabilizer_tableau[[p[0], self.qubit_num * 2]] = 0;
            }

            // 値を格納
            register.set(self.stabilizer_tableau[[p[0], self.qubit_num * 2]]);
        } else {
            // 測定結果が決定的のとき
            let mut temp: Array1<u8> = Array::zeros(self.qubit_num * 2 + 1);

            self.stabilizer_tableau
                .slice(s![..self.qubit_num, a])
                .iter()
                .enumerate()
                .filter(|(_, &i)| i == 1)
                .for_each(|(i, _)| self.row_sum_temp(i + self.qubit_num, &mut temp));

            // 値を格納
            register.set(temp[self.qubit_num * 2]);
        }

        if self.rng.gen::<f64>() < error_rate {
            register.set(register.get() ^ 1);
            if cfg!(debug_assertions) {
                println!("measurement error: {}", a);
            }
            true
        } else {
            false
        }
    }

    /// measurement
    fn measurement_to_zero(&mut self, a: usize) {
        let p = self
            .stabilizer_tableau
            .slice(s![self.qubit_num.., a])
            .iter()
            .enumerate()
            .filter(|(_, &x)| x == 1)
            .map(|(i, _)| i + self.qubit_num)
            .collect::<Vec<usize>>();

        // 一つでもXpa = 1のとき、結果はランダム
        if !p.is_empty() {
            let p_destabilizer = self
                .stabilizer_tableau
                .slice(s![..self.qubit_num, a])
                .iter()
                .enumerate()
                .filter(|(_, &x)| x == 1)
                .map(|(i, _)| i)
                .collect::<Vec<usize>>();

            p_destabilizer.iter().for_each(|&i| self.row_sum(i, p[0]));
            p.iter().skip(1).for_each(|&i| self.row_sum(i, p[0]));

            // (p[0] - qubit_num) 行目をp[0]行目に置換
            let (mut q_n, mut q) = self
                .stabilizer_tableau
                .multi_slice_mut((s![p[0] - self.qubit_num, ..], s![p[0], ..]));

            q_n.assign(&q);
            for i in q.iter_mut() {
                *i = 0;
            }
            q[self.qubit_num + a] = 1;

            // 必ず0にセット(固有値1)
            self.stabilizer_tableau[[p[0], self.qubit_num * 2]] = 0;
        } else {
            // 測定結果が決定的のとき
            let mut temp: Array1<u8> = Array::zeros(self.qubit_num * 2 + 1);

            self.stabilizer_tableau
                .slice(s![..self.qubit_num, a])
                .iter()
                .enumerate()
                .filter(|(_, &i)| i == 1)
                .for_each(|(i, _)| self.row_sum_temp(i + self.qubit_num, &mut temp));
        }
    }

    /// measurement and reset
    fn measurement_and_reset(&mut self, a: usize, register: &Rc<Cell<u8>>, error_rate: f64) {
        let measurement_error = self.measurement(a, register, error_rate);
        if ((register.get() == 1) && !measurement_error)
            || ((register.get() == 0) && measurement_error)
        {
            self.x(a);
        }
    }

    /// insert depolarizing noise
    fn depolarizing(&mut self, a: usize, p: f64) {
        if self.rng.gen::<f64>() < p {
            // insert noise
            match self.rng.gen::<f64>() {
                x if (0.0..1.0 / 3.0).contains(&x) => {
                    // self.z(a);
                    if cfg!(debug_assertions) {
                        println!("z error: {}", a)
                    }
                } // Z error
                x if (1.0 / 3.0..2.0 / 3.0).contains(&x) => {
                    self.x(a);
                    if cfg!(debug_assertions) {
                        println!("x error: {}", a)
                    }
                } // X error
                x if (2.0 / 3.0..1.0).contains(&x) => {
                    // Y error
                    // self.x(a);
                    // self.z(a);
                    if cfg!(debug_assertions) {
                        println!("y error: {}", a)
                    }
                }
                _ => panic!("rng must be 0.0..1.0"),
            }
        }
    }

    fn reset(&mut self) {
        let size = self.qubit_num * 2;
        self.stabilizer_tableau = concatenate![Axis(1), Array::eye(size), Array::zeros((size, 1))];
    }
}

impl SimulatorInterface for CHPSimulator {
    /// add CNOT gate
    fn add_cx(&mut self, a: usize, b: usize) {
        self.dispatcher.push(Operation::CX(a, b));
    }

    /// add Hadamard gate
    fn add_h(&mut self, a: usize) {
        self.dispatcher.push(Operation::H(a));
    }

    /// add S gate (Phase gate)
    fn add_s(&mut self, a: usize) {
        self.dispatcher.push(Operation::S(a));
    }

    /// add X gate
    fn add_x(&mut self, a: usize) {
        self.dispatcher.push(Operation::X(a))
    }

    /// add Z gate
    fn add_z(&mut self, a: usize) {
        self.dispatcher.push(Operation::Z(a))
    }

    /// add measurement
    fn add_measurement(&mut self, a: usize, register: Rc<Cell<u8>>, error_rate: f64) {
        self.dispatcher.push(Operation::M(a, register, error_rate));
    }

    fn add_measurement_to_zero(&mut self, a: usize) {
        self.dispatcher.push(Operation::MToZero(a));
    }

    /// add measurement_and_reset
    fn add_measurement_and_reset(&mut self, a: usize, register: Rc<Cell<u8>>, error_rate: f64) {
        self.dispatcher.push(Operation::MR(a, register, error_rate));
    }

    fn add_noise(&mut self, a: usize, noise_type: NoiseType) {
        match noise_type {
            NoiseType::Depolarizing(p) => self.dispatcher.push(Operation::Depolarizing(a, p)),
        }
    }

    /// Reset stabilizer tableau
    fn reset(&mut self) {
        self.core.reset();
    }

    /// measurement direct
    fn measurement(&mut self, a: usize, register: Rc<Cell<u8>>, error_rate: f64) {
        self.core.measurement(a, &register, error_rate);
    }

    /// run circuit
    fn run(&mut self) {
        let Self { core, dispatcher } = self;

        for op in dispatcher.operations().iter() {
            match op {
                Operation::CX(a, b) => core.cx(*a, *b),
                Operation::H(a) => core.h(*a),
                Operation::Depolarizing(a, p) => core.depolarizing(*a, *p),
                Operation::MR(a, register, error_rate) => {
                    core.measurement_and_reset(*a, register, *error_rate)
                }
                Operation::M(a, register, error_rate) => {
                    core.measurement(*a, register, *error_rate);
                }
                Operation::MToZero(a) => core.measurement_to_zero(*a),
                Operation::S(a) => core.s(*a),
                Operation::X(a) => core.x(*a),
                Operation::Z(a) => core.z(*a),
            }
            //println!("{:?}", op);
        }
    }
}