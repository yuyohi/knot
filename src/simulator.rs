use self::chp_simulator::CHPSimulator;
use crate::noise::noise_model::NoiseType;
use std::cell::Cell;
use std::rc::Rc;

pub mod chp_simulator;
pub mod core;
pub mod frame;

/// シミュレータの外部からもアクセスできるオペレーション
pub trait SimulatorInterface {
    /// add CNOT gate
    fn add_cx(&mut self, a: usize, b: usize);

    /// add Hadamard gate
    fn add_h(&mut self, a: usize);

    /// add S gate (Phase gate)
    fn add_s(&mut self, a: usize);

    /// add X gate
    fn add_x(&mut self, a: usize);

    /// add Z gate
    fn add_z(&mut self, a: usize);

    /// add measurement
    fn add_measurement(&mut self, a: usize, register: Rc<Cell<u8>>, error_rate: f64);

    // add measurement as once
    //fn add_measurement_at_once(&mut self, a: Vec<usize>, register: &mut Array3<u8>);

    /// add measurement coercion to zero
    fn add_measurement_to_zero(&mut self, a: usize);

    /// add measurement_and_reset
    fn add_measurement_and_reset(&mut self, a: usize, register: Rc<Cell<u8>>, error_rate: f64);

    /// add noise
    fn add_noise(&mut self, a: usize, noise_type: NoiseType);

    /// add Reset stabilizer tableau
    fn reset(&mut self);

    /// add run circuit
    fn run(&mut self);

    /// direct measurement
    fn measurement(&mut self, a: usize, register: Rc<Cell<u8>>, error_rate: f64);

}

pub enum SimulatorWrapper {
    CHPSimulator(CHPSimulator),
}

impl SimulatorInterface for SimulatorWrapper {
    fn add_cx(&mut self, a: usize, b: usize) {
        match *self {
            SimulatorWrapper::CHPSimulator(ref mut sim) => sim.add_cx(a, b),
        };
    }

    fn add_h(&mut self, a: usize) {
        match *self {
            SimulatorWrapper::CHPSimulator(ref mut sim) => sim.add_h(a),
        };
    }

    fn add_s(&mut self, a: usize) {
        match *self {
            SimulatorWrapper::CHPSimulator(ref mut sim) => sim.add_s(a),
        };
    }

    fn add_x(&mut self, a: usize) {
        match *self {
            SimulatorWrapper::CHPSimulator(ref mut sim) => sim.add_x(a),
        };
    }

    fn add_z(&mut self, a: usize) {
        match *self {
            SimulatorWrapper::CHPSimulator(ref mut sim) => sim.add_z(a),
        };
    }

    fn add_measurement(&mut self, a: usize, register: Rc<Cell<u8>>, error_rate: f64) {
        match *self {
            SimulatorWrapper::CHPSimulator(ref mut sim) => sim.add_measurement(a, register, error_rate),
        }
    }

    fn add_measurement_to_zero(&mut self, a: usize) {
        match *self {
            SimulatorWrapper::CHPSimulator(ref mut sim) => sim.add_measurement_to_zero(a),
        }
    }

    fn add_measurement_and_reset(&mut self, a: usize, register: Rc<Cell<u8>>, error_rate: f64) {
        match *self {
            SimulatorWrapper::CHPSimulator(ref mut sim) => sim.add_measurement_and_reset(a, register, error_rate),
        }
    }

    fn measurement(&mut self, a: usize, register: Rc<Cell<u8>>, error_rate: f64) {
        match *self {
            SimulatorWrapper::CHPSimulator(ref mut sim) => sim.measurement(a, register, error_rate),
        }
    }

    fn add_noise(&mut self, a: usize, noise_type: NoiseType) {
        match *self {
            SimulatorWrapper::CHPSimulator(ref mut sim) => sim.add_noise(a, noise_type),
        }
    }

    fn reset(&mut self) {
        match *self {
            SimulatorWrapper::CHPSimulator(ref mut sim) => sim.reset(),
        };
    }

    fn run(&mut self) {
        match *self {
            SimulatorWrapper::CHPSimulator(ref mut sim) => sim.run(),
        }
    }
}

#[derive(Debug)]
pub enum Operation {
    CX(usize, usize),
    H(usize),
    S(usize),
    X(usize),
    Z(usize),
    M(usize, Rc<Cell<u8>>, f64),
    MR(usize, Rc<Cell<u8>>, f64),
    MToZero(usize),
    Depolarizing(usize, f64),
    //MAll(char)
}

pub enum Type {
    CHPSimulator,
}
