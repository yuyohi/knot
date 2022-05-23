pub struct Stabilizer {
    pub ancilla: (i32, i32),
    pauli_product: Vec<Option<(i32, i32)>>,
}

impl Stabilizer {
    pub fn new(ancilla: (i32, i32), pauli_product: Vec<Option<(i32, i32)>>) -> Self {
        Self {
            ancilla,
            pauli_product
        }
    }

    pub fn pauli_product(&self) -> &Vec<Option<(i32, i32)>> {
        &self.pauli_product
    }
}