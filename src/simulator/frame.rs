use ndarray::prelude::*;

pub struct PauliFrame {
    z_frame: Array2<u8>,
    x_frame: Array2<u8>,
}

impl PauliFrame {
    pub fn new_rotated_surface_code(distance: usize) -> Self {
        let z_frame = Array::zeros((distance, distance));
        let x_frame = Array::zeros((distance, distance));

        Self {
            z_frame,
            x_frame,
        }
    }

    /// return z frame
    pub fn z_frame(&self) -> ArrayView2<u8> {
        self.z_frame.view()
    }

    /// return x frame
    pub fn x_frame(&self) -> ArrayView2<u8> {
        self.x_frame.view()
    }

    /// return mutable z frame 
    pub fn z_frame_mut(&mut self) -> ArrayViewMut2<u8> {
        self.z_frame.view_mut()
    }

    /// return mutable x frame 
    pub fn x_frame_mut(&mut self) -> ArrayViewMut2<u8> {
        self.x_frame.view_mut()
    }

    /// reset frame
    pub fn reset(&mut self) {
        let shape = (self.z_frame.shape()[0], self.z_frame.shape()[1]);
        self.z_frame = Array::zeros(shape);
        self.x_frame = Array::zeros(shape);
    }
    
}