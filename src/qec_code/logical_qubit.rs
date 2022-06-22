pub trait LogicalQubit {
    // logical h gate
    fn logical_h(&mut self);

    // logical x gate
    fn logical_x(&mut self);

    // logical z gate
    fn logical_z(&mut self);
}
