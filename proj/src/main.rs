// This is in addition to computational algebra a personal learning experience with rust.

pub mod algebraic_structure;
pub mod tests;
pub mod integer_computations;


fn main() {
    let loops = 10;
    let naive_square_points: usize = 20;
    let naive_points: usize = 100;
    let square_points: usize = 50;
    tests::plot_timing_naive_square(naive_square_points, loops).expect("Should not fail");
    tests::plot_timing_naive(naive_points, loops).expect("Should not fail");
    tests::plot_timing_square(square_points, loops).expect("Should not fail");
}
