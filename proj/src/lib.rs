// This is in addition to computational algebra a personal learning experience with rust.
pub mod algebraic_structure;
// pub mod tests;
pub mod integers;
// pub mod integer_computations;
// pub mod prime;
// pub mod prime_statistics;
pub mod random;


// fn options_handler(args: Vec<String>) -> bool {
//     let mut used = false;
//     for arg in &args[1..] {
//         if arg == "generate-test-lists" {
//             prime::generate_primes();
//             prime::generate_non_primes();
//             used = true;
//         } else if arg == "generate-small-primes" {
//             let n = 100000;
//             prime::generate_small_primes(n);
//             used = true;
//         } else if arg == "proj1" {
//             let loops = 10;
//             let naive_square_points: usize = 20;
//             let naive_points: usize = 50;
//             let square_points: usize = 50;
//             tests::plot_timing_naive_square(naive_square_points, loops).expect("Should not fail");
//             tests::plot_timing_naive(naive_points, loops).expect("Should not fail");
//             tests::plot_timing_square(square_points, loops).expect("Should not fail");
//             used = true;
//         } else if arg == "proj2" {
//             continue;
//         } else {
//             println!("Argument {} is not valid.", arg);
//         }
//     }
//     used
// }
//
//
// fn main() {
//     if !options_handler(env::args().collect()) {
//         main_without_arg();
//     }
// }
//
//
// fn main_without_arg() {
//     prime_statistics::time_finding_primes_trial_division_vs_sieving();
// }
