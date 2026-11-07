//! Print Earth's season lengths. `cargo run -p fms-planet --example seasons`

fn main() {
    let lengths = fms_planet::Planet::earth().season_lengths_days();
    println!("Earth northern spring/summer/autumn/winter (days): {lengths:.2?}");
}
