use knot::qec_code::rotated_surface_code::RotatedSurfaceCode;
use knot::qubit_network::ErrorDistribution;
use colored::*;

fn main() {
    let loop_num: u64 = 10000;
    let distance = 5;
    let seed = 10;
    let (mean, std_dev) = (0.06, 0.005);
    let distribution = ErrorDistribution::TruncNormal {
        mean,
        std_dev,
        seed,
    };
    // let distribution = ErrorDistribution::Equal(0.01);
    let mut code = RotatedSurfaceCode::new(distance, distance, distribution, 0.01, seed);

    code.initialize();
    code.syndrome_measurement();

    let mut error_num = 0;
    let mut abnormal = 0;

    for i in 0..loop_num {
        /* 
        let distribution = ErrorDistribution::TruncNormal {
            mean,
            std_dev,
            seed: i,
        };
        let mut code = RotatedSurfaceCode::new(distance, distance, distribution, 0.01, seed);

        code.initialize();
        code.syndrome_measurement(); */
        code.reset();
        code.run();
        code.decode_mwpm(distance);

        let ans = code.logical_value();

        if ans == 1 {
            error_num += 1;
        }

        if ans == u8::MAX {
            abnormal += 1;
        }
        println!("ans = {}, loop {}", ans, i);
        if ans != 0 {
            println!("{}", "#########################################################################################\nerror\n#########################################################################################".red());
            break;
        }
        println!("");
    }

    println!("abnormal: {}", abnormal);
    println!("error_num: {}", error_num);
    println!(
        "error rate: {}",
        (error_num + abnormal) as f32 / loop_num as f32
    );
}