use knot::qec_code::rotated_surface_code::RotatedSurfaceCode;
use knot::qubit_network::ErrorDistribution;
use indicatif::ProgressBar;

fn main() {
    let loop_num = 10000;
    let distance = [5, 7];
    let error_rate = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08];
    let std_dev = 0.005;
    let seed = 1;

    let loop_time = distance.len() * error_rate.len();
    let mut count = 1;

    let mut result = vec![Vec::new(); distance.len()];

    for (&d, r) in distance.iter().zip(result.iter_mut()) {
        for &p in error_rate.iter() {
            println!("Progress {}/{}", count, loop_time);

            let mut error_num = 0;

            let bar = ProgressBar::new(loop_num);

            for i in 0..loop_num {
                let distribution = ErrorDistribution::TruncNormal {
                    mean: p,
                    std_dev,
                    seed: i,
                };
                let mut code = RotatedSurfaceCode::new(d, d, distribution, p, i);

                code.initialize();
                code.syndrome_measurement();
                code.run();
                code.decode_mwpm(d);

                let ans = code.logical_value();

                if ans != 0 {
                    error_num += 1;
                }

                bar.inc(1)
            }

            r.push(error_num as f32 / loop_num as f32);

            bar.finish();

            count += 1;
        }
    }

    println!("{:?}", result);
}