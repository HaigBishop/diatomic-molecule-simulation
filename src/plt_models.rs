use plotters::prelude::*;

// Hydrogen parameters in atomic units (au)
const H2_K_AU: f32 = 3.665358E-01;      // Force constant (harmonic & Morse)
const H2_D_AU: f32 = 1.818446E-01;      // Dissociation energy (Morse)
const H2_ALPHA_AU: f32 = 1.003894E+00;  // Bond strength (Morse)

// Harmonic potential: V(x) = (1/2) k x²
fn harmonic_potential(x: f32) -> f32 {
    0.5 * H2_K_AU * x * x
}

// Morse potential: V(x) = D (1 - exp(-α x))²
// Here x is defined as the displacement from equilibrium (with x = 0 at the minimum).
fn morse_potential(x: f32) -> f32 {
    H2_D_AU * (1.0 - (-H2_ALPHA_AU * x).exp()).powi(2)
}


fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Create output directory if it doesn't exist
    std::fs::create_dir_all("./output/morse_vs_harm")?;
    
    // Create a drawing area for the plot (PNG file)
    let root = BitMapBackend::new("./output/morse_vs_harm/morse_vs_harm.png", (2560, 1440)).into_drawing_area();
    root.fill(&WHITE)?;

    // Define a range for displacements (in atomic units, au).
    let x_min = -2.0;
    let x_max = 6.0;
    let n_points = 500;
    let dx = (x_max - x_min) / (n_points as f32);

    // Generate data points for both potentials.
    let mut harmonic_data = Vec::with_capacity(n_points + 1);
    let mut morse_data = Vec::with_capacity(n_points + 1);
    for i in 0..=n_points {
        let x = x_min + i as f32 * dx;
        harmonic_data.push((x, harmonic_potential(x)));
        morse_data.push((x, morse_potential(x)));
    }

    // Define the vertical range.
    // Since the harmonic potential grows quadratically, we base the y_max on its value at x_max,
    // adding a 10% margin.
    let y_max = harmonic_potential(x_max) * 1.1;

    // Build the chart.
    let mut chart = ChartBuilder::on(&root)
        // .caption("H2: Morse Potential vs Harmonic Oscillator", ("sans-serif", 100))
        .margin(80)
        .x_label_area_size(100)
        .y_label_area_size(140)
        .build_cartesian_2d(x_min..x_max, 0.0..y_max)?;

    chart.configure_mesh()
        .x_desc("Displacement (Å)")
        .y_desc("Potential Energy (au)")
        .x_label_style(("sans-serif", 40))
        .y_label_style(("sans-serif", 40))
        .draw()?;

    // Draw the harmonic potential curve (in blue)
    chart.draw_series(LineSeries::new(
        harmonic_data,
        BLUE.stroke_width(3),
    ))?
    .label("Harmonic")
    .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLUE));

    // Draw the Morse potential curve (in red)
    chart.draw_series(LineSeries::new(
        morse_data,
        RED.stroke_width(3),
    ))?
    .label("Morse")
    .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));

    // Configure and draw the legend with a larger font size.
    chart.configure_series_labels()
        .border_style(&BLACK)
        .label_font(("sans-serif", 40))
        .draw()?;

    // Save the result.
    root.present()?;
    println!("Plot has been saved to './output/morse_vs_harm/morse_vs_harm.png'");

    Ok(())
}
