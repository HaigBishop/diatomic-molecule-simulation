use std::fs::File;
use std::io::{self, BufRead};
use std::path::Path;
use plotters::prelude::*;

// List of directories and their corresponding dt values
const DIRECTORIES: &[(&str, f32)] = &[
    ("harmonic_osc_dt0_temp298_diff", 0.5),
    ("harmonic_osc_dt1_temp298_diff", 1.0),
    ("harmonic_osc_dt2_temp298_diff", 2.0),
    ("harmonic_osc_dt4_temp298_diff", 4.0),
    ("harmonic_osc_dt5_temp298_diff", 5.0),
    ("harmonic_osc_dt8_temp298_diff", 8.0),
    ("harmonic_osc_dt10_temp298_diff", 10.0),
];

// Plot size and style constants
const PLOT_SIZE: (u32, u32) = (1200, 800);
const FONT_SIZE: u32 = 30;

fn read_max_displacement(directory: &str) -> Result<f32, Box<dyn std::error::Error>> {
    let filepath = format!("output/{}/simulation_recording.csv", directory);
    let file = File::open(filepath)?;
    let reader = io::BufReader::new(file);
    
    let mut max_displacement = 0.0f32;
    let mut is_header = true;
    
    for line in reader.lines() {
        let line = line?;
        
        // Skip the header row
        if is_header {
            is_header = false;
            continue;
        }
        
        // Parse the CSV data
        let fields: Vec<&str> = line.split(',').collect();
        if fields.len() >= 2 {
            if let Ok(displacement) = fields[1].parse::<f32>() {
                max_displacement = max_displacement.max(displacement.abs());
            }
        }
    }
    
    Ok(max_displacement)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Collect dt values and their corresponding max displacements
    let mut dt_displacement_data = Vec::new();
    
    println!("Reading maximum displacement values from CSV files...");
    for (dir, dt) in DIRECTORIES {
        match read_max_displacement(dir) {
            Ok(max_disp) => {
                println!("Directory: {}, Δt: {}, Max Displacement: {}", dir, dt, max_disp);
                dt_displacement_data.push((*dt, max_disp));
            },
            Err(e) => {
                eprintln!("Error reading from directory {}: {}", dir, e);
            }
        }
    }
    
    // Sort data by dt for proper line plotting
    dt_displacement_data.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    
    // Create the output directory if it doesn't exist
    let output_dir = "output/dt_analysis";
    if !Path::new(output_dir).exists() {
        std::fs::create_dir_all(output_dir)?;
    }
    
    let output_file = format!("{}/max_displacement_vs_dt.png", output_dir);
    
    // Set up the plot
    let root = BitMapBackend::new(&output_file, PLOT_SIZE).into_drawing_area();
    root.fill(&WHITE)?;
    
    // Determine the x and y axis ranges
    let dt_min = dt_displacement_data.iter().map(|(dt, _)| *dt).fold(f32::INFINITY, f32::min);
    let dt_max = dt_displacement_data.iter().map(|(dt, _)| *dt).fold(f32::NEG_INFINITY, f32::max);
    let disp_min = 0.0; // Typically start at zero for error measurements
    let disp_max = dt_displacement_data.iter().map(|(_, disp)| *disp).fold(f32::NEG_INFINITY, f32::max) * 1.1;
    
    // Create the chart
    let mut chart = ChartBuilder::on(&root)
        .margin(80)
        .x_label_area_size(60)
        .y_label_area_size(120)
        .build_cartesian_2d(dt_min..dt_max, disp_min..disp_max)?;
    
    chart.configure_mesh()
        .x_desc("Time Step (Δt) (au)")
        .y_desc("Maximum Displacement Error (Å)")
        .x_label_style(("sans-serif", FONT_SIZE).into_font())
        .y_label_style(("sans-serif", FONT_SIZE).into_font())
        .axis_desc_style(("sans-serif", FONT_SIZE).into_font())
        .draw()?;
    
    // Draw line connecting the data points
    chart.draw_series(LineSeries::new(
        dt_displacement_data.iter().map(|(dt, disp)| (*dt, *disp)),
        RED.stroke_width(3),
    ))?
    .label("Maximum Displacement Error")
    .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], RED.stroke_width(3)));
    
    // Add data points
    chart.draw_series(
        dt_displacement_data.iter().map(|(dt, disp)| {
            Circle::new((*dt, *disp), 5, BLUE.filled())
        })
    )?;
    
    // Configure and draw the legend
    chart.configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .label_font(("sans-serif", FONT_SIZE).into_font())
        .position(SeriesLabelPosition::UpperLeft)
        .draw()?;
    
    // Save the result
    root.present()?;
    println!("Plot has been saved to '{}'", output_file);
    
    Ok(())
} 