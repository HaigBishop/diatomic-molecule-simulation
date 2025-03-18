use std::fs::File;
use std::io::{self, BufRead};
use std::path::Path;
use plotters::prelude::*;

// List of directory pairs and their corresponding temperatures
const DIRECTORIES: &[(f32, &str, &str)] = &[
    (100.0, "morse_potential_dt1_temp100_sim", "harmonic_osc_dt1_temp100_sim"),
    (200.0, "morse_potential_dt1_temp200_sim", "harmonic_osc_dt1_temp200_sim"),
    (300.0, "morse_potential_dt1_temp300_sim", "harmonic_osc_dt1_temp300_sim"),
    (500.0, "morse_potential_dt1_temp500_sim", "harmonic_osc_dt1_temp500_sim"),
    (1000.0, "morse_potential_dt1_temp1000_sim", "harmonic_osc_dt1_temp1000_sim"),
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
    // Collect temperature values and corresponding max displacements for both models
    let mut temp_morse_displacement_data = Vec::new();
    let mut temp_harmonic_displacement_data = Vec::new();
    
    println!("Reading maximum displacement values from CSV files...");
    for (temp, morse_dir, harmonic_dir) in DIRECTORIES {
        // Process Morse potential data
        match read_max_displacement(morse_dir) {
            Ok(max_disp) => {
                println!("Temperature: {}, Morse Directory: {}, Max Displacement: {}", temp, morse_dir, max_disp);
                temp_morse_displacement_data.push((*temp, max_disp));
            },
            Err(e) => {
                eprintln!("Error reading from directory {}: {}", morse_dir, e);
            }
        }
        
        // Process Harmonic oscillator data
        match read_max_displacement(harmonic_dir) {
            Ok(max_disp) => {
                println!("Temperature: {}, Harmonic Directory: {}, Max Displacement: {}", temp, harmonic_dir, max_disp);
                temp_harmonic_displacement_data.push((*temp, max_disp));
            },
            Err(e) => {
                eprintln!("Error reading from directory {}: {}", harmonic_dir, e);
            }
        }
    }
    
    // Sort data by temperature for proper line plotting
    temp_morse_displacement_data.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    temp_harmonic_displacement_data.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    
    // Create the output directory if it doesn't exist
    let output_dir = "output/temp_analysis";
    if !Path::new(output_dir).exists() {
        std::fs::create_dir_all(output_dir)?;
    }
    
    let output_file = format!("{}/max_displacement_vs_temp.png", output_dir);
    
    // Set up the plot
    let root = BitMapBackend::new(&output_file, PLOT_SIZE).into_drawing_area();
    root.fill(&WHITE)?;
    
    // Determine the x and y axis ranges
    let temp_min = DIRECTORIES.iter().map(|(temp, _, _)| *temp).fold(f32::INFINITY, f32::min);
    let temp_max = DIRECTORIES.iter().map(|(temp, _, _)| *temp).fold(f32::NEG_INFINITY, f32::max);
    
    let disp_values: Vec<f32> = temp_morse_displacement_data.iter().map(|(_, disp)| *disp)
        .chain(temp_harmonic_displacement_data.iter().map(|(_, disp)| *disp))
        .collect();
    let disp_min = 0.0; // Typically start at zero for displacement measurements
    let disp_max = disp_values.iter().fold(f32::NEG_INFINITY, |max, &val| max.max(val)) * 1.1;
    
    // Create the chart
    let mut chart = ChartBuilder::on(&root)
        .margin(80)
        .x_label_area_size(60)
        .y_label_area_size(120)
        .build_cartesian_2d(temp_min..temp_max, disp_min..disp_max)?;
    
    chart.configure_mesh()
        .x_desc("Temperature (K)")
        .y_desc("Maximum Displacement (Å)")
        .x_label_style(("sans-serif", FONT_SIZE).into_font())
        .y_label_style(("sans-serif", FONT_SIZE).into_font())
        .axis_desc_style(("sans-serif", FONT_SIZE).into_font())
        .draw()?;
    
    // Draw line for Morse potential data
    chart.draw_series(LineSeries::new(
        temp_morse_displacement_data.iter().map(|(temp, disp)| (*temp, *disp)),
        RED.stroke_width(3),
    ))?
    .label("Morse Potential")
    .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], RED.stroke_width(3)));
    
    // Add data points for Morse potential
    chart.draw_series(
        temp_morse_displacement_data.iter().map(|(temp, disp)| {
            Circle::new((*temp, *disp), 5, RED.filled())
        })
    )?;
    
    // Draw line for Harmonic oscillator data
    chart.draw_series(LineSeries::new(
        temp_harmonic_displacement_data.iter().map(|(temp, disp)| (*temp, *disp)),
        BLUE.stroke_width(3),
    ))?
    .label("Harmonic Oscillator")
    .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], BLUE.stroke_width(3)));
    
    // Add data points for Harmonic oscillator
    chart.draw_series(
        temp_harmonic_displacement_data.iter().map(|(temp, disp)| {
            Circle::new((*temp, *disp), 5, BLUE.filled())
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