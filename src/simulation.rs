use std::fs::OpenOptions;
use std::io::Write;
use plotters::prelude::*;

/// Helper function for formatting x-axis tick labels as integers.
fn integer_tick_formatter(v: &f32) -> String {
    format!("{:.0}", v)
}

#[derive(Clone)]
pub struct SimulationState {
    pub time: f32,
    pub displacement: f32,
    pub force: f32,
    pub acceleration: f32,
    pub velocity: f32,
    pub kinetic_e: f32,
    pub potential_e: f32,
    pub total_e: f32,
}

const PLOT_SIZE: (u32, u32) = (1200, 800);

pub fn write_history(history: &Vec<SimulationState>, exp_name: &str) {
    let file_path = format!("output/{}/simulation_recording.csv", exp_name);
    
    // Only try to delete if the file exists
    if std::path::Path::new(&file_path).exists() {
        if let Err(e) = std::fs::remove_file(&file_path) {
            eprintln!("Failed to delete file: {}", e);
        }
    }

    // Open the file in append mode, create it if it doesn't exist
    let mut file = match OpenOptions::new()
        .append(true)
        .create(true)
        .open(file_path) {
            Ok(file) => file,
            Err(e) => {
                eprintln!("Failed to open file: {}", e);
                return;
            }
        };

    // Write the header to the file
    if let Err(e) = writeln!(
        file, "Time,Displacement,Force,Acceleration,Velocity,Kinetic Energy,Potential Energy,Total Energy"
    ) {
        eprintln!("Failed to write to file: {}", e);
    }

    // Write the history to the file in CSV format
    for state in history {
        // Write the state to the file in CSV format
        if let Err(e) = writeln!(
            file, "{},{},{},{},{},{},{},{}",
            state.time, state.displacement, state.force, state.acceleration, 
            state.velocity, state.kinetic_e, state.potential_e, state.total_e
        ) {
            eprintln!("Failed to write to file: {}", e);
        }
    }
}

pub fn print_state(state: &SimulationState) {
    println!(
        "Time: {}, Displacement: {}, Total Energy: {}",
        state.time, state.displacement, state.total_e
    );
}

pub fn plot_history(history: &Vec<SimulationState>, y_col: &str, filepath: &str, _source_name: &str) {
    // Create the output directory if it doesn't exist
    let output_dir = "output";
    if !std::path::Path::new(output_dir).exists() {
        std::fs::create_dir(output_dir).unwrap();
    }

    let y_values: Vec<f32> = match y_col {
        "displacement" => history.iter().map(|s| s.displacement).collect(),
        "force" => history.iter().map(|s| s.force).collect(),
        "acceleration" => history.iter().map(|s| s.acceleration).collect(),
        "velocity" => history.iter().map(|s| s.velocity).collect(),
        "kinetic_e" => history.iter().map(|s| s.kinetic_e).collect(),
        "potential_e" => history.iter().map(|s| s.potential_e).collect(),
        "total_e" => history.iter().map(|s| s.total_e).collect(),
        _ => panic!("Invalid y_col value"),
    };

    let y_min = y_values.iter().cloned().fold(f32::INFINITY, f32::min);
    let y_max = y_values.iter().cloned().fold(f32::NEG_INFINITY, f32::max);
    
    // Add 10% padding to the y-axis range
    let y_range = y_max - y_min;
    let y_padding = y_range * 0.1;
    let y_min = y_min - y_padding;
    let y_max = y_max + y_padding;

    // Compute common scaling factor so that tick labels have at most one digit before the decimal place
    let y_abs_max = y_min.abs().max(y_max.abs());
    let y_scale = if y_abs_max == 0.0 { 1.0 } else { 10f32.powf(y_abs_max.log10().floor()) };

    // Modify y-axis label to include the magnitude multiplier if it's not 1 (i.e. e0)
    let base_label = match y_col {
        "displacement" => ("Displacement", "Å"),
        "force" => ("Force", "au"),
        "acceleration" => ("Acceleration", "au"),
        "velocity" => ("Velocity", "au"),
        "kinetic_e" => ("Kinetic Energy", "au"),
        "potential_e" => ("Potential Energy", "au"),
        "total_e" => ("Total Energy", "au"),
        _ => (y_col, ""),
    };

    let y_label = if (y_scale - 1.0).abs() < std::f32::EPSILON {
        format!("{} ({})", base_label.0, base_label.1)
    } else {
        let exponent = y_scale.log10() as i32;
        format!("{} (10^{} {})", base_label.0, exponent, base_label.1)
    };

    let root_area = BitMapBackend::new(filepath, PLOT_SIZE)
        .into_drawing_area();
    root_area.fill(&WHITE).unwrap();

    let mut chart = ChartBuilder::on(&root_area)
        // .caption(format!("{} History ({})", source_name, y_col.to_string().chars().next().unwrap().to_uppercase().collect::<String>() + &y_col[1..]), ("sans-serif", 25).into_font())
        .margin_left(40)
        .margin_right(30)
        .margin_bottom(30)
        .x_label_area_size(85)
        .y_label_area_size(120)
        .build_cartesian_2d(
            0.0..history.last().unwrap().time,
            y_min..y_max,
        )
        .unwrap();

    chart.configure_mesh()
        .x_desc("Time (au)")
        .y_desc(y_label)
        .x_label_formatter(&integer_tick_formatter)
        .y_label_formatter(&|v: &f32| format!("{:.1}", *v / y_scale))
        .x_labels(10)
        .y_labels(10)
        .x_label_style(("sans-serif", 30).into_font())
        .y_label_style(("sans-serif", 30).into_font())
        .draw()
        .unwrap();

    chart
        .draw_series(LineSeries::new(
            history.iter().zip(y_values.iter()).map(|(s, y)| (s.time, *y)),
            &RED,
        ))
        .unwrap()
        .label(y_col.to_string().chars().next().unwrap().to_uppercase().collect::<String>() + &y_col[1..])
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));

    chart
        .configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .label_font(("sans-serif", 30).into_font())
        .position(SeriesLabelPosition::UpperRight)
        .draw()
        .unwrap();
}

fn plot_energy_history(history: &Vec<SimulationState>, filepath: &str, _source_name: &str) {
    let root_area = BitMapBackend::new(filepath, PLOT_SIZE)
        .into_drawing_area();
    root_area.fill(&WHITE).unwrap();

    // Find min and max energy values across all energy types
    let e_min = history.iter()
        .flat_map(|s| vec![s.kinetic_e, s.potential_e, s.total_e])
        .fold(f32::INFINITY, f32::min);
    let e_max = history.iter()
        .flat_map(|s| vec![s.kinetic_e, s.potential_e, s.total_e])
        .fold(f32::NEG_INFINITY, f32::max);

    // Add 10% padding to the y-axis range
    let e_range = e_max - e_min;
    let e_padding = e_range * 0.1;
    let e_min = e_min - e_padding;
    let e_max = e_max + e_padding;

    // Compute scaling factor for energy y-axis
    let y_abs_max = e_min.abs().max(e_max.abs());
    let y_scale = if y_abs_max == 0.0 { 1.0 } else { 10f32.powf(y_abs_max.log10().floor()) };

    let y_label = if (y_scale - 1.0).abs() < std::f32::EPSILON {
        "Energy (au)".to_string()
    } else {
        let exponent = y_scale.log10() as i32;
        format!("Energy (10^{} au)", exponent)
    };

    let mut chart = ChartBuilder::on(&root_area)
        // .caption(format!("{} Energy History", source_name), ("sans-serif", 25).into_font())
        .margin_left(40)
        .margin_right(30)
        .margin_bottom(30)
        .x_label_area_size(85)
        .y_label_area_size(120)
        .build_cartesian_2d(
            0.0..history.last().unwrap().time,
            e_min..e_max,
        )
        .unwrap();

    chart.configure_mesh()
        .x_desc("Time (au)")
        .y_desc(y_label)
        .x_label_formatter(&integer_tick_formatter)
        .y_label_formatter(&|v: &f32| format!("{:.1}", *v / y_scale))
        .x_labels(10)
        .y_labels(10)
        .x_label_style(("sans-serif", 30).into_font())
        .y_label_style(("sans-serif", 30).into_font())
        .draw()
        .unwrap();

    // Draw kinetic energy in red
    chart
        .draw_series(LineSeries::new(
            history.iter().map(|s| (s.time, s.kinetic_e)),
            &RED,
        ))
        .unwrap()
        .label("Kinetic Energy")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));

    // Draw potential energy in blue
    chart
        .draw_series(LineSeries::new(
            history.iter().map(|s| (s.time, s.potential_e)),
            &BLUE,
        ))
        .unwrap()
        .label("Potential Energy")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLUE));

    // Draw total energy in green
    chart
        .draw_series(LineSeries::new(
            history.iter().map(|s| (s.time, s.total_e)),
            &GREEN,
        ))
        .unwrap()
        .label("Total Energy")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &GREEN));

    chart
        .configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .label_font(("sans-serif", 30).into_font())
        .position(SeriesLabelPosition::UpperRight)
        .draw()
        .unwrap();
}

pub fn plot_histories(history: &Vec<SimulationState>, exp_name: &str, image_type: &str, source_name: &str) {
    // Create the experiment directory inside output
    let exp_dir = format!("output/{}", exp_name);
    std::fs::create_dir_all(&exp_dir).unwrap();

    // Plot these columns
    plot_history(history, "displacement", 
        &format!("{}/displacement.{}", exp_dir, image_type), source_name);
    plot_history(history, "force", 
        &format!("{}/force.{}", exp_dir, image_type), source_name);
    plot_history(history, "acceleration", 
        &format!("{}/acceleration.{}", exp_dir, image_type), source_name);
    plot_history(history, "velocity", 
        &format!("{}/velocity.{}", exp_dir, image_type), source_name);
    plot_history(history, "total_e", 
        &format!("{}/total_e.{}", exp_dir, image_type), source_name);

    // Plot energy history
    plot_energy_history(history, 
        &format!("{}/energies.{}", exp_dir, image_type), source_name);
}

pub fn plot_histories_overlayed(history_1: &Vec<SimulationState>, history_2: &Vec<SimulationState>, 
                               exp_name: &str, image_type: &str, source_name: &str) {
    // Create the experiment directory inside output
    let exp_dir = format!("output/{}", exp_name);
    std::fs::create_dir_all(&exp_dir).unwrap();

    // Plot these columns overlayed
    plot_history_overlayed(history_1, history_2, "displacement", 
        &format!("{}/displacement.{}", exp_dir, image_type), source_name);
    plot_history_overlayed(history_1, history_2, "force", 
        &format!("{}/force.{}", exp_dir, image_type), source_name);
    plot_history_overlayed(history_1, history_2, "acceleration", 
        &format!("{}/acceleration.{}", exp_dir, image_type), source_name);
    plot_history_overlayed(history_1, history_2, "velocity", 
        &format!("{}/velocity.{}", exp_dir, image_type), source_name);
    
    // Plot energy histories overlayed
    plot_energy_histories_overlayed(history_1, history_2, 
        &format!("{}/energies.{}", exp_dir, image_type), source_name);
}

pub fn get_difference_of_histories(history_1: &Vec<SimulationState>, history_2: &Vec<SimulationState>) -> Vec<SimulationState> {
    let mut difference_of_histories: Vec<SimulationState> = Vec::new();
    for i in 0..history_1.len() {
        let difference_of_states = SimulationState {
            time: history_1[i].time,  // Keep time unchanged
            displacement: history_1[i].displacement - history_2[i].displacement,
            force: history_1[i].force - history_2[i].force,
            acceleration: history_1[i].acceleration - history_2[i].acceleration,
            velocity: history_1[i].velocity - history_2[i].velocity,
            kinetic_e: history_1[i].kinetic_e - history_2[i].kinetic_e,
            potential_e: history_1[i].potential_e - history_2[i].potential_e,
            total_e: history_1[i].total_e - history_2[i].total_e,
        };
        difference_of_histories.push(difference_of_states);
    }
    difference_of_histories
}

pub fn plot_history_overlayed(history_1: &Vec<SimulationState>, history_2: &Vec<SimulationState>, y_col: &str, filepath: &str, _source_name: &str) {
    // Create the output directory if it doesn't exist
    let output_dir = "output";
    if !std::path::Path::new(output_dir).exists() {
        std::fs::create_dir(output_dir).unwrap();
    }

    // Get y values for both histories
    let get_y_values = |history: &Vec<SimulationState>| -> Vec<f32> {
        match y_col {
            "displacement" => history.iter().map(|s| s.displacement).collect(),
            "force" => history.iter().map(|s| s.force).collect(),
            "acceleration" => history.iter().map(|s| s.acceleration).collect(),
            "velocity" => history.iter().map(|s| s.velocity).collect(),
            "kinetic_e" => history.iter().map(|s| s.kinetic_e).collect(),
            "potential_e" => history.iter().map(|s| s.potential_e).collect(),
            "total_e" => history.iter().map(|s| s.total_e).collect(),
            _ => panic!("Invalid y_col value"),
        }
    };

    let y_values_1 = get_y_values(history_1);
    let y_values_2 = get_y_values(history_2);

    // Find overall min and max
    let y_min = y_values_1.iter()
        .chain(y_values_2.iter())
        .cloned()
        .fold(f32::INFINITY, f32::min);
    let y_max = y_values_1.iter()
        .chain(y_values_2.iter())
        .cloned()
        .fold(f32::NEG_INFINITY, f32::max);
    
    // Add 10% padding to the y-axis range
    let y_range = y_max - y_min;
    let y_padding = y_range * 0.1;
    let y_min = y_min - y_padding;
    let y_max = y_max + y_padding;

    // Compute common scaling factor for the y-axis
    let combined_y_abs_max = y_min.abs().max(y_max.abs());
    let y_scale = if combined_y_abs_max == 0.0 { 1.0 } else { 10f32.powf(combined_y_abs_max.log10().floor()) };

    let base_label = match y_col {
        "displacement" => ("Displacement", "Å"),
        "force" => ("Force", "au"),
        "acceleration" => ("Acceleration", "au"),
        "velocity" => ("Velocity", "au"),
        "kinetic_e" => ("Kinetic Energy", "au"),
        "potential_e" => ("Potential Energy", "au"),
        "total_e" => ("Total Energy", "au"),
        _ => (y_col, ""),
    };

    let y_label = if (y_scale - 1.0).abs() < std::f32::EPSILON {
        format!("{} ({})", base_label.0, base_label.1)
    } else {
        let exponent = y_scale.log10() as i32;
        format!("{} (10^{} {})", base_label.0, exponent, base_label.1)
    };

    let root_area = BitMapBackend::new(filepath, PLOT_SIZE)
        .into_drawing_area();
    root_area.fill(&WHITE).unwrap();

    let mut chart = ChartBuilder::on(&root_area)
        // .caption(format!("{} History ({})", source_name, y_col.to_string().chars().next().unwrap().to_uppercase().collect::<String>() + &y_col[1..]), ("sans-serif", 25).into_font())
        .margin_left(40)
        .margin_right(30)
        .margin_bottom(30)
        .x_label_area_size(85)
        .y_label_area_size(120)
        .build_cartesian_2d(
            0.0..history_1.last().unwrap().time,
            y_min..y_max,
        )
        .unwrap();

    chart.configure_mesh()
        .x_desc("Time (au)")
        .y_desc(y_label)
        .x_label_formatter(&integer_tick_formatter)
        .y_label_formatter(&|v: &f32| format!("{:.1}", *v / y_scale))
        .x_labels(10)
        .y_labels(10)
        .x_label_style(("sans-serif", 30).into_font())
        .y_label_style(("sans-serif", 30).into_font())
        .draw()
        .unwrap();

    // Draw first history in red
    chart
        .draw_series(LineSeries::new(
            history_1.iter().zip(y_values_1.iter()).map(|(s, y)| (s.time, *y)),
            &RED,
        ))
        .unwrap()
        .label(format!("{} (Morse)", y_col.to_string().chars().next().unwrap().to_uppercase().collect::<String>() + &y_col[1..]))
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));

    // Draw second history in blue
    chart
        .draw_series(LineSeries::new(
            history_2.iter().zip(y_values_2.iter()).map(|(s, y)| (s.time, *y)),
            &BLUE,
        ))
        .unwrap()
        .label(format!("{} (Harmonic)", y_col.to_string().chars().next().unwrap().to_uppercase().collect::<String>() + &y_col[1..]))
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLUE));

    chart
        .configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .label_font(("sans-serif", 30).into_font())
        .position(SeriesLabelPosition::UpperRight)
        .draw()
        .unwrap();
}

fn plot_energy_histories_overlayed(history_1: &Vec<SimulationState>, history_2: &Vec<SimulationState>, filepath: &str, _source_name: &str) {
    let root_area = BitMapBackend::new(filepath, PLOT_SIZE)
        .into_drawing_area();
    root_area.fill(&WHITE).unwrap();

    // Find min and max energy values across all energy types and both histories
    let e_min = history_1.iter()
        .chain(history_2.iter())
        .flat_map(|s| vec![s.kinetic_e, s.potential_e, s.total_e])
        .fold(f32::INFINITY, f32::min);
    let e_max = history_1.iter()
        .chain(history_2.iter())
        .flat_map(|s| vec![s.kinetic_e, s.potential_e, s.total_e])
        .fold(f32::NEG_INFINITY, f32::max);

    // Add 10% padding to the y-axis range
    let e_range = e_max - e_min;
    let e_padding = e_range * 0.1;
    let e_min = e_min - e_padding;
    let e_max = e_max + e_padding;

    // Compute scaling factor for energy y-axis
    let y_abs_max = e_min.abs().max(e_max.abs());
    let y_scale = if y_abs_max == 0.0 { 1.0 } else { 10f32.powf(y_abs_max.log10().floor()) };

    let y_label = if (y_scale - 1.0).abs() < std::f32::EPSILON {
        "Energy (au)".to_string()
    } else {
        let exponent = y_scale.log10() as i32;
        format!("Energy (10^{} au)", exponent)
    };

    let mut chart = ChartBuilder::on(&root_area)
        // .caption(format!("{} Energy History", source_name), ("sans-serif", 25).into_font())
        .margin_left(40)
        .margin_right(30)
        .margin_bottom(30)
        .x_label_area_size(85)
        .y_label_area_size(120)
        .build_cartesian_2d(
            0.0..history_1.last().unwrap().time,
            e_min..e_max,
        )
        .unwrap();

    chart.configure_mesh()
        .x_desc("Time (au)")
        .y_desc(y_label)
        .x_label_formatter(&integer_tick_formatter)
        .y_label_formatter(&|v: &f32| format!("{:.1}", *v / y_scale))
        .x_labels(10)
        .y_labels(10)
        .x_label_style(("sans-serif", 30).into_font())
        .y_label_style(("sans-serif", 30).into_font())
        .draw()
        .unwrap();

    // Draw first history energies (Morse)
    chart
        .draw_series(LineSeries::new(
            history_1.iter().map(|s| (s.time, s.kinetic_e)),
            &RED,
        ))
        .unwrap()
        .label("Kinetic Energy (Morse)")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED));

    chart
        .draw_series(LineSeries::new(
            history_1.iter().map(|s| (s.time, s.potential_e)),
            &BLUE,
        ))
        .unwrap()
        .label("Potential Energy (Morse)")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLUE));

    chart
        .draw_series(LineSeries::new(
            history_1.iter().map(|s| (s.time, s.total_e)),
            &GREEN,
        ))
        .unwrap()
        .label("Total Energy (Morse)")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &GREEN));

    // Draw second history energies (Harmonic) with lighter colors
    chart
        .draw_series(LineSeries::new(
            history_2.iter().map(|s| (s.time, s.kinetic_e)),
            &RED.mix(0.5),
        ))
        .unwrap()
        .label("Kinetic Energy (Harmonic)")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &RED.mix(0.5)));

    chart
        .draw_series(LineSeries::new(
            history_2.iter().map(|s| (s.time, s.potential_e)),
            &BLUE.mix(0.5),
        ))
        .unwrap()
        .label("Potential Energy (Harmonic)")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &BLUE.mix(0.5)));

    chart
        .draw_series(LineSeries::new(
            history_2.iter().map(|s| (s.time, s.total_e)),
            &GREEN.mix(0.5),
        ))
        .unwrap()
        .label("Total Energy (Harmonic)")
        .legend(|(x, y)| PathElement::new(vec![(x, y), (x + 20, y)], &GREEN.mix(0.5)));

    chart
        .configure_series_labels()
        .background_style(&WHITE.mix(0.8))
        .border_style(&BLACK)
        .label_font(("sans-serif", 30).into_font())
        .position(SeriesLabelPosition::UpperRight)
        .draw()
        .unwrap();
}



