/*
Simple atomic simulations in Rust.
 - Here, displacement is defined as the difference between the current bond length and the bond length at the equilibrium position.
 - Units used in simulation and calculation are atomic units (au).
*/

// Import the simulation module
mod simulation;
use simulation::{SimulationState, write_history, print_state, 
    plot_histories, get_difference_of_histories, plot_histories_overlayed};

// Conversion factors and constants
const KB: f32 = 1.3806488E-23;
const A0_TO_M: f32 = 5.2917721092E-11;

// H2 physical constants in atomic units (au)
const M_AU: f32 = 9.114400E+02;      // Mass (harmonic & morse)
const K_AU: f32 = 3.665358E-01;      // Force constant (harmonic & morse)
const D_AU: f32 = 1.818446E-01;      // Dissociation energy (morse)
const ALPHA_AU: f32 = 1.003894E+00;      // Bond strength (morse)

// H2 physical constants in atomic units (SI)
const K_SI: f32 = 5.706570E+02;
const D_SI: f32 = 7.928147E-19;
const ALPHA_SI: f32 = 1.897085E+10;

// Time steps
const DT_AU: f32 = 20.0;
const EXP_LEN_AU: f32 = 3200.0;
const N_STEPS: i32 = (EXP_LEN_AU / DT_AU) as i32;
const PRINT_FREQ: i32 = 200;

// Temperature in Kelvin
const TEMP_K: f32 = 100.0; // 100.0 200.0 298.15 500.0 1000.0


fn init_harmonic_osc(r0_a0: f32) -> SimulationState {
    SimulationState {
        time: 0.0,
        displacement: r0_a0,
        force: -K_AU * r0_a0,
        acceleration: -K_AU * r0_a0 / M_AU,
        velocity: 0.0,
        kinetic_e: 0.0,
        potential_e: 0.5 * K_AU * r0_a0.powi(2),
        total_e: 0.5 * K_AU * r0_a0.powi(2),
    }
}

fn init_morse_potential(r0_a0: f32) -> SimulationState {
    let exp_alpha_r0 = f32::exp(-ALPHA_AU * r0_a0);
    let init_force = -2.0 * D_AU * ALPHA_AU * exp_alpha_r0 * (1.0 - exp_alpha_r0);
    let exp_alpha_r0_sq = (1.0 - exp_alpha_r0).powi(2);
    SimulationState {
        time: 0.0,
        displacement: r0_a0,
        force: init_force,
        acceleration: init_force / M_AU,
        velocity: 0.0,
        kinetic_e: 0.0,
        potential_e: D_AU * exp_alpha_r0_sq,
        total_e: D_AU * exp_alpha_r0_sq,
    }
}

fn simulate_harmonic_osc_state(state: &mut SimulationState) {
    state.time += DT_AU;
    let old_acceleration = state.acceleration;
    state.displacement += state.velocity * DT_AU + 0.5 * old_acceleration * DT_AU * DT_AU;
    state.force = -K_AU * state.displacement;
    state.acceleration = state.force / M_AU;
    state.velocity += 0.5 * (old_acceleration + state.acceleration) * DT_AU;
    state.kinetic_e = 0.5 * M_AU * state.velocity.powi(2);
    state.potential_e = 0.5 * K_AU * state.displacement.powi(2);
    state.total_e = state.kinetic_e + state.potential_e;
}

fn simulate_morse_potential_state(state: &mut SimulationState) {
    state.time += DT_AU;
    let old_acceleration = state.acceleration;
    state.displacement += state.velocity * DT_AU + 0.5 * old_acceleration * DT_AU * DT_AU; // I assume this is correct
    let exp_alpha_r0 = f32::exp(-ALPHA_AU * state.displacement);
    let next_force = -2.0 * D_AU * ALPHA_AU * exp_alpha_r0 * (1.0 - exp_alpha_r0);
    let exp_alpha_r0_sq = (1.0 - exp_alpha_r0).powi(2);
    state.force = next_force;
    state.acceleration = next_force / M_AU;
    state.velocity += 0.5 * (old_acceleration + state.acceleration) * DT_AU; // I assume this is correct
    state.kinetic_e = 0.5 * M_AU * state.velocity.powi(2);
    state.potential_e = D_AU * exp_alpha_r0_sq;
    state.total_e = state.kinetic_e + state.potential_e;
    
}

fn calculate_harmonic_osc_state(state: &mut SimulationState, r0_a0: f32) {
    state.time += DT_AU;
    let omega = (K_AU / M_AU).sqrt();
    let omega_time_cos = (omega * state.time).cos();
    let omega_time_sin = (omega * state.time).sin();
    state.displacement = r0_a0 * omega_time_cos;
    state.force = -K_AU * state.displacement;
    state.acceleration = -r0_a0 * omega * omega * omega_time_cos;
    state.velocity = -r0_a0 * omega * omega_time_sin;
    state.kinetic_e = 0.5 * K_AU * r0_a0 * r0_a0 * omega_time_sin * omega_time_sin;
    state.potential_e = 0.5 * K_AU * r0_a0 * r0_a0 * omega_time_cos * omega_time_cos;
    state.total_e = state.kinetic_e + state.potential_e;
}

fn main() {
    // Generate experiment names (output directories)
    let exp_name_base = format!("dt{}_temp{}", DT_AU as i32, TEMP_K as i32);
    let exp_name_harm_sim = format!("harmonic_osc_{}_sim", exp_name_base);
    let exp_name_harm_calc = format!("harmonic_osc_{}_calc", exp_name_base);
    let exp_name_harm_diff = format!("harmonic_osc_{}_diff", exp_name_base);
    let exp_name_morse_sim = format!("morse_potential_{}_sim", exp_name_base);
    let exp_name_morse_overlay = format!("morse_potential_{}_overlay", exp_name_base);

    println!("Calculation Initial Displacements...");
    // Calculate the initial displacement in SI
    let r0_si_harm: f32 = ((2.0 * KB * TEMP_K) / K_SI).sqrt();
    let r0_si_morse: f32 = (1.0 - (K_SI * r0_si_harm * r0_si_harm / (2.0 * D_SI)).sqrt()).ln() / (-ALPHA_SI);
    // Convert the initial displacement to AU
    let r0_a0_harm: f32 = r0_si_harm / A0_TO_M; // Initial bond length displacement  (harmonic)
    let r0_a0_morse: f32 = r0_si_morse / A0_TO_M;  // Initial bond length displacement  (morse)
    println!("r0_a0_harm: {}", r0_a0_harm);
    println!("r0_a0_morse: {}", r0_a0_morse);

    println!("\nStarting Simulation...");

    // Harmonic Oscillator Simulation //////////////////////////////////
    // Initialize the simulation system
    let mut state = init_harmonic_osc(r0_a0_harm);
    // Create a vector to store the history of states
    let mut harm_sim_history: Vec<SimulationState> = Vec::new();
    // For N_STEPS steps
    for sim_iter in 0..N_STEPS {
        // Simulate the next state
        simulate_harmonic_osc_state(&mut state);
        // Store the current state in the harm_sim_history
        harm_sim_history.push(state.clone());
        // Print state
        if sim_iter % PRINT_FREQ == 0 {
            print_state(&state);
        }
    }
    // Plot all histories in the harmonic_osc directory
    plot_histories(&harm_sim_history, &exp_name_harm_sim, "png");
    // Write harm_sim_history to file
    write_history(&harm_sim_history, &exp_name_harm_sim);


    // Harmonic Oscillator Calculation //////////////////////////////////
    // Initialize the predicted simulation system
    let mut state = init_harmonic_osc(r0_a0_harm);
    // Create a vector to store the history of states
    let mut harm_calc_history: Vec<SimulationState> = Vec::new();
    // For N_STEPS steps
    for sim_iter in 0..N_STEPS {
        // Calculate the next state
        calculate_harmonic_osc_state(&mut state, r0_a0_harm);
        // Store the current state in the harm_calc_history
        harm_calc_history.push(state.clone());
        // Print state
        if sim_iter % PRINT_FREQ == 0 {
            print_state(&state);
        }
    }
    // Plot all histories in the harmonic_osc directory
    plot_histories(&harm_calc_history, &exp_name_harm_calc, "png");
    // Write harm_calc_history to file
    write_history(&harm_calc_history, &exp_name_harm_calc);



    // Harmonic Oscillator Difference //////////////////////////////////
    let difference_of_harm_histories = get_difference_of_histories(&harm_sim_history, &harm_calc_history);
    plot_histories(&difference_of_harm_histories, &exp_name_harm_diff, "png");
    write_history(&difference_of_harm_histories, &exp_name_harm_diff);



    // Morse Potential Simulation //////////////////////////////////
    // Initialize the simulation system
    let mut state = init_morse_potential(r0_a0_morse);
    // Create a vector to store the history of states
    let mut morse_sim_history: Vec<SimulationState> = Vec::new();
    // For N_STEPS steps
    for sim_iter in 0..N_STEPS {
        // Simulate the next state
        simulate_morse_potential_state(&mut state);
        // Store the current state in the morse_sim_history
        morse_sim_history.push(state.clone());
        // Print state
        if sim_iter % PRINT_FREQ == 0 {
            print_state(&state);
        }
    }
    // Plot all histories in the morse_potential directory
    plot_histories(&morse_sim_history, &exp_name_morse_sim, "png");
    // Plot the morse history overlayed on the harmonic history
    plot_histories_overlayed(&morse_sim_history, &harm_sim_history, &exp_name_morse_overlay, "png");
    // Write morse_sim_history to file
    write_history(&morse_sim_history, &exp_name_morse_sim);




    println!("Done!");
}


