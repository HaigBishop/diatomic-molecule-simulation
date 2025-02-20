/*
Simple atomic simulations in Rust.
Displacement is defined as the difference between the current bond length and the equilibrium bond length.
*/

// Import the simulation module
mod simulation;
use simulation::{SimulationState, write_history, print_state, plot_histories};

// Physical constants in atomic units (au)
const M_AU: f32 = 9.114400E+02;      // Mass
const K_AU: f32 = 3.665358E-01;      // Force constant
const R0_A0: f32 = 7.177700E-02;     // Initial bond length displacement 

// Time steps
const DT_AU: f32 = 2.0;
const EXP_LEN_AU: f32 = 3200.0;
const N_STEPS: i32 = (EXP_LEN_AU / DT_AU) as i32;
const PRINT_FREQ: i32 = 200;

// Experiment name (output directory)
const EXP_NAME_SIM: &str = "harmonic_osc_dt2_sim";
const EXP_NAME_CALC: &str = "harmonic_osc_dt2_calc";
const EXP_NAME_DIFF: &str = "harmonic_osc_dt2_diff";

fn init_harmonic_osc() -> SimulationState {
    SimulationState {
        time: 0.0,
        displacement: R0_A0,
        force: -K_AU * R0_A0,
        acceleration: -K_AU * R0_A0 / M_AU,
        velocity: 0.0,
        kinetic_e: 0.0,
        potential_e: 0.5 * K_AU * R0_A0.powi(2),
        total_e: 0.5 * K_AU * R0_A0.powi(2),
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

fn calculate_harmonic_osc_state(state: &mut SimulationState) {
    state.time += DT_AU;
    let omega = (K_AU / M_AU).sqrt();
    let omega_time_cos = (omega * state.time).cos();
    let omega_time_sin = (omega * state.time).sin();
    state.displacement = R0_A0 * omega_time_cos;
    state.force = -K_AU * state.displacement;
    state.acceleration = -R0_A0 * omega * omega * omega_time_cos;
    state.velocity = -R0_A0 * omega * omega_time_sin;
    state.kinetic_e = 0.5 * K_AU * R0_A0 * R0_A0 * omega_time_sin * omega_time_sin;
    state.potential_e = 0.5 * K_AU * R0_A0 * R0_A0 * omega_time_cos * omega_time_cos;
    state.total_e = state.kinetic_e + state.potential_e;
}

fn get_difference_of_histories(history_1: &Vec<SimulationState>, history_2: &Vec<SimulationState>) -> Vec<SimulationState> {
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

fn main() {
    // Harmonic Oscillator Simulation //////////////////////////////////
    // Initialize the simulation system
    let mut state = init_harmonic_osc();
    // Create a vector to store the history of states
    let mut sim_history: Vec<SimulationState> = Vec::new();
    // For N_STEPS steps
    for sim_iter in 0..N_STEPS {
        // Simulate the next state
        simulate_harmonic_osc_state(&mut state);
        // Store the current state in the sim_history
        sim_history.push(state.clone());
        // Print state
        if sim_iter % PRINT_FREQ == 0 {
            print_state(&state);
        }
    }
    // Plot all histories in the harmonic_osc directory
    plot_histories(&sim_history, EXP_NAME_SIM, "png");
    // Write sim_history to file
    write_history(&sim_history, EXP_NAME_SIM);


    // Harmonic Oscillator Calculation //////////////////////////////////
    // Initialize the predicted simulation system
    let mut state = init_harmonic_osc();
    // Create a vector to store the history of states
    let mut calc_history: Vec<SimulationState> = Vec::new();
    // For N_STEPS steps
    for sim_iter in 0..N_STEPS {
        // Calculate the next state
        calculate_harmonic_osc_state(&mut state);
        // Store the current state in the calc_history
        calc_history.push(state.clone());
        // Print state
        if sim_iter % PRINT_FREQ == 0 {
            print_state(&state);
        }
    }
    // Plot all histories in the harmonic_osc directory
    plot_histories(&calc_history, EXP_NAME_CALC, "png");
    // Write calc_history to file
    write_history(&calc_history, EXP_NAME_CALC);


    // Harmonic Oscillator Difference //////////////////////////////////
    let difference_of_histories = get_difference_of_histories(&sim_history, &calc_history);
    plot_histories(&difference_of_histories, EXP_NAME_DIFF, "png");
    write_history(&difference_of_histories, EXP_NAME_DIFF);

    println!("Done!");
}


