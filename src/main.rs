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

// H2 physical constants (au)
const H2_M_AU: f32 = 9.114400E+02;      // Mass (harmonic & morse)
const H2_K_AU: f32 = 3.665358E-01;      // Force constant (harmonic & morse)
const H2_D_AU: f32 = 1.818446E-01;      // Dissociation energy (morse)
const H2_ALPHA_AU: f32 = 1.003894E+00;      // Bond strength (morse)
// H2 physical constants (SI)
const H2_K_SI: f32 = 5.706570E+02;
const H2_D_SI: f32 = 7.928147E-19;
const H2_ALPHA_SI: f32 = 1.897085E+10;

// Hg2 physical constants (au)
const HG2_M_AU: f32 = 1.840841E+05;
const HG2_K_AU: f32 = 1.374407E-03;
const HG2_RSTR_AU: f32 = 6.952302E+00;
const HG2_EPS_AU: f32 = 1.845314E-03;
// Hg2 physical constants (SI)
const HG2_K_SI: f32 = 2.139865E+00;

// Ar2 physical constants (au)
const AR2_M_AU: f32 = 3.641021E+04;
const AR2_K_AU: f32 = 3.232914E-04;
const AR2_RSTR_AU: f32 = 7.107260E+00;
const AR2_EPS_AU: f32 = 4.536240E-04;
// Ar2 physical constants (SI)
const AR2_K_SI: f32 = 5.033442E-01;



// Time steps
const DT_AU: f32 = 10.0;
const EXP_LEN_AU: f32 = 300000.0;
const N_STEPS: i32 = (EXP_LEN_AU / DT_AU) as i32;
const PRINT_FREQ: i32 = 200;

// Temperature in Kelvin
const TEMP_K: f32 = 50.0; // 100.0 200.0 298.15 500.0 1000.0


fn init_harmonic_osc(r0_a0: f32, k_au: f32, m_au: f32) -> SimulationState {
    SimulationState {
        time: 0.0,
        displacement: r0_a0,
        force: -k_au * r0_a0,
        acceleration: -k_au * r0_a0 / m_au,
        velocity: 0.0,
        kinetic_e: 0.0,
        potential_e: 0.5 * k_au * r0_a0.powi(2),
        total_e: 0.5 * k_au * r0_a0.powi(2),
    }
}

fn init_morse_potential(r0_a0: f32, d_au: f32, alpha_au: f32, m_au: f32) -> SimulationState {
    let exp_alpha_r0 = f32::exp(-alpha_au * r0_a0);
    let init_force = -2.0 * d_au * alpha_au * exp_alpha_r0 * (1.0 - exp_alpha_r0);
    let exp_alpha_r0_sq = (1.0 - exp_alpha_r0).powi(2);
    SimulationState {
        time: 0.0,
        displacement: r0_a0,
        force: init_force,
        acceleration: init_force / m_au,
        velocity: 0.0,
        kinetic_e: 0.0,
        potential_e: d_au * exp_alpha_r0_sq,
        total_e: d_au * exp_alpha_r0_sq,
    }
}

fn simulate_harmonic_osc_state(state: &mut SimulationState, k_au: f32, m_au: f32) {
    state.time += DT_AU;
    let old_acceleration = state.acceleration;
    state.displacement += state.velocity * DT_AU + 0.5 * old_acceleration * DT_AU * DT_AU;
    state.force = -k_au * state.displacement;
    state.acceleration = state.force / m_au;
    state.velocity += 0.5 * (old_acceleration + state.acceleration) * DT_AU;
    state.kinetic_e = 0.5 * m_au * state.velocity.powi(2);
    state.potential_e = 0.5 * k_au * state.displacement.powi(2);
    state.total_e = state.kinetic_e + state.potential_e;
}

fn simulate_morse_potential_state(state: &mut SimulationState, d_au: f32, alpha_au: f32, m_au: f32) {
    state.time += DT_AU;
    let old_acceleration = state.acceleration;
    state.displacement += state.velocity * DT_AU + 0.5 * old_acceleration * DT_AU * DT_AU;
    let exp_alpha_r0 = f32::exp(-alpha_au * state.displacement);
    let next_force = -2.0 * d_au * alpha_au * exp_alpha_r0 * (1.0 - exp_alpha_r0);
    let exp_alpha_r0_sq = (1.0 - exp_alpha_r0).powi(2);
    state.force = next_force;
    state.acceleration = next_force / m_au;
    state.velocity += 0.5 * (old_acceleration + state.acceleration) * DT_AU;
    state.kinetic_e = 0.5 * m_au * state.velocity.powi(2);
    state.potential_e = d_au * exp_alpha_r0_sq;
    state.total_e = state.kinetic_e + state.potential_e;
}

fn calculate_harmonic_osc_state(state: &mut SimulationState, r0_a0: f32, k_au: f32, m_au: f32) {
    state.time += DT_AU;
    let omega = (k_au / m_au).sqrt();
    let omega_time_cos = (omega * state.time).cos();
    let omega_time_sin = (omega * state.time).sin();
    state.displacement = r0_a0 * omega_time_cos;
    state.force = -k_au * state.displacement;
    state.acceleration = -r0_a0 * omega * omega * omega_time_cos;
    state.velocity = -r0_a0 * omega * omega_time_sin;
    state.kinetic_e = 0.5 * k_au * r0_a0 * r0_a0 * omega_time_sin * omega_time_sin;
    state.potential_e = 0.5 * k_au * r0_a0 * r0_a0 * omega_time_cos * omega_time_cos;
    state.total_e = state.kinetic_e + state.potential_e;
}

fn init_lennard_jones(m_au: f32, rstar: f32, eps: f32, r0_lj: f32) -> SimulationState {
    let rstar_over = rstar / (r0_lj + rstar);
    let init_force = (12.0 / (r0_lj + rstar)) * eps * (rstar_over.powi(12) - rstar_over.powi(6));
    SimulationState {
        time: 0.0,
        displacement: r0_lj,
        force: init_force,
        acceleration: init_force / m_au,
        velocity: 0.0,
        kinetic_e: 0.0,
        potential_e: eps * (rstar_over.powi(12) - 2.0 * rstar_over.powi(6) + 1.0),
        total_e: eps * (rstar_over.powi(12) - 2.0 * rstar_over.powi(6) + 1.0),
    }
}

fn simulate_lennard_jones_state(state: &mut SimulationState, m_au: f32, rstar: f32, eps: f32) {
    state.time += DT_AU;
    let old_acceleration = state.acceleration;
    state.displacement += state.velocity * DT_AU + 0.5 * old_acceleration * DT_AU * DT_AU;
    let rstar_over = rstar / (state.displacement + rstar);
    state.force  = (12.0 / (state.displacement + rstar)) * eps * (rstar_over.powi(12) - rstar_over.powi(6));
    state.acceleration = state.force / m_au;
    state.velocity += 0.5 * (old_acceleration + state.acceleration) * DT_AU;
    state.kinetic_e = 0.5 * m_au * state.velocity.powi(2);
    state.potential_e = eps * (rstar_over.powi(12) - 2.0 * rstar_over.powi(6) + 1.0);
    state.total_e = state.kinetic_e + state.potential_e;
}


fn main() {

    // Set Up //////////////////////////////////
    // Generate experiment names (output directories)
    let exp_name_base = format!("dt{}_temp{}", DT_AU as i32, TEMP_K as i32);
    let exp_name_harm_sim = format!("harmonic_osc_{}_sim", exp_name_base);
    let exp_name_harm_calc = format!("harmonic_osc_{}_calc", exp_name_base);
    let exp_name_harm_diff = format!("harmonic_osc_{}_diff", exp_name_base);
    let exp_name_morse_sim = format!("morse_potential_{}_sim", exp_name_base);
    let exp_name_morse_overlay = format!("morse_potential_{}_overlay", exp_name_base);
    let exp_name_hg2_lj = format!("hg2_lj_{}_sim", exp_name_base);
    let exp_name_ar2_lj = format!("ar2_lj_{}_sim", exp_name_base);

    println!("Calculation Initial Displacements...");
    // Calculate the initial displacement in SI
    let r0_si_harm_h2: f32 = ((2.0 * KB * TEMP_K) / H2_K_SI).sqrt();
    let r0_si_morse: f32 = (1.0 - (H2_K_SI * r0_si_harm_h2 * r0_si_harm_h2 / (2.0 * H2_D_SI)).sqrt()).ln() / (-H2_ALPHA_SI);
    let r0_si_harm_hg2: f32 = ((2.0 * KB * TEMP_K) / HG2_K_SI).sqrt();
    let r0_si_harm_ar2: f32 = ((2.0 * KB * TEMP_K) / AR2_K_SI).sqrt();
    // Convert the initial displacement to AU
    let r0_a0_harm_h2: f32 = r0_si_harm_h2 / A0_TO_M;
    let r0_a0_morse_h2: f32 = r0_si_morse / A0_TO_M;
    let r0_a0_harm_hg2: f32 = r0_si_harm_hg2 / A0_TO_M;
    let r0_a0_harm_ar2: f32 = r0_si_harm_ar2 / A0_TO_M;
    // let r0_a0_lj_hg2: f32 = HG2_RSTR_AU * (((2.0 * HG2_EPS_AU).powf(1.0 / 12.0) * ((HG2_K_AU).sqrt() * r0_a0_harm_hg2 + (2.0 * HG2_EPS_AU).sqrt()).powf(-1.0 / 6.0)) - 1.0);
    // let r0_a0_lj_ar2: f32 = AR2_RSTR_AU * (((2.0 * AR2_EPS_AU).powf(1.0 / 12.0) * ((AR2_K_AU).sqrt() * r0_a0_harm_ar2 + (2.0 * AR2_EPS_AU).sqrt()).powf(-1.0 / 6.0)) - 1.0);
    let r0_a0_lj_hg2 = HG2_RSTR_AU * (((1.0 + (HG2_K_AU * r0_a0_harm_hg2.powi(2) / HG2_EPS_AU).sqrt()).powf(-1.0 / 6.0)) - 1.0);
    let r0_a0_lj_ar2 = AR2_RSTR_AU * (((1.0 + (AR2_K_AU * r0_a0_harm_ar2.powi(2) / AR2_EPS_AU).sqrt()).powf(-1.0 / 6.0)) - 1.0);

    println!("r0_si_harm_hg2: {}", r0_si_harm_hg2);
    println!("r0_si_harm_ar2: {}", r0_si_harm_ar2);
    println!("r0_a0_harm_hg2: {}", r0_a0_harm_hg2);
    println!("r0_a0_harm_ar2: {}", r0_a0_harm_ar2);
    println!("r0_a0_lj_hg2: {}", r0_a0_lj_hg2);
    println!("r0_a0_lj_ar2: {}", r0_a0_lj_ar2);

    println!("\nStarting Simulation...");


    // H2 Harmonic Oscillator Simulation //////////////////////////////////
    // Initialize the simulation system
    let mut state = init_harmonic_osc(r0_a0_harm_h2, H2_K_AU, H2_M_AU);
    // Create a vector to store the history of states
    let mut harm_sim_history: Vec<SimulationState> = Vec::new();
    // For N_STEPS steps
    for sim_iter in 0..N_STEPS {
        // Simulate the next state
        simulate_harmonic_osc_state(&mut state, H2_K_AU, H2_M_AU);
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


    // H2 Harmonic Oscillator Calculation //////////////////////////////////
    // Initialize the predicted simulation system
    let mut state = init_harmonic_osc(r0_a0_harm_h2, H2_K_AU, H2_M_AU);
    // Create a vector to store the history of states
    let mut harm_calc_history: Vec<SimulationState> = Vec::new();
    // For N_STEPS steps
    for sim_iter in 0..N_STEPS {
        // Calculate the next state
        calculate_harmonic_osc_state(&mut state, r0_a0_harm_h2, H2_K_AU, H2_M_AU);
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



    // H2 Harmonic Oscillator Difference //////////////////////////////////
    let difference_of_harm_histories = get_difference_of_histories(&harm_sim_history, &harm_calc_history);
    plot_histories(&difference_of_harm_histories, &exp_name_harm_diff, "png");
    write_history(&difference_of_harm_histories, &exp_name_harm_diff);



    // H2 Morse Potential Simulation //////////////////////////////////
    // Initialize the simulation system
    let mut state = init_morse_potential(r0_a0_morse_h2, H2_D_AU, H2_ALPHA_AU, H2_M_AU);
    // Create a vector to store the history of states
    let mut morse_sim_history: Vec<SimulationState> = Vec::new();
    // For N_STEPS steps
    for sim_iter in 0..N_STEPS {
        // Simulate the next state
        simulate_morse_potential_state(&mut state, H2_D_AU, H2_ALPHA_AU, H2_M_AU);
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



    // Hg2 Lennard-Jones Simulation //////////////////////////////////
    println!("\nStarting Hg2 Lennard–Jones Simulation...");
    let mut state = init_lennard_jones(HG2_M_AU, HG2_RSTR_AU, HG2_EPS_AU, r0_a0_lj_hg2);
    let mut hg2_lj_history: Vec<SimulationState> = Vec::new();
    for sim_iter in 0..N_STEPS {
        simulate_lennard_jones_state(&mut state, HG2_M_AU, HG2_RSTR_AU, HG2_EPS_AU);
        hg2_lj_history.push(state.clone());
        if sim_iter % PRINT_FREQ == 0 {
            print_state(&state);
        }
    }
    plot_histories(&hg2_lj_history, &exp_name_hg2_lj, "png");
    write_history(&hg2_lj_history, &exp_name_hg2_lj);



    // Ar2 Lennard-Jones Simulation //////////////////////////////////
    println!("\nStarting Ar2 Lennard–Jones Simulation...");
    let mut state = init_lennard_jones(AR2_M_AU, AR2_RSTR_AU, AR2_EPS_AU, r0_a0_lj_ar2);
    let mut ar2_lj_history: Vec<SimulationState> = Vec::new();
    for sim_iter in 0..N_STEPS {
        simulate_lennard_jones_state(&mut state, AR2_M_AU, AR2_RSTR_AU, AR2_EPS_AU);
        ar2_lj_history.push(state.clone());
        if sim_iter % PRINT_FREQ == 0 {
            print_state(&state);
        }
    }
    plot_histories(&ar2_lj_history, &exp_name_ar2_lj, "png");
    write_history(&ar2_lj_history, &exp_name_ar2_lj);

    println!("Done!");
}


