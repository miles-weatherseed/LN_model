mod functions;
use rand::Rng;
use std::env;

// Declaring the constants for the LN model in the global scope. Defining as floats will help with SA later on
const DEFAULT_RADIUS: f64 = 500.0;
const DEFAULT_CONTACT_RADIUS: f64 = 20.0;
const T_CELL_DENSITY: f64 = 1e6;
const DEFAULT_AGSPEC_FREQ: f64 = 1e-5;
const TOTAL_TNUM: f64 = DEFAULT_AGSPEC_FREQ*T_CELL_DENSITY*4.0*0.33333*3.1415926*DEFAULT_RADIUS*DEFAULT_RADIUS*DEFAULT_RADIUS*1e-9;
const T_VELOCITY_MEAN: f64 = 10.06354;
const T_VELOCITY_STDEV: f64 = 0.5;
const _T_GAMMA_SHAPE: f64 = 2.28567314561;
const _T_GAMMA_SCALE: f64 = 4.40287799651;
const T_FREE_PATH_MEAN: f64 = 25.0;
const T_FREE_PATH_STDEV: f64 = 3.0;
const NUM_ANTIGEN_IN_CONTACT_AREA: u32 = 500;
const NUM_ANTIGEN_ON_DC: u32 = 100000;
const TCELL_ACTIVATION_THRESHOLD: u32 = 20;
const FIRST_DC_ARRIVAL: f64 = 1080.0;
const DC_ARRIVAL_DURATION: f64 = 5.0;//360.0;
const _NUMBER_DCS: u32 = 720;
const DC_FLUX: f64 = 0.03333; // per minute
const DEFAULT_ANTIGEN_DECAY_RATE: f64 = 0.001;
const DEFAULT_COGNATE_RATIO_MEAN_AT_DERMIS: f64 = 0.2;
const _ANTIGEN_OOB_TOLERANCE: f64 = 0.1*_NUMBER_DCS as f64;
const _POSITION_OOB_TOLERANCE: f64 = 100.0*_NUMBER_DCS as f64;
const MINIMUM_CONTACT_GRID_MEMORY_REDUCTION_FACTOR: f64 = 0.25;
const _TIME_TOLERANCE: f64 = 1e-6;

const NUM_REPEATS: u16 = 20;
const DEFAULT_TIMESTEP: f64 = 0.05;
const NUM_TIME_MEASUREMENTS: u32 = 500;
const NUM_TIMESTEPS_2DAY: f64 = 2880.0 / DEFAULT_TIMESTEP;
const NUM_TIMESTEPS_28DAY: f64 = 28.0 * 1440.0 / DEFAULT_TIMESTEP;

// cross-presentation parameters


// the main function which coordinates all the activities of the model
fn main() {
    let args: Vec<String> = env::args().collect();
    let ag_on_arrival:f64 = args[1].parse().unwrap();
    let off_rate:f64 = args[2].parse().unwrap();
    let contact_radius_sq:f64 = args[3].parse().unwrap();
    let dc_flux: f64 = args[4].parse().unwrap();
    let t_cell_velocity: f64 = args[5].parse().unwrap();
    let t_free_path: f64 = args[6].parse().unwrap();
    let t_act_threshold: f64 = args[7].parse().unwrap();

    return simulation(NUM_TIMESTEPS_28DAY.round() as u32, DEFAULT_TIMESTEP, dc_flux, DEFAULT_RADIUS, contact_radius_sq.sqrt(), TOTAL_TNUM.round() as u32, ag_on_arrival, off_rate, t_cell_velocity, t_free_path, t_act_threshold.round() as u32)
}

fn simulation(num_time_steps: u32, time_step: f64, dc_flux:f64, radius: f64, contact_radius: f64, num_t_cells: u32, ag_on_arrival:f64, off_rate:f64, t_cell_velocity:f64, t_free_path:f64, t_act_threshold:u32) {

    // creating the variables we need
    let mut total_num_activated: u32 = 0;
    let mut _total_num_unique_interactions: u32 = 0;
    let mut total_num_interactions: u32 = 0;
    let cog_ag_on_arrival: f64 = ag_on_arrival;
    let num_d_cells: u32 = (2880.0*dc_flux) as u32;
    
    // we must first make containers of dCells and tCells
    let mut d_cell_list = Vec::new();
    let mut t_cell_list = Vec::new();
    for i in 0..num_d_cells {
        let some = functions::functions::DCell {x: -1.0, y : -1.0, z : -1.0, cog_ag_ratio: 0.0, prob_activation: 0.0, time_antigen_count_last_updated: 0.0};
        d_cell_list.push(some);
    }
    for i in 0..num_t_cells {
        let some = functions::functions::TCell {x: -1.0, y : -1.0, z : -1.0, vx: 0.0, vy: 0.0, vz: 0.0, failed_dc_contact: vec![-1], interactions: 0};
        t_cell_list.push(some);
    }
    
    // create a discrete grid to search space for T cells looking for DCs. Make this grid iteratively more coarse until the memory usage is less than 1GB
    let mut occupied_positions_array_size = 100000000;
    let mut contract_grid_reduction_factor = 2.0;
    let mut num_positions: u32 = 0;

    while (occupied_positions_array_size*10*4) as f64 > 2e9 {
        if contract_grid_reduction_factor < MINIMUM_CONTACT_GRID_MEMORY_REDUCTION_FACTOR{
            return;
        }
        else {
            contract_grid_reduction_factor /= 2.0;
        }
        num_positions = (2.0*contract_grid_reduction_factor*radius / contact_radius+1.0) as u32;
        occupied_positions_array_size = num_positions*num_positions*num_positions;
    }
    let cell_side: f64 = contact_radius/contract_grid_reduction_factor;

    // allocating memory
    let mut occupied_positions: Vec<Vec<u32>> = Vec::with_capacity(occupied_positions_array_size as usize);
    for i in 0..occupied_positions_array_size {
        occupied_positions.insert(i as usize, Vec::new())
    }
    assert_eq!(occupied_positions.len(), occupied_positions_array_size as usize);

    let mut preoccupied_positions: Vec<Vec<u32>> = Vec::with_capacity(occupied_positions_array_size as usize);
    for i in 0..occupied_positions_array_size {
        preoccupied_positions.insert(i as usize, Vec::new())
    }
    assert_eq!(occupied_positions.len(), occupied_positions_array_size as usize);

    // more containers for T cells
    let mut cell_movement_order: Vec<u32> = Vec::with_capacity(num_t_cells as usize);
    let mut free_path_remaining: Vec<f64> = vec![0.0; num_t_cells as usize];
    let _time_steps_between_dc_arrival: f64 = (DC_ARRIVAL_DURATION/time_step)/(num_d_cells - 1) as f64;

    for repeat in 0..NUM_REPEATS {
        // reinitialise number of activated t-cells and the positions of DCs
        let mut d_cells_present: u32 = 0;
        let mut this_num_activated: u32 = 0;
        let mut _this_num_unique_interactions: u32 = 0;
        let mut this_num_interactions: u32 = 0;
        let mut d_travelled_x = vec![0.0; num_t_cells as usize];
        let mut d_travelled_y = vec![0.0; num_t_cells as usize];
        let mut d_travelled_z = vec![0.0; num_t_cells as usize];

        if repeat > 0 {
            for i in 0..occupied_positions_array_size {
                occupied_positions[i as usize].clear();
                preoccupied_positions[i as usize].clear();
            }
        };
        cell_movement_order.resize(num_t_cells as usize, 0);
        for i in 0..num_t_cells {cell_movement_order[i as usize] = i};

        // place D Cells on lattice
        for i in 0..num_d_cells {
            let posn  = functions::functions::place_d_cells(i, radius, cell_side, num_positions, &mut preoccupied_positions, &d_cell_list, contact_radius);
            d_cell_list[i as usize].x = posn.0;
            d_cell_list[i as usize].y = posn.1;
            d_cell_list[i as usize].z = posn.2;
            d_cell_list[i as usize].cog_ag_ratio = cog_ag_on_arrival;
        }
        // place T Cells on a lattice
        for i in 0..num_t_cells {
            let (xnew, ynew, znew)  = functions::functions::place_t_cells(radius, cell_side, num_positions, &preoccupied_positions, &d_cell_list, contact_radius);
            t_cell_list[i as usize].x = xnew;
            t_cell_list[i as usize].y = ynew;
            t_cell_list[i as usize].z = znew;
            t_cell_list[i as usize].failed_dc_contact = vec![-1];
            free_path_remaining[i as usize] = 0.0;
        }

        for t in 0..num_time_steps {

            //for i in 0..occupied_positions_array_size {occupied_positions[i as usize].clear()};
            
            // spawn new DCs
            // while t as f64*time_step < DC_ARRIVAL_DURATION && t as f64>= d_cells_present as f64*time_steps_betwee_dc_arrival{
            //     let d_cell_num: u32 = d_cells_present;
            //     d_cell_list[d_cell_num as usize].time_antigen_count_last_updated = 0.0;
            //     d_cell_list[d_cell_num as usize].cog_ag_ratio = cog_ag_on_arrival;
            //     functions::functions::place_dc_on_discrete_grid(d_cell_list[d_cell_num as usize].x, d_cell_list[d_cell_num as usize].y, d_cell_list[d_cell_num as usize].z, d_cell_num, radius, cell_side, &mut occupied_positions, num_positions);
            //     d_cells_present += 1;
            // }

            // spawning function for flux model
            while t as f64*time_step*dc_flux >= d_cells_present as f64 { // we haven't got as many DCs as we should, better add one   
            let dloc = d_cells_present % num_d_cells; // dloc is the position in the d_cell_list we will be adding/removing
            if d_cells_present > num_d_cells {
                functions::functions::remove_dc_from_discrete_grid(d_cell_list[dloc as usize].x, d_cell_list[dloc as usize].y, d_cell_list[dloc as usize].z, dloc, radius, cell_side, &mut occupied_positions, num_positions); // get rid of old DC that was here!
                let posn  = functions::functions::place_d_cells(dloc, radius, cell_side, num_positions, &mut occupied_positions, &d_cell_list, contact_radius);
                d_cell_list[dloc as usize].x = posn.0;
                d_cell_list[dloc as usize].y = posn.1;
                d_cell_list[dloc as usize].z = posn.2;
            }
            d_cell_list[dloc as usize].time_antigen_count_last_updated = t as f64 * time_step;
            d_cell_list[dloc as usize].cog_ag_ratio = cog_ag_on_arrival;
            functions::functions::place_dc_on_discrete_grid(d_cell_list[dloc as usize].x, d_cell_list[dloc as usize].y, d_cell_list[dloc as usize].z, dloc, radius, cell_side, &mut occupied_positions, num_positions);
            d_cells_present += 1;
            }

            // move our t cells and apply boundary conditions
            let mut pull_backs = 0;
            //if t%500==0{println!("{}, {}, {}", cell_movement_order[0], cell_movement_order[1], cell_movement_order[2])};
            for cell in 0..cell_movement_order.len() {
                let t_cell_num = cell_movement_order[cell - pull_backs];
                let mut x: f64 = t_cell_list[t_cell_num as usize].x; let mut y: f64 = t_cell_list[t_cell_num as usize].y; let mut z: f64 = t_cell_list[t_cell_num as usize].z;
                let vx: f64 = t_cell_list[t_cell_num as usize].vx; let vy: f64 = t_cell_list[t_cell_num as usize].vy; let vz: f64 = t_cell_list[t_cell_num as usize].vz;
                let dx: f64 = vx*time_step; let dy: f64 = vy*time_step; let dz: f64 = vz*time_step;
                d_travelled_x[t_cell_num as usize] += (dx*dx).sqrt();
                d_travelled_z[t_cell_num as usize] += (dz*dz).sqrt();
                d_travelled_y[t_cell_num as usize] += (dy*dy).sqrt();
                let (mut newx, mut newy, mut newz) = (x + dx, y + dy, z + dz);
                free_path_remaining[t_cell_num as usize] -= functions::functions::magnitude(dx, dy, dz);
                
                // see if we've left sphere or reached end of mean free path
                if functions::functions::outside_sphere(newx, newy, newz, radius) {
                    // regenerate velocity to move t cell away from sphere surface
                    let vmag = functions::functions::genrand_t_velocity(t_cell_velocity, T_VELOCITY_STDEV);
                    let theta = functions::functions::genrand_uniform_pos().acos();
                    let phi = 6.283185307*functions::functions::genrand_uniform_pos();
                    let magsq = functions::functions::magnitude(x, y, z);
                    x /= magsq; y /= magsq; z /= magsq;
                    let (v0, v1, v2): (f64, f64, f64);
                    v0 = -vmag*x; v1 = -vmag*y; v2 = -vmag*z;
                    let (v0, v1, v2) = functions::functions::arbitrary_axis_rotation(z, y, -x, v0, v1, v2, theta);
                    let (v0, v1, v2) = functions::functions::arbitrary_axis_rotation(x, y, z, v0, v1, v2, phi);
                    //println!("{}", functions::functions::magnitude(v0, v1, v2));
                    t_cell_list[t_cell_num as usize].vx = v0; t_cell_list[t_cell_num as usize].vy = v1; t_cell_list[t_cell_num as usize].vz = v2;
                    free_path_remaining[t_cell_num as usize] = functions::functions::genrand_t_freepath(t_free_path, T_FREE_PATH_STDEV);
            
                    // finally, reverse previous step we tried to make and take a new one
                    newx += v0*time_step - dx; newy += v1*time_step - dy; newz += v2*time_step - dz;
                }
                else if free_path_remaining[t_cell_num as usize] <= 0.0 {
                    let vmag = functions::functions::genrand_t_velocity(t_cell_velocity, T_VELOCITY_STDEV);
                    let theta = functions::functions::genrand_uniform_pos().acos();
                    let phi = 6.283185307*functions::functions::genrand_uniform_pos();
                    t_cell_list[t_cell_num as usize].vx = vmag*theta.sin()*phi.cos(); 
                    t_cell_list[t_cell_num as usize].vy = vmag*theta.sin()*phi.sin(); 
                    t_cell_list[t_cell_num as usize].vz = vmag*theta.cos();
                    //println!("{}", functions::functions::magnitude(t_cell_list[t_cell_num as usize].vx, t_cell_list[t_cell_num as usize].vy, t_cell_list[t_cell_num as usize].vz));
                    //println!("{}, {}, {}", t_cell_list[t_cell_num as usize].vx, t_cell_list[t_cell_num as usize].vy, t_cell_list[t_cell_num as usize].vz);
                    free_path_remaining[t_cell_num as usize] = functions::functions::genrand_t_freepath(t_free_path, T_FREE_PATH_STDEV);
                }

                // check for DCs nearby to activate t cell and either update position or remove from simulation
                let (d_cell_num, in_contact) = functions::functions::check_contact_with_dendrites_t(newx, newy, newz, radius, cell_side, num_positions, &occupied_positions, &d_cell_list, contact_radius, &t_cell_list[t_cell_num as usize].failed_dc_contact);
                if in_contact {
                    d_cell_list[d_cell_num as usize].cog_ag_ratio *= f64::exp(-((t as f64 * time_step) - d_cell_list[d_cell_num as usize].time_antigen_count_last_updated) * off_rate);
                    let (d_coord_x, d_coord_y, d_coord_z) = functions::functions::set_coordinates(d_cell_list[d_cell_num as usize].x, d_cell_list[d_cell_num as usize].y, d_cell_list[d_cell_num as usize].z, radius, cell_side);
                    //println!("{}, {}", t_cell_num, d_cell_num);
                    this_num_interactions += 1;
                    //println!("Time is {}, Ratio is {}", t as f64*time_step, d_cell_list[d_cell_num as usize].cog_ag_ratio);
                    if rand::random::<f64>() < functions::functions::get_prob_activation(NUM_ANTIGEN_IN_CONTACT_AREA, NUM_ANTIGEN_ON_DC, t_act_threshold, d_cell_list[d_cell_num as usize].cog_ag_ratio) {
                        // activation! We no longer need to track this T cell. Erase it and pull the iterator back one
                        assert_eq!(cell_movement_order.remove(cell - pull_backs), t_cell_num);
                        pull_backs += 1; // we have to pull back the iterator so we don't skip any cells
                        // instead we choose to respawn this t cell somewhere else to keep the density in the paracortex at the same level
                        //let (xnew, ynew, znew)  = functions::functions::place_t_cells(radius, cell_side, num_positions, &occupied_positions, &d_cell_list, contact_radius);
                        //t_cell_list[t_cell_num as usize].x = xnew;
                        //t_cell_list[t_cell_num as usize].y = ynew;
                        //t_cell_list[t_cell_num as usize].z = znew;
                        //t_cell_list[t_cell_num as usize].failed_dc_contact = vec![-1];
                        //free_path_remaining[t_cell_num as usize] = 0.0;
                        this_num_activated +=1 ;
                    }
                    else {
                        // activation failed - mark this pair so that they won't repeatedly try to interact
                        t_cell_list[t_cell_num as usize].failed_dc_contact.push(d_cell_num);
                        // regenerate velocity to move t cell away from DC after interaction
                        let vmag = functions::functions::genrand_t_velocity(t_cell_velocity, T_VELOCITY_STDEV);
                        let theta = functions::functions::genrand_uniform_pos().acos();
                        let phi = 6.283185307*functions::functions::genrand_uniform_pos();
                        x = newx - d_cell_list[d_cell_num as usize].x;
                        y = newy - d_cell_list[d_cell_num as usize].y;
                        z = newz - d_cell_list[d_cell_num as usize].z;
                        let magsq = functions::functions::magnitude(x, y, z);
                        x /= magsq; y /= magsq; z /= magsq;
                        let (v0, v1, v2): (f64, f64, f64);
                        v0 = -vmag*x; v1 = -vmag*y; v2 = -vmag*z;
                        let (v0, v1, v2) = functions::functions::arbitrary_axis_rotation(z, y, -x, v0, v1, v2, theta);
                        let (v0, v1, v2) = functions::functions::arbitrary_axis_rotation(x, y, z, v0, v1, v2, phi);
                        t_cell_list[t_cell_num as usize].vx = v0; t_cell_list[t_cell_num as usize].vy = v1; t_cell_list[t_cell_num as usize].vz = v2;
                        free_path_remaining[t_cell_num as usize] = functions::functions::genrand_t_freepath(t_free_path, T_FREE_PATH_STDEV);
                    }
                }
                t_cell_list[t_cell_num as usize].x = newx; t_cell_list[t_cell_num as usize].y = newy; t_cell_list[t_cell_num as usize].z = newz; 
            } // end of t cell loop
        } // end of time loop
        total_num_activated += this_num_activated;
        total_num_interactions += this_num_interactions;
        //println!("{}, {}", this_num_interactions, this_num_activated);
    } // end of repeat loop
    println!("{}    {}", total_num_interactions, total_num_activated);
    //println!("Total number of activations: {}", total_num_activated);
} // end of function
// working perfectly with gaussian velocities as of 28/1/21

/* TO ADD:_
[x] flux of DCs arriving -> also need DCs leaving at a floating rate
[] cross-presentation model -> run once at beginning, should then have arrival values of all 5 peptides
*/