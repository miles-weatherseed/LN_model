use rand::Rng;

pub mod functions {
    use rand::StdRng;
    use rand::distributions::{Distribution, Gamma, Normal};
    use statrs::distribution::{Binomial, Univariate, Distribution as OtherDistribution};
    pub struct DCell {
        pub x: f64,
        pub y: f64,
        pub z: f64,
        pub cog_ag_ratio: f64,
        pub prob_activation: f64,
        pub time_antigen_count_last_updated: f64,
    }
    pub struct TCell {
        pub x: f64,
        pub y: f64,
        pub z: f64,
        pub vx: f64,
        pub vy: f64,
        pub vz: f64,
        pub failed_dc_contact: Vec<i32>,
        pub interactions: u32,
    }
    pub fn add_one(x: u32) -> u32{
        x + 1
    }
    pub fn magnitude(a: f64, b:f64, c:f64) -> f64 {
        f64::sqrt(f64::powf(a, 2.0) + f64::powf(b, 2.0) + f64::powf(c, 2.0))
    }
    pub fn sq_magnitude(a: f64, b:f64, c:f64) -> f64 {
        f64::powf(a, 2.0) + f64::powf(b, 2.0) + f64::powf(c, 2.0)
    }
    pub fn generate_coords(radius: f64) ->  (f64, f64, f64){
        (radius*(2.0*rand::random::<f64>() - 1.0), radius*(2.0*rand::random::<f64>() - 1.0), radius*(2.0*rand::random::<f64>() - 1.0))
    }
    pub fn outside_sphere(x: f64, y: f64, z: f64, radius: f64) -> bool {
        magnitude(x, y, z) > radius
    }
    pub fn set_coordinates(x: f64, y: f64, z: f64, radius: f64, cell_side: f64) -> (i32, i32, i32){
        (((x + radius)/cell_side) as i32, ((y + radius)/cell_side) as i32, ((z + radius)/cell_side) as i32)
    }
    pub fn threed_to_oned(coord_x: i32, coord_y: i32, coord_z: i32, num_positions: u32) -> i32 {
        return num_positions as i32*num_positions as i32*coord_x + num_positions as i32*coord_y + coord_z
    }
    pub fn check_contact_with_dendrites_t(x: f64, y: f64, z: f64, radius: f64, cell_side: f64, num_positions: u32, occupied_positions: &Vec<Vec<u32>>, d_cell_list: &Vec<DCell>, contact_radius: f64, failed_interaction_id: &Vec<i32>) -> (i32, bool){
        let (coord_x, coord_y, coord_z) = set_coordinates(x, y, z, radius, cell_side);
        let ref nearby_dcs = &occupied_positions[threed_to_oned(coord_x, coord_y, coord_z, num_positions) as usize];
        let mut d_cell_num: i32 = -1; let mut in_contact: bool = false;
        for &nb_dc in *nearby_dcs {
            if failed_interaction_id.iter().any(|&i| i== nb_dc as i32) {continue;}
            let dendx: f64 = d_cell_list[nb_dc as usize].x; let dendy: f64 = d_cell_list[nb_dc as usize].y; let dendz: f64 = d_cell_list[nb_dc as usize].z;
            if sq_magnitude(x - dendx, y - dendy, z - dendz) <= contact_radius*contact_radius {
                in_contact = true; d_cell_num = nb_dc as i32; break;
            }
        }
        return (d_cell_num, in_contact)
    }
    pub fn check_contact_with_dendrites_d(x: f64, y: f64, z: f64, radius: f64, cell_side: f64, num_positions: u32, occupied_positions: &Vec<Vec<u32>>, d_cell_list: &Vec<DCell>, contact_radius: f64) -> (i32, bool) {
        let (coord_x, coord_y, coord_z) = set_coordinates(x, y, z, radius, cell_side);
        let ref nearby_dcs = &occupied_positions[threed_to_oned(coord_x, coord_y, coord_z, num_positions) as usize];
        let mut d_cell_num: i32 = -1; let mut in_contact: bool = false;
        for &nb_dc in *nearby_dcs {
            let dendx: f64 = d_cell_list[nb_dc as usize].x; let dendy: f64 = d_cell_list[nb_dc as usize].y; let dendz: f64 = d_cell_list[nb_dc as usize].z;
            if sq_magnitude(x - dendx, y - dendy, z - dendz) <= contact_radius*contact_radius {
                in_contact = true; d_cell_num = nb_dc as i32; break;
            }
        }
        return (d_cell_num, in_contact)
    }
    pub fn place_t_cells(radius: f64, cell_side:f64, num_positions: u32, preoccupied_positions: &Vec<Vec<u32>>, d_cell_list: &Vec<DCell>, contact_radius: f64) -> (f64, f64, f64){
        let mut in_contact: bool = false;
        let mut d_cell_num: i32 = -1;
        let mut x: f64; let mut y: f64; let mut z: f64;
        let (x, y, z) = loop {
            let (x, y, z) = generate_coords(radius);
            if outside_sphere(x, y, z, radius) {
                continue;
            }
            let (d_cell_num, in_contact) = check_contact_with_dendrites_t(x, y, z, radius, cell_side, num_positions, preoccupied_positions, d_cell_list, contact_radius, &vec![-1]);
            if in_contact == false {break (x, y, z);}
        };
        return (x, y, z);
    }
    pub fn place_d_cells(mdci: u32, radius: f64, cell_side:f64, num_positions: u32, preoccupied_positions: &mut Vec<Vec<u32>>, d_cell_list: &Vec<DCell>, contact_radius: f64) -> (f64, f64, f64){
        let mut in_contact: bool = false;
        let mut d_cell_num: i32 = -1;
        let mut x: f64; let mut y: f64; let mut z: f64;
        let (x, y, z) = loop {
            let (x, y, z) = generate_coords(radius);
            if outside_sphere(x, y, z, radius) {
                continue;
            }
            let (d_cell_num, in_contact) = check_contact_with_dendrites_d(x, y, z, radius, cell_side, num_positions, preoccupied_positions, d_cell_list, contact_radius);
            if in_contact == false {break (x, y, z);}
        };
        place_dc_on_discrete_grid(x, y, z, mdci, radius, cell_side, preoccupied_positions, num_positions);
        return (x, y, z);
    }
    pub fn place_dc_on_discrete_grid(x: f64, y: f64, z: f64, mdci: u32, radius: f64, cell_side: f64, occupied_positions: &mut Vec<Vec<u32>>, num_positions: u32) {
        let (coord_x, coord_y, coord_z) = set_coordinates(x, y, z, radius, cell_side);
        for p in 0..2 {
            if coord_x + p - 1 >= num_positions as i32{break;}
            else if coord_x + p - 1 < 0 {continue;}
            for q in 0..2 {
                if coord_y + q - 1 >= num_positions as i32{break;}
                else if coord_y + q - 1 < 0 {continue;}
                for r in 0..2 {
                    if coord_z + r - 1 >= num_positions as i32{break;}
                    else if coord_z + r - 1 < 0 {continue;}
                    occupied_positions[threed_to_oned(coord_x + p - 1, coord_y + q - 1, coord_z + r - 1, num_positions) as usize].push(mdci);
                }
            }
        }

    }
    pub fn remove_dc_from_discrete_grid(x: f64, y: f64, z: f64, mdci: u32, radius: f64, cell_side: f64, occupied_positions: &mut Vec<Vec<u32>>, num_positions: u32) {
        let (coord_x, coord_y, coord_z) = set_coordinates(x, y, z, radius, cell_side);
        for p in 0..2 {
            if coord_x + p - 1 >= num_positions as i32{break;}
            else if coord_x + p - 1 < 0 {continue;}
            for q in 0..2 {
                if coord_y + q - 1 >= num_positions as i32{break;}
                else if coord_y + q - 1 < 0 {continue;}
                for r in 0..2 {
                    if coord_z + r - 1 >= num_positions as i32{break;}
                    else if coord_z + r - 1 < 0 {continue;}
                    let index = occupied_positions[threed_to_oned(coord_x + p - 1, coord_y + q - 1, coord_z + r - 1, num_positions) as usize].iter().position(|x| *x == mdci).unwrap();
                    occupied_positions[threed_to_oned(coord_x + p - 1, coord_y + q - 1, coord_z + r - 1, num_positions) as usize].remove(index);
                }
            }
        }

    }
    pub fn genrand_t_velocity(mean: f64, stdev: f64) -> f64 {
        //let gamma = Gamma::new(shape, scale);
        let gamma = Normal::new(mean, stdev);
        let v = gamma.sample(&mut rand::thread_rng());
        v
    }
    pub fn genrand_uniform_pos() -> f64 {
        rand::random::<f64>()
    }
    pub fn arbitrary_axis_rotation(ax_x: f64, ax_y: f64, ax_z: f64, vec_x: f64, vec_y: f64, vec_z: f64, angle: f64) -> (f64, f64, f64) {
        let (mut res_x, mut res_y, mut res_z): (f64, f64, f64);
        let mut rotation_matrix = (angle.cos() + ax_x*ax_x*(1.0 - angle.cos()), ax_x*ax_y*(1.0 - angle.cos()) - ax_z*angle.sin(), ax_x*ax_z*(1.0 - angle.cos()) + ax_y*angle.sin(),
        ax_x*ax_y*(1.0 - angle.cos()) + ax_z*angle.sin(), angle.cos() + ax_y*ax_y*(1.0 - angle.cos()), ax_y*ax_z*(1.0 - angle.cos()) - ax_x*angle.sin(),
        ax_x*ax_z*(1.0 - angle.cos()) - ax_y*angle.sin(), ax_y*ax_z*(1.0 - angle.cos()) + ax_x*angle.sin(), angle.cos() + ax_z*ax_z*(1.0 - angle.cos()));
        res_x = vec_x*&rotation_matrix.0 + vec_y*&rotation_matrix.1 + vec_z*&rotation_matrix.2;
        res_y = vec_x*&rotation_matrix.3 + vec_y*&rotation_matrix.4 + vec_z*&rotation_matrix.5;
        res_z = vec_x*&rotation_matrix.6 + vec_y*&rotation_matrix.7 + vec_z*&rotation_matrix.8;
        return (res_x, res_y, res_z);
    }
    pub fn genrand_t_freepath(mean: f64, stdev: f64) -> f64 {
        let normal = Normal::new(mean, stdev);
        let v = normal.sample(&mut rand::thread_rng());
        v
    }
    pub fn get_prob_activation(mhc_in_contact: u32, total_mhc: u32, threshold: u32, cog_ag_ratio: f64) -> f64 {
        let binom = Binomial::new(cog_ag_ratio, mhc_in_contact as u64).unwrap();
        return 1.0 - binom.cdf((threshold - 1) as f64)
    }
    pub fn add_vectors(a: Vec<f64>, b: Vec<f64>) -> Vec<f64> {
        let n = a.len() as u32;
        let mut c = vec![0.0; n as usize];
        for i in 0..n {
            c[i as usize] = a[i as usize] + b[i as usize];
        }
        println!("{:?}", c);
        c
    }
}