pub struct Atom {
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub vx: f64,
    pub vy: f64,
    pub vz: f64,
    pub fx: f64,
    pub fy: f64,
    pub fz: f64,
}

impl Atom {
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Atom { x, y, z, vx: 0.0, vy: 0.0, vz: 0.0, fx: 0.0, fy: 0.0, fz: 0.0 }
    }

    pub fn distance(&self, other: &Atom, box_size: f64) -> f64 {
        let mut dx = self.x - other.x;
        let mut dy = self.y - other.y;
        let mut dz = self.z - other.z;

        // Apply minimum image convention for PBC
        if dx > box_size / 2.0 { dx -= box_size; }
        if dx < -box_size / 2.0 { dx += box_size; }
        if dy > box_size / 2.0 { dy -= box_size; }
        if dy < -box_size / 2.0 { dy += box_size; }
        if dz > box_size / 2.0 { dz -= box_size; }
        if dz < -box_size / 2.0 { dz += box_size; }

        (dx * dx + dy * dy + dz * dz).sqrt()
    }
}


pub struct Simulation {
    atoms: Vec<Atom>,
    box_size: f64,
    dt: f64,
    cutoff_radius: f64,
    neighbor_list: Vec<Vec<usize>>,
}

impl Simulation {
    pub fn new(atoms: Vec<Atom>, box_size: f64, dt: f64, cutoff_radius: f64) -> Self {
        let neighbor_list = vec![Vec::new(); atoms.len()];
        Simulation {
            atoms,
            box_size,
            dt,
            cutoff_radius,
            neighbor_list,
        }
    }

    pub fn compute_neighbor_list(&mut self) {
        let n = self.atoms.len();
        for i in 0..n {
            self.neighbor_list[i].clear();
        }

        for i in 0..n {
            for j in i + 1..n {
                let dist = self.atoms[i].distance(&self.atoms[j], self.box_size);
                if dist < self.cutoff_radius {
                    self.neighbor_list[i].push(j);
                    self.neighbor_list[j].push(i);
                }
            }
        }
    }

    pub fn compute_forces(&mut self) {
        for atom in &mut self.atoms {
            atom.fx = 0.0;
            atom.fy = 0.0;
            atom.fz = 0.0;
        }

        let n = self.atoms.len();
        for i in 0..n {
            for &j in &self.neighbor_list[i] {
                let dist = self.atoms[i].distance(&self.atoms[j], self.box_size);
                if dist == 0.0 {
                    continue;
                }
                let force = 1.0 / dist.powi(2); // Inverse square law for simplification
                let fx = force * (self.atoms[j].x - self.atoms[i].x) / dist;
                let fy = force * (self.atoms[j].y - self.atoms[i].y) / dist;
                let fz = force * (self.atoms[j].z - self.atoms[i].z) / dist;

                self.atoms[i].fx += fx;
                self.atoms[i].fy += fy;
                self.atoms[i].fz += fz;

                self.atoms[j].fx -= fx;
                self.atoms[j].fy -= fy;
                self.atoms[j].fz -= fz;
            }
        }
    }

    pub fn verlet_step(&mut self) {
        let dt2 = self.dt * self.dt;

        // Store old forces
        let old_forces: Vec<(f64, f64, f64)> = self.atoms.iter().map(|atom| (atom.fx, atom.fy, atom.fz)).collect();

        for atom in &mut self.atoms {
            // Update positions based on current velocities and forces
            atom.x += atom.vx * self.dt + 0.5 * atom.fx * dt2;
            atom.y += atom.vy * self.dt + 0.5 * atom.fy * dt2;
            atom.z += atom.vz * self.dt + 0.5 * atom.fz * dt2;

            // Apply periodic boundary conditions
            if atom.x < 0.0 { atom.x += self.box_size; }
            if atom.x >= self.box_size { atom.x -= self.box_size; }
            if atom.y < 0.0 { atom.y += self.box_size; }
            if atom.y >= self.box_size { atom.y -= self.box_size; }
            if atom.z < 0.0 { atom.z += self.box_size; }
            if atom.z >= self.box_size { atom.z -= self.box_size; }
        }

        // Compute new neighbor list
        self.compute_neighbor_list();

        // Compute new forces
        self.compute_forces();

        for (atom, old_force) in self.atoms.iter_mut().zip(old_forces.iter()) {
            // Update velocities based on average of old and new forces
            atom.vx += 0.5 * (old_force.0 + atom.fx) * self.dt;
            atom.vy += 0.5 * (old_force.1 + atom.fy) * self.dt;
            atom.vz += 0.5 * (old_force.2 + atom.fz) * self.dt;
        }
    }

    pub fn print_state(&self) {
        for (i, atom) in self.atoms.iter().enumerate() {
            println!("Atom {}: position ({:.2}, {:.2}, {:.2}), velocity ({:.2}, {:.2}, {:.2})", 
                     i, atom.x, atom.y, atom.z, atom.vx, atom.vy, atom.vz);
        }
    }
}
