use rand::Rng;
mod simulation;
use simulation::{Simulation, Atom};

fn main() {
    let box_size = 5.0;
    let num_atoms = 11; // 10 carbon atoms + 1 lithium atom
    let dt = 0.01;  // Time step
    let cutoff_radius = 3.5; // Cutoff radius for neighbor list

    // Create lithium atom at the center of the box
    let mut atoms = Vec::new();
    atoms.push(Atom::new(box_size / 2.0, box_size / 2.0, box_size / 2.0));

    // Create carbon atoms with random positions
    let mut rng = rand::thread_rng();
    for _ in 0..num_atoms-1 {
        let x = rng.gen_range(0.0..box_size);
        let y = rng.gen_range(0.0..box_size);
        let z = rng.gen_range(0.0..box_size);
        atoms.push(Atom::new(x, y, z));
    }

    // Initialize the simulation
    let mut sim = Simulation::new(atoms, box_size, dt, cutoff_radius);

    // Run the simulation
    for _ in 0..100 {
        sim.verlet_step();
        sim.print_state();
    }
}
