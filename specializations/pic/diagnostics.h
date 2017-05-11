#ifndef flecsi_sp_pic_diagnostics_h
#define flecsi_sp_pic_diagnostics_h

namespace flecsi {
    namespace sp {
        namespace pic {
            namespace pic {
                namespace diagnostics {

                    double calculate_kentic_energy(double part_ux, double part_uy, double part_uz, double part_m) {

                        //#ifdef PER_PARTICLE_WEIGHT
                        //double mass = part_weight * part_m;
                        //#else
                        double mass = part_m;
                        //#endif

                        double gamma = sqrt(part_ux * part_ux + part_uy * part_uy + part_uz * part_uz + 1.0);
                        double kinetic_energy = mass * (gamma - 1.0);

                        return kinetic_energy;

                    }

                    double particle_energy(mesh_t& m, species_t& sp)
                    {
                        double total_particle_energy = 0.0;

                        auto particles_accesor = get_particle_accessor(m, sp.key);

                        for ( auto c : m.cells() ) {

                            auto& cell_particles = particles_accesor[c];

                            for (size_t i = 0; i < cell_particles.block_number+1; i++)
                            {

                                for (int v = 0; v < PARTICLE_BLOCK_SIZE; v++)
                                {
                                    real_t ux = cell_particles.get_ux(i, v);
                                    real_t uy = cell_particles.get_uy(i, v);
                                    real_t uz = cell_particles.get_uz(i, v);

                                    total_particle_energy += calculate_kentic_energy(ux,uy,uz, sp.m);
                                }
                            }
                        }
                        return total_particle_energy;
                    }

                    // Function to calculate the total energy based on fields and particles
                    void energy_diagnostic()
                    {
                        // In vpic this calls energy_f
                        // And then calls energy_p for each species
                    }

                }
        }
    }
}

#endif  // guard
