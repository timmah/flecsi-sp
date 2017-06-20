#ifndef flecsi_sp_pic_initial_conditions_h
#define flecsi_sp_pic_initial_conditions_h

#include <iostream>
#include <fstream>

#include <flecsi-sp/pic/types.h>
#include <flecsi-sp/pic/helpers.h>
#include <flecsi-sp/pic/species.h>
#include <flecsi-sp/pic/particles.h>

// TODO: this could also init species...

namespace flecsi {
    namespace sp {
        namespace pic {

            // Type imports
            using species_t = types::species_t;
            using particle_list_t = types::particle_list_t;
            using species_list_t = types::species_list_t;

            // What's the best way to implement this?
            // Want something which can be trivially user extendible
            class initial_conditions_t {

                public:
                    virtual inline int common_function() = 0;
                    virtual void process_species(
                            mesh_t& m,
                            species_t& sp
                    ) = 0;

                // Prevent instantiation of base class
                protected:
                    species_list_t& species;

                    initial_conditions_t (
                            species_list_t& species_in
                    ) :
                        species(species_in)
                    {
                        // empty
                    }

                    initial_conditions_t (const initial_conditions_t&) = delete;
                    initial_conditions_t & operator = (const initial_conditions_t&) = delete;

                    void add_particle(
                            species_t& sp,
                            particle_list_t& cell_particles,
                            size_t cell_id,
                            size_t nppc,
                            real_t x,
                            real_t y,
                            real_t z
                    )
                    {

                        std::array<real_t,3> velocity = sp.initial_velocity;

                        real_t ux = velocity[0];
                        real_t uy = velocity[1];
                        real_t uz = velocity[2];

                        real_t w = sp.m;

                        cell_particles.add_particle(x, y, z, cell_id, ux, uy, uz, w);
                        sp.num_particles++;
                    }

            };

            class two_stream : public initial_conditions_t
            {

                public:
                    two_stream(
                                species_list_t& species
                            ) : initial_conditions_t(
                                species
                            )
                    {
                        real_t q = 1.0;
                        real_t m = 1.0;

                        // Two identical species
                        species.push_back( species_t(q,m) );
                        species.push_back( species_t(q,m) );

                        // Set keys for species selector
                        species[0].key = Species_Keys::ELECTRON;
                        species[1].key = Species_Keys::NEGATIVE;

                        // Two stream
                        species[0].set_initial_velocity(0,0,1);
                        species[1].set_initial_velocity(0,0,-1);
                    }

                    inline int common_function() final
                    {
                        return 1;
                    }

                    void process_species(
                            mesh_t& m,
                            species_t& sp
                    )
                    {

                        auto particles_accesor = get_particle_accessor(m, sp.key);

                        for ( auto c : m.cells() )
                        {
                            bool proccess_this_cell = false;

                            auto& cell_particles = particles_accesor[c];
                            auto v = m.vertices(c)[0];
                            auto coords = v->coordinates();

                            // Just the local bit? So needs to be nx? not NX_global?
                            size_t width = Parameters::instance().nx;
                            size_t height = Parameters::instance().ny;

                            size_t cell_id = coords_to_1d(coords, width, height);

                            // Lets try put it a third in, going in Y

                            size_t cell_min = 0;
                            size_t cell_max = 0;

                            real_t x = random_real(0,1);
                            real_t y = random_real(0,1);
                            real_t z = random_real(0,1);

                            // Negative
                            if (sp.initial_velocity[2] == -1)
                            {
                                // Try and place in the top third row
                                cell_min = height/3.0;
                                cell_max = cell_min+(width*height);
                                y = random_real(1.0/3.0, 1.0/3.0+0.02);
                            }
                            // Positive
                            else {
                                cell_min = 2*(height/3.0);
                                cell_max = cell_min+(width*height);
                                y = random_real(2.0/3.0, 2.0/3.0+0.02);
                            }

                            if ( (cell_id > cell_min) && (cell_id < cell_max)) {
                                proccess_this_cell = true;
                            }

                            if (proccess_this_cell) {
                                size_t nppc = Parameters::instance().NPPC;
                                for (size_t i = 0; i < nppc; i++)
                                {
                                    add_particle(sp, cell_particles, cell_id, nppc, x, y, z);
                                }
                            }
                        }
                    }
            };

        } // namespace pic
    } // namespace sp
} // namespace flecsi

#endif // initial_conditions
