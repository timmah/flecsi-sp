#ifndef flecsi_sp_pic_initial_conditions_h
#define flecsi_sp_pic_initial_conditions_h

#include <iostream>
#include <fstream>

#include <flecsi-sp/pic/types.h>
#include <flecsi-sp/pic/helpers.h>
#include <flecsi-sp/pic/species.h>
#include <flecsi-sp/pic/particles.h>

namespace flecsi {
    namespace sp {
        namespace pic {

            // Type imports
            using species_t = types::species_t;
            using particle_list_t = types::particle_list_t;

            // What's the best way to implement this?
            // Want something which can be trivially user extendible 
            class initial_conditions_t {


                public:
                    virtual inline int common_function() = 0;

                // Prevent instantiation of base class
                protected:
                    size_t width = 0;
                    size_t height = 0;
                    size_t nppc = 0;
                    // TODO: Make construtor

                    initial_conditions_t () {}
                    initial_conditions_t (const initial_conditions_t&) = delete;
                    initial_conditions_t & operator = (const initial_conditions_t&) = delete;
            };

            class two_stream : public initial_conditions_t
            {

                inline int common_function() final
                {
                    return 1;
                }

                void process_cell(mesh_t& m, species_t& sp, auto particles_accesor)
                {
                    for ( auto c : m.cells() )
                    {
                        auto v = m.vertices(c)[0]; // Try and grab the bottom corner of this cell
                        auto coord = v->coordinates();

                        auto& cell_particles = particles_accesor[c];

                        size_t width = initial_conditions_t::width;
                        size_t height = initial_conditions_t::height;

                        size_t cell_id = coords_to_1d(coord, width, height);

                        for (size_t i = 0; i < nppc; i++)
                        {
                            add_particle(sp, cell_particles, cell_id, nppc);
                        }
                    }
                }

                void add_particle(
                        species_t& sp,
                        particle_list_t& cell_particles,
                        size_t cell_id,
                        size_t nppc
                )
                {

                    real_t x = 1; //random_real( x_min, x_max );
                    real_t y = 1; //random_real( y_min, y_max );
                    real_t z = 1; //random_real( z_min, z_max );
                    // TODO: Set

                    std::array<real_t,3> velocity = sp.initial_velocity;

                    real_t ux = velocity[0];
                    real_t uy = velocity[1];
                    real_t uz = velocity[2];

                    real_t w = sp.m;

                    cell_particles.add_particle(x, y, z, cell_id, ux, uy, uz, w);
                }
            };

        } // namespace pic
    } // namespace sp
} // namespace flecsi

#endif // initial_conditions
