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
                    size_t width = 0;
                    size_t height = 0;
                    size_t nppc = 0;
                    // TODO: Make construtor

                    initial_conditions_t (
                            size_t width_in,
                            size_t height_in,
                            size_t nppc_in
                    ) :
                        width(width_in),
                        height(height_in),
                        nppc(nppc_in)
                    {
                        // empty
                    }

                    initial_conditions_t (const initial_conditions_t&) = delete;
                    initial_conditions_t & operator = (const initial_conditions_t&) = delete;
            };

            class two_stream : public initial_conditions_t
            {

                public:
                    two_stream(
                            size_t width_in,
                            size_t height_in,
                            size_t nppc_in
                            ) : initial_conditions_t(
                                width_in,
                                height_in,
                                nppc_in
                                )
                {
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
                        //auto v = m.vertices(c)[0];
                        //auto coord = v->coordinates();

                        //auto& cell_particles = particles_accesor[c];

                        auto particles_accesor = get_particle_accessor(m, sp.key);

                        for ( auto c : m.cells() )
                        {
                            bool proccess_this_cell = false;

                            auto& cell_particles = particles_accesor[c];
                            auto v = m.vertices(c)[0];
                            auto coords = v->coordinates();

                            size_t width = initial_conditions_t::width;
                            size_t height = initial_conditions_t::height;

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
                                cell_min = height/3.0;
                                cell_max = cell_min+(width*height);
                                y = random_real(1.0/3.0, 1.0/3.0+0.02);
                            }
                            // Positive
                            else {
                                cell_min = height/3.0;
                                cell_max = cell_min+(width*height);
                                y = random_real(2.0/3.0, 2.0/3.0+0.02);
                            }


                            if ( (cell_id > cell_min) && (cell_id < cell_max)) {
                                proccess_this_cell = true;
                            }

                            if (proccess_this_cell) {

                                for (size_t i = 0; i < nppc; i++)
                                {
                                    add_particle(sp, cell_particles, cell_id, nppc, x, y, z);
                                }
                            }
                        }
                    }

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

        } // namespace pic
    } // namespace sp
} // namespace flecsi

#endif // initial_conditions
