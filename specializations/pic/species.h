#ifndef flecsi_sp_pic_species_h
#define flecsi_sp_pic_species_h

#include "flecsi-sp/pic/entity_types.h"

namespace flecsi {
    namespace sp {
        namespace pic {

            static constexpr size_t num_dimensions = pic_config_t::num_dimensions;

            template <class real_t> class species_
            {
                public:

                    // Species properties
                    real_t q;
                    real_t m;
                    std::array<real_t, num_dimensions> initial_velocity;
                    size_t key;

                    size_t num_particles = 0;

                    species_(real_t q_i, real_t m_i, size_t num_particles_i = 0)
                    {
                        q = q_i;
                        m = m_i;
                        num_particles = num_particles_i;
                    }

                    void set_initial_velocity(real_t x, real_t y, real_t z)
                    {
                        initial_velocity = { x, y, z };
                    }
            };

        } // namespace pic
    } // namespace sp
} // namespace flecsi

#endif // flecsi_sp_pic_species_h
