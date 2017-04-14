#ifndef flecsi_sp_pic_species_h
#define flecsi_sp_pic_species_h

namespace flecsi {
    namespace sp {
        namespace pic {

            template <class real_t> class species_
            {
                public: 

                    // Species properties
                    real_t q;
                    real_t m;

                    size_t num_particles = 0;

                    species_(real_t q_i, real_t m_i, size_t num_particles_i = 0)
                    {
                        q = q_i;
                        m = m_i;
                        num_particles = num_particles_i;
                    }
            };

        } // namespace pic
    } // namespace sp
} // namespace flecsi

#endif // flecsi_sp_pic_species_h
