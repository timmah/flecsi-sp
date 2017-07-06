#ifndef flecsi_sp_pic_particles_h
#define flecsi_sp_pic_particles_h

namespace flecsi {
    namespace sp {
        namespace pic {

            // Let's go right for AoSoA and see if it breaks anything...
            enum { PARTICLE_BLOCK_SIZE = 16 };

            template <class real_t> class particle_
            {
                public:
                    // TODO: Need to attribute_aligned on each array
                    real_t  dx[PARTICLE_BLOCK_SIZE];
                    real_t  dy[PARTICLE_BLOCK_SIZE];
                    real_t  dz[PARTICLE_BLOCK_SIZE];
                    int32_t i[PARTICLE_BLOCK_SIZE];
                    real_t ux[PARTICLE_BLOCK_SIZE];
                    real_t uy[PARTICLE_BLOCK_SIZE];
                    real_t uz[PARTICLE_BLOCK_SIZE];
                    real_t  w[PARTICLE_BLOCK_SIZE];
                    int32_t count = 0;
            };


            template <class real_t> class particle_list_
            {
                using particle_t = particle_<real_t>;
                public:

                    size_t block_number = 0;
                    static const size_t num_blocks = 4;

                    particle_t block[num_blocks];

                    /// General Methods
                    void add_particle(
                            real_t x,
                            real_t y,
                            real_t z,
                            int ii,
                            real_t ux,
                            real_t uy,
                            real_t uz,
                            real_t w
                    )
                    {

                        // If the current block is full
                        if (block[block_number].count >= PARTICLE_BLOCK_SIZE)
                        {
                            // Move onto the next block
                            block_number++;
                        }

                        // It should never be overfull
                        assert(block[block_number].count < PARTICLE_BLOCK_SIZE);

                        int i = block_number;
                        int v = block[block_number].count;

                        set_x( i, v, x );
                        set_y( i, v, y );
                        set_z( i, v, z );

                        set_ux( i, v, ux );
                        set_uy( i, v, uy );
                        set_uz( i, v, uz );

                        set_w( i, v, w  );
                        set_i( i, v, ii );

                        block[block_number].count++;

                        // This means we over flowed our particle store...
                            // Could replace this with a refusal to store
                            // the particle?
                        //assert(block_number < num_blocks);
                    }

                    /// Setters
                    inline void set_x(size_t i, size_t v, real_t val)
                    {
                        block[i].dx[v] = val;
                    }
                    inline void set_y(size_t i, size_t v, real_t val)
                    {
                        block[i].dy[v] = val;
                    }
                    inline void set_z(size_t i, size_t v, real_t val)
                    {
                        block[i].dz[v] = val;
                    }

                    inline void set_ux(size_t i, size_t v, real_t val)
                    {
                        block[i].ux[v] = val;
                    }
                    inline void set_uy(size_t i, size_t v, real_t val)
                    {
                        block[i].uy[v] = val;
                    }
                    inline void set_uz(size_t i, size_t v, real_t val)
                    {
                        block[i].uz[v] = val;
                    }
                    inline void set_w(size_t i, size_t v, real_t val)
                    {
                        block[i].w[v] = val;
                    }
                    inline void set_i(size_t i, size_t v, int val)
                    {
                        block[i].i[v] = val;
                    }

                    /// Getters
                    inline real_t get_x(size_t i, size_t v)
                    {
                        return block[i].dx[v];
                    }
                    inline real_t get_y(size_t i, size_t v)
                    {
                        return block[i].dy[v];
                    }
                    inline real_t get_z(size_t i, size_t v)
                    {
                        return block[i].dz[v];
                    }

                    inline real_t get_ux(size_t i, size_t v)
                    {
                        return block[i].ux[v];
                    }
                    inline real_t get_uy(size_t i, size_t v)
                    {
                        return block[i].uy[v];
                    }
                    inline real_t get_uz(size_t i, size_t v)
                    {
                        return block[i].uz[v];
                    }
                    inline real_t get_w(size_t i, size_t v)
                    {
                        return block[i].w[v];
                    }

            };

        } // namespace pic
    } // namespace sp
} // namespace flecsi

#endif // flecsi_sp_pic_particles_h
