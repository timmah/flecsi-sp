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
                    /*
                    particle_() {}
                    particle_(const particle_& p) {
                        std::cout << "particle copy constructor" << std::endl;
                    }
                    ~particle_() { }
                    */

                    real_t  dx[PARTICLE_BLOCK_SIZE];
                    real_t  dy[PARTICLE_BLOCK_SIZE];
                    real_t  dz[PARTICLE_BLOCK_SIZE];
                    int32_t i[PARTICLE_BLOCK_SIZE];
                    real_t ux[PARTICLE_BLOCK_SIZE];
                    real_t uy[PARTICLE_BLOCK_SIZE];
                    real_t uz[PARTICLE_BLOCK_SIZE];
                    real_t  w[PARTICLE_BLOCK_SIZE];
                    int32_t count = 0;
                    //TODO: Need to attribute_aligned on each array
            };


            template <class real_t> class particle_list_
            {
                using particle_t = particle_<real_t>;
                public:

                    size_t block_number = 0;
                    static const size_t num_blocks = 4;

                    particle_t block[num_blocks];

                    // We don't want to do anything with complex object types
                    /*
                    particle_list_()
                    {
                        std::cout << "Default Constructor " << num_blocks << std::endl;
                    }

                    particle_list_(int num_blocks)
                    {
                        block = new particle_t[num_blocks];
                        std::cout << "Constructor " << num_blocks << std::endl;
                        this->num_blocks = num_blocks;
                    }
                    //  TODO: Add destructor
                    ~particle_list_()
                    {
                        //delete block;
                        std::cout << "Destructor" << std::endl;
                    }

                    // Copy consturctor
                    particle_list_(const particle_list_& c)
                    {
                        // TODO: Think about correct way to implement this copy Constructor
                        std::cout << "Copy Constructor" << std::endl;
                        block = new particle_t[c.num_blocks];
                        for (int i = 0; i < num_blocks; i++)
                        {
                            block[i] = c.block[i];
                        }

                        block_number = c.block_number;
                        num_blocks = c.num_blocks;
                    }

                    // Assignment contructor
                    particle_list_& operator= (const particle_list_& p)
                    {
                        std::cout << "Assignment Constructor" << std::endl;
                    }

                    */
                    // TODO: This probably shouldn't be here

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
