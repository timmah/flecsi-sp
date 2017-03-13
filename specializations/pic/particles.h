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
                public:
                    particle_<real_t>* block;
                    
                    // TODO: This really shouldn't be here
                    int block_number = 0;

                    /// General Methods
                    inline int add_particle(
                            real_t x, 
                            real_t y, 
                            real_t z, 
                            int i, 
                            real_t ux, 
                            real_t uy, 
                            real_t uz, 
                            real_t w
                    )
                    {
                        //block->count++;
                        assert(block->count < PARTICLE_BLOCK_SIZE);
                    }

                    /// Setters
                    inline void set_particle_x(size_t i, size_t v, real_t val)
                    {
                        block[i].dx[v] = val;
                    }
                    inline void set_particle_y(size_t i, size_t v, real_t val)
                    {
                        block[i].dy[v] = val;
                    }
                    inline void set_particle_z(size_t i, size_t v, real_t val)
                    {
                        block[i].dz[v] = val;
                    }

                    inline void set_particle_ux(size_t i, size_t v, real_t val)
                    {
                        block[i].ux[v] = val;
                    }
                    inline void set_particle_uy(size_t i, size_t v, real_t val)
                    {
                        block[i].uy[v] = val;
                    }
                    inline void set_particle_uz(size_t i, size_t v, real_t val)
                    {
                        block[i].uz[v] = val;
                    }

                    /// Getters
                    inline real_t get_particle_x(size_t i, size_t v)
                    {
                        return block[i].dx[v];
                    }
                    inline real_t get_particle_y(size_t i, size_t v)
                    {
                        return block[i].dy[v];
                    }
                    inline real_t get_particle_z(size_t i, size_t v)
                    {
                        return block[i].dz[v];
                    }

                    inline real_t get_particle_ux(size_t i, size_t v)
                    {
                        return block[i].ux[v];
                    }
                    inline real_t get_particle_uy(size_t i, size_t v)
                    {
                        return block[i].uy[v];
                    }
                    inline real_t get_particle_uz(size_t i, size_t v)
                    {
                        return block[i].uz[v];
                    }

            };

        } // namespace pic
    } // namespace sp
} // namespace flecsi

#endif // flecsi_sp_pic_particles_h
