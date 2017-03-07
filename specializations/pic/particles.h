namespace flecsi {
namespace sp {
namespace pic {

enum { PARTICLE_BLOCK_SIZE = 16 };

typedef struct particle
{
    float  dx[PARTICLE_BLOCK_SIZE];  
    float  dy[PARTICLE_BLOCK_SIZE];
    float  dz[PARTICLE_BLOCK_SIZE];
    int32_t i[PARTICLE_BLOCK_SIZE];  
    float ux[PARTICLE_BLOCK_SIZE];  
    float uy[PARTICLE_BLOCK_SIZE];
    float uz[PARTICLE_BLOCK_SIZE];
    float  w[PARTICLE_BLOCK_SIZE];   
    //TODO: Need to attribute_aligned on each array
} particle_t;

template <class real_t> class particle_list_
{ 
    public:

#ifdef AoSoA
        particle_t* block;
        inline real_t load_particle_x(size_t i, size_t v)
        {
            return block[i].dx[v];
        }
        inline real_t load_particle_y(size_t i, size_t v)
        {
            return block[i].dy[v];
        }
        inline real_t load_particle_z(size_t i, size_t v)
        {
            return block[i].dz[v];
        }

        inline real_t load_particle_ux(size_t i, size_t v)
        {
            return block[i].ux[v];
        }
        inline real_t load_particle_uy(size_t i, size_t v)
        {
            return block[i].uy[v];
        }
        inline real_t load_particle_uz(size_t i, size_t v)
        {
            return block[i].uz[v];
        }
#else  // AoS

        particle_t* block;
        inline real_t load_particle_x(size_t i, size_t v)
        {
            return block[i].dx;
        }
        inline real_t load_particle_y(size_t i, size_t v)
        {
            return block[i].dy;
        }
        inline real_t load_particle_z(size_t i, size_t v)
        {
            return block[i].dz;
        }

        inline real_t load_particle_ux(size_t i, size_t v)
        {
            return block[i].ux;
        }
        inline real_t load_particle_uy(size_t i, size_t v)
        {
            return block[i].uy;
        }
        inline real_t load_particle_uz(size_t i, size_t v)
        {
            return block[i].uz;
        }
#endif

};

} // namespace pic
} // namespace sp
} // namespace flecsi
