// TODO: Document this file 
//

#include <flecsi-sp/pic/particles.h>

#ifndef flecsi_sp_pic_boundary_h
#define flecsi_sp_pic_boundary_h


namespace flecsi {
    namespace sp {
        namespace pic {


            // TODO: Do I need to make a distinction between "internal" and "physical"?
            // Assume everything is physical for now and trust in the flecsi magic for cell data?
            enum Boundary_type {
                left, // X
                right, // X
                top, // Y
                bottom, // Y
                front, // Z 
                back // Z
            };

            // TODO: Template the class based on this?
            //
            // TODO: Passing particle_list_t AND real_t doesn't make much sense
            template <class particle_list_t, class real_t> class BoundaryStrategy
            {
                public: 
                    //using particle_list_t = particle_list_<real_t>;
                    //using particle_t = particle_<real_t>;

                    // TODO: It's not (yet) clear to me how to make this valid
                        // to avoid the re-specification of the list type (incase
                        // the real_ts don't match etc)
                    //using typename particle_list_<real_t>::particle_t; 

                    // TODO: Is there a zero-overhead way of doing this but enforcing the child to implement? 
                    // Does it matter assuming you can't have an instance of the base class?

                    // TODO: Is this trying to deal with inter_mpi boundary, or physical, or both?
                    virtual inline int process_particle(particle_list_t& p, size_t i, size_t v) = 0;

                // Prevent instantiation of base class
                protected:
                    BoundaryStrategy() {}
                    BoundaryStrategy(const BoundaryStrategy&) = delete;
                    BoundaryStrategy& operator = (const BoundaryStrategy&) = delete;
            };


            template <class particle_list_t, class real_t> 
                class ReflectiveBoundary : 
                    public BoundaryStrategy<particle_list_t, real_t>
            {
                //using typename BoundaryStrategy<real_t>::particle_t; 

                // TODO: What data interface does this need? 
                // Is just the particle information sufficient to decide where to send it next ?
                // Is this more complex than it initially seems given you may have to
                // "move()" the particle if it does something like bound off?  
                // It could take the object and directly query the data, or we could try pass the Particle to it?
                // Passing a particle would marry us to a data layout of AoS..
                
                inline int process_particle(
                        particle_list_t& p, 
                        size_t i, 
                        size_t v
                ) final
                {

                    const real_t local_cell_max = 1.0;

                    real_t x = p.get_x(i, v);
                    real_t y = p.get_y(i, v);
                    real_t z = p.get_z(i, v);

                    real_t ux = p.get_ux(i, v);
                    real_t uy = p.get_uy(i, v);
                    real_t uz = p.get_uz(i, v);
                    
                    // Invert velocity 
                    ux *= -1;
                    uy *= -1;
                    uz *= -1;

                    // Change position after bounding on boundary
                    
                    // TODO: This needs to know the boundary position and
                        // therefor may need additional data?
                    // TODO: This also needs to know how to move the particles, which may mean we dulicate code/logic
                        // We could already assume it got moved? 
                        // If {xyz} are relative, it would just be 1-that?
                    x = local_cell_max - x;
                    y = local_cell_max - y;
                    z = local_cell_max - z;

                    p.set_x(x, i, v);
                    p.set_y(y, i, v);
                    p.set_z(z, i, v);

                    p.set_ux(ux, i, v);
                    p.set_uy(uy, i, v);
                    p.set_uz(uz, i, v);
                }
            
            };

            /*
            class ReflectiveBoundary : public BoundaryStrategy
            {
                inline int process_particle(particle_t& p) final
                {
                }
            };

            class PeriodicBoundary : public BoundaryStrategy
            {
                inline int process_particle(particle_t& p) final
                {
                }
            };
            */

        } // namespace pic
    } // namespace sp
} // namespace flecsi

#endif // flecsi_sp_pic_boundary_h
