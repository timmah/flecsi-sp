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
            template <class particle_list_t> class BoundaryCondition
            {
                public:

                    real_t x_min = 0.0;
                    real_t y_min = 0.0;
                    real_t z_min = 0.0;
                    real_t x_max = 0.0;
                    real_t y_max = 0.0;
                    real_t z_max = 0.0;

                    void set_x_min(real_t val) { x_min = val; }
                    void set_y_min(real_t val) { y_min = val; }
                    void set_z_min(real_t val) { z_min = val; }

                    void set_x_max(real_t val) { x_max = val; }
                    void set_y_max(real_t val) { y_max = val; }
                    void set_z_max(real_t val) { z_max = val; }

                    real_t get_x_min() { return x_min; }
                    real_t get_y_min() { return y_min; }
                    real_t get_z_min() { return z_min; }

                    real_t get_x_max() { return x_max; }
                    real_t get_y_max() { return y_max; }
                    real_t get_z_max() { return z_max; }

                    //using particle_list_t = particle_list_<real_t>;
                    //using particle_t = particle_<real_t>;

                    // TODO: It's not (yet) clear to me how to make this valid
                    // to avoid the re-specification of the list type (incase
                    // the real_ts don't match etc)
                    //using typename particle_list_<real_t>::particle_t;

                    // TODO: Is there a zero-overhead way of doing this but enforcing the child to implement?
                    // Does it matter assuming you can't have an instance of the base class?

                    // TODO: Is this trying to deal with inter_mpi boundary, or physical, or both?
                    // TODO: These should be inline?

                    //virtual int process_particle(particle_list_t& p, size_t i, size_t v) = 0;

                    inline int process_particle(
                            particle_list_t& cell_particles,
                            size_t i,
                            size_t v
                            )
                    {
                        real_t x = cell_particles.get_x(i, v);
                        real_t y = cell_particles.get_y(i, v);
                        real_t z = cell_particles.get_z(i, v);

                        // Out by X
                        if (x > BoundaryCondition<particle_list_t>::get_x_max())
                        {
                            process_particle_x(
                                cell_particles,
                                i,
                                v,
                                BoundaryCondition<particle_list_t>::get_x_max(),
                                BoundaryCondition<particle_list_t>::get_x_min()
                            );
                        }
                        else if (x < BoundaryCondition<particle_list_t>::get_x_min()) {
                            process_particle_x(
                                cell_particles,
                                i,
                                v,
                                BoundaryCondition<particle_list_t>::get_x_min(),
                                BoundaryCondition<particle_list_t>::get_x_max()
                            );
                        }

                        // Out by Y
                        if (y > BoundaryCondition<particle_list_t>::get_y_max()) {
                            process_particle_y(
                                cell_particles,
                                i,
                                v,
                                BoundaryCondition<particle_list_t>::get_y_max(),
                                BoundaryCondition<particle_list_t>::get_y_min()
                            );
                        }
                        else if (y < BoundaryCondition<particle_list_t>::get_y_min()) {
                            process_particle_y(
                                cell_particles,
                                i,
                                v,
                                BoundaryCondition<particle_list_t>::get_y_min(),
                                BoundaryCondition<particle_list_t>::get_y_max()
                            );
                        }

                        if (z > BoundaryCondition<particle_list_t>::get_z_max()) {
                            process_particle_z(
                                cell_particles,
                                i,
                                v,
                                BoundaryCondition<particle_list_t>::get_z_max(),
                                BoundaryCondition<particle_list_t>::get_z_min()
                            );
                        }
                        else if (z < BoundaryCondition<particle_list_t>::get_z_min()) {
                            process_particle_z(
                                cell_particles,
                                i,
                                v,
                                BoundaryCondition<particle_list_t>::get_z_min(),
                                BoundaryCondition<particle_list_t>::get_z_max()
                            );
                        }
                    }
                    virtual int process_particle_x(
                            particle_list_t& p,
                            size_t i,
                            size_t v,
                            real_t extent,
                            real_t extent_far) = 0;
                    virtual int process_particle_y(
                            particle_list_t& p,
                            size_t i,
                            size_t v,
                            real_t extent,
                            real_t extent_far) = 0;
                    virtual int process_particle_z(
                            particle_list_t& p,
                            size_t i,
                            size_t v,
                            real_t extent,
                            real_t extent_far) = 0;

                    // Prevent instantiation of base class
                protected:
                    BoundaryCondition() {}
                    BoundaryCondition(const BoundaryCondition&) = delete;
                    BoundaryCondition& operator = (const BoundaryCondition&) = delete;

                    BoundaryCondition(
                            real_t x_min_in,
                            real_t x_max_in,
                            real_t y_min_in,
                            real_t y_max_in,
                            real_t z_min_in,
                            real_t z_max_in
                    ) :
                        x_min(x_min_in),
                        x_max(x_max_in),
                        y_min(y_min_in),
                        y_max(y_max_in),
                        z_min(z_min_in),
                        z_max(z_max_in)
                {
                    // empty
                }
            };


            template <class particle_list_t>
                class ReflectiveBoundary :
                    public BoundaryCondition<particle_list_t>
            {
                //using typename BoundaryConditionparticle_t;

                // TODO: What data interface does this need?
                // Is just the particle information sufficient to decide where to send it next ?
                // Is this more complex than it initially seems given you may have to
                // "move()" the particle if it does something like bound off?
                // It could take the object and directly query the data, or we could try pass the Particle to it?
                // Passing a particle would marry us to a data layout of AoS..

                public:
                    ReflectiveBoundary(
                            real_t x_min_in,
                            real_t x_max_in,
                            real_t y_min_in,
                            real_t y_max_in,
                            real_t z_min_in,
                            real_t z_max_in
                            ) :
                        BoundaryCondition<particle_list_t>(
                            x_min_in,
                            x_max_in,
                            y_min_in,
                            y_max_in,
                            z_min_in,
                            z_max_in
                    )
                    {} // Empty

                    inline real_t update_velocity(real_t v)
                    {
                        return v*-1.0;
                    }
                    inline real_t update_position(real_t v, real_t boundary)
                    {
                        return (boundary - (v-boundary));
                    }

                    inline int process_particle_x(
                            particle_list_t& p,
                            size_t i,
                            size_t v,
                            real_t extent,
                            real_t extent_far
                            ) final
                    {

                        // TODO: Can I get away with hard coding this? (if we store
                        // relative offsets into the cell)

                        // Invert velocity
                        real_t ux = p.get_ux(i, v);
                        ux = update_velocity(ux);

                        real_t x = p.get_x(i, v);
                        x = update_position(x, extent);

                        p.set_x(i, v, x);
                        p.set_ux(i, v, ux);

                        return 1;
                    }

                    inline int process_particle_y(
                            particle_list_t& p,
                            size_t i,
                            size_t v,
                            real_t extent,
                            real_t extent_far
                            ) final
                    {

                        // Invert velocity
                        real_t uy = p.get_uy(i, v);
                        uy = update_velocity(uy);

                        // Change position after bounding on boundary
                        real_t y = p.get_y(i, v);
                        y = update_position(y, extent);

                        p.set_y(i, v, y);
                        p.set_uy(i, v, uy);
                    }
                    inline int process_particle_z(
                            particle_list_t& p,
                            size_t i,
                            size_t v,
                            real_t extent,
                            real_t extent_far
                            ) final
                    {

                        // Invert velocity
                        real_t uz = p.get_uz(i, v);
                        uz = update_velocity(uz);

                        // Change position after bounding on boundarz
                        real_t z = p.get_z(i, v);
                        z = update_position(z, extent);

                        p.set_z(i, v, z);
                        p.set_uz(i, v, uz);
                    }
            };

            template <class particle_list_t>
                class PeriodicBoundary :
                    public BoundaryCondition<particle_list_t>
            {
                public:
                    PeriodicBoundary(
                            real_t x_min_in,
                            real_t x_max_in,
                            real_t y_min_in,
                            real_t y_max_in,
                            real_t z_min_in,
                            real_t z_max_in
                        ) : BoundaryCondition<particle_list_t>(
                            x_min_in,
                            x_max_in,
                            y_min_in,
                            y_max_in,
                            z_min_in,
                            z_max_in
                    )
                    {} // Empty

                    inline real_t update_velocity(real_t v)
                    {
                        return v;
                    }
                    inline real_t update_position(
                            real_t v,
                            real_t boundary,
                            real_t additional // The other boundary
                    )
                    {
                        return v - boundary + additional;
                    }

                    inline int process_particle_x(
                            particle_list_t& p,
                            size_t i,
                            size_t v,
                            real_t extent,
                            real_t extent_far
                            ) final
                    {

                        // TODO: Can I get away with hard coding this? (if we store
                        // relative offsets into the cell)

                        // Invert velocity
                        real_t ux = p.get_ux(i, v);
                        ux = update_velocity(ux);

                        real_t x = p.get_x(i, v);
                        x = update_position(x, extent, extent_far);

                        p.set_x(i, v, x);
                        p.set_ux(i, v, ux);

                        return 1;
                    }

                    inline int process_particle_y(
                            particle_list_t& p,
                            size_t i,
                            size_t v,
                            real_t extent,
                            real_t extent_far
                            ) final
                    {

                        // Invert velocity
                        real_t uy = p.get_uy(i, v);
                        uy = update_velocity(uy);

                        // Change position after bounding on boundary
                        real_t y = p.get_y(i, v);
                        y = update_position(y, extent, extent_far);

                        p.set_y(i, v, y);
                        p.set_uy(i, v, uy);
                    }
                    inline int process_particle_z(
                            particle_list_t& p,
                            size_t i,
                            size_t v,
                            real_t extent,
                            real_t extent_far
                            ) final
                    {

                        // Invert velocity
                        real_t uz = p.get_uz(i, v);
                        uz = update_velocity(uz);

                        // Change position after bounding on boundarz
                        real_t z = p.get_z(i, v);
                        z = update_position(z, extent, extent_far);

                        p.set_z(i, v, z);
                        p.set_uz(i, v, uz);
                    }
            };


        } // namespace pic
    } // namespace sp
} // namespace flecsi

#endif // flecsi_sp_pic_boundary_h
