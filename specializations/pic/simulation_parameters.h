#ifndef flecsi_sp_pic_simulation_paramaters_h
#define flecsi_sp_pic_simulation_paramaters_h

// TODO: Document this
// This class could capture all of the "required" (or "assumed") fields needed
// for a PIC simulation using this specialization.
// Users of this specialization could then subclass this class (or even the
// specialization if they want to write one) to extend the support it offers?
//
// Much of this can be populated (or calculated) during input deck reading. 
// It could also be used as a mechanism to provide sensible default values?
//
// This can be extended by users to add additional fields?
//
// Perhaps this should be a singleton? If it was a static class people couldn't extend?

namespace flecsi {
    namespace sp {
        namespace pic {

            template <class real_t> class Parameters_
            {

                // TODO: Setters and getters
                public: 
                    //Ben way of doing things
                    static Parameters_& instance()
                    {
                        static Parameters_ instance_;
                        return instance_;
                    }

                    /*
                    // Design patterns book way (favoring lazy initialization)
                    static Parameters_* instance_;

                    static Parameters_& instance()
                    {
                        // TODO: move this init out to a function?
                        if (!instance_)
                        {
                            instance_ = new Parameters_;
                        }
                        return instance_;
                    }
                    */

                    // Consts
                    const real_t mu = 4.0 * M_PI * 1.0e-7; // permeability of free space
                    const real_t c = 299792458; // Speed of light
                    const real_t eps = 1.0 / (c * c * mu); // permittivity of free space

                    size_t num_species = 0;

                    size_t NX_global = 64;
                    size_t NY_global = 64;
                    size_t NZ_global = 64;

                    size_t nx = NX_global;
                    size_t ny = NY_global;
                    size_t nz = NZ_global;

                    size_t NPPC = 32;

                    double dt = 0.1; // TODO: units?
                    int num_steps = 10;

                    real_t len_x_global = 1.0;
                    real_t len_y_global = 1.0;
                    real_t len_z_global = 1.0;

                    real_t len_x = len_x_global;
                    real_t len_y = len_y_global;
                    real_t len_z = len_z_global;

                    // Assume domain starts at [0,0,0] and goes to [len,len,len]
                    real_t local_x_min = 0.0;
                    real_t local_y_min = 0.0;
                    real_t local_z_min = 0.0;

                    real_t local_x_max = len_x;
                    real_t local_y_max = len_y;
                    real_t local_z_max = len_z;

                    real_t dx = len_x/nx;
                    real_t dy = len_y/ny;
                    real_t dz = len_z/nz;
            };

            //int Parameters_::instance2() { }

        } // namespace pic
    } // namespace sp
} // namespace flecsi

/*
flecsi::sp::pic::Parameters_<> flecsi::sp::pic::Parameters_<>::instance2()
{
    return instance_;
}
*/

#endif // flecsi_sp_pic_simulation_paramaters_h
