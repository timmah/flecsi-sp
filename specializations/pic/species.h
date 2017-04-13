#ifndef flecsi_sp_pic_species_h
#define flecsi_sp_pic_species_h

namespace flecsi {
    namespace sp {
        namespace pic {

            template <class real_t> class species_
            {
                public: 

                    // Species properties
                    real_t charge;
                    real_t weight;
            };



        } // namespace pic
    } // namespace sp
} // namespace flecsi

#endif // flecsi_sp_pic_species_h

