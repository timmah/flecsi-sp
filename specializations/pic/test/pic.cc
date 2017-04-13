/*~-------------------------------------------------------------------------~~*
 * Copyright (c) 2014 Los Alamos National Security, LLC
 * All rights reserved.
 *~-------------------------------------------------------------------------~~*/

// Note: This currently doesn't test anything of note, and is more a sanity
// check to ensure the files all compile

#include <cinchtest.h>
#include <vector>

#include <flecsi/topology/mesh_topology.h>
#include "flecsi/data/data.h"
#include "flecsi/utils/common.h"
#include "flecsi-sp/pic/mesh.h"

using namespace flecsi;
using namespace flecsi::sp::pic;
using namespace flecsi::topology;
using vertex_t = mesh_t::vertex_t;
//using real_t = float;

class pic_t
  : public ::testing::Test
{
protected:

  static constexpr size_t N = 3;
  mesh_t m;

  void
  SetUp()
  {

    std::vector<vertex_t *> vs;

    for(size_t k(0); k<N+1; ++k) {
      for(size_t j(0); j<N+1; ++j) {
        for(size_t i(0); i<N+1; ++i) {
          // TODO: Check this logic
          bool is_domain_boundary = i==0 || j==0 || i==(N-1) || j==(N-1) || 
            k == 0 || k == (N-1); 
          vs.push_back(
              m.make_vertex({double(i), double(j), double(k)},
              is_domain_boundary ? entity_type_t::domain_boundary :
                entity_type_t::unknown
          ));
        } // for
      } // for
    } // for

    size_t width = N+1;
    size_t depth = N+1;

    for(size_t k(0); k<N; ++k) {
      for(size_t j(0); j<N; ++j) {
        for(size_t i(0); i<N; ++i) {
          bool is_domain_boundary = i==0 || j==0 || i==(N-1) || j==(N-1) || 
            k == 0 || k == (N-1); 
          // Indexing x + y*WIDTH + Z*WIDTH*DEPTH
          m.make_cell({
              vs[ (k    *width*depth) + (j    *width) + i ], // {0, 0, 0} 
              vs[ (k    *width*depth) + (j    *width) + (i+1) ], // {0, 0, 1} 
              vs[ (k    *width*depth) + ((j+1)*width) + (i+1) ], // {0, 1, 1} 
              vs[ (k    *width*depth) + ((j+1)*width) + i ], // {0, 1, 0} 
              vs[ ((k+1)*width*depth) + ((j+1)*width) + i ], // {1, 1, 0} 
              vs[ ((k+1)*width*depth) + (j    *width) + i ], // {1, 0, 0} 
              vs[ ((k+1)*width*depth) + (j    *width) + (i+1) ], // {1, 0, 1} 
              vs[ ((k+1)*width*depth) + ((j+1)*width) + (i+1) ], // {1, 1, 1} 
          }, is_domain_boundary ? entity_type_t::domain_boundary :
            entity_type_t::unknown);
        } // for
      } // for
    }
    

    m.init();

  } // SetUp

  virtual void TearDown() {}

}; // class pic_t

TEST_F(pic_t, sanity) {

  // In tree_topology they define:
    // An entity vector
      //  using entity_t = typename Policy::entity_t;
      //  using entity_vector_t = std::vector<entity_t*>;
      
    // An entity space using index
      //  using entity_space_t = index_space<entity_t*, true, true, false>;
    
    // Then store entitys by
    //   entity_space_t entities_;
    //
    //   entity_t* make_entity(Args&&... args)
    //   {
    //     auto ent = new entity_t(std::forward<Args>(args)...);
    //     entity_id_t id = entities_.size();
    //     ent->set_id_(id);
    //     entities_.push_back(ent);
    //     return ent;
    //   }
    //  
    //  For the simple_entry stuff:
        //  index_space_ is;
        //  is.push_back({id++, entry});
          //   Where: size_t entry = itr->entry; (size_t?!)
    
    using entity_t = particle_<real_t>; //typename Policy::entity_t;
    using entity_vector_t = std::vector<entity_t*>;

    using simple_t = topology::simple_entry<size_t>;
    topology::index_space<simple_t, false, true, false> particles_is;  

    // TODO: How do I register data to a index space?
    // This isn't (currently) supported, index into the particle array using
    // the index space directly
    
    //register_data(m, particles, electrons, particle_t, global, particles_is); 



    // this is just lots of attempts to do something useful ..


    //register_data(m, hydro, pressure, double, global, 1);
    //using entity_space_t = topology::index_space<entity_t*, true, true, false>;
    // TODO: This should use a simple_entry?

    //entity_space_t entities_; // this is particles_

    //particle_t entry = new Particle();

    //register_data(...);
    //using simple_t = simple_entry<size_t>;
    //topology::index_space<simple_t, false, true, false> particles_;  

    //particles_.push_back(&entry);

    //std::cout << m.particles_[0].index_space_id() << std::endl;

    //entity_vector_t ev;
    //m.particles_.push_back({id++, entry});
    //m.particles_.push_back(&p);

    // TODO: This is where I'm going to draft out the algorithm  I want to use.
    //  Hopefully doing it this way will help expose what I need access to
   
    // TODO: How does it know how many particles I want?
    // Is it just backed by a vector which I can push_back onto?
    //topology::index_space<
    //
      //topology::simple_entry<particle_t*>, true, true, false> particles_;
      //
      //

  /*
  //From ./data/serial/sparse.h
    index_space_ is;

    size_t id = 0;
    while(itr != end){
      is.push_back({id++, itr->entry});
      ++itr;
    }

    return is; 
  */

  //////// Plan ////////
  //
  // Define a grid. 
  //
  //
  ///// METHOD /////
  // Leap frog mehtod. Define velocity and position at different times, and do:
  //
  // 1D (?) Finite difference
  // F_old = m( (v_new - v_old) / delta t ) 
  // V_new = (x_new - x_old) / delta t 
  //
  // But they're at different times. Need to push v back to -1/2(delta t) 
  // using F at t=0. Also need to update kinetic and potential/field energies 
  //
  // Field Equations:
  // F = F_e + F_B 
  // F = qE + q(v X B)  (cross product?)
  // E and B are done at the particle, and must be interpolated.


} // TEST

/*----------------------------------------------------------------------------*
 * Cinch test Macros
 *
 *  ==== I/O ====
 *  CINCH_CAPTURE()              : Insertion stream for capturing output.
 *                                 Captured output can be written or
 *                                 compared using the macros below.
 *
 *    EXAMPLE:
 *      CINCH_CAPTURE() << "My value equals: " << myvalue << std::endl;
 *
 *  CINCH_COMPARE_BLESSED(file); : Compare captured output with
 *                                 contents of a blessed file.
 *
 *  CINCH_WRITE(file);           : Write captured output to file. 
 *
 *  CINCH_ASSERT(ASSERTION, ...) : Call Google test macro and automatically
 *                                 dump captured output (from CINCH_CAPTURE)
 *                                 on failure.
 *
 *  CINCH_EXPECT(ASSERTION, ...) : Call Google test macro and automatically
 *                                 dump captured output (from CINCH_CAPTURE)
 *                                 on failure.
 *
 * Google Test Macros
 *
 * Basic Assertions:
 *
 *  ==== Fatal ====             ==== Non-Fatal ====
 *  ASSERT_TRUE(condition);     EXPECT_TRUE(condition)
 *  ASSERT_FALSE(condition);    EXPECT_FALSE(condition)
 *
 * Binary Comparison:
 *
 *  ==== Fatal ====             ==== Non-Fatal ====
 *  ASSERT_EQ(val1, val2);      EXPECT_EQ(val1, val2)
 *  ASSERT_NE(val1, val2);      EXPECT_NE(val1, val2)
 *  ASSERT_LT(val1, val2);      EXPECT_LT(val1, val2)
 *  ASSERT_LE(val1, val2);      EXPECT_LE(val1, val2)
 *  ASSERT_GT(val1, val2);      EXPECT_GT(val1, val2)
 *  ASSERT_GE(val1, val2);      EXPECT_GE(val1, val2)
 *
 * String Comparison:
 *
 *  ==== Fatal ====                     ==== Non-Fatal ====
 *  ASSERT_STREQ(expected, actual);     EXPECT_STREQ(expected, actual)
 *  ASSERT_STRNE(expected, actual);     EXPECT_STRNE(expected, actual)
 *  ASSERT_STRCASEEQ(expected, actual); EXPECT_STRCASEEQ(expected, actual)
 *  ASSERT_STRCASENE(expected, actual); EXPECT_STRCASENE(expected, actual)
 *----------------------------------------------------------------------------*/

/*~------------------------------------------------------------------------~--*
 * Formatting options for vim.
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~------------------------------------------------------------------------~--*/
