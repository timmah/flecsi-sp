/*~--------------------------------------------------------------------------~*
 *
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~--------------------------------------------------------------------------~*/

#ifndef driver_h
#define driver_h

///
//
///

#include <iostream>

// Files from flecsi
#include <flecsi/data/data.h>

// Files from flecsi-sp
#include <flecsi-sp/minimal/minimal_mesh.h>

using namespace flecsi;
using namespace flecsi::data;
using namespace flecsi::sp;
using namespace flecsi::sp::minimal;

using vertex_t = minimal_mesh_t::vertex_t;
using real_t = float;

static constexpr size_t MESH_SIZE = 8;

///
// Initialize a basic mesh.
///
void init_mesh(minimal_mesh_t & m) {
  std::vector<vertex_t *> vs;

    for(size_t j(0); j<MESH_SIZE+1; ++j) {
      for(size_t i(0); i<MESH_SIZE+1; ++i) {
        vs.push_back(m.make_vertex({double(i), double(j)}));
      } // for
    } // for

    size_t width = N+1;

    for(size_t j(0); j<MESH_SIZE; ++j) {
      for(size_t i(0); i<MESH_SIZE; ++i) {
        m.make_cell({
          vs[ i    + ( j    * width)],
          vs[(i+1) + ( j    * width)],
          vs[(i+1) + ((j+1) * width)],
          vs[ i    + ((j+1) * width)]
        });
      } // for
    } // for

    m.init();

} // init_mesh

void register_data(minimal_mesh_t& m)
{
  // Need to register data for:
    // Fields
    // Particles

  // Format for register_data(
    // Mesh
    // Data store (?)  (name space?)
    // Name (key)
    // Type 
    // Sparse/Dense
    // ? 
    // Attachment site
  // );

  // DEPRICATED?
  // Format for register_state( Mesh, Name (key), Attachment Site, Data_type,
    // persistent, );

  // Register Particles 
    // 
  //register_data(m, sovler, particles, real_t, dense, 2, vertices);
  
  // Fields
  register_data(m, sovler, fields, real_t, dense, 2, vertices);


}

void driver(int argc, char ** argv) {
	minimal_mesh_t m;

	// Initialize the mesh
	//
	// For real programs, mesh initialization will be handled through
	// an I/O object.
	init_mesh(m);

  register_data(m);



  auto u = get_accessor(m, solver, unknowns, real_t, dense, 0);

  for(auto v: m.vertices()) {
    u[v] = 0.0;
  } // for

} // driver

#endif // driver_h

/*~-------------------------------------------------------------------------~-*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
