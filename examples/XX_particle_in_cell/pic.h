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
#include <flecsi-sp/pic/mesh.h>
#include <flecsi-sp/pic/entity_types.h>

using namespace flecsi;
using namespace flecsi::data;
using namespace flecsi::sp;
using namespace flecsi::sp::pic;

// TODO: This should be PIC?
using mesh_t = pic_mesh_t;
using vertex_t = pic_vertex_t;
using namespace flecsi::sp;
using real_t = float;

static constexpr size_t MESH_SIZE = 8;

///
// Initialize a basic mesh.
///
void init_mesh(mesh_t & m, size_t N) {
  std::vector<vertex_t *> vs;

  for(size_t k(0); k<MESH_SIZE+1; ++k) {
    for(size_t j(0); j<MESH_SIZE+1; ++j) {
      for(size_t i(0); i<MESH_SIZE+1; ++i) {
        // TODO: Think about the order of these paramaters 
        vs.push_back(m.make_vertex({double(i), double(j), double(k)}));
      } // for
    } // for
  } // for

    size_t width = MESH_SIZE+1;
    size_t height = MESH_SIZE+1;

    for(size_t k(0); k<MESH_SIZE; ++k) {
      for(size_t j(0); j<MESH_SIZE; ++j) {
        for(size_t i(0); i<MESH_SIZE; ++i) {
          // x + y*WIDTH + Z*WIDTH*DEPTH
          m.make_cell(
              { // [x,y,z]
                vs[ (i+0) + ((j+0)*width) + ((k+0)*width*height)], // [0,0,0]
                vs[ (i+1) + ((j+0)*width) + ((k+0)*width*height)], // [1,0,0]
                vs[ (i+1) + ((j+1)*width) + ((k+0)*width*height)], // [1,1,0]
                vs[ (i+0) + ((j+1)*width) + ((k+0)*width*height)], // [0,1,0]
                vs[ (i+0) + ((j+1)*width) + ((k+1)*width*height)], // [0,1,1]
                vs[ (i+0) + ((j+0)*width) + ((k+1)*width*height)], // [0,0,1]
                vs[ (i+0) + ((j+1)*width) + ((k+1)*width*height)], // [0,1,1]
                vs[ (i+1) + ((j+1)*width) + ((k+1)*width*height)] // [1,1,1]
              },
              cell_type_t::unknown
          );       


        } // for
      } // for
    } // for

    m.init();

} // init_mesh

void register_pic_data(mesh_t& m)
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
    // TODO: Ben how do I do this?
    // I get that it's ~:
    //register_data(, particles, electrons, particle_t, dense, 1, is);
    // But I'm struggling find docs for the specifics 
  
  // Fields
  register_data(m, fields, jx, real_t, dense, 1, vertices);
  // TODO: Etc


}

void driver(int argc, char ** argv) {
	mesh_t m;
  const size_t N = 64;

	// Initialize the mesh
	//
	// For real programs, mesh initialization will be handled through
	// an I/O object.
	init_mesh(m, N);

  register_pic_data(m);



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
