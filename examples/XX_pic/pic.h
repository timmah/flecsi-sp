/*~--------------------------------------------------------------------------~*
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

using namespace flecsi;
using namespace flecsi::data;

void init_mesh(mesh_t* m, size_t nx, size_t ny, size_t nz) {

  std::vector<vertex_t *> vs;

  for(size_t k(0); k<nz+1; ++k) {
    for(size_t j(0); j<ny+1; ++j) {
      for(size_t i(0); i<nx+1; ++i) {

        vs.push_back(m.make_vertex({double(i), double(j), double(k)}));

      } // for
    } // for
  } // for

  size_t width = nx+1;
  size_t depth = ny+1;

  for(size_t k(0); k<nz; ++k) { // z
    for(size_t j(0); j<ny; ++j) { // y
      for(size_t i(0); i<nx; ++i) { // x

        bool is_domain_boundary = i==0 || j==0 || i==(nx-1) || j==(ny-1) || 
          k == 0 || k == (nz-1); 

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
            }, is_domain_boundary ? cell_type_t::domain_boundary :
            cell_type_t::unknown);
      } // for
    } // for
  } // for

  m.init();
}

void init_simulation() {
}

// Test main PIC kernel implementation
void kernel() {


  // Equations of motion
  // Leap frog method. 
    // TODO: Does this change in 3d?
    // Integrate force and velocity independently 
    //
    //  m * (dv)/(dt) = F 
    //      (dx)/(dt) = V
    //
    //  Can do this as a finite difference. 
    //
    //  m * (v_new - v_old) / (delta t) = F_old
    //      (x_new - x_old) / (delta t) = v_new
    //
    // Advances v_t to v_(t+delta t) and x_t to x_(t+delta_t)
      // v and x however are at difference times 
     
    // TODO: Push this back
    // Push v back to v_(t-(delta t / 2)) 
     
    // F has two parts 
    // F = F_electric + F_magnetic 
    // F = qE + q(v x B) 
      // Cross product implies a rotation?
      // E and B calculated at particle 
      // Interpolate E and B from grid to the particle
      
    // Field Equations
    
    // Yee Grid
      // Replaces Maxwell's equations with a set of finite difference equations
      
  // FDTD (Finite-Difference Time-Domain)


}





void driver(int argc, char ** argv) {

  // Declare mesh
  mesh_t m;

  // TODO: Read input deck
  size_t NX = 64;
  size_t NY = 64;
  size_t NZ = 64;
  size_t NPPC = 32;
  double dt = 0.1; // TODO: units?

  // Init mesh
  init_mesh(m, NX, NY, NZ);

  // Init Simulation (generic particles etc)
  init_simulation();

  // Register data
  register_data(m, solver, unknowns, double, dense, 2, vertices);
  register_data(m, solver, p, double, sparse, 2, vertices, 3);


  // Perform main loop
  for (size_t i = 0; i < num_steps; i++) 
  {
    kernel();
  }

  /*
  auto ap = get_mutator(m, solver, p, double, sparse, 0, 3);

  for(auto v: m.vertices()) 
  {
    for(size_t i(0); i<3; ++i) 
    {
      ap(v, i) = 1.0;
    } // for
  } // for

  auto u = get_accessor(m, solver, unknowns, double, dense, 0);

  for(auto v: m.vertices()) {
    std::cout << v->coordinates() << std::endl;
    u[v] = 0.0;
  } // for
  */

} // driver

#endif // driver_h

/*~-------------------------------------------------------------------------~-*
 * Formatting options
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
