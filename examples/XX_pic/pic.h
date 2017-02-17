/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~--------------------------------------------------------------------------~*/

#ifndef driver_h
#define driver_h

/// TODO ////
//
// boundary style things? periodic? 
// data layout? method to load particle properties?
//
/// TODO /// 

// Includes 
#include <iostream>

// Files from flecsi
#include <flecsi/data/data.h>

// Files from flecsi-sp
#include <flecsi-sp/pic/mesh.h>

// Namespaces
using namespace flecsi;
using namespace flecsi::data;
using namespace flecsi::sp;
using namespace flecsi::sp::pic;

using mesh_t = pic_mesh_t;
using vertex_t = pic_types_t::vertex_t;

using real_t = double;


// Constants
#define AoS // TODO: Need other types
#define NDIM 3 // Number of dimensions, i.e 3D 
// This can be pulled from the PIC class itself 

// Types
using dim_array_t = std::array<real_t,NDIM>;

const real_t q = 1; // TODO: Give these values
const real_t m = 1; // TODO: Give these values
const real_t mu = 4.0 * M_PI * 1.0e-7; // permeability of free space
const real_t c = 299792458; // Speed of light   
const real_t eps = 1.0 / (c * c * mu); // permittivity of free space

// TODO: Read input deck
size_t NX = 64;
size_t NY = 64;
size_t NZ = 64;
size_t NPPC = 32;
double dt = 0.1; // TODO: units?
int num_steps = 10;

real_t len_x = 1.0;
real_t len_y = 1.0;
real_t len_z = 1.0;

real_t dx = len_x/NX;
real_t dy = len_y/NY;
real_t dz = len_z/NZ;

// TODO:" Most of the uses for this can really be 1d attached to cells/vertexs
using three_d_real_vector_t = std::vector< std::vector< std::vector< real_t > >>;
// TEMP DATA
//
// TODO: Remove these add register them with flecsi. 
// These are here just so it compiles 
three_d_real_vector_t Ex;
three_d_real_vector_t Ey;
three_d_real_vector_t Ez;

three_d_real_vector_t Bx;
three_d_real_vector_t By;
three_d_real_vector_t Bz;

three_d_real_vector_t Jx;
three_d_real_vector_t Jy;
three_d_real_vector_t Jz;
// Methods
void init_mesh(mesh_t& m, size_t nx, size_t ny, size_t nz) {

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

/* // TODO: Don't do this as **
void load_particle_properties(
            real_t** dx, 
            real_t** dy, 
            real_t** dz, 
            real_t** ux, 
            real_t** uy, 
            real_t** uz
            )
{
#ifdef AoS 
  *dx = &particles[species][i].dx;
  *dy = &particles[species][i].dx;
  *dz = &particles[species][i].dz;

  *ux = &particles[species][i].ux;
  *uy = &particles[species][i].uy;
  *uz = &particles[species][i].uz;
#endif
}*/

dim_array_t cross_product( dim_array_t a, dim_array_t b) {

  dim_array_t result; 

  result[0] = a[1]*b[2] - a[2]*b[1];
  result[1] = a[2]*b[0] - a[0]*b[2];
  result[2] = a[0]*b[1] - a[1]*b[0];

  return result; 

}

void init_simulation() {
}

void field_solve(real_t dt) {
  
  // TODO: Double check these are the right values for EM (not ES)

  // General Idea, taken directly from _the_ Yee paper
  //
  // Electric Field 
    //  Take equation for D
    //  Dx_n - Dx_n-1 / dt = (Hz - Hz) / dy - (hy - hy) / dz + Jx 
    //  Rearrange:
    //  Dx_n = ((Hz - Hz) / dy - (hy - hy) / dz + Jx) *dt + Dx_n-1
    //  Swap out D and B:
    //  eps*Ex_n = ( 1/mu * ((Bz - Bz) / dy - (By - By) / dz) + Jx) * dt + eps*Ex_n-1
    //  Alternatively:
    //  Ex_n = (( 1/mu * ((Bz - Bz) / dy - (By - By) / dz) + Jx) * dt)/eps + Ex_n-1
    
  // Magnetic Field 
    // (Bx_h - Bx_-h) / dt  = (Ey - Ey) / dz - (Ez - Ez) / dy
    // Rearrange:
    // Bx = ((Ey - Ey) / dz - (Ez - Ez) / dy) * dt + Bx_-h
  
  // Note: Care needs to be taken as the sign changes from x to y to z 
  
  for (int k = 0; k < NZ; k++) {
    for (int j = 0; j < NY; j++) {
      for (int i = 0; i < NX; i++) {

        Bx[k][j][i] = Bx[k][j][i] + dt * ( (Ey[k+1][j][i] - Ey[k][j][i]) / dz - 
            (Ez[k][j+1][i] - Ez[k][j][i]) / dy );

        By[k][j][i] = By[k][j][i] + dt * ( (Ez[k+1][j][i] - Ez[k][j][i]) / dx - 
            (Ex[k][j+1][i] - Ex[k][j][i]) / dz );

        Bz[k][j][i] = Bz[k][j][i] + dt * ( (Ex[k+1][j][i] - Ex[k][j][i]) / dy - 
            (Ey[k][j+1][i] - Ey[k][j][i]) / dx );

        // TODO: Check these indexs
        Ex[k][j][i] = ( (1/(mu*eps)) * ( ((Bz[k][j][i] - Bz[k][j-1][i] ) / dy ) +
              ((By[k][j][i] - By[k][j][i-1]) / dz)) + Jx[k][j][i]) * dt + Ex[k][j][i];

        Ey[k][j][i] = ( (1/(mu*eps)) * ( ((Bx[k][j][i] - Bx[k][j-1][i] ) / dz ) +
              ((Bz[k][j][i] - Bz[k][j][i-1]) / dx)) + Jy[k][j][i]) * dt + Ey[k][j][i];

        Ez[k][j][i] = ( (1/(mu*eps)) * ( ((By[k][j][i] - By[k][j-1][i] ) / dx ) +
              ((Bx[k][j][i] - Bx[k][j][i-1]) / dy)) + Jz[k][j][i]) * dt + Ez[k][j][i];

      }
    }
  }


}
void particle_move() { 

  // mult by delta_t and divide by delta_x 
  // Swap in F/m = qE/m 
    // v = v + (F * delta_t) / m 
    // x = x + v * delta_t
    
  // Gives:
    // (v*delta_t)/delta_x = (v*delta_t / delta_x) + (q/m)*((E * delta_t^2) / delta_x)

  // KE = m/2 * v_old * v_new 
  //
  
  /* TODO: Implement this once we have a particle store
  dx += ux*dt;
  dy += uy*dt;
  dz += uz*dt;
  */
}

void update_velocity(real_t dt) {
  // Use the Boris method to update the velocity and rotate (P62 in Birdsall) 
  // Example of this can also be found here 
  // https://www.particleincell.com/wp-content/uploads/2011/07/ParticleIntegrator.java 
 
  // Equations:
  
  // v_old = v- - qE/m * delta_t / 2
    // => v- = v + qE/m * delta_t / 2
    
  // v_new = v+ + qE/m * delta_t / 2
  
  // v' = v- + v- X t 
  // v+ = v- + v' X s 
  
  // |v-|^2 = |v+|^2
  // s = 2t / (1+t^2) 
  // t = qB /m * delta_t / 2 
  //
  // TODO: move these data declarations
  dim_array_t t;
  dim_array_t E;
  dim_array_t B;
  dim_array_t v;
  dim_array_t s;
  dim_array_t v_plus;
  dim_array_t v_minus;
  dim_array_t v_prime;
 
  real_t t_squared = 0.0;

  // Calculate t and |t^2|
  for (size_t i = 0; i < NDIM; i++) 
  {
    t[i] = q/m * B[i] * 0.5 * dt;
    t_squared += t[i]*t[i];
  }

  // Calculate s
  for (size_t i = 0; i < NDIM; i++) 
  {
    s[i] = 2*t[i] / (1+t_squared);
  }

  // Calculate v-
  for (size_t i = 0; i < NDIM; i++) 
  {
    v_minus[i] = v[i] + q/m * E[i] * 0.5 * dt;
  }

  // Calculate v'
  dim_array_t vt_cross_product = cross_product( v_minus, t);
  for (size_t i = 0; i < NDIM; i++) 
  {
    v_prime[i] = v_minus[i] + vt_cross_product[i];
  }

  // Calculate v+
  dim_array_t vs_cross_product = cross_product( v_prime, s);
  for (size_t i = 0; i < NDIM; i++) 
  {
    v_plus[i] = v_minus[i] + vs_cross_product[i];
  }

  // Calculate v_new
  for (size_t i = 0; i < NDIM; i++) 
  {
    v[i] = v_plus[i] + q/m * E[i] * 0.5 * dt;
  }

  // TODO: Can hoist that q/m E dt/2?
  
}



void particle_push() {

}

void fields_half() {} 
void fields_final() {} 

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
  
  fields_half();

  particle_push();

  fields_final();

}


void particle_initialization() {

}



void driver(int argc, char ** argv) {

  // Declare mesh
  mesh_t m;


  // Init mesh
  init_mesh(m, NX, NY, NZ);

  // Init Simulation (generic particles etc)
  init_simulation();

  // Register data
  //register_data(m, solver, unknowns, double, dense, 1, m.vertices);
  
  //register_data(m, hydro, materials, double, dense, 1, m.cells);
  register_data(m, solver, unknowns, double, dense, 1, vertices);        

  //register_data(m, solver, p, double, sparse, 2, vertices, 3);
  
  // Experiment with data
  auto u = get_accessor(m, solver, unknowns, double, dense, 0);                    
                                                                                     
  std::cout << "top ent " << m.topology::mesh_topology_t<pic_types_t>::num_entities(0) << std::endl;
  for(auto v: m.vertices()) {                                                      
    std::cout << v->coordinates() << std::endl;                                    
    u[v] = 0.0;                                                                    
  } // for                  


  particle_initialization();

  // Do an initial velocity move of dt/2 backwards to enable the leap frogging
  update_velocity(-1 * (dt/2));

  // Perform main loop
  for (size_t i = 0; i < num_steps; i++) 
  {
    //kernel();
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
