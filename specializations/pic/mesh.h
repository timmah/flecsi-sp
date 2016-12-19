/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*/

#ifndef felcsisp_pic_mesh_h
#define felcsisp_pic_mesh_h

#include "flecsi-sp/geometry/point.h"
#include "flecsi-sp/pic/types.h"

///
// \file pic_mesh.h
// \authors bergen
// \date Initial file creation: Oct 19, 2016
///

namespace flecsi {
namespace sp {

///
// Use this namespace to expose enumerations and types.
///
namespace pic {

  enum pic_index_spaces_t : size_t {
    vertices,
    cells
  }; // enum pic_index_spaces_t

  enum pic_cell_index_spaces_t : size_t {
    interior,
    boundary
  }; // enum pic_cell_index_spaces_t

} // namespace pic

///
// \class pic_mesh_t pic_mesh.h
// \brief pic_mesh_t provides...
///
class pic_mesh_t
  : public topology::mesh_topology_t<pic_types_t>
{
public:

  using base_t = topology::mesh_topology_t<pic_types_t>;

  static constexpr size_t dimension = pic_config_t::num_dimensions;
  
  using vertex_t = pic_types_t::vertex_t;
  using cell_t = pic_types_t::cell_t;

  using point_t = pic_config_t::point_t;

  /// Default constructor
  pic_mesh_t()
    : base_t() {}

  /// Copy constructor (disabled)
  pic_mesh_t(const pic_mesh_t &) = delete;

  /// Assignment operator (disabled)
  pic_mesh_t & operator = (const pic_mesh_t &) = delete;

  /// Destructor
   ~pic_mesh_t() {}

  ///
  // Initialize the mesh.
  ///
  void
  init()
  {
    // Initialize domain 0 of the mesh topology.
    base_t::init<0>();

    // Use a predicate function to create the interior cells
    // index space
    interior_cells_ =
      base_t::entities<dimension, 0>().filter(is_interior);

    // Use a predicate function to create the domain boundary cells
    // index space
    boundary_cells_ =
      base_t::entities<dimension, 0>().filter(is_domain_boundary);
  } // init

  ///
  // Add a vertex to the mesh topology.
  ///
  vertex_t *
  make_vertex(
    const point_t & pos
  )
  {
    auto v = base_t::make<vertex_t>(*this, pos);
    base_t::add_entity<0, 0>(v);
    return v;
  } // make_vertex

  ///
  // Add a cell to the mesh topology
  ///
  cell_t *
  make_cell(
    const std::initializer_list<vertex_t *> & vertices,
    cell_type_t type
  )
  {
    auto c = base_t::make<cell_t>(*this, type);
    base_t::add_entity<dimension, 0>(c);
    base_t::init_entity<0, dimension, 0>(c, vertices);
    return c;
  } // make_cell

  ///
  //
  ///
  size_t
  indices(
    size_t index_space_id
  ) override
  {
    switch(index_space_id) {
      case pic::vertices:
        return base_t::num_entities(0);
      case pic::cells:
        return base_t::num_entities(dimension);
      default:
        assert(false && "unknown index space");
    } // switch
  } // indices

  ///
  //
  ///
  auto
  vertices()
  {
    return base_t::entities<0, 0>();
  } // vertices

  ///
  //
  ///
  template<
    typename E
  >
  auto
  vertices(
    E * e
  )
  {
    return base_t::entities<0, 0>(e);
  } // vertices

  ///
  //
  ///
  auto
  cells()
  {
    return base_t::entities<dimension, 0>();
  } // cells

  ///
  //
  ///
  template<
    typename E
  >
  auto
  cells(E * e)
  {
    return base_t::entities<dimension, 0>(e);
  } // cells

  auto
  cells(
    size_t is
  )
  {
    switch(is) {
      case pic::interior:
        return interior_cells_;
      case pic::boundary:
        return boundary_cells_;
      default:
        assert(false && "unknown index space");
    } // switch
  } // cells

private:

  ///
  // Predicate function to create index space for accessing
  // domain boundary cells.
  ///
  static
  bool
  is_domain_boundary(
    cell_t * c
  )
  {
    return c->type() == cell_type_t::domain_boundary;
  } // is_domain_boundary

  ///
  // Predicate function to create index space for accessing
  // interior cells.
  ///
  static
  bool
  is_interior(
    cell_t * c
  )
  {
    return !is_domain_boundary(c);
  } // is_interior

  //topology::index_space<
    //particle_t, true, true, false> particles_;

  topology::index_space<
    topology::domain_entity<0, cell_t>, false, true, false> interior_cells_;

  topology::index_space<
    topology::domain_entity<0, cell_t>, false, true, false> boundary_cells_;

}; // class pic_mesh_t

} // namespace sp
} // namespace flecsi

#endif // felcsisp_pic_mesh_h

/*~-------------------------------------------------------------------------~-*
 * Formatting options for vim.
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
