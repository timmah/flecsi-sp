/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2015 Los Alamos National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*/

#ifndef flecsi_sp_pic_entity_types_h
#define flecsi_sp_pic_entity_types_h

#include <flecsi/topology/mesh_types.h>

#include "flecsi-sp/pic/config.h"
#include "flecsi-sp/pic/particles.h"

///
// \file pic_entity_types.h
// \authors bergen bird
// \date Initial file creation: Oct 19, 2016
///

namespace flecsi {
  namespace sp {
    namespace pic {

      enum class entity_type_t : size_t {
        unknown,
        domain_boundary
      }; // enum entity_type_t


      ///
      // \struct pic_vertex_t
      // \brief FIXME
      ///
      struct pic_vertex_t
        : public topology::mesh_entity_t<0,
        pic_config_t::num_domains>
      {
        using point_t = pic_config_t::point_t;

        ///
        // Constructor.
        //
        // \param mesh FIXME
        ///
        pic_vertex_t(
            topology::mesh_topology_base_t & mesh,
            const point_t & coordinates,
            entity_type_t type
            )
          : mesh_(mesh), coordinates_(coordinates), type_(type)
        {
        }

        const point_t &
          coordinates()
          const
          {
            return coordinates_;
          } // coordinates

        entity_type_t
          type()
          {
            return type_;
          } // type

        private:

        topology::mesh_topology_base_t & mesh_;
        point_t coordinates_;
        entity_type_t type_;

      }; // class pic_vertex_t

      /*
         enum class cell_type_t : size_t {
         unknown,
         domain_boundary
         }; // enum cell_type_t
         */

      ///
      // \struct pic_cell_t
      // \breif FIXME
      ///
      struct pic_cell_t
        : public topology::mesh_entity_t<FLECSI_MESH_DIMENSION,
        pic_config_t::num_domains>
      {
        using real_t = pic_config_t::real_t;

        ///
        // Constructor.
        //
        // \param mesh FIXME
        ///
        pic_cell_t(
            topology::mesh_topology_base_t & mesh, 
            entity_type_t type
            )
          : mesh_(mesh), type_(type)
        {
        }

        std::vector<size_t>
          create_entities(
              id_t cell_id,
              size_t dim,
              topology::domain_connectivity<FLECSI_MESH_DIMENSION> & c,
              id_t * e
              )
          {
            // FIXME
            return {};
          } // create_entities

        real_t
          volume()
          {
            // FIXME
            return 1.0;
          } // volume

        entity_type_t
          type()
          {
            return type_;
          } // type

        private:

        topology::mesh_topology_base_t & mesh_;
        entity_type_t type_;

      }; // class pic_cell_t


      /* Moved to particles.h
      // TODO :Where is the correct place to define particles
      struct particle_t {

      using id_t = size_t;
      using real_t = pic_config_t::real_t;

      real_t dx, dy, dz; 
      real_t ux, uy, uz; 

      real_t w;          
      int32_t i;        
      }; // particle_t

*/


    } // namespace pic
  } // namespace sp
} // namespace flecsi

#endif // flecsi_sp_pic_entity_types_h

/*~-------------------------------------------------------------------------~-*
 * Formatting options for vim.
 * vim: set tabstop=2 shiftwidth=2 expandtab :
 *~-------------------------------------------------------------------------~-*/
