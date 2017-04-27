#ifndef flecsi_sp_pic_visualization_h
#define flecsi_sp_pic_visualization_h

#include <iostream>
#include <fstream>

namespace flecsi {
    namespace sp {
        namespace pic {
            class Visualizer {

                public:
                    std::ofstream vis_file;

                    void write_header(size_t total_num_particles, size_t step) {

                        std::stringstream sstm;

                        sstm << "vis" << step << ".vtk";
                        std::string file_name = sstm.str();

                        vis_file.open(file_name);

                        vis_file << "# vtk DataFile Version 2.0" << std::endl;
                        vis_file << "Unstructured Grid Example" << std::endl;
                        vis_file << "ASCII" << std::endl;
                        vis_file << "" << std::endl;
                        vis_file << "DATASET UNSTRUCTURED_GRID" << std::endl;

                        vis_file << "POINTS " << total_num_particles << " float" << std::endl;
                    }

                    void write_particles(auto particles_accesor, size_t num_particles, mesh_t& m)
                    {

                        size_t write_count = 0;
                        for ( auto c : m.cells() ) {

                            auto& cell_particles = particles_accesor[c];

                            // TODO: Is there a way abstract this loop structure with the current particle structure
                            // TODO: this may need to be "active ppc" or similar
                            //
                            // Only iterate over used blocks
                            for (size_t i = 0; i < cell_particles.block_number+1; i++)
                            {

                                for (size_t v = 0; v < PARTICLE_BLOCK_SIZE; v++)
                                {
                                    // TODO: Does this need masking for the empty unfilled blocks

                                    real_t x = cell_particles.get_x(i, v);
                                    real_t y = cell_particles.get_y(i, v);
                                    real_t z = cell_particles.get_z(i, v);

                                    // TODO: Find a better way to filter particles out, USE MASK
                                    if (x == 0 && y == 0 && z == 0) continue;

                                    vis_file << x << " " << y << " " << z << std::endl;
                                    write_count++;
                                }
                            }
                        }

                        /*
                           for (unsigned int sn = 0; sn < species.size(); sn++) {

                           int particle_count = species[sn].num_particles;
                           for (int particle_number = 0; particle_number < particle_count; particle_number++) {

                           double* part_x_in;
                           double* part_y_in;

                           load_particle_position(sn, particle_number, &part_x_in, &part_y_in);

                           vis_file << *part_x_in << " " << *part_y_in << " 0.0" << std::endl;
                           }
                           }
                           */

                        vis_file << "CELL_TYPES " << num_particles << std::endl;

                        for (size_t p = 0; p < num_particles; p++)
                        {
                            vis_file << "1" << std::endl;
                        }

                        vis_file << "POINT_DATA " << num_particles << std::endl;
                        vis_file << "SCALARS weight float 1"  << std::endl;
                        vis_file << "LOOKUP_TABLE default" << std::endl;

                        for ( auto c : m.cells() ) {

                            auto& cell_particles = particles_accesor[c];
                            for (size_t i = 0; i < cell_particles.block_number+1; i++)
                            {
                                for (size_t v = 0; v < PARTICLE_BLOCK_SIZE; v++)
                                {
                                    real_t x = cell_particles.get_x(i, v);
                                    real_t y = cell_particles.get_y(i, v);
                                    real_t z = cell_particles.get_z(i, v);

                                    real_t w = cell_particles.get_w(i, v);

                                    // TODO: Find a better way to filter particles out, USE MASK
                                    if (x == 0 && y == 0 && z == 0) continue;

                                    vis_file << w << std::endl;
                                }
                            }
                        }

                        /*
                    //0.0 1.0 2.0 3.0 4.0 5.0 6.0 7.0
                    //for (unsigned int species_number = 0; species_number < species.size(); species_number++) {
                        //Species& this_species = get_species(species_number);
                        //int particle_count = this_species.get_num_particles();
                        for (size_t particle_number = 0; particle_number < num_particles; particle_number++) {
                            // TODO: Figure out what quantities are useful to output.
                            double weight = 1.1;
                            vis_file << weight << std::endl;
                        }
                    //}
                    */

                    }
                    void finalize()
                    {
                        vis_file.close();
                    }

        };
    } // namespace pic
} // namespace sp
} // namespace flecsi

#endif // Visualizer
