#ifndef flecsi_sp_pic_helpers_h
#define flecsi_sp_pic_helpers_h

#include "flecsi-sp/pic/types.h"

#if ENABLE_DEBUG
  #define logger std::cout << "LOG:" << __FILE__ << ":" << __LINE__ << " \t :: \t " 
#else
#  define logger while(0) std::cout
#endif /* ENABLE_DEBUG */

using dim_array_t = std::array<real_t,NDIM>;
//using a = flesci::sp::dim_array_t;

////// HELPER METHODS //////
//
/** 
 * @brief Helper function to generate a random real between two bounds.
 * Most useful in particle placement
 * 
 * @param min Lower bound of the number to generate
 * @param max Upper bound of the number to generate
 * 
 * @return Random number between bounds
 */
double random_real(real_t min, real_t max)
{
  double r = (double)rand() / (double)RAND_MAX;
  return min + r * (max - min);
}
/** 
 * @brief Helper function to calculate the cross produce of two arrays
 * 
 * @param a First array
 * @param b Second array
 * 
 * @return Result of cross product calculation
 */
dim_array_t cross_product(dim_array_t a, dim_array_t b) 
{
  dim_array_t result; 

  result[0] = a[1]*b[2] - a[2]*b[1];
  result[1] = a[2]*b[0] - a[0]*b[2];
  result[2] = a[0]*b[1] - a[1]*b[0];

  return result; 
}

////// END HELPERS ////// 
#endif // flecsi_sp_pic_helpers_h
