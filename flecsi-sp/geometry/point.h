/*~-------------------------------------------------------------------------~~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~-------------------------------------------------------------------------~~*/
////////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief Provides a dimensioned array which functions as a vector.
////////////////////////////////////////////////////////////////////////////////
#pragma once

// user includes

#include "flecsi-sp/math/vector.h"

#include <flecsi/utils/common.h>

#include <cmath>

namespace flecsi{
namespace sp {
namespace geometry{

////////////////////////////////////////////////////////////////////////////////
//!  \brief The dimensioned_array type provides a general base for defining
//!  contiguous array types that have a specific dimension.
//!
//!  \tparam T The type of the array, e.g., P.O.D. type.
//!  \tparam D The dimension of the array, i.e., the number of elements
//!    to be stored in the array.
////////////////////////////////////////////////////////////////////////////////
template <typename T, std::size_t D> 
using point = math::vector<T,D>;


///
// \function distance
///
template <typename T, size_t D>
T distance(const point<T, D> & a, const point<T, D> & b)
{
  T sum(0);
  for (size_t d(0); d < D; ++d) {
    sum += flecsi::utils::square(a[d] - b[d]);
  } // for

  return std::sqrt(sum);
} // distance
  
} // namespace geometry
} // namespace sp
} // namespace flecsi
