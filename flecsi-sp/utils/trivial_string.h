/*~-------------------------------------------------------------------------~~*
 * Copyright (c) 2016 Los Alamos National Laboratory, LLC
 * All rights reserved
 *~-------------------------------------------------------------------------~~*/
////////////////////////////////////////////////////////////////////////////////
/// \file
/// \brief Some commonly used utilities.
////////////////////////////////////////////////////////////////////////////////
#pragma once

// user includes
#include "flecsi-sp/common/constants.h"

// system includes
#include <array>
#include <cstring>
#include <string>

namespace flecsi {
namespace sp {
namespace utils {

////////////////////////////////////////////////////////////////////////////////
//! \brief A constexpr, trivially copyable string.
////////////////////////////////////////////////////////////////////////////////
class trivial_string_t {

public:

  using size_type = std::size_t;
  using value_type = char;
  using pointer = value_type *;
  using const_pointer = const value_type *;

  //! \brief return the value at a specific position `n`
  auto & operator[] ( size_type n )
  { 
    assert( n <= size() );
    return str_[n]; 
  }

  //! \brief return the value at a specific position `n`
  const auto & operator[] ( size_type n ) const
  { 
    assert( n <= size() );
    return str_[n]; 
  }

  //! \brief Return the size of the string.
  static constexpr size_type size() { return size_; }

  //! \brief Return a pointer to a c-string.
  const_pointer c_str() const { return str_.data(); }

  //! \brief Return a regular string.
  std::string str() const { return str_.data(); }

  //! \brief Return a pointer to the data.
  const_pointer data() const { return str_.data(); }
  pointer data() { return str_.data(); }

private:

  //! the length of the string
  static constexpr size_type size_ = common::max_char_length;

  //! a trivially copyable character array
  std::array<value_type, size_> str_;

};

///////////////////////////////////////////////////////////////////////////////
//! \brief Convert a std::string to a character array
///////////////////////////////////////////////////////////////////////////////
static auto to_trivial_string( const std::string & str )
{
  trivial_string_t tmp;
  strcpy( tmp.data(), str.c_str() );
  return tmp;
}

} // namespace
} // namespace
} // namespace
