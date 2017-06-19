/*
 * bitstring.hpp
 *
 *  Created on: Feb 25, 2017
 *      Author: xiangzhang
 */
#ifndef BITSTRING_HPP
#define BITSTRING_HPP

#include <vector>
#include <limits>

namespace phsim {

struct BitString {

  using Byte = unsigned char;
  using Word = unsigned long int;
  using Size = unsigned long int;
  using Data = std::vector<Word>;

  //WORDSIZE is 64-bit
  static const Byte WORDSIZE = 8 * sizeof(Word);

  static const Word WORDMAX = std::numeric_limits<Word>::max();

  BitString() = default;
  BitString(Data &&, Size);
  BitString(const Data &, Size);

  Data data;
  Size length = 0;
  std::vector <Data> AdjTable;
};

//inline BitString::BitString(Data &&d, Size l)
 //   : data(std::move(d)), length(l) { }

inline BitString::BitString(const Data &d, Size l)
    : data(d), length(l) { }
}

#endif
