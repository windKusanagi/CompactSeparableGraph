/*
 * bitstring_writer.hpp
 *
 *  Created on: Feb 25, 2017
 *      Author: xiangzhang
 */

#ifndef BITSTRING_WRITER_HPP
#define BITSTRING_WRITER_HPP

#include "bitstring.h"

namespace phsim {

class BitStringWriter {

 public:
  explicit BitStringWriter(BitString &);
  ~BitStringWriter();

  void writeBit(bool);
  void writeBits(BitString::Byte, BitString::Word);



  BitString::Word __buffer;

  //a variable to calculate how many bits are stored in the buffer
  //if the length of buffer reaches wordsize, empty buffer and this variable
  //push the buffer into vector<Word>
  BitString::Byte __dirty_bits;

  //vector<Word> for the actual storage of data.
  BitString &__data;
};

//Initalize the BitString element inside writer, if length % WORDSIZE != 0 that means
//the last word of BitString is not full, so store the tail into buffer and pop it from
//the vector<Word>
inline BitStringWriter::BitStringWriter(BitString &data)
    : __data(data) {
  if ((__dirty_bits = __data.length % BitString::WORDSIZE)) {
    __buffer = __data.data.back();     // back() operation returns a reference of the last element in vector
    __data.data.pop_back();
  }
  else {
    __buffer = 0UL;
  }
}

//Before deconstruction, if the last word is not fully used, just push it back to the
//vector<Word>
inline BitStringWriter::~BitStringWriter() {
  if (__dirty_bits)
    __data.data.push_back(__buffer);
}

//write one single bit to BitString, first enlarge the length of bitarray by 1,
//if buffer is full, push it into data array and empty dirtybits and buffer
//else
inline void BitStringWriter::writeBit(bool bit) {
  ++__data.length;
  //if buffer is full, empty it and push it back to data at first
  if (__dirty_bits == BitString::WORDSIZE) {
    __data.data.push_back(__buffer);
    __buffer = 0UL;
    __dirty_bits = 0;
  }
  ++__dirty_bits;
  //unless the bit is 1, the buffer should not be changed because original buffer is full of 0s.
  //here bits are stored from left to right.
  if (bit)
    __buffer |= 1UL << (BitString::WORDSIZE - __dirty_bits);
}

//write bits into BitString,
inline void BitStringWriter::writeBits(BitString::Byte num_bits, BitString::Word bits) {
  bits &= BitString::WORDMAX >> (BitString::WORDSIZE - num_bits);
  //std::cout << " HERE "<< (int)std::numeric_limits<BitString::Word>::max() << " , " << (int)BitString::WORDSIZE <<"," << (int)bits << std::endl;
  __data.length += num_bits;
  __dirty_bits += num_bits;
  //after inserting the new bits, if buffer is overwhelming, first fill the buffer with
  //corresponding bits and push it into data
  if (__dirty_bits > BitString::WORDSIZE) {
    __buffer |= bits >> (__dirty_bits - BitString::WORDSIZE);
    __data.data.push_back(__buffer);
    __buffer = 0UL;
    __dirty_bits -= BitString::WORDSIZE;
  }
  __buffer |= bits << (BitString::WORDSIZE - __dirty_bits);
}
}

#endif
