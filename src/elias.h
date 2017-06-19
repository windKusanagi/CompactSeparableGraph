/*
 * elias.hpp
 *
 *  Created on: Feb 25, 2017
 *      Author: xiangzhang
 */


#ifndef ELIAS_HPP
#define ELIAS_HPP

#include <iostream>

#include "bitstring_reader.h"
#include "bitstring_writer.h"

namespace phsim {


namespace {

BitString::Byte bitLength(BitString::Word word) {
  BitString::Byte length = 0;
  while (word > 0)
    word >>= 1, ++length;
  //std::cout<< (int)length <<std::endl;;
  return length;
}

}

inline void gammaEncode(BitStringWriter &writer, BitString::Word word , unsigned int &curBitNum) {
  BitString::Byte num_bits = bitLength(word);
  writer.writeBits(num_bits - 1, 0UL);
  writer.writeBits(num_bits, word);
  curBitNum += 2* num_bits -1;
}

inline BitString::Word gammaDecode(BitStringReader &reader) {
  BitString::Byte num_bits = 0;
  while (!reader.readBit()){
    ++num_bits;
  }
  //std::cout<<std::endl;
  //std::cout<<"in gammaDecode required bits is "<<(int)num_bits<<std::endl;
  if (num_bits)
    return (1UL << num_bits) | reader.readBits(num_bits);
  else
    return 1UL;
}

inline  BitString::Word gammaDecodeWithOffset(BitStringReader &reader, int offset) {
  BitString::Byte num_bits = 0;
  for (int i=0; i<offset; i++) reader.readBit();
  while (!reader.readBit()){
    ++num_bits;
  }
  //std::cout<<"in gammaDecode required bits is "<<(int)num_bits<<std::endl<<std::endl;
  if (num_bits)
    return (1UL << num_bits) | reader.readBits(num_bits);
  else
    return 1UL;
}

inline void omegaEncode(BitStringWriter &writer, BitString::Word word) {
  BitString::Byte result_bits = 1;
  BitString::Word high_bits = 0;
  BitString::Word low_bits = 0;
  while (word > 1) {
    BitString::Byte num_bits = bitLength(word);
    if (result_bits >= BitString::WORDSIZE) {
      high_bits |= word << (result_bits - BitString::WORDSIZE);
      result_bits += num_bits;
    }
    else {
      low_bits |= word << result_bits;
      result_bits += num_bits;
      if (result_bits > BitString::WORDSIZE)
        high_bits = word >> (BitString::WORDSIZE + num_bits - result_bits);
    }
    word = num_bits - 1;
  }
  if (result_bits > BitString::WORDSIZE)
    writer.writeBits(result_bits - BitString::WORDSIZE, high_bits);
  writer.writeBits(std::min(result_bits, BitString::WORDSIZE), low_bits);
}

inline BitString::Word omegaDecode(BitStringReader &reader) {
  BitString::Word result = 1;
  while (reader.readBit()) {
    result = (1 << result) | reader.readBits(result);
  }
  return result;
}

inline void deltaEncode(BitStringWriter &writer, BitString::Word word) {
  BitString::Word num_bits = bitLength(word);
  BitString::Word num_num_bits = bitLength(num_bits);
  writer.writeBits(num_num_bits - 1, 0);
  writer.writeBits(num_num_bits, num_bits);
  writer.writeBits(num_bits - 1, word);
}

inline BitString::Word deltaDecode(BitStringReader &reader) {
  BitString::Byte num_bits = 0;
  while (!reader.readBit())
    ++num_bits;
  if (num_bits) {
    num_bits = ((1 << num_bits) | reader.readBits(num_bits)) - 1UL;
    return (1UL << num_bits) | reader.readBits(num_bits);
  }
  else {
    return 1UL;
  }
}
}

#endif
