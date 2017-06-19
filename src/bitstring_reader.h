/*
 * bitstring_reader.hpp
 *
 *  Created on: Feb 25, 2017
 *      Author: xiangzhang
 */

#ifndef BITSTRING_READER_HPP
#define BITSTRING_READER_HPP

#include <algorithm>
#include <iostream>
#include "bitstring.h"

namespace phsim {

class BitStringReader {

 public:
  explicit BitStringReader(BitString &);
  explicit BitStringReader(BitString &, int, int , int);     // BitStringReader constructor to read data from a specific start position by n length

  bool eof();
  BitString::Size bitsLeft();

  bool readBit();
  BitString::Word readBits(BitString::Byte);

 private:
  void fillBuffer();

  BitString::Word __buffer;
  BitString::Word __buffered_bits;
  BitString::Data::iterator __data;
  BitString::Size __bits_left;
};

inline BitStringReader::BitStringReader(BitString &data)
    : __buffer(0), __buffered_bits(0),
      __data(data.data.begin()),
      __bits_left(data.length) {
			//std::cout<< "data.length is"<<data.length<<std::endl;
	}

inline BitStringReader::BitStringReader(BitString &data, int vectorIndex, int length, int vOffset)
    : __buffer(0), __buffered_bits(0),
      __data(data.data.begin()+=vectorIndex),
      __bits_left(length + vOffset) {
			//std::cout<< "data.length is: "<<data.length<<std::endl;
	}
//if no bits left, then it is the end of file
inline bool BitStringReader::eof() {
  return (__bits_left == 0);
}

//return how many bits are left
inline BitString::Size BitStringReader::bitsLeft() {
  return __bits_left;
}

inline bool BitStringReader::readBit() {
  fillBuffer();
  --__bits_left;
  return __buffer & (1UL << --__buffered_bits);//let 1 left shift the number of buffered bits and make and operation with the buffer
};

inline BitString::Word BitStringReader::readBits(BitString::Byte num_bits) {
  fillBuffer();
  //std::cout<<"inside readBits:__buffered_bits: "<<__buffered_bits<<" num_bits: "<<(int)num_bits<<" __bits_left:"<<__bits_left<<std::endl;
  if (__buffered_bits >= num_bits) {//if the num of buffered bits is larger than required bits.
    __bits_left -= num_bits;        //bits left should reduce corresponding number
	//std::cout<<"inside readBits:__bits_left:"<<__bits_left<<std::endl;
    __buffered_bits -= num_bits;    //buffered bits are consumed.
    return (__buffer >> __buffered_bits) & (BitString::WORDMAX >> (BitString::WORDSIZE - num_bits));
  }
  else {//if cached bits is not enough for the requirement
    BitString::Word high_bits = __buffer & (BitString::WORDMAX >> (BitString::WORDSIZE - __buffered_bits));
    num_bits -= __buffered_bits;
    __bits_left -= __buffered_bits;
	//std::cout<<"__bits_left in second half "<<__bits_left<<std::endl;
    __buffered_bits = 0;
    fillBuffer();
	//std::cout<<"__bits_left after fillbuffer "<<__bits_left<<std::endl;
	//std::cout<<"__bits_left - num_bits "<<__bits_left<<" "<<(int)num_bits<<std::endl;
    __bits_left -= num_bits;
	//std::cout<<"__bits_left after -=num_bits "<<__bits_left<<std::endl;
    __buffered_bits -= num_bits;
    return (high_bits << num_bits)
        | ((__buffer >> __buffered_bits) & (BitString::WORDMAX >> (BitString::WORDSIZE - num_bits)));
  }
};

inline void BitStringReader::fillBuffer() {
  //if buffered bits is empty
  //std::cout<<"inside fillBuffer: bufferedbits:"<<__buffered_bits<<std::endl;
  if (!__buffered_bits) {
    __buffer = *(__data++);   //current buffer is a memeber of the data vector
    __buffered_bits = std::min<BitString::Size>(BitString::WORDSIZE, __bits_left);//return smaller values between WORDSIZE and bits left.
	//std::cout<<"|new "<<__buffered_bits<<" |"<<std::endl;
    __buffer >>= BitString::WORDSIZE - __buffered_bits;//right shift the num of bits that are not buffered.
  }
};
}

#endif
