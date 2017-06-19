/*
 * decodeLookupTable.h
 *
 *  Created on: Apr 28, 2017
 *      Author: xiangzhang
 */

#ifndef DECODELOOKUPTABLE_H_
#define DECODELOOKUPTABLE_H_

#include <iostream>
#include <bitset>
#include <vector>
#include <math.h>
#include <string>

class decodeLookupTable{

public:
	decodeLookupTable();
	std::bitset<32> BuildTableEntry ( std::bitset<16> );
	void SetNumOfBitsForNext (std::bitset<32> &, int );

	std::vector< std::bitset<32> > lookupTable;  // Component of each enrty in this table:
												 // From left to right ( reverse to the natural index of a bit vector )
												 // bit index 0-4 (5 bits) to indicate how many bits are left for the next 16 bits
												 // bit index 5-6 (2 bits) to indicate the number of encoded numbers can be decoded in this 16 bits
												 // bit index 7-30 (24 bits) to store the decoded numbers respectively, each of them takes 8-bits space
												 // bit index 31 (1 bit) meaningless
};

decodeLookupTable::decodeLookupTable(){

	for ( int i=0; i< pow(2,16) ; i++){
		std::bitset<16> bs ( pow(2,10)-i);
		std::bitset<32> rst = BuildTableEntry(bs);
		this->lookupTable.push_back(rst);
	}
}

std::bitset<32> decodeLookupTable::BuildTableEntry ( std::bitset<16> bs0){

	std::bitset<32> rst;

	int counter = 0;
	int pos = 15;
	int rstPos = 17;

	int oneCounter = 0;
	for( int i = 0; i <15 ;i++){
		if (bs0[pos -i] == true){
			oneCounter ++ ;
		}
	}
	if ( oneCounter > 8){
		return rst;
	}

	for( int i = 0; i< 15; i++){
		if (bs0[pos - i] == false){
			counter ++;
		}else{
			break;
		}
	}

	if ( counter >= 8){   // All the bits in this 16 bits in part of the next gamma code, rst [31 : 27 ] = "10000"
		rst.set(31);

	}else{

		rst.set(25); // Indicate there are at least one number can be decoded in this 16 bits

		pos = pos -(2*counter);
		//cout << pos << endl;
		for( int i = 0; i < counter+1 ; i++){
			if (bs0[pos+i] == true){
				rst.set(rstPos+i);
			}
		}
		pos -= 1;
		int temp = counter;
		counter = 0;
		rstPos -= 8;

		if (bs0[pos] == true){
			rst.reset();
			return rst;
		}

		for( int i = 0; i< pos; i++){
			if (bs0[pos - i] == false){
				counter ++;
			}else{
				break;
			}
		}

		if ( pos < 2*counter +1 ){
			SetNumOfBitsForNext ( rst, pos+1);   // Set number of bits which is part of next bits

		}else{
			rst.reset(25); // Indicate there are at least two numbers can be decoded in this 16 bits
			rst.set(26);

			pos = pos -(2*counter);
			//cout << pos << endl;
			for( int i = 0; i < counter+1 ; i++){
				if (bs0[pos+i] == true){
					rst.set(rstPos+i);
				}
			}
			pos -= 1;

			if (bs0[pos] == true){
				rst.reset();
				return rst;
			}

			temp = counter;
			counter = 0;
			rstPos -= 8;
			for( int i = 0; i< pos; i++){
				if (bs0[pos - i] == false){
					counter ++;
				}else{
					break;
				}
			}
			if ( pos < 2*counter +1 ){
				SetNumOfBitsForNext ( rst, pos+1);   // Set number of bits which is part of next bits

			}else{
				rst.set(25); // Indicate there are three numbers can be decoded in this 16 bits
				pos = pos -(2*counter);
				for( int i = 0; i < counter+1 ; i++){
					if (bs0[pos+i] == true){
						rst.set(rstPos+i);
					}
				}
				SetNumOfBitsForNext ( rst, pos);

			}
		}
	}
	return rst;

}

void decodeLookupTable::SetNumOfBitsForNext (std::bitset<32> &bits, int pos){
	if (pos==1){
		bits.set(27);
	}
	if (pos==2){
		bits.set(28);
	}
	if (pos==3){
		bits.set(27);
		bits.set(28);
	}
	if (pos==4){
		bits.set(29);
	}
	if (pos==5){
		bits.set(27);
		bits.set(29);
	}
	if (pos==6){
		bits.set(28);
		bits.set(29);
	}
	if (pos==7){
		bits.set(27);
		bits.set(28);
		bits.set(29);
	}
	if (pos==8){
		bits.set(30);
	}
	if (pos==9){
		bits.set(27);
		bits.set(30);
	}
	if (pos==10){
		bits.set(28);
		bits.set(30);
	}
	if (pos==11){
		bits.set(27);
		bits.set(28);
		bits.set(30);
	}
	if (pos==12){
		bits.set(29);
		bits.set(30);
	}
	if (pos==13){
		bits.set(27);
		bits.set(29);
		bits.set(30);
	}
}


#endif /* DECODELOOKUPTABLE_H_ */
