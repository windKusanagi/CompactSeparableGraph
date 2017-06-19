/*
 * IndexStructure.h
 *
 *  Created on: Apr 13, 2017
 *      Author: xiangzhang
 */


#ifndef INDEXSTRUCTURE_H
#define INDEXSTRUCTURE_H

#include <iostream>
#include <vector>
#include <string>
#include <bitset>
#include <math.h>
#include <sdsl/bit_vectors.hpp>
using namespace sdsl;

class IndexStructure{


public:
	IndexStructure(unsigned int nV);
	void countDITotalBits();
	void buildSemiDIndex();
	int integerBitLength(int num);
	void buildIndirectIndex ();
	void buildrrrDI(int);
	void buildsdbDI(int);

	// Direct Index Class Member
	std::vector <unsigned int> directIndex;
	unsigned long directIndexTotalBit;
	unsigned int vNum;



	// Semi DirectIndex Class Member
	unsigned int semiDIflag;
	std::vector <int> semiDIFW;
	std::vector < std::vector <std::bitset<10>  >  > semiDISW;
	std::vector <unsigned int * > semiDISWSub;


	// IndirectIndex Class Member
	unsigned blockNum;
	unsigned blockVNum;
	unsigned subBlockLength;
	unsigned subBlockCounter;

	//std::vector <std::bitset<14> > bitVofEachBlock;  // Number of vertex in each block on 4elt.graph
	//std::vector <std::bitset<18> > bitVofEachBlock;	 // Number of vertex in each block on feocean.graph
	//std::vector <std::bitset<18> > bitVofEachBlock;	 // Number of vertex in each block on m14b.graph

	std::vector < bit_vector > bitVofEachBlock;  // Number of vertex in each block on 4elt.graph

	//std::vector < std::vector< std::vector<unsigned int>  >  > indrectIndex;
	std::vector < std::vector<unsigned int>   > indrectIndex;
	std::vector < std::vector< unsigned int >  >subBlockVNum;
	std::vector < rrr_vector<> > rrrBVofEachBlock;
	rrr_vector<> rrr_direct_index;
	sd_vector<> sdb_vector;
};


IndexStructure::IndexStructure( unsigned int nV ){

	this->directIndexTotalBit = 0;
	this->semiDIflag = 0;
	this->vNum = nV;
	this->blockVNum=0;
	this->blockNum=0;
	this->subBlockLength = 0;
	this->subBlockCounter = 0;

}

void IndexStructure::countDITotalBits(){
	unsigned long sum=0;
	for (int i=0; i<this->directIndex.size() ;i++){
		sum+= integerBitLength(this->directIndex[i]);
	}
	this->directIndexTotalBit = sum;
}

void IndexStructure::buildSemiDIndex(){
	std::vector <int> tempOffsetSet;
	for (int i=0; i<this->directIndex.size();i++){
		if (i%4 == 0){
			this->semiDIFW.push_back(this->directIndex[i]);
		}else{
			tempOffsetSet.push_back(this->directIndex[i]);
		}
	}
	for (int i=0 ; i<tempOffsetSet.size()/3; i++){

		int offset1 = tempOffsetSet[i*3]   - this->semiDIFW[i];
		int offset2 = tempOffsetSet[i*3+1] - this->semiDIFW[i];
		int offset3 = tempOffsetSet[i*3+2] - this->semiDIFW[i];
		int len1 = integerBitLength( offset1 );
		int len2 = integerBitLength( offset2 );
		int len3 = integerBitLength( offset3 );

		std::vector< bitset<10> > secWord;
		if ( len3 <= 11){
			std::bitset<10> n1(offset1);
			std::bitset<10> n2(offset2);
			std::bitset<10> n3(offset3);
			secWord.push_back(n1);
			secWord.push_back(n2);
			secWord.push_back(n3);
			this->semiDISW.push_back(secWord);
			this->semiDIflag = (i+1)*4;
		}else{
			unsigned int secSubWord[3] = {tempOffsetSet[i*3], tempOffsetSet[i*3+1],tempOffsetSet[i*3+2]};
			this->semiDISWSub.push_back(secSubWord);
		}
	}
	int rest = tempOffsetSet.size()%3 ;
	unsigned int *tail = new unsigned int [rest];
	for (int i=0; i<rest; i++){
		tail[i] = tempOffsetSet[tempOffsetSet.size()/3*3 + i];
	}
	this->semiDISWSub.push_back(tail);


}
int IndexStructure::integerBitLength(int num) {
	int length = 0;
	while (num > 0)
		num >>= 1, ++length;
	return length;
}
void IndexStructure::buildIndirectIndex(){
	this->blockVNum = std::ceil(std::log2(this->vNum));
	this->blockNum = std::ceil(this->vNum*1.0/this->blockVNum);
	this->subBlockLength = 16* this->blockVNum;

	std::vector <unsigned int> vIndex;   // 	only need to store the address of vertex which is indicated as 1 in bitVector of each block
	std::vector<unsigned int>  vNumInEachSubBlock;
	//std::bitset<14> bitVPerBlock;   // 4elt.graph
	//std::bitset<18> bitVPerBlock;	 // feocean.graph
	//std::bitset<18> bitVPerBlock;	 // m14b.graph
	bit_vector bitVPerBlock = bit_vector(this->blockVNum, 0);

	for (int i = 0; i < this->blockNum -1 ; i++) {
	//for (int i = 0; i < 3 ; i++) {
		//cout << "------------------------------------------------------------------------------------"<<endl;
		unsigned int offsetCountforSubblock = 0;
		unsigned int vertexCounter = 0;
		//bitVPerBlock.reset();
		//bitVPerBlock.set(0);
		bitVPerBlock = bit_vector(this->blockVNum,0);
		bitVPerBlock[0] = true;
		vIndex.push_back(this->directIndex[i * this->blockVNum]);
		unsigned int lastPos = this->directIndex[i * this->blockVNum];
//		if (i == 92){
//			cout << "Block 92 :  lastPos in  "<< i * this->blockVNum << " index: "  << this->directIndex[i * this->blockVNum]<<endl;;
//		}
		for (int j = 0; j < this->blockVNum; j++) {
			//cout << "j : " << j <<endl;

			vertexCounter += 1;
			offsetCountforSubblock = this->directIndex[i * this->blockVNum + j] - lastPos;
//			if(i==92){
//				cout << "current vertex index : " << this->directIndex[i * this->blockVNum + j] <<endl;
//				cout << "current offset : " << offsetCountforSubblock <<endl;;
//			}
			if (j < this->blockVNum-1){
				//cout << offsetCountforSubblock << "     " << lastPos << endl;
				if (offsetCountforSubblock >= this->subBlockLength  ) {
					offsetCountforSubblock = 0;
					vNumInEachSubBlock.push_back(vertexCounter);
					vertexCounter = 0;
					//bitVPerBlock.set(j+1);
					bitVPerBlock[j+1] = true;
					vIndex.push_back(this->directIndex[i * this->blockVNum + j+1]);
					lastPos = this->directIndex[i * this->blockVNum + j+1];
					this->subBlockCounter += 1;
				}
			}else{
				if (bitVPerBlock[bitVPerBlock.size()-1] == true){

					vNumInEachSubBlock.push_back(1);
					this->subBlockCounter += 1;
				}else{
					vNumInEachSubBlock.push_back(vertexCounter);
					this->subBlockCounter += 1;
				}
			}
		}
//		if (i == 92){
//			cout << bitVPerBlock<<endl;
//		}
		this->subBlockVNum.push_back(vNumInEachSubBlock);
		rrr_vector<> rrr_bv = bitVPerBlock;
		this->rrrBVofEachBlock.push_back(rrr_bv);
		this->bitVofEachBlock.push_back(bitVPerBlock);
		//cout << "ccc   " <<  vIndex.size() <<endl;
		//cout << vIndex.size()<<endl;
		this->indrectIndex.push_back(vIndex);
		vIndex.clear();
	}

}
void IndexStructure::buildrrrDI( int codeLength){
	bit_vector b = bit_vector(codeLength, 0);
	for (int i=0; i<this->vNum; i++){
		b[this->directIndex[i]]=true;
	}
	this->rrr_direct_index = b;
}
void IndexStructure::buildsdbDI( int codeLength){
	bit_vector b = bit_vector(codeLength, 0);
	for (int i=0; i<this->vNum; i++){
		b[this->directIndex[i]]=true;
	}
	//this->sdb_vector(b);
	sd_vector<> s_v(b);
	this->sdb_vector = s_v;
}

#endif
