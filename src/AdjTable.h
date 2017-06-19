/*
 * AdjTable.h
 *
 *  Created on: Dec 28, 2016
 *      Author: xiangzhang
 */
#ifndef ADJTABLE_H_
#define ADJTABLE_H_

#include <cstdlib>
#include <iostream>
#include <string>

using namespace std;

struct AdjEntry{
		std::string GammaCode;

		AdjEntry* next;

};

class AdjTable{

public:
	AdjTable();

	void AddItem(std::string Gcode, int index);

	void PrintTable();
	static const int tableSize = 10;

	AdjEntry* Table[tableSize];

};

AdjTable::AdjTable(){
	for(int i=0 ; i<tableSize ; i++){
		Table[i] = new AdjEntry;
		Table[i]->GammaCode = "empty";
		//HashTable[i]->drink = "empty";
		Table[i]->next = NULL;
	}

}

void AdjTable::AddItem(string Gcode, int index){


	if(Table[index]->GammaCode == "empty"){

		Table[index]->GammaCode = Gcode;

	}else{
		AdjEntry* Ptr = Table[index];

		AdjEntry* n = new AdjEntry;
		n->GammaCode = Gcode;
		n->next = NULL;

		while(Ptr->next != NULL){
			Ptr = Ptr->next;
		}
		Ptr->next = n;
	}
}


#endif



