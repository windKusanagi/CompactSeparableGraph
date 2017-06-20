/*
 * NeighbourListNonSorted.h
 *
 *  Created on: Dec 28, 2016
 *      Author: xiangzhang
 */


#ifndef NEIGHBOURLISTNONSORTED_H_
#define NEIGHBOURLISTNONSORTED_H_



#include <cstdlib>
#include <iostream>
#include <string>
#include <algorithm>
#include <string>

using namespace std;

struct NeighbourEntry{
		int VertexId;

		NeighbourEntry* next;

};

class NeighbourListNonSorted{

public:
	NeighbourListNonSorted();

	static const int tableSize = 15606; // table size of 4elt graph
	//static const int tableSize = 214765; // table size of m14b graph
	//static const int tableSize = 143437; // table size of feocean graph
	//static const int tableSize = 15; // table size of hyper-small mesh graph
	//static const int tableSize = 144649;    //table size  of 144.graph
	//static const int tableSize = 110971; //table size of 598a graph
	//static const int tableSize = 28924; //table size of bcsstk30 graph
	//static const int tableSize = 35588;  //table size of bcsstk31 graph
	//static const int tableSize = 45087;  //table size of fe_body graph
	//static const int tableSize = 156317;  //table size of wave graph

	void addItem( int Vid, int index);
	void sortItself();
	void printTable();
	void printItemsInIndex(int index);
	int numberOfItemsInIndex(int index);


	int** sortedNeighbourTable;

	NeighbourEntry* Table[tableSize];

};
NeighbourListNonSorted::NeighbourListNonSorted(){
	for(int i=0 ; i<tableSize ; i++){
		Table[i] = new NeighbourEntry;
		Table[i]->VertexId = -1;
		Table[i]->next = NULL;
	}
	sortedNeighbourTable = NULL;
}

void NeighbourListNonSorted::addItem(int vid, int index){


	if(Table[index]->VertexId == -1){

		Table[index]->VertexId = vid;

	}else{
		NeighbourEntry* Ptr = Table[index];

		NeighbourEntry* n = new NeighbourEntry;
		n->VertexId = vid;
		n->next = NULL;

		while(Ptr->next != NULL){
			Ptr = Ptr->next;
		}
		Ptr->next = n;
	}
}

void NeighbourListNonSorted::printItemsInIndex(int index){
	NeighbourEntry* Ptr = Table[index];
	if(Ptr->VertexId == -1){
		cout<<"index = "<< index <<" is empty "<<endl;
	}else{
		cout<<"index "<< index << " contains the following items\n" ;
		while(Ptr != NULL){
			cout << "------------\n";
			cout <<Ptr->VertexId <<endl;

			Ptr = Ptr->next;
		}
	}
}
int NeighbourListNonSorted::numberOfItemsInIndex(int index){
	int count = 0;
	if(Table[index]->VertexId == -1){
		return count;
	}else{
		count++;
		NeighbourEntry* Ptr = Table[index];
		while(Ptr->next != NULL){
			count++;
			Ptr = Ptr->next;
		}
	}
	return count;
}

void NeighbourListNonSorted::sortItself(){

	sortedNeighbourTable = new int*[tableSize];

	for(int i=0; i<tableSize ; i++){
		sortedNeighbourTable[i] = new int [ numberOfItemsInIndex(i)];
	}

	for(int i=0; i<tableSize ; i++){

		int tempArr [numberOfItemsInIndex(i)];
		NeighbourEntry* Ptr = Table[i];
		for(int j=0; j<numberOfItemsInIndex(i); j++){
			//cout << Ptr->VertexId<< " * " ;
			tempArr [j] = Ptr->VertexId;
			if(Ptr->next != NULL){

				Ptr = Ptr->next;
			}
		}
		sort(tempArr, tempArr + numberOfItemsInIndex(i) );
		for (int k=0; k<numberOfItemsInIndex(i) ; k++){
			sortedNeighbourTable[i][k] = tempArr[k];
		}
	}

}
#endif
