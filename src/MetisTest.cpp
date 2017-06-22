/*
 * MetisTest.cpp
 *
 *  Created on: Nov 28, 2016
 *      Author: xiangzhang
 */

#include <cstddef>
#include <iostream>
#include <metis.h>
#include <sstream>
#include <vector>
#include <set>
#include <math.h>
#include <string>
#include <algorithm>
#include <fstream>
#include <sdsl/bit_vectors.hpp>
#include <chrono>
#include <ctime>

#include "MetisNode.h"
#include "AdjTable.h"
#include "bitstring.h"
#include "NeighbourListNonSorted.h"
#include "bitstring_reader.h"
#include "bitstring_writer.h"
#include "elias.h"
#include "IndexStructure.h"

using namespace std;
using namespace sdsl;

struct BstNode {
	vector <string> EdgeSparators;
	BstNode* left;
	BstNode* right;
	int ** mappingOrg;
	int VertexId;
	int NewVertexId;
	int flag; // "1" is internal node and "0" is leaf node.
	int vertexNum; // number of vertex under the current node
	//set <int> leaves; // set of leaves under current node
	int * ordered_array;
	int array_len;
	//set<int> N1;  // used in child flipping section
	//set<int> N2;

};
// Function to create a new Node in heap
BstNode* getNewNode(vector <string> EdgeSparator , int num) {
	BstNode* newNode = new BstNode();
	newNode->EdgeSparators = EdgeSparator;
	newNode->vertexNum = num;
	newNode->left = newNode->right = NULL;
	newNode->mappingOrg = NULL;
	newNode->flag = 1;
	newNode->VertexId = -1;
	newNode->NewVertexId = -1;
	newNode->ordered_array = NULL;
	newNode->array_len = 0;
	return newNode;
}

BstNode* getLeafNode( int id){
	BstNode* newNode = new BstNode();
	newNode->flag = 0;
	newNode->left = newNode->right = NULL;
	newNode->VertexId = id;
	newNode->NewVertexId = -1;
	newNode->mappingOrg = NULL;
	newNode->vertexNum = 1;
	//newNode->leaves.insert(id);
	newNode->ordered_array = new int[1];
	newNode->ordered_array[0]  = id;
	newNode->array_len = 1;
	return newNode;
}

BstNode * root = NULL;
int InorderTraversalVID = 0;
int *mapAfterRenumber = NULL;
int *mappingForRestore = NULL;
vector< vector<string> > adjTable;  //Store the gamma code of
NeighbourListNonSorted nlns; // initialize the neighbour list without sorting


vector< vector<phsim::BitString::Word> > inputWV;
vector< phsim::BitString::Word > directDecodeResult;
vector<int> offsetIndex;

//idx_t nV = 15;

idx_t nV = 15606;  // nV of 4elt.graph
idx_t nE = 45878;  // nE of 4elt.graph

//idx_t nV = 214765;  // nV of m14b.graph
//idx_t nE = 1679018;  // nE of m14b.graph

//idx_t nV = 143437;  // nV of feocean.graph
//idx_t nE = 819186;  // 2*nE of feocean.graph


//idx_t nV = 144649;    // nV of 144.graph
//idx_t nE = 1074393;   // nE of 144.graph

//idx_t nV = 110971;    // nV of 598a.graph
//idx_t nE = 741934;   // nE of 598a.graph

//idx_t nV = 28924;    // nV of bcsstk30.graph
//idx_t nE = 1007284;   // nE of bcsstk30.graph

//idx_t nV = 35588;    // nV of bcsstk31.graph
//idx_t nE = 572914;   // nE of bcsstk31.graph

//idx_t nV = 45087;    // nV of fe_body.graph
//idx_t nE = 163734;   // nE of fe_body.graph

//idx_t nV = 156317;    // nV of wave.graph
//idx_t nE = 1059331;   // nE of wave.graph


idx_t * Inputadj = new int [91756];  // size of adj of 4elt graph
idx_t * InputXadj = new int [15607]; // nV of 4elt + 1

//idx_t * Inputadj = new int [3358036];  // size of adj of m14b graph
//idx_t * InputXadj = new int [214766]; // nV of m14b + 1

//idx_t * Inputadj = new int [819186];  // size of adj of feocean graph
//idx_t * InputXadj = new int [143438]; // nV of feocean + 1

//idx_t * Inputadj = new int [2148786];  // size of adj of 144 graph
//idx_t * InputXadj = new int [144650]; // nV of 144 + 1

//idx_t * Inputadj = new int [1483868];  // size of adj of 598a graph
//idx_t * InputXadj = new int [110972]; // nV of 598a + 1

//idx_t * Inputadj = new int [2014568];  // size of adj of bcsstk30 graph
//idx_t * InputXadj = new int [28925]; // nV of bcsstk30 + 1

//idx_t * Inputadj = new int [1145828];  // size of adj of bcsstk31 graph
//idx_t * InputXadj = new int [35589]; // nV of bcsstk31 + 1

//idx_t * Inputadj = new int [327468];  // size of adj of fe_body graph
//idx_t * InputXadj = new int [45088]; // nV of fe_body + 1

//idx_t * Inputadj = new int [2118662];  // size of adj of wave graph
//idx_t * InputXadj = new int [156318]; // nV of wave + 1



//idx_t nVertices = 15;
//idx_t nEdges    = 22;
//idx_t InputXadj[16] = {0, 2, 5, 8, 11, 13, 16, 20, 24, 28,31, 33, 36, 39, 42, 44};
//idx_t Inputadj[44] = {1,5,0, 2, 6, 1, 3, 7, 2, 4, 8, 3, 9, 0, 6, 10, 1, 5, 7, 11, 2, 6, 8, 12, 3, 7, 9, 13, 4, 8, 14, 5,
//			11, 6, 10, 12, 7, 11, 13, 8, 12, 14, 9, 13};


bool *negativeDiff = new bool [nV];

IndexStructure indexStructure (nV) ; // initialize the index structure which includes direct, semi-direct and indirect index structure

int metisUsage();
string IntToString1 (int a);
BstNode * generateTree ( idx_t nVertex , idx_t *xadj, idx_t *adj , int ** mapOrg);
BstNode * generateLeaves (idx_t nVertex , idx_t *xadj, idx_t *adj, int ** mapOrg);
void inorderTraversal(BstNode* x);
void preorderTraversal(BstNode* x);
void renumberByInorderTraversal(BstNode* x);
void findNeighbour(BstNode * x);
void formAdjTable(NeighbourListNonSorted x );
void readInput();

void postOrderTraCountLeaf(BstNode* x);
void preorderTraversalCheckLeaf(BstNode* x);
void childFlipping (BstNode* node, vector< int* > &, vector< int*>& , vector<int>&, vector<int>&);
//void childFlipping1( BstNode* node, vector< set<int> > & , vector< set<int> > &);
void directDecodeAll(phsim::BitString );
int* mergeArray(int *, int *, int, int);

vector <int> getAdjForOneVUsingDI (phsim::BitString, unsigned int);
vector <int> getAdjForOneVUsingSemiDI (phsim::BitString, unsigned int);
vector <int> getAdjForOneVUsingrrrDI (phsim::BitString, unsigned int);
vector <int> getAdjForOneVUsingIDI (phsim::BitString, unsigned int);
vector <int> getAdjForOneVUsingsdbDI (phsim::BitString, unsigned int);

vector <int> getAdjForOneVUsingDIWithD (phsim::BitString, unsigned int);
vector <int> getAdjForOneVUsingSemiDIWithD (phsim::BitString, unsigned int);
vector <int> getAdjForOneVUsingrrrDIWithD (phsim::BitString, unsigned int);
vector <int> getAdjForOneVUsingsdbDIWithD (phsim::BitString, unsigned int);

int getDegreeUsingDI(phsim::BitString, unsigned int, bool);
int getDegreeUsingSemiDI(phsim::BitString, unsigned int, bool);
int getDegreeUsingIDI (phsim::BitString, unsigned int);
int getDegreeUsingrrrDI (phsim::BitString, unsigned int, bool);
int getDegreeUsingsdbDI (phsim::BitString, unsigned int, bool);

int compareints (const void * a, const void * b){
  return ( *(int*)a - *(int*)b );
}
void testUnit (phsim::BitString);
void testUnitWithD (phsim::BitString, int );
void DFS(phsim::BitString);


int main() {



	readInput();
	cout << "Read input completed" <<endl;
    int ** FirstRound = NULL;
    mapAfterRenumber = new int[nV];
    mappingForRestore = new int[nV];

    std::chrono::time_point<std::chrono::system_clock> start1, end1;
    start1 = std::chrono::system_clock::now();

    //generateTree(nVertices , InputXadj , Inputadj , FirstRound );
    generateTree(nV , InputXadj , Inputadj , FirstRound );
    cout << "Generating separator tree completed" <<endl;
    //preorderTraversal(root);

    end1 = std::chrono::system_clock::now();

    std::chrono::duration<double> elapsed_seconds = end1-start1;
    std::time_t end_time1 = std::chrono::system_clock::to_time_t(end1);

    std::cout << "Generating tree structure finished computation at " << std::ctime(&end_time1)
              << "elapsed time: " << elapsed_seconds.count() << "s\n";

    postOrderTraCountLeaf(root);
    cout << "Making leaf set completed !" <<endl;
    //preorderTraversal(root);

    //preorderTraversalCheckLeaf(root);

    std::chrono::time_point<std::chrono::system_clock> start2, end2;
    start2 = std::chrono::system_clock::now();

    vector<int*> n1;
    vector<int*> n2;
    vector<int > n1size;
    vector<int > n2size;
    //buildN1N2(root);
    childFlipping(root , n1, n2 , n1size , n2size);

    cout << "Flipping child completed!" <<endl;

    end2 = std::chrono::system_clock::now();

    std::chrono::duration<double> elapsed_seconds2 = end2-start2;
    std::time_t end_time2 = std::chrono::system_clock::to_time_t(end2);

    std::cout << "Child flipping finished computation at " << std::ctime(&end_time2)
              << "elapsed time: " << elapsed_seconds2.count() << "s\n";



    renumberByInorderTraversal(root);
    cout << "Renumber completed ! " <<  endl;

    findNeighbour(root);
    nlns.sortItself();
    formAdjTable(nlns);
    cout<< "Forming adj table completed ! " <<endl;


    //preorderTraversal(root);

    // add the number of neighbors for each vertex at the beginning of each adj
    //vector< int > degree;


    phsim::BitString data;
	{
		phsim::BitStringWriter writer(data);
		unsigned int curBitN = 0;
		for (int i=0; i< inputWV.size() ; i++){
			indexStructure.directIndex.push_back(curBitN);
			//degree.push_back(inputWV[i].size());
			phsim::gammaEncode(writer ,inputWV[i].size() ,curBitN );
			for(auto x:inputWV[i]){
				phsim::gammaEncode(writer ,x ,curBitN );
			}

		}
	}


//	int degreeBitCount = 0;
//	for ( int i=0; i<degree.size();i++){
//		degreeBitCount+= 2* floor(log2(degree[i])) ;
//		degreeBitCount ++;
//	}
//
//	cout << degreeBitCount <<endl;

	// without the number of neighbors for each vertex at the beginning of each adj

//    phsim::BitString data;
//	{
//		phsim::BitStringWriter writer(data);
//		unsigned int curBitN = 0;
//		for (int i=0; i< inputWV.size() ; i++){
//			indexStructure.directIndex.push_back(curBitN);
//			for(auto x:inputWV[i]){
//				phsim::gammaEncode(writer ,x ,curBitN );
//			}
//
//		}
//	}

    cout<< "Gamma Encoding completed ! " <<endl;
    cout<< data.data.size() <<endl;



    //cout << degreeBitCount<<endl;;
    //indexStructure.countDITotalBits();
    cout << "Building Direct Index completed !" << endl;
//    //cout << index.directIndexTotalBit << endl;
//
    indexStructure.buildSemiDIndex();
    cout << "Building Semi Direct Index completed ! " <<endl;

    //cout << index.semiDIflag <<endl;
    //cout << "first word size  : " << index.semiDIFW.size() << endl;
    //cout << "second word fit size  :  " << index.semiDISW.size() <<endl;
    //cout << "second word not fit size  : " << index.semiDISWSub.size() <<endl;


    indexStructure.buildIndirectIndex();
    cout << "Building Indirect Index completed ! " <<endl;
    indexStructure.buildrrrDI(data.length);
    cout << "Building rrr direct Index completed ! " <<endl;
    indexStructure.buildsdbDI(data.length);
    cout << "Building sbd direct Index completed ! " <<endl;



    //DFS(data);

//    directDecodeAll(data);
//    cout << "Directly decode the dataset completed! " <<endl;

//    cout << indexStructure.subBlockCounter <<endl;
//    cout << indexStructure.blockVNum << endl;
//    cout << indexStructure.blockNum << endl;


//    cout << size_in_bytes(indexStructure.rrr_direct_index) <<endl;
//    cout << size_in_bytes(indexStructure.sdb_vector) <<endl;
//
//


//    cout << " vid '1401' :" <<mappingForRestore[1401]<<endl;
//    cout << "Using directIndex:" <<endl;
//    vector<int> test = getAdjForOneVUsingDIWithD(data , 1401);
//    for (int i=0; i<test.size(); i++){
//    	cout << test[i] << endl;
//    }
//

//    cout << "-----------------------------------------------------------------------" <<endl;
//    cout << " vid '1401' :" <<mappingForRestore[1401]<<endl;
//    cout << "Using directIndex:" <<endl;
//    vector<int> test = getAdjForOneVUsingDI(data , 1401);
//    for (int i=0; i<test.size(); i++){
//    	cout << test[i] << endl;
//    }
    // 35588
//    cout << "Using SemiDIndex :" <<endl;
//    for ( int i=0 ;i< nV-1 ;i++){
//    	getAdjForOneVUsingSemiDI(data , i);
//    }
//    cout << "complete" <<endl;
//    vector<int> test1 = getAdjForOneVUsingSemiDI(data , 1401);
//	for (int i=0; i<test1.size(); i++){
//		cout << test1[i] << endl;
//	}
//    cout << "Using rrr DIndex :" <<endl;
//    vector<int> test2 = getAdjForOneVUsingrrrDI(data , 1401);
//	for (int i=0; i<test2.size(); i++){
//		cout << test2[i] << endl;
//	}
//    cout << "Using sdb DIndex :" <<endl;
//    vector<int> test3 = getAdjForOneVUsingsdbDI(data , 1401);
//	for (int i=0; i<test3.size(); i++){
//		cout << test3[i] << endl;
//	}
//    cout << " vid '4' :" <<mappingForRestore[4] <<endl;
//    cout << "Using indirect DIndex :"<<endl;
//    vector<int> test4 = getAdjForOneVUsingIDI (data, 1401);
//	for (int i=0; i<test4.size(); i++){
//		cout << test4[i] << endl;
//	}
//	cout << "----------------------------------------------------------------------" <<endl;


//	for( int i=0; i<20 ; i++){
//		cout<<i<< "  " <<getDegreeUsingDI(data, i, false)<<endl;
//	}
//
//	for( int i=0; i<20 ; i++){
//		cout<<i<< "  " <<getDegreeUsingSemiDI(data, i, false)<<endl;
//	}
//	for( int i=0; i<20 ; i++){
//		cout<<i<< "  " <<getDegreeUsingrrrDI(data, i, false)<<endl;
//	}
//	for( int i=0; i<20 ; i++){
//		cout<<i<< "  " <<getDegreeUsingsdbDI(data, i, false)<<endl;
//	}


//    cout << indexStructure.blockNum <<endl;
//    cout << indexStructure.subBlockCounter <<endl;
//    int max = 0;
//    for ( int i=0; i<indexStructure.blockNum -1  ; i++ ){
//    	//cout << indexStructure.bitVofEachBlock[i] << "  "  << indexStructure.indrectIndex[i].size() <<endl;
//    	//cout << max << "   " << indexStructure.indrectIndex[i].size()<<endl;
//    	if (indexStructure.indrectIndex[i].size() > max){
//    		max = indexStructure.indrectIndex[i].size();
//    	}
//    }
//    cout <<"max"<< max <<endl;
    //15606
    cout << "Using Indirect : " <<endl;
    for ( int i =0 ; i< 15590 ;i++){
    	getAdjForOneVUsingIDI(data , i);
    }
    cout<< "complete!~" <<endl;

    vector<int> res1 = getAdjForOneVUsingIDI(data, 13000);
    cout << " vid '9000' :" <<mappingForRestore[13000]<<endl;
    for (int i=0; i<res1.size();i++){
    	cout << res1[i] << "  " ;
    }
    //testUnit(data);
    //testUnitWithD(data, 0);





//    cout << " vid '1801' :" <<mappingForRestore[1801]<<endl;
//    cout << getDegreeUsingDI(data, 1801, true)<<endl;
//    cout << getDegreeUsingIDI(data, 1801)<<endl;



	delete FirstRound;
	delete mapAfterRenumber;
	delete root;

	FirstRound = NULL;
	mapAfterRenumber = NULL;
	root = NULL;

    return 0;
}

string IntToString1 (int a){
    ostringstream temp;
    temp<<a;
    return temp.str();
}




void readInput(){
	string line, line2;
	//ifstream myfile1 ("/home/xiangzhang/MetisTestGraph/m14b/m14bXadj.txt");  // Reading data of m14b graph
	//ifstream myfile2 ("/home/xiangzhang/MetisTestGraph/m14b/m14bAdj.txt");
	//ifstream myfile1 ("/home/xiangzhang/MetisTestGraph/feocean/feoceanXadj.txt");  // Reading data of feocean
	//ifstream myfile2 ("/home/xiangzhang/MetisTestGraph/feocean/feoceanAdj.txt");
	ifstream myfile1 ("/home/xiangzhang/MetisTestGraph/4elt_xadj.txt");  // Reading data of 4elt
	ifstream myfile2 ("/home/xiangzhang/MetisTestGraph/4elt_adj.txt");
	//ifstream myfile1 ("/home/xiangzhang/MetisTestGraph/144/144_xadj.txt");  // Reading data of 144.graph
	//ifstream myfile2 ("/home/xiangzhang/MetisTestGraph/144/144_adj.txt");
	//ifstream myfile1 ("/home/xiangzhang/MetisTestGraph/598a/598a_xadj.txt");  // Reading data of 598a.graph
	//ifstream myfile2 ("/home/xiangzhang/MetisTestGraph/598a/598a_adj.txt");
	//ifstream myfile1 ("/home/xiangzhang/MetisTestGraph/bcsstk30/bcsstk30_xadj.txt");  // Reading data of bcsstk30.graph
	//ifstream myfile2 ("/home/xiangzhang/MetisTestGraph/bcsstk30/bcsstk30_adj.txt");
	//ifstream myfile1 ("/home/xiangzhang/MetisTestGraph/bcsstk31/bcsstk31_xadj.txt");  // Reading data of bcsstk31.graph
	//ifstream myfile2 ("/home/xiangzhang/MetisTestGraph/bcsstk31/bcsstk31_adj.txt");
	//ifstream myfile1 ("/home/xiangzhang/MetisTestGraph/fe_body/fe_body_xadj.txt");  // Reading data of febody.graph
	//ifstream myfile2 ("/home/xiangzhang/MetisTestGraph/fe_body/fe_body_adj.txt");
	//ifstream myfile1 ("/home/xiangzhang/MetisTestGraph/wave/wave_xadj.txt");  // Reading data of wave.graph
	//ifstream myfile2 ("/home/xiangzhang/MetisTestGraph/wave/wave_adj.txt");


    //ifstream myfile1 ("/raid6/workspace/zichu/xiang/macsproject/graph_set/m14b/m14bXadj.txt");  // Reading data of m14b graph
    //ifstream myfile2 ("/raid6/workspace/zichu/xiang/macsproject/graph_set/m14b/m14bAdj.txt");
    //ifstream myfile1 ("/raid6/workspace/zichu/xiang/macsproject/graph_set/feocean/feoceanXadj.txt");  // Reading data of feocean
    //ifstream myfile2 ("/raid6/workspace/zichu/xiang/macsproject/graph_set/feocean/feoceanAdj.txt");
    //ifstream myfile1 ("/raid6/workspace/zichu/xiang/macsproject/graph_set/4elt/4elt_xadj.txt");  // Reading data of 4elt
    //ifstream myfile2 ("/raid6/workspace/zichu/xiang/macsproject/graph_set/4elt/4elt_adj.txt");
    //ifstream myfile1 ("/raid6/workspace/zichu/xiang/macsproject/graph_set/144/144_xadj.txt");  // Reading data of 144.graph
    //ifstream myfile2 ("/raid6/workspace/zichu/xiang/macsproject/graph_set/144/144_adj.txt");
    //ifstream myfile1 ("/raid6/workspace/zichu/xiang/macsproject/graph_set/598a/598a_xadj.txt");  // Reading data of 598a.graph
    //ifstream myfile2 ("/raid6/workspace/zichu/xiang/macsproject/graph_set/598a/598a_adj.txt");
    //ifstream myfile1 ("/raid6/workspace/zichu/xiang/macsproject/graph_set/bcsstk30/bcsstk30_xadj.txt");  // Reading data of bcsstk30.graph
    //ifstream myfile2 ("/raid6/workspace/zichu/xiang/macsproject/graph_set/bcsstk30/bcsstk30_adj.txt");
    //ifstream myfile1 ("/raid6/workspace/zichu/xiang/macsproject/graph_set/bcsstk31/bcsstk31_xadj.txt");  // Reading data of bcsstk31.graph
    //ifstream myfile2 ("/raid6/workspace/zichu/xiang/macsproject/graph_set/bcsstk31/bcsstk31_adj.txt");
    //ifstream myfile1 ("/raid6/workspace/zichu/xiang/macsproject/graph_set/fe_body/fe_body_xadj.txt");  // Reading data of febody.graph
    //ifstream myfile2 ("/raid6/workspace/zichu/xiang/macsproject/graph_set/fe_body/fe_body_adj.txt");
	//ifstream myfile1 ("/raid6/workspace/zichu/xiang/macsproject/graph_set/wave/wave_xadj.txt");  // Reading data of bcsstk31.graph
	//ifstream myfile2 ("/raid6/workspace/zichu/xiang/macsproject/graph_set/wave/wave_adj.txt");



	if (myfile1.is_open()){
		int pos = 0;
		while ( getline (myfile1,line) )
		{
			InputXadj[pos] = stoi(line);
			pos ++ ;
		}
		myfile1.close();
		//cout << pos <<endl;

	}

	else cout << "Unable to open file";

	if (myfile2.is_open()){
		int pos2 = 0;
		while ( getline (myfile2,line2) )
		{
			Inputadj[pos2] = stoi(line2);
			pos2 ++ ;
		}
		myfile2.close();
		//cout << pos2<<endl;
	}
	else cout << "Unable to open file";


}

void testUnit( phsim::BitString data ){
	std::chrono::time_point<std::chrono::system_clock> start1, end1;
	std::chrono::duration<double> elapsed_seconds;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	start1 = std::chrono::system_clock::now();

	for( int i=0; i<100000 ; i++){
		getDegreeUsingDI(data, 1401, false);

	}

	end1 = std::chrono::system_clock::now();
	elapsed_seconds = end1-start1;
	std::cout << "DI elapsed time: " << elapsed_seconds.count() << "s\n";

	start1 = std::chrono::system_clock::now();

	for( int i=0; i<100000 ; i++){
		getDegreeUsingSemiDI(data , 1401, false);

	}

	end1 = std::chrono::system_clock::now();
	elapsed_seconds = end1-start1;
	std::cout << "SEMI DI elapsed time: " << elapsed_seconds.count() << "s\n";

	start1 = std::chrono::system_clock::now();

	for( int i=0; i<100000 ; i++){

		getDegreeUsingrrrDI(data, 1401, false);

	}

	end1 = std::chrono::system_clock::now();
	elapsed_seconds = end1-start1;
	std::cout << "RRR elapsed time: " << elapsed_seconds.count() << "s\n";

	start1 = std::chrono::system_clock::now();

	for( int i=0; i<100000 ; i++){
		getDegreeUsingsdbDI(data, 1401, false);
	}

	end1 = std::chrono::system_clock::now();
	elapsed_seconds = end1-start1;
	std::cout << "SDB elapsed time: " << elapsed_seconds.count() << "s\n";

	cout << "----------------------------------------------------" <<endl;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	start1 = std::chrono::system_clock::now();

	for ( int i=0; i< 100000 ; i++){
		getAdjForOneVUsingDI(data , 1401);

	}

	end1 = std::chrono::system_clock::now();
	elapsed_seconds = end1-start1;
	std::cout << "DI N-L elapsed time: " << elapsed_seconds.count() << "s\n";

	start1 = std::chrono::system_clock::now();

	for ( int i=0; i< 100000 ; i++){

		getAdjForOneVUsingSemiDI(data , 1401);

	}

	end1 = std::chrono::system_clock::now();
	elapsed_seconds = end1-start1;
	std::cout << "SEMI DI N-L elapsed time: " << elapsed_seconds.count() << "s\n";

	start1 = std::chrono::system_clock::now();

	for ( int i=0; i< 100000 ; i++){

		getAdjForOneVUsingrrrDI(data , 1401);

	}

	end1 = std::chrono::system_clock::now();
	elapsed_seconds = end1-start1;
	std::cout << "RRR N-L elapsed time: " << elapsed_seconds.count() << "s\n";

	start1 = std::chrono::system_clock::now();

	for ( int i=0; i< 100000 ; i++){

		getAdjForOneVUsingsdbDI(data , 1401);

	}

	end1 = std::chrono::system_clock::now();
	elapsed_seconds = end1-start1;
	std::cout << "SDB N-L elapsed time: " << elapsed_seconds.count() << "s\n";

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



}

void testUnitWithD (phsim::BitString data , int vid){

	std::chrono::time_point<std::chrono::system_clock> start1, end1;
	std::chrono::duration<double> elapsed_seconds;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	start1 = std::chrono::system_clock::now();

	for( int i=0; i<100000 ; i++){
		getDegreeUsingDI(data, vid, true);

	}

	end1 = std::chrono::system_clock::now();
	elapsed_seconds = end1-start1;
	std::cout << "DI elapsed time: " << elapsed_seconds.count() << "s\n";

	start1 = std::chrono::system_clock::now();

	for( int i=0; i<100000 ; i++){
		getDegreeUsingSemiDI(data , vid, true);

	}

	end1 = std::chrono::system_clock::now();
	elapsed_seconds = end1-start1;
	std::cout << "SEMI DI elapsed time: " << elapsed_seconds.count() << "s\n";

	start1 = std::chrono::system_clock::now();

	for( int i=0; i<100000 ; i++){

		getDegreeUsingrrrDI(data, vid, true);

	}

	end1 = std::chrono::system_clock::now();
	elapsed_seconds = end1-start1;
	std::cout << "RRR elapsed time: " << elapsed_seconds.count() << "s\n";

	start1 = std::chrono::system_clock::now();

	for( int i=0; i<100000 ; i++){
		getDegreeUsingsdbDI(data, vid, true);
	}

	end1 = std::chrono::system_clock::now();
	elapsed_seconds = end1-start1;
	std::cout << "SDB elapsed time: " << elapsed_seconds.count() << "s\n";


	start1 = std::chrono::system_clock::now();

	for( int i=0; i<100000 ; i++){
		getDegreeUsingIDI(data, vid);
	}

	end1 = std::chrono::system_clock::now();
	elapsed_seconds = end1-start1;
	std::cout << "IDI elapsed time: " << elapsed_seconds.count() << "s\n";

	cout << "----------------------------------------------------" <<endl;
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	start1 = std::chrono::system_clock::now();

	for ( int i=0; i< 100000 ; i++){
		getAdjForOneVUsingDIWithD(data , vid);

	}

	end1 = std::chrono::system_clock::now();
	elapsed_seconds = end1-start1;
	std::cout << "DI N-L elapsed time: " << elapsed_seconds.count() << "s\n";


	start1 = std::chrono::system_clock::now();

	for ( int i=0; i< 100000 ; i++){

		getAdjForOneVUsingSemiDIWithD(data , vid);

	}

	end1 = std::chrono::system_clock::now();
	elapsed_seconds = end1-start1;
	std::cout << "SEMI DI N-L elapsed time: " << elapsed_seconds.count() << "s\n";

	start1 = std::chrono::system_clock::now();

	for ( int i=0; i< 100000 ; i++){

		getAdjForOneVUsingrrrDIWithD(data , vid);

	}

	end1 = std::chrono::system_clock::now();
	elapsed_seconds = end1-start1;
	std::cout << "RRR N-L elapsed time: " << elapsed_seconds.count() << "s\n";

	start1 = std::chrono::system_clock::now();

	for ( int i=0; i< 100000 ; i++){

		getAdjForOneVUsingsdbDIWithD(data , vid);

	}

	end1 = std::chrono::system_clock::now();
	elapsed_seconds = end1-start1;
	std::cout << "SDB N-L elapsed time: " << elapsed_seconds.count() << "s\n";


	start1 = std::chrono::system_clock::now();

	for ( int i=0; i< 100000 ; i++){
		getAdjForOneVUsingIDI(data , vid);

	}

	end1 = std::chrono::system_clock::now();
	elapsed_seconds = end1-start1;
	std::cout << "IDI N-L elapsed time: " << elapsed_seconds.count() << "s\n";

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



}



BstNode * generateLeaves (idx_t nVertex , idx_t *xadjancy, idx_t *adjancy , int ** mapOrg){
	//BstNode * topNode = NULL;

	if( nVertex == 3){
		//case1: xadj = {0,1,3,4} adj= {1,0,2,1}
		if ( xadjancy[0]==0 && xadjancy[1]==1 && xadjancy[2]==3 && xadjancy[3]== 4 ){
			//string str1 = "0,1";
			//string str2 = "1,2";
			/*
			if(mapOrg[0][0]==23 ||mapOrg[1][0]==23 ||mapOrg[2][0]==23){
				cout << "case1"<<endl;
			}*/

			string str1 = IntToString1(mapOrg[0][0])+ "," + IntToString1(mapOrg[1][0]);
			string str2 = IntToString1(mapOrg[1][0])+ "," + IntToString1(mapOrg[2][0]);
			vector<string> v1;
			vector<string> v2;
			v1.push_back(str1);
			v2.push_back(str2);
			BstNode * topNode = getNewNode(v1,3);
			BstNode * leftChild = getNewNode(v2,2);
			BstNode * leafOne = getLeafNode(mapOrg[0][0]);
			BstNode * leafTwo = getLeafNode(mapOrg[1][0]);
			BstNode * leafThree = getLeafNode(mapOrg[2][0]);
			topNode->left = leftChild;
			topNode->right = leafOne;
			leftChild->left = leafTwo;
			leftChild->right = leafThree;

			return topNode;

		}else if(xadjancy[0]==0 && xadjancy[1]==1 && xadjancy[2]==2 && xadjancy[3]== 4){ //case2: xadj = {0,1,2,4} adj = {2,2,0,1}

			//string str1 = "0,2";
			//string str2 = "1,2";
			/*
			if(mapOrg[0][0]==23 ||mapOrg[1][0]==23 ||mapOrg[2][0]==23){
				cout << "case2"<<endl;
			}*/


			string str1 = IntToString1(mapOrg[0][0])+ "," + IntToString1(mapOrg[2][0]);
			string str2 = IntToString1(mapOrg[1][0])+ "," + IntToString1(mapOrg[2][0]);
			vector<string> v1;
			vector<string> v2;
			v1.push_back(str1);
			v2.push_back(str2);
			BstNode * topNode = getNewNode(v1,3);
			BstNode * leftChild = getNewNode(v2,2);
			BstNode * leafOne = getLeafNode(mapOrg[0][0]);
			BstNode * leafTwo = getLeafNode(mapOrg[1][0]);
			BstNode * leafThree = getLeafNode(mapOrg[2][0]);
			topNode->left = leftChild;
			topNode->right = leafOne;
			leftChild->left = leafTwo;
			leftChild->right = leafThree;

			return topNode;

		}else if(xadjancy[0]==0 && xadjancy[1]==2 && xadjancy[2]==3 && xadjancy[3]== 4){ //case3: xadj = {0,2,3,4} adj = {1,2,0,0}
			//string str1 = "0,1";
			//string str2 = "0,2";
			/*
			if(mapOrg[0][0]==23 ||mapOrg[1][0]==23 ||mapOrg[2][0]==23){
				cout << "case3"<<endl;
			}*/
			string str1 = IntToString1(mapOrg[0][0])+ "," + IntToString1(mapOrg[1][0]);
			string str2 = IntToString1(mapOrg[0][0])+ "," + IntToString1(mapOrg[2][0]);
			vector<string> v1;
			vector<string> v2;
			v1.push_back(str1);
			v2.push_back(str2);
			BstNode * topNode = getNewNode(v1,3);
			BstNode * leftChild = getNewNode(v2,2);
			BstNode * leafOne = getLeafNode(mapOrg[1][0]);
			BstNode * leafTwo = getLeafNode(mapOrg[0][0]);
			BstNode * leafThree = getLeafNode(mapOrg[2][0]);
			topNode->right = leafOne;
			topNode->left = leftChild;
			leftChild->left = leafTwo;
			leftChild->right = leafThree;

			return topNode;
		}else if(xadjancy[0]==0 && xadjancy[1]==2 && xadjancy[2]==4 && xadjancy[3]== 6){ //case4: xadj = {0,2,4,6} adj = {1,2,0,2,0,1}
			//string str1 = "0,1";
			//string str2 = "0,2";
			//string str3 = "1,2";
			/*
			if(mapOrg[0][0]==23 ||mapOrg[1][0]==23 ||mapOrg[2][0]==23){
				cout << "case4"<<endl;
			}*/
			string str1 = IntToString1(mapOrg[0][0])+ "," + IntToString1(mapOrg[1][0]);
			string str2 = IntToString1(mapOrg[0][0])+ "," + IntToString1(mapOrg[2][0]);
			string str3 = IntToString1(mapOrg[1][0])+ "," + IntToString1(mapOrg[2][0]);
			vector<string> v1;
			vector<string> v2;
			v1.push_back(str1);
			v1.push_back(str2);
			v2.push_back(str3);
			BstNode * topNode = getNewNode(v1,3);
			BstNode * leftChild = getNewNode(v2,2);
			BstNode * leafOne = getLeafNode(mapOrg[0][0]);
			BstNode * leafTwo = getLeafNode(mapOrg[1][0]);
			BstNode * leafThree = getLeafNode(mapOrg[2][0]);
			topNode->right = leafOne;
			topNode->left = leftChild;
			leftChild->left = leafTwo;
			leftChild->right = leafThree;

			return topNode;
		}else{

			if( xadjancy[0]==0 && xadjancy[1]==1 && xadjancy[2]==1 && xadjancy[3]== 2 && adjancy[0]==2 && adjancy[1] == 0){ //case5: xadj = {0,1,1,2} adj = {2,0}
				/*
				if(mapOrg[0][0]==23 ||mapOrg[1][0]==23 ||mapOrg[2][0]==23){
					cout << "case5"<<endl;
				}*/
				vector<string> v1;
				vector<string> v2;
				string str0 = "NoE";
				v1.push_back(str0) ; // Means that no edge separator here
				string str1 = IntToString1(mapOrg[0][0])+ "," + IntToString1(mapOrg[2][0]);
				v2.push_back(str1);
				BstNode * topNode = getNewNode(v1,3);
				BstNode * leftChild = getNewNode(v2,2);
				BstNode * leafOne = getLeafNode(mapOrg[1][0]);
				BstNode * leafTwo = getLeafNode(mapOrg[0][0]);
				BstNode * leafThree = getLeafNode(mapOrg[2][0]);
				topNode->right = leafOne;
				topNode->left = leftChild;
				leftChild->left = leafTwo;
				leftChild->right = leafThree;

				//return NULL;
				return topNode;

			}else if(xadjancy[0]==0 && xadjancy[1]==0 && xadjancy[2]==0 && xadjancy[3]== 0){ //case6: xadj = {0,0,0,0}
				/*
				if(mapOrg[0][0]==23 ||mapOrg[1][0]==23 ||mapOrg[2][0]==23){
					cout << "case6"<<endl;
				}*/
				string str0 = "NoE";
				vector<string> v1;
				vector<string> v2;
				v1.push_back(str0) ; // Means that no edge separator here
				v2.push_back(str0);
				BstNode * topNode = getNewNode(v1,3);
				BstNode * leftChild = getNewNode(v2,2);
				BstNode * leafOne = getLeafNode(mapOrg[0][0]);
				BstNode * leafTwo = getLeafNode(mapOrg[1][0]);
				BstNode * leafThree = getLeafNode(mapOrg[2][0]);
				topNode->right = leafOne;
				topNode->left = leftChild;
				leftChild->left = leafTwo;
				leftChild->right = leafThree;

				//return NULL;
				return topNode;

			}else if(xadjancy[0]==0 && xadjancy[1]==0 && xadjancy[2]==1 && xadjancy[3]== 2 && adjancy[0]==2 && adjancy[1] == 1){ //case7: xadj = {0,0,1,2} adj = {2,1}
				/*
				if(mapOrg[0][0]==23 ||mapOrg[1][0]==23 ||mapOrg[2][0]==23){
					cout << "case7"<<endl;
				}*/
				string str0 = "NoE";
				vector<string> v1;
				vector<string> v2;
				v1.push_back(str0) ; // Means that no edge separator here
				string str1 = IntToString1(mapOrg[1][0])+ "," + IntToString1(mapOrg[2][0]);
				v2.push_back(str1);
				BstNode * topNode = getNewNode(v1,3);
				BstNode * leftChild = getNewNode(v2,2);
				BstNode * leafOne = getLeafNode(mapOrg[0][0]);
				BstNode * leafTwo = getLeafNode(mapOrg[1][0]);
				BstNode * leafThree = getLeafNode(mapOrg[2][0]);
				topNode->right = leafOne;
				topNode->left = leftChild;
				leftChild->left = leafTwo;
				leftChild->right = leafThree;

				//return NULL;
				return topNode;

			}else if(xadjancy[0]==0 && xadjancy[1]==1 && xadjancy[2]==2 && xadjancy[3]== 2 && adjancy[0]==1 && adjancy[1] == 0){ //case8: xadj = {0,1,2,2} adj = {1,0}
				/*
				if(mapOrg[0][0]==23 ||mapOrg[1][0]==23 ||mapOrg[2][0]==23){
					cout << "case8"<<endl;
				}*/

				vector<string> v1;
				vector<string> v2;
				string str0 = "NoE";
				v1.push_back(str0) ; // Means that no edge separator here
				string str1 = IntToString1(mapOrg[0][0])+ "," + IntToString1(mapOrg[1][0]);
				v2.push_back(str1);
				BstNode * topNode = getNewNode(v1,3);
				BstNode * leftChild = getNewNode(v2,2);
				BstNode * leafOne = getLeafNode(mapOrg[2][0]);
				BstNode * leafTwo = getLeafNode(mapOrg[0][0]);
				BstNode * leafThree = getLeafNode(mapOrg[1][0]);
				topNode->right = leafOne;
				topNode->left = leftChild;
				leftChild->left = leafTwo;
				leftChild->right = leafThree;

				//return NULL;

				return topNode;

			}else{//caseSpecial: Unknown case
				cout<< "Unknown case!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"<<endl;
				cout<< xadjancy[0]<< " " <<xadjancy[1]<< " " << xadjancy[2] << " "  << xadjancy [3] <<endl;
				cout<< adjancy[0]<< " " << adjancy[1] << endl<<endl;
				return NULL;
			}
		}

	}else if (nVertex == 2){
		//cout << "nVertex == 2:  "<<endl;
		/*
		for (int i = 0 ; i < 2 ; i++){
			cout << mapOrg[i][0] << " : " <<mapOrg[i][1]<<endl;
		}*/
		//cout<< adjancy[0]<< " " << adjancy[1] <<endl;
		//string str1 = "0,1";
		string str1 = IntToString1(mapOrg[0][0])+ "," + IntToString1(mapOrg[1][0]);
		vector<string> v1;
		v1.push_back(str1);
		BstNode *leaf = getNewNode(v1,2);
		BstNode * leafOne = getLeafNode(mapOrg[0][0]);
		BstNode * leafTwo = getLeafNode(mapOrg[1][0]);
		leaf->left = leafOne;
		leaf->right = leafTwo;
		return leaf;

	}else{
		cout << "HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH"<<endl;
		return NULL;
	}
}

BstNode * generateTree ( idx_t nVertex , idx_t *xadj, idx_t *adj , int **mapOrg){

	//termination
	if(nVertex <= 3){

		BstNode * leafNode = NULL;
		leafNode = generateLeaves (nVertex , xadj, adj, mapOrg );
		return leafNode;
		//return NULL;
	}

    BstNode * curr;
	//define the root
	if(nVertex == nV){
		MetisNode mynode( nVertex ,xadj,adj);
		//mynode.printMapping();
		//mynode.printAdj();
		root = getNewNode (mynode.EdgeSeparator , nV);
		curr =  root;
		curr->left = generateTree (mynode.NumberOfVL , mynode.xal , mynode.adjAftMLC , mynode.mapOrgLC);
		curr->right = generateTree (mynode.NumberOfVR , mynode.xar , mynode.adjAftMRC , mynode.mapOrgRC);
	}
	else{
		MetisNode mynode( nVertex ,xadj,adj, mapOrg);

		curr = getNewNode (mynode.EdgeSeparator , nVertex);

		curr->left = generateTree (mynode.NumberOfVL , mynode.xal , mynode.adjAftMLC , mynode.mapOrgLC);
		curr->right = generateTree (mynode.NumberOfVR , mynode.xar , mynode.adjAftMRC , mynode.mapOrgRC);
	}

	return curr;

}



void preorderTraversal(BstNode* x){
	if(x != NULL ){
		if(x->flag == 1){


			cout << "Edge Separators: "<<endl;
			for( int i = 0; i< x->EdgeSparators.size();i++){
				cout << x->EdgeSparators[i] << endl;
			}

			preorderTraversal(x->left);

			preorderTraversal(x->right);
		}else{
			cout << "Node contains only one Vertex: " <<endl;
			cout<<"Vertex id: " << x->VertexId <<endl;
			cout<<"New Vertex id: " << x->NewVertexId <<endl;
		}

	}
}




void inorderTraversal(BstNode* x){
	if(x != NULL){
		if(x->flag == 1){
			inorderTraversal(x->left);
			for( int i = 0; i< x->EdgeSparators.size();i++){
				cout <<  x->EdgeSparators[i] << endl;
			}
			inorderTraversal(x->right);

		}else{
			cout << "Node contains only one Vertex: " <<endl;
			cout<<"Vertex id: " << x->VertexId <<endl;
			cout<<"New Vertex id: " << x->NewVertexId <<endl;
		}
	}
}

void renumberByInorderTraversal(BstNode* x){

	if(x != NULL){
		renumberByInorderTraversal(x->left);

		if(x->flag == 0){

			x->NewVertexId = InorderTraversalVID;

			mapAfterRenumber[x->VertexId] = InorderTraversalVID;
			mappingForRestore[InorderTraversalVID] = x->VertexId;
			InorderTraversalVID++;

			//cout << x->VertexId <<endl;

		}

		renumberByInorderTraversal(x->right);
	}

}

void findNeighbour (BstNode *x){
	if(x != NULL){
		if(x->flag == 1){
			if(x->EdgeSparators.size() != 0 ){
				if( x->EdgeSparators[0] != "NoE"){
					for( int i=0; i<x->EdgeSparators.size(); i++){
						int vertexA , vertexB = -1;
						int index = x->EdgeSparators[i].find(",");
						vertexA = stoi(x->EdgeSparators[i].substr(0,index));
						vertexB = stoi(x->EdgeSparators[i].substr(index+1,x->EdgeSparators[i].size()));
						int mapA = mapAfterRenumber[vertexA];
						int mapB = mapAfterRenumber[vertexB];
						nlns.addItem(mapB , mapA );
						nlns.addItem(mapA , mapB );
					}
				}
			}
			findNeighbour(x->left);
			findNeighbour(x->right);
		}
	}
}

void formAdjTable(NeighbourListNonSorted NLS ){
    for(int i=0; i<nV; i++){
    	negativeDiff[i] = false;
    }

	for(int i = 0; i<NLS.tableSize ;i++){

		vector<phsim::BitString::Word> curVertexAdjList;
    	int v0 = i;
    	int v1 = NLS.sortedNeighbourTable[i][0];

    	string gammaCode = "";

		if( v0 > v1){ // Which means v0-v1 is less than 0
			//gammaCode = EncodingEliasGamma( v0-v1);
			negativeDiff[i] = true;
			curVertexAdjList.push_back((int)(v0-v1));

		}else{
			//gammaCode = EncodingEliasGamma( v1-v0);

			curVertexAdjList.push_back((int)(v1-v0));

		}
		//AdjList.push_back(gammaCode);

		for ( int j=0 ; j<NLS.numberOfItemsInIndex(i) -1 ; j++){
			v0 = v1;
			v1 = NLS.sortedNeighbourTable[i][j+1];

			string gammaCode1 = "";
			curVertexAdjList.push_back((int)(v1-v0));

    	}
		offsetIndex.push_back(NLS.numberOfItemsInIndex(i));
		inputWV.push_back(curVertexAdjList);
		//adjTable.push_back(AdjList);
		//cout << endl;
    }
    //cout << input.size();

}


void postOrderTraCountLeaf (BstNode* x){
	if(x != NULL ){
		if( x->flag == 1){
			if ( x->left->flag == 1 && x->right->flag == 1){
				postOrderTraCountLeaf( x->left);
				postOrderTraCountLeaf( x->right);

				x->array_len = x->left->array_len + x->right->array_len;
				x->ordered_array = new int [x->array_len];

				x->ordered_array = mergeArray (x->left->ordered_array , x->right->ordered_array,
													x->left->array_len , x->right->array_len);
			}

			if (x->left->flag == 0 && x->right->flag == 1){
				postOrderTraCountLeaf( x->right);
				x->array_len = x->left->array_len + x->right->array_len;
				x->ordered_array = new int [x->array_len];

				x->ordered_array = mergeArray (x->left->ordered_array , x->right->ordered_array,
												x->left->array_len , x->right->array_len);
			}

			if (x->left->flag == 1 and x->right->flag == 0 ){
				postOrderTraCountLeaf( x->left);
				x->array_len = x->left->array_len + x->right->array_len;
				x->ordered_array = new int [x->array_len];

				x->ordered_array = mergeArray (x->left->ordered_array , x->right->ordered_array,
									x->left->array_len , x->right->array_len);

			}

			if ( x->left->flag == 0 and x->right->flag == 0){
				x->array_len = x->left->array_len + x->right->array_len;
				x->ordered_array = new int [x->array_len];

				x->ordered_array = mergeArray (x->left->ordered_array , x->right->ordered_array,
									x->left->array_len , x->right->array_len);
			}
		}
	}
}



void preorderTraversalCheckLeaf(BstNode* x){
	if(x != NULL ){
		if(x->flag == 1){
			for (int i=0; i<x->array_len ;i++){
				cout << x->ordered_array [i] << " ";
			}
			cout << endl;

			preorderTraversalCheckLeaf(x->left);

			preorderTraversalCheckLeaf(x->right);
		}

	}
}

int* mergeArray( int* a , int* b, int len_a, int len_b ){
	int new_len = len_a + len_b;
	int* new_arr = new int [new_len];

	int i=0, j=0, index = 0;
	while ( i<len_a && j < len_b ){
		if (a[i] < b[j]){
			new_arr[index] = a[i];
			index++;
			i++;
		}else{
			new_arr[index] = b[j];
			index++;
			j++;
		}
	}
	while ( i < len_a){
		new_arr[index] = a[i];
		index++;
		i++;
	}
	while ( j < len_b){
		new_arr[index] = b[j];
		index++;
		j++;
	}
	return new_arr;
}


void childFlipping( BstNode * x , vector< int * >&n1, vector< int * >&n2 , vector<int> &n1_size , vector<int> &n2_size){
	if ( x->left != NULL && x->right != NULL){

		int count1L = 0;
		int count1R = 0;
		int count2R = 0;
		int count2L = 0;

		// if both children are not leaf node
		if ( x->left->flag != 0 && x->right->flag != 0 ){
			//cout << n1.size() << "   "<<n2.size() << endl;

			if( !n1.size() == 0 || !n2.size()==0 ){  // when current node is not the root

				if ( n1.size()==0 ){

					for(int i=0; i < x->left->array_len; i++){

						int index1 = x->left->ordered_array[i];
						int index2 = index1+1;
						int neighbor_num = InputXadj[index2]-InputXadj[index1];

						for( int j=InputXadj[index1]; j<InputXadj[index1]+neighbor_num; j++){
							int anotherV = Inputadj[j];

							for ( int k=0; k<n2.size();k++){
//								if (n2[k].count(anotherV)){
//									count2L ++;
//								}
							  int * pItem;
							  pItem = (int*) bsearch (&anotherV, n2[k], n2_size[k], sizeof (int), compareints);
							  if (pItem!=NULL)
								  count2L ++;
							}
						}

					}

					for(int i=0; i < x->right->array_len; i++){

						int index1 = x->right->ordered_array[i];
						int index2 = index1+1;
						int neighbor_num = InputXadj[index2]-InputXadj[index1];

						for( int j=InputXadj[index1]; j<InputXadj[index1]+neighbor_num; j++){
							int anotherV = Inputadj[j];

							for ( int k=0; k<n2.size();k++){
//								if (n2[k].count(anotherV)){
//									count2R ++;
//								}
								int * pItem;
								pItem = (int*) bsearch (&anotherV, n2[k], n2_size[k], sizeof (int), compareints);
								if (pItem!=NULL)
									count2R ++;
							}
						}
					}

				}else if( n2.size()==0){

					for(int i=0; i < x->left->array_len; i++){

						int index1 = x->left->ordered_array[i];
						int index2 = index1+1;
						int neighbor_num = InputXadj[index2]-InputXadj[index1];

						for( int j=InputXadj[index1]; j<InputXadj[index1]+neighbor_num; j++){
							int anotherV = Inputadj[j];

							for ( int k=0; k<n1.size();k++){
//								if (n1[k].count(anotherV)){
//									count1L ++;
//								}
								int * pItem;
								pItem = (int*) bsearch (&anotherV, n1[k], n1_size[k], sizeof (int), compareints);
								if (pItem!=NULL)
									count1L ++;
							}
						}
					}

					for(int i=0; i < x->right->array_len; i++){

						int index1 = x->right->ordered_array[i];
						int index2 = index1+1;
						int neighbor_num = InputXadj[index2]-InputXadj[index1];

						for( int j=InputXadj[index1]; j<InputXadj[index1]+neighbor_num; j++){
							int anotherV = Inputadj[j];

							for ( int k=0; k<n1.size();k++){
//								if (n1[k].count(anotherV)){
//									count1R ++;
//								}
								int * pItem;
								pItem = (int*) bsearch (&anotherV, n1[k], n1_size[k], sizeof (int), compareints);
								if (pItem!=NULL)
									count1R ++;
							}
						}
					}

				}else{


					for(int i=0; i< x->left->array_len; i++){

						int index1 = x->left->ordered_array[i];
						int index2 = index1+1;
						int neighbor_num = InputXadj[index2]-InputXadj[index1];

						for( int j=InputXadj[index1]; j<InputXadj[index1]+neighbor_num; j++){
							int anotherV = Inputadj[j];

							for ( int k=0; k<n1.size();k++){
//								if (n1[k].count(anotherV)){
//									count1L ++;
//								}
								int * pItem;
								pItem = (int*) bsearch (&anotherV, n1[k], n1_size[k], sizeof (int), compareints);
								if (pItem!=NULL)
									count1L ++;
							}

							for ( int m=0; m<n2.size();m++){
//								if (n2[m].count(anotherV)){
//									count2L ++;
//								}
								int * pItem;
								pItem = (int*) bsearch (&anotherV, n2[m], n2_size[m], sizeof (int), compareints);
								if (pItem!=NULL)
									count2L ++;
							}
						}
					}

					for(int i=0; i< x->right->array_len; i++){

						int index1 = x->right->ordered_array[i];
						int index2 = index1+1;
						int neighbor_num = InputXadj[index2]-InputXadj[index1];

						for( int j=InputXadj[index1]; j<InputXadj[index1]+neighbor_num; j++){
							int anotherV = Inputadj[j];
							for ( int k=0; k<n1.size();k++){
//								if (n1[k].count(anotherV)){
//									count1R ++;
//								}
								int * pItem;
								pItem = (int*) bsearch (&anotherV, n1[k], n1_size[k], sizeof (int), compareints);
								if (pItem!=NULL)
									count1R ++;
							}
							for ( int m=0; m<n2.size();m++){
//								if (n2[m].count(anotherV)){
//									count2R ++;
//								}
								int * pItem;
								pItem = (int*) bsearch (&anotherV, n2[m], n2_size[m], sizeof (int), compareints);
								if (pItem!=NULL)
									count2R ++;
							}
						}
					}
				}
				// Flip the child nodes if En1nl + En2nr < Enln2 + En1nr

				//cout << countL1 << "  " << count2R << "  " <<countL2 << "  "<< count1R <<endl;

                if ((count1L + count2R) < (count2L + count1R)){
                        BstNode * temp = NULL;
                        temp = x->left;
                        x->left = x->right;
                        x->right = temp;
                }
			}

			n2.push_back(x->right->ordered_array);
			n2_size.push_back(x->right->array_len);
			childFlipping( x->left , n1, n2, n1_size , n2_size);
			n2.pop_back();
			n2_size.pop_back();

			n1.push_back(x->left->ordered_array);
			n1_size.push_back(x->left->array_len);
			childFlipping( x->right, n1, n2, n1_size, n2_size);
			n1.pop_back();
			n1_size.pop_back();

		}

		// left child is a leaf while not right child
		if ( x->left->flag == 0 && x->right->flag == 1 ){

			for(int i=0; i< x->left->array_len; i++){

				int index1 = x->left->ordered_array[i];
				int index2 = index1+1;
				int neighbor_num = InputXadj[index2]-InputXadj[index1];

				for( int j=InputXadj[index1]; j<InputXadj[index1]+neighbor_num; j++){
					int anotherV = Inputadj[j];

					for ( int k=0; k<n1.size();k++){
						//if (n1[k].count(anotherV)){
						//	count1L ++;
						//}
						int * pItem;
						pItem = (int*) bsearch (&anotherV, n1[k], n1_size[k], sizeof (int), compareints);
						if (pItem!=NULL)
							count1L ++;
					}

					for ( int m=0; m<n2.size();m++){
						//if (n2[m].count(anotherV)){
						//	count2L ++;
						//}
						int * pItem;
						pItem = (int*) bsearch (&anotherV, n2[m], n2_size[m], sizeof (int), compareints);
						if (pItem!=NULL)
							count2L ++;
					}
				}
			}

			for(int i=0; i< x->right->array_len; i++){

				int index1 = x->right->ordered_array[i];
				int index2 = index1+1;
				int neighbor_num = InputXadj[index2]-InputXadj[index1];

				for( int j=InputXadj[index1]; j<InputXadj[index1]+neighbor_num; j++){
					int anotherV = Inputadj[j];
					for ( int k=0; k<n1.size();k++){
						//if (n1[k].count(anotherV)){
						//	count1R ++;
						//}
						int * pItem;
						pItem = (int*) bsearch (&anotherV, n1[k], n1_size[k], sizeof (int), compareints);
						if (pItem!=NULL)
							count1R ++;
					}
					for ( int m=0; m<n2.size();m++){
						//if (n2[m].count(anotherV)){
						//	count2R ++;
						//}
						int * pItem;
						pItem = (int*) bsearch (&anotherV, n2[m], n2_size[m], sizeof (int), compareints);
						if (pItem!=NULL)
							count2R ++;
					}
				}
			}

			if ((count1L + count2R) < (count2L + count1R)){
                    BstNode * temp = NULL;
                    temp = x->left;
                    x->left = x->right;
                    x->right = temp;

        			n2.push_back(x->right->ordered_array);
        			n2_size.push_back(x->right->array_len);
        			childFlipping( x->left , n1, n2, n1_size , n2_size);  // now right child is the leaf
        			n2.pop_back();
        			n2_size.pop_back();

			}else{

				n1.push_back(x->left->ordered_array);
				n1_size.push_back(x->left->array_len);
				childFlipping( x->right, n1, n2, n1_size, n2_size);
				n1.pop_back();
				n1_size.pop_back();
            }
		}
		// if the right child is a leaf
		if ( x->left->flag == 1 && x->right->flag == 0 ){

			for(int i=0; i< x->left->array_len; i++){

				int index1 = x->left->ordered_array[i];
				int index2 = index1+1;
				int neighbor_num = InputXadj[index2]-InputXadj[index1];

				for( int j=InputXadj[index1]; j<InputXadj[index1]+neighbor_num; j++){
					int anotherV = Inputadj[j];

					for ( int k=0; k<n1.size();k++){
						//if (n1[k].count(anotherV)){
						//	count1L ++;
						//}
						int * pItem;
						pItem = (int*) bsearch (&anotherV, n1[k], n1_size[k], sizeof (int), compareints);
						if (pItem!=NULL)
							count1L ++;
					}

					for ( int m=0; m<n2.size();m++){
						//if (n2[m].count(anotherV)){
						//	count2L ++;
						//}
						int * pItem;
						pItem = (int*) bsearch (&anotherV, n2[m], n2_size[m], sizeof (int), compareints);
						if (pItem!=NULL)
							count2L ++;
					}
				}
			}

			for(int i=0; i< x->right->array_len; i++){

				int index1 = x->right->ordered_array[i];
				int index2 = index1+1;
				int neighbor_num = InputXadj[index2]-InputXadj[index1];

				for( int j=InputXadj[index1]; j<InputXadj[index1]+neighbor_num; j++){
					int anotherV = Inputadj[j];
					for ( int k=0; k<n1.size();k++){
						//if (n1[k].count(anotherV)){
						//	count1R ++;
						//}
						int * pItem;
						pItem = (int*) bsearch (&anotherV, n1[k], n1_size[k], sizeof (int), compareints);
						if (pItem!=NULL)
							count1R ++;
					}
					for ( int m=0; m<n2.size();m++){
						//if (n2[m].count(anotherV)){
						//	count2R ++;
						//}
						int * pItem;
						pItem = (int*) bsearch (&anotherV, n2[m], n2_size[m], sizeof (int), compareints);
						if (pItem!=NULL)
							count2R ++;
					}
				}
			}

			if ((count1L + count2R) < (count2L + count1R)){
				BstNode * temp = NULL;
				temp = x->left;
				x->left = x->right;
				x->right = temp;


				n1.push_back(x->left->ordered_array);
				n1_size.push_back(x->left->array_len);
				childFlipping( x->right, n1, n2, n1_size, n2_size);   // now left child is the leaf
				n1.pop_back();
				n1_size.pop_back();


			}else{

    			n2.push_back(x->right->ordered_array);
    			n2_size.push_back(x->right->array_len);
    			childFlipping( x->left , n1, n2, n1_size , n2_size);
    			n2.pop_back();
    			n2_size.pop_back();
            }

		}

		// if both children are leaf node
		if ( x->left->flag == 0 && x->right->flag == 0 ){

			for(int i=0; i< x->left->array_len; i++){

				int index1 = x->left->ordered_array[i];
				int index2 = index1+1;
				int neighbor_num = InputXadj[index2]-InputXadj[index1];

				for( int j=InputXadj[index1]; j<InputXadj[index1]+neighbor_num; j++){
					int anotherV = Inputadj[j];

					for ( int k=0; k<n1.size();k++){
						//if (n1[k].count(anotherV)){
						//	count1L ++;
						//}
						int * pItem;
						pItem = (int*) bsearch (&anotherV, n1[k], n1_size[k], sizeof (int), compareints);
						if (pItem!=NULL)
							count1L ++;
					}

					for ( int m=0; m<n2.size();m++){
						//if (n2[m].count(anotherV)){
						//	count2L ++;
						//}
						int * pItem;
						pItem = (int*) bsearch (&anotherV, n2[m], n2_size[m], sizeof (int), compareints);
						if (pItem!=NULL)
							count2L ++;
					}
				}
			}

			for(int i=0; i< x->right->array_len; i++){

				int index1 = x->right->ordered_array[i];
				int index2 = index1+1;
				int neighbor_num = InputXadj[index2]-InputXadj[index1];

				for( int j=InputXadj[index1]; j<InputXadj[index1]+neighbor_num; j++){
					int anotherV = Inputadj[j];
					for ( int k=0; k<n1.size();k++){
						//if (n1[k].count(anotherV)){
						//	count1R ++;
						//}
						int * pItem;
						pItem = (int*) bsearch (&anotherV, n1[k], n1_size[k], sizeof (int), compareints);
						if (pItem!=NULL)
							count1R ++;
					}
					for ( int m=0; m<n2.size();m++){
						//if (n2[m].count(anotherV)){
						//	count2R ++;
						//}
						int * pItem;
						pItem = (int*) bsearch (&anotherV, n2[m], n2_size[m], sizeof (int), compareints);
						if (pItem!=NULL)
							count2R ++;
					}
				}
			}

			if ((count1L + count2R) < (count2L + count1R)){
				BstNode * temp = NULL;
				temp = x->left;
				x->left = x->right;
				x->right = temp;

			}

		}

	}
}



//
//void childFlipping1( BstNode * x , vector< set<int> >&n1, vector< set<int> >&n2){
//	if ( x->left != NULL && x->right != NULL){
//		if ( x->left->flag != 0 && x->right->flag != 0 ){
//			//cout << n1.size() << "   "<<n2.size() << endl;
//
//			int count1L = 0;
//			int count1R = 0;
//			int count2R = 0;
//			int count2L = 0;
//
//			if( !n1.size() == 0 || !n2.size()==0 ){  // when current node is not the root
//
//
//
//				if ( n1.size()==0 ){
//
//					set<int>::iterator cur2L = x->left->leaves.begin();
//					for(int i=0; i<x->left->leaves.size(); i++){
//
//						int index1 = *cur2L;
//						int index2 = index1+1;
//						int neighbor_num = InputXadj[index2]-InputXadj[index1];
//
//						for( int j=InputXadj[index1]; j<InputXadj[index1]+neighbor_num; j++){
//							int anotherV = Inputadj[j];
//
//							for ( int k=0; k<n2.size();k++){
//								if (n2[k].count(anotherV)){
//									count2L ++;
//								}
//							}
//						}
//						advance(cur2L,1);
//					}
//
//					set<int>::iterator cur2R = x->right->leaves.begin();
//					for(int i=0; i<x->right->leaves.size(); i++){
//
//						int index1 = *cur2R;
//						int index2 = index1+1;
//						int neighbor_num = InputXadj[index2]-InputXadj[index1];
//
//
//						for( int j=InputXadj[index1]; j<InputXadj[index1]+neighbor_num; j++){
//							int anotherV = Inputadj[j];
//							for ( int k=0; k<n2.size();k++){
//								if (n2[k].count(anotherV)){
//									count2R ++;
//								}
//							}
//						}
//						advance(cur2R,1);
//					}
//
//
//				}else if( n2.size()==0){
//
//					set<int>::iterator cur1L = x->left->leaves.begin();
//					for(int i=0; i<x->left->leaves.size(); i++){
//
//						int index1 = *cur1L;
//						int index2 = index1+1;
//						int neighbor_num = InputXadj[index2]-InputXadj[index1];
//
//						for( int j=InputXadj[index1]; j<InputXadj[index1]+neighbor_num; j++){
//							int anotherV = Inputadj[j];
//
//							for ( int k=0; k<n1.size();k++){
//								if (n1[k].count(anotherV)){
//									count1L ++;
//								}
//							}
//
//						}
//						advance(cur1L,1);
//					}
//
//					set<int>::iterator cur1R = x->right->leaves.begin();
//					for(int i=0; i<x->right->leaves.size(); i++){
//
//						int index1 = *cur1R;
//						int index2 = index1+1;
//						int neighbor_num = InputXadj[index2]-InputXadj[index1];
//
//						for( int j=InputXadj[index1]; j<InputXadj[index1]+neighbor_num; j++){
//							int anotherV = Inputadj[j];
//							for ( int k=0; k<n1.size();k++){
//								if (n1[k].count(anotherV)){
//									count1R ++;
//								}
//							}
//						}
//						advance(cur1R,1);
//					}
//
//				}else{
//
//					set<int>::iterator curL12 = x->left->leaves.begin();
//					for(int i=0; i< x->left->leaves.size(); i++){
//
//						int index1 = *curL12;
//						int index2 = index1+1;
//						int neighbor_num = InputXadj[index2]-InputXadj[index1];
//
//						for( int j=InputXadj[index1]; j<InputXadj[index1]+neighbor_num; j++){
//							int anotherV = Inputadj[j];
//							for ( int k=0; k<n1.size();k++){
//								if (n1[k].count(anotherV)){
//									count1L ++;
//								}
//							}
//							for ( int m=0; m<n2.size();m++){
//								if (n2[m].count(anotherV)){
//									count2L ++;
//								}
//							}
//
//						}
//						advance(curL12,1);
//					}
//
//					set<int>::iterator curR12 = x->right->leaves.begin();
//					for(int i=0; i< x->right->leaves.size(); i++){
//
//						int index1 = *curR12;
//						int index2 = index1+1;
//						int neighbor_num = InputXadj[index2]-InputXadj[index1];
//
//						for( int j=InputXadj[index1]; j<InputXadj[index1]+neighbor_num; j++){
//							int anotherV = Inputadj[j];
//							for ( int k=0; k<n1.size();k++){
//								if (n1[k].count(anotherV)){
//									count1R ++;
//								}
//							}
//							for ( int m=0; m<n2.size();m++){
//								if (n2[m].count(anotherV)){
//									count2R ++;
//								}
//							}
//						}
//						advance(curR12,1);
//					}
//				}
//				// Flip the child nodes if En1nl + En2nr < Enln2 + En1nr
//
//				//cout << countL1 << "  " << count2R << "  " <<countL2 << "  "<< count1R <<endl;
//
//                if ((count1L + count2R) < (count2L + count1R)){
//                        BstNode * temp = NULL;
//                        temp = x->left;
//                        x->left = x->right;
//                        x->right = temp;
//                }
//			}
//
//			n2.push_back(x->right->leaves);
//			childFlipping1( x->left , n1, n2);
//			n2.pop_back();
//
//			n1.push_back(x->left->leaves);
//			childFlipping1( x->right, n1, n2);
//			n1.pop_back();
//
//		}
//	}
//}


void directDecodeAll(phsim::BitString data){
	phsim::BitStringReader reader(data);
	while(!reader.eof()){
		directDecodeResult.push_back(phsim::gammaDecode(reader));
	}
}



int getDegreeUsingDI(phsim::BitString data, unsigned int vid, bool flag_degree){

	int degree = 0;

	if (flag_degree == false){
		int vectorIndex = 0;
		int offsetInV = 0;
		int dataLength = 0;

		// Parameters Preparation
		vectorIndex = indexStructure.directIndex[vid]/64;
		offsetInV = indexStructure.directIndex[vid]%64;
		if ( vid != indexStructure.directIndex.size()-1){
			dataLength = indexStructure.directIndex[vid+1] - indexStructure.directIndex[vid];
		}else{
			dataLength = data.length - indexStructure.directIndex[vid];
		}

		phsim::BitStringReader reader(data, vectorIndex, dataLength , offsetInV); // reader.bits_left = dataLength + offsetInV
		int flag=0;
		while(!reader.eof()){
			if (flag==0){
				phsim::gammaDecodeWithOffset(reader , offsetInV);
				flag++;
				degree++;
			}else{
				phsim::gammaDecode(reader);
				degree++;
			}
		}

	}else{

		int vectorIndex = 0;
		int offsetInV = 0;
		int dataLength = 0;

		unsigned long temp = 0;
		// Parameters Preparation
		vectorIndex = indexStructure.directIndex[vid]/64;
		offsetInV = indexStructure.directIndex[vid]%64;
		if ( vid != indexStructure.directIndex.size()-1){
			dataLength = indexStructure.directIndex[vid+1] - indexStructure.directIndex[vid];
		}else{
			dataLength = data.length - indexStructure.directIndex[vid];
		}

		phsim::BitStringReader reader(data, vectorIndex, dataLength , offsetInV); // reader.bits_left = dataLength + offsetInV
		int flag=0;
		while(!reader.eof()){
			if (flag==0){
				degree = (unsigned long)phsim::gammaDecodeWithOffset(reader , offsetInV);
				break;
			}
		}
	}
	return degree;
}

int getDegreeUsingSemiDI(phsim::BitString data, unsigned int vid, bool flag_degree){

	int degree = 0;

	if (flag_degree == false){
		int vectorIndex = 0;
		int offsetInV = 0;
		int dataLength = 0;

		// Parameters Preparation
		int order = vid/4;
		if ( vid%4 == 0){      // The first word in semiDI is enough
			vectorIndex = indexStructure.semiDIFW[order]/64;
			offsetInV = indexStructure.semiDIFW[order]%64;
			dataLength = indexStructure.semiDISW[order][0].to_ulong();
			// Decoding
			phsim::BitStringReader reader(data, vectorIndex, dataLength , offsetInV); // reader.bits_left = dataLength + offsetInV
			int flag=0;
			while(!reader.eof()){
				if (flag==0){
					phsim::gammaDecodeWithOffset(reader , offsetInV);
					flag++;
					degree++;
				}else{
					phsim::gammaDecode(reader);
					degree++;
				}
			}
		}else{
			if ( vid%4 == 3){
				unsigned int baseAddress = indexStructure.semiDIFW[order];
				unsigned int realAddress = baseAddress + indexStructure.semiDISW[order][2].to_ulong();
				vectorIndex = realAddress/64;
				offsetInV = realAddress%64;
				dataLength = indexStructure.semiDIFW[order+1] - (indexStructure.semiDIFW[order] + indexStructure.semiDISW[order][2].to_ulong()) ;
			}else{
				unsigned int baseAddress = indexStructure.semiDIFW[order];
				unsigned int realAddress = baseAddress + indexStructure.semiDISW[order][vid % 4 -1].to_ulong();
				vectorIndex = realAddress/64;
				offsetInV = realAddress%64;
				dataLength = indexStructure.semiDISW[order][vid % 4].to_ulong() - indexStructure.semiDISW[order][vid % 4-1].to_ulong() ;
			}
			// Decoding
			phsim::BitStringReader reader(data, vectorIndex, dataLength , offsetInV); // reader.bits_left = dataLength + offsetInV
			int flag=0;
			while(!reader.eof()){
				if (flag==0){
					phsim::gammaDecodeWithOffset(reader , offsetInV);
					flag++;
					degree++;
				}else{
					phsim::gammaDecode(reader);
					degree++;
				}
			}
		}

	}else{
		int vectorIndex = 0;
		int offsetInV = 0;
		int dataLength = 0;



		// Parameters Preparation
		int order = vid/4;
		if ( vid%4 == 0){      // The first word in semiDI is enough
			vectorIndex = indexStructure.semiDIFW[order]/64;
			offsetInV = indexStructure.semiDIFW[order]%64;
			dataLength = indexStructure.semiDISW[order][0].to_ulong();
			// Decoding
			phsim::BitStringReader reader(data, vectorIndex, dataLength , offsetInV); // reader.bits_left = dataLength + offsetInV
			int flag=0;
			//while(!reader.eof()){
				//if (flag==0){
					degree = (unsigned long)phsim::gammaDecodeWithOffset(reader , offsetInV);
					//break;

				//}
			//}
		}else{
			if ( vid%4 == 3){
				unsigned int baseAddress = indexStructure.semiDIFW[order];
				unsigned int realAddress = baseAddress + indexStructure.semiDISW[order][2].to_ulong();
				vectorIndex = realAddress/64;
				offsetInV = realAddress%64;
				dataLength = indexStructure.semiDIFW[order+1] - (indexStructure.semiDIFW[order] + indexStructure.semiDISW[order][2].to_ulong()) ;
			}else{
				unsigned int baseAddress = indexStructure.semiDIFW[order];
				unsigned int realAddress = baseAddress + indexStructure.semiDISW[order][vid % 4 -1].to_ulong();
				vectorIndex = realAddress/64;
				offsetInV = realAddress%64;
				dataLength = indexStructure.semiDISW[order][vid % 4].to_ulong() - indexStructure.semiDISW[order][vid % 4-1].to_ulong() ;
			}
			// Decoding
			phsim::BitStringReader reader(data, vectorIndex, dataLength , offsetInV); // reader.bits_left = dataLength + offsetInV
			int flag=0;
			//while(!reader.eof()){
				//if (flag==0){
					degree = phsim::gammaDecodeWithOffset(reader , offsetInV);
					//break;
				//}
			//}
		}
	}

	return degree;
}
int getDegreeUsingIDI (phsim::BitString data, unsigned int vid ){
	int degree = 0;

	vector< int > adjList;
	vector< int > temp;
	int block_num = vid / indexStructure.blockVNum;
	int in_block_index = vid % indexStructure.blockVNum;

	int vectorIndex = 0;
	int offsetInV = 0;
	int dataLength = 0;

	rrr_vector<>::select_1_type rrr_s(&indexStructure.rrrBVofEachBlock[block_num]);
	rrr_vector<>::rank_1_type rrr_r(&indexStructure.rrrBVofEachBlock[block_num]);


	if ( rrr_r(indexStructure.blockVNum) == 1 ){

		dataLength = indexStructure.indrectIndex[block_num+1][0] - indexStructure.indrectIndex[block_num][0];
		vectorIndex = indexStructure.indrectIndex[block_num][0]/64;
		offsetInV = indexStructure.indrectIndex[block_num][0]%64;

	}else if ( rrr_r(indexStructure.blockVNum) == 2 ){

		if ( in_block_index < rrr_s(2)){
			dataLength = indexStructure.indrectIndex[block_num][1] - indexStructure.indrectIndex[block_num][0];
			vectorIndex = indexStructure.indrectIndex[block_num][0]/64;
			offsetInV = indexStructure.indrectIndex[block_num][0]%64;
		}else{
			dataLength = indexStructure.indrectIndex[block_num+1][0] - indexStructure.indrectIndex[block_num][1];
			vectorIndex = indexStructure.indrectIndex[block_num][1]/64;
			offsetInV = indexStructure.indrectIndex[block_num][1]%64;
		}

	}else if ( rrr_r(indexStructure.blockVNum) == 3 ){

		if ( in_block_index < rrr_s(2)){
			dataLength = indexStructure.indrectIndex[block_num][1] - indexStructure.indrectIndex[block_num][0];
			vectorIndex = indexStructure.indrectIndex[block_num][0]/64;
			offsetInV = indexStructure.indrectIndex[block_num][0]%64;

		}else if (in_block_index < rrr_s(3)){

			dataLength = indexStructure.indrectIndex[block_num][2] - indexStructure.indrectIndex[block_num][1];
			vectorIndex = indexStructure.indrectIndex[block_num][1]/64;
			offsetInV = indexStructure.indrectIndex[block_num][1]%64;

		}else{

			dataLength = indexStructure.indrectIndex[block_num+1][0] - indexStructure.indrectIndex[block_num][2];
			vectorIndex = indexStructure.indrectIndex[block_num][2]/64;
			offsetInV = indexStructure.indrectIndex[block_num][2]%64;
		}

	}else if ( rrr_r(indexStructure.blockVNum) == 4 ){

		if ( in_block_index < rrr_s(2)){

			dataLength = indexStructure.indrectIndex[block_num][1] - indexStructure.indrectIndex[block_num][0];
			vectorIndex = indexStructure.indrectIndex[block_num][0]/64;
			offsetInV = indexStructure.indrectIndex[block_num][0]%64;

		}else if (in_block_index < rrr_s(3)){

			dataLength = indexStructure.indrectIndex[block_num][2] - indexStructure.indrectIndex[block_num][1];
			vectorIndex = indexStructure.indrectIndex[block_num][1]/64;
			offsetInV = indexStructure.indrectIndex[block_num][1]%64;

		}else if (in_block_index < rrr_s(4)){

			dataLength = indexStructure.indrectIndex[block_num][3] - indexStructure.indrectIndex[block_num][2];
			vectorIndex = indexStructure.indrectIndex[block_num][2]/64;
			offsetInV = indexStructure.indrectIndex[block_num][2]%64;

		}else{
			dataLength = indexStructure.indrectIndex[block_num+1][0] - indexStructure.indrectIndex[block_num][3];
			vectorIndex = indexStructure.indrectIndex[block_num][3]/64;
			offsetInV = indexStructure.indrectIndex[block_num][3]%64;
		}

	}else{

	}

//
	phsim::BitStringReader reader(data, vectorIndex, dataLength , offsetInV); // reader.bits_left = dataLength + offsetInV
	int flag=0;
	while(!reader.eof()){
		if (flag==0){
			temp.push_back(phsim::gammaDecodeWithOffset(reader , offsetInV));
			flag++;
		}else{
			temp.push_back(phsim::gammaDecode(reader));
		}

	}

	if (in_block_index==0){
		degree = temp[0];


	}else{

		int start_index = 0;
		int degree_num = temp[0];
		start_index+=degree_num+1;
		int next_pos = 0;
		for (int i=1; i< in_block_index ; i++ ){
			next_pos += degree_num+1;
			degree_num = temp[next_pos];
			start_index += degree_num+1;
		}
		degree = temp[start_index];
	}

	return degree;
}

int getDegreeUsingrrrDI (phsim::BitString data, unsigned int vid, bool flag_degree){

	int degree=0;
	if (flag_degree == false){
		// Parameters initialization
		int vectorIndex = 0;
		int offsetInV = 0;
		int dataLength = 0;

		rrr_vector<>::select_1_type rrr_s(&indexStructure.rrr_direct_index);
		// Parameters Preparation
		vectorIndex = rrr_s(vid+1)/64;
		offsetInV = rrr_s(vid+1)%64;
		if ( vid != indexStructure.directIndex.size()-1){
			dataLength = rrr_s(vid+2) - rrr_s(vid+1);
		}else{
			dataLength = data.length - rrr_s(vid+1);
		}

		phsim::BitStringReader reader(data, vectorIndex, dataLength , offsetInV); // reader.bits_left = dataLength + offsetInV
		int flag=0;
		while(!reader.eof()){
			if (flag==0){
				phsim::gammaDecodeWithOffset(reader , offsetInV);
				flag++;
				degree++;
			}else{
				phsim::gammaDecode(reader);
				degree++;
			}

		}
	}else{
		// Parameters initialization
		int vectorIndex = 0;
		int offsetInV = 0;
		int dataLength = 0;


		rrr_vector<>::select_1_type rrr_s(&indexStructure.rrr_direct_index);
		// Parameters Preparation
		vectorIndex = rrr_s(vid+1)/64;
		offsetInV = rrr_s(vid+1)%64;
		if ( vid != indexStructure.directIndex.size()-1){
			dataLength = rrr_s(vid+2) - rrr_s(vid+1);
		}else{
			dataLength = data.length - rrr_s(vid+1);
		}

		phsim::BitStringReader reader(data, vectorIndex, dataLength , offsetInV); // reader.bits_left = dataLength + offsetInV
		int flag=0;
		//while(!reader.eof()){
			//if (flag==0){
				degree = (unsigned long)phsim::gammaDecodeWithOffset(reader , offsetInV);
			//	break;
			//}
		//}

	}
	return degree;
}
int getDegreeUsingsdbDI (phsim::BitString data, unsigned int vid, bool flag_degree){

	int degree=0;
	if (flag_degree == false){
		// Parameters initialization
		int vectorIndex = 0;
		int offsetInV = 0;
		int dataLength = 0;

		sd_vector<>::select_1_type sdb_s(&indexStructure.sdb_vector);
		// Parameters Preparation
		vectorIndex = sdb_s(vid+1)/64;
		offsetInV = sdb_s(vid+1)%64;
		if ( vid != indexStructure.directIndex.size()-1){
			dataLength = sdb_s(vid+2) - sdb_s(vid+1);
		}else{
			dataLength = data.length - sdb_s(vid+1);
		}

		phsim::BitStringReader reader(data, vectorIndex, dataLength , offsetInV); // reader.bits_left = dataLength + offsetInV
		int flag=0;
		while(!reader.eof()){
			if (flag==0){
				phsim::gammaDecodeWithOffset(reader , offsetInV);
				flag++;
				degree++;
			}else{
				phsim::gammaDecode(reader);
				degree++;
			}

		}
	}else{
		// Parameters initialization
		int vectorIndex = 0;
		int offsetInV = 0;
		int dataLength = 0;


		sd_vector<>::select_1_type sdb_s(&indexStructure.sdb_vector);
		// Parameters Preparation
		vectorIndex = sdb_s(vid+1)/64;
		offsetInV = sdb_s(vid+1)%64;
		if ( vid != indexStructure.directIndex.size()-1){
			dataLength = sdb_s(vid+2) - sdb_s(vid+1);
		}else{
			dataLength = data.length - sdb_s(vid+1);
		}

		phsim::BitStringReader reader(data, vectorIndex, dataLength , offsetInV); // reader.bits_left = dataLength + offsetInV
		int flag=0;
		//while(!reader.eof()){
			//if (flag==0){
				degree = (unsigned long)phsim::gammaDecodeWithOffset(reader , offsetInV);
			//	break;
			//}
		//}

	}
	return degree;
}


vector <int> getAdjForOneVUsingDI (phsim::BitString data, unsigned int vid  ){
	// Parameters initialization
	int vectorIndex = 0;
	int offsetInV = 0;
	int dataLength = 0;
	vector <int> adjList;
	vector <int> temp;
	// Parameters Preparation
	vectorIndex = indexStructure.directIndex[vid]/64;
	offsetInV = indexStructure.directIndex[vid]%64;
	if ( vid != indexStructure.directIndex.size()-1){
		dataLength = indexStructure.directIndex[vid+1] - indexStructure.directIndex[vid];
	}else{
		dataLength = data.length - indexStructure.directIndex[vid];
	}

	phsim::BitStringReader reader(data, vectorIndex, dataLength , offsetInV); // reader.bits_left = dataLength + offsetInV
	int flag=0;
	while(!reader.eof()){
		if (flag==0){
			temp.push_back(phsim::gammaDecodeWithOffset(reader , offsetInV));
			flag++;
		}else{
			temp.push_back(phsim::gammaDecode(reader));
		}

	}
	int base = 0;
	if ( negativeDiff[vid] == false){
		base = vid + temp[0];
	}else{
		base = vid - temp[0];
	}

	adjList.push_back(mappingForRestore[base]);
	for( int i=1; i<temp.size(); i++){
		int id = base + temp[i];
		adjList.push_back(mappingForRestore[id]);
		base = id;
	}

	return adjList;
}
vector <int> getAdjForOneVUsingDIWithD (phsim::BitString data, unsigned int vid){
	// Parameters initialization
	int vectorIndex = 0;
	int offsetInV = 0;
	int dataLength = 0;
	vector <int> adjList;
	vector <int> temp;
	// Parameters Preparation
	vectorIndex = indexStructure.directIndex[vid]/64;
	offsetInV = indexStructure.directIndex[vid]%64;
	if ( vid != indexStructure.directIndex.size()-1){
		dataLength = indexStructure.directIndex[vid+1] - indexStructure.directIndex[vid];
	}else{
		dataLength = data.length - indexStructure.directIndex[vid];
	}

	phsim::BitStringReader reader(data, vectorIndex, dataLength , offsetInV); // reader.bits_left = dataLength + offsetInV
	int flag=0;
	while(!reader.eof()){
		if (flag==0){
			temp.push_back(phsim::gammaDecodeWithOffset(reader , offsetInV));
			flag++;
		}else{
			temp.push_back(phsim::gammaDecode(reader));
		}

	}
	int base = 0;
	if ( negativeDiff[vid] == false){
		base = vid + temp[1];
	}else{
		base = vid - temp[1];
	}
	adjList.push_back(mappingForRestore[base]);
	for( int i=2; i<temp.size(); i++){
		int id = base + temp[i];
		adjList.push_back(mappingForRestore[id]);
		base = id;
	}

	return adjList;
}

vector <int> getAdjForOneVUsingSemiDI (phsim::BitString data, unsigned int vid){
	// Parameters initialization
	int vectorIndex = 0;
	int offsetInV = 0;
	int dataLength = 0;
	vector <int> adjList;
	vector <int> temp;
	// Parameters Preparation
	int order = vid/4;
	if ( vid%4 == 0){      // The first word in semiDI is enough
		vectorIndex = indexStructure.semiDIFW[order]/64;
		offsetInV = indexStructure.semiDIFW[order]%64;
		dataLength = indexStructure.semiDISW[order][0].to_ulong();
		// Decoding
		phsim::BitStringReader reader(data, vectorIndex, dataLength , offsetInV); // reader.bits_left = dataLength + offsetInV
		int flag=0;
		while(!reader.eof()){
			if (flag==0){
				temp.push_back(phsim::gammaDecodeWithOffset(reader , offsetInV));
				flag++;
			}else{
				temp.push_back(phsim::gammaDecode(reader));
			}
		}
	}else{
		if ( vid%4 == 3){
			unsigned int baseAddress = indexStructure.semiDIFW[order];
			unsigned int realAddress = baseAddress + indexStructure.semiDISW[order][2].to_ulong();
			vectorIndex = realAddress/64;
			offsetInV = realAddress%64;
			dataLength = indexStructure.semiDIFW[order+1] - (indexStructure.semiDIFW[order] + indexStructure.semiDISW[order][2].to_ulong()) ;
		}else{
			unsigned int baseAddress = indexStructure.semiDIFW[order];
			unsigned int realAddress = baseAddress + indexStructure.semiDISW[order][vid % 4 -1].to_ulong();
			vectorIndex = realAddress/64;
			offsetInV = realAddress%64;
			dataLength = indexStructure.semiDISW[order][vid % 4].to_ulong() - indexStructure.semiDISW[order][vid % 4-1].to_ulong() ;
		}
		// Decoding
		phsim::BitStringReader reader(data, vectorIndex, dataLength , offsetInV); // reader.bits_left = dataLength + offsetInV
		int flag=0;
		while(!reader.eof()){
			if (flag==0){
				temp.push_back(phsim::gammaDecodeWithOffset(reader , offsetInV));
				flag++;
			}else{
				temp.push_back(phsim::gammaDecode(reader));
			}
		}
	}
	// Building adjacency list for that vertex
	int base = 0;
	if ( negativeDiff[vid] == false){
		base = vid + temp[0];
	}else{
		base = vid - temp[0];
	}
	adjList.push_back(mappingForRestore[base]);
	for( int i=1; i<temp.size(); i++){
		int id = base + temp[i];
		adjList.push_back(mappingForRestore[id]);
		base = id;
	}

	return adjList;
}

vector <int> getAdjForOneVUsingSemiDIWithD (phsim::BitString data, unsigned int vid){
	// Parameters initialization
	int vectorIndex = 0;
	int offsetInV = 0;
	int dataLength = 0;
	vector <int> adjList;
	vector <int> temp;
	// Parameters Preparation
	int order = vid/4;
	if ( vid%4 == 0){      // The first word in semiDI is enough
		vectorIndex = indexStructure.semiDIFW[order]/64;
		offsetInV = indexStructure.semiDIFW[order]%64;
		dataLength = indexStructure.semiDISW[order][0].to_ulong();
		// Decoding
		phsim::BitStringReader reader(data, vectorIndex, dataLength , offsetInV); // reader.bits_left = dataLength + offsetInV
		int flag=0;
		while(!reader.eof()){
			if (flag==0){
				temp.push_back(phsim::gammaDecodeWithOffset(reader , offsetInV));
				flag++;
			}else{
				temp.push_back(phsim::gammaDecode(reader));
			}
		}
	}else{
		if ( vid%4 == 3){
			unsigned int baseAddress = indexStructure.semiDIFW[order];
			unsigned int realAddress = baseAddress + indexStructure.semiDISW[order][2].to_ulong();
			vectorIndex = realAddress/64;
			offsetInV = realAddress%64;
			dataLength = indexStructure.semiDIFW[order+1] - (indexStructure.semiDIFW[order] + indexStructure.semiDISW[order][2].to_ulong()) ;
		}else{
			unsigned int baseAddress = indexStructure.semiDIFW[order];
			unsigned int realAddress = baseAddress + indexStructure.semiDISW[order][vid % 4 -1].to_ulong();
			vectorIndex = realAddress/64;
			offsetInV = realAddress%64;
			dataLength = indexStructure.semiDISW[order][vid % 4].to_ulong() - indexStructure.semiDISW[order][vid % 4-1].to_ulong() ;
		}
		// Decoding
		phsim::BitStringReader reader(data, vectorIndex, dataLength , offsetInV); // reader.bits_left = dataLength + offsetInV
		int flag=0;
		while(!reader.eof()){
			if (flag==0){
				temp.push_back(phsim::gammaDecodeWithOffset(reader , offsetInV));
				flag++;
			}else{
				temp.push_back(phsim::gammaDecode(reader));
			}
		}
	}
	// Building adjacency list for that vertex
	int base = 0;
	if ( negativeDiff[vid] == false){
		base = vid + temp[1];
	}else{
		base = vid - temp[1];
	}
	adjList.push_back(mappingForRestore[base]);
	for( int i=2; i<temp.size(); i++){
		int id = base + temp[i];
		adjList.push_back(mappingForRestore[id]);
		base = id;
	}

	return adjList;
}

vector <int> getAdjForOneVUsingrrrDI (phsim::BitString data, unsigned int vid){
	// Parameters initialization
	int vectorIndex = 0;
	int offsetInV = 0;
	int dataLength = 0;
	vector <int> adjList;
	vector <int> temp;
	rrr_vector<>::select_1_type rrr_s(&indexStructure.rrr_direct_index);
	// Parameters Preparation
	vectorIndex = rrr_s(vid+1)/64;
	offsetInV = rrr_s(vid+1)%64;
	if ( vid != indexStructure.directIndex.size()-1){
		dataLength = rrr_s(vid+2) - rrr_s(vid+1);
	}else{
		dataLength = data.length - rrr_s(vid+1);
	}

	phsim::BitStringReader reader(data, vectorIndex, dataLength , offsetInV); // reader.bits_left = dataLength + offsetInV
	int flag=0;
	while(!reader.eof()){
		if (flag==0){
			temp.push_back(phsim::gammaDecodeWithOffset(reader , offsetInV));
			flag++;
		}else{
			temp.push_back(phsim::gammaDecode(reader));
		}

	}
	int base = 0;
	if ( negativeDiff[vid] == false){
		base = vid + temp[0];
	}else{
		base = vid - temp[0];
	}
	adjList.push_back(mappingForRestore[base]);
	for( int i=1; i<temp.size(); i++){
		int id = base + temp[i];
		adjList.push_back(mappingForRestore[id]);
		base = id;
	}
	return adjList;
}
vector <int> getAdjForOneVUsingrrrDIWithD (phsim::BitString data, unsigned int vid){
	// Parameters initialization
	int vectorIndex = 0;
	int offsetInV = 0;
	int dataLength = 0;
	vector <int> adjList;
	vector <int> temp;
	rrr_vector<>::select_1_type rrr_s(&indexStructure.rrr_direct_index);
	// Parameters Preparation
	vectorIndex = rrr_s(vid+1)/64;
	offsetInV = rrr_s(vid+1)%64;
	if ( vid != indexStructure.directIndex.size()-1){
		dataLength = rrr_s(vid+2) - rrr_s(vid+1);
	}else{
		dataLength = data.length - rrr_s(vid+1);
	}

	phsim::BitStringReader reader(data, vectorIndex, dataLength , offsetInV); // reader.bits_left = dataLength + offsetInV
	int flag=0;
	while(!reader.eof()){
		if (flag==0){
			temp.push_back(phsim::gammaDecodeWithOffset(reader , offsetInV));
			flag++;
		}else{
			temp.push_back(phsim::gammaDecode(reader));
		}

	}
	int base = 0;
	if ( negativeDiff[vid] == false){
		base = vid + temp[1];
	}else{
		base = vid - temp[1];
	}
	adjList.push_back(mappingForRestore[base]);
	for( int i=2; i<temp.size(); i++){
		int id = base + temp[i];
		adjList.push_back(mappingForRestore[id]);
		base = id;
	}
	return adjList;
}
vector <int> getAdjForOneVUsingsdbDI (phsim::BitString data, unsigned int vid){
	// Parameters initialization
	int vectorIndex = 0;
	int offsetInV = 0;
	int dataLength = 0;
	vector <int> adjList;
	vector <int> temp;
	sd_vector<>::select_1_type sdb_s(&indexStructure.sdb_vector);
	// Parameters Preparation
	vectorIndex = sdb_s(vid+1)/64;
	offsetInV = sdb_s(vid+1)%64;

	if ( vid != indexStructure.directIndex.size()-1){
		dataLength = sdb_s(vid+2) - sdb_s(vid+1);
	}else{
		dataLength = data.length - sdb_s(vid+1);
	}

	phsim::BitStringReader reader(data, vectorIndex, dataLength , offsetInV); // reader.bits_left = dataLength + offsetInV
	int flag=0;
	while(!reader.eof()){
		if (flag==0){
			temp.push_back(phsim::gammaDecodeWithOffset(reader , offsetInV));
			flag++;
		}else{
			temp.push_back(phsim::gammaDecode(reader));
		}

	}
	int base = 0;
	if ( negativeDiff[vid] == false){
		base = vid + temp[0];
	}else{
		base = vid - temp[0];
	}
	adjList.push_back(mappingForRestore[base]);
	for( int i=1; i<temp.size(); i++){
		int id = base + temp[i];
		adjList.push_back(mappingForRestore[id]);
		base = id;
	}
	return adjList;
}
vector <int> getAdjForOneVUsingsdbDIWithD (phsim::BitString data, unsigned int vid){
	// Parameters initialization
	int vectorIndex = 0;
	int offsetInV = 0;
	int dataLength = 0;
	vector <int> adjList;
	vector <int> temp;
	sd_vector<>::select_1_type sdb_s(&indexStructure.sdb_vector);
	// Parameters Preparation
	vectorIndex = sdb_s(vid+1)/64;
	offsetInV = sdb_s(vid+1)%64;

	if ( vid != indexStructure.directIndex.size()-1){
		dataLength = sdb_s(vid+2) - sdb_s(vid+1);
	}else{
		dataLength = data.length - sdb_s(vid+1);
	}

	phsim::BitStringReader reader(data, vectorIndex, dataLength , offsetInV); // reader.bits_left = dataLength + offsetInV
	int flag=0;
	while(!reader.eof()){
		if (flag==0){
			temp.push_back(phsim::gammaDecodeWithOffset(reader , offsetInV));
			flag++;
		}else{
			temp.push_back(phsim::gammaDecode(reader));
		}

	}
	int base = 0;
	if ( negativeDiff[vid] == false){
		base = vid + temp[1];
	}else{
		base = vid - temp[1];
	}
	adjList.push_back(mappingForRestore[base]);
	for( int i=2; i<temp.size(); i++){
		int id = base + temp[i];
		adjList.push_back(mappingForRestore[id]);
		base = id;
	}
	return adjList;
}

vector <int> getAdjForOneVUsingIDI (phsim::BitString data,unsigned int vid){

	vector< int > adjList;
	vector< int > temp;
	int block_num = vid / indexStructure.blockVNum;
	int in_block_index = vid % indexStructure.blockVNum;

	int vectorIndex = 0;
	int offsetInV = 0;
	int dataLength = 0;

	int offset_in_subblock = 0;
	int case_flag = 0;

	rrr_vector<>::select_1_type rrr_select(&indexStructure.rrrBVofEachBlock[block_num]);
	rrr_vector<>::rank_1_type rrr_rank(&indexStructure.rrrBVofEachBlock[block_num]);


	if ( indexStructure.indrectIndex[block_num].size() == 1 ){

		dataLength = indexStructure.indrectIndex[block_num+1][0] - indexStructure.indrectIndex[block_num][0];
		vectorIndex = indexStructure.indrectIndex[block_num][0]/64;
		offsetInV = indexStructure.indrectIndex[block_num][0]%64;
		offset_in_subblock = in_block_index;
		//case_flag = 1;

	}else if ( indexStructure.indrectIndex[block_num].size() == 2 ){

		if ( in_block_index < rrr_select(2)){
			dataLength = indexStructure.indrectIndex[block_num][1] - indexStructure.indrectIndex[block_num][0];
			vectorIndex = indexStructure.indrectIndex[block_num][0]/64;
			offsetInV = indexStructure.indrectIndex[block_num][0]%64;
			offset_in_subblock = in_block_index;
			//case_flag = 2;

		}else{
			dataLength = indexStructure.indrectIndex[block_num+1][0] - indexStructure.indrectIndex[block_num][1];
			vectorIndex = indexStructure.indrectIndex[block_num][1]/64;
			offsetInV = indexStructure.indrectIndex[block_num][1]%64;
			offset_in_subblock = in_block_index - rrr_select(2);
			//case_flag = 2;
		}

	}else if ( indexStructure.indrectIndex[block_num].size() == 3 ){

		if ( in_block_index < rrr_select(2)){
			dataLength = indexStructure.indrectIndex[block_num][1] - indexStructure.indrectIndex[block_num][0];
			vectorIndex = indexStructure.indrectIndex[block_num][0]/64;
			offsetInV = indexStructure.indrectIndex[block_num][0]%64;
			offset_in_subblock = in_block_index;
			//case_flag = 3;

		}else if (in_block_index < rrr_select(3)){

			dataLength = indexStructure.indrectIndex[block_num][2] - indexStructure.indrectIndex[block_num][1];
			vectorIndex = indexStructure.indrectIndex[block_num][1]/64;
			offsetInV = indexStructure.indrectIndex[block_num][1]%64;
			offset_in_subblock = in_block_index - rrr_select(2);
			//case_flag = 3;
		}else{

			dataLength = indexStructure.indrectIndex[block_num+1][0] - indexStructure.indrectIndex[block_num][2];
			vectorIndex = indexStructure.indrectIndex[block_num][2]/64;
			offsetInV = indexStructure.indrectIndex[block_num][2]%64;
			offset_in_subblock = in_block_index - rrr_select(3);
			//case_flag = 3;
		}

	}else if ( indexStructure.indrectIndex[block_num].size() == 4 ){

		if ( in_block_index < rrr_select(2)){

			dataLength = indexStructure.indrectIndex[block_num][1] - indexStructure.indrectIndex[block_num][0];
			vectorIndex = indexStructure.indrectIndex[block_num][0]/64;
			offsetInV = indexStructure.indrectIndex[block_num][0]%64;
			offset_in_subblock = in_block_index;
			//case_flag = 4;

		}else if (in_block_index < rrr_select(3)){

			dataLength = indexStructure.indrectIndex[block_num][2] - indexStructure.indrectIndex[block_num][1];
			vectorIndex = indexStructure.indrectIndex[block_num][1]/64;
			offsetInV = indexStructure.indrectIndex[block_num][1]%64;
			offset_in_subblock = in_block_index - rrr_select(2);
			//case_flag = 4;

		}else if (in_block_index < rrr_select(4)){

			dataLength = indexStructure.indrectIndex[block_num][3] - indexStructure.indrectIndex[block_num][2];
			vectorIndex = indexStructure.indrectIndex[block_num][2]/64;
			offsetInV = indexStructure.indrectIndex[block_num][2]%64;
			offset_in_subblock = in_block_index - rrr_select(3);
			//case_flag = 4;
		}else{
			dataLength = indexStructure.indrectIndex[block_num+1][0] - indexStructure.indrectIndex[block_num][3];
			vectorIndex = indexStructure.indrectIndex[block_num][3]/64;
			offsetInV = indexStructure.indrectIndex[block_num][3]%64;
			offset_in_subblock = in_block_index - rrr_select(4);
			//case_flag = 4;
		}

	}else{
		//cout << "here" <<endl;
	}
	//cout << "111" <<endl;
//
	phsim::BitStringReader reader(data, vectorIndex, dataLength , offsetInV); // reader.bits_left = dataLength + offsetInV
	int flag=0;
	while(!reader.eof()){
		if (flag==0){
			temp.push_back(phsim::gammaDecodeWithOffset(reader , offsetInV));
			flag++;
		}else{
			temp.push_back(phsim::gammaDecode(reader));
		}

	}


	if (offset_in_subblock == 0){
		int degree = temp[0];
		int base = 0;

		if ( negativeDiff[vid] == false){
			base = vid + temp[1];
		}else{
			base = vid - temp[1];
		}
		//cout << "degree : " <<degree << endl;

		adjList.push_back(mappingForRestore[base]);
		for( int i=2; i<degree+1; i++){
			int id = base + temp[i];
			adjList.push_back(mappingForRestore[id]);
			base = id;
		}
	}else{
		int start_index = 0;
		int degree_num = temp[0];
		start_index+=degree_num+1;
		int next_pos = 0;
		for (int i=1; i< offset_in_subblock ; i++ ){
			next_pos += degree_num+1;
			degree_num = temp[next_pos];
			start_index += degree_num+1;
		}

		int degree = temp[start_index];
		int base = 0;

		if ( negativeDiff[vid] == false){
			base = vid + temp[start_index + 1];
		}else{
			base = vid - temp[start_index + 1];
		}
		//cout << "degree : " <<degree << endl;

		adjList.push_back(mappingForRestore[base]);
		for( int i=2; i<degree+1; i++){
			int id = base + temp[start_index + i];
			adjList.push_back(mappingForRestore[id]);
			base = id;
		}

	}

	return adjList;
}


void DFS(phsim::BitString data){
	vector<int> thisLoop;
	vector<int> nextLoop;
	bool *visited = new bool [nV];

	for(int i=0; i<nV ;i++){
		visited[i] = false;
	}

	std::chrono::time_point<std::chrono::system_clock> start1, end1;
	std::chrono::duration<double> elapsed_seconds;
	//start1 = std::chrono::system_clock::now();

////////////////////////////////////////////////////////////////////////////////////////////////
	visited[0] = true;
	//visited[nV-1] = true;
	thisLoop.push_back(0);
	vector<int> temp_set = getAdjForOneVUsingDI(data, 0);
	for( int i=0; i<temp_set.size();i++){
		nextLoop.push_back(temp_set[i]);
		visited[temp_set[i]] = true;
	}


	start1 = std::chrono::system_clock::now();
//	for ( int i=0; i< nextLoop.size();i++){
//		cout << nextLoop[i] << " ";
//	}
//	while (nextLoop.size()){
//		//cout << nextLoop.size()<< " ";
//		thisLoop = nextLoop;
//		nextLoop.clear();
//		for ( int i=0; i<thisLoop.size() ;i++){
//			temp_set = getAdjForOneVUsingDI(data, thisLoop[i], false);
//			for ( int j=0; j<temp_set.size() ;j++){
//				if (visited[temp_set[j]] == false){
//					nextLoop.push_back(temp_set[j]);
//					visited[temp_set[j]] = true;
//				}
//			}
//		}
//	}

	while (nextLoop.size()){
		//cout << nextLoop.size()<<"  ";
		thisLoop = nextLoop;
		nextLoop.clear();
		for ( int i=0; i<thisLoop.size() ;i++){

			int vectorIndex = 0;
			int offsetInV = 0;
			int dataLength = 0;
			vector <int> adjList;
			//vector <int> temp;
			// Parameters Preparation
			vectorIndex = indexStructure.directIndex[thisLoop[i]]/64;
			offsetInV = indexStructure.directIndex[thisLoop[i]]%64;
			if ( thisLoop[i] != indexStructure.directIndex.size()-1){
				dataLength = indexStructure.directIndex[thisLoop[i]+1] - indexStructure.directIndex[thisLoop[i]];
			}else{
				dataLength = data.length - indexStructure.directIndex[thisLoop[i]];
			}

			phsim::BitStringReader reader(data, vectorIndex, dataLength , offsetInV); // reader.bits_left = dataLength + offsetInV
			int flag=0;
//			while(!reader.eof()){
//				if (flag==0){
//					temp.push_back(phsim::gammaDecodeWithOffset(reader , offsetInV));
//					flag++;
//				}else{
//					temp.push_back(phsim::gammaDecode(reader));
//				}
//			}
			int base = 0;
			int diff = 0;
			while(!reader.eof()){
				if (flag==0){
					diff = phsim::gammaDecodeWithOffset(reader , offsetInV);
					flag++;
					if ( negativeDiff[thisLoop[i]] == false){
						base = thisLoop[i] + diff;
					}else{
						base = thisLoop[i] - diff;
					}

					if (visited[base] == false){
						nextLoop.push_back(base);
						visited[base] = true;
					}

				}else{
					diff = phsim::gammaDecode(reader);
					int id = base + diff;
					if (visited[id] == false){
						nextLoop.push_back(id);
						visited[id] = true;
					}
					base = id;
				}
			}
//			int base = 0;
//			if ( negativeDiff[thisLoop[i]] == false){
//				base = thisLoop[i] + temp[0];
//			}else{
//				base = thisLoop[i] - temp[0];
//			}
//
//			if (visited[base] == false){
//				nextLoop.push_back(base);
//				visited[base] = true;
//			}
//			//adjList.push_back(mappingForRestore[base]);
//			for( int i=1; i<temp.size(); i++){
//				int id = base + temp[i];
//
//				if (visited[id] == false){
//					nextLoop.push_back(id);
//					visited[id] = true;
//				}
//				base = id;
//			}
		}


	}


////////////////////////////////////////////////////////////////////////////////////////////////


	end1 = std::chrono::system_clock::now();
	elapsed_seconds = end1-start1;
	std::cout << "DFS (Direct) elapsed time: " << elapsed_seconds.count() << "s\n";

	thisLoop.clear();
	nextLoop.clear();
	start1 = std::chrono::system_clock::now();

//////////////////////////////////////////////////////////////////////////////////////////

	cout << endl;

	bool *visited2 = new bool [nV];
	for(int i=0; i<nV ;i++){
		visited2[i] = false;
	}
	thisLoop.push_back(0);
	visited2[0] = true;

	for( int i=0; i<InputXadj[1];i++){
		nextLoop.push_back(Inputadj[i]);
		visited2[Inputadj[i]] = true;
	}
//	for ( int i=0; i< nextLoop.size();i++){
//		cout << nextLoop[i] << " ";
//	}
	while (nextLoop.size()){
		//cout << nextLoop.size()<< " ";
		thisLoop = nextLoop;
		nextLoop.clear();
		for ( int i=0; i<thisLoop.size() ;i++){

			for ( int j=0; j< InputXadj[ thisLoop[i]+1 ] - InputXadj[thisLoop[i]] ;j++){

				if (visited2[Inputadj[InputXadj[thisLoop[i]]+j]] == false){

					nextLoop.push_back( Inputadj[InputXadj[thisLoop[i]]+j] );
					visited2[ Inputadj[InputXadj[thisLoop[i]]+j] ] = true;

				}

			}
		}
	}

//////////////////////////////////////////////////////////////////////////////////////
	end1 = std::chrono::system_clock::now();
	elapsed_seconds = end1-start1;
	std::cout << "DFS (array-based adjacency list) elapsed time: " << elapsed_seconds.count() << "s\n";

	thisLoop.clear();
	nextLoop.clear();
	cout << endl;

///////////////////////////////////////////////////////////////////////////////////////
	bool *visited3 = new bool [nV];

	for(int i=0; i<nV ;i++){
		visited3[i] = false;
	}
	visited3[0] = true;
	visited3[nV-1] = true;
	visited3[nV-2] = true;
	visited3[nV-3] = true;

	thisLoop.push_back(0);
	temp_set.clear();
	temp_set = getAdjForOneVUsingSemiDI(data, 0);
	for( int i=0; i<temp_set.size();i++){
		nextLoop.push_back(temp_set[i]);
		visited3[temp_set[i]] = true;
	}

	//cout << "here" <<endl;
	start1 = std::chrono::system_clock::now();
	//int count = 0;
	while (nextLoop.size()){
	//while (count <=44){
		//cout << nextLoop.size() <<"  ";
		thisLoop = nextLoop;
		nextLoop.clear();
		for ( int i=0; i<thisLoop.size() ;i++){
			if ( thisLoop[i] == nV-1 || thisLoop[i] == nV-2 || thisLoop[i] == nV-3 ){

			}else{
				int vectorIndex = 0;
				int offsetInV = 0;
				int dataLength = 0;
				vector <int> adjList;
				vector <int> temp;
				// Parameters Preparation
				int order = thisLoop[i]/4;
				if ( thisLoop[i]%4 == 0){      // The first word in semiDI is enough
					vectorIndex = indexStructure.semiDIFW[order]/64;
					offsetInV = indexStructure.semiDIFW[order]%64;
					dataLength = indexStructure.semiDISW[order][0].to_ulong();
					// Decoding
					phsim::BitStringReader reader(data, vectorIndex, dataLength , offsetInV); // reader.bits_left = dataLength + offsetInV
					int flag=0;
//					while(!reader.eof()){
//						if (flag==0){
//							temp.push_back(phsim::gammaDecodeWithOffset(reader , offsetInV));
//							flag++;
//						}else{
//							temp.push_back(phsim::gammaDecode(reader));
//						}
//					}
					int base = 0;
					int diff = 0;
					while(!reader.eof()){
						if (flag==0){
							diff = phsim::gammaDecodeWithOffset(reader , offsetInV);
							flag++;
							if ( negativeDiff[thisLoop[i]] == false){
								base = thisLoop[i] + diff;
							}else{
								base = thisLoop[i] - diff;
							}

							if (visited3[base] == false){
								nextLoop.push_back(base);
								visited3[base] = true;
							}

						}else{
							diff = phsim::gammaDecode(reader);
							int id = base + diff;
							if (visited3[id] == false){
								nextLoop.push_back(id);
								visited3[id] = true;
							}
							base = id;
						}
					}
				}else{
					if ( thisLoop[i]%4 == 3){
						unsigned int baseAddress = indexStructure.semiDIFW[order];
						unsigned int realAddress = baseAddress + indexStructure.semiDISW[order][2].to_ulong();
						vectorIndex = realAddress/64;
						offsetInV = realAddress%64;
						dataLength = indexStructure.semiDIFW[order+1] - (indexStructure.semiDIFW[order] + indexStructure.semiDISW[order][2].to_ulong()) ;
					}else{
						unsigned int baseAddress = indexStructure.semiDIFW[order];
						unsigned int realAddress = baseAddress + indexStructure.semiDISW[order][thisLoop[i] % 4 -1].to_ulong();
						vectorIndex = realAddress/64;
						offsetInV = realAddress%64;
						dataLength = indexStructure.semiDISW[order][thisLoop[i] % 4].to_ulong() - indexStructure.semiDISW[order][thisLoop[i] % 4-1].to_ulong() ;
					}
					// Decoding
					phsim::BitStringReader reader(data, vectorIndex, dataLength , offsetInV); // reader.bits_left = dataLength + offsetInV
					int flag=0;
//					while(!reader.eof()){
//						if (flag==0){
//							temp.push_back(phsim::gammaDecodeWithOffset(reader , offsetInV));
//							flag++;
//						}else{
//							temp.push_back(phsim::gammaDecode(reader));
//						}
//					}
					int base = 0;
					int diff = 0;
					while(!reader.eof()){
						if (flag==0){
							diff = phsim::gammaDecodeWithOffset(reader , offsetInV);
							flag++;
							if ( negativeDiff[thisLoop[i]] == false){
								base = thisLoop[i] + diff;
							}else{
								base = thisLoop[i] - diff;
							}

							if (visited3[base] == false){
								nextLoop.push_back(base);
								visited3[base] = true;
							}

						}else{
							diff = phsim::gammaDecode(reader);
							int id = base + diff;
							if (visited3[id] == false){
								nextLoop.push_back(id);
								visited3[id] = true;
							}
							base = id;
						}
					}
				}

			}
		}

	}

////////////////////////////////////////////////////////////////////////////////////////////////


		end1 = std::chrono::system_clock::now();
		elapsed_seconds = end1-start1;
		std::cout << "DFS (Semi-Direct) elapsed time: " << elapsed_seconds.count() << "s\n";
		thisLoop.clear();
		nextLoop.clear();
		cout << endl;

///////////////////////////////////////////////////////////////////////////////////////////////

		bool *visited4 = new bool [nV];

		for(int i=0; i<nV ;i++){
			visited4[i] = false;
		}
		visited4[0] = true;
		visited4[nV-1] = true;
		visited4[nV-2] = true;
		thisLoop.push_back(0);
		temp_set.clear();
		temp_set = getAdjForOneVUsingrrrDI(data, 0);
		for( int i=0; i<temp_set.size();i++){
			nextLoop.push_back(temp_set[i]);
			visited4[temp_set[i]] = true;
		}
		start1 = std::chrono::system_clock::now();
		//int count = 0;
		while (nextLoop.size()){
			thisLoop = nextLoop;
			nextLoop.clear();
			for ( int i=0; i<thisLoop.size() ;i++){
				int vectorIndex = 0;
				int offsetInV = 0;
				int dataLength = 0;
				vector <int> adjList;
				vector <int> temp;
				rrr_vector<>::select_1_type rrr_s(&indexStructure.rrr_direct_index);
				// Parameters Preparation
				vectorIndex = rrr_s(thisLoop[i]+1)/64;
				offsetInV = rrr_s(thisLoop[i]+1)%64;
				if ( thisLoop[i] != indexStructure.directIndex.size()-1){
					dataLength = rrr_s(thisLoop[i]+2) - rrr_s(thisLoop[i]+1);
				}else{
					dataLength = data.length - rrr_s(thisLoop[i]+1);
				}

				phsim::BitStringReader reader(data, vectorIndex, dataLength , offsetInV); // reader.bits_left = dataLength + offsetInV
				int flag=0;
				int base = 0;
				int diff = 0;
				while(!reader.eof()){
					if (flag==0){
						diff = phsim::gammaDecodeWithOffset(reader , offsetInV);
						flag++;
						if ( negativeDiff[thisLoop[i]] == false){
							base = thisLoop[i] + diff;
						}else{
							base = thisLoop[i] - diff;
						}

						if (visited4[base] == false){
							nextLoop.push_back(base);
							visited4[base] = true;
						}

					}else{
						diff = phsim::gammaDecode(reader);
						int id = base + diff;
						if (visited4[id] == false){
							nextLoop.push_back(id);
							visited4[id] = true;
						}
						base = id;
					}
				}
			}

		}
////////////////////////////////////////////////////////////////////////////////////////////////


		end1 = std::chrono::system_clock::now();
		elapsed_seconds = end1-start1;
		std::cout << "DFS (rrr) elapsed time: " << elapsed_seconds.count() << "s\n";
		thisLoop.clear();
		nextLoop.clear();
		cout << endl;
///////////////////////////////////////////////////////////////////////////////////////////////

		bool *visited5 = new bool [nV];

		for(int i=0; i<nV ;i++){
			visited5[i] = false;
		}
		visited5[0] = true;
		visited5[nV-1] = true;
		visited5[nV-2] = true;
		thisLoop.push_back(0);
		temp_set.clear();
		temp_set = getAdjForOneVUsingsdbDI(data, 0);
		for( int i=0; i<temp_set.size();i++){
			nextLoop.push_back(temp_set[i]);
			visited5[temp_set[i]] = true;
		}
		start1 = std::chrono::system_clock::now();
		//int count = 0;
		while (nextLoop.size()){
			thisLoop = nextLoop;
			nextLoop.clear();
			for ( int i=0; i<thisLoop.size() ;i++){
				int vectorIndex = 0;
				int offsetInV = 0;
				int dataLength = 0;
				vector <int> adjList;
				vector <int> temp;
				sd_vector<>::select_1_type sdb_s(&indexStructure.sdb_vector);
				// Parameters Preparation
				vectorIndex = sdb_s(thisLoop[i]+1)/64;
				offsetInV = sdb_s(thisLoop[i]+1)%64;

				if ( thisLoop[i] != indexStructure.directIndex.size()-1){
					dataLength = sdb_s(thisLoop[i]+2) - sdb_s(thisLoop[i]+1);
				}else{
					dataLength = data.length - sdb_s(thisLoop[i]+1);
				}

				phsim::BitStringReader reader(data, vectorIndex, dataLength , offsetInV); // reader.bits_left = dataLength + offsetInV
				int flag=0;
				int base = 0;
				int diff = 0;
				while(!reader.eof()){
					if (flag==0){
						diff = phsim::gammaDecodeWithOffset(reader , offsetInV);
						flag++;
						if ( negativeDiff[thisLoop[i]] == false){
							base = thisLoop[i] + diff;
						}else{
							base = thisLoop[i] - diff;
						}

						if (visited5[base] == false){
							nextLoop.push_back(base);
							visited5[base] = true;
						}

					}else{
						diff = phsim::gammaDecode(reader);
						int id = base + diff;
						if (visited5[id] == false){
							nextLoop.push_back(id);
							visited5[id] = true;
						}
						base = id;
					}
				}
			}

		}
////////////////////////////////////////////////////////////////////////////////////////////////

		end1 = std::chrono::system_clock::now();
		elapsed_seconds = end1-start1;
		std::cout << "DFS (sd vector) elapsed time: " << elapsed_seconds.count() << "s\n";
		thisLoop.clear();
		nextLoop.clear();
		cout << endl;
}


int metisUsage(){

	//
	//	idx_t nVertices = 15;
	//	idx_t nEdges    = 22;
	//    idx_t xadj[nVertices+1] = {0, 2, 5, 8, 11, 13, 16, 20, 24, 28,31, 33, 36, 39, 42, 44};
	//    idx_t adjncy[2*nEdges] = {1,5,0, 2, 6, 1, 3, 7, 2, 4, 8, 3, 9, 0, 6, 10, 1, 5, 7, 11, 2, 6, 8, 12, 3, 7, 9, 13, 4, 8, 14, 5,
	//        		11, 6, 10, 12, 7, 11, 13, 8, 12, 14, 9, 13};

    idx_t nVertices = 6;
    idx_t nEdges    = 5;
    idx_t nWeights  = 1;
    idx_t nParts    = 2;

    idx_t objval;
    idx_t part[nVertices];


    // Indexes of starting points in adjacent array
    idx_t xadj[nVertices+1] = {0,2,4,5,8,9,10};

    // Adjacent vertices in consecutive index order
    idx_t adjncy[2 * nEdges] = {2,3,3,4,0,5,0,1,1,3};

    // Weights of vertices
    // if all weights are equal then can be set to NULL
    idx_t vwgt[nVertices * nWeights];


    int ret = METIS_PartGraphRecursive(&nVertices,& nWeights, xadj, adjncy,
     				       NULL, NULL, NULL, &nParts, NULL,
     				       NULL, NULL, &objval, part);

    /*
    int ret = METIS_PartGraphKway(&nVertices,& nWeights, xadj, adjncy,
				       NULL, NULL, NULL, &nParts, NULL,
				       NULL, NULL, &objval, part);
	*/
    //std::cout << ret << std::endl;
    cout<< "------------------------------------------" <<endl;
    for(unsigned part_i = 0; part_i < nVertices; part_i++){
	std::cout << part_i << " " << part[part_i] << std::endl;
    }
    cout<< "------------------------------------------" <<endl;
    MetisNode mynode( nVertices ,xadj,adjncy);



    //			cout << "L N1  : " ;
    //			for (const auto &x : nextN1L){
    //				cout << x << "  ";
    //			}
    //			cout << endl;
    //			cout << "L N2  : " ;
    //			for (const auto &x : nextN2L){
    //				cout << x << "  ";
    //			}
    //			cout << endl;
    //			cout << "R N1  : " ;
    //			for (const auto &x : nextN1R){
    //				cout << x << "  ";
    //			}
    //			cout << endl;
    //			cout << "R N2  : " ;
    //			for (const auto &x : nextN2R){
    //				cout << x << "  ";
    //			}
    //			cout<< endl;;
    //			cout <<"-----------------------------------" << endl;






    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();




    end = std::chrono::system_clock::now();

    std::chrono::duration<double> elapsed_seconds = end-start;
    std::time_t end_time = std::chrono::system_clock::to_time_t(end);

    std::cout << "finished computation at " << std::ctime(&end_time)
              << "elapsed time: " << elapsed_seconds.count() << "s\n";



    return 0;
}


