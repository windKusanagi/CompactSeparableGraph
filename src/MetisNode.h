/*
 * MetisNode.h
 *
 *  Created on: Nov 28, 2016
 *      Author: xiangzhang
 */

#ifndef METISNODE_H_
#define METISNODE_H_

#include <metis.h>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>

using namespace std;

class MetisNode{

public:
	MetisNode(idx_t NumOfVertice,idx_t xadj[] ,idx_t adj[]);
	MetisNode(idx_t nVertices , idx_t xadj[], idx_t adjncy[] , int **mapOrg );
	void printItself();
	void printMapping();
	void printXadj();
	void printAdj();

	idx_t NumberOfVL;
	idx_t NumberOfVR;
	std::vector <std::string> EdgeSeparator;

    idx_t * adjl;
    idx_t * xal;
    idx_t * adjr;
    idx_t * xar;
    //int ** mapLC;
    //int ** mapRC;
    int ** mapOrgLC;
    int ** mapOrgRC;
    idx_t *adjAftMLC;
	idx_t *adjAftMRC;

};

string IntToString (int a)
{
    ostringstream temp;
    temp<<a;
    return temp.str();
}

MetisNode::MetisNode(idx_t nVertices , idx_t xadj[], idx_t adjncy[] ){

    idx_t nWeights  = 1;
    idx_t nParts    = 2;

    idx_t objval; //number of edge separators
    idx_t part[nVertices];

    int ret = METIS_PartGraphRecursive(&nVertices,& nWeights, xadj, adjncy,
     				       NULL, NULL, NULL, &nParts, NULL,
     				       NULL, NULL, &objval, part);
/*
    int ret = METIS_PartGraphKway(&nVertices,& nWeights, xadj, adjncy,
				       NULL, NULL, NULL, &nParts, NULL,
				       NULL, NULL, &objval, part);
													*/
    //cout << " step 1 complete" <<endl;

    idx_t nVerticesLC = NULL;
    idx_t nVerticesRC = NULL;
    idx_t nEdgesLC = NULL;
    idx_t nEdgesRC = NULL;

    /*
    for(unsigned part_i = 0; part_i < nVertices; part_i++){
    	std::cout << part_i << " " << part[part_i] << std::endl;
    }*/

    //int edges [objval*2];

    for( int i = 0 ; i < nVertices; i++){
       	if (part[i] == 0){
       		nVerticesLC ++;

       		for(int j = xadj[i]; j < xadj[i+1]; j++){

       			if(part[adjncy[j]] == 1){

       				string temp = IntToString(i)+","+ IntToString(adjncy[j]);
       				//cout << temp <<endl;
       				this->EdgeSeparator.push_back(temp);

       				//cout<<"{" << i <<" , "<< adjncy[j] <<"}" << endl;
       			}
       		}
       	}
    }

    nVerticesRC = nVertices - nVerticesLC;

    //cout<< " step 2 complete ! " <<endl;

    // building mapping table of original vertex id for both children
    this->mapOrgLC = new int*[nVerticesLC];
    for(int i = 0; i < nVerticesLC; i++){
       this->mapOrgLC[i] = new int[2];
    }

    this->mapOrgRC = new int*[nVerticesRC];
    for(int i = 0; i < nVerticesRC; i++){
       this->mapOrgRC[i] = new int[2];
    }

    for( int i = 0 , flagL = 0 , flagR = 0 ; i < nVertices; i++){
       	if (part[i] == 0){
       		this->mapOrgLC[flagL][0] = i;
       		this->mapOrgLC[flagL][1] = flagL;
       		flagL++;

       	}else{
       		this->mapOrgRC[flagR][0] = i;
       		this->mapOrgRC[flagR][1] = flagR;
       		flagR++;
       	}
    }

    //cout << "step3 complete !" <<endl;


    vector <int> adjncyl;
    vector <int> xadjl;
    vector <int> adjncyr;
    vector <int> xadjr;

    xadjl.push_back(0);
    xadjr.push_back(0);
       //cout << nVerticesLC<< "," << nVerticesRC<< endl;


    int flagl=0;
    int flagr=0;
    for(int i = 0; i < nVertices; i++){   		//Generate adjncy and xadj for both children
    	if(part[i]==0){

    		for (int j = xadj[i]; j < xadj[i+1]; j++){
       			if(part[adjncy[j]] == 0){
       				adjncyl.push_back(adjncy[j]);
       				flagl++;
       			}
       		}
       		xadjl.push_back(flagl);

       	}else{

       		for (int j = xadj[i]; j < xadj[i+1]; j++){
       			if(part[adjncy[j]] == 1){
       				adjncyr.push_back(adjncy[j]);
       				flagr++;
       			}
       		}
       		xadjr.push_back(flagr);

       	}
    }

    //cout << " step 4 complete ! " <<endl;

    //Mapping for vertices id of left and right child
    int ** mappingLC = new int*[nVerticesLC];
    for(int i = 0; i < nVerticesLC; i++)
    	mappingLC[i] = new int[2];

    int ** mappingRC = new int*[nVerticesRC];
    for(int i = 0; i < nVerticesRC; i++)
       	mappingRC[i] = new int[2];

    for(int i = 0; i < nVerticesLC; i++){ //initialization for LC
       	mappingLC[i][0] = -1;
       	mappingLC[i][1] = i;

    }
    for(int i = 0; i < nVerticesRC; i++){ //initialization for RC
     	mappingRC[i][0] = -1;
       	mappingRC[i][1] = i;

    }
    for(int i = 0, j = 0; i < nVertices; i++){ //building mapping
       	if (part[i] == 0){
       		mappingLC[j][0] =i;
       		j++;
       	}

    }
    for(int i = 0, j = 0; i < nVertices; i++){ //building mapping
       	if (part[i] == 1){
       		mappingRC[j][0] =i;
       		j++;
       	}

    }


    int * trickMap  =  new int [nVertices];

    for ( int i=0; i<nVerticesLC ; i++){
    	trickMap[mappingLC[i][0]] = mappingLC[i][1];
    }
    for ( int i=0; i<nVerticesRC ; i++){
    	trickMap[mappingRC[i][0]] = mappingRC[i][1];
    }



    //cout << " step 5 complete !" <<endl;


    idx_t *adjAfterMLC = new int [adjncyl.size()];
    idx_t *adjAfterMRC = new int [adjncyr.size()];

    //cout << " step 5.1 complete !" <<endl;



    for (int i = 0; i < adjncyl.size(); i++){
//    	jlc = 0;
//    	while (jlc < nVerticesLC && flaglc == 0){
//       		if( mappingLC[jlc][0] == adjncyl[i] ){
//       			//cout << i <<endl;
//       			adjAfterMLC[i] = mappingLC[jlc][1];
//       			flaglc = 1;
//       		}
//       		jlc ++ ;
//    	}
//    	flaglc = 0;
    	adjAfterMLC[i] = trickMap[adjncyl[i]];
    }
    //cout << " step 5.2 complete !" <<endl;



    for (int i = 0; i < adjncyr.size(); i++){
//    	jrc = 0;
//    	while (jrc < nVerticesRC && flagrc == 0){
//       		if( mappingRC[jrc][0] == adjncyr[i] ){
//       			//cout << i << endl;
//       			adjAfterMRC[i] = mappingRC[jrc][1];
//       			flagrc = 1;
//       		}
//       		jrc ++ ;
//       	}
//    	flagrc = 0;
    	adjAfterMRC[i] = trickMap[adjncyr[i]];
    }

    //cout << " step 6 complete ! " <<endl;

    this->NumberOfVL = nVerticesLC;
    this->NumberOfVR = nVerticesRC;

    this->adjl = new int [adjncyl.size()];
       for(int i=0; i<adjncyl.size();i++){
    	   this->adjl[i] = adjncyl[i];
    }
    this->adjr = new int [adjncyr.size()];
       for(int i=0; i<adjncyr.size();i++){
    	   this->adjr[i] = adjncyr[i];
    }

    this->xal = new int [xadjl.size()];
       for(int i=0; i<xadjl.size();i++){
    	   this->xal[i] = xadjl[i];
    }

    this->xar = new int [xadjr.size()];
       for(int i=0; i<xadjr.size();i++){
    	   this->xar[i] = xadjr[i];
    }



    this->adjAftMLC = new int [adjncyl.size()];
    for(int i = 0; i < adjncyl.size(); i++){

        	this->adjAftMLC[i] = adjAfterMLC[i];

    }
    this->adjAftMRC = new int [adjncyr.size()];
    for(int i = 0; i < adjncyr.size(); i++){

        	this->adjAftMRC[i] = adjAfterMRC[i];

    }

    //cout << " step 8 complete ! " <<endl;

	delete mappingLC;
	delete mappingRC;
	delete adjAfterMLC;
	delete adjAfterMRC;
	delete trickMap;
	mappingLC = NULL;
	mappingRC = NULL;
	adjAfterMLC = NULL;
	adjAfterMRC = NULL;
	trickMap = NULL;
}

MetisNode::MetisNode(idx_t nVertices , idx_t xadj[], idx_t adjncy[] , int **mapOrg ){

    idx_t nWeights  = 1;
    idx_t nParts    = 2;

    idx_t objval; //number of edge separators
    idx_t part[nVertices];


    int ret = METIS_PartGraphRecursive(&nVertices,& nWeights, xadj, adjncy,
     				       NULL, NULL, NULL, &nParts, NULL,
     				       NULL, NULL, &objval, part);

/*
    int ret = METIS_PartGraphKway(&nVertices,& nWeights, xadj, adjncy,
				       NULL, NULL, NULL, &nParts, NULL,
				       NULL, NULL, &objval, part);*/



    idx_t nVerticesLC = NULL;
    idx_t nVerticesRC = NULL;
    idx_t nEdgesLC = NULL;
    idx_t nEdgesRC = NULL;

    /*
    for(unsigned part_i = 0; part_i < nVertices; part_i++){
    	std::cout << part_i << " " << part[part_i] << std::endl;
    }*/

    //int edges [objval*2];

    vector < vector<int> > tempv;
    vector < vector<int> > tempVCP;
    vector <string> EdgeSPR;

    for( int i = 0 ; i < nVertices; i++){
       	if (part[i] == 0){
       		nVerticesLC ++;

       		for(int j = xadj[i]; j < xadj[i+1]; j++){

       			if(part[adjncy[j]] == 1){

       				vector <int> tempVid;
       				tempVid.push_back(i);
       				tempVid.push_back(adjncy[j]);

       				//string temp = IntToString(i)+","+ IntToString(adjncy[j]);
       				//cout << temp <<endl;
       				//tempv.push_back(temp);
       				tempv.push_back(tempVid);

       				//cout<<"{" << i <<" , "<< adjncy[j] <<"}" << endl;
       			}
       		}
       	}
    }
    tempVCP = tempv;

    nVerticesRC = nVertices - nVerticesLC;

    for(int i = 0; i < tempv.size() ; i++){

    	for(int k = 0; k<nVertices; k++){
    		if(tempv[i][0]== mapOrg[k][1]){
   				tempVCP[i][0]= mapOrg[k][0];
   			}

   		}
    	for(int k = 0; k<nVertices; k++){
    		if(tempv[i][1]== mapOrg[k][1]){
    			//cout <<"see?  "<< tempv[i][1] <<" " << mapOrg[k][1]<< "  " << mapOrg[k][0] <<endl;
   				tempVCP[i][1]= mapOrg[k][0];
   			}

   		}

    }

    for( int i=0; i<tempVCP.size();i++){
    	string temp = IntToString(tempVCP[i][0])+","+ IntToString(tempVCP[i][1]);
    	EdgeSPR.push_back(temp);
    }
    this->EdgeSeparator = EdgeSPR;


    // building mapping table of original vertex id for both children
    this->mapOrgLC = new int*[nVerticesLC];
    int ** tempMOLC = new int*[nVerticesLC];
    for(int i = 0; i < nVerticesLC; i++){
       this->mapOrgLC[i] = new int[2];
       tempMOLC[i] = new int[2];
    }

    this->mapOrgRC = new int*[nVerticesRC];
    int ** tempMORC = new int*[nVerticesRC];
    for(int i = 0; i < nVerticesRC; i++){
       this->mapOrgRC[i] = new int[2];
       tempMORC[i] = new int[2];
    }


    for( int i = 0 , flagL = 0 , flagR = 0 ; i < nVertices; i++){
       	if (part[i] == 0){
       		tempMOLC[flagL][0] = i;
       		this->mapOrgLC[flagL][1] = flagL;
       		tempMOLC[flagL][1] = flagL;
       		flagL++;

       	}else{
       		tempMORC[flagR][0] = i;
       		this->mapOrgRC[flagR][1] = flagR;
       		tempMORC[flagR][1] = flagR;
       		flagR++;
       	}
    }



    for( int i=0; i< nVerticesLC ; i++){
    	for( int j=0; j< nVertices ; j++){
        	if(tempMOLC[i][0] == mapOrg[j][1]){
        		this->mapOrgLC[i][0]= mapOrg[j][0];
        	}
    	}

    }
    for( int i=0; i< nVerticesRC ; i++){
    	for( int j=0; j< nVertices ; j++){
        	if(tempMORC[i][0] == mapOrg[j][1]){
        		this->mapOrgRC[i][0]= mapOrg[j][0];
        	}
    	}

    }





    vector <int> adjncyl;
    vector <int> xadjl;
    vector <int> adjncyr;
    vector <int> xadjr;

    xadjl.push_back(0);
    xadjr.push_back(0);
       //cout << nVerticesLC<< "," << nVerticesRC<< endl;


    int flagl=0;
    int flagr=0;
    for(int i = 0; i < nVertices; i++){   		//Generate adjncy and xadj for both children
    	if(part[i]==0){

    		for (int j = xadj[i]; j < xadj[i+1]; j++){
       			if(part[adjncy[j]] == 0){
       				adjncyl.push_back(adjncy[j]);
       				flagl++;
       			}
       		}
       		xadjl.push_back(flagl);

       	}else{

       		for (int j = xadj[i]; j < xadj[i+1]; j++){
       			if(part[adjncy[j]] == 1){
       				adjncyr.push_back(adjncy[j]);
       				flagr++;
       			}
       		}
       		xadjr.push_back(flagr);

       	}
    }

    //Mapping for vertices id of left and right child
    int ** mappingLC = new int*[nVerticesLC];
    for(int i = 0; i < nVerticesLC; i++)
    	mappingLC[i] = new int[2];

    int ** mappingRC = new int*[nVerticesRC];
    for(int i = 0; i < nVerticesRC; i++)
       	mappingRC[i] = new int[2];

    for(int i = 0; i < nVerticesLC; i++){ //initialization for LC
       	mappingLC[i][0] = -1;
       	mappingLC[i][1] = i;

    }
    for(int i = 0; i < nVerticesRC; i++){ //initialization for RC
     	mappingRC[i][0] = -1;
       	mappingRC[i][1] = i;

    }
    for(int i = 0, j = 0; i < nVertices; i++){ //building mapping
       	if (part[i] == 0){
       		mappingLC[j][0] =i;
       		j++;
       	}

    }
    for(int i = 0, j = 0; i < nVertices; i++){ //building mapping
       	if (part[i] == 1){
       		mappingRC[j][0] =i;
       		j++;
       	}

    }


    int * trickMap  =  new int [nVertices];

    for ( int i=0; i<nVerticesLC ; i++){
    	trickMap[mappingLC[i][0]] = mappingLC[i][1];
    }
    for ( int i=0; i<nVerticesRC ; i++){
    	trickMap[mappingRC[i][0]] = mappingRC[i][1];
    }



    idx_t *adjAfterMLC = new int [adjncyl.size()];
    idx_t *adjAfterMRC = new int [adjncyr.size()];




    for (int i = 0; i < adjncyl.size(); i++){

    	adjAfterMLC[i] = trickMap[adjncyl[i]];
    }
    //cout << " step 5.2 complete !" <<endl;


    for (int i = 0; i < adjncyr.size(); i++){

    	adjAfterMRC[i] = trickMap[adjncyr[i]];
    }
    //cout << " step 5.2 complete !" <<endl;




    this->NumberOfVL = nVerticesLC;
    this->NumberOfVR = nVerticesRC;

    this->adjl = new int [adjncyl.size()];
       for(int i=0; i<adjncyl.size();i++){
    	   this->adjl[i] = adjncyl[i];
    }
    this->adjr = new int [adjncyr.size()];
       for(int i=0; i<adjncyr.size();i++){
    	   this->adjr[i] = adjncyr[i];
    }

    this->xal = new int [xadjl.size()];
       for(int i=0; i<xadjl.size();i++){
    	   this->xal[i] = xadjl[i];
    }

    this->xar = new int [xadjr.size()];
       for(int i=0; i<xadjr.size();i++){
    	   this->xar[i] = xadjr[i];
    }



//    this->mapLC = new int*[nVerticesLC];
//    for(int i = 0; i < nVerticesLC; i++){
//    	   this->mapLC[i] = new int[2];
//    }
//    for(int i = 0; i < nVerticesLC; i++){ //initialization for LC
//       	this->mapLC[i][0] = -1;
//       	this->mapLC[i][1] = i;
//
//    }
//
//    for(int i = 0; i < nVerticesLC; i++){
//
//	   for(int j = 0; j<2 ; j++){
//		   this->mapLC[i][j] = mappingLC[i][j];
//	   }
//    }
//
//
//
//    this->mapRC = new int*[nVerticesRC];
//    for(int i = 0; i < nVerticesRC; i++){
//       this->mapRC[i] = new int[2];
//    }
//    for(int i = 0; i < nVerticesRC; i++){ //initialization for LC
//       	this->mapRC[i][0] = -1;
//       	this->mapRC[i][1] = i;
//
//    }
//
//    for(int i = 0; i < nVerticesRC; i++){
//    	for(int j = 0; j<2 ; j++){
//        	this->mapRC[i][j] = mappingRC[i][j];
//    	}
//
//    }

    this->adjAftMLC = new int [adjncyl.size()];
    for(int i = 0; i < adjncyl.size(); i++){

        	this->adjAftMLC[i] = adjAfterMLC[i];

    }
    this->adjAftMRC = new int [adjncyr.size()];
    for(int i = 0; i < adjncyr.size(); i++){

        	this->adjAftMRC[i] = adjAfterMRC[i];

    }


	delete mappingLC;
	delete mappingRC;
	delete adjAfterMLC;
	delete adjAfterMRC;
	delete tempMOLC;
	delete tempMORC;
	mappingLC = NULL;
	mappingRC = NULL;
	adjAfterMLC = NULL;
	adjAfterMRC = NULL;
	tempMOLC = NULL;
	tempMORC = NULL;

}

void MetisNode::printItself(){
	cout<<"\n" <<"print: "<<endl;
    for (int i =0; i < this->EdgeSeparator.size();i++){
    	cout << this->EdgeSeparator[i] <<" ";
    }
    cout <<endl;
}

void MetisNode::printMapping(){

//	cout<<"\n\n" << "print mapping for left child: " <<endl;
//	for (int i = 0 ; i<this->NumberOfVL; i++){
//		cout << this->mapLC [i][0] << "," << this->mapLC [i][1] << "    ";
//	}
	cout<<"\n" << "print  ORGNIAL mapping for left child: " <<endl;
	for (int i = 0 ; i<this->NumberOfVL; i++){
		cout << this->mapOrgLC [i][0] << "," << this->mapOrgLC [i][1] << "    ";
	}



//	cout<<"\n\n" << "print mapping for right child: " <<endl;
//	for (int i = 0 ; i<this->NumberOfVR; i++){
//		cout << this->mapRC [i][0] << ","  << this->mapRC [i][1] << "    ";
//	}
	cout<<"\n" << "print  ORGNIAL mapping for right child: " <<endl;
	for (int i = 0 ; i<this->NumberOfVR; i++){
		cout << this->mapOrgRC [i][0] << "," << this->mapOrgRC [i][1] << "    ";
	}

	cout <<endl;
}
void MetisNode::printXadj(){
	cout<< "\nnVertex --- left : "<< this->NumberOfVL << "  right : " << this->NumberOfVR ;
	cout<<"\n" << "print Xadj for left child: " <<endl;
	for (int i = 0 ; i<this->NumberOfVL+1; i++){
		cout << this->xal[i] << "," ;
	}

	cout<<"\n" << "print Xadj for right child: " <<endl;
	for (int i = 0 ; i<this->NumberOfVR+1; i++){
		cout << this->xar [i] << "," ;
	}
	cout <<endl;
}
void MetisNode::printAdj(){
//	cout<<"\n" << "print adj for left child: " <<endl;
//	for (int i = 0 ; i<(this->xal[this->NumberOfVL]); i++){
//		cout << this->adjl[i] << endl ;
//	}

	cout<<"\n" << "print adj for right child: " <<endl;
	for (int i = 0 ; i<(this->xar[this->NumberOfVR]); i++){
		cout << this->adjr[i] << endl; ;
	}


//
//	cout<<"\n" << "print adj after mapping for left child: " <<endl;
//	for (int i = 0 ; i<(this->xal[this->NumberOfVL]); i++){
//		cout << this->adjAftMLC[i] << endl ;
//	}

	cout<<"\n" << "print adj after mapping for right child: " <<endl;
	for (int i = 0 ; i<(this->xar[this->NumberOfVR]); i++){
		cout << this->adjAftMRC[i] << endl; ;
	}
	cout <<endl;
}

#endif /* METISNODE_H_ */
