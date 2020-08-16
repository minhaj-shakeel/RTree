#include<iostream>
#include "file_manager.h"
#include "errors.h"
#include "constants.h"
#include<cstring>
#include <fstream>
#include <limits.h>

using namespace std;

FileManager fm;

struct Node{
    int ID;   // 1 int
    char *mbr; // 2*d int
    int parentID; // 1 int
    char *childID; // maxCap int                
    char *childData; // 2*d*maxCap int
    
    Node(int d, int maxCap){
        mbr = (char*)malloc(2*d*INT_SIZE);
        childID = (char*)malloc(maxCap*INT_SIZE);
        childData = (char*)malloc(2*d*maxCap*INT_SIZE);
    }
};

void BulkLoad(string fileName,int N, int d, int maxCap);
void AssignParent(string fileName, int start, int end, int d, int maxCap);

void ReadNode(char *dataIn, Node *node);

void TestFunction(string fileName, int d, int maxCap){
    FileHandler fh = fm.OpenFile("treeData.txt");
    PageHandler ph;
    for (int i=0; i< 111; i++){
        ph = fh.PageAt(i);
        char *data = ph.GetData();

        Node *newNode = new Node(d,maxCap);
        memcpy(newNode, data, sizeof(Node));
        cout << newNode->ID  << " " << newNode->parentID << endl;
        for (int j = 0; j < maxCap; j++){
            cout << *(int *)(data+j*INT_SIZE) << " ";
        }
        cout << endl;
        fh.UnpinPage(i);
    }
}


int main(int argc, char* argv[]){
    string queryFile = argv[1];
    int maxCap = atoi(argv[2]);
    int d = atoi(argv[3]);
    string outputFile = argv[4];

    fstream inFile,outFile;
    inFile.open(queryFile,ios::in);
    outFile.open(outputFile,ios::out);


    string cmd;
    while(inFile >> cmd){
        if (cmd=="BULKLOAD"){
            int N;
            string fileName;
            inFile >> fileName;
            inFile >> N; 
            //BulkLoad(fileName, N, d, maxCap);
            TestFunction(fileName, d, maxCap);
            outFile << cmd;
        } else if (cmd=="INSERT"){
            int PointArray[d];
            for (int i=0; i<d ;i++){
                inFile>>PointArray[i];
            }
            outFile << cmd;
        } else if (cmd=="QUERY"){
            int PointArray[d];
            for (int i=0; i<d ;i++){
                inFile>>PointArray[i];
            }
            outFile << cmd;
        } else{
            cout << "INVALID COMMAND";
        }
        outFile<<"\n\n\n";
    }
    inFile.close();
    outFile.close();


}

void BulkLoad(string fileName, int N, int d, int maxCap){
    FileHandler fhIn, fhOut;
    PageHandler phIn, phOut;
    char *dataIn, *dataOut;
    fhIn = fm.OpenFile("Testcases/TC_1/sortedData_2_10_100.txt");
    fhOut = fm.CreateFile("treeData.txt");
    
    int pointSize = INT_SIZE*d;
    int pointPerPage = PAGE_CONTENT_SIZE/pointSize;
    int resVal = INT_MIN;
    //Read the N points one by one
    for (int i=0; i<N; i++){
        //Calculate the page Number and Offset
        int pageNum =  (i*pointSize)/PAGE_CONTENT_SIZE;
        int offset = (i%pointPerPage)*pointSize;
        
        phIn = fhIn.PageAt(pageNum);
        dataIn = phIn.GetData();
        char* pointStart = dataIn+offset;

        Node* newNode = new Node(d,maxCap);
        newNode->ID = i;
        
        // Store residual Value in first d intgers
        for(int j =0; j<d; j++){
            memcpy(newNode->mbr+j*INT_SIZE,&resVal,INT_SIZE);
        }
        //store d points from the file into mbr
        memcpy(newNode->mbr+d*INT_SIZE,pointStart,d*INT_SIZE);
        //cout << *(int *)pointStart << endl;
        //set all the child ID's to 1
        for(int j=0; j<d; j++){
           *(newNode->childID+j*INT_SIZE) = -1;
        }

        //store residual Value in the child Data
        for(int j=0; j<2*d*maxCap; j++){
            memcpy(newNode->childData+j*INT_SIZE,&resVal,INT_SIZE);
        }
        
        //Store the every Node into the New Page
        phOut = fhOut.NewPage();
        dataOut = phOut.GetData();
        memcpy(dataOut,newNode,sizeof(Node));

        //Mark dirty and unpin the page
        int outPageNum = phOut.GetPageNum();
        fhOut.MarkDirty(outPageNum);
        fhOut.UnpinPage(outPageNum);
    }
    fhOut.FlushPages();
    fm.CloseFile(fhOut);
    fm.CloseFile(fhIn);  

    AssignParent(fileName, 0, N, d, maxCap);  
}

void AssignParent(string fileName, int start, int end, int d, int maxCap){
    cout << "AssignParent "<< start << " " << end << endl;
    FileHandler fhIn, fhOut;
    PageHandler phIn, phOut;
    char *dataIn, *dataOut;
    fhIn = fm.OpenFile("treeData.txt");
    fhOut = fm.OpenFile("treeData.txt");
    if (end-start==1)
        return;
    int childNodes = end-start;
    int numParents = childNodes/maxCap+((childNodes%maxCap)!=0);
    for(int i=0 ;i<numParents; i++){
        Node *parentNode = new Node(d,maxCap);
        parentNode->ID = end+i;
        //store Child Data into Parent Node
        for(int j=0; j<maxCap; j++){
            int childID = start+i*maxCap+j;
            phIn = fhIn.PageAt(childID);
            dataIn = phIn.GetData();

            //Read the ChildNode
            Node *childNode = new Node(d,maxCap);
            memcpy(childNode,dataIn,sizeof(Node));
            childNode->parentID = parentNode->ID;
            memcpy(dataIn,childNode,sizeof(Node));
            fhIn.MarkDirty(childID);
            fhIn.UnpinPage(childID);

            //Copy Child Data
            *(int *)(parentNode->childID+j)=childNode->ID;
            memcpy(parentNode->childData+2*d*j*INT_SIZE,childNode->mbr,2*d*INT_SIZE); 
            
        }

        //Store the Parent's Data
        phOut = fhOut.NewPage();
        int parentPageNum = phOut.GetPageNum();
        dataOut = phOut.GetData();
        memcpy(dataOut,parentNode,sizeof(Node));
        fhOut.MarkDirty(parentPageNum);
        fhOut.UnpinPage(parentPageNum);
 
    } 
    fhOut.FlushPages();
    AssignParent(fileName, end, end+numParents, d, maxCap);

}



void ReadNode(char *dataIn, Node *node,int d,int maxCap){
    node->ID = *(int *)dataIn;
    memcpy(node->mbr,dataIn+INT_SIZE,2*d*INT_SIZE);
    node->parentID =*(int *)(dataIn+(2*d+1)*INT_SIZE);
    memcpy(node->childID,dataIn+(2*d+1)*INT_SIZE,maxCap*INT_SIZE);
    memcpy(node->childData,dataIn+(2*d+1+maxCap)*INT_SIZE,2*d*maxCap*INT_SIZE);
}