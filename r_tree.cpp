#include<iostream>
#include "file_manager.h"
#include "errors.h"
#include "constants.h"
#include<cstring>
#include <fstream>
#include <limits.h>

using namespace std;

FileManager fm;
string tempFile = "treeData.txt";

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
    int getChildID(int i);
    char *getChildData(int i, int d);
    bool isLeaf(int maxCap);
};

int Node::getChildID(int i){
    return *(int *)(childID+i*INT_SIZE);
}

char* Node::getChildData(int i, int d){
    return childData+i*2*d*INT_SIZE;
}

bool isInside(char *point, char *mbr, int d){
    for(int i=0; i<d; i++){
        int minD = *(int *)(mbr+i*INT_SIZE);
        int maxD = *(int *)(mbr+(i+d)*INT_SIZE);
        int pointD = *(int *)(point+i*INT_SIZE);
        //cout << minD << " " << pointD 
        if ((pointD < minD) | (pointD > maxD)){
            return false;
        } 
    }
    return true;
}

bool isEqual(char *point1, char *point2, int d){
    for(int i=0; i<d; i++){
        if (*(int *)(point1+i*INT_SIZE)!= *(int *)(point2+i*INT_SIZE))
            return false;
    }
    return true;
}



bool Node::isLeaf(int maxCap){
    for (int i=0 ;i<maxCap ;i++){
        if (getChildID(i) >= 0)
            return false;
    }
    return true;
}

int BulkLoad(char* fileName,int N, int d, int maxCap);
int AssignParent(string fileName, int start, int end, int d, int maxCap);
bool isPresent(char *point, int nodeID, int d, int maxCap);
void TestFunction(string fileName, int d, int maxCap, int N){
    FileHandler fh = fm.OpenFile("treeData.txt");
    PageHandler ph;
    for (int i=N; i<N+10; i++){
        ph = fh.PageAt(i);
        char *data = ph.GetData();

        Node *newNode = new Node(d,maxCap);
        memcpy(newNode, data, sizeof(Node));
        
        for(int k=0; k<maxCap; k++){
            for(int j=0 ;j< 2*d; j++){
                cout << *(int *)(newNode->childData+(2*d*k+j)*INT_SIZE) << " ";
            }
            cout << endl;
        }
        cout << endl;
        
        for (int j = 0; j < 2*d; j++){
            cout << *(int *)(newNode->mbr+(j)*INT_SIZE) << " ";
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
    int rootIndex;
    while(inFile >> cmd){
        if (cmd=="BULKLOAD"){
            int N;
            string fName;
            inFile >> fName;
            inFile >> N; 
            
            char *fileName = (char*)malloc(fName.size());
            memcpy(fileName, &fName,fName.size());
            
            rootIndex = BulkLoad(fileName, N, d, maxCap);
            cout << "root= " << rootIndex << endl; 
            int array[2];
            array[0]=585167411;
            array[1]=2791691; 
            
            cout << isPresent((char *)array, rootIndex, d, maxCap);
            //TestFunction(fileName, d, maxCap, N);
            outFile << cmd;
        } else if (cmd=="INSERT"){
            int PointArray[d];
            for (int i=0; i<d ;i++){
                inFile>>PointArray[i];
            }
            bool result = isPresent((char *)PointArray,rootIndex, d, maxCap);
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
    fm.DestroyFile("treeData.txt");


}

int BulkLoad(char* fileName, int N, int d, int maxCap){
    FileHandler fhIn, fhOut;
    PageHandler phIn, phOut;
    char *dataIn, *dataOut;
    fhIn = fm.OpenFile("Testcases/TC_1/sortedData_2_10_100.txt");
    //fhOut = fm.CreateFile((char *)&tempFile);
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
        
        // Store residual Value in first d integers
        for(int j =0; j<d; j++){
            memcpy(newNode->mbr+j*INT_SIZE,&resVal,INT_SIZE);
        }
        //store d points from the file into mbr
        memcpy(newNode->mbr+d*INT_SIZE,pointStart,d*INT_SIZE);
        //cout << *(int *)pointStart << endl;
        //set all the child ID's to -1
        for(int j=0; j<maxCap; j++){
           *(int *)(newNode->childID+j*INT_SIZE) = -1;
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

    return AssignParent(fileName, 0, N, d, maxCap);  
}

int AssignParent(string fileName, int start, int end, int d, int maxCap){
    int resVal = INT_MIN;
    FileHandler fhIn, fhOut;
    PageHandler phIn, phOut;
    char *dataIn, *dataOut;
    fhIn = fm.OpenFile("treeData.txt");
    fhOut = fm.OpenFile("treeData.txt");
    if (end-start==1)
        return start;
    int childNodes = end-start;
    int numParents = childNodes/maxCap+((childNodes%maxCap)!=0);
    int o;
    start==0? o=d: o=0; 
    for(int i=0 ;i<numParents; i++){
        Node *parentNode = new Node(d,maxCap);
        parentNode->ID = end+i;
        
        int minPoint[d];
        int maxPoint[d];
        for(int j=0; j<d ;j++){
            minPoint[j]=INT_MAX;
            maxPoint[j]=INT_MIN;
        }
        
        //store Child Data into Parent Node
        int lastNodeChilds = maxCap;
        if (i==numParents-1){
            lastNodeChilds=childNodes-(numParents-1)*maxCap;
        }
        for(int j=0; j<min(maxCap,lastNodeChilds); j++){
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

            for(int k=0; k<d ;k++){
                minPoint[k]=min(minPoint[k],*(int *)(childNode->mbr+(k+o)*INT_SIZE));
                maxPoint[k]=max(maxPoint[k],*(int *)(childNode->mbr+(k+d)*INT_SIZE));
            }
            
            //Copy Child Data
            *(int *)(parentNode->childID+j*INT_SIZE)=childNode->ID;
            memcpy(parentNode->childData+2*d*j*INT_SIZE,childNode->mbr,2*d*INT_SIZE); 
            
        }
        
        //set all the child ID's to -1 for the remainign childs
        for(int j=lastNodeChilds; j<maxCap; j++){
           *(int *)(parentNode->childID+j*INT_SIZE) = -1;
        }

        //store residual Value in the child Data for the remaining childs
        for(int j=2*d*lastNodeChilds; j<2*d*maxCap; j++){
            memcpy(parentNode->childData+j*INT_SIZE,&resVal,INT_SIZE);
        }
        
        //Store the Parent's mbr
        memcpy(parentNode->mbr,minPoint,d*INT_SIZE);
        memcpy(parentNode->mbr+d*INT_SIZE,maxPoint,d*INT_SIZE);

        //Store the Parent's Data
        phOut = fhOut.NewPage();
        int parentPageNum = phOut.GetPageNum();
        dataOut = phOut.GetData();
        memcpy(dataOut,parentNode,sizeof(Node));
        fhOut.MarkDirty(parentPageNum);
        fhOut.UnpinPage(parentPageNum);
 
    } 
    fhOut.FlushPages();
    return AssignParent(fileName, end, end+numParents, d, maxCap);

}

bool isPresent(char *point, int nodeID, int d, int maxCap){
    FileHandler fh = fm.OpenFile("treeData.txt");
    PageHandler ph = fh.PageAt(nodeID);
    char* data = ph.GetData();
    struct Node *nodeData = new Node(d, maxCap);
    memcpy(nodeData,data,sizeof(Node));
    fh.UnpinPage(nodeID);
    // If the Node is LeafNode
    if (nodeData->isLeaf(maxCap)){
        char *pointData = (nodeData->mbr+d*INT_SIZE);
        return isEqual(point, pointData, d);
    } 
    
    // If the Node is non LeafNode
    else{
        bool ans = false;
        for(int i=0; i<maxCap; i++){
            int childID = nodeData->getChildID(i);
            if (childID <0)
                ans = ans|false;
            else{
                char *childMBR = nodeData->getChildData(i, d);
                if (!isInside(point, childMBR, d)){
                    ans = ans|false;
                }
                else{
                    ans = ans|isPresent(point, childID, d, maxCap);
                }
            }
        }
        return ans; 
    }
}
