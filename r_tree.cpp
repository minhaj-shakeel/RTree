#include<iostream>
#include "file_manager.h"
#include "errors.h"
#include "constants.h"
#include<cstring>
#include <fstream>
#include <float.h>
#include <limits.h>

using namespace std;

FileManager fm;
# define INT_SIZE  sizeof(int)

struct Node{
    int ID;   // 1 int
    int *mbr; // 2*d int
    bool LeafNode;  // 1 int
    int parentID; // 1 int
    int childCnt; // 1 int
    int *childID; // maxCap+1 int                
    int *childData; // 2*d*(maxCap+1) int    Giving Extra Space for One Hyper Rectangle 
    
    Node(int d, int maxCap){
        childCnt=0;
        parentID=-1;
        mbr = (int*)malloc(2*d*INT_SIZE);
        childID = (int*)malloc((maxCap+1)*INT_SIZE);
        childData = (int*)malloc(2*d*(maxCap+1)*INT_SIZE);
    }
    int getChildID(int i);
    int* getChildData(int i, int d);
    int selectChild(int *point,int d, int maxCap);
    bool isLeaf();
    void CalculateMBR(int d, int maxCap);
    bool addChildData(int *Point, int cid, int d);
};

class Rtree{
        int cnt, height;
    public:
        int d, maxCap, root;
        string treeFile;
        Rtree(int dim,int maxChild){
            root=-1;
            height=-1;
            cnt=0;
            d=dim;
            maxCap=maxChild;
            treeFile = "treeData.txt";
        }

        int GetNewId();
        int IncreaseCnt();
        int BulkLoad(const char* fileName,int N);
        int AssignParent(int start, int end);
        bool isPresent(int *point, int nodeID);
};


void PrintNode(Node *newNode, int d, int maxCap);
bool isEqual(int *point1, int *point2, int d);
bool isInside(int *point, int *mbr, int d);

int main(int argc, char* argv[]){
    string queryFile = argv[1];
    int maxCap = atoi(argv[2]);
    int d = atoi(argv[3]);
    string outputFile = argv[4];

    Rtree rtree(d,maxCap);
    
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
            
            char fileName[fName.length()+1];
            strcpy(fileName,fName.c_str());

            rootIndex = rtree.BulkLoad(fileName, N);
            
            outFile << cmd;
            outFile<<"\n\n\n";
        } else if (cmd=="INSERT"){
            int PointArray[d];
            for (int i=0; i<d ;i++){
                inFile>>PointArray[i];
            }
    
            outFile << cmd;
            outFile<<"\n\n\n";
        } else if (cmd=="QUERY"){
        
            int PointArray[d];
            for (int i=0; i<d ;i++){
                inFile>>PointArray[i];
            }

            bool result = rtree.isPresent(PointArray,rootIndex);
            result==true?outFile << "TRUE":outFile << "FALSE";
            outFile<<"\n\n\n";
        } else{
            cout << "INVALID COMMAND";
        }
        
    }
    inFile.close();
    outFile.close();
    fm.DestroyFile("treeData.txt");
}

int Rtree::BulkLoad(const char* fileName,int N){
    FileHandler fhIn, fhOut;
    PageHandler phIn, phOut;
    char *dataIn, *dataOut;
    fhIn = fm.OpenFile(fileName);
    fhOut = fm.CreateFile("treeData.txt");
    int pointSize = INT_SIZE*d;
    int pointPerPage = PAGE_CONTENT_SIZE/pointSize;
    int NumNodes = N/maxCap+((N%maxCap)!=0);
    
    
    int resVal[d];
    for(int i=0; i<d; i++){
        resVal[i]=INT_MIN;
    }
    //Read the N points one by one
    for (int i=0; i<NumNodes; i++){
        // Create A New Node, make it leaf Node, and then Start Adding Points
        Node *newNode = new Node(d,maxCap);
        newNode->ID = GetNewId();
        IncreaseCnt();
        newNode->LeafNode=true;
        
        // Initialize Space for calculating MBR in running hand
        int minPoint[d];
        int maxPoint[d];
        for(int j=0; j<d ;j++){
            minPoint[j]=INT_MAX;
            maxPoint[j]=INT_MIN;
        }
        
        // reset all the child IDs to -1
        for(int j=0; j<maxCap+1; j++){
            *(newNode->childID+j)=-1;
        }

        int remaining = N-i*maxCap;
        for(int j=0; j<min(maxCap,remaining); j++){  
            
            //Calculate the page Number and Offset
            int pageNum =  ((i*maxCap+j)*pointSize)/PAGE_CONTENT_SIZE;
            int offset =   ((i*maxCap+j)%pointPerPage)*d;
            
            phIn = fhIn.PageAt(pageNum);
            dataIn = phIn.GetData();
            int* pointStart = (int *)(dataIn)+offset;
            *(newNode->childID+j)=0;
            newNode->childCnt++;
            int *storeStart = newNode->getChildData(j,d);
            memcpy(storeStart,resVal,pointSize);
            memcpy(storeStart+d,pointStart,pointSize);

            for(int k=0; k<d ;k++){
                minPoint[k]=min(minPoint[k],*(pointStart+k));
                maxPoint[k]=max(maxPoint[k],*(pointStart+k));
            }
            fhIn.UnpinPage(pageNum);
        }   

        // Store the calculated MBR
        memcpy(newNode->mbr,minPoint,pointSize);
        memcpy(newNode->mbr+d,maxPoint,pointSize);

        // Store the parent node into the New Page 
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
        return AssignParent(0,NumNodes);
}




bool Node::addChildData(int *Point, int cid, int d){
    int *position = getChildData(childCnt, d);
    childID[childCnt]=cid;
    if (LeafNode){
        memcpy(position+d,Point,d*INT_SIZE);
        childCnt++;
        return true;
    }
    else{
        memcpy(position, Point, 2*d*INT_SIZE);
        childCnt++;
        return true;
    }
    return true;
}


bool isInside(int *point, int *mbr, int d){
    for(int i=0; i<d; i++){
        int minD = *(mbr+i);
        int maxD = *(mbr+(i+d));
        int pointD = *(point+i);
        //cout << minD << " " << pointD 
        if ((pointD < minD) | (pointD > maxD)){
            return false;
        } 
    }
    return true;
}




bool Node::isLeaf(){
   return LeafNode;
}



void TestFunction(string fileName, int d, int maxCap, int N){
    FileHandler fh = fm.OpenFile("treeData.txt");
    PageHandler ph;
    for (int i=0; i<111; i++){
        ph = fh.PageAt(i);
        char *data = ph.GetData();

        Node *newNode = new Node(d,maxCap);
        memcpy(newNode, data, sizeof(Node));
        
        for(int k=0; k<maxCap; k++){
            for(int j=0 ;j< 2*d; j++){
                cout << *(int *)(newNode->childData+2*d*k+j) << " ";
            }
            cout << endl;
        }
        cout << endl;
        
        for (int j = 0; j < 2*d; j++){
            cout << *(int *)(newNode->mbr+j) << " ";
        }
        cout << endl;
        fh.UnpinPage(i);
    }
}




bool Rtree::isPresent(int *point, int nodeID){
    FileHandler fh = fm.OpenFile("treeData.txt");
    PageHandler ph = fh.PageAt(nodeID);
    char* data = ph.GetData();
    struct Node *nodeData = new Node(d, maxCap);
    memcpy(nodeData,data,sizeof(Node));
    fh.UnpinPage(nodeID);
    // If the Node is LeafNode
    if (nodeData->isLeaf()){
        for (int i=0; i<maxCap; i++){
            int childID = nodeData->getChildID(i);
            if (childID==-1)
                continue;
            int *pointData = nodeData->getChildData(i,d)+d;
            if (isEqual(point, pointData, d)){
                fm.CloseFile(fh);
                return true;
            }
        }
        fm.CloseFile(fh);
        return false;
    } 
    
    // If the Node is non LeafNode
    else{
        bool ans = false;
        for(int i=0; i<maxCap; i++){
            int childID = nodeData->getChildID(i);
            if (childID <0)
                ans = ans|false;
            else{
                int *childMBR = nodeData->getChildData(i, d);
                if (!isInside(point, childMBR, d)){
                    ans = ans|false;
                }
                else{
                    fm.CloseFile(fh);
                    ans = ans|isPresent(point, childID);
                }
            }
        }
        fm.CloseFile(fh);
        return ans; 
    }
}

int Rtree::AssignParent(int start, int end){
    if (end-start==1){
        root = start;
        return root;
    }
    
    FileHandler fhIn, fhOut;
    PageHandler phIn, phOut;
    char *dataIn, *dataOut;
    fhIn = fm.OpenFile("treeData.txt");
    fhOut = fm.OpenFile("treeData.txt");
    
    int childNodes = end-start;
    int numParents = childNodes/maxCap+((childNodes%maxCap)!=0);

    int resVal[d];
    for(int i=0; i<d; i++){
        resVal[i]=INT_MIN;
    }

    for(int i=0 ;i<numParents; i++){
        Node *parentNode = new Node(d,maxCap);
        parentNode->ID = GetNewId();
        IncreaseCnt();
        parentNode->LeafNode = false;
        
        int minPoint[d];
        int maxPoint[d];
        for(int j=0; j<d ;j++){
            minPoint[j]=INT_MAX;
            maxPoint[j]=INT_MIN;
        }

        // Store Child's Data into Parent Node
        int remaining = childNodes-i*maxCap;
        for(int j=0; j<min(maxCap,remaining); j++ ){
            int childID = start+i*maxCap+j;
            phIn = fhIn.PageAt(childID);
            dataIn = phIn.GetData();

            // Read The Child Node
            Node *childNode = new Node(d,maxCap);
            memcpy(childNode,dataIn,sizeof(Node));
            childNode->parentID = parentNode->ID;
            memcpy(dataIn,childNode,sizeof(Node));
            fhIn.MarkDirty(childID);
            fhIn.UnpinPage(childID);
            fhIn.FlushPage(childID);

            for(int k=0; k<d ;k++){
                minPoint[k]=min(minPoint[k],*(childNode->mbr+k));
                maxPoint[k]=max(maxPoint[k],*(childNode->mbr+k+d));
            }
            
            //Copy Child Data and Increase the Child Count
            *(int *)(parentNode->childID+j)=childNode->ID;
            memcpy(parentNode->childData+2*d*j,childNode->mbr,2*d*INT_SIZE); 
            parentNode->childCnt++;
        }

        //set all the child ID's to -1 for the remainign childs
        for(int j=remaining; j<maxCap; j++){
           *(parentNode->childID+j) = -1;
           int *cdata = parentNode->getChildData(j,d);
           memcpy(cdata,resVal,d*INT_SIZE);

        }
    
        //Store the Parent's mbr
        memcpy(parentNode->mbr,minPoint,d*INT_SIZE);
        memcpy(parentNode->mbr+d,maxPoint,d*INT_SIZE);

        //Store the Parent's Data
        phOut = fhOut.NewPage();
        int parentPageNum = phOut.GetPageNum();
        dataOut = phOut.GetData();
        memcpy(dataOut,parentNode,sizeof(Node));
        fhOut.MarkDirty(parentPageNum);
        fhOut.UnpinPage(parentPageNum);
    }
    fhOut.FlushPages();
    fm.CloseFile(fhIn);
    fm.CloseFile(fhOut);
    return AssignParent(end, end+numParents);
}

int Rtree::GetNewId(){
    return cnt;
}

int Rtree::IncreaseCnt(){
    cnt++;
    return cnt;
}

int Node::getChildID(int i){
    return *(childID+i);
}

int* Node::getChildData(int i, int d){
    return childData+i*2*d;
}

void PrintNode(Node *newNode, int d, int maxCap){
    
    cout << "Node id " << newNode->ID << endl;
    for (int j = 0; j < 2*d; j++){
            cout << *(int *)(newNode->mbr+j) << " ";
    }
    cout << endl;
    cout << endl;
    for(int k=0; k<maxCap+1; k++){
            cout << (newNode->childID)[k] << " " ;
            for(int j=0 ;j< 2*d; j++){
                cout << *(int *)(newNode->childData+2*d*k+j) << " ";
            }
            cout << endl;
        }
        cout << endl;

}

bool isEqual(int *point1, int *point2, int d){
    for(int i=0; i<d; i++){
        if (*(point1+i)!= *(point2+i))
            return false;
    }
    return true;
}
