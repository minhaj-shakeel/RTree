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
string tempFile = "treeData.txt";

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
    int *getChildData(int i, int d);
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
        int insert(int *point, int rn);
        int SplitNode(int nodeId);
        int updateMBR(int nodeId);
        int updateMBR(int nodeId, int childId);
};

int Rtree::updateMBR(int nodeId){
    FileHandler fh = fm.OpenFile("treeData.txt");
    PageHandler ph = fh.PageAt(nodeId);
    Node *node = new Node(d, maxCap);
    char *data = ph.GetData();
    memcpy(node, data, sizeof(Node));
    node->CalculateMBR(d, maxCap);
    fh.MarkDirty(nodeId);
    fh.UnpinPage(nodeId);
    fh.FlushPage(nodeId);

    // Make Changes to The Parent
    int parentId = node->parentID;
    
    // If There is no Valid ParentId, then
    // we are at root, return as it is.
    if (parentId < 0 )
        return nodeId;

    cout << "mbr test" << endl;
    // Else, read the Parent Node, update the child's MBR 
    // and call updateMBR recursively.
    ph = fh.PageAt(parentId);
    Node *parent = new Node(d, maxCap);
    char *parentData = ph.GetData();
    memcpy(parent,parentData, sizeof(Node));

    for (int i=0; i<maxCap; i++){
        if (parent->getChildID(i) == nodeId){
            int *newPos = parent->getChildData(i, d);
            memcpy(newPos, node->mbr, 2*d*INT_MAX);
        }
    }

    fh.MarkDirty(parentId);
    fh.UnpinPage(parentId);
    fh.FlushPage(parentId);
    fm.CloseFile(fh);
    return updateMBR(parentId);



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
double distance(int *point1, int *point2, int d){
    double dist = 0;
    for(int i=0; i<d; i++)
        dist+=(point2[i]-point1[i])*(point2[i]-point1[i]);
    return dist;
}


double Wastage(int *point1, int *point2, int d){
    double wastage = 0;
    double area1, area2, full;
    area1 = area2 = full = 1;
    for(int i =0; i< d; i++){
        area1*=(point1[d+i]-point1[i]);
        area1*=(point2[d+i]-point2[i]);
        int maxVal = max(point1[d+i],point2[d+i]);
        int minVal = min(point1[i],point2[i]);
        full*=(maxVal-minVal);
    }
    return full-(area1+area2);
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


int Rtree::SplitNode(int nodeId){
    FileHandler fh = fm.OpenFile("treeData.txt");
    PageHandler ph = fh.PageAt(nodeId);
    char *data = ph.GetData();
    
    // Read the overflowed Node
    Node *currNode = new Node(d, maxCap);
    memcpy(currNode, data, sizeof(Node));

    // Create 2 New Nodes, Assign the id of
    // overflowed node to one Node, and fresh
    // new Id to the other one, make both
    // the leaf nodes.
    Node *newNode1 = new Node(d, maxCap);
    Node *newNode2 = new Node(d, maxCap);
    newNode1->ID = currNode->ID;
    newNode2->ID = GetNewId();
    IncreaseCnt();
   
    newNode1->LeafNode = currNode->isLeaf();
    newNode2->LeafNode = currNode->isLeaf();
    for(int i=0; i<maxCap; i++){
        newNode1->childID[i]=-1;
        newNode2->childID[i]=-1;
    }
    
    
    // find Index e1, e2 using Quadratic Split
    int e1, e2;
    if (currNode->isLeaf()) {
        double maxDist = DBL_MIN;
        for(int i=0; i<maxCap+1; i++){
            int *point1 = currNode->getChildData(i,d);
            for(int j=i+1; j<maxCap+1; j++){
                int *point2 = currNode->getChildData(j,d);
                double currDist = distance(point1+d,point2+d,d);
                if (currDist>maxDist){
                    e1=i;
                    e2=j;
                    maxDist=currDist;
                }
            }
        }
    }
    
    else{
        double maxWastage = DBL_MIN;
        for(int i=0; i<maxCap+1; i++){
            int *point1 = currNode->getChildData(i,d);
            for(int j=i+1; j<maxCap+1; j++){
                int *point2 = currNode->getChildData(j,d);
                double wastage = Wastage(point1,point2,d);
                if (wastage>maxWastage)
                e1=i; 
                e2=j;
                maxWastage = wastage;
            }
        }
    }


    int *point1 = currNode->getChildData(e1,d);
    int *point2 = currNode->getChildData(e2,d);

    int *point1Str = newNode1->getChildData(0,d);
    int *point2Str = newNode2->getChildData(0,d);


    int m = maxCap/2+maxCap%2;
    // Assign these Two Seed Points to Diff Nodes
    if (currNode->isLeaf()){
    newNode1->addChildData(point1+d, currNode->getChildID(e1),d);
    newNode2->addChildData(point2+d, currNode->getChildID(e2),d);
    }
    else {
    newNode1->addChildData(point1, currNode->getChildID(e1),d);
    newNode2->addChildData(point2, currNode->getChildID(e2),d);
    }
    
    // Assign the remaining Points between Nodes
    for( int i =0; i <maxCap+1; i++){
        int *currPoint = currNode->getChildData(i,d);
        // If index is one of seed child's
        if (i == e1 | i== e2)
            continue;

        // Else Check for minimum distance and insert the Point
        // If node is Leaf, use the Distance criteria else use
        // the wastage of area as criteria
        if (currNode->isLeaf()){

            // If 1st Node is half Full, fill int the Next 
            if (newNode1->childCnt == m){
                newNode2->addChildData(currPoint+d, 0, d);
                continue;
            }

            // If 2nd Node is half Full, fill in the first
            if (newNode2->childCnt == m){
                newNode1->addChildData(currPoint+d, 0, d);
                continue;
            }


            double dist1 = distance(currPoint+d,point1+d,d);
            double dist2 = distance(currPoint+d,point2+d,d);

            if (dist1<dist2){
                newNode1->addChildData(currPoint+d, 0, d);
            } else{
                newNode2->addChildData(currPoint+d,0, d);
            }
        }
        else{
            // Case of Non Leaf Node
            int insertedId = currNode->getChildID(i);

            // If 1st Node is half Full, fill int the Next 
            if (newNode1->childCnt == m){
                newNode2->addChildData(currPoint, insertedId, d);
                continue;
            }

            // If 2nd Node is half Full, fill in the first
            if (newNode2->childCnt == m){
                newNode1->addChildData(currPoint, insertedId, d);
                continue;
            }


            double wastage1 = distance(currPoint,point1,d);
            double wastage2 = distance(currPoint,point2,d);

            if (wastage1<wastage2){
                newNode1->addChildData(currPoint, insertedId, d);
            } else{
                newNode2->addChildData(currPoint,insertedId, d);
            }
        }

    }

    



    // Calculate The New MBR for both the nodes;
    newNode1->CalculateMBR(d,maxCap);
    newNode2->CalculateMBR(d,maxCap);

    //Copy Node1 to the destination of overflowed Node
    memcpy(data, newNode1,sizeof(Node));
    fh.MarkDirty(nodeId);
    fh.UnpinPage(nodeId);
    fh.FlushPage(nodeId);

    // Create a New Page and Store the Second Node into that
    ph = fh.NewPage();
    data = ph.GetData();
    memcpy(data,newNode2, sizeof(Node));
    fh.MarkDirty(ph.GetPageNum());
    fh.UnpinPage(ph.GetPageNum());
    fh.FlushPage(ph.GetPageNum());

    int parentId = currNode->parentID;
    
    // Split is done at the root level
    if (parentId < 0){
        Node *newParent = new Node(d, maxCap);
        newParent->parentID=-1;
        newParent->addChildData(newNode1->mbr, newNode1->ID,d);
        newParent->addChildData(newNode2->mbr, newNode2->ID,d);
        // Assign a New Id
        newParent->ID = GetNewId();
        IncreaseCnt();
        ph = fh.NewPage();
        root = newParent->ID;
        char *rootData = ph.GetData();
        memcpy(rootData, newParent, sizeof(Node));
        fh.MarkDirty(ph.GetPageNum());
        fh.UnpinPage(ph.GetPageNum());
        fh.FlushPage(ph.GetPageNum());
        fh.FlushPages();
        return updateMBR(newParent->ID);
    }

    // Splitting is done at non root level
    
    Node *Parent = new Node(d, maxCap);
    ph = fh.PageAt(parentId);
    char *pdata = ph.GetData();
    memcpy(Parent, pdata, sizeof(Node));

    // Update the MBR of new Node 1 which is just 
    // a replacement for overflowed Node
    for(int i=0; i<maxCap; i++){
        if (Parent->getChildID(i) == currNode->ID){
            int *prevPoint = Parent->getChildData(i,d);
            memcpy(prevPoint,newNode1->mbr,2*d*INT_SIZE); 
        }
    }


    // If the Parent has not maxCap children
    // then store the new Node's MBR into the 
    // empty child Space and call updateMBR
    if (Parent->childCnt<maxCap){
        Parent->addChildData(newNode2->mbr, newNode2->ID, d);
        memcpy(pdata, Parent, sizeof(Node));
        fh.FlushPages();
        return updateMBR(Parent->ID);
    }
    

    
    // This node will also result in the overflow
    // Paste the MBR of 2nd newly formed node in the end
    // and then recursively call the split node
    else{
        
        int *resSpace = Parent->getChildData(maxCap, d);
        Parent->addChildData(newNode2->mbr, newNode2->ID, d);
        //memcpy(resSpace, newNode2->mbr, 2*d*INT_SIZE);
        memcpy(pdata, Parent, sizeof(Node));
        fh.FlushPages();
        return SplitNode(Parent->ID);

    }


}

void Node::CalculateMBR(int d, int maxCap){
    int minPoint[d];
    int maxPoint[d];
    for(int j=0; j<d ;j++){
        minPoint[j]=INT_MAX;
        maxPoint[j]=INT_MIN;
    }

    int o=0;
    if (LeafNode == true)
        o=d;

    for (int i=0; i<maxCap; i++){
        // If there doesn't exist Point for that index
        if (getChildID(i)<0)
            continue;

        int *pointStart = getChildData(i,d);
        // Update Min Max Value
        for(int j=0; j<d ;j++){
                minPoint[j]=min(minPoint[j],*(pointStart+j+o));
                maxPoint[j]=max(maxPoint[j],*(pointStart+j+d));
        }
    }

    // Copy Calculated Value into the MBR field
    memcpy(mbr,minPoint,d*INT_SIZE);
    memcpy(mbr+d,maxPoint,d*INT_SIZE);
}



double Area(int *MinPoint, int *MaxPoint, int d){
    double area = 1;
    for(int i =0; i<d; i++){
        area=area*(MaxPoint[i]-MinPoint[i]);
    }
    return area;
}


// Calculate the Minimum Area Enlargement if A point is to be included into the Hyper Rectangle
double MinAreaEnlargement(int *point, int* MinPoint, int* MaxPoint, int d){
    for(int i=0; i< d; i++){
        if (point[i] < MinPoint[i])
            MinPoint[i]=point[i];
        if (point[i] > MaxPoint[i])
            MaxPoint[i]=point[i];
    }
    return Area(MinPoint,MaxPoint,d);
}

// Iterative Function To Calculate The Area 

int Node::selectChild(int *point, int d, int maxCap){
    // Todo :- Cross Check the data type for Storing area
   
    // Initialize the temporary Variables
    int index=-1;
    double minAreaEnlarge = DBL_MAX;
    double minArea = DBL_MAX;
    
    // Iterate over all the possible children
    for(int i=0; i <maxCap;i++){
        
        int cid = getChildID(i);
        // This is not a valid Child (this place is empty)
        if (cid < 0)
            continue;

        // get the child's data and calculate the Area Enlargement and Area
        // of the hyper rectangle for the comparison purpose
        int *cdata = getChildData(i, d);
        double cAreaEnlarge = MinAreaEnlargement(point,cdata,cdata+d,d);
        double cArea = Area(cdata, cdata+d, d);
        
        // Ignore the Child
        if (cAreaEnlarge> minAreaEnlarge)
            continue;
        
        // Set this is current child with Minimum Enlargement Area
        else if (cAreaEnlarge<minAreaEnlarge){
            index=i;
            minAreaEnlarge=cAreaEnlarge;
            minArea = cArea;
        } 

        // If the Area Enlargement is equal to current minimum 
        // then compare the Area of the hyper rectangle
        else {
            // If the actual area is less than that of previously
            // stored child, update the index and the minArea,
            // Minimum Area Enlargement is still the same 
            if (cArea < minArea){
                index=i;
                minArea = cArea;
            }

            // We don't update index if Area is equal to currently stored Minimum 
            // area as in that case the node with lower index is preffered which is 
            // by default the case here.
        }
    }
    return index;
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

bool isEqual(int *point1, int *point2, int d){
    for(int i=0; i<d; i++){
        if (*(point1+i)!= *(point2+i))
            return false;
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
            const char *fName;
            inFile >> fName;
            inFile >> N; 
            
            //FileHandler dummy = fm.OpenFile(fileName);
            //memcpy(fileName, &fName,fName.size());
        
            rootIndex = rtree.BulkLoad(fName, N);
            cout << "root= " << rootIndex << endl; 
            
            outFile << cmd;
            outFile<<"\n\n\n";
        } else if (cmd=="INSERT"){
            int PointArray[d];
            for (int i=0; i<d ;i++){
                inFile>>PointArray[i];
                //cout << PointArray[i] << " " ;
            }
            //cout << endl;

          //rtree.insert(PointArray,rtree.root);
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
    cout << fileName << endl;
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
        
        Node *newNode = new Node(d,maxCap);
        newNode->ID = GetNewId();
        IncreaseCnt();
        newNode->LeafNode=true;
        
        int minPoint[d];
        int maxPoint[d];
        for(int j=0; j<d ;j++){
            minPoint[j]=INT_MAX;
            maxPoint[j]=INT_MIN;
        }
        
        // reset all the child IDs to -1
        for(int j=0; j<maxCap; j++){
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


int Rtree::AssignParent(int start, int end){
    FileHandler fhIn, fhOut;
    PageHandler phIn, phOut;
    char *dataIn, *dataOut;
    fhIn = fm.OpenFile("treeData.txt");
    fhOut = fm.OpenFile("treeData.txt");
    if (end-start==1){
        root = start;
        return root;
    }
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
        parentNode->parentID = -1;
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
    return AssignParent(end, end+numParents);
}


int Rtree::insert(int *point, int rn){
    FileHandler fh = fm.OpenFile("treeData.txt");
    PageHandler ph = fh.PageAt(rn);
    char *data = ph.GetData();
    struct Node*nodeData = new Node(d, maxCap);
    memcpy(nodeData, data, sizeof(Node));

    
    if (nodeData->isLeaf()){
        // We are at leaf Node, Lets figure out Insertion 
        
        // If the node have less than maxCap points
        // then store point directly
        if (nodeData->childCnt < maxCap){
            nodeData->addChildData(point,0,d);
            
            // update the MBR, MBR could have changed only due to 
            // the addition of last point
            int *MBR = nodeData->mbr;
            for(int i=0; i<d; i++){
                MBR[i] = min(MBR[i],point[i]);
                MBR[i+d] = max(MBR[i+d], point[i]);                
            }    


            // Copy back the Node and flush Page
            memcpy(data, nodeData, sizeof(Node));
            fh.MarkDirty(rn);
            fh.UnpinPage(rn);
            fh.FlushPage(rn);

            // Recursively Update The MBR
            return updateMBR(nodeData->parentID, nodeData->ID);
        }

        // If there were maxCap children Already, then 
        // we need to figure out overflow
       
        
        // There was No extra Space in the Leaf Node, 
        // We need to copy the point in the residual space
        // and then Recursively Split the Nodes
        int *resSpace = nodeData->getChildData(maxCap,d);
        (nodeData->childID)[maxCap]=0;
        memcpy(resSpace+d, point, d*INT_SIZE);
        // Copy back the Node and flush Page
        memcpy(data, nodeData, sizeof(Node));
        fh.MarkDirty(rn);
        fh.UnpinPage(rn);
        fh.FlushPage(rn);

        return SplitNode(rn);
    }
    
    // If the node is not a leaf Node, the select the appropriate Child
    // and then recursively call the insert 
    else{
        int childID = nodeData->selectChild(point, d, maxCap);
        fh.FlushPages();
        return insert(point, childID);
    }

    return 1;
}

int Rtree::updateMBR(int nodeId, int childId){
    if (nodeId < 0)
        return childId;
    
    FileHandler fh = fm.OpenFile("treeData.txt");
    PageHandler ph = fh.PageAt(childId);
    Node *childNode = new Node(d, maxCap);
    char *data = ph.GetData();
    memcpy(childNode, data, sizeof(Node));
    fh.UnpinPage(childId);

    Node *mainNode = new Node(d, maxCap);
    ph = fh.PageAt(nodeId);
    data = ph.GetData();
    memcpy(mainNode, data, sizeof(Node)); 
    for (int i=0 ;i< maxCap; i++){
        // Select the proper index of Child under the Parent Node
        if (mainNode->getChildID(i) == childId){
            // Store the Child's New MBR under the Parent's Node
            int *childDataPos = mainNode->getChildData(i,d); 
            memcpy(childDataPos, childNode->mbr, 2*d*INT_SIZE);

            // Now update the Parent's MBR
            int *parentMBR = mainNode->mbr;
            for(int j=0; j<d; j++){
                parentMBR[i] = min(parentMBR[i], childDataPos[i]);
                parentMBR[i+d] = max(parentMBR[i+d], childDataPos[i+d]);
            }
            // Store the updated the Parent's Node
            memcpy(data, mainNode, sizeof(Node));
        }
    }

    fh.MarkDirty(nodeId);
    fh.UnpinPage(nodeId);
    fh.FlushPage(nodeId);

    // Current Parent Node is the Root Node
    if (mainNode->parentID == 0)
        return nodeId;

    return updateMBR(mainNode->parentID, nodeId);
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
            if (isEqual(point, pointData, d))
                return true;
        }
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
                    ans = ans|isPresent(point, childID);
                }
            }
        }
        return ans; 
    }
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