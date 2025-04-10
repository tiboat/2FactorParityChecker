#include <bits/stdc++.h>


// Author: Jorik Jooken

using namespace std;

vector< set<int> > graph;
vector< vector<int> > adjMatrix;
int n;
int requiredGirth;
int nbGroupElements;
vector< vector<int> > groupMultiplicationTable;

vector<int> orderVec;
vector< pair<int, int> > remainingEdges;
vector<int> choicesMade;

void writeAdjMatrix()
{
    cout << "Lift of Gray graph of girth at least " << requiredGirth << ":" << endl;
    for(int i=0; i<n; i++)
    {
        for(int j=0; j<n; j++)
        {
            cout << adjMatrix[i][j] << " ";
        }
        cout << endl;
    }
}

int calcShortestPathErasedEdge(int u, int v)
{
    vector<int> dist(n,1e9);
    dist[u]=0;
    queue<int> q;
    q.push(u);
    while(!q.empty())
    {
        int now=q.front();
        q.pop();
        for(int neigh : graph[now])
        {
            if(now==u && neigh==v) continue;
            if(now==v && neigh==u) continue;
            if(dist[neigh]<1e9) continue;
            dist[neigh]=dist[now]+1;
            q.push(neigh);
            if(neigh==v) return dist[neigh];
        }
    }
    return dist[v];
}

int calcGirthUpperBound()
{
    int ret=1e9;
    for(int u=0; u<n; u++)
    {
        for(int v : graph[u])
        {
            if(v<u) continue;
            int dist=calcShortestPathErasedEdge(u,v);
            ret=min(ret,dist+1);
            if(ret<requiredGirth)
            {
                return ret;
            }
        }
    }
    return ret;
}

set<int> determineSuitableElements(int edgeIndex)
{
    set<int> suitableElements;
    for(int topElGroup=0; topElGroup<nbGroupElements; topElGroup++)
    {
        bool ok=true;

        bool loops=false;
        // top
        for(int el=0; el<nbGroupElements; el++)
        {
            int newEl=groupMultiplicationTable[el][topElGroup];
            int revEl=groupMultiplicationTable[newEl][topElGroup];
            int idx1=el*54+remainingEdges[edgeIndex].first;
            int idx2=newEl*54+remainingEdges[edgeIndex].second;
            adjMatrix[idx1][idx2]=adjMatrix[idx2][idx1]=1;
            graph[idx1].insert(idx2);
            graph[idx2].insert(idx1);
        }
        int currGirth=-1;
        if(!loops) currGirth=calcGirthUpperBound();
        if(currGirth<requiredGirth || loops)
        {
            ok=false;
        }
        if(ok) suitableElements.insert(topElGroup);
        // undo top
        for(int el=0; el<nbGroupElements; el++)
        {
            int newEl=groupMultiplicationTable[el][topElGroup];
            int idx1=el*54+remainingEdges[edgeIndex].first;
            int idx2=newEl*54+remainingEdges[edgeIndex].second;
            adjMatrix[idx1][idx2]=adjMatrix[idx2][idx1]=0;
            graph[idx1].erase(idx2);
            graph[idx2].erase(idx1);
        }
    }
    
    return suitableElements;
}

vector< set<int> > suitableElementsPerEdge;

void recursivelyAssign(int edgeIndex)
{
    if(edgeIndex==remainingEdges.size())
    {
        writeAdjMatrix();
        return;
    }
    for(int topElGroup : suitableElementsPerEdge[edgeIndex])
    {
        bool ok=true;

        bool loops=false;
        // top
        for(int el=0; el<nbGroupElements; el++)
        {
            int newEl=groupMultiplicationTable[el][topElGroup];
            int revEl=groupMultiplicationTable[newEl][topElGroup];
            int idx1=el*54+remainingEdges[edgeIndex].first;
            int idx2=newEl*54+remainingEdges[edgeIndex].second;
            adjMatrix[idx1][idx2]=adjMatrix[idx2][idx1]=1;
            graph[idx1].insert(idx2);
            graph[idx2].insert(idx1);
        }
        int currGirth=-1;
        if(!loops) currGirth=calcGirthUpperBound();
        if(currGirth<requiredGirth || loops)
        {
            ok=false;
        }
        if(ok)
        {
            choicesMade.push_back(topElGroup);
            recursivelyAssign(edgeIndex+1);
            choicesMade.pop_back();
        }
        
        // undo top
        for(int el=0; el<nbGroupElements; el++)
        {
            int newEl=groupMultiplicationTable[el][topElGroup];
            int idx1=el*54+remainingEdges[edgeIndex].first;
            int idx2=newEl*54+remainingEdges[edgeIndex].second;
            adjMatrix[idx1][idx2]=adjMatrix[idx2][idx1]=0;
            graph[idx1].erase(idx2);
            graph[idx2].erase(idx1);
        }
    }
}


// lifts of Gray graph
// us??????????B??o@C?GG@A?GO?B??CG?AG??S??@G??B????_???G???A????G???A?????O???@????A????A????O??????G????@??????_@????_A????B?????B??????CG?????a??????Q?????AG??????@_?????@_???????g??????B??O?????AOO?????@_A??????A`??????@OC?????c??_????@_??
int main(int argc, char *argv[])
{
    if(argc!=2)
    {
        fprintf(stderr,"Wrong number of command line arguments!\nExpected girthRequired\n");
    }
    requiredGirth=atoi(argv[1]);

    ios::sync_with_stdio(false);
    cin.tie(0);
    string line;
    while(getline(cin,line))
    {
        nbGroupElements=atoi(line.c_str());

        vector<int> foo(nbGroupElements,-1);
        groupMultiplicationTable.assign(nbGroupElements,foo);
        for(int i=0; i<nbGroupElements; i++)
        {
            getline(cin,line);
            stringstream ss;
            ss << line;
            for(int j=0; j<nbGroupElements; j++)
            {
                int tmp;
                ss >> tmp;
                groupMultiplicationTable[i][j]=tmp-1;
            }
        }
        if(nbGroupElements!=3) continue;

        n=nbGroupElements*54;
        adjMatrix.assign(n, vector<int>(n,0));
        set<int> emp;
        graph.assign(n,emp);


        orderVec.push_back(0);
        orderVec.push_back(2);
        orderVec.push_back(53);
        orderVec.push_back(38);
        orderVec.push_back(28);        
        orderVec.push_back(42);
        orderVec.push_back(49);
        orderVec.push_back(1);
        orderVec.push_back(48);
        orderVec.push_back(43);
        orderVec.push_back(29);
        orderVec.push_back(45);
        orderVec.push_back(50);
        orderVec.push_back(47);
        orderVec.push_back(34);
        orderVec.push_back(37);
        orderVec.push_back(24);
        orderVec.push_back(36);

	    orderVec.push_back(35);
        orderVec.push_back(22);
        orderVec.push_back(11);
        orderVec.push_back(15);
        orderVec.push_back(26);        
        orderVec.push_back(41);
        orderVec.push_back(30);
        orderVec.push_back(20);
        orderVec.push_back(6);
        orderVec.push_back(19);
        orderVec.push_back(31);
        orderVec.push_back(44);
        orderVec.push_back(51);
        orderVec.push_back(46);
        orderVec.push_back(33);
        orderVec.push_back(13);
        orderVec.push_back(7);
        orderVec.push_back(21);

	    orderVec.push_back(9);
        orderVec.push_back(23);
        orderVec.push_back(10);
        orderVec.push_back(16);
        orderVec.push_back(5);        
        orderVec.push_back(12);
        orderVec.push_back(4);
        orderVec.push_back(17);
        orderVec.push_back(27);
        orderVec.push_back(40);
        orderVec.push_back(32);
        orderVec.push_back(18);
        orderVec.push_back(8);
        orderVec.push_back(14);
        orderVec.push_back(25);
        orderVec.push_back(39);
        orderVec.push_back(52);
        orderVec.push_back(3);
        // tree gets identity
        // 5 edges in spanning tree
        for(int el=0; el<nbGroupElements; el++)
        {
            for(int i=0; i<=52; i++)
            {
                int idx1=el*54+orderVec[i];
                int idx2=el*54+orderVec[i+1];
                adjMatrix[idx1][idx2]=adjMatrix[idx2][idx1]=1;
                graph[idx1].insert(idx2);
                graph[idx2].insert(idx1);
            }
        }
        remainingEdges.clear();
        remainingEdges.push_back(make_pair(0,1));
	    remainingEdges.push_back(make_pair(0,3));
	    remainingEdges.push_back(make_pair(2,51));
	    remainingEdges.push_back(make_pair(3,50));
	    remainingEdges.push_back(make_pair(4,14));
        remainingEdges.push_back(make_pair(5,15));
        remainingEdges.push_back(make_pair(6,13));
        remainingEdges.push_back(make_pair(7,18));
        remainingEdges.push_back(make_pair(8,22));
        remainingEdges.push_back(make_pair(9,17));
        remainingEdges.push_back(make_pair(10,20));
        remainingEdges.push_back(make_pair(11,19));
        remainingEdges.push_back(make_pair(12,24));
        remainingEdges.push_back(make_pair(16,28));
        remainingEdges.push_back(make_pair(21,29));
        remainingEdges.push_back(make_pair(23,34));
        remainingEdges.push_back(make_pair(25,43));
        remainingEdges.push_back(make_pair(26,39));
        remainingEdges.push_back(make_pair(27,38));
        remainingEdges.push_back(make_pair(30,45));
        remainingEdges.push_back(make_pair(31,42));
        remainingEdges.push_back(make_pair(32,44));
        remainingEdges.push_back(make_pair(33,47));
        remainingEdges.push_back(make_pair(35,46));
        remainingEdges.push_back(make_pair(36,52));
        remainingEdges.push_back(make_pair(37,53));
        remainingEdges.push_back(make_pair(40,48));
        remainingEdges.push_back(make_pair(41,49));

        // determine suitable elements per edge
        suitableElementsPerEdge.clear();
        for(int i=0; i<remainingEdges.size(); i++)
        {
            set<int> suitableElements=determineSuitableElements(i);
            suitableElementsPerEdge.push_back(suitableElements);
        }

        // recursively try to assign group elements to the remaining edges            
        recursivelyAssign(0);
    }
    return 0;
}
