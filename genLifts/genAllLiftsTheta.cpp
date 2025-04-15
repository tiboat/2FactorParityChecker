#include <bits/stdc++.h>

// Author: Jorik Jooken

using namespace std;

vector< set<int> > graph;
vector< vector<int> > adjMatrix;
int n;
int requiredGirth;
int nbGroupElements;
vector< vector<int> > groupMultiplicationTable;

void writeAdjMatrix()
{
    cout << "n=" << n << endl;
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

vector< pair<int, int> > remainingEdges;

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
            int idx1=el*2+remainingEdges[edgeIndex].first;
            int idx2=newEl*2+remainingEdges[edgeIndex].second;
            if(idx1==idx2 || adjMatrix[idx1][idx2]==1)
            {
                loops=true;
                break;
            }
            adjMatrix[idx1][idx2]=adjMatrix[idx2][idx1]=1;
            graph[idx1].insert(idx2);
            graph[idx2].insert(idx1);
        }
        if(loops) continue;
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
            int idx1=el*2+remainingEdges[edgeIndex].first;
            int idx2=newEl*2+remainingEdges[edgeIndex].second;
            adjMatrix[idx1][idx2]=adjMatrix[idx2][idx1]=0;
            graph[idx1].erase(idx2);
            graph[idx2].erase(idx1);
        }
    }
    // cerr << "nbGroupElements and suitableElements.size(): " << nbGroupElements << " " << suitableElements.size() << endl;
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
            int idx1=el*2+remainingEdges[edgeIndex].first;
            int idx2=newEl*2+remainingEdges[edgeIndex].second;
            if(idx1==idx2 || adjMatrix[idx1][idx2]==1)
            {
                loops=true;
                break;
            }
            adjMatrix[idx1][idx2]=adjMatrix[idx2][idx1]=1;
            graph[idx1].insert(idx2);
            graph[idx2].insert(idx1);
        }
        if(loops) continue;
        int currGirth=-1;
        if(!loops) currGirth=calcGirthUpperBound();
        if(currGirth<requiredGirth || loops)
        {
            ok=false;
        }
        if(ok)
        {
            recursivelyAssign(edgeIndex+1);
        }
        // undo top
        for(int el=0; el<nbGroupElements; el++)
        {
            int newEl=groupMultiplicationTable[el][topElGroup];
            int idx1=el*2+remainingEdges[edgeIndex].first;
            int idx2=newEl*2+remainingEdges[edgeIndex].second;
            adjMatrix[idx1][idx2]=adjMatrix[idx2][idx1]=0;
            graph[idx1].erase(idx2);
            graph[idx2].erase(idx1);
        }
    }
}


// lifts of Theta graph: cubic graph with 2 vertices and 3 parallel edges
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
        
        // cerr << "nbGroupElements: " << nbGroupElements << endl;

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
        if(nbGroupElements!=15) continue;

        n=nbGroupElements*2;
        adjMatrix.assign(n, vector<int>(n,0));
        set<int> emp;
        graph.assign(n,emp);

        vector<int> orderVec;

        orderVec.push_back(0);
        orderVec.push_back(1);
        // 1 edge in spanning tree
        for(int el=0; el<nbGroupElements; el++)
        {
            for(int i=0; i<=0; i++)
            {
                int idx1=el*2+orderVec[i];
                int idx2=el*2+orderVec[i+1];
                adjMatrix[idx1][idx2]=adjMatrix[idx2][idx1]=1;
                graph[idx1].insert(idx2);
                graph[idx2].insert(idx1);
            }
        }

        remainingEdges.clear();
        remainingEdges.push_back(make_pair(0,1));
        remainingEdges.push_back(make_pair(0,1));
        
        // cerr << orderVec.size() << " " << remainingEdges.size() << endl;

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
