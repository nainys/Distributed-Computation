#include <stdio.h>
#include <bits/stdc++.h>
#include <mpi.h>

using namespace std;
typedef struct Graph
{
	int V,E;
	int *edges;
} Graph;

Graph * createGraph(int V,int E)
{
	Graph* graph = new Graph;
	graph->V = V;
	graph->E = E;
	graph->edges = new int[E*3];
	return graph;
}

Graph * readGraph(string fileName)
{
	ifstream fin(fileName.c_str());
	int V,E;
	fin>>V>>E;
	Graph *g = createGraph(V,E);
	int src,dest,weight;
	for(int i = 0;i < E;i++)
	{
		fin>>src>>dest>>weight;
		g->edges[i * 3 + 0] = src-1;
		g->edges[i * 3 + 1] = dest-1;
		g->edges[i * 3 + 2] = weight;
	}
	fin.close();
	return g;
}

int find(int parent[], int i) 
{ 
    if (parent[i] == -1) 
        return i; 
    return find(parent, parent[i]); 
} 
   
void Union(int parent[], int x, int y) 
{ 
    int xset = find(parent, x); 
    int yset = find(parent, y);
    if(xset!=yset){ 
       parent[xset] = yset; 
    } 
} 


void mstBoruvka(Graph *graph,int *subsets)
{
	int rank,size;
	int V,E;
	int *edgeArray;
	int MSTweight = 0;
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	int vertexPerProc;
	int numTrees;
	
	if(rank == 0)
	{
		V = graph->V;
		E = graph->E;
	}

	// Printing Number of Edges in the MST
	if(rank == 0)
		cout<<V-1<<endl;
	
	// Broadcasting number of vertices and edges
	MPI_Bcast(&V,1,MPI_INT,0,MPI_COMM_WORLD);
	vertexPerProc = V/size;
	MPI_Bcast(&E,1,MPI_INT,0,MPI_COMM_WORLD);

	edgeArray = new int[E*3];
	if(rank == 0)
	{
		for(int i = 0;i<E*3;i++)
		{
			edgeArray[i] = graph->edges[i];
		}
	}

	//  Broadcasting Edges
	MPI_Bcast(edgeArray,E*3,MPI_INT,0,MPI_COMM_WORLD);

	// Vertex Partitioning
	int vLow,vHigh;
	vLow = rank*vertexPerProc;
	vHigh = ((rank+1)*vertexPerProc) - 1;
	numTrees = V;

	int* unionFind = new int[V];
	int *cheapest = new int[V];
	int* cheapestReceive;
	int* cheapestFinal;
	if(rank == 0)
	{
		for(int i=0;i<V;i++)
			unionFind[i] = subsets[i];
		cheapestReceive = new int[size*V];
		cheapestFinal = new int[V];
	}

	while(numTrees > 1)
	{	
		// Broadcasting Union Find Data Structure
		MPI_Bcast(unionFind,V,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
		for(int i = 0;i < V;i++)
			cheapest[i] = -1;

		// Finding Cheapest Edges for respective vertices
		for(int i = 0;i<E;i++)
		{
			int src,dest,weight;
			src = edgeArray[i * 3 + 0];
			dest = edgeArray[i * 3 + 1];
			weight = edgeArray[i * 3 + 2];
			// cout<<"\nTrying an Edge at Rank: "<<rank<<" with Source: "<<src<<" and Dest: "<<dest<<endl;
			// cout<<"\nvLow: "<<vLow<<" and vHigh: "<<vHigh<<endl;
			if((src >= vLow && src <= vHigh) || (dest >= vLow && dest <= vHigh))
			{
				int set1 = find(unionFind, src);
				int set2 = find(unionFind, dest);

				if(set1 == set2)
					continue;

				else
				{

					if(cheapest[set1] == -1 || edgeArray[cheapest[set1] * 3 + 2] > edgeArray[i*3 + 2])
						cheapest[set1] = i;

					if(cheapest[set2] == -1 || edgeArray[cheapest[set2]*3 + 2] > edgeArray[i*3 + 2])
						cheapest[set2] = i;
					
				}
			}
		}

		MPI_Barrier(MPI_COMM_WORLD);

		// cout<<"Cheapest array on Rank"<<rank<<endl;
		// for(int i=0;i<V;i++)
		// 	cout<<cheapest[i]<<" ";
		// cout<<endl;
		
		// Gathering Cheapest edges at the root
		MPI_Gather(cheapest, V, MPI_INT, cheapestReceive, V, MPI_INT, 0, MPI_COMM_WORLD);

		// if(rank == 0)
		// {
		// 	cout<<"Cheapest Receive Array"<<endl;
		// 	for(int i=0;i<size*V;i++)
		// 	{
		// 		cout<<cheapestReceive[i]<<" ";
		// 	}
		// 	cout<<endl;
		// }
		
		if(rank == 0)
		{
			for(int i=0;i<V;i++)
			{
				cheapestFinal[i] = -1;
				for(int j=0;j<size;j++)
				{
					int index = cheapestReceive[i+V*j];

					if(index == -1)
						continue;
					if(cheapestFinal[i] == -1 || edgeArray[cheapestFinal[i]*3 + 2] > edgeArray[index*3 + 2])
						cheapestFinal[i] = index;
				}
			}
		}

		// Adding cheapest edges to MST
		if(rank == 0)
		{
			for (int i=0; i<V; i++) 
			{ 
				// cout<<"Vertex - "<<i<<endl;
				if (cheapestFinal[i] != -1) 
				{ 
					int set1 = find(unionFind, edgeArray[cheapestFinal[i]*3]); 
					int set2 = find(unionFind, edgeArray[cheapestFinal[i]*3 + 1]); 
					// cout<<"Set 1: "<<set1<<" and Set 2: "<<set2<<endl;
					if (set1 == set2) 
						continue; 
					MSTweight += edgeArray[cheapestFinal[i]*3 + 2]; 
					// cout<<"MSTWeight: "<<MSTweight<<endl;
					// cout<<"Edge "<<edgeArray[cheapestReceive[i]*3]<<"-"<<edgeArray[cheapestReceive[i]*3 + 1]<<" with weight - "<<edgeArray[cheapestReceive[i]*3 + 2]<<" included in MST"<<endl; 
					cout<<edgeArray[cheapestFinal[i]*3]+1<<" "<<edgeArray[cheapestFinal[i]*3 + 1]+1<<endl;
					// Do a union of set1 and set2 and decrease number of trees 
					Union(unionFind, set1, set2); 
					// cout<<"\nPrinting UF Structure: ";
					// for(int i = 0;i<V;i++)
					// {
					// 	cout<<" "<<unionFind[i];
					// }
					// cout<<endl;
					numTrees--; 
				} 
			} 
		}
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Bcast(&numTrees,1,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Barrier(MPI_COMM_WORLD);
	}
	if(rank == 0)
		cout<<MSTweight<<endl;
}


int main(int argc, char** argv)
{
	int rank,size;
	Graph *graph;
	int *subsets;
	int V,E;

	string fileName(argv[1]);
	MPI_Init(NULL,NULL);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	if(rank == 0)
	{
		graph = readGraph(fileName);
		V = graph->V;
		E = graph->E;
		subsets = new int[V];
		for(int i = 0;i < V;i++)
			subsets[i] = -1;
	}
	MPI_Barrier(MPI_COMM_WORLD);
	mstBoruvka(graph,subsets);
	MPI_Finalize();
}
