#ifndef INFINITIVE_VALUE
#define INFINITIVE_VALUE 1000000
#endif

#ifndef _C_GRAPH_LIB
#define _C_GRAPH_LIB

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "libfdr/jrb.h"
#include "libfdr/dllist.h"

#define DIRECT_GRAPH 1
#define UNDIRECT_GRAPH 0

typedef struct
{
	JRB edges;
	JRB vertices;
	int type;
} Graph;

typedef Dllist Queue;
typedef Dllist Stack;

Graph createGraph();
int numVertex(Graph g);
int addVertex(Graph graph, int id, char *name);
char *getVertex(Graph graph, int id);
int addEdge(Graph graph, int v1, int v2, double weight);
int hasEdge(Graph graph, int v1, int v2);
double getEdgeValue(Graph graph, int v1, int v2);
int indegree(Graph graph, int v, int *output);
int outdegree(Graph graph, int v, int *output);
void dropGraph(Graph graph);
double Dijkstra(Graph graph, int s, int t, int *path, int *length);
void MST(Graph); // IN RA FILE CAY KHUNG NHO NHAT
void Enqueue(Queue *q, Jval jval);
int Dequeue(Queue *q);
int isInQueue(Queue q, int i);
void Put(Stack *S, Jval jval);
int Pop(Stack *S);
int DAG(Graph graph);
void TSort(Graph g, int output[], int *n);
void DFS(Graph g, int start);
void BFS(Graph g, int start);
void Coloring(Graph g);
int numConnectedComponent(Graph g);
Graph createGraphReverse(Graph g);
void previsit(int v, int *countClock);
void postvisit(int v, int *countClock, JRB post);
int numStrongConnectComponent(Graph g);
void DFS_Edge_Reverse(Graph gReverse, int numV, int *countClock, JRB post);
int DFSForNumStrongConnect(Graph g, int numV, JRB post);
Graph importFile(char *fileName, int typeGraph);
void print2Dot(Graph g, char *fileName);

Graph importFile(char *fileName, int typeGraph)
{
	Graph g;
	int numV, numEgde, v1, v2;
	double w;
	FILE *f = fopen(fileName, "r");
	if (f == NULL)
	{
		printf("Can't open file %s\n", fileName);
		return g;
	}

	fscanf(f, "%d%*c", &numV);
	if (numV > 10000)
	{
		printf("Support Max Vertex: 10000!\n");
		fclose(f);
		return g;
	}
	if ((typeGraph != UNDIRECT_GRAPH) && (typeGraph != DIRECT_GRAPH))
	{
		printf("uncorrect type graph!\n");
		fclose(f);
		return g;
	}

	g = createGraph(typeGraph);
	fscanf(f, "%d%*c", &numEgde);
	while (!feof(f))
	{
		fscanf(f, "%d %d %lf%*c", &v1, &v2, &w);
		addEdge(g, v1, v2, w);
		addVertex(g, v1, "");
		addVertex(g, v2, "");
	}

	return g;
}

void print2Dot(Graph g, char *fileName)
{

	JRB Node, tmp, ptr;
	int v1, v2;
	int numV = numVertex(g);
	double w;

	FILE *fp;
	fp = fopen(fileName, "w");
	if (g.type == DIRECT_GRAPH)
		fprintf(fp, "digraph dothi\n{\n");
	else
		fprintf(fp, "graph dothi\n{\n");

	// IN RA các đỉnh
	jrb_traverse(Node, g.vertices)
		fprintf(fp, "\t%d\n", jval_i(Node->key));

	Graph printed;
	if (g.type == UNDIRECT_GRAPH)
	{
		printed = createGraph(UNDIRECT_GRAPH);
	}

	//Vẽ các cạnh
	jrb_traverse(Node, g.edges)
	{
		tmp = (JRB)jval_v(Node->val);
		jrb_traverse(ptr, tmp)
		{
			v1 = jval_i(Node->key);
			v2 = jval_i(ptr->key);
			w = jval_d(ptr->val);
			if (g.type == DIRECT_GRAPH)
				fprintf(fp, "\t%d -> %d [label=\"%lf\"]\n", v1, v2, w);
			else
			{
				//kiểm tra xem đã có cạnh v1 -- v2 trong graph printed hay chua
				if (!hasEdge(printed, v1, v2))
				{
					// Thêm cạnh đã in vào graph printed
					addEdge(printed, v1, v2, 0);
					fprintf(fp, "\t%d -- %d [label=\"%lf\"]\n", v1, v2, w);
				}
			}
		}
	}
	fprintf(fp, "}");

	if (g.type == UNDIRECT_GRAPH)
		dropGraph(printed);
	fclose(fp);
}

Graph createGraph(int type)
{
	if ((type != DIRECT_GRAPH) && (type != UNDIRECT_GRAPH))
	{
		printf("Please write correct type graph!!!\n");
		exit(0);
	}
	Graph g;
	g.edges = make_jrb();
	g.vertices = make_jrb();
	g.type = type;
	return g;
}

// Đếm số đỉnh của đò thị
int numVertex(Graph g)
{
	int numV = 0;
	JRB ptr;
	jrb_traverse(ptr, g.vertices)
		numV++;
	return numV;
}

int addVertex(Graph g, int id, char *name)
{
	JRB node = jrb_find_int(g.vertices, id);
	if (node == NULL)
	{
		jrb_insert_int(g.vertices, id, new_jval_s(strdup(name)));
		return 1;
	}
	return 0;
}

char *getVertex(Graph g, int id)
{
	JRB node = jrb_find_int(g.vertices, id);
	if (node == NULL)
		return NULL;
	else
		return jval_s(node->val);
}

int addEdge(Graph graph, int v1, int v2, double weight)
{
	JRB node, tree;
	if (getEdgeValue(graph, v1, v2) == INFINITIVE_VALUE)
	{
		node = jrb_find_int(graph.edges, v1);
		if (node == NULL)
		{
			tree = make_jrb();
			jrb_insert_int(graph.edges, v1, new_jval_v((JRB)tree));
		}
		else
		{
			tree = (JRB)jval_v(node->val);
		}
		jrb_insert_int(tree, v2, new_jval_d(weight));

		if (graph.type == UNDIRECT_GRAPH)
		{
			if (getEdgeValue(graph, v2, v1) == INFINITIVE_VALUE)
			{
				node = jrb_find_int(graph.edges, v2);
				if (node == NULL)
				{
					tree = make_jrb();
					jrb_insert_int(graph.edges, v2, new_jval_v((JRB)tree));
				}
				else
				{
					tree = (JRB)jval_v(node->val);
				}
				jrb_insert_int(tree, v1, new_jval_d(weight));
			}
		}
		return 1;
	}
	else return 0;
}

int hasEdge(Graph graph, int v1, int v2)
{
	JRB Node = jrb_find_int(graph.edges, v1);
	if (Node != NULL)
	{
		JRB tmp = jrb_find_int((JRB)jval_v(Node->val), v2);
		if (tmp != NULL)
		{
			return 1;
		}
		else
		{
			return 0;
		}
	}
}

double getEdgeValue(Graph graph, int v1, int v2)
{
	JRB node, tree;
	node = jrb_find_int(graph.edges, v1);
	if (node == NULL)
		return INFINITIVE_VALUE;
	tree = (JRB)jval_v(node->val);
	node = jrb_find_int(tree, v2);
	if (node == NULL)
		return INFINITIVE_VALUE;
	else
		return jval_d(node->val);
}

int indegree(Graph graph, int v, int *output)
{
	JRB tree, node;
	int total = 0;
	jrb_traverse(node, graph.edges)
	{
		tree = (JRB)jval_v(node->val);
		if (jrb_find_int(tree, v))
		{
			output[total] = jval_i(node->key);
			total++;
		}
	}
	return total;
}

int outdegree(Graph graph, int v, int *output)
{
	JRB tree, node;
	int total;
	node = jrb_find_int(graph.edges, v);
	if (node == NULL)
		return 0;
	tree = (JRB)jval_v(node->val);
	total = 0;
	jrb_traverse(node, tree)
	{
		output[total] = jval_i(node->key);
		total++;
	}
	return total;
}

double Dijkstra(Graph g, int s, int t, int *path, int *length)
{
	// Đếm số đỉnh của đồ thị
	int numV = numVertex(g);

	if (numV == 0)
		return INFINITIVE_VALUE;

	// isSelectMin dùng để xem đỉnh đó đã được chọn chưa (chọn min)
	// beenInPQ dùng để kiểm tra xem đỉnh đó đã ở trong hàng đợi hay chưa
	JRB isSelectMin = make_jrb(), beenInPQ = make_jrb(), distance = make_jrb(), previous = make_jrb(), distanceNode, previousNode, isSelectMinNode, isSelectMinU, isSelectMinV, beenInPQNode, beenInPQV, vertex, sNode, distanceU, distanceV, previousV;
	double min_dist, w, total;
	int min, u, tmp[numV], n, output[numV], v, idVertex;
	Dllist ptr, queue, node;

	jrb_traverse(vertex, g.vertices)
	{
		idVertex = jval_i(vertex->key);
		jrb_insert_int(isSelectMin, idVertex, new_jval_i(0));
		jrb_insert_int(beenInPQ, idVertex, new_jval_i(0));
		jrb_insert_int(distance, idVertex, new_jval_d(INFINITIVE_VALUE));
		jrb_insert_int(previous, idVertex, new_jval_i(-1));
	}

	// Dijkstra algorithm
	sNode = jrb_find_int(distance, s);
	sNode->val.d = 0;

	previousNode = jrb_find_int(previous, s);
	previousNode->val.i = s;

	queue = new_dllist();
	dll_append(queue, new_jval_i(s));

	while (!dll_empty(queue))
	{
		// get u from the priority queue
		min_dist = INFINITIVE_VALUE;
		dll_traverse(ptr, queue)
		{
			u = jval_i(ptr->val);
			distanceU = jrb_find_int(distance, u);
			if (min_dist > distanceU->val.d)
			{
				min_dist = distanceU->val.d;
				min = u;
				node = ptr;
			}
		}
		dll_delete_node(node);
		u = min;
		beenInPQNode = jrb_find_int(beenInPQ, u);
		beenInPQNode->val.i = 0;

		if (u == t)
			break;

		isSelectMinU = jrb_find_int(isSelectMin, u);
		isSelectMinU->val.i = 1;
		n = outdegree(g, u, output);
		for (int i = 0; i < n; i++)
		{
			v = output[i];
			isSelectMinV = jrb_find_int(isSelectMin, v);
			if (isSelectMinV->val.i == 0)
			{
				w = getEdgeValue(g, u, v);
				distanceU = jrb_find_int(distance, u);
				distanceV = jrb_find_int(distance, v);
				if (distanceV->val.d > distanceU->val.d + w)
				{
					distanceV->val.d = distanceU->val.d + w;
					previousV = jrb_find_int(previous, v);
					previousV->val.i = u;
				}
				beenInPQV = jrb_find_int(beenInPQ, v);
				if (beenInPQV->val.i == 0)
				{
					dll_append(queue, new_jval_i(v));
					beenInPQV->val.i = 1;
				}
			}
		}
	}

	distanceNode = jrb_find_int(distance, t);
	total = distanceNode->val.d;
	if (total != INFINITIVE_VALUE)
	{
		tmp[0] = t;
		n = 1;
		while (t != s)
		{
			previousNode = jrb_find_int(previous, t);
			t = previousNode->val.i;
			tmp[n++] = t;
		}
		for (int i = n - 1; i >= 0; i--)
			path[n - i - 1] = tmp[i];
		*length = n;
	}

	jrb_free_tree(isSelectMin);
	jrb_free_tree(beenInPQ);
	jrb_free_tree(distance);
	jrb_free_tree(previous);
	free_dllist(queue);

	return total;
}

void Enqueue(Queue *q, Jval jval)
{
	dll_append(*q, jval);
}

int Dequeue(Queue *q)
{
	Queue tmp = dll_first(*q);
	int v = jval_i(tmp->val);
	dll_delete_node(tmp);
	return v;
}

void Put(Stack *S, Jval jval)
{
	dll_append(*S, jval);
}

int Pop(Stack *S)
{
	Stack tmp = dll_last(*S);
	int v = jval_i(tmp->val);
	dll_delete_node(tmp);
	return v;
}

int Top(Stack S)
{
	if (!dll_empty(S))
	{
		return jval_i(dll_last(S)->val);
	}
	else
		return -1;
}

// Duyệt graph theo BFS nếu trong quá trình duyệt mà v == start có nghĩa đồ thị có chu trình trả vê là 0
// Duyệt hết graph mà không thấy v == start thì trả về 1 (có nghĩa không có chu trình)
int DAG(Graph graph)
{
	// Đếm số đỉnh của đồ thị
	int numV = numVertex(graph);

	if (numV == 0)
		return -1;

	Queue q = new_dllist();
	JRB visited = make_jrb(), vertex, visitedU, visitedV;
	int u, v;
	int output[numV];
	int re, start;

	jrb_traverse(vertex, graph.vertices)
		jrb_insert_int(visited, jval_i(vertex->key), new_jval_i(0));

	jrb_traverse(vertex, graph.vertices)
	{
		Enqueue(&q, vertex->key);
		start = jval_i(vertex->key);

		while (!dll_empty(q))
		{
			u = Dequeue(&q);
			visitedU = jrb_find_int(visited, u);
			if (visitedU->val.i == 0)
			{
				visitedU->val.i = 1;
				re = outdegree(graph, u, output);
				for (int i = 0; i < re; i++)
				{
					v = output[i];
					visitedV = jrb_find_int(visited, v);
					if (v == start)
					{
						return 0;
					}
					if (visitedV->val.i == 0)
					{
						Enqueue(&q, new_jval_i(v));
					}
				}
			}
		}
	}

	jrb_free_tree(visited);

	return 1;
}

// Sắp xếp Topo
void TSort(Graph g, int output[], int *n)
{
	// Đếm số đỉnh của đồ thị
	int numV = numVertex(g);

	if (numV == 0)
		return;

	if (!DAG(g))
	{
		printf("Detected Cycle\n");
		return;
	}

	Queue q = new_dllist();
	*n = 0;
	int v, out[numV], idVertex;
	JRB indeg = make_jrb(), vertex, indegNode;

	jrb_traverse(vertex, g.vertices)
	{
		idVertex = jval_i(vertex->key);
		jrb_insert_int(indeg, idVertex, new_jval_i(0));
	}

	jrb_traverse(vertex, g.vertices)
	{
		idVertex = jval_i(vertex->key);
		indegNode = jrb_find_int(indeg, idVertex);
		indegNode->val.i += indegree(g, idVertex, out);
		if (indegNode->val.i == 0)
		{
			Enqueue(&q, new_jval_i(idVertex));
		}
	}

	while (!dll_empty(q))
	{
		v = Dequeue(&q);
		output[(*n)++] = v;
		int k = outdegree(g, v, out);
		for (int i = 0; i < k; i++)
		{
			indegNode = jrb_find_int(indeg, out[i]);
			indegNode->val.i--;
			if (indegNode->val.i == 0)
				Enqueue(&q, new_jval_i(out[i]));
		}
	}

	free_dllist(q);

	jrb_free_tree(indeg);
}

// DFS
void DFS(Graph g, int start)
{
	// Đếm số đỉnh của đồ thị
	int numV = numVertex(g);

	if (numV == 0)
		return;

	Stack S = new_dllist();
	JRB visited = make_jrb(), vertex, visitedU, visitedNode, visitedV;
	int u, v, output[numV], n = 0, idVertex;

	jrb_traverse(vertex, g.vertices)
	{
		idVertex = jval_i(vertex->key);
		jrb_insert_int(visited, idVertex, new_jval_i(0));
	}

	// Kiểm tra đã thăm hết đỉnh chưa
	jrb_traverse(vertex, g.vertices)
	{
		idVertex = jval_i(vertex->key);
		visitedNode = jrb_find_int(visited, idVertex);
		if (jval_i(visitedNode->val) == 0)
		{
			Put(&S, new_jval_i(start));
			while (!dll_empty(S))
			{
				u = Pop(&S);
				visitedU = jrb_find_int(visited, u);
				if (visitedU->val.i == 0)
				{
					printf("%d ", u);
					visitedU->val.i = 1;
					n = outdegree(g, u, output);

					for (int i = n - 1; i >= 0; i--)
					{
						v = output[i];
						visitedV = jrb_find_int(visited, v);
						if (visitedV->val.i == 0)
						{
							Put(&S, new_jval_i(v));
						}
					}
				}
			}
		}
	}

	jrb_free_tree(visited);
	free_dllist(S);
}

//BFS
void BFS(Graph g, int start)
{
	// Đếm số đỉnh của đồ thị
	int numV = numVertex(g);

	if (numV == 0)
		return;

	Queue q = new_dllist();
	JRB visited = make_jrb(), vertex, visitedU, visitedNode, visitedV;
	int u, v, output[numV], n = 0, idVertex;

	jrb_traverse(vertex, g.vertices)
	{
		idVertex = jval_i(vertex->key);
		jrb_insert_int(visited, idVertex, new_jval_i(0));
	}

	jrb_traverse(vertex, g.vertices)
	{
		idVertex = jval_i(vertex->key);
		visitedNode = jrb_find_int(visited, idVertex);
		if (jval_i(visitedNode->val) == 0)
		{
			Enqueue(&q, new_jval_i(start));

			while (!dll_empty(q))
			{
				u = Dequeue(&q);
				visitedU = jrb_find_int(visited, u);
				if (visitedU->val.i == 0)
				{
					printf("%d ", u);
					visitedU->val.i = 1;
					n = outdegree(g, u, output);

					for (int i = 0; i < n; i++)
					{
						v = output[i];
						visitedV = jrb_find_int(visited, v);
						if (jval_i(visitedV->val) == 0)
						{
							Enqueue(&q, new_jval_i(v));
						}
					}
				}
			}
		}
	}

	jrb_free_tree(visited);
	free_dllist(q);
}

int isInQueue(Queue q, int i)
{
	Dllist ptr;
	dll_traverse(ptr, q)
	{
		if (ptr->val.i == i)
		{
			return 1;
		}
	}

	return 0;
}

// Sử dụng thuật toán Prim
// Thuật toán Prim chỉ áp dụng với đồ thị vô hướng liên thông
void MST(Graph g)
{
	if (g.type != UNDIRECT_GRAPH)
	{
		printf("Only for undirect graph!!!\n");
		return;
	}
	// Đếm số đỉnh của đồ thị
	int numV = numVertex(g);

	if (numV == 0)
		return;

	// Khai báo các biến cần thiết
	JRB isInMST = make_jrb(), isInMSTNode, vertex, distance = make_jrb(), distanceNode, distanceU, distanceV, previous = make_jrb(), previousNode, previousV, isInMSTU, isInMSTV;
	double min_dist, w;
	int min, u, n, output[numV], v;
	Dllist ptr, queue, node;

	// gán visit các đỉnh bằng 0 (có nghĩa là chưa thăm)
	jrb_traverse(vertex, g.vertices)
	{
		int idVertex = jval_i(vertex->key);
		// GÁN CHƯA THĂM CHO CÁC ĐỈNH
		jrb_insert_int(isInMST, idVertex, new_jval_i(0));
		// Khỏi tạo chi phí ban đầu đến các đỉnh là vô cùng
		jrb_insert_int(distance, idVertex, new_jval_d(INFINITIVE_VALUE));

		jrb_insert_int(previous, idVertex, new_jval_d(-1));
	}

	// Khỏi tạo chi phí ban đầu đến các đỉnh là vô cùng

	// Chọn một đỉnh đâu tiên trong cây JRB làm đỉnh bắt đầu
	int s = jval_i(jrb_first(g.edges)->key);

	// Gán chi phí đến đỉnh s là 0
	distanceNode = jrb_find_int(distance, s);
	distanceNode->val.d = 0;

	// Gán cha của s là s
	previousNode = jrb_find_int(previous, s);
	previousNode->val.i = s;

	// Make queue
	queue = new_dllist();

	// Append s
	dll_append(queue, new_jval_i(s));

	while (!dll_empty(queue))
	{
		// get u from the priority queue
		// Chọn đỉnh có chi phí nhỏ nhất
		min_dist = INFINITIVE_VALUE;
		dll_traverse(ptr, queue)
		{
			u = jval_i(ptr->val);
			distanceNode = jrb_find_int(distance, u);
			if (min_dist > distanceNode->val.d)
			{
				min_dist = distanceNode->val.d;
				min = u;
				node = ptr;
			}
		}

		dll_delete_node(node);
		u = min;

		// Đánh đấu là đã thăm (hay đã cho vào tập X)
		isInMSTU = jrb_find_int(isInMST, u);
		isInMSTU->val.i = 1;

		// Lấy các đỉnh kề của đỉnh u vào mảng output
		n = outdegree(g, u, output);

		// Lặp qua các đỉnh kê của u
		for (int i = 0; i < n; i++)
		{
			// Lấy ra ID của đỉnh kề thứ i của u
			v = output[i];

			// Lấy giá trị cạnh u->v
			w = getEdgeValue(g, u, v);

			// Gán lại chi phí đường đi tư u -> v
			distanceNode = jrb_find_int(distance, v);
			isInMSTV = jrb_find_int(isInMST, v);
			previousV = jrb_find_int(previous, v);
			if (distanceNode->val.d > w && isInMSTV->val.i == 0)
			{
				distanceNode->val.d = w;
				previousV->val.i = u;
				// Thêm v vào hàng đợi
				if (isInQueue(queue, v) == 0)
					dll_append(queue, new_jval_i(v));
			}
		}
	}

	// Export to dot file
	FILE *fp = fopen("MST.dot", "w");
	fprintf(fp, "graph MST\n{\n");

	// printed đc dùng để in ra file .dot
	Graph printed = createGraph(UNDIRECT_GRAPH);
	JRB Node, temp, ptr2;
	int v1, v2;

	// IN RA các đỉnh
	jrb_traverse(Node, g.vertices)
		fprintf(fp, "\t%d\n", jval_i(Node->key));

	// IN RA MST (Cây khung nhỏ nhất)
	jrb_traverse(ptr2, previous)
	{
		v1 = ptr2->key.i;
		v2 = ptr2->val.i;
		if (v1 == 0 && v2 == 0)
			continue;
		w = getEdgeValue(g, v1, v2);
		fprintf(fp, "\t%d -- %d [color=\"#ff0000\", label=\"%lf\", penwidth=3]\n", v2, v1, w);
		addEdge(printed, v1, v2, 0);
	}

	// In ra các cạnh còn lại
	jrb_traverse(Node, g.edges)
	{
		temp = (JRB)jval_v(Node->val);
		jrb_traverse(ptr2, temp)
		{
			v1 = jval_i(Node->key);
			v2 = jval_i(ptr2->key);

			// Kiểm tra xem đã có cạnh v1 -- v2 trong printed hay chưa (Kiểm tra xem đã in cạnh v1 -- v2 ra file .dot chưa)
			if (!hasEdge(printed, v1, v2))
			{
				w = jval_d(ptr2->val);

				// Thêm cạnh đã in vào graph printed
				addEdge(printed, v1, v2, 0);
				fprintf(fp, "\t%d -- %d [label=\"%lf\"]\n", v1, v2, w);
			}
		}
	}

	fprintf(fp, "}");
	fclose(fp);

	// Free
	jrb_free_tree(isInMST);
	free_dllist(queue);
	jrb_free_tree(distance);
	jrb_free_tree(previous);

	dropGraph(printed);
}

//To mau do thi
//Input: Do thi g
//Output: file .dot
void Coloring(Graph g)
{
	int m, n = 0, a, b, adj[500], nadj;
	int i = 0, *A, *B, *C, j, v1, v2;
	double w;
	JRB ptr, node, tmp;
	FILE *fp;

	//Đếm số đỉnh
	jrb_traverse(ptr, g.vertices)
		n++;

	A = (int *)calloc(n, sizeof(int));
	B = (int *)malloc(n * sizeof(int));
	C = (int *)malloc(n * sizeof(int));

	i = 0;
	jrb_traverse(ptr, g.vertices)
		C[i++] = jval_i(ptr->key);

	//Tìm số bậc của đỉnh, lưu vào A
	for (i = 0; i < n; ++i)
	{
		nadj = outdegree(g, C[i], adj);
		A[i] += nadj;
		for (j = 0; j < nadj; ++j)
			if (hasEdge(g, adj[j], C[i]) == 1)
				A[i]--;
		A[i] += indegree(g, C[i], adj);
	}

	//Copy mảng A sang mảng B
	for (i = 0; i < n; ++i)
		B[i] = A[i];
	//Tìm thứ tự đỉnh xếp theo bậc từ cao đến thấp, lưu vào A
	for (a = 0; a < n; ++a)
	{
		i = a;

		//Tìm đỉnh bậc cao nhất
		for (b = 0; b < n; b++)
			if (B[b] > B[i])
				i = b;

		B[i] = 0;
		A[a] = i;
	}

	B[A[0]] = 0; //Đỉnh đầu tiên có màu 0
	a = 1;

	//Tô màu các đỉnh còn lại
	for (i = 1; i < n; ++i)
	{
		B[A[i]] = 0; //Màu mặc định = 0

		for (b = 0; b < i; ++b)
			if ((hasEdge(g, C[A[i]], C[A[b]]) || hasEdge(g, C[A[b]], C[A[i]])) && B[A[i]] == B[A[b]])
			{
				//Nếu 2 đỉnh của 1 cạnh có cùng màu
				B[A[i]]++;
				b = -1;
				if (B[A[i]] == a)
					a++;
			}
	}

	fp = fopen("dothitomau.dot", "w");
	if (g.type == DIRECT_GRAPH)
		fprintf(fp, "digraph dothi\n{\n");
	else
		fprintf(fp, "graph dothi\n{\n");

	//Tô màu các đỉnh
	for (i = 0; i < n; ++i)
		fprintf(fp, "\t%d [fillcolor=\".%d .3 1.0\", style=filled];\n", C[i], B[i]);

	Graph printed;
	if (g.type == UNDIRECT_GRAPH)
	{
		printed = createGraph(UNDIRECT_GRAPH);
	}

	//Vẽ các cạnh
	jrb_traverse(node, g.edges)
	{
		tmp = (JRB)jval_v(node->val);
		jrb_traverse(ptr, tmp)
		{
			v1 = jval_i(node->key);
			v2 = jval_i(ptr->key);
			w = jval_d(ptr->val);
			if (g.type == DIRECT_GRAPH)
				fprintf(fp, "\t%d -> %d [label=\"%g\"]\n", v1, v2, w);
			else
			{
				//kiểm tra xem đã có cạnh v1 -- v2 trong graph printed hay chua
				if (!hasEdge(printed, v1, v2))
				{
					// Thêm cạnh đã in vào graph printed
					addEdge(printed, v1, v2, 0);
					fprintf(fp, "\t%d -- %d [label=\"%g\"]\n", v1, v2, w);
				}
			}
		}
	}
	fprintf(fp, "}");

	if (g.type == UNDIRECT_GRAPH)
		dropGraph(printed);

	free(A);
	free(B);
	free(C);
	fclose(fp);
	return;
}

// Đếm số thành phần liên thông của đồ thị
int numConnectedComponent(Graph g)
{
	// Đếm số đỉnh của đồ thị
	int numV = numVertex(g);

	if (numV == 0)
		return 0;

	Stack S = new_dllist();
	JRB visited = make_jrb(), visitedNode, vertex, visitedU, visitedV;
	int u, v;
	int adj[numV];
	int nFind = 0;
	int cc = 0;

	jrb_traverse(vertex, g.vertices)
	{
		int idVertex = jval_i(vertex->key);
		// GÁN CHƯA THĂM CHO CÁC ĐỈNH
		jrb_insert_int(visited, idVertex, new_jval_i(0));
	}

	jrb_traverse(vertex, g.vertices)
	{
		visitedNode = jrb_find_int(visited, jval_i(vertex->key));
		if (jval_i(visitedNode->val) == 0)
		{
			cc = cc + 1;
			Put(&S, vertex->key);

			while (!dll_empty(S))
			{
				u = Pop(&S);
				visitedU = jrb_find_int(visited, u);
				if (visitedU->val.i == 0)
				{
					visitedU->val.i = 1;
					nFind = outdegree(g, u, adj);
					for (int i = 0; i < nFind; i++)
					{
						v = adj[i];
						visitedV = jrb_find_int(visited, v);
						if (visitedV->val.i == 0)
						{
							Put(&S, new_jval_i(v));
						}
					}
				}
			}
		}
	}

	jrb_free_tree(visited);
	free_dllist(S);

	return cc;
}

// Tạo đồ thị cạnh ngược từ đồ thị g
Graph createGraphReverse(Graph g)
{
	Graph gReverse = createGraph(DIRECT_GRAPH);
	JRB node, tree, ptr;
	int v1, v2;
	double w;

	jrb_traverse(node, g.vertices)
		addVertex(gReverse, jval_i(node->key), "");

	jrb_traverse(node, g.edges)
	{
		tree = (JRB)jval_v(node->val);
		jrb_traverse(ptr, tree)
		{
			v1 = jval_i(node->key);
			v2 = jval_i(ptr->key);
			w = jval_d(ptr->val);
			addEdge(gReverse, v2, v1, w);
		}
	}

	return gReverse;
}

void previsit(int v, int *countClock)
{
	++(*countClock);
}

void postvisit(int v, int *countClock, JRB post)
{
	(*countClock)++;
	jrb_insert_int(post, *countClock, new_jval_i(v));
}

// DFS trên đồ thị ngược
void DFS_Edge_Reverse(Graph gReverse, int numV, int *countClock, JRB post)
{
	Stack S = new_dllist();
	JRB visited = make_jrb(), isInStack = make_jrb(), vertex, visitedS, isInStackNode, visitedU, visitedV;
	int u, v, start, adj[numV], nFind = 0;

	jrb_traverse(vertex, gReverse.vertices)
	{
		int idVertex = jval_i(vertex->key);
		// GÁN CHƯA THĂM CHO CÁC ĐỈNH
		jrb_insert_int(visited, idVertex, new_jval_i(0));

		// Chưa ở trong Stack
		jrb_insert_int(isInStack, idVertex, new_jval_i(0));
	}

	jrb_traverse(vertex, gReverse.vertices)
	{
		start = jval_i(vertex->key);
		visitedS = jrb_find_int(visited, start);
		if (visitedS->val.i == 0)
		{
			Put(&S, new_jval_i(start));
			isInStackNode = jrb_find_int(isInStack, start);
			isInStackNode->val.i = 1;

			while (!dll_empty(S))
			{
				u = Top(S);
				visitedU = jrb_find_int(visited, u);
				if (visitedU->val.i == 1)
				{
					int available = 0;
					nFind = outdegree(gReverse, u, adj);

					for (int i = 0; i < nFind; i++)
					{
						v = adj[i];
						visitedV = jrb_find_int(visited, v);
						isInStackNode = jrb_find_int(isInStack, v);
						if ((visitedV->val.i == 0) && (isInStackNode->val.i == 0))
						{
							Put(&S, new_jval_i(v));
							isInStackNode->val.i = 1;
							available = 1;
							break;
						}
					}
					if (available == 0)
					{
						postvisit(u, countClock, post);

						isInStackNode = jrb_find_int(isInStack, u);
						isInStackNode->val.i = 0;
						Pop(&S);
					}
				}
				else
				{
					visitedU->val.i = 1;
					previsit(u, countClock);
					nFind = outdegree(gReverse, u, adj);

					for (int i = 0; i < nFind; i++)
					{
						v = adj[i];
						visitedV = jrb_find_int(visited, v);
						isInStackNode = jrb_find_int(isInStack, v);
						if ((visitedV->val.i == 0) && (isInStackNode->val.i == 0))
						{
							Put(&S, new_jval_i(v));
							isInStackNode->val.i = 1;
							break;
						}
					}
				}
			}
		}
	}

	jrb_free_tree(visited);
	jrb_free_tree(isInStack);
	free_dllist(S);
}

// DFS trên đồ thị g
int DFSForNumStrongConnect(Graph g, int numV, JRB post)
{
	Stack S = new_dllist();
	JRB visited = make_jrb(), visitedJ, visitedU, visitedV, vertex, ptr;
	int u, v, adj[numV], nFind = 0, cc = 0, idVertex;

	jrb_traverse(vertex, g.vertices)
	{
		int idVertex = jval_i(vertex->key);
		// GÁN CHƯA THĂM CHO CÁC ĐỈNH
		jrb_insert_int(visited, idVertex, new_jval_i(0));
	}

	jrb_rtraverse(ptr, post)
	{
		int j = jval_i(ptr->val);
		visitedJ = jrb_find_int(visited, j);
		if (visitedJ->val.i == 0)
		{
			cc = cc + 1;
			Put(&S, new_jval_i(j));

			while (!dll_empty(S))
			{
				u = Pop(&S);
				visitedU = jrb_find_int(visited, u);
				if (visitedU->val.i == 0)
				{
					visitedU->val.i = 1;
					nFind = outdegree(g, u, adj);

					for (int i = nFind - 1; i >= 0; i--)
					{
						v = adj[i];
						visitedV = jrb_find_int(visited, v);
						if (visitedV->val.i == 0)
						{
							Put(&S, new_jval_i(v));
						}
					}
				}
			}
		}
	}

	jrb_free_tree(visited);
	free_dllist(S);

	return cc;
}

// Số thành phần liên thông mạnh
int numStrongConnectComponent(Graph g)
{
	// Đếm số đỉnh của đồ thị
	int numV = numVertex(g);

	if (numV == 0)
		return 0;

	int countClock = 0, re = 0;
	JRB post = make_jrb();

	Graph gReverse = createGraphReverse(g);

	DFS_Edge_Reverse(gReverse, numV, &countClock, post);
	re = DFSForNumStrongConnect(g, numV, post);

	dropGraph(gReverse);
	jrb_free_tree(post);

	return re;
}

void dropGraph(Graph graph)
{
	JRB node, tree;
	jrb_traverse(node, graph.edges)
	{
		tree = (JRB)jval_v(node->val);
		jrb_free_tree(tree);
	}

	jrb_free_tree(graph.edges);
	jrb_free_tree(graph.vertices);
}

#endif