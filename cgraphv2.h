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
void Enqueue(Queue *q, Jval jval);
int Dequeue(Queue *q);
int isInQueue(Queue q, int i);
void Put(Stack *S, Jval jval);
int Pop(Stack *S);
int DAG(Graph graph);
void TSort(Graph g, int output[], int *n);
void DFS(Graph g, int start);
void BFS(Graph g, int start);


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