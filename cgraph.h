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
void addVertex(Graph graph, int id, char *name);
char *getVertex(Graph graph, int id);
void addEdge(Graph graph, int v1, int v2, double weight);
int hasEdge(Graph graph, int v1, int v2);
double getEdgeValue(Graph graph, int v1, int v2);
int indegree(Graph graph, int v, int *output);
int outdegree(Graph graph, int v, int *output);
void dropGraph(Graph graph);
void MST(Graph); // IN RA FILE CAY KHUNG NHO NHAT
void Enqueue(Queue *q, Jval jval);
int Dequeue(Queue *q);
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
void hamBFS(Graph g, int source, int dest);
Dllist outgoingVertices(Graph g, int v);
void hamDFS(Graph g, int source, int dest);
Dllist dll_min(Dllist list, JRB cost, JRB visited);
double shortestPath(Graph g, int source, int dest);
Dllist findPath(JRB parent, int source, int dest);
void topologicalSort(Graph g);
void recursionDFS(Graph g, int source, JRB visited, Dllist stack);
void TPLT(Graph g);
void recursionDFSdao(Graph g, int source, JRB visited, Dllist stack);
Dllist incomingVertices(Graph g, int v);
void SCC(Graph g);
int isInQueue(Queue q, int i);


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

void addVertex(Graph g, int id, char *name)
{
    JRB node = jrb_find_int(g.vertices, id);
    if (node == NULL) // only add new vertex
        jrb_insert_int(g.vertices, id, new_jval_s(strdup(name)));
}

char *getVertex(Graph g, int id)
{
    JRB node = jrb_find_int(g.vertices, id);
    if (node == NULL)
        return NULL;
    else
        return jval_s(node->val);
}

void addEdge(Graph graph, int v1, int v2, double weight)
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
    }

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
    int visited[numV];
    int u, v;
    int output[numV];
    int re, start;
    JRB vertex;

    jrb_traverse(vertex, graph.vertices)
    {
        memset(visited, 0, sizeof(visited));
        Enqueue(&q, vertex->key);
        start = jval_i(vertex->key);

        while (!dll_empty(q))
        {
            u = Dequeue(&q);
            if (!visited[u])
            {
                visited[u] = 1;
                re = outdegree(graph, u, output);
                for (int i = 0; i < re; i++)
                {
                    v = output[i];
                    if (v == start)
                    {
                        return 0;
                    }
                    if (visited[v] == 0)
                    {
                        Enqueue(&q, new_jval_i(v));
                    }
                }
            }
        }
    }
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
    JRB tmp;
    int v;
    int out[numV];
    int *indeg = (int *)calloc(numV, sizeof(int));
    jrb_traverse(tmp, g.vertices)
    {
        indeg[jval_i(tmp->key)] += indegree(g, jval_i(tmp->key), out);
        if (indeg[jval_i(tmp->key)] == 0)
        {
            Enqueue(&q, tmp->key);
        }
    }

    while (!dll_empty(q))
    {
        v = Dequeue(&q);
        output[(*n)++] = v;
        int k = outdegree(g, v, out);
        for (int i = 0; i < k; i++)
        {
            indeg[out[i]]--;
            if (indeg[out[i]] == 0)
                Enqueue(&q, new_jval_i(out[i]));
        }
    }

    free_dllist(q);

    free(indeg);
}

// DFS
void DFS(Graph g, int start)
{
    // Đếm số đỉnh của đồ thị
    int numV = numVertex(g);

    if (numV == 0)
        return;

    Stack S = new_dllist();
    int *visited = (int *)calloc(numV, sizeof(int));
    int u;
    int output[numV];
    int n = 0;

    // Kiểm tra đã thăm hết đỉnh chưa
    for (int j = 0; j < numV; j++)
    {
        if (visited[j] == 0)
        {

            Put(&S, new_jval_i(start));

            while (!dll_empty(S))
            {
                u = Pop(&S);
                if (visited[u] == 0)
                {
                    printf("%d ", u);
                    visited[u] = 1;
                    n = outdegree(g, u, output);

                    for (int i = 0; i < n; i++)
                    {
                        int v = output[i];
                        if (visited[v] == 0)
                        {
                            Put(&S, new_jval_i(v));
                        }
                    }
                }
            }
        }
    }

    free(visited);
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
    int *visited = (int *)calloc(numV, sizeof(int));
    int u;
    int output[numV];
    int n = 0;

    // Dùng vòng lặp for kiểm tra xem đã visit hết đỉnh chưa
    for (int j = 0; j < numV; j++)
    {
        if (visited[j] == 0)
        {
            Enqueue(&q, new_jval_i(start));

            while (!dll_empty(q))
            {
                u = Dequeue(&q);
                if (visited[u] == 0)
                {
                    printf("%d ", u);
                    visited[u] = 1;
                    n = outdegree(g, u, output);

                    for (int i = 0; i < n; i++)
                    {
                        int v = output[i];
                        if (visited[v] == 0)
                        {
                            Enqueue(&q, new_jval_i(v));
                        }
                    }
                }
            }
        }
    }

    free(visited);
    free_dllist(q);
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
    int *visited = (int *)calloc(numV, sizeof(int));
    int u, v;
    int adj[numV];
    int nFind = 0;
    int cc = 0;

    for (int j = 0; j < numV; j++)
    {
        if (visited[j] == 0)
        {
            cc = cc + 1;
            Put(&S, new_jval_i(j));

            while (!dll_empty(S))
            {
                u = Pop(&S);
                if (visited[u] == 0)
                {
                    visited[u] = 1;
                    nFind = outdegree(g, u, adj);
                    for (int i = 0; i < nFind; i++)
                    {
                        v = adj[i];
                        if (visited[v] == 0)
                        {
                            Put(&S, new_jval_i(v));
                        }
                    }
                }
            }
        }
    }

    free(visited);
    free_dllist(S);

    return cc;
}

// Tạo đồ thị cạnh ngược từ đồ thị g
Graph createGraphReverse(Graph g)
{
    Graph gReverse = createGraph(DIRECT_GRAPH);
    JRB node, temp, ptr;
    int v1, v2;
    double w;

    jrb_traverse(node, g.edges)
    {
        temp = (JRB)jval_v(node->val);
        jrb_traverse(ptr, temp)
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
    /* preArr[v] =  */ ++(*countClock);
}

void postvisit(int v, int *countClock, JRB post)
{
    (*countClock)++;
    jrb_insert_int(post, *countClock, new_jval_i(v));
    // postArr[v] = clock;
}

// DFS trên đồ thị ngược
void DFS_Edge_Reverse(Graph gReverse, int numV, int *countClock, JRB post)
{
    Stack S = new_dllist();
    int *visited = (int *)calloc(numV, sizeof(int));
    int *isInStack = (int *)calloc(numV, sizeof(int));
    int u, v;
    int adj[numV];
    int nFind = 0;

    for (int j = 0; j < numV; j++)
    {
        if (visited[j] == 0)
        {
            Put(&S, new_jval_i(j));
            isInStack[j] = 1;

            while (!dll_empty(S))
            {
                u = Top(S);
                if (visited[u] == 1)
                {
                    int available = 0;
                    nFind = outdegree(gReverse, u, adj);

                    for (int i = 0; i < nFind; i++)
                    {
                        v = adj[i];
                        if (visited[v] == 0 && isInStack[v] == 0)
                        {
                            Put(&S, new_jval_i(v));
                            isInStack[v] = 1;
                            available = 1;
                            break;
                        }
                    }
                    if (available == 0)
                    {
                        // printf("%d\n", u);
                        postvisit(u, countClock, post);
                        // printf("%d\n", clock);
                        Pop(&S);
                    }
                }
                else
                {
                    visited[u] = 1;
                    // printf("%d\n", u);

                    previsit(u, countClock);
                    // printf("%d\n", clock);
                    nFind = outdegree(gReverse, u, adj);

                    for (int i = 0; i < nFind; i++)
                    {
                        v = adj[i];
                        if (visited[v] == 0 && isInStack[v] == 0)
                        {
                            Put(&S, new_jval_i(v));
                            isInStack[v] = 1;
                            break;
                        }
                    }
                }
            }
            // printf("------------------\n");
        }
    }

    // printf("PostLast: %d\n", clock);

    free(visited);
    free(isInStack);
    free_dllist(S);
}

// DFS trên đồ thị g
int DFSForNumStrongConnect(Graph g, int numV, JRB post)
{
    Stack S = new_dllist();
    int *visited = (int *)calloc(numV, sizeof(int));
    int u, v;
    int adj[numV];
    int nFind = 0;
    int cc = 0;
    JRB tmp;

    jrb_rtraverse(tmp, post)
    {
        // printf("%d\n", key);
        int j = jval_i(tmp->val);
        if (visited[j] == 0)
        {
            cc = cc + 1;
            Put(&S, new_jval_i(j));

            while (!dll_empty(S))
            {
                u = Pop(&S);
                if (visited[u] == 0)
                {
                    visited[u] = 1;
                    // printf("%d\n", u);
                    nFind = outdegree(g, u, adj);

                    for (int i = 0; i < nFind; i++)
                    {
                        v = adj[i];
                        if (visited[v] == 0)
                        {
                            Put(&S, new_jval_i(v));
                        }
                    }
                }
            }
            // printf("------------------\n");
        }
    }

    free(visited);
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

void hamBFS(Graph g, int source, int dest)
{
    JRB visited = make_jrb();
    Dllist q = new_dllist();
    Dllist n;
    Dllist adjs;
    int i, u;

    dll_append(q, new_jval_i(source));
    while (!dll_empty(q))
    {
        n = dll_first(q);
        u = jval_i(dll_val(n));
        dll_delete_node(n);
        if (jrb_find_int(visited, u) == NULL)
        {
            printf("%d ", u);
            jrb_insert_int(visited, u, new_jval_i(1));
            adjs = outgoingVertices(g, u);

            dll_traverse(n, adjs)
            {
                if (jrb_find_int(visited, jval_i(dll_val(n))) == NULL)
                    dll_append(q, dll_val(n));
                if (jval_i(dll_val(n)) == dest)
                {
                    free_dllist(adjs);
                    jrb_free_tree(visited);
                    free_dllist(q);
                    return;
                }
            }
            free_dllist(adjs);
        }
    }
    printf("\n");
    jrb_free_tree(visited);
    free_dllist(q);
}

void hamDFS(Graph g, int source, int dest)
{
    JRB visited = make_jrb();
    Dllist q = new_dllist();
    Dllist n;
    Dllist adjs;
    int i, u;

    dll_append(q, new_jval_i(source));
    while (!dll_empty(q))
    {
        n = dll_last(q);
        u = jval_i(dll_val(n));
        dll_delete_node(n);
        if (jrb_find_int(visited, u) == NULL)
        {
            printf("%d ", u);
            jrb_insert_int(visited, u, new_jval_i(1));
            adjs = outgoingVertices(g, u);

            dll_rtraverse(n, adjs)
            {
                if (jrb_find_int(visited, jval_i(dll_val(n))) == NULL)
                    dll_append(q, dll_val(n));
                if (jval_i(dll_val(n)) == dest)
                {
                    free_dllist(adjs);
                    jrb_free_tree(visited);
                    free_dllist(q);
                    return;
                }
            }
            free_dllist(adjs);
        }
    }
    printf("\n");
    jrb_free_tree(visited);
    free_dllist(q);
}

Dllist outgoingVertices(Graph g, int v)
{
    JRB node, l, i;
    Dllist output = new_dllist();

    if ((node = jrb_find_int(g.edges, v)) == NULL)
        return output;

    l = (JRB)jval_v(node->val);
    jrb_traverse(i, l)
    {
        dll_append(output, i->key);
    }
    return output;
}
Dllist incomingVertices(Graph g, int v)
{
    JRB node, l, i;
    Dllist output = new_dllist();

    jrb_traverse(node, g.edges)
    {
        l = (JRB)jval_v(node->val);
        if (jrb_find_int(l, v) != NULL)
            dll_append(output, node->key);
    }
    return output;
}
double shortestPath(Graph g, int source, int dest)
{
    JRB visited = make_jrb();
    JRB cost = make_jrb();
    JRB parent = make_jrb();
    JRB node;
    Dllist q = new_dllist();
    Dllist path;
    Dllist n;
    Dllist adjs;
    int u;
    double k, i, j;
    dll_append(q, new_jval_i(source));
    jrb_traverse(node, g.vertices)
    {
        jrb_insert_int(cost, jval_i(node->key), new_jval_d(INFINITIVE_VALUE));
    }
    jrb_find_int(cost, source)->val = new_jval_d(0);
    while (1)
    {
        n = dll_min(q, cost, visited);
        if (n == NULL)
        {
            printf("No road from %d to %d\n", source, dest);
            jrb_free_tree(visited);
            jrb_free_tree(cost);
            jrb_free_tree(parent);
            free_dllist(q);
            return INFINITIVE_VALUE;
        }
        u = jval_i(dll_val(n));
        dll_delete_node(n);
        if (u == dest)
            break;
        k = jval_d(jrb_find_int(cost, u)->val); //gia tri cua cost u
        jrb_insert_int(visited, u, new_jval_i(1));
        adjs = outgoingVertices(g, u);
        dll_traverse(n, adjs)
        {
            if (jrb_find_int(visited, jval_i(dll_val(n))) == NULL)
            {
                node = jrb_find_int(cost, jval_i(dll_val(n))); //gia tri cua cost node
                i = getEdgeValue(g, u, jval_i(dll_val(n)));
                if (jval_d(node->val) > i + k)
                {
                    node->val = new_jval_d(i + k);

                    node = jrb_find_int(parent, jval_i(dll_val(n)));
                    if (node == NULL)
                    {
                        jrb_insert_int(parent, jval_i(dll_val(n)), new_jval_i(u));
                        node = jrb_find_int(parent, jval_i(dll_val(n)));
                    }
                    else
                        node->val = new_jval_i(u);
                }
                dll_append(q, dll_val(n));
            }
        }

        free_dllist(adjs);
    }
    j = jval_d(jrb_find_int(cost, dest)->val);

    printf("Shortest path length is %lf\nRoad: ", j);
    path = findPath(parent, source, dest);
    dll_traverse(n, path)
        printf("%d ", jval_i(dll_val(n)));
    printf("\n");
    jrb_free_tree(visited);
    jrb_free_tree(cost);
    jrb_free_tree(parent);
    free_dllist(q);
    free_dllist(path);
    return j;
}

Dllist dll_min(Dllist list, JRB cost, JRB visited)
{
    Dllist node;
    JRB node_g, node_m;
    double min = INFINITIVE_VALUE;
    jrb_traverse(node_g, cost)
    {
        if (jval_d(node_g->val) < min && jrb_find_int(visited, jval_i(node_g->key)) == NULL)
        {
            min = jval_d(node_g->val);
            node_m = node_g;
        }
    }
    dll_traverse(node, list)
    {
        if (jval_i(dll_val(node)) == jval_i(node_m->key))
            return node;
    }
    if (min == INFINITIVE_VALUE)
        return NULL;
}

Dllist findPath(JRB parent, int source, int dest)
{
    JRB node = jrb_find_int(parent, dest);
    Dllist path = new_dllist();
    dll_prepend(path, new_jval_i(dest));
    while (jval_i(node->val) != source)
    {
        dll_prepend(path, node->val);
        node = jrb_find_int(parent, jval_i(node->val));
    }
    dll_prepend(path, new_jval_i(source));
    return path;
}

void recursionDFS(Graph g, int source, JRB visited, Dllist stack)
{
    Dllist queue, node;
    if (jrb_find_int(visited, source) == NULL)
    {
        jrb_insert_int(visited, source, new_jval_i(1));
        queue = outgoingVertices(g, source);
        dll_rtraverse(node, queue)
            recursionDFS(g, jval_i(dll_val(node)), visited, stack);
        dll_prepend(stack, new_jval_i(source));

        free_dllist(queue);
    }
}

void topologicalSort(Graph g)
{
    JRB visited = make_jrb(), node;
    Dllist stack = new_dllist(), nodedll;
    if (!DAG(g))
    {
        printf("Do thi co chu trinh!\n");
        return;
    }
    jrb_traverse(node, g.vertices)
    {
        recursionDFS(g, jval_i(node->key), visited, stack);
    }
    dll_traverse(nodedll, stack)
    {
        printf("%d ", jval_i(dll_val(nodedll)));
    }
    printf("\n");
    free_dllist(stack);
    jrb_free_tree(visited);
}

void TPLT(Graph g)
{
    JRB visited = make_jrb(), node;
    Dllist stack = new_dllist();
    int count = 0;
    jrb_traverse(node, g.vertices)
    {
        if (jrb_find_int(visited,jval_i(node->key)) == NULL)
            count++;
        recursionDFS(g, jval_i(node->key), visited, stack);
    }
    printf("%d\n", count);
    free_dllist(stack);
    jrb_free_tree(visited);
}

void recursionDFSdao(Graph g, int source, JRB visited, Dllist stack)
{
    Dllist queue, node;
    if (jrb_find_int(visited, source) == NULL)
    {
        jrb_insert_int(visited, source, new_jval_i(1));
        queue = incomingVertices(g, source);
        dll_rtraverse(node, queue)
            recursionDFSdao(g, jval_i(dll_val(node)), visited, stack);
        dll_prepend(stack, new_jval_i(source));

        free_dllist(queue);
    }
}

void SCC(Graph g)
{
    Dllist stack, nodedll, stacktemp;
    JRB node, visited;
    stack = new_dllist();
    stacktemp = new_dllist();
    int scc = 0;
    visited = make_jrb();
    jrb_traverse(node, g.vertices)
    {
        recursionDFSdao(g, jval_i(node->key), visited, stack);
    }
    jrb_free_tree(visited);
    visited = make_jrb();

    dll_traverse(nodedll, stack)
    {
        if (jrb_find_int(visited, jval_i(dll_val(nodedll))) == NULL)
        {
            recursionDFS(g, jval_i(dll_val(nodedll)),visited, stacktemp);
            scc++;
        }
    }
    printf("%d\n",scc);
    free_dllist(stack);
    free_dllist(stacktemp);
    jrb_free_tree(visited);
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
#endif