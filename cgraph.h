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

typedef struct
{
  JRB edges;
  JRB vertices;
} Graph;

typedef Dllist Queue;
typedef Dllist Stack;

Graph createGraph();
int numVertex(Graph g);
void addVertex(Graph graph, int id, char *name);
char *getVertex(Graph graph, int id);
void addEdge(Graph graph, int v1, int v2, double weight);
double getEdgeValue(Graph graph, int v1, int v2);
int indegree(Graph graph, int v, int *output);
int outdegree(Graph graph, int v, int *output);
void dropGraph(Graph graph);
double shortestPath(Graph graph, int s, int t, int *path, int *length);
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

Graph createGraph()
{
  Graph g;
  g.edges = make_jrb();
  g.vertices = make_jrb();
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

// Dijkstra
double shortestPath(Graph g, int s, int t, int *path, int *length)
{
  // Đếm số đỉnh của đồ thị
  int numV = numVertex(g);

  if (numV == 0)
    return;

  // isSelectMin dùng để xem đỉnh đó đã được chọn chưa (chọn min)
  int *isSelectMin = (int *)calloc(numV, sizeof(int));

  // beenInPQ dùng để kiểm tra xem đỉnh đó đã ở trong hàng đợi hay chưa
  int *beenInPQ = (int *)calloc(numV, sizeof(int));
  double distance[numV], min_dist, w, total;
  int min, u;
  int previous[numV], tmp[numV];
  int n, output[numV], v;
  Dllist ptr, queue, node;

  // Dijkstra algorithm
  for (int i = 0; i < numV; i++)
    distance[i] = INFINITIVE_VALUE;
  distance[s] = 0;
  previous[s] = s;

  queue = new_dllist();
  dll_append(queue, new_jval_i(s));

  while (!dll_empty(queue))
  {
    // get u from the priority queue
    min_dist = INFINITIVE_VALUE;
    dll_traverse(ptr, queue)
    {
      u = jval_i(ptr->val);
      if (min_dist > distance[u])
      {
        min_dist = distance[u];
        min = u;
        node = ptr;
      }
    }
    dll_delete_node(node);
    u = min;
    beenInPQ[u] = 0;

    if (u == t)
      break;

    isSelectMin[u] = 1;
    n = outdegree(g, u, output);
    for (int i = 0; i < n; i++)
    {
      v = output[i];
      if (!isSelectMin[v])
      {
        w = getEdgeValue(g, u, v);
        if (distance[v] > distance[u] + w)
        {
          distance[v] = distance[u] + w;
          previous[v] = u;
        }
        if (!beenInPQ[v])
        {
          dll_append(queue, new_jval_i(v));
          beenInPQ[v] = 1;
        }
      }
    }
  }

  free(isSelectMin);
  free(beenInPQ);
  free(queue);

  total = distance[t];
  if (total != INFINITIVE_VALUE)
  {
    tmp[0] = t;
    n = 1;
    while (t != s)
    {
      t = previous[t];
      tmp[n++] = t;
    }
    for (int i = n - 1; i >= 0; i--)
      path[n - i - 1] = tmp[i];
    *length = n;
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
    return jval_i(dll_first(S)->val);
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
    return;

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
  
  if(!DAG(g)){
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
  // Đếm số đỉnh của đồ thị
  int numV = numVertex(g);

  if (numV == 0)
    return;

  // Khai báo các biến cần thiết
  int *visited = (int *)calloc(numV, sizeof(int));
  double distance[numV], min_dist, w, lastID;
  int min, u;
  int tmp[numV], previous[numV], path[numV];
  int n, output[numV], v, t;
  Dllist ptr, queue, node;

  // Khỏi tạo chi phí ban đầu đến các đỉnh là vô cùng
  for (int i = 0; i < numV; i++)
    distance[i] = INFINITIVE_VALUE;

  // Chọn một đỉnh đâu tiên trong cây JRB làm đỉnh bắt đầu
  int s = jval_i(jrb_first(g.edges)->key);

  // Gán chi phí đến đỉnh s là 0
  distance[s] = 0;

  // Gán cha của s là s
  previous[s] = s;

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
      if (min_dist > distance[u])
      {
        min_dist = distance[u];
        min = u;
        node = ptr;
      }
    }

    // Make empty queue
    while (!dll_empty(queue))
      dll_delete_node(dll_first(queue));

    u = min;

    // lastID dùng để lưu đỉnh cuối cùng thao tác
    lastID = u;

    // Đánh đấu là đã thăm (hay đã cho vào tập X)
    visited[u] = 1;

    // Lấy các đỉnh kề của đỉnh u vào mảng output
    n = outdegree(g, u, output);

    // Lặp qua các đỉnh kê của u
    for (int i = 0; i < n; i++)
    {
      // Lấy ra ID của đỉnh kề thứ i của u
      v = output[i];

      // Nếu chưa đưa vào tập X
      if (!visited[v])
      {
        // Lấy giá trị cạnh u->v
        w = getEdgeValue(g, u, v);
        // Gán lại chi phí đường đi tư u -> v
        distance[v] = w;

        // Gán cha của v là u
        previous[v] = u;

        // Thêm v vào hàng đợi
        dll_append(queue, new_jval_i(v));
      }
    }
  }

  // Free
  free(visited);
  free(queue);

  // Truy ngược lại previous để tìm MST
  t = lastID;
  tmp[0] = t;
  n = 1;
  while (t != s)
  {
    t = previous[t];
    tmp[n++] = t;
  }

  // Chuyển cuối lên đâu
  for (int i = n - 1; i >= 0; i--)
    path[n - i - 1] = tmp[i];

  // Export to dot file
  FILE *fp = fopen("MST.dot", "w");
  fprintf(fp, "graph MST\n{\n");

  // printed đc dùng để in ra file .dot
  Graph printed = createGraph();
  JRB Node, temp, ptr2;
  int v1, v2;
  jrb_traverse(Node, g.vertices)
      fprintf(fp, "\t%d\n", jval_i(Node->key));

  // IN RA MST (Cây khung nhỏ nhất)
  for (int i = 0; i < n - 1; i++)
  {
    v1 = path[i];
    v2 = path[i + 1];
    w = getEdgeValue(g, v1, v2);

    fprintf(fp, "\t%d -- %d [color=\"#ff0000\", label=\"%g\", penwidth=3]\n", v1, v2, w);
    addEdge(printed, v1, v2, 0);
    addEdge(printed, v2, v1, 0);
  }

  // In ra các cạnh còn lại
  jrb_traverse(Node, g.edges)
  {
    temp = (JRB)jval_v(Node->val);
    jrb_traverse(ptr2, temp)
    {
      v1 = jval_i(Node->key);
      v2 = jval_i(ptr2->key);
      if (!hasEdge(printed, v1, v2) && !hasEdge(printed, v1, v2))
      {
        w = jval_d(ptr2->val);
        addEdge(printed, v1, v2, 0);
        addEdge(printed, v2, v1, 0);
        fprintf(fp, "\t%d -- %d [label=\"%g\"]\n", v1, v2, w);
      }
    }
  }

  fprintf(fp, "}");
  fclose(fp);

  dropGraph(printed);
}

//To mau do thi
//Input: Do thi g
//Output: file .dot
void Coloring(Graph g)
{
  int m, n = 0, a, b, adj[500], nadj;
  int i, *A, *B, j;
  JRB ptr, node, tmp;
  FILE *fp;

  //Đếm số đỉnh
  jrb_traverse(ptr, g.vertices)
      n++;

  A = (int *)calloc(n, sizeof(int));
  B = (int *)malloc(n * sizeof(int));

  //Tìm số bậc của đỉnh, lưu vào A
  for (i = 0; i < n; ++i)
  {
    nadj = outdegree(g, i, adj);
    A[i] += nadj;
    for (j = 0; j < nadj; ++j)
      if (hasEdge(g, adj[j], i) == 1)
        A[i]--;
    A[i] += indegree(g, i, adj);
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
      if ((hasEdge(g, A[i], A[b]) || hasEdge(g, A[b], A[i])) && B[A[i]] == B[A[b]])
      {
        //Nếu 2 đỉnh của 1 cạnh có cùng màu
        B[A[i]]++;
        b = -1;
        if (B[A[i]] == a)
          a++;
      }
  }

  fp = fopen("dothitomau.dot", "w");
  fprintf(fp, "digraph dothi\n{\n");

  //Tô màu các đỉnh
  for (i = 0; i < n; ++i)
    fprintf(fp, "%d [fillcolor=\".%d .3 1.0\", style=filled];\n", i, B[i]);

  //Vẽ các cạnh
  jrb_traverse(node, g.edges)
  {
    tmp = (JRB)jval_v(node->val);
    jrb_traverse(ptr, tmp)
        fprintf(fp, "%d -> %d\n", jval_i(node->key), jval_i(ptr->key));
  }
  fprintf(fp, "}");

  free(A);
  free(B);
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
  Graph gReverse = createGraph();
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

#endif