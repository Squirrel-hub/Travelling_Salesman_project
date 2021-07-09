/*
Name - Samyak Chakrabarty , Roll No. - 19EC10084
Name - Kaizer Rahman , Roll No. - 19IE10044
Topic - Algorithms 1 Term Project
*/

#include <iostream>
#include <math.h>
using namespace std;

typedef struct POINT
{
    int id;
    int x;
    int y;
    int q;
} point;

typedef struct GRAPH
{
    int v;
    double **adj;
} graph;

struct AVLTreeNode
{
    struct AVLTreeNode *left;
    int weight;
    struct AVLTreeNode *right;
    point *P;
    int length;
    int height;
    int cap;
};

// Global Variable for Array of Points that constitute the tour
point *tour;
// Counter for previous variable
int t;
// Distance covered by the truck
double dis = 0;
// Pointer to previous point crossed
point Prev;

// Eucledian Distance Norm
double distance(point A, point B)
{
    double D = pow((A.x - B.x), 2) + pow((A.y - B.y), 2);
    double dist = pow(D, 0.5);
    return dist;
}

graph *make_Graph(int V, point *P)
{
    graph *G = (graph *)malloc(sizeof(graph));
    G->v = V + 1;

    int maxx_dim = V;
    G->adj = new double *[maxx_dim + 1];
    for (int i = 0; i < maxx_dim + 1; i++)
        G->adj[i] = new double[maxx_dim + 1];
    for (int i = 0; i < V + 1; i++)
        for (int j = 0; j < V + 1; j++)
            G->adj[i][j] = distance(P[i], P[j]);

    return G;
}

int Max(int a,int b)
{
    if(a>b)
        return a;
    else
        return b;
}

/* modified AVL Tree Functions */
int Height(struct AVLTreeNode *root)
{
    if(!root)
        return -1;
    else
        return root->height;
}

struct AVLTreeNode *SingleRotateLeft(struct AVLTreeNode *X)
{
    struct AVLTreeNode *W=X->left;
    X->left = W->right;
    W->right=X;
    X->height=Max(Height(X->left),Height(X->right))+1;
    W->height=Max(Height(W->left),X->height)+1;
    return W;
}

struct AVLTreeNode *SingleRotateRight(struct AVLTreeNode *W)
{
    struct AVLTreeNode *X=W->right;
    W->right=X->left;
    X->left=W;
    W->height=Max(Height(W->right),Height(W->left))+1;
    X->height=Max(Height(X->right),W->height)+1;
    return X;
}

struct AVLTreeNode *DoubleRotatewithLeft(struct AVLTreeNode *Z)
{
    Z->left=SingleRotateRight(Z->left);
    return SingleRotateLeft(Z);
}

struct AVLTreeNode *DoubleRotatewithRight(struct AVLTreeNode *Z)
{
    Z->right=SingleRotateLeft(Z->right);
    return SingleRotateRight(Z);
}

struct AVLTreeNode *Insert(struct AVLTreeNode *root,struct AVLTreeNode *parent,point P,point Source,int Max_capacity,int n)
{
    if(!root)
    {
        root=(struct AVLTreeNode *)malloc(sizeof(struct AVLTreeNode));
        if(!root)
        {
            cout<<"\nMemory Error";
            return NULL;
        }
        else
        {
            root->weight=P.q;
            root->height=0;
            root->left=root->right=NULL;
            root->P=(point *)malloc(n*sizeof(point));
            root->P[0]=Source;
            root->P[1]=P;
            root->length=2;
            root->cap=n;
        }
    }
    else if(root->weight+P.q<=Max_capacity)
    {
        root->P[root->length]=P;
        root->length++;
        root->weight+=P.q;
    }
    else
    {
        root->left=Insert(root->left,root,P,Source,Max_capacity,n);
        if((Height(root->left)-Height(root->right))==2)
        {
            if(P.q<root->left->weight)
                root=SingleRotateLeft(root);
            else
                root=DoubleRotatewithLeft(root);
        }
    }
    root->height=Max(Height(root->left),Height(root->right))+1;
    return root;
}

// Minimum edge detection for Prim's algorithm
int minKey(int key[], bool mstSet[], int V)
{
    int min = INT_MAX, min_index;

    for (int v = 0; v < V; v++)
        if (mstSet[v] == false && key[v] < min)
            min = key[v], min_index = v;

    return min_index;
}

// Pre-order walk
void preorder(int parent[], graph *G, point *P, int root, int visit[])
{
    double min_dist = 0;
    tour[t++] = P[root];

    visit[root] = 1;
    for (int i = 1; i < G->v; i++)
    {

        if (parent[i] == root && visit[i] == 0)
            preorder(parent, G, P, i, visit);

        if (i == root && visit[parent[i]] == 0)
            preorder(parent, G, P, parent[i], visit);
    }
}

// In-order walk for final printing
void inorder(struct AVLTreeNode *T)
{
    if(T == NULL) return;

    for(int i=0;i<T->length;i++)
    {
        cout << T->P[i].id << " ";
        dis += distance(Prev,T->P[i]);
        Prev = T->P[i];
    }

    inorder(T->left);
    inorder(T->right);

}

// MST based 2-approximate TSP
void approx_TSP(graph *G, point *P)
{
    int V = G->v;
    int parent[V];
    int key[V];
    bool mstSet[V];
    for (int i = 0; i < V; i++)
        key[i] = INT_MAX, mstSet[i] = false;
    key[0] = 0;
    parent[0] = -1;
    for (int count = 0; count < V - 1; count++)
    {
        int u = minKey(key, mstSet, V);
        mstSet[u] = true;
        for (int v = 0; v < V; v++)
            if (G->adj[u][v] && mstSet[v] == false && G->adj[u][v] < key[v])
                parent[v] = u, key[v] = G->adj[u][v];
    }
    int visit[V];
    for (int i = 0; i < V; i++)
        visit[i] = 0;

    tour = (point *)calloc(V + 1, sizeof(point));
    t = 0;
    preorder(parent, G, P, 0, visit);

}

// 2-approximate bin packing using First fit allocation
AVLTreeNode *FFA(point *P, int C, int V)
{
    struct AVLTreeNode *root=NULL;
    for(int i=1;i<V+1;i++)
        root=Insert(root,NULL,P[tour[i].id],P[tour[0].id],C,V+1);

    Prev = root->P[0]; 
    inorder(root);
    cout << P[0].id << " ";
    dis += distance(Prev,P[0]);

    return root;        
}

// Driver Code
int main()
{
    int N, M, C;
    cout << "Write N and C: ";
    cin >> N >> C;
    point *P = (point *)malloc((N + 1) * sizeof(point));
    cout << "Write Pickup Point Coordinate: ";
    cin >> P[0].x >> P[0].y;
    P[0].id = 0;
    cout << "Write Delivery Locations and Weights:\n";
    for (int j = 1; j < N + 1; j++)
    {
        P[j].id = j;
        cin >> P[j].x >> P[j].y >> P[j].q;
    }

    graph *G = make_Graph(N, P);
    approx_TSP(G, P);
    cout << "Optimal Tour: \n";
    AVLTreeNode *root= FFA(P,C,N);
    
    cout << "Distance: " << dis;
    
}