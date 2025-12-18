#include "ArcadiaEngine.h"
#include <algorithm>
#include <queue>
#include <cstring>
#include <numeric>
#include <climits>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <string>
#include <iostream>
#include <map>
#include <set>

#define C first
#define U second.first
#define V second.second


using namespace std;

typedef pair<long long, pair<int,int>> edge;
// =========================================================
// PART A: DATA STRUCTURES (Concrete Implementations)
// =========================================================

// --- 1. PlayerTable (Double Hashing) ---

class ConcretePlayerTable : public PlayerTable {
private:

    static const int TABLE_SIZE = 101;

    struct Player{
        int playerID;
        string name;
        bool occupied;

        Player() : playerID(-1) , name("") , occupied(false) {}
    };

    Player table[TABLE_SIZE];

    // Primary Hash Function
    int hash1(int key){
        return key % TABLE_SIZE;
    }

    // Secondary Hash Function
    int hash2(int key){
        return 1 + (key % (TABLE_SIZE - 1));
    }

public:
    ConcretePlayerTable() {

    }

    void insert(int playerID, string name) override {

        int idx1 = hash1(playerID);
        int idx2 = hash2(playerID);

        for(int i = 0 ; i < TABLE_SIZE; i++){
            int idx = (idx1 + i * idx2) % TABLE_SIZE; // Formula of Double Hashing
            if(!table[idx].occupied){
                table[idx].playerID = playerID;
                table[idx].name = name;
                table[idx].occupied = true;
                return;
            }
        }
        cout << "The Table is Full.\n";
    }

    string search(int playerID) override {

        int idx1 = hash1(playerID);
        int idx2 = hash2(playerID);

        for(int i = 0;i < TABLE_SIZE ; i++){
            int idx = (idx1 + i * idx2) % TABLE_SIZE; // Formula of Double Hashing

            if(table[idx].playerID == playerID){
                return table[idx].name; // found
            }
        }

        return "Not Found!\n"; // not found
    }
};

// --- 2. Leaderboard (Skip List) ---

class ConcreteLeaderboard : public Leaderboard {
private:
    static const int MAX_LEVEL = 16;
    const float prob = 0.5;
    int currentLevel;

    struct Node {
        int playerID;
        int score;
        vector<Node*> forward;

        Node(int id, int sc, int level) {
            playerID = id;
            score = sc;
            forward.resize(level + 1, nullptr);
        }
    };

    Node* sentinal;

    int randomLevel() {
        int level = 0;
        while (((float)rand() / RAND_MAX) < prob && level < MAX_LEVEL)
            level++;
        return level;
    }

    // ORDER: score DESC, playerID ASC
    bool comesBefore(Node* node, int score, int id) {
        if (node->score != score)
            return node->score > score;
        return node->playerID < id;
    }

public:
    ConcreteLeaderboard() {
        currentLevel = 0;
        sentinal = new Node(-1, INT_MAX, MAX_LEVEL);
    }

    void removePlayer(int playerID) override {
        Node* current = sentinal;
        vector<Node*> update(MAX_LEVEL + 1);

        for (int i = currentLevel; i >= 0; i--) {
            while (current->forward[i] &&
                   current->forward[i]->playerID != playerID) {
                current = current->forward[i];
            }
            update[i] = current;
        }

        current = current->forward[0];

        if (current && current->playerID == playerID) {
            for (int i = 0; i <= currentLevel; i++) {
                if (update[i]->forward[i] != current) break;
                update[i]->forward[i] = current->forward[i];
            }
            delete current;

            while (currentLevel > 0 &&
                   sentinal->forward[currentLevel] == nullptr)
                currentLevel--;
        }
    }

    void addScore(int playerID, int score) override {
        Node* current = sentinal;
        vector<Node*> update(MAX_LEVEL + 1);

        for (int i = currentLevel; i >= 0; i--) {
            while (current->forward[i] &&
                   comesBefore(current->forward[i], score, playerID)) {
                current = current->forward[i];
            }
            update[i] = current;
        }

        current = current->forward[0];

        // Player exists → update score
        if (current && current->playerID == playerID) {
            int newScore = current->score + score;
            removePlayer(playerID);
            addScore(playerID, newScore);
            return;
        }

        int level = randomLevel();
        if (level > currentLevel) {
            for (int i = currentLevel + 1; i <= level; i++)
                update[i] = sentinal;
            currentLevel = level;
        }

        Node* newNode = new Node(playerID, score, level);
        for (int i = 0; i <= level; i++) {
            newNode->forward[i] = update[i]->forward[i];
            update[i]->forward[i] = newNode;
        }
    }

    vector<int> getTopN(int n) override {
        vector<int> result;
        Node* current = sentinal->forward[0];
        while (current && n--) {
            result.push_back(current->playerID);
            current = current->forward[0];
        }
        return result;
    }
};

// --- 3. AuctionTree (Red-Black Tree) ---

enum Color { RED, BLACK };

// ---------------- NODE ----------------
class Node {
public:
    int itemID;
    int price;
    Color color;
    Node* left;
    Node* right;
    Node* parent;

    Node(int id, int price) {
        this->itemID = id;
        this->price = price;
        color = RED;
        left = right = parent = nullptr;
    }
};

// ---------------- RED-BLACK TREE ----------------
class RBTree {
private:
    Node* root;
    Node* NIL;

    void leftRotate(Node* x) {
        Node* y = x->right;
        x->right = y->left;
        if (y->left != NIL) y->left->parent = x;
        y->parent = x->parent;
        if (x->parent == nullptr) root = y;
        else if (x == x->parent->left) x->parent->left = y;
        else x->parent->right = y;
        y->left = x;
        x->parent = y;
    }

    void rightRotate(Node* x) {
        Node* y = x->left;
        x->left = y->right;
        if (y->right != NIL) y->right->parent = x;
        y->parent = x->parent;
        if (x->parent == nullptr) root = y;
        else if (x == x->parent->right) x->parent->right = y;
        else x->parent->left = y;
        y->right = x;
        x->parent = y;
    }

    void bstInsert(Node* z) {
        Node* y = nullptr;
        Node* x = root;
        while (x != NIL) {
            y = x;
            if (z->itemID < x->itemID) x = x->left;
            else x = x->right;
        }
        z->parent = y;
        if (y == nullptr) root = z;
        else if (z->itemID < y->itemID) y->left = z;
        else y->right = z;
        z->left = z->right = NIL;
    }

    void fixInsert(Node* z) {
        while (z->parent != nullptr && z->parent->color == RED) {
            Node* parent = z->parent;
            Node* grandparent = parent->parent;

            if (parent == grandparent->left) {
                Node* uncle = grandparent->right;
                if (uncle->color == RED) {
                    parent->color = BLACK;
                    uncle->color = BLACK;
                    grandparent->color = RED;
                    z = grandparent;
                }
                else {
                    if (z == parent->right) {
                        z = parent;
                        leftRotate(z);
                    }
                    parent->color = BLACK;
                    grandparent->color = RED;
                    rightRotate(grandparent);
                }
            }
            else {
                Node* uncle = grandparent->left;
                if (uncle->color == RED) {
                    parent->color = BLACK;
                    uncle->color = BLACK;
                    grandparent->color = RED;
                    z = grandparent;
                }
                else {
                    if (z == parent->left) {
                        z = parent;
                        rightRotate(z);
                    }
                    parent->color = BLACK;
                    grandparent->color = RED;
                    leftRotate(grandparent);
                }
            }
        }
        root->color = BLACK;
    }

    void transplant(Node* u, Node* v) {
        if (u->parent == nullptr) root = v;
        else if (u == u->parent->left) u->parent->left = v;
        else u->parent->right = v;
        v->parent = u->parent;
    }

    Node* minimum(Node* node) {
        while (node->left != NIL) node = node->left;
        return node;
    }

    void bstDelete(Node* z) {
        Node* y = z;
        Node* x;
        Color yOriginalColor = y->color;

        if (z->left == NIL) {
            x = z->right;
            transplant(z, z->right);
        }
        else if (z->right == NIL) {
            x = z->left;
            transplant(z, z->left);
        }
        else {
            y = minimum(z->right);
            yOriginalColor = y->color;
            x = y->right;
            if (y->parent == z) x->parent = y;
            else {
                transplant(y, y->right);
                y->right = z->right;
                y->right->parent = y;
            }
            transplant(z, y);
            y->left = z->left;
            y->left->parent = y;
            y->color = z->color;
        }
        delete z;
        if (yOriginalColor == BLACK) fixDelete(x);
    }

    void fixDelete(Node* x) {
        while (x != root && x->color == BLACK) {
            if (x == x->parent->left) {
                Node* w = x->parent->right;
                if (w->color == RED) {
                    w->color = BLACK;
                    x->parent->color = RED;
                    leftRotate(x->parent);
                    w = x->parent->right;
                }
                if (w->left->color == BLACK && w->right->color == BLACK) {
                    w->color = RED;
                    x = x->parent;
                }
                else {
                    if (w->right->color == BLACK) {
                        w->left->color = BLACK;
                        w->color = RED;
                        rightRotate(w);
                        w = x->parent->right;
                    }
                    w->color = x->parent->color;
                    x->parent->color = BLACK;
                    w->right->color = BLACK;
                    leftRotate(x->parent);
                    x = root;
                }
            }
            else {
                Node* w = x->parent->left;
                if (w->color == RED) {
                    w->color = BLACK;
                    x->parent->color = RED;
                    rightRotate(x->parent);
                    w = x->parent->left;
                }
                if (w->left->color == BLACK && w->right->color == BLACK) {
                    w->color = RED;
                    x = x->parent;
                }
                else {
                    if (w->left->color == BLACK) {
                        w->right->color = BLACK;
                        w->color = RED;
                        leftRotate(w);
                        w = x->parent->left;
                    }
                    w->color = x->parent->color;
                    x->parent->color = BLACK;
                    w->left->color = BLACK;
                    rightRotate(x->parent);
                    x = root;
                }
            }
        }
        x->color = BLACK;
    }

public:
    RBTree() {
        NIL = new Node(0, 0);
        NIL->color = BLACK;
        NIL->left = NIL->right = NIL->parent = nullptr;
        root = NIL;
    }

    void insert(int itemID, int price) {
        Node* z = new Node(itemID, price);
        z->left = z->right = NIL;
        bstInsert(z);
        fixInsert(z);
    }

    void remove(int itemID) {
        Node* current = root;
        while (current != NIL) {
            if (itemID == current->itemID) break;
            else if (itemID < current->itemID) current = current->left;
            else current = current->right;
        }
        if (current != NIL) bstDelete(current);
    }
};

class ConcreteAuctionTree : public AuctionTree {
private:
    RBTree tree;  // Using RBTree internally
public:
    ConcreteAuctionTree() {}

    void insertItem(int itemID, int price) override {
        tree.insert(itemID, price);
    }

    void deleteItem(int itemID) override {
        tree.remove(itemID);
    }
};

// =========================================================
// PART B: INVENTORY SYSTEM (Dynamic Programming)
// =========================================================

int InventorySystem::optimizeLootSplit(int n, vector<int>& coins) {
    int total = 0;
    for (int i : coins) total += i;
    int half = total / 2;
    vector<vector<bool>> table(n, vector<bool>(half + 1, false));

    for (int i = 0; i < n; i++) table[i][0] = true;
    if (coins[0] <= half) table[0][coins[0]] = true;

    for (int i = 1; i < n; i++) {
        for (int j = 1; j <= half; j++) {
            table[i][j] = table[i - 1][j];
            if (j >= coins[i])
                table[i][j] = table[i][j] || table[i - 1][j - coins[i]];
        }
    }

    int firstPartitionSum = 0;
    for (int j = half; j >= 0; j--) {
        if (table[n - 1][j]) {
            firstPartitionSum = j;
            break;
        }
    }

    return total - 2 * firstPartitionSum;
}

int InventorySystem::maximizeCarryValue(int capacity, vector<pair<int, int>>& items) {

    // n is the number of items
    // the first index in each item is the weight if the item
    // the second index in each item is the value of the item

    int n = items.size();
    int weights[n+1];
    int values[n+1];

    //filling the weights and values array with the passed data
    weights[0] = 0;
    values[0] = 0;
    for (int i = 1; i <= n; ++i) {
        weights[i] = items[i - 1].first;
        values[i] = items[i - 1].second;
    }


    //filling each cell in the matrix with the prober data
    int matrix[n+1][capacity+1];

    // initialize first row and column as zeros

    for(int i = 0; i <= n;i++){
        matrix[0][i] = 0;
    }
    for(int i = 0; i <= n;i++){
        matrix[i][0] = 0;
    }

    for (int i = 1; i <= n; ++i) {
        for (int j = 1; j <= capacity; ++j) {
            if (i==1 || j==1)
                matrix[i][j] = 0;
            else if (weights[i] > j)
                matrix[i][j] = matrix[i-1][j];
            else
                matrix[i][j] = max(matrix[i-1][j] , matrix[i-1][j-weights[i]] + values[i]);
        }
    }

    return matrix[n][capacity];

}

long long InventorySystem::countStringPossibilities(string s) {
    // TODO: Implement string decoding DP
    // Rules: "uu" can be decoded as "w" or "uu"
    //        "nn" can be decoded as "m" or "nn"
    // Count total possible decodings
    return 0;
}

// =========================================================
// PART C: WORLD NAVIGATOR (Graphs)
// =========================================================

bool WorldNavigator::pathExists(int n, vector<vector<int>>& edges, int source, int dest) {
    // BFS
    // define the Adjacency list
    vector<vector<int>> adj(n);

    for(auto &e: edges){
        int u = e[0],v = e[1];
        // bidirectional Graph
        adj[u].push_back(v);
        adj[v].push_back(u);
    }

    vector <int> visited(n,false);
    visited[source] = true;
    queue<int> q;
    q.push(source);

    while(!q.empty()){ // nodes not visited
        int u = q.front();
        q.pop();

        if(u == dest){
            return true;
        }

        for(int v : adj[u]){ // neighbours of u
            if(!visited[v]){
                visited[v] = true;
                q.push(v);
            }
        }
    }

    return false;
}

// =======================
// Union-Find (Disjoint Set)
// =======================
vector<int> parent, rankSet;

void makeSet(int n) {
    parent.resize(n);
    rankSet.assign(n, 0);
    for (int i = 0; i < n; i++)
        parent[i] = i;
}

int findSet(int u) {
    if (parent[u] == u)
        return u;
    return parent[u] = findSet(parent[u]); // path compression
}

void unite(int u, int v) {
    u = findSet(u);
    v = findSet(v);
    if (u == v) return;

    if (rankSet[u] < rankSet[v])
        parent[u] = v;
    else if (rankSet[u] > rankSet[v])
        parent[v] = u;
    else {
        parent[v] = u;
        rankSet[u]++;
    }
}

bool sameSet(int u, int v) {
    return findSet(u) == findSet(v);
}


long long WorldNavigator::minBribeCost(int n, int m, long long goldRate, long long silverRate,
                                       vector<vector<int>>& roadData) {


    vector<edge> edgeList;

    for (int i = 0; i < m; i++) {
        int u = roadData[i][0];
        int v = roadData[i][1];
        int gold = roadData[i][2];
        int silver = roadData[i][3];


        long long cost = (long long)gold * goldRate + (long long)silver * silverRate;

        edgeList.push_back({cost, {u, v}});
    }

    // Kruskal Algorithm using DSU
    sort(edgeList.begin(), edgeList.end());
    makeSet(n);

    long long mstCost = 0;
    int edgesUsed = 0;

    for (auto &e : edgeList) {
        int u = e.second.first;
        int v = e.second.second;
        long long cost = e.first;

        if (!sameSet(u, v)) {
            unite(u, v);
            mstCost += cost;
            edgesUsed++;
        }

        if (edgesUsed == n - 1)
            break;
    }

    if (edgesUsed != n - 1) return -1;

    return mstCost;
}



string WorldNavigator::sumMinDistancesBinary(int n, vector<vector<int>>& roads) {

    const long long inf = 0x3f3f3f3f;
    long long dist[n][n];
    memset(dist,inf,sizeof dist);

    for(int i = 0 ; i < n ;i++){
        dist[i][i] = 0;

    }

    // store roads

    for(auto &r : roads){
        int u = r[0];
        int v = r[1];
        long long w = r[2];

        dist[u][v] = min(dist[u][v], w);
        dist[v][u] = min(dist[v][u], w);
    }


    // Floyd Warshall Algorithm O(V^3)
    for (int k = 0; k < n; k++) {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (dist[i][k] < inf && dist[k][j] < inf) {
                    dist[i][j] = min(dist[i][j], dist[i][k] + dist[k][j]);
                }
            }
        }
    }


    long long sum = 0;
    for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
            if (dist[i][j] < inf)
                sum += dist[i][j];
        }
    }

    if (sum == 0) return "0";

    string binary = "";
    while (sum > 0) {
        binary = char('0' + (sum % 2)) + binary;
        sum /= 2;
    }

    return binary;

}

// =========================================================
// PART D: SERVER KERNEL (Greedy)
// =========================================================

int ServerKernel::minIntervals(vector<char>& tasks, int n) {

    priority_queue<int> Maxheap;
    queue<pair<int, int>> waiting;

    int freq[26] = { 0 };
    int timer = 0;


    for (char c : tasks) {
        freq[c - 'A']++;
    }


    for (int i = 0; i < 26; i++) {
        if (freq[i] > 0)
            Maxheap.push(freq[i]);
    }

    while (!Maxheap.empty() || !waiting.empty()) {

        if (!waiting.empty() && waiting.front().second == timer) {
            Maxheap.push(waiting.front().first);
            waiting.pop();
        }
        if (!Maxheap.empty()) {
            int temp = Maxheap.top();
            Maxheap.pop();
            temp--;
            if (temp > 0) {
                waiting.push({ temp, timer + n + 1 });
            }
        }
        timer++;
    }

    return timer;
}

// =========================================================
// FACTORY FUNCTIONS (Required for Testing)
// =========================================================

extern "C" {
    PlayerTable* createPlayerTable() {
        return new ConcretePlayerTable();
    }

    Leaderboard* createLeaderboard() {
        return new ConcreteLeaderboard();
    }

    AuctionTree* createAuctionTree() {
        return new ConcreteAuctionTree();
    }
}
