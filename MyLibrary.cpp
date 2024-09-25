#include <bits/stdc++.h>
using namespace std;

// Vector implementation
class DynamicArray
{
private:
    int size_;
    int max_;
    int *arrayholder;

public:
    DynamicArray(int n = 10)
    {
        this->size_ = 0;
        this->max_ = n;
        this->arrayholder = new int[n];
    }

    ~DynamicArray()
    {
        delete[] this->arrayholder;
    }

    int size()
    {
        return this->size_;
    }

    int &operator[](int i)
    {
        assert(i < this->size_);
        return this->arrayholder[i];
    }

    int push(int n)
    {
        if (this->max_ < this->size_ + 1)
        {
            this->max_ *= 2;
            int *tmp_ = new int[this->max_];
            for (int i = 0; i < this->size_; i++)
            {
                tmp_[i] = this->arrayholder[i];
            }
            delete[] this->arrayholder;
            this->arrayholder = tmp_;
            this->arrayholder[this->size_] = n;
            this->size_ += 1;
        }
        else
        {
            this->arrayholder[this->size_] = n;
            this->size_ += 1;
        }
    }
};

// Linkedlist implementation
template <typename T>
class Node
{
public:
    T value;
    Node *next;
    Node *prev;

    Node(T value)
    {
        this->value = value;
    }
};

template <typename T>
class LinkedList
{
private:
    int size_;
    Node<T> *head = NULL;
    Node<T> *tail = NULL;
    Node<T> *itr = NULL;

public:
    LinkedList()
    {
        this->size_ = 0;
    }

    void append(T value)
    {
        if (this->head == NULL)
        {
            this->head = new Node<T>(value);
            this->tail = this->head;
        }
        else
        {
            this->tail->next = new Node<T>(value);
            this->tail->next->prev = this->tail;
            this->tail = this->tail->next;
        }
        this->size_ += 1;
    }

    void prepend(T value)
    {
        if (this->head == NULL)
        {
            this->head = new Node<T>(value);
            this->tail = this->head;
        }
        else
        {
            this->head->prev = new Node<T>(value);
            this->head->prev->next = this->head;
            this->head = this->head->prev;
        }
        this->size_ += 1;
    }

    Node<T> *iterate()
    {
        if (this->itr == NULL)
        {
            this->itr = this->head;
        }
        else
        {
            this->itr = this->itr->next;
        }
        return this->itr;
    }

    T ptr()
    {
        return this->itr->value;
    }

    void resetIterator()
    {
        this->itr = NULL;
    }
};

// Hashset implementation
struct bucket
{
    vector<int> v;

    bool has(int val)
    {
        return find(v.begin(), v.end(), val) != v.end();
    }

    void push(int val)
    {
        v.push_back(val);
    }

    void erase(int val)
    {
        v.erase(find(v.begin(), v.end(), val));
    }
};

class Hashset
{
private:
    bucket *container;
    int bucket_count = 0;
    int count = 0;

    int hash(int val)
    {
        int key = val % bucket_count;
        return key;
    }

public:
    Hashset(int size = 100)
    {
        this->bucket_count = size;
        this->container = new bucket[size];
    }

    void insert(int val)
    {
        int key = hash(val);
        if (this->container[key].has(val))
        {
            return;
        }
        this->container[key].push(val);
        count++;
    }

    void remove(int val)
    {
        int key = hash(val);
        if (this->container[key].has(val))
        {
            this->container[key].erase(val);
            count--;
        }
    }

    bool find(int val)
    {
        int key = hash(val);
        return this->container[key].has(val);
    }

    int size()
    {
        return count;
    }
};

// Hashmap implementation
struct bucket
{
    int v[1000];

    bucket()
    {
        memset(v, 0, sizeof(v));
    }

    bool has(int val)
    {
        return v[val] != 0;
    }

    void push(int val)
    {
        v[val]++;
    }

    void erase(int val)
    {
        v[val] = 0;
    }

    int value(int val)
    {
        return v[val];
    }
};

class Hashmap
{
private:
    bucket *container;
    int bucket_count = 0;
    int count = 0;

    int hash(int val)
    {
        int key = val % bucket_count;
        return key;
    }

public:
    Hashmap(int size = 100)
    {
        this->bucket_count = size;
        this->container = new bucket[size];
    }

    void insert(int val)
    {
        int key = hash(val);
        if (this->container[key].has(val))
        {
            this->container[key].push(val);
            return;
        }
        this->container[key].push(val);
        count++;
    }

    void remove(int val)
    {
        int key = hash(val);
        if (this->container[key].has(val))
        {
            this->container[key].erase(val);
            count--;
        }
    }

    bool find(int val)
    {
        int key = hash(val);
        return this->container[key].has(val);
    }

    int size()
    {
        return count;
    }

    int operator[](int val)
    {
        int key = hash(val);
        if (this->container[key].has(val) == 0)
        {
            this->container[key].push(val);
        }
        else
        {
            this->container[key].push(val);
        }
        return this->container[key].value(val);
    }
};

// BST implementation
struct Node
{
    int data;
    struct Node *left = NULL;
    struct Node *right = NULL;

    Node(int val = 0)
    {
        this->data = val;
    }
};

class BST
{
private:
    Node *root;

    void inorder(Node *node)
    {
        if (!node)
        {
            return;
        }
        inorder(node->left);
        cout << node->data << " ";
        inorder(node->right);
    }

    Node *insert(Node *&root, int val)
    {
        if (!root)
        {
            root = new Node(val);
        }
        else
        {
            if (val < root->data)
            {
                root->left = insert(root->left, val);
            }
            else if (val > root->data)
            {
                root->right = insert(root->right, val);
            }
        }
        return root;
    }

    bool find(Node *root, int val)
    {
        if (!root)
        {
            return 0;
        }
        else if (root->data == val)
        {
            return 1;
        }
        else if (root->data > val)
        {
            find(root->left, val);
        }
        else if (root->data < val)
        {
            find(root->right, val);
        }
    }

    Node *findMin(Node *node)
    {
        if (!node)
        {
            return NULL;
        }
        else if (node->left == NULL)
        {
            return node;
        }
        else
        {
            return findMin(node->left);
        }
    }

    Node *erase(Node *&root, int val)
    {
        Node *temp;
        if (root == NULL)
        {
            return NULL;
        }
        else if (root->data > val)
        {
            root->left = erase(root->left, val);
        }
        else if (root->data < val)
        {
            root->right = erase(root->right, val);
        }
        else if (root->left && root->right)
        {
            temp = findMin(root->right);
            root->data = temp->data;
            root->right = erase(root->right, root->data);
        }
        else
        {
            temp = root;
            if (root->left == NULL)
            {
                root = root->right;
            }
            else if (root->right == NULL)
            {
                root = root->left;
            }
            delete temp;
        }
        return root;
    }

public:
    BST()
    {
        root = NULL;
    }
    void insert(int val)
    {
        insert(root, val);
    }

    void traverse()
    {
        inorder(root);
        cout << "\n";
    }

    void search(int val)
    {
        cout << find(root, val) << "\n";
    }

    void erase(int val)
    {
        erase(root, val);
    }
};

// Heap (max) implementation
class Heap
{
private:
    vector<int> v;
    int ptr = 0;

public:
    Heap()
    {
        v.push_back(-1);
    }

    void insert(int val)
    {
        ptr += 1;
        v.push_back(val);
        int i = ptr;
        while (i > 1)
        {
            int parent = i / 2;
            if (v[parent] < v[i])
            {
                swap(v[parent], v[i]);
                i = parent;
            }
            else
            {
                return;
            }
        }
    }

    void pop()
    {
        if (v.size() <= 1)
        {
            return;
        }
        v[1] = v[v.size() - 1];
        v.pop_back();
        ptr -= 1;
        int i = 1;
        while (i < ptr)
        {
            int left = v[i * 2];
            int right = v[i * 2 + 1];
            int larger = left > right ? i * 2 : i * 2 + 1;
            if (v[i] < v[larger])
            {
                swap(v[i], v[larger]);
                i = larger;
            }
            else
            {
                return;
            }
        }
    }

    int top()
    {
        if (v.size() <= 1)
        {
            return -1;
        }
        return v[1];
    }

    int size()
    {
        return v.size() - 1;
    }
};

int main()
{

    // DynamicArray arr(11);
    // for (int i = 0; i < 220; i++)
    // {
    //     arr.push(i);
    // }
    // for (int i = 0; i < arr.size(); i++)
    // {
    //     cout << arr[i] << " ";
    // }
    // cout << "\n";

    // LinkedList<int> ll;
    // for (int i = 1; i <= 10; i++)
    // {
    //     ll.append(i);
    //     ll.prepend(10 * i);
    // }
    // while (ll.iterate() != NULL)
    // {
    //     cout << ll.ptr() << " ";
    // }
    // cout << "\n";

    Hashset set = Hashset(10);
    set.insert(1);
    set.insert(2);
    set.insert(10);
    set.insert(11);
    set.insert(55);
    set.insert(99);
    cout << set.size() << " ";
    set.remove(2);
    cout << set.find(2) << "\n";

    return 0;
}