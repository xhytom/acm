class SegTree {
private:
    struct Node {
        Node () : left_(nullptr), right_(nullptr), val_(0), lazy_(0) {}
        int val_;
        int lazy_;
        Node* left_;
        Node* right_;
    };

public:
    Node* root_;
    SegTree() { root_ = new Node(); }
    ~SegTree() {}

    // 更新区间值
    void upDate(Node* curNode, int curLeft, int curRight, int upDateLeft, int upDateRight, int addVal) {
        if (upDateLeft <= curLeft && upDateRight >= curRight) {
            // 如果需要更新的区间[upDateLeft, upDateRight] 包含了 当前这个区间[curLeft, curRight] 
            // 那么暂存一下更新的值
            // 等到什么时候用到孩子结点了，再把更新的值发放给孩子
            curNode->val_ += addVal * (curRight - curLeft + 1);
            curNode->lazy_ += addVal;
            return;
        }

        // 到这里说明要用到左右孩子了
        // 因此，要用pushDown函数把懒标签的值传递下去
        int mid = (curLeft + curRight) / 2;
        pushDown(curNode, mid - curLeft + 1, curRight - mid);

        // 说明在[curLeft, curRight]中，
        if (upDateLeft <= mid) { 
            upDate(curNode->left_, curLeft, mid, upDateLeft, upDateRight, addVal);
        } 
        if (upDateRight > mid) {
            upDate(curNode->right_, mid + 1, curRight, upDateLeft, upDateRight, addVal);
        }

        // 更新了子节点还需要更新现在的结点
        pushUp(curNode);
    }
    

    // 把结点curNode的懒标记分发给左右孩子 然后自己的懒标记清零
    void pushDown(Node* curNode, int leftChildNum, int rightChildNum) {
        if (curNode->left_ == nullptr) curNode->left_ = new Node;
        if (curNode->right_ == nullptr) curNode->right_ = new Node;

        if (curNode->lazy_ == 0) return;

        curNode->left_->val_ += curNode->lazy_ * leftChildNum;
        curNode->left_->lazy_ += curNode->lazy_;

        curNode->right_->val_ += curNode->lazy_ * rightChildNum;
        curNode->right_->lazy_ += curNode->lazy_;

        curNode->lazy_ = 0;

        // 注意不需要递归再继续下推懒标签 
        // 每次只需要推一层即可
    }

    // 一般是子节点因为要被用到了，所以需要更新值 因此也要同时更新父节点的值
    void pushUp(Node* curNode) {
        curNode->val_ = curNode->left_->val_ + curNode->right_->val_;
    }

    // 查询
    int query(Node* curNode, int curLeft, int curRight, int queryLeft, int queryRight) {
        if (queryLeft <= curLeft && queryRight >= curRight) {
            return curNode->val_;
        }
        // 用到左右结点力 先下推！
        int mid = (curLeft + curRight) / 2;
        pushDown(curNode, mid - curLeft + 1, curRight - mid);

        int curSum = 0;
        if (queryLeft <= mid) curSum += query(curNode->left_, curLeft, mid, queryLeft, queryRight);
        if (queryRight > mid) curSum += query(curNode->right_, mid + 1, curRight, queryLeft, queryRight);

        return curSum;
    }
};