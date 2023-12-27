#ifndef RTREE_SRC_H_INCLUDED
#define RTREE_SRC_H_INCLUDED

#include "Geometry.h"
#include "Tree.h"

#include <algorithm>
#include <array>
#include <queue>
#include <string>
#include <vector>

#include "CMakeIn.h"

namespace hw6 {

    /// <summary>
    /// </summary>
    /// <param name="M">Maximum number of children of each nodes.</param>
    template <uint8_t M> class RNode {
    private:
        RNode<M>* parent = nullptr;
        std::array<RNode<M>*, M> children = []() {
            decltype(children) ret;
            for (decltype(M) i = 0; i < M; ++i)
                ret[i] = nullptr;
            return ret;
        }();
        decltype(M) childrenNum = 0;
        Envelope bbox;
        std::vector<Feature> features;

    public:
        RNode() = delete;
        RNode(const Envelope& box) : bbox(box) {}

        inline bool isLeafNode() { return childrenNum == 0; }

        inline const Envelope& getEnvelope() { return bbox; }

        inline RNode<M>* getParent() { return parent; }

        inline void setEnvelope(const Envelope& box) { bbox = box; }

        inline RNode<M>* getChildNode(size_t i) {
            return i < childrenNum ? children[i] : nullptr;
        }

        inline const RNode<M>* getChildNode(size_t i) const {
            return i < childrenNum ? children[i] : nullptr;
        }

        inline decltype(M) getChildNum() const { return childrenNum; }

        inline size_t getFeatureNum() const { return features.size(); }

        inline const Feature& getFeature(size_t i) const { return features[i]; }

        inline const std::vector<Feature>& getFeatures() const { return features; }

        inline void add(const Feature& f) { features.push_back(f); }

        inline void add(RNode<M>* child) {
            children[childrenNum] = child;
            child->parent = this;
            ++childrenNum;
        }

        inline void remove(const Feature& f) {
            auto where = [&]() {
                for (auto itr = features.begin(); itr != features.end(); ++itr)
                    if (itr->getName() == f.getName())
                        return itr;
            }();
            features.erase(where);
            if (features.empty())
                features.shrink_to_fit(); // free memory unused but allocated
        }

        inline void remove(RNode<M>* child) {
            for (decltype(M) i = 0; i < childrenNum; ++i)
                if (children[i] == child) {
                    --childrenNum;
                    std::swap(children[i], children[childrenNum]);
                    children[childrenNum] = nullptr;
                    break;
                }
        }

        inline Feature popBackFeature() {
            auto ret = features.back();
            features.pop_back();
            return ret;
        }

        inline RNode<M>* popBackChildNode() {
            --childrenNum;
            auto ret = children[childrenNum];
            children[childrenNum] = nullptr;
            return ret;
        }

        void countNode(int& interiorNum, int& leafNum) {
            if (isLeafNode()) {
                ++leafNum;
            }
            else {
                ++interiorNum;
                for (decltype(M) i = 0; i < childrenNum; ++i)
                    children[i]->countNode(interiorNum, leafNum);
            }
        }

        int countHeight(int height) {
            ++height;
            if (!isLeafNode()) {
                int cur = height;
                for (decltype(M) i = 0; i < childrenNum; ++i)
                    height = max(height, children[i]->countHeight(cur));
            }
            return height;
        }

        inline void draw() {
            if (isLeafNode()) {
                bbox.draw();
            }
            else
                for (decltype(M) i = 0; i < childrenNum; ++i)
                    children[i]->draw();
        }

        void rangeQuery(const Envelope& rect, std::vector<Feature>& features) {
            // TODO
            //新添加的代码开始
            if (this->childrenNum == 0) {
                for (int j = 0; j < this->getFeatureNum(); j++) {
                    features.push_back(this->getFeature(j));
                }
            }
            for (int i = 0; i < this->childrenNum; i++) {//遍历R树
                if (this->getChildNode(i)->getEnvelope().intersect(rect)) {//MBR与Rect相交
                    if (this->getChildNode(i)->isLeafNode()) {//如果是叶节点将要素加入候选集
                        for (int j = 0; j < this->getChildNode(i)->getFeatureNum(); j++) {
                            features.push_back(this->getChildNode(i)->getFeature(j));
                        }
                    }
                    else {
                        this->children[i]->rangeQuery(rect, features);//不是则向下搜索
                    }
                }
            }
            //新添加的代码结束
        }
        //新添加的代码
        void UnionAllEnvelope() {//将结点的MBR更新
            Envelope Enve(DBL_MAX, -DBL_MAX, DBL_MAX, -DBL_MAX);
            if (!this->isLeafNode()) {
                for (int i = 0; i < this->getChildNum(); i++)
                    Enve = Enve.unionEnvelope(this->getChildNode(i)->getEnvelope());
            }
            else {
                for (int i = 0; i < this->getFeatureNum(); i++)
                    Enve = Enve.unionEnvelope(this->getFeature(i).getEnvelope());
            }
            this->setEnvelope(Enve);
        }

        RNode<M>* pointInLeafNode(double x, double y) {
            // TODO
            if (isLeafNode()) {//叶节点则结束
                return this;
            }
            else {
                for (int i = 0; i < this->childrenNum; i++) {
                    if (children[i]->bbox.contain(x, y)) {
                        return this->getChildNode(i)->pointInLeafNode(x, y);
                    }
                }
            }
            return nullptr;
        }
        //新添加的代码——调整R树的结构
        RNode<M>* AdjustTree(RNode<M>* New) {
            //如果不用调整，则递归至根节点返回
            if (New == nullptr) {
                if (parent == nullptr)
                    return this;
                return parent->AdjustTree(nullptr);
            }
            //如果是根节点
            if (this->getParent() == nullptr) {
                Envelope e(DBL_MAX, -DBL_MAX, DBL_MAX, -DBL_MAX);
                RNode<M>* NewRoot = new RNode<M>(e);//创建新的根节点
                NewRoot->add(New); NewRoot->add(this);//新的根节点将旧结点和增加的结点作为子节点
                this->parent = NewRoot;
                New->parent = NewRoot;
                New->UnionAllEnvelope(); this->UnionAllEnvelope();
                NewRoot->UnionAllEnvelope();
                return NewRoot;
            }
            if (parent->childrenNum < M) {//如果增加的结点加入后不需要继续分裂
                parent->add(New);
                parent->UnionAllEnvelope();
                return parent->AdjustTree(nullptr);
            }
            else {//需要继续分裂
                RNode<M>* OldP = parent;
                RNode<M>* NewP = parent->split(New);//分裂出新的子节点
                OldP->UnionAllEnvelope(); NewP->UnionAllEnvelope();
                return OldP->AdjustTree(NewP);//继续调整
            }
            return nullptr;
        }

        RNode<M>* split(RNode<M>* Added) {//将结点分裂
            vector<Envelope> Enve; Enve.clear();//判断结点分裂标准的包围盒集
            if (!this->isLeafNode()) {//如果不是叶节点，则用其子节点分裂
                for (int i = 0; i < this->getChildNum(); i++) {
                    Enve.push_back(this->getChildNode(i)->getEnvelope());
                }
                Enve.push_back(Added->getEnvelope());
            }
            else {//反之，用要素分裂
                for (int i = 0; i < this->getFeatureNum(); i++) {
                    Enve.push_back(this->getFeature(i).getEnvelope());
                }
            }
            vector<vector<int>> TwoSet(2);
            int Seed1, Seed2; double Max = -INFINITY;
            //计算D值
            for (int i = 0; i < Enve.size(); i++) {
                for (int j = i; j < Enve.size(); j++) {
                    double D = Enve[i].unionEnvelope(Enve[j]).getArea() - Enve[i].getArea() - Enve[j].getArea();
                    if (D > Max) {//最大D值所在对即种子值
                        Max = D; Seed1 = i; Seed2 = j;
                    }
                }
            }
            TwoSet[0].push_back(Seed1); TwoSet[1].push_back(Seed2);
            //分配至两组
            for (int i = 0; i < Enve.size(); i++) {
                if (i == Seed1 || i == Seed2)continue;
                double d1 = Enve[i].unionEnvelope(Enve[Seed1]).getArea() - Enve[i].getArea() - Enve[Seed1].getArea(),
                    d2 = Enve[i].unionEnvelope(Enve[Seed2]).getArea() - Enve[i].getArea() - Enve[Seed2].getArea();
                if (d1 < d2) {//d值越小越紧密
                    TwoSet[0].push_back(i);
                }
                else {
                    TwoSet[1].push_back(i);
                }
            }
            if (this->isLeafNode()) {//叶节点对要素集进行分裂
                vector<Feature> Temp; Temp.clear();
                for (int i = 0; i < this->getFeatureNum(); i++) {
                    Temp.push_back(this->getFeature(i));
                }
                while (this->getFeatureNum() >= 1) {
                    this->popBackFeature();
                }
                Envelope e(DBL_MAX, -DBL_MAX, DBL_MAX, -DBL_MAX);
                RNode<M>* NewNode = new RNode<M>(e);
                this->setEnvelope(e);
                for (int x : TwoSet[0]) {
                    NewNode->add(Temp[x]);
                    NewNode->setEnvelope(NewNode->getEnvelope().unionEnvelope(Temp[x].getEnvelope()));
                }
                for (int x : TwoSet[1]) {
                    this->add(Temp[x]);
                    this->setEnvelope(this->getEnvelope().unionEnvelope(Temp[x].getEnvelope()));
                }
                return NewNode;
            }
            else {//非叶子结点，对子节点进行分裂
                Envelope e(DBL_MAX, -DBL_MAX, DBL_MAX, -DBL_MAX);
                RNode<M>* NewNode = new RNode<M>(e);
                RNode<M>* NewNode_1 = new RNode<M>(e);
                for (int x : TwoSet[0]) {
                    if (x == M) {
                        NewNode->add(Added);
                        NewNode->setEnvelope(NewNode->getEnvelope().unionEnvelope(Added->getEnvelope()));
                    }
                    else {
                        NewNode->add(this->getChildNode(x));
                        NewNode->setEnvelope(NewNode->getEnvelope().unionEnvelope(this->getChildNode(x)->getEnvelope()));
                    }
                }
                for (int x : TwoSet[1]) {
                    if (x == M) {
                        NewNode_1->add(Added);
                        NewNode_1->setEnvelope(NewNode_1->getEnvelope().unionEnvelope(Added->getEnvelope()));
                    }
                    else {
                        NewNode_1->add(this->getChildNode(x));
                        NewNode_1->setEnvelope(NewNode_1->getEnvelope().unionEnvelope(this->getChildNode(x)->getEnvelope()));
                    }
                }
                while (this->getChildNum() >= 1) {
                    this->popBackChildNode();
                }
                this->setEnvelope(e);
                while (NewNode_1->childrenNum >= 1) {
                    RNode<M>* temp = NewNode_1->popBackChildNode();
                    this->add(temp);
                    this->setEnvelope(this->getEnvelope().unionEnvelope(temp->getEnvelope()));
                }
                free(NewNode_1);
                return NewNode;
            }
            return nullptr;
        }
    };

    template <uint8_t M> class RTree : public Tree {
    private:
        // Vars here MAY be useful, but it's ok to delete them {
        constexpr static auto m = M >> 1;
        constexpr static auto M_minus_m = M - m;
        constexpr static double EQ_ERR = 0.0000000005;
        // }

        RNode<M>* root = nullptr;

    public:
        RTree() : Tree(M) { static_assert(M >= 4); }
        ~RTree() {
            if (root != nullptr)
                delete root;
            root = nullptr;
        }

        virtual void setCapacity(int capacity) override {
            // DO NOTHING, since capacity is immutable in R tree
        }
        virtual bool constructTree(const std::vector<Feature>& features) override {
            // TODO
            if (features.empty())return false;//如果要素集为空，创建失败
            Envelope e(DBL_MAX, -DBL_MAX, DBL_MAX, -DBL_MAX);
            root = new RNode<M>(e);
            for (Feature f : features) {
                this->root = this->Insert(f);//几何要素逐个插入
            }
            bbox = this->root->getEnvelope();//更新包围盒
            //bbox = Envelope(-74.1, -73.8, 40.6, 40.8);
            return true;
        }
        virtual RNode<M>* Insert(const Feature feature) {//插入几何要素
            RNode<M>* InsertNode = this->FindNode(feature);//找到增加几何要素后，包围盒面积增加最小的叶节点
            if (InsertNode->getFeatureNum() > M) {//几何要素大于M时
                RNode<M>* NewNode = InsertNode->split(nullptr);//分裂
                return InsertNode->AdjustTree(NewNode);//调整树结构
            }
            else {
                return this->root;//返回根节点
            }
        }

        virtual RNode<M>* FindNode(const Feature feature) {
            RNode<M>* p = this->root;
            p->setEnvelope(p->getEnvelope().unionEnvelope(feature.getEnvelope()));//更新包围盒
            while (!p->isLeafNode()) {//知道找到叶节点为止
                double MinArea = INFINITY;
                int record = -1;
                for (int i = 0; i < p->getChildNum(); i++) {//找出增量最小的结点
                    double Diff = p->getChildNode(i)->getEnvelope().unionEnvelope(feature.getEnvelope()).getArea() - p->getChildNode(i)->getEnvelope().getArea();
                    if (MinArea > Diff) {
                        MinArea = Diff;
                        record = i;
                    }
                }
                p = p->getChildNode(record);//向下查找
                p->setEnvelope(p->getEnvelope().unionEnvelope(feature.getEnvelope()));//更新包围盒
            }
            p->add(feature);
            return p;
        }
        virtual void countNode(int& interiorNum, int& leafNum) override {
            interiorNum = leafNum = 0;
            if (root != nullptr)
                root->countNode(interiorNum, leafNum);
        }

        virtual void countHeight(int& height) override {
            height = 0;
            if (root != nullptr)
                height = root->countHeight(height);
        }

        virtual void rangeQuery(const Envelope& rect,
            std::vector<Feature>& features) override {
            features.clear();
            if (root != nullptr)
                root->rangeQuery(rect, features);
        }

        virtual bool NNQuery(double x, double y,
            std::vector<Feature>& features) override {
            features.clear();
            // TODO
            if (!this->root || !(this->root->getEnvelope().contain(x, y)))
                return false;
            double MinMaxDist = std::fmax(root->getEnvelope().getWidth(), root->getEnvelope().getHeight());//最大最小距离初始化
            RNode<M>* p = pointInLeafNode(x, y);
            if (p != nullptr) {
                for (int i = 0; i < p->getFeatures().size(); i++) {
                    MinMaxDist = std::fmin(MinMaxDist, p->getFeature(i).maxDistance2Envelope(x, y));
                }
            }
            Envelope Rect = Envelope(x - MinMaxDist, x + MinMaxDist, y - MinMaxDist, y + MinMaxDist);//构建范围查询所用的较小范围
            this->rangeQuery(Rect, features);
            return true;
        }

        RNode<M>* pointInLeafNode(double x, double y) {
            if (root != nullptr)
                return root->pointInLeafNode(x, y);
            else
                return nullptr;
        }

        virtual void draw() override {
            if (root != nullptr)
                root->draw();
        }

    public:
        static void test(int t);

        static void analyse();
    };

} // namespace hw6

#endif // !RTREE_SRC_H_INCLUDED
