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
            //����ӵĴ��뿪ʼ
            if (this->childrenNum == 0) {
                for (int j = 0; j < this->getFeatureNum(); j++) {
                    features.push_back(this->getFeature(j));
                }
            }
            for (int i = 0; i < this->childrenNum; i++) {//����R��
                if (this->getChildNode(i)->getEnvelope().intersect(rect)) {//MBR��Rect�ཻ
                    if (this->getChildNode(i)->isLeafNode()) {//�����Ҷ�ڵ㽫Ҫ�ؼ����ѡ��
                        for (int j = 0; j < this->getChildNode(i)->getFeatureNum(); j++) {
                            features.push_back(this->getChildNode(i)->getFeature(j));
                        }
                    }
                    else {
                        this->children[i]->rangeQuery(rect, features);//��������������
                    }
                }
            }
            //����ӵĴ������
        }
        //����ӵĴ���
        void UnionAllEnvelope() {//������MBR����
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
            if (isLeafNode()) {//Ҷ�ڵ������
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
        //����ӵĴ��롪������R���Ľṹ
        RNode<M>* AdjustTree(RNode<M>* New) {
            //������õ�������ݹ������ڵ㷵��
            if (New == nullptr) {
                if (parent == nullptr)
                    return this;
                return parent->AdjustTree(nullptr);
            }
            //����Ǹ��ڵ�
            if (this->getParent() == nullptr) {
                Envelope e(DBL_MAX, -DBL_MAX, DBL_MAX, -DBL_MAX);
                RNode<M>* NewRoot = new RNode<M>(e);//�����µĸ��ڵ�
                NewRoot->add(New); NewRoot->add(this);//�µĸ��ڵ㽫�ɽ������ӵĽ����Ϊ�ӽڵ�
                this->parent = NewRoot;
                New->parent = NewRoot;
                New->UnionAllEnvelope(); this->UnionAllEnvelope();
                NewRoot->UnionAllEnvelope();
                return NewRoot;
            }
            if (parent->childrenNum < M) {//������ӵĽ��������Ҫ��������
                parent->add(New);
                parent->UnionAllEnvelope();
                return parent->AdjustTree(nullptr);
            }
            else {//��Ҫ��������
                RNode<M>* OldP = parent;
                RNode<M>* NewP = parent->split(New);//���ѳ��µ��ӽڵ�
                OldP->UnionAllEnvelope(); NewP->UnionAllEnvelope();
                return OldP->AdjustTree(NewP);//��������
            }
            return nullptr;
        }

        RNode<M>* split(RNode<M>* Added) {//��������
            vector<Envelope> Enve; Enve.clear();//�жϽ����ѱ�׼�İ�Χ�м�
            if (!this->isLeafNode()) {//�������Ҷ�ڵ㣬�������ӽڵ����
                for (int i = 0; i < this->getChildNum(); i++) {
                    Enve.push_back(this->getChildNode(i)->getEnvelope());
                }
                Enve.push_back(Added->getEnvelope());
            }
            else {//��֮����Ҫ�ط���
                for (int i = 0; i < this->getFeatureNum(); i++) {
                    Enve.push_back(this->getFeature(i).getEnvelope());
                }
            }
            vector<vector<int>> TwoSet(2);
            int Seed1, Seed2; double Max = -INFINITY;
            //����Dֵ
            for (int i = 0; i < Enve.size(); i++) {
                for (int j = i; j < Enve.size(); j++) {
                    double D = Enve[i].unionEnvelope(Enve[j]).getArea() - Enve[i].getArea() - Enve[j].getArea();
                    if (D > Max) {//���Dֵ���ڶԼ�����ֵ
                        Max = D; Seed1 = i; Seed2 = j;
                    }
                }
            }
            TwoSet[0].push_back(Seed1); TwoSet[1].push_back(Seed2);
            //����������
            for (int i = 0; i < Enve.size(); i++) {
                if (i == Seed1 || i == Seed2)continue;
                double d1 = Enve[i].unionEnvelope(Enve[Seed1]).getArea() - Enve[i].getArea() - Enve[Seed1].getArea(),
                    d2 = Enve[i].unionEnvelope(Enve[Seed2]).getArea() - Enve[i].getArea() - Enve[Seed2].getArea();
                if (d1 < d2) {//dֵԽСԽ����
                    TwoSet[0].push_back(i);
                }
                else {
                    TwoSet[1].push_back(i);
                }
            }
            if (this->isLeafNode()) {//Ҷ�ڵ��Ҫ�ؼ����з���
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
            else {//��Ҷ�ӽ�㣬���ӽڵ���з���
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
            if (features.empty())return false;//���Ҫ�ؼ�Ϊ�գ�����ʧ��
            Envelope e(DBL_MAX, -DBL_MAX, DBL_MAX, -DBL_MAX);
            root = new RNode<M>(e);
            for (Feature f : features) {
                this->root = this->Insert(f);//����Ҫ���������
            }
            bbox = this->root->getEnvelope();//���°�Χ��
            //bbox = Envelope(-74.1, -73.8, 40.6, 40.8);
            return true;
        }
        virtual RNode<M>* Insert(const Feature feature) {//���뼸��Ҫ��
            RNode<M>* InsertNode = this->FindNode(feature);//�ҵ����Ӽ���Ҫ�غ󣬰�Χ�����������С��Ҷ�ڵ�
            if (InsertNode->getFeatureNum() > M) {//����Ҫ�ش���Mʱ
                RNode<M>* NewNode = InsertNode->split(nullptr);//����
                return InsertNode->AdjustTree(NewNode);//�������ṹ
            }
            else {
                return this->root;//���ظ��ڵ�
            }
        }

        virtual RNode<M>* FindNode(const Feature feature) {
            RNode<M>* p = this->root;
            p->setEnvelope(p->getEnvelope().unionEnvelope(feature.getEnvelope()));//���°�Χ��
            while (!p->isLeafNode()) {//֪���ҵ�Ҷ�ڵ�Ϊֹ
                double MinArea = INFINITY;
                int record = -1;
                for (int i = 0; i < p->getChildNum(); i++) {//�ҳ�������С�Ľ��
                    double Diff = p->getChildNode(i)->getEnvelope().unionEnvelope(feature.getEnvelope()).getArea() - p->getChildNode(i)->getEnvelope().getArea();
                    if (MinArea > Diff) {
                        MinArea = Diff;
                        record = i;
                    }
                }
                p = p->getChildNode(record);//���²���
                p->setEnvelope(p->getEnvelope().unionEnvelope(feature.getEnvelope()));//���°�Χ��
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
            double MinMaxDist = std::fmax(root->getEnvelope().getWidth(), root->getEnvelope().getHeight());//�����С�����ʼ��
            RNode<M>* p = pointInLeafNode(x, y);
            if (p != nullptr) {
                for (int i = 0; i < p->getFeatures().size(); i++) {
                    MinMaxDist = std::fmin(MinMaxDist, p->getFeature(i).maxDistance2Envelope(x, y));
                }
            }
            Envelope Rect = Envelope(x - MinMaxDist, x + MinMaxDist, y - MinMaxDist, y + MinMaxDist);//������Χ��ѯ���õĽ�С��Χ
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
