#include "QuadTree.h"
#include <set>

namespace hw6 {

/*
 * QuadNode
 */
void QuadNode::split(size_t capacity) {
    for (int i = 0; i < 4; ++i) {
        delete children[i];
        children[i] = nullptr;
    }

    // Task construction
    // TODO

    Envelope e[4];
    double midX = (bbox.getMinX() + bbox.getMaxX()) / 2;
    double midY = (bbox.getMinY() + bbox.getMaxY()) / 2;
    e[0] = Envelope(bbox.getMinX(), midX, bbox.getMinY(), midY);
    e[1] = Envelope(midX, bbox.getMaxX(), bbox.getMinY(), midY);
    e[2] = Envelope(bbox.getMinX(), midX, midY, bbox.getMaxY());
    e[3] = Envelope(midX, bbox.getMaxX(), midY, bbox.getMaxY());

    children[0] = new QuadNode(e[0]);
    children[1] = new QuadNode(e[1]);
    children[2] = new QuadNode(e[2]);
    children[3] = new QuadNode(e[3]);
    for (int i = 0; i < features.size(); i++)
    {
        if (children[0]->getEnvelope().intersect(features[i].getEnvelope()))
        {
            children[0]->add(features[i]);
        }
        if (children[1]->getEnvelope().intersect(features[i].getEnvelope()))
        {
            children[1]->add(features[i]);
        }
        if (children[2]->getEnvelope().intersect(features[i].getEnvelope()))
        {
            children[2]->add(features[i]);
        }
        if (children[3]->getEnvelope().intersect(features[i].getEnvelope()))
        {
            children[3]->add(features[i]);
        }
    }
    features.clear();

    for (int i = 0; i < 4; i++)
    {
        if (children[i]->getFeatureNum() > capacity)
            children[i]->split(capacity);
    }
}

void QuadNode::countNode(int &interiorNum, int &leafNum) {
    if (isLeafNode()) {
        ++leafNum;
    } else {
        ++interiorNum;
        for (int i = 0; i < 4; ++i)
            children[i]->countNode(interiorNum, leafNum);
    }
}

int QuadNode::countHeight(int height) {
    ++height;
    if (!isLeafNode()) {
        int cur = height;
        for (int i = 0; i < 4; ++i) {
            height = std::max(height, children[i]->countHeight(cur));
        }
    }
    return height;
}

void QuadNode::rangeQuery(const Envelope &rect,
                          std::vector<Feature> &features) {
    if (!bbox.intersect(rect))
        return;

    // Task range query
    // TODO

    for (int i = 0; i < 4; i++) {
        Envelope e = children[i]->bbox;
        if (bbox.intersect(children[i]->bbox)) {
            if (children[i]->isLeafNode()) {
                //store feature								
                for (int j = 0; j < children[i]->getFeatureNum(); j++) {
                    Feature f = children[i]->getFeature(j);
                    features.push_back(f);
                }

            }
            else {
                children[i]->rangeQuery(rect, features);
            }

        }
    }
}

QuadNode *QuadNode::pointInLeafNode(double x, double y) {
    // Task NN query
    // TODO

    if (isLeafNode()) {
        return this;
    }
    else {
        for (int i = 0; i < 4; i++) {
            if (children[i]->bbox.contain(x, y)) {
                return children[i]->pointInLeafNode(x, y);
            }
        }
    }
    return nullptr;
}

void QuadNode::draw() {
    if (isLeafNode()) {
        bbox.draw();
    } else {
        for (int i = 0; i < 4; ++i)
            children[i]->draw();
    }
}

/*
 * QuadTree
 */
bool QuadTree::constructTree(const std::vector<Feature> &features) {
    if (features.empty())
        return false;

    // Task construction
    // TODO
    Envelope ori_env(features[0].getEnvelope());//相当于空包围盒
    for (int i = 0; i < features.size(); i++)
    {
        ori_env = ori_env.unionEnvelope(features[i].getEnvelope());//合并所有要素的包围盒
    }
    root = new QuadNode(ori_env);
    root->add(features);
    if (root->getFeatureNum() > capacity)
        root->split(capacity);

    bbox = ori_env;
    //bbox = Envelope(-74.1, -73.8, 40.6, 40.8); // 注意此行代码需要更新为features的包围盒，或根节点的包围盒
    return true;
}

void QuadTree::countNode(int &interiorNum, int &leafNum) {
    interiorNum = 0;
    leafNum = 0;
    if (root)
        root->countNode(interiorNum, leafNum);
}

void QuadTree::countHeight(int &height) {
    height = 0;
    if (root)
        height = root->countHeight(0);
}

void QuadTree::rangeQuery(const Envelope &rect,
                          std::vector<Feature> &features) {
    features.clear();

    // Task range query
    // TODO
    if (rect.intersect(bbox))root->rangeQuery(rect, features);

    // filter step (选择查询区域与几何对象包围盒相交的几何对象)

    // 注意四叉树区域查询仅返回候选集，精炼步在hw6的rangeQuery中完成
}

bool QuadTree::NNQuery(double x, double y, std::vector<Feature> &features) {

    if (!root || !(root->getEnvelope().contain(x, y)))
        return false;

    // Task NN query
    // TODO

    // filter step
    // (使用maxDistance2Envelope函数，获得查询点到几何对象包围盒的最短的最大距离，然后区域查询获得候选集)

    const Envelope &envelope = root->getEnvelope();
    double minDist = std::max(envelope.getWidth(), envelope.getHeight());


    QuadNode* n = root->pointInLeafNode(x, y);
    for (int i = 0; i < n->getFeatureNum(); i++) {
        minDist = std::min(minDist, n->getFeature(i).maxDistance2Envelope(x, y));
    }
    Envelope rect = Envelope(x - minDist, x + minDist, y - minDist, y + minDist);
    rangeQuery(rect, features);

    // 注意四叉树邻近查询仅返回候选集，精炼步在hw6的NNQuery中完成

    return true;
}

void QuadTree::draw() {
    if (root)
        root->draw();
}

} // namespace hw6
