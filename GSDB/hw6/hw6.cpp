// hw6.cpp : �������̨Ӧ�ó������ڵ㡣
//

#include "Common.h"
#include "Geometry.h"
#include "shapelib/shapefil.h"

#include "CMakeIn.h" // ����ע��

#include <gl/freeglut.h> // GLUT��ͷ�ļ�

#include <cstdio>
#include <ctime>
#include <iostream>
#include <list>
#include <map>
#include <memory>
#include <vector>
#include<set>


#ifdef USE_RTREE
#include "RTree.h"
using TreeTy = hw6::RTree<8>;
#else
#include "QuadTree.h"
using TreeTy = hw6::QuadTree;
#endif

using namespace std;

int screenWidth = 640;
int screenHeight = 480;

double pointSize = 2.0;

int mode;

vector<hw6::Feature> features;
vector<hw6::Feature> roads;
vector<hw6::Feature> polygons;//New

bool showRoad = true;
bool showPoly = true;//New

unique_ptr<hw6::Tree> pointTree;
unique_ptr<hw6::Tree> roadTree;
unique_ptr<hw6::Tree> PolyTree;//New

bool showTree = false;

hw6::Feature nearestFeature;

bool firstPoint = true;
hw6::Point corner[2];
hw6::Envelope selectedRect;
vector<hw6::Feature> selectedFeatures;

/*
 * shapefile�ļ���name��geometry���Զ�ȡ
 */
vector<string> readName(const char* filename) {
    DBFHandle file = DBFOpen(filename, "r");

    vector<string> res;
    int cct = DBFGetRecordCount(file);
    res.reserve(cct);
    for (int i = 0; i < cct; ++i) {
        string a = DBFReadStringAttribute(file, i, 0);
        res.push_back(a);
    }

    DBFClose(file);

    return res;
}

vector<hw6::Geometry*> readGeom(const char* filename) {
    SHPHandle file = SHPOpen(filename, "r");

    int pnEntities, pnShapeType;
    double padfMinBound[4], padfMaxBound[4];
    SHPGetInfo(file, &pnEntities, &pnShapeType, padfMinBound, padfMaxBound);

    vector<hw6::Point> points;
    vector<hw6::Geometry*> geoms;
    geoms.reserve(pnEntities);
    switch (pnShapeType) {
    case SHPT_POINT:
        for (int i = 0; i < pnEntities; ++i) {
            SHPObject* pt = SHPReadObject(file, i);
            geoms.push_back(new hw6::Point(pt->padfY[0], pt->padfX[0]));
            SHPDestroyObject(pt);
        }
        break;

    case SHPT_ARC:
        for (int i = 0; i < pnEntities; ++i) {
            points.clear();
            SHPObject* pt = SHPReadObject(file, i);
            for (int j = 0; j < pt->nVertices; ++j) {
                points.push_back(hw6::Point(pt->padfY[j], pt->padfX[j]));
            }
            SHPDestroyObject(pt);
            geoms.push_back(new hw6::LineString(points));
        }
        break;

    case SHPT_POLYGON:
        for (int i = 0; i < pnEntities; ++i) {
            points.clear();
            SHPObject* pt = SHPReadObject(file, i);
            for (int j = 0; j < pt->nVertices; ++j) {
                points.push_back(hw6::Point(pt->padfY[j], pt->padfX[j]));
            }
            SHPDestroyObject(pt);
            hw6::LineString line(points);
            hw6::Polygon* poly = new hw6::Polygon(line);
            geoms.push_back(new hw6::Polygon(line));
        }
        break;
    }

    SHPClose(file);
    return geoms;
}

void transformValue(double& res, const char* format = "%.2lf") {
    char buf[20];
    sprintf(buf, format, res);
    sscanf(buf, "%lf", &res);
}

void wrongMessage(hw6::Envelope e1, hw6::Envelope e2, bool cal) {
    cout << "Your answer is " << cal << " for test between ";
    e1.print();
    cout << " and ";
    e2.print();
    cout << ", but the answer is " << !cal << endl;
}

void wrongMessage(const hw6::Point& pt1, const hw6::Point& pt2, double dis,
    double res) {
    cout << "Your answer is " << dis << " for test ";
    pt1.print();
    cout << " and ";
    pt2.print();
    cout << ", but the answer is " << res << endl;
}

void wrongMessage(hw6::Envelope e1, hw6::Envelope e2, hw6::Envelope cal,
    hw6::Envelope res) {
    cout << "Your answer is ";
    cal.print();
    cout << " for test between ";
    e1.print();
    cout << " and ";
    e2.print();
    cout << ", but the answer is ";
    res.print();
    cout << endl;
}

/*
 * ���������Ϣ
 */
void printGeom(vector<hw6::Geometry*>& geom) {
    cout << "Geometry:" << endl;
    for (vector<hw6::Geometry*>::iterator it = geom.begin(); it != geom.end();
        ++it) {
        (*it)->print();
    }
}

/*
 * ɾ��������Ϣ
 */
void deleteGeom(vector<hw6::Geometry*>& geom) {
    for (vector<hw6::Geometry*>::iterator it = geom.begin(); it != geom.end();
        ++it) {
        delete* it;
        *it = NULL;
    }
    geom.clear();
}

/*
 * ��ȡŦԼ��·����
 */
void loadRoadData() {
    vector<hw6::Geometry*> geom = readGeom(PROJ_SRC_DIR "/data/highway");

    roads.clear();
    for (size_t i = 0; i < geom.size(); ++i)
        roads.push_back(hw6::Feature(to_string(i), geom[i]));

    cout << "road number: " << geom.size() << endl;
    roadTree->setCapacity(20);
    roadTree->constructTree(roads);
}

/*
 * ��ȡ��������� //New
 */
void loadPolygonData() {
    vector<hw6::Geometry*> geom = readGeom(PROJ_SRC_DIR "/data/polygon");

    polygons.clear();
    for (size_t i = 0; i < geom.size(); ++i)
        polygons.push_back(hw6::Feature(to_string(i), geom[i]));

    cout << "Polygon number: " << geom.size() << endl;
    PolyTree->setCapacity(1);
    PolyTree->constructTree(polygons);
}

/*
 * ��ȡŦԼ���г����޵�����
 */
void loadStationData() {
    vector<hw6::Geometry*> geom = readGeom(PROJ_SRC_DIR "/data/station");
    vector<string> name = readName(PROJ_SRC_DIR "/data/station");

    features.clear();
    for (size_t i = 0; i < geom.size(); ++i)
        features.push_back(hw6::Feature(name[i], geom[i]));

    cout << "station number: " << geom.size() << endl;
    pointTree->setCapacity(5);
    pointTree->constructTree(features);
}

/*
 * ��ȡŦԼ���⳵�򳵵�����
 */
void loadTaxiData() {
    vector<hw6::Geometry*> geom = readGeom(PROJ_SRC_DIR "/data/taxi");
    vector<string> name = readName(PROJ_SRC_DIR "/data/taxi");

    features.clear();
    for (size_t i = 0; i < geom.size(); ++i)
        features.push_back(hw6::Feature(name[i], geom[i]));

    cout << "taxi number: " << geom.size() << endl;
    pointTree->setCapacity(100);
    pointTree->constructTree(features);
}

/*
 * �����ѯ
 */

void rangeQuery() {
    vector<hw6::Feature> candidateFeatures;

    // filter step (ʹ���Ĳ�����ò�ѯ����ͼ���������Χ���ཻ�ĺ�ѡ����
    if (mode == RANGEPOINT)
        pointTree->rangeQuery(selectedRect, candidateFeatures);
    else if (mode == RANGELINE)
        roadTree->rangeQuery(selectedRect, candidateFeatures);
    else if (mode == RANGEPOLY)//New
        PolyTree->rangeQuery(selectedRect, candidateFeatures);//New
    // refine step (��ȷ�ж�ʱ����Ҫȥ�أ������ѯ����ͼ��ζ�����ظ�����)
    // TODO
    vector<hw6::Feature>feature1;
    set<hw6::Feature>featureMap;
    for (hw6::Feature f : candidateFeatures) {
        if (featureMap.find(f) == featureMap.end()) {
            if (f.getGeom()->intersects(selectedRect))
                feature1.push_back(f);
            featureMap.insert(f);
        }
    }
    selectedFeatures = feature1;

    //else selectedFeatures = candidateFeatures;
}

/*
 * �ڽ���ѯ
 */
void NNQuery(hw6::Point p) {
    vector<hw6::Feature> candidateFeatures;

    // filter step (ʹ���Ĳ�����þ���Ͻ��ļ���������ѡ��)
    if (mode == NNPOINT)
        pointTree->NNQuery(p.getX(), p.getY(), candidateFeatures);
    else if (mode == NNLINE)
        roadTree->NNQuery(p.getX(), p.getY(), candidateFeatures);
    else if (mode == NNPOLY)//New
        PolyTree->NNQuery(p.getX(), p.getY(), candidateFeatures);//New

    // refine step (��ȷ�����ѯ���뼸�ζ���ľ���)
    // TODO
    if (!candidateFeatures.empty()) {
        size_t nearestFeatureIndex = 0;
        double minDist = INFINITY;
        for (size_t j = 0; j < candidateFeatures.size(); ++j) {
            if (candidateFeatures[j].distance(p.getX(), p.getY()) < minDist) {
                minDist = candidateFeatures[j].distance(p.getX(), p.getY());
                nearestFeatureIndex = j;
            }
        }
        nearestFeature = candidateFeatures[nearestFeatureIndex];
    }

}

/*
 * ����Ļ����ת������������
 */
void transfromPt(hw6::Point& pt) {
    const hw6::Envelope bbox = pointTree->getEnvelope();//***PolyTree2pointTree
    double width = bbox.getMaxX() - bbox.getMinX() + 0.002;
    double height = bbox.getMaxY() - bbox.getMinY() + 0.002;

    double x = pt.getX() * width / screenWidth + bbox.getMinX() - 0.001;
    double y = pt.getY() * height / screenHeight + bbox.getMinY() - 0.001;

    x = max(bbox.getMinX(), x);
    x = min(bbox.getMaxX(), x);
    y = max(bbox.getMinY(), y);
    y = min(bbox.getMaxY(), y);
    pt = hw6::Point(x, y);
}

/*
 * ���ƴ���
 */
void display() {
    // glClearColor(241 / 255.0, 238 / 255.0, 232 / 255.0, 0.0);
    glClearColor(1.0, 1.0, 1.0, 0.0);
    glClear(GL_COLOR_BUFFER_BIT);

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    const hw6::Envelope bbox = pointTree->getEnvelope();//***PolyTree2pointTree
    gluOrtho2D(bbox.getMinX() - 0.001, bbox.getMaxX() + 0.001,
        bbox.getMinY() - 0.001, bbox.getMaxY() + 0.001);

    // ��·����
    if (showRoad) {
        glColor3d(252 / 255.0, 214 / 255.0, 164 / 255.0);
        for (size_t i = 0; i < roads.size(); ++i)
            roads[i].draw();
    }
    // ����λ���
    if (showPoly) {//New
        glColor3d(0.0, 146 / 255.0, 247 / 255.0);//New
        for (size_t i = 0; i < polygons.size(); ++i)//New
            polygons[i].draw();//New
    }//New
    // �����
    if (!(mode == RANGELINE || mode == NNLINE)) {
        glPointSize((float)pointSize);
        glColor3d(0.0, 146 / 255.0, 247 / 255.0);
        for (size_t i = 0; i < features.size(); ++i)
            features[i].draw();
    }

    // �Ĳ�������
    if (showTree) {
        glColor3d(0.0, 146 / 255.0, 247 / 255.0);
        if (mode == RANGELINE || mode == NNLINE)//New
            roadTree->draw();
        else  if (mode == RANGEPOINT || mode == NNPOINT)
            pointTree->draw();
        else
            PolyTree->draw();
    }

    // �������������
    if (mode == NNPOINT) {
        glPointSize(5.0);
        glColor3d(0.9, 0.0, 0.0);
        nearestFeature.draw();
    }

    // ����������·�Ͷ���λ���
    if (mode == NNLINE || mode == NNPOLY) {//New
        glLineWidth(3.0);
        glColor3d(0.9, 0.0, 0.0);
        nearestFeature.draw();
        glLineWidth(1.0);
    }

    // ����ѡ�����
    if (mode == RANGEPOINT || mode == RANGELINE || mode == RANGEPOLY) {//New
        glColor3d(0.0, 0.0, 0.0);
        selectedRect.draw();
        glColor3d(1.0, 0.0, 0.0);
        for (size_t i = 0; i < selectedFeatures.size(); ++i)
            selectedFeatures[i].draw();
    }

    glFlush();
    glutSwapBuffers();
}

/*
 * ���ͼ��̽���
 */
void mouse(int button, int state, int x, int y) {
    if (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN) {
        if (mode == RANGEPOINT || mode == RANGELINE || mode == RANGEPOLY) {//New
            if (firstPoint) {
                selectedFeatures.clear();
                corner[0] = hw6::Point(x, screenHeight - y);
                transfromPt(corner[0]);
            }
            else {
                corner[1] = hw6::Point(x, screenHeight - y);
                transfromPt(corner[1]);
                selectedRect =
                    hw6::Envelope(min(corner[0].getX(), corner[1].getX()),
                        max(corner[0].getX(), corner[1].getX()),
                        min(corner[0].getY(), corner[1].getY()),
                        max(corner[0].getY(), corner[1].getY()));
                rangeQuery();
            }
            firstPoint = !firstPoint;
            glutPostRedisplay();
        }
    }
}

void passiveMotion(int x, int y) {
    corner[1] = hw6::Point(x, screenHeight - y);

    if ((mode == RANGEPOINT || mode == RANGELINE || mode == RANGEPOLY) && !firstPoint) {//New
        corner[1] = hw6::Point(x, screenHeight - y);
        transfromPt(corner[1]);
        selectedRect = hw6::Envelope(min(corner[0].getX(), corner[1].getX()),
            max(corner[0].getX(), corner[1].getX()),
            min(corner[0].getY(), corner[1].getY()),
            max(corner[0].getY(), corner[1].getY()));
        rangeQuery();

        glutPostRedisplay();
    }
    else if (mode == NNPOINT || mode == NNLINE || mode == NNPOLY) {//New
        hw6::Point p(x, screenHeight - y);
        transfromPt(p);
        NNQuery(p);

        glutPostRedisplay();
    }
}

void changeSize(int w, int h) {
    screenWidth = w;
    screenHeight = h;
    glViewport(0, 0, w, h);
    glutPostRedisplay();
}

void processNormalKeys(unsigned char key, int x, int y) {
    switch (key) {
    case 27:
        exit(0);
        break;
    case 'N':
        mode = NNLINE;
        break;
    case 'n':
        mode = NNPOINT;
        break;
    case 'o':
        mode = NNPOLY;
        break;
    case 'S':
        mode = RANGELINE;
        firstPoint = true;
        break;
    case 's':
        mode = RANGEPOINT;
        firstPoint = true;
        break;
    case 'O'://New
        mode = RANGEPOLY;//New
        firstPoint = true;//New
        break;//New
    case 'B':
    case 'b':
        loadStationData();
        mode = Default;
        break;
    case 'T':
    case 't':
        loadTaxiData();
        mode = Default;
        break;
    case 'R':
    case 'r':
        showRoad = !showRoad;
        break;
    case 'Q':
    case 'q':
        showTree = !showTree;
        break;
    case '+':
        pointSize *= 1.1;
        break;
    case '-':
        pointSize /= 1.1;
        break;
    case '1':
    case '2':
    case '3':
    case '4':
    case '5':
    case '6':
    case '7':
    case '8':
        TreeTy::test(key - '0');
        break;
    default:
        mode = Default;
        break;
    }
    glutPostRedisplay();
}

int main(int argc, char* argv[]) {
    cout << "Key Usage:\n"
        << "  S  : range search for roads\n"
        << "  s  : range search for stations\n"
        << "  O  : range search for polygons\n"//New
        << "  N  : nearest road search\n"
        << "  n  : nearest station search\n"
        << "  o  : nearest polygon search\n"//New
        << "  B/b: Bicycle data\n"
        << "  T/t: Taxi data\n"
        << "  R/r: show Road\n"
        << "  Q/q: show Tree\n"
        << "  +  : increase point size\n"
        << "  -  : decrease point size\n"
        << "  1  : Test Envelope contain, interset and union\n"
        << "  2  : Test distance between Point and LineString\n"
        << "  3  : Test distance between Point and Polygon\n"
        << "  4  : Test tree construction\n"
        << "  5  : Test (your option here)\n"
        << "  8  : Tree performance analysis\n"
        << "  ESC: quit\n"
        << endl;

    pointTree = make_unique<TreeTy>();
    roadTree = make_unique<TreeTy>();
    PolyTree = make_unique<TreeTy>();//New

    loadRoadData();
    loadPolygonData();//New
    loadStationData();

    glutInit(&argc, argv);
    glutInitWindowSize(screenWidth, screenHeight);
    glutInitWindowPosition(0, 0);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE);
    glutCreateWindow("New York");

    glutMouseFunc(mouse);
    glutDisplayFunc(display);
    glutPassiveMotionFunc(passiveMotion);
    glutReshapeFunc(changeSize);
    glutKeyboardFunc(processNormalKeys);

    glutMainLoop();

    return 0;
}
