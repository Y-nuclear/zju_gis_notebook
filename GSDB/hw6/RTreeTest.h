#ifndef RTREE_TEST_H_INCLUDED
#define RTREE_TEST_H_INCLUDED

#include "Common.h"
#include "RTreeSrc.h"

#include "CMakeIn.h"

using namespace hw6;

extern int mode;
extern std::vector<Geometry *> readGeom(const char *filename);
extern std::vector<std::string> readName(const char *filename);
extern void transformValue(double &res, const char *format);
extern void wrongMessage(Envelope e1, Envelope e2, bool cal);
extern void wrongMessage(const Point &pt1, const Point &pt2, double dis,
                         double res);
extern void wrongMessage(Envelope e1, Envelope e2, Envelope cal, Envelope res);

namespace hw6 {

template <uint8_t M> void RTree<M>::test(int t) {
    using namespace std;

    std::cout << "*********************Start*********************" << std::endl;
    if (t == TEST1) {
        std::cout << "²âÊÔ1: Envelope Contain, Intersect, and Union" << endl;

        int failedCase = 0;
        Envelope e1(-1, 1, -1, 1);
        vector<Envelope> tests;
        tests.push_back(Envelope(-0.5, 0.5, -0.5, 0.5));
        tests.push_back(Envelope(-0.5, 0.5, 0.5, 1.5));
        tests.push_back(Envelope(0.5, 1.5, -0.5, 0.5));
        tests.push_back(Envelope(-1.5, -0.5, -1.5, -0.5));
        tests.push_back(Envelope(-2, -1, -0.5, 0.5));
        tests.push_back(Envelope(1, 1.5, 1, 1.5));
        tests.push_back(Envelope(-2, -1.5, -0.5, 0.5));
        tests.push_back(Envelope(-0.5, 0.5, 1.5, 2));
        tests.push_back(Envelope(-2, -1.5, 0.5, 1.5));
        tests.push_back(Envelope(0.5, 1.5, 1.5, 2));

        for (size_t i = 0; i < tests.size(); ++i) {
            if (e1.contain(tests[i]) != (i == 0)) {
                failedCase += 1;
                wrongMessage(e1, tests[i], (i != 0));
            }
            if (tests[i].contain(e1) == true) {
                failedCase += 1;
                wrongMessage(tests[i], e1, true);
            }
        }
        cout << "Envelope Contain: " << tests.size() * 2 - failedCase << " / "
             << tests.size() * 2 << " tests are passed" << endl;

        failedCase = 0;
        for (size_t i = 0; i < tests.size(); ++i) {
            if (e1.intersect(tests[i]) != (i < 6)) {
                failedCase += 1;
                wrongMessage(e1, tests[i], (i < 6));
            }
            if (tests[i].intersect(e1) != (i < 6)) {
                failedCase += 1;
                wrongMessage(tests[i], e1, (i < 6));
            }
        }
        cout << "Envelope Intersect: " << tests.size() * 2 - failedCase << " / "
             << tests.size() * 2 << " tests are passed" << endl;

        failedCase = 0;
        vector<Envelope> results;
        results.push_back(Envelope(-1, 1, -1, 1));
        results.push_back(Envelope(-1, 1, -1, 1.5));
        results.push_back(Envelope(-1, 1.5, -1, 1));
        results.push_back(Envelope(-1.5, 1, -1.5, 1));
        results.push_back(Envelope(-2, 1, -1, 1));
        results.push_back(Envelope(-1, 1.5, -1, 1.5));
        results.push_back(Envelope(-2, 1, -1, 1));
        results.push_back(Envelope(-1, 1, -1, 2));
        results.push_back(Envelope(-2, 1, -1, 1.5));
        results.push_back(Envelope(-1, 1.5, -1, 2));
        for (size_t i = 0; i < tests.size(); ++i) {
            if (e1.unionEnvelope(tests[i]) != results[i]) {
                failedCase += 1;
                wrongMessage(e1, tests[i], e1.unionEnvelope(tests[i]),
                             results[i]);
            }
            if (tests[i].unionEnvelope(e1) != results[i]) {
                failedCase += 1;
                wrongMessage(tests[i], e1, e1.unionEnvelope(tests[i]),
                             results[i]);
            }
        }
        cout << "Envelope Union: " << tests.size() * 2 - failedCase << " / "
             << tests.size() * 2 << " tests are passed" << endl;
    } else if (t == TEST2) {
        cout << "²âÊÔ2: Distance between Point and LineString" << endl;

        vector<Point> points;
        points.push_back(Point(0, 0));
        points.push_back(Point(10, 10));
        LineString line(points);

        points.push_back(Point(-10, -10));
        points.push_back(Point(20, 20));
        points.push_back(Point(5, 5));
        points.push_back(Point(10, 0));
        points.push_back(Point(10, -10));
        points.push_back(Point(0, 10));
        points.push_back(Point(0, 20));
        points.push_back(Point(20, 0));

        double dists[] = {0,       0,       14.1421, 14.1421, 0,
                          7.07107, 14.1421, 7.07107, 14.1421, 14.1421};

        int failedCase = 0;
        for (size_t i = 0; i < points.size(); ++i) {
            double dist = points[i].distance(&line);
            if (fabs(dist - dists[i]) > 0.0001) {
                failedCase += 1;
                cout << "Your answer is " << dist << " for test between ";
                line.print();
                cout << " and ";
                points[i].print();
                cout << ", but the answer is " << dists[i] << endl;
            }
        }
        cout << "Distance between Point and LineString: "
             << points.size() - failedCase << " / " << points.size()
             << " tests are passed" << endl;
    } else if (t == TEST3) {
        cout << "²âÊÔ3: Distance between Point and Polygon" << endl;

        vector<Point> points;
        points.push_back(Point(5, 0));
        points.push_back(Point(3, 6));
        points.push_back(Point(2, 4));
        points.push_back(Point(-2, 4));
        points.push_back(Point(-3, 5));
        points.push_back(Point(-5, 0));
        points.push_back(Point(0, -3));
        points.push_back(Point(5, 0));
        LineString line(points);
        Polygon poly(line);

        points.clear();
        points.push_back(Point(5, 4));
        points.push_back(Point(3, 4));
        points.push_back(Point(0, 4));
        points.push_back(Point(-3, 4));
        points.push_back(Point(-5, 4));
        points.push_back(Point(5, 5));
        points.push_back(Point(3, 5));
        points.push_back(Point(0, 5));
        points.push_back(Point(-3, 5));
        points.push_back(Point(0, 0));

        double dists[] = {1.26491, 0, 0, 0, 1.48556, 1.58114, 0, 1, 0, 0};

        int failedCase = 0;
        for (size_t i = 0; i < points.size(); ++i) {
            double dist = points[i].distance(&poly);
            if (fabs(dist - dists[i]) > 0.00001) {
                failedCase += 1;
                cout << "Your answer is " << dist << " for test between ";
                poly.print();
                cout << " and ";
                points[i].print();
                cout << ", but the answer is " << dists[i] << endl;
            }
        }
        cout << "Distance between Point and Polygon: "
             << points.size() - failedCase << " / " << points.size()
             << " tests are passed" << endl;
    } else if (t == TEST4) {
        cout << "²âÊÔ4: RTree Construction" << endl;
        int ncase, cct;
        ncase = cct = 2;

        RTree<8> rtree;
        vector<Geometry *> geom = readGeom(PROJ_SRC_DIR "/data/station");
        vector<Feature> features;

        for (size_t i = 0; i < geom.size(); ++i)
            features.push_back(Feature("", geom[i]));

        rtree.constructTree(features);

        int height, interiorNum, leafNum;
        rtree.countHeight(height);
        rtree.countNode(interiorNum, leafNum);

        if (!(height == 4 && interiorNum == 18 && leafNum == 86)) {
            cout << "Case 1: "
                 << "Your answer is height: " << height
                 << ", interiorNum: " << interiorNum << ", leafNum: " << leafNum
                 << ". One possible answer is height: 4, interiorNum: "
                    "18, "
                    "leafNum: 86\n";
            --cct;
        }

        features.clear();
        for (size_t i = 0; i < geom.size(); ++i)
            delete geom[i];
        geom.clear();

        vector<Geometry *> geom2 = readGeom(PROJ_SRC_DIR "/data/highway");
        vector<Feature> features2;
        RTree<8> rtree2;

        for (size_t i = 0; i < geom2.size(); ++i)
            features2.push_back(Feature("", geom2[i]));

        rtree2.constructTree(features2);

        int height2, interiorNum2, leafNum2;
        rtree2.countHeight(height2);
        rtree2.countNode(interiorNum2, leafNum2);

        if (!(height2 == 6 && interiorNum2 == 484 && leafNum2 == 2305)) {
            cout << "Case 2: "
                 << "Your answer is height: " << height2
                 << ", interiorNum: " << interiorNum2
                 << ", leafNum: " << leafNum2
                 << ". One possible answer is height: 6, interiorNum: "
                    "484, leafNum: 2305\n";
            --cct;
        }

        features2.clear();
        for (size_t i = 0; i < geom2.size(); ++i)
            delete geom2[i];
        geom2.clear();

        //cout << "RTree Construction: " << cct << " / " << ncase
        //     << " tests are passed" << endl;
    }else if (t == TEST5) {
        cout << "²âÊÔ5: Distance between LineString and LineString && LineString and Polygon" << endl;

        vector<Point> points;
        points.push_back(Point(0, 0));
        points.push_back(Point(5, 0));
        LineString line(points);
        vector<LineString> lines;
        points.clear();
        points.push_back(Point(0, 1));
        points.push_back(Point(1, 1));
        lines.push_back(LineString(points));
        points.push_back(Point(1, 0));
        lines.push_back(LineString(points));
        points.clear();
        points.push_back(Point(0, 1));
        points.push_back(Point(2, -1));
        lines.push_back(LineString(points));
        points.clear();
        points.push_back(Point(0, 1));
        points.push_back(Point(2, 0.5));
        points.push_back(Point(4, 2));
        lines.push_back(LineString(points));
        points.clear();
        points.push_back(Point(6, 1));
        points.push_back(Point(8, 8));
        lines.push_back(LineString(points));
        points.clear();
        points.push_back(Point(4, 1));
        points.push_back(Point(6, 0));
        lines.push_back(LineString(points));
        double dists[] = { 1, 0, 0, 0.5,1.4142135,0.447213};

        int failedCase = 0;
        for (size_t i = 0; i < lines.size(); ++i) {
            double dist = lines[i].distance(&line);
            if (fabs(dist - dists[i]) > 0.00001) {
                failedCase += 1;
                cout << "Your answer is " << dist << " for test between ";
                line.print();
                cout << " and ";
                lines[i].print();
                cout << ", but the answer is " << dists[i] << endl;
            }
        }
        cout << "Distance between LineString and LineString: "
            << lines.size() - failedCase << " / " << lines.size()
            << " tests are passed" << endl;
        vector<Polygon> poly;
        points.clear();
        lines.clear();
        points.push_back(Point(0, 1)); points.push_back(Point(1, 1));
        points.push_back(Point(1, 2)); points.push_back(Point(0, 1));
        lines.push_back(LineString(points));
        poly.push_back(Polygon(lines[0]));
        points.clear();
        lines.clear();
        points.push_back(Point(0, 1)); points.push_back(Point(1, 1));
        points.push_back(Point(1, 0)); points.push_back(Point(0, 1));
        lines.push_back(LineString(points));
        poly.push_back(Polygon(lines[0]));
        points.clear();
        lines.clear();
        points.push_back(Point(0, 1)); points.push_back(Point(1, 1));
        points.push_back(Point(2, -1)); points.push_back(Point(0, 1));
        lines.push_back(LineString(points));
        poly.push_back(Polygon(lines[0]));
        points.clear();
        lines.clear();
        points.push_back(Point(0, 1)); points.push_back(Point(1, 1));
        points.push_back(Point(2, 0.5)); points.push_back(Point(0, 1));
        lines.push_back(LineString(points));
        poly.push_back(Polygon(lines[0]));
        points.clear();
        lines.clear();
        points.push_back(Point(6, 1)); points.push_back(Point(8, 8));
        points.push_back(Point(8, 4)); points.push_back(Point(6, 1));
        lines.push_back(LineString(points));
        poly.push_back(Polygon(lines[0]));
        points.clear();
        lines.clear();
        points.push_back(Point(7, 0)); points.push_back(Point(5, 2));
        points.push_back(Point(8, 4)); points.push_back(Point(7, 0));
        lines.push_back(LineString(points));
        poly.push_back(Polygon(lines[0]));
        double dists1[] = { 1, 0, 0, 0.5,1.4142135,1.4142135 };

        failedCase = 0;
        for (size_t i = 0; i < poly.size(); ++i) {
            double dist = line.distance(&poly[i]);
            if (fabs(dist - dists1[i]) > 0.00001) {
                failedCase += 1;
                cout << "Your answer is " << dist << " for test between ";
                line.print();
                cout << " and ";
                poly[i].print();
                cout << ", but the answer is " << dists1[i] << endl;
            }
        }
        cout << "Distance between LineString and Polygon: "
            << poly.size() - failedCase << " / " << poly.size()
            << " tests are passed" << endl;
        Envelope E(2, 8, -2, 6);
        double dists2[] = { 0, 0, 1, 1,1,1 };

        failedCase = 0;
        for (size_t i = 0; i < poly.size(); ++i) {
            bool dist = E.intersectPoly(&poly[i]);
            if (dist != dists2[i]) {
                failedCase += 1;
                cout << "Your answer is " << dist << " for test between ";
                E.print();
                cout << " and ";
                poly[i].print();
                cout << ", but the answer is " << dists2[i] << endl;
            }
        }
        cout << "Distance between LineString and Polygon: "
            << poly.size() - failedCase << " / " << poly.size()
            << " tests are passed" << endl;
    }
    else if (t == TEST8) {
        cout << "²âÊÔ8: RTreeAnalysis" << endl;
        analyse();
    }

    cout << "**********************End**********************" << endl;
}

template <uint8_t I, uint8_t Last, uint8_t Step>
void forConstCapAnalyseRTree(const std::vector<Feature> &features) {
    if constexpr (I <= Last) {
        RTree<I> rtree;
        clock_t start_time, end_time;
        start_time = clock();
        rtree.constructTree(features);
        end_time = clock();
        int height, interiorNum, leafNum;
        rtree.countHeight(height);
        rtree.countNode(interiorNum, leafNum);
        cout << int(I) << "cons#Capacity " << int(I) << endl;
        printf("Height: %d \tInterior node number: %d \tLeaf node number: %d\n", height, interiorNum, leafNum);
        cout << "Construction Time:" << (double)(end_time - start_time) / 1000.0 << "s" << endl;
        // TODO
        double x, y;
        vector<Feature> candidateFeatures;
        start_time = clock();
        for (int i = 0; i < 100000; ++i) {
            x = -((rand() % 225) / 10000.0 + 73.9812);
            y = (rand() % 239) / 10000.0 + 40.7247;
            candidateFeatures.clear();

            rtree.NNQuery(x, y, candidateFeatures);
            // refine step
            // TODO

            Feature nearestFeature;
            if (!candidateFeatures.empty()) {
                size_t nearestFeatureIndex = 0;
                double minDist = INFINITY;
                for (size_t j = 0; j < candidateFeatures.size(); ++j) {
                    if (candidateFeatures[j].distance(x, y) < minDist) {
                        minDist = candidateFeatures[j].distance(x, y);
                        nearestFeatureIndex = j;
                    }
                }
                nearestFeature = candidateFeatures[nearestFeatureIndex];
            }
        }
        end_time = clock();


        cout << "NNQuery Time:" << (double)(end_time - start_time) / 1000.0 << "s" << endl;
        forConstCapAnalyseRTree<I + Step, Last, Step>(features);
    }
}

template <uint8_t M> void RTree<M>::analyse() {
    using namespace std;

    vector<Feature> features;
    vector<Geometry *> geom = readGeom(PROJ_SRC_DIR "/data/taxi");
    vector<string> name = readName(PROJ_SRC_DIR "/data/taxi");

    features.clear();
    features.reserve(geom.size());
    for (size_t i = 0; i < geom.size(); ++i)
        features.push_back(Feature(name[i], geom[i]));

    cout << "taxi number: " << geom.size() << endl;

    srand(time(nullptr));

    forConstCapAnalyseRTree<70, 200, 10>(features);
}


} // namespace hw6

#endif // !RTREE_TEST_H_INCLUDED