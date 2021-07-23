#pragma once
#include "ofMain.h"

float distance_pow(ofVec2f p1, ofVec2f p2);
float distance(ofVec2f p1, ofVec2f p2);
float scalarOfVector(ofVec2f vec);
float scalarOfVector_pow(ofVec2f vec);
float CrossOfVector(ofVec2f vec1, ofVec2f vec2);
float calcAngle(ofVec2f middlePoint, ofVec2f pointA, ofVec2f pointB);
float calcAngle(ofVec2f vecA, ofVec2f vecB, float scalarA, float scalarB);
float mapValue(float value, float min, float max);
float mapValue(float value, float min1, float max1, float min2, float max2);
float calcInternalDivision(float value, float x1, float y1, float x2, float y2);
int GetDigit(int num);
float calcDispersion(vector<ofVec2f> points);
ofVec2f calcAvgPos(vector<ofVec2f> points);




ofVec3f multMatToVec(ofMatrix3x3 _mat, ofVec3f _vec);
ofMatrix3x3 translateMat(float _x, float _y);
ofMatrix3x3 rotateMat(float _angle);
ofMatrix3x3 scaleMat(float _x, float _y);
ofMatrix3x3 affinebyTRS(ofMatrix3x3 _matT, ofMatrix3x3 _matR, ofMatrix3x3 _matS);

/*contour*/
void drawContour(vector<ofVec2f> _contour);
void drawAllContours(vector<vector<ofVec2f>> _contours);
void changeContoursScale(vector<vector<ofVec2f>> _contours, vector<vector<ofVec2f>>* _outcontours, int w, int h, int w_out, int h_out);
void changeContourScale(vector<ofVec2f> _contour, vector<ofVec2f>* _outcontour, int w, int h, int w_out, int h_out);
void contourToDence(vector<ofVec2f> _contour, vector<ofVec2f>* denceContour, vector<int>* indexes, float d);
void edgeToLaplacian(vector<ofVec2f> _edges, vector<ofVec2f>* Laplacian, int _lapNum);
void edgeToLaplacianWithIndex(vector<ofVec2f> _edges, vector<ofVec2f>* Laplacian, vector<int>& _indexes, int _lapNum);
/*contour*/

bool vectorFinder(std::vector<int> vec, int number);


ofVec3f calcDistanceLineToPoint(ofVec2f _p0, ofVec2f _p1, ofVec2f _p2);
ofVec3f calcLaplacianLineToPoint(ofVec2f _p0, ofVec2f _p1, ofVec2f _p2, ofVec2f _p3, ofVec2f _p4, ofVec2f _preLaplacian, ofVec4f elements);
ofVec3f calcLaplacianLineToPointInSmallArea(ofVec2f _p0, ofVec2f _p1, ofVec2f _p2, ofVec2f _p3, ofVec2f _p4, ofVec2f _preLaplacian, float _detectArea);
ofVec2f detectMostNearpoint_edgevector(ofVec2f _point, vector<vector<ofVec2f>> _edgeVectors);
ofVec3f detectMostNearLaplacianpoint_edgevector(ofVec2f _point, vector<vector<ofVec2f>> _edgeVectors, ofVec2f _point1, ofVec2f _point2, ofVec2f _preLaplacian, ofVec4f elements);
void calcLaplacian(vector<ofVec2f> _contour, vector<ofVec2f>* _laplacians);
ofVec4f calcSimilarityLaplacian(ofVec2f v1, ofVec2f v2, ofVec2f vv1, ofVec2f vv2);

ofVec2f calcIntersectionPoint(ofVec2f _p1, ofVec2f _p2, ofVec2f _p3, ofVec2f _p4);