#include "basicFunc.h"
float distance_pow(ofVec2f p1, ofVec2f p2) {
    float d;
    d = pow((p1.x - p2.x), 2) + pow((p1.y - p2.y), 2);
    return d;
}
float distance(ofVec2f p1, ofVec2f p2) {
    float d;
    d = sqrt(pow((p1.x - p2.x), 2) + pow((p1.y - p2.y), 2));
    return d;
}

float scalarOfVector(ofVec2f vec) {
    float scalar = pow(vec.x, 2) + pow(vec.y, 2);
    scalar=sqrt(scalar);
    return scalar;
}
float scalarOfVector_pow(ofVec2f vec) {
    float scalar = sqrt(pow(vec.x, 2) + pow(vec.y, 2));
    return scalar;
}
float CrossOfVector(ofVec2f vec1, ofVec2f vec2){
    float cross=vec1.x*vec2.y-vec1.y*vec2.x;
	return cross;
}
float calcAngle(ofVec2f middlePoint, ofVec2f pointA, ofVec2f pointB){
    ofVec2f vecA=pointA-middlePoint;
    ofVec2f vecB=pointB-middlePoint;
    float scalarA=scalarOfVector(vecA);
    float scalarB=scalarOfVector(vecB);
    float cosAngle=vecA.dot(vecB)/(scalarA*scalarB);
    float angle=acos(cosAngle);
    float cross=CrossOfVector(vecA, vecB);
    if (cosAngle>=1) { angle=0; }
    if (cross<0) { angle=(2*PI)-angle; }
    return angle;
}
float calcAngle(ofVec2f vecA, ofVec2f vecB, float scalarA, float scalarB){
    float cosAngle=vecA.dot(vecB)/(scalarA*scalarB);
    float angle=acos(cosAngle);
    float cross=CrossOfVector(vecA, vecB);
    if (cosAngle>=1) { angle=0; }
    if (cross<0) { angle=(2*PI)-angle; }
    return angle;

}

float mapValue(float value, float min, float max){
    float result;
    result=(value-min)/(max-min);
    return result;
}
float mapValue(float value, float min1, float max1, float min2, float max2){
    float result;
    float normalized=mapValue(value, min1, max1);
    result=normalized*(max2-min2)+min2;
    return result;
}

float calcInternalDivision(float value, float x1, float y1, float x2, float y2) {
	float a = (y2 - y1) / (x2 - x1);
	float b = y1 - (a*x1);
	float answer = a * value + b;
	return answer;
}

int GetDigit(int num) {
	int digit = 0;
	while (num != 0) {
		num /= 10;
		digit++;
	}
	return digit;
}


float calcDispersion(vector<ofVec2f> points){
    ofVec2f sumPos=ofVec2f(0,0);
    for (int i=0; i<points.size(); i++) {
        sumPos+=points[i];
    }
    ofVec2f avgPos=sumPos/points.size();
    
    float sumDiff=0;
    for (int i=0; i<points.size(); i++) {
        float diff=distance_pow(points[i], avgPos);
        sumDiff+=diff;
    }
    float avgDiff=sqrt(sumDiff/(points.size()+0.01f));
    
    return avgDiff;
}
ofVec2f calcAvgPos(vector<ofVec2f> points){
    ofVec2f sumPos=ofVec2f(0,0);
    for (int i=0; i<points.size(); i++) {
        sumPos+=points[i];
    }
    ofVec2f avgPos=sumPos/points.size();
    return avgPos;
}




ofVec3f multMatToVec(ofMatrix3x3 _mat, ofVec3f _vec) {
    ofVec3f v;
    float x, y, z;
    for (int i = 0; i < 3; i++) {
        x = _mat.a*_vec.x + _mat.b*_vec.y + _mat.c*_vec.z;
        y = _mat.d*_vec.x + _mat.e*_vec.y + _mat.f*_vec.z;
        z = _mat.g*_vec.x + _mat.h*_vec.y + _mat.i*_vec.z;
    }
    v = ofVec3f(x, y, z);
    return v;
}
ofMatrix3x3 translateMat(float _x, float _y) {
    ofMatrix3x3 mat;
    mat.a = 1; mat.b = 0; mat.c = _x;
    mat.d = 0; mat.e = 1; mat.f = _y;
    mat.g = 0; mat.h = 0; mat.i = 1;
    return mat;
}
ofMatrix3x3 rotateMat(float _angle) {
    ofMatrix3x3 mat;
    mat.a = cos(_angle); mat.b = -sin(_angle);  mat.c = 0;
    mat.d = sin(_angle); mat.e = cos(_angle);   mat.f = 0;
    mat.g = 0;           mat.h = 0;             mat.i = 1;
    return mat;
}
ofMatrix3x3 scaleMat(float _x, float _y) {
    ofMatrix3x3 mat;
    mat.a = _x; mat.b = 0;  mat.c = 0;
    mat.d = 0;  mat.e = _y; mat.f = 0;
    mat.g = 0;  mat.h = 0;  mat.i = 1;
    return mat;
}
ofMatrix3x3 affinebyTRS(ofMatrix3x3 _matT, ofMatrix3x3 _matR, ofMatrix3x3 _matS) {
    ofMatrix3x3 mat;
    //mat=_matR.operator*(_matS);
    mat = _matT.operator*(_matR.operator*(_matS));
    return mat;
}

void drawContour(vector<ofVec2f> _contour) {
	ofBeginShape();
	for (int i = 0; i < _contour.size(); i++) {
		ofVertex(_contour[i].x, _contour[i].y);
	}
	ofVertex(_contour[0].x, _contour[0].y);
	ofEndShape();
}


void drawAllContours(vector<vector<ofVec2f>> _contours){
    for(int i=0; i<_contours.size(); i++){
        ofBeginShape();
        for (int j=0; j<_contours[i].size(); j++) {
            ofVertex(_contours[i][j].x, _contours[i][j].y);
        }
        ofVertex(_contours[i][0].x, _contours[i][0].y);
        ofEndShape();
    }
}

void changeContoursScale(vector<vector<ofVec2f>> _contours, vector<vector<ofVec2f>>* _outcontours, int w, int h, int w_out, int h_out){
    vector<ofVec2f> outcontour;
    for (int i = 0; i < _contours.size(); i++) {
        outcontour.clear();
        for (int j = 0; j < _contours[i].size(); j++) {
            outcontour.push_back(ofVec2f(_contours[i][j].x*w_out / w, _contours[i][j].y*h_out / h));
        }
        _outcontours->push_back(outcontour);
    }
}
void changeContourScale(vector<ofVec2f> _contour, vector<ofVec2f>* _outcontour, int w, int h, int w_out, int h_out) {
	_outcontour->clear();
	for (int i = 0; i < _contour.size(); i++) {
		_outcontour->push_back(ofVec2f(_contour[i].x*w_out / w, _contour[i].y*h_out / h));
	}
}
void contourToDence(vector<ofVec2f> _contour, vector<ofVec2f>* denceContour, vector<int>* indexes, float d) {
	vector<ofVec2f> contour = _contour;
	ofVec2f pos;
	denceContour->clear();
	indexes->clear();
	for (int i = 0; i < contour.size(); i++) {
		int ia = i + 1;
		if (ia > contour.size()-1) { ia = 0; }
		float dist = distance(contour[ia], contour[i]);
		ofVec2f vec = contour[ia] - contour[i];
		float dist_i = 0;
		pos = ofVec2f(0, 0);
		while (dist_i < 1) {
			if (dist_i == 0) {
				indexes->push_back(denceContour->size());
			}
			pos = contour[i] + (dist_i*vec);
			denceContour->push_back(pos);
			dist_i += (d / dist);
		}
	}
}
void edgeToLaplacian(vector<ofVec2f> _edges, vector<ofVec2f>* Laplacian, int _lapNum) {
	int lapNum = _lapNum;
	vector<ofVec2f> edges = _edges;
	Laplacian->clear();
	for (int i = 0; i < edges.size(); i++) {
		ofVec2f sumVec = ofVec2f(0, 0);
		for (int j = -(lapNum / 2); j <= (lapNum / 2); j++) {
			if (j != 0) {
				int id = i - j;
				if (id < 0) { id += edges.size(); }
				if (id > edges.size() - 1) { id -= edges.size(); }
				sumVec += edges[id];
			}
		}
		ofVec2f avgVec = sumVec / lapNum;
		ofVec2f Lap = edges[i] - avgVec;
		Laplacian->push_back(Lap);
	}
}
void edgeToLaplacianWithIndex(vector<ofVec2f> _edges, vector<ofVec2f>* Laplacian, vector<int>& _indexes, int _lapNum) {
	int lapNum = _lapNum;
	vector<ofVec2f> edges = _edges;
	vector<int> indexes = _indexes;
	Laplacian->clear();
	for (int i = 0; i < indexes.size(); i++) {
		ofVec2f sumVec = ofVec2f(0, 0);
		for (int j = -(lapNum / 2); j <= (lapNum / 2); j++) {
			if (j != 0) {
				int id = indexes[i] - j;
				if (id < 0) { id += edges.size(); }
				if (id > edges.size() - 1) { id -= edges.size(); }
				sumVec += edges[id];
			}
		}
		ofVec2f avgVec = sumVec / lapNum;
		ofVec2f Lap = edges[indexes[i]] - avgVec;
		Laplacian->push_back(Lap);
	}
}



bool vectorFinder(std::vector<int> vec, int number) {
	bool Find=false;
    auto itr = std::find(vec.begin(), vec.end(), number);
    size_t index = std::distance( vec.begin(), itr );
    if (index != vec.size()) { // 発見できたとき
		Find = true;
    }
    else { // 発見できなかったとき
		Find = false;
    }
	return Find;
}


ofVec3f calcDistanceLineToPoint(ofVec2f _p0, ofVec2f _p1, ofVec2f _p2) {
	ofVec2f p0 = _p0;
	ofVec2f p1 = _p1;
	ofVec2f p2 = _p2;

	float a = p2.x - p1.x;
	float b = p2.y - p1.y;
	
	float a1 = p1.x - p0.x;
	float b1 = p1.y - p0.y;
	float AB = (a*a) + (b*b);

	float t;
	float d;
	ofVec2f crossPoint;
	t = -((a*a1) + (b*b1)) / AB;
	if (t < 0) {
		crossPoint = p1;
		d = (a1*a1) + (b1*b1);
	}else if(t > 1) {
		crossPoint = p2;
		float a2 = p2.x - p0.x;
		float b2 = p2.y - p0.y;
		d = (a2*a2) + (b2*b2);
	}
	else {
		crossPoint = ofVec2f(a*t + p1.x, b*t + p1.y);
		d = ((a*b1) - (b * a1))*((a*b1) - (b * a1)) / AB;
	}

	ofVec3f PandDist = ofVec3f(crossPoint.x, crossPoint.y, d);

	return PandDist;
}

ofVec3f calcLaplacianLineToPoint(ofVec2f _p0, ofVec2f _p1, ofVec2f _p2, ofVec2f _p3, ofVec2f _p4, ofVec2f _preLaplacian, ofVec4f elements) {
	ofVec2f p0 = _p0;
	ofVec2f p1 = _p1;
	ofVec2f p2 = _p2;
	ofVec2f p3 = _p3;
	ofVec2f p4 = _p4;

	ofVec2f preLaplacian = _preLaplacian;

	float a = p2.x - p1.x;
	float b = p2.y - p1.y;

	float a1 = p1.x - (((p3.x + p4.x) / 2) + (preLaplacian.x*elements.x) + (preLaplacian.y*elements.y));
	float b1 = p1.y - (((p3.y + p4.y) / 2) + (-preLaplacian.x*elements.y) + (preLaplacian.y*elements.x));
	float AB = (a*a) + (b*b);

	float t;
	float d;
	ofVec2f crossPoint;
	t = -((a*a1) + (b*b1)) / AB;
	if (t < 0) {
		crossPoint = p1;
		d = (a1*a1) + (b1*b1);
	}
	else if (t > 1) {
		crossPoint = p2;
		float a2 = p2.x - (((p3.x + p4.x) / 2) + (preLaplacian.x*elements.x) + (preLaplacian.y*elements.y));
		float b2 = p2.y - (((p3.y + p4.y) / 2) + (-preLaplacian.x*elements.y) + (preLaplacian.y*elements.x));
		d = (a2*a2) + (b2*b2);
	}
	else {
		crossPoint = ofVec2f(a*t + p1.x, b*t + p1.y);
		d = ((a*b1) - (b * a1))*((a*b1) - (b * a1)) / AB;
	}
	ofVec3f PandDist = ofVec3f(crossPoint.x, crossPoint.y, d);

	return PandDist;
}

ofVec3f calcLaplacianLineToPointInSmallArea(ofVec2f _p0, ofVec2f _p1, ofVec2f _p2, ofVec2f _p3, ofVec2f _p4, ofVec2f _preLaplacian, float _detectArea){
	ofVec2f p0 = _p0;
	ofVec2f p1 = _p1;
	ofVec2f p2 = _p2;
	ofVec2f p3 = _p3;
	ofVec2f p4 = _p4;
	float length = _detectArea;

	ofVec2f preLaplacian = _preLaplacian;

	float a = p2.x - p1.x;
	float b = p2.y - p1.y;

	float x1 = p1.x - p0.x;
	float y1 = p1.y - p0.y;
	float AB = (a*a) + (b*b);

	float D = 4 * ( (pow(length, 2)*AB)  -  pow((a*y1) - (b*x1), 2) );
	
	ofVec3f PandDist;
	if (D >= 0) {
		//printf("%f, ", D);
		float min_t = ( -((a*x1) + (b*y1))  -  sqrt(D) ) / AB;
		float max_t = (-((a*x1) + (b*y1))  +  sqrt(D)) / AB;

		float a1 = p1.x - (((p3.x + p4.x) / 2) + preLaplacian.x);
		float b1 = p1.y - (((p3.y + p4.y) / 2) + preLaplacian.y);

		float t;
		float d;
		ofVec2f crossPoint;
		t = -((a*a1) + (b*b1)) / AB;
		if (t < min_t) {
			//crossPoint = p1;
			crossPoint = ofVec2f(a*min_t + p1.x, b*min_t + p1.y);
			//d = (a1*a1) + (b1*b1);
			d = pow(a*min_t + a1, 2) + pow(b*min_t + b1, 2);
		}
		else if (t > max_t) {
			//crossPoint = p2;
			crossPoint = ofVec2f(a*max_t + p1.x, b*max_t + p1.y);
			float a2 = p2.x - (((p3.x + p4.x) / 2) + preLaplacian.x);
			float b2 = p2.y - (((p3.y + p4.y) / 2) + preLaplacian.y);
			//d = (a2*a2) + (b2*b2);
			d = pow(a*max_t + a1, 2) + pow(b*max_t + b1, 2);
		}
		else {
			crossPoint = ofVec2f(a*t + p1.x, b*t + p1.y);
			d = ((a*b1) - (b * a1))*((a*b1) - (b * a1)) / AB;
		}
		PandDist = ofVec3f(crossPoint.x, crossPoint.y, d);
	}
	else {
		PandDist = ofVec3f(p0.x, p0.y, 100000);
	}

	return PandDist;
}

ofVec2f detectMostNearpoint_edgevector(ofVec2f _point, vector<vector<ofVec2f>> _edgeVectors) {
	ofVec2f point = _point;
	vector<vector<ofVec2f>> edgeVectors = _edgeVectors;

	ofVec2f mostNearPoint;
	ofVec3f PandDist;
	ofVec2f nearPoint;
	float mostNearDist = 1000000;
	for (int i = 0; i < edgeVectors.size(); i++) {
		for (int j = 0; j < edgeVectors[i].size(); j++) {
			int ja = j + 1;
			if (ja > edgeVectors[i].size() - 1) { ja = 0; }
			PandDist = calcDistanceLineToPoint(point, edgeVectors[i][j], edgeVectors[i][ja]);
			nearPoint = ofVec2f(PandDist.x, PandDist.y);
			float dist = PandDist.z;
			if (dist < mostNearDist) {
				mostNearDist = dist;
				mostNearPoint = nearPoint;
			}
		}
	}

	return mostNearPoint;
}

ofVec3f detectMostNearLaplacianpoint_edgevector(ofVec2f _point, vector<vector<ofVec2f>> _edgeVectors, ofVec2f _point1, ofVec2f _point2, ofVec2f _preLaplacian, ofVec4f elements) {
	ofVec2f point = _point;
	ofVec2f point1 = _point1;
	ofVec2f point2 = _point2;

	vector<vector<ofVec2f>> edgeVectors = _edgeVectors;
	ofVec2f preLaplacian = _preLaplacian;

	ofVec2f mostNearPoint=point;
	ofVec3f mostNearPandD;
	ofVec3f PandDist;
	ofVec2f nearPoint;
	float mostNearDist = 1000000;
	for (int i = 0; i < edgeVectors.size(); i++) {
		for (int j = 0; j < edgeVectors[i].size(); j++) {
			int ja = j + 1;
			if (ja > edgeVectors[i].size() - 1) { ja = 0; }
			PandDist = calcLaplacianLineToPoint(point, edgeVectors[i][j], edgeVectors[i][ja], point1, point2, preLaplacian, elements);
			//PandDist = calcLaplacianLineToPointInSmallArea(point, edgeVectors[i][j], edgeVectors[i][ja], point1, point2, preLaplacian, 3);

			nearPoint = ofVec2f(PandDist.x, PandDist.y);
			float distPP = distance_pow(nearPoint, point);
			//if (distPP < 0.5) {
				float dist = PandDist.z;
				if (dist < mostNearDist) {
					mostNearDist = dist;
					mostNearPoint = nearPoint;
				}
			//}
		}
	}
	//printf("%f, ", mostNearDist);
	mostNearPandD = ofVec3f(mostNearPoint.x, mostNearPoint.y, mostNearDist);

	return mostNearPandD;
}


void calcLaplacian(vector<ofVec2f> _contour, vector<ofVec2f>* _laplacians) {
	vector<ofVec2f> contour = _contour;
	ofVec2f avgVec;
	ofVec2f laplacian;
	for (int i = 0; i < contour.size(); i++) {
		int ia = i + 1; int ib = i - 1;
		if (ia > contour.size() - 1) { ia = 0; }
		if (ib < 0) { ib = contour.size() - 1; }
		avgVec = (contour[ia] + contour[ib]) / 2;
		laplacian = contour[i] - avgVec;
		_laplacians->push_back(laplacian);
	}
}

ofVec4f calcSimilarityLaplacian(ofVec2f v1, ofVec2f v2, ofVec2f vv1, ofVec2f vv2) {
	float a, w, tx, ty;
	ofVec2f v12 = v1 - v2;
	ofVec2f vv12 = vv1 - vv2;

	float d = pow(v12.x, 2) + pow(v12.y, 2);
	a = (vv12.x*v12.x + vv12.y*v12.y) / d;
	w = (vv12.x*v12.y - vv12.y*v12.x) / d;
	tx = vv1.x - (a*v1.x) - (w*v1.y);
	ty = vv1.y - (a*v1.y) + (w*v1.x);

	ofVec4f elements = ofVec4f(a, w, tx, ty);
	return elements;
}


ofVec2f calcIntersectionPoint(ofVec2f _p1, ofVec2f _p2, ofVec2f _p3, ofVec2f _p4) {
	ofVec2f p1, p2, p3, p4;
	p1 = _p1; p2 = _p2; p3 = _p3; p4 = _p4;
	ofVec2f intersectPoint = ofVec2f(0, 0);

	float f3, f4;
	f3 = ((p2.x - p1.x)*(p3.y - p1.y)) - ((p2.y - p1.y)*(p3.x - p1.x));
	f4 = ((p2.x - p1.x)*(p4.y - p1.y)) - ((p2.y - p1.y)*(p4.x - p1.x));

	bool isIntersect = false;
	if (f3*f4 > 0) { isIntersect = true; }

	if (isIntersect) {
		float detA, t;
		detA = ((p1.x - p2.x)*(p4.x - p3.x)) - ((p4.x - p3.x)*(p1.y - p2.y));
		t = (((p4.y - p3.y)*(p4.x - p2.x)) + ((p3.x - p4.x)*(p4.y - p2.y))) / detA;

		if (0 < t&&t < 1) {
			intersectPoint.x = (t*p1.x) + ((1 - t)*p2.x);
			intersectPoint.y = (t*p1.y) + ((1 - t)*p2.y);
		}
		else {
			isIntersect = false;
		}
	}

	return intersectPoint;
}