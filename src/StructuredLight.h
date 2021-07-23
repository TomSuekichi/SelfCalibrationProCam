#pragma once
#include "ofMain.h"
#include "ofxBasicFunc.h"
#include "ofxCVfunc.h"

class StructuredLight {
public:
	void setup(int _camW, int _camH);
	void update(int _count, ofImage& _baseImg);
	void draw();
	void drawGrayCode();
	vector<ofVec2f> camPosSet;
	vector<ofVec2f> proPosSet;
	int isCalculateRun;

private:
	void generateGrayCode(vector<int>* grayValue, int value, int num_bit);
	int camW, camH;
	const int proW = 1280;
	const int proH = 720;

	int num_bit_x, num_bit_y;
	vector<vector<int>> grayCodeSet_x, grayCodeSet_y;
	vector<ofFbo> grayCodeFboSet;
	vector<ofFbo> grayCodeNegaFboSet;
	vector<ofImage> savedImgSet;
	vector<ofImage> savedGrayImgSet;
	vector<ofImage> savedNegaImgSet;

	bool isSaved;
	int count;
	int preCount;
	ofFbo correspondFbo;
	bool isSetupFinished;
	ofImage baseImg;

	ofShader shader;

	ofFloatPixels correspondPixels;

};

