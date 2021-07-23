#pragma once
#include "ofMain.h"
#include "ofxBasicFunc.h"
#include "ofxCVfunc.h"
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/SVD>
using namespace Eigen;


class RadialFundMatrix {
public:
	void setup(vector<ofVec2f>& _camPosSet, vector<ofVec2f>& _proPosSet);
	void calcLoop();
	void drawTest(int _num);

private:
	vector<ofVec2f> camPosSet, proPosSet;
	vector<int> randomNums;
	vector<ofVec4f> randomCamPosSet;
	vector<ofVec4f> randomProPosSet;
	MatrixXd calcfundMat(vector<ofVec4f> _camPosSet, vector<ofVec4f> _proPosSet);
	ofVec2f calcInitialDistortion(MatrixXd _radialFundMat, float _ac, float _bc, float _ap, float _bp);
	MatrixXd calcInitialFundMat(MatrixXd _radialFundMat, float _ac, float _bc, float _ap, float _bp, float _dc, float _dp);
	MatrixXd calcInitialEssentialMat(MatrixXd _fundMat, float _fc, float _fp, float _pc, float qc, float _pp, float _qp);
	ofVec3f computeAnglesFromMatrix(MatrixXd _MatR);

	int camW = 2048;
	int camH = 1088;
	int proW = 1280;
	int proH = 720;
	float ap, bp, ac, bc;
	float dc, dp;
	float fc, fp;

	MatrixXd RadialFundMat;
	MatrixXd FundMat;
	MatrixXd EssentialMat;
	MatrixXd optimizedRadialFundMat;
	int inlierNum;
	int maxInlierNum;

	vector<ofVec2f> inlierCamPos;
	vector<ofVec2f> inlierProPos;
	vector<ofVec2f> preInlierCamPos;
	vector<ofVec2f> preInlierProPos;

};

