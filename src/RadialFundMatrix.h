#pragma once
#include "ofMain.h"
#include "ofxBasicFunc.h"
#include "ofxCVfunc.h"
#include <Eigen/Core>
#include <Eigen/LU>
#include <Eigen/SVD>
#include "Eigen/Dense"
#include "unsupported/Eigen/NonLinearOptimization"
#include "unsupported/Eigen/NumericalDiff"
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

	void optimizedLM(MatrixXd _MatF, ofVec2f _pp, ofVec2f _pc, double _dp, double _dc);

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
	MatrixXd MatKc;
	MatrixXd MatKp;

	int inlierNum;
	int maxInlierNum;

	vector<ofVec2f> inlierCamPos;
	vector<ofVec2f> inlierProPos;
	vector<ofVec2f> preInlierCamPos;
	vector<ofVec2f> preInlierProPos;


	template<typename _Scalar, int NX = Dynamic, int NY = Dynamic>
	struct Functor
	{
		typedef _Scalar Scalar;
		enum {
			InputsAtCompileTime = NX,
			ValuesAtCompileTime = NY
		};
		typedef Matrix<Scalar, InputsAtCompileTime, 1> InputType;
		typedef Matrix<Scalar, ValuesAtCompileTime, 1> ValueType;
		typedef Matrix<Scalar, ValuesAtCompileTime, InputsAtCompileTime> JacobianType;
	};

	struct misra1a_functor : Functor<double>
	{
		misra1a_functor(int inputs, int values, double *_camX, double *_camY, double *_proX, double *_proY)
			: inputs_(inputs), values_(values), _camX(_camX), _camY(_camY), _proX(_proX), _proY(_proY) {}

		double *_camX;
		double *_camY;
		double *_proX;
		double *_proY;
		int operator()(const VectorXd& b, VectorXd& fvec) const
		{
			for (int i = 0; i < values_; ++i) {
				//write objection function in this area
				//[example] fvec[i] = b[0] * (1.0 - exp(-b[1] * x[i])) - y[i];
				double Cr, Cf, Cp;

				//prepare
				MatrixXd _MatF(3, 3);
				_MatF(0, 0) = b[0]; _MatF(0, 1) = b[1]; _MatF(0, 2) = b[2];
				_MatF(1, 0) = b[3]; _MatF(1, 1) = b[4]; _MatF(1, 2) = b[5];
				_MatF(2, 0) = b[6]; _MatF(2, 1) = b[7]; _MatF(2, 2) = b[8];
				double _ap, _bp, _ac, _bc, _dp, _dc;
				_ap = b[9];
				_bp = b[10];
				_ac = b[11];
				_bc = b[12];
				_dp = b[13];
				_dc = b[14];

				MatrixXd _MatDp(3, 4); MatrixXd _MatDc(3, 4);
				_MatDc(0, 0) = _dc * _ac; _MatDc(0, 1) = 1 - (2 * _dc * pow(_ac, 2)); _MatDc(0, 2) = -2 * _dc * _ac*_bc;          _MatDc(0, 3) = _dc * _ac*(pow(_ac, 2) + pow(_bc, 2));
				_MatDc(1, 0) = _dc * _bc; _MatDc(1, 1) = -2 * _dc * _ac;              _MatDc(1, 2) = 1 - (2 * _dc * pow(_bc, 2)); _MatDc(1, 3) = _dc * _bc*(pow(_ac, 2) + pow(_bc, 2));
				_MatDc(2, 0) = _dc;       _MatDc(2, 1) = -2 * _dc * _ac;              _MatDc(2, 2) = -2 * _dc *_bc;               _MatDc(2, 3) = 1 + _dc * (pow(_ac, 2) + pow(_bc, 2));

				_MatDp(0, 0) = _dp * _ap; _MatDp(0, 1) = 1 - (2 * _dp * pow(_ap, 2)); _MatDp(0, 2) = -2 * _dp * _ap*_bp;          _MatDp(0, 3) = _dp * _ap*(pow(_ap, 2) + pow(_bp, 2));
				_MatDp(1, 0) = _dp * _bp; _MatDp(1, 1) = -2 * _dp * _ap;              _MatDp(1, 2) = 1 - (2 * _dp * pow(_bp, 2)); _MatDp(1, 3) = _dp * _bp*(pow(_ap, 2) + pow(_bp, 2));
				_MatDp(2, 0) = _dp;       _MatDp(2, 1) = -2 * _dp * _ap;              _MatDp(2, 2) = -2 * _dp *_bp;               _MatDp(2, 3) = 1 + _dp * (pow(_ap, 2) + pow(_bp, 2));

				MatrixXd _MatR(4, 4);
			}
			return 0;
		}
		/*
		  int df(const VectorXd& b, MatrixXd& fjac)
		  {
			  for (int i = 0; i < values_; ++i) {
			fjac(i, 0) = (1.0 - exp(-b[1]*x[i]));
			fjac(i, 1) = (b[0]*x[i] * exp(-b[1]*x[i]));
			  }
			  return 0;
		  }
		*/
		const int inputs_;
		const int values_;
		int inputs() const { return inputs_; }
		int values() const { return values_; }
	};
};

