#include "RadialFundMatrix.h"

void RadialFundMatrix::setup(vector<ofVec2f>& _camPosSet, vector<ofVec2f>& _proPosSet) {
	camPosSet = _camPosSet;
	proPosSet = _proPosSet;

	ac = camW / 2;
	bc = camH / 2;
	ap = proW / 2;
	bp = proH / 2;
	fc = sqrt(pow(camW, 2) + pow(camH, 2));
	fp = sqrt(pow(proW, 2) + pow(proH, 2));
}

void RadialFundMatrix::drawTest(int _num) {
	ofSetColor(255);
	ofSetLineWidth(1);
	//ofDrawLine(inlierCamPos[_num] / 5, inlierProPos[_num] * ofVec2f(float(camW) / proW, float(camH) / proH) / 5 + ofVec2f(camW / 5, 0));

	for (int i = 0; i < inlierCamPos.size(); i++) {
		if (i % 15000 == 0) {
			//ofSetHSVColor(i * 360 / inlierCamPos.size(), 100, 100);
			ofDrawLine(inlierCamPos[i] / 5, inlierProPos[i] * ofVec2f(float(camW) / proW, float(camH) / proH) / 5 );
		}
	}
	ofSetColor(255);
}

void RadialFundMatrix::calcLoop() {
	int numCorrespond = camPosSet.size();
	float minDsum = 100000000;
	maxInlierNum = 0;
	preInlierCamPos.clear();
	preInlierProPos.clear();
	for (int k = 0; k < 100; k++) {
		randomNums.clear();
		bool isNum = false;
		inlierNum = 0;
		while (randomNums.size() < 15) {
			isNum = false;
			int num = ofRandom(0, camPosSet.size());
			for (int i = 0; i < randomNums.size(); i++) {
				if (num == randomNums[i]) {
					isNum = true;
					break;
				}
			}
			if (isNum == false) {
				randomNums.push_back(num);
			}
		}
		randomCamPosSet.clear();
		randomProPosSet.clear();
		for (int i = 0; i < randomNums.size(); i++) {
			float x = camPosSet[randomNums[i]].x;
			float y = camPosSet[randomNums[i]].y;
			float u = proPosSet[randomNums[i]].x;
			float v = proPosSet[randomNums[i]].y;
			ofVec4f campos = ofVec4f((x*x) + (y*y), x, y, 1);
			ofVec4f propos = ofVec4f((u*u) + (v*v), u, v, 1);
			//printf("%f, %f, %f, %f\n", x, y, u, v);
			randomCamPosSet.push_back(campos);
			randomProPosSet.push_back(propos);
		}
		RadialFundMat = calcfundMat(randomCamPosSet, randomProPosSet);
		//float dSum = 0;
		for (int i = 0; i < camPosSet.size(); i++) {
			float x = camPosSet[i].x;
			float y = camPosSet[i].y;
			float u = proPosSet[i].x;
			float v = proPosSet[i].y;
			Vector2d campos;
			campos << (x*x) + (y*y), x, y, 1;
			Vector2d propos;
			propos << (u*u) + (v*v), u, v, 1;
			float d = propos.transpose() * RadialFundMat*campos;
			if (pow(d, 2) < 100000) {
				inlierNum += 1;
				preInlierCamPos.push_back(camPosSet[i]);
				preInlierProPos.push_back(proPosSet[i]);
			}
			//dSum += pow(d, 2);
			//printf("%f, ", d);
		}
		//printf("%d", inlierNum);
		if (maxInlierNum < inlierNum) {
			maxInlierNum = inlierNum;
			optimizedRadialFundMat = RadialFundMat;
			inlierCamPos = preInlierCamPos;
			inlierProPos = preInlierProPos;
		}
		//dSum /= camPosSet.size();
		//if (dSum < minDsum) {
		//	minDsum = dSum;
		//	optimizedRadialFundMat = RadialFundMat;
		//}
		//printf("%f", dSum);
		
		
	}
	//printf("\n");
	std::cout << optimizedRadialFundMat << std::endl;
	ofVec2f dd = calcInitialDistortion(optimizedRadialFundMat, ac, bc, ap, bp);
	dc = dd.x;
	dp = dd.y;
	FundMat = calcInitialFundMat(optimizedRadialFundMat, ac, bc, ap, bp, dc, dp);
	EssentialMat = calcInitialEssentialMat(FundMat, fc, fp, ac, bc, ap, bp);
	JacobiSVD<MatrixXd> svdE(EssentialMat, ComputeThinU | ComputeThinV);
	printf("sing\n");
	std::cout << svdE.singularValues() << std::endl;
	MatrixXd matD;
	Vector4d diag;
	diag << 1, 1, 0;
	matD = diag.asDiagonal();
	EssentialMat = svdE.matrixU()*matD*(svdE.matrixV().transpose());
	std::cout << EssentialMat << std::endl;

	//5
	FundMat = MatKp.transpose()*EssentialMat*MatKc;
	std::cout << FundMat << std::endl;

	//MatrixXd MatZ = MatrixXd::Zero(3, 3);
	//MatrixXd MatW = MatrixXd::Zero(3, 3);
	//MatZ(0, 1) = 1; MatZ(1, 0) = -1;
	//MatW(0, 1) = -1; MatW(1, 0) = 1; MatW(2, 2) = 1;
	//VectorXd e3(3);
	//e3 << 0, 0, 1;
	//printf("\n");
	//std::cout << e3 << std::endl;

	//printf("\n");
	//std::cout << svdE.matrixU() << std::endl;
	//VectorXd t1(3);
	//t1 << svdE.matrixU()(0, 2), svdE.matrixU()(1, 2), svdE.matrixU()(2, 2);
	//MatrixXd R1 = svdE.matrixU()*MatW*svdE.matrixV().transpose();
	//MatrixXd R2 = svdE.matrixU()*MatW.transpose()*svdE.matrixV().transpose();
	//MatrixXd MatR = R1.transpose()*R2;
	////MatrixXd t2 = -svdE.matrixU()*MatZ*svdE.matrixU().transpose();
	//printf("\n");
	//std::cout << t1 << std::endl;
	////std::cout << t2 << std::endl;
	//printf("\n");
	//std::cout << R1 << std::endl;
	//std::cout << R2 << std::endl;
	//std::cout << MatR << std::endl;

	//ofVec3f Angle1 = computeAnglesFromMatrix(R1);
	//Angle1 = Angle1 * 360 / (2 * PI);
	//printf("%f, %f, %f\n", Angle1.x, Angle1.y, Angle1.z);
	//ofVec3f Angle2 = computeAnglesFromMatrix(R2);
	//Angle2 = Angle2 * 360 / (2 * PI);
	//printf("%f, %f, %f\n", Angle2.x, Angle2.y, Angle2.z);
	//ofVec3f Angle = computeAnglesFromMatrix(MatR);
	//Angle = Angle * 360 / (2 * PI);
	//printf("%f, %f, %f\n", Angle.x, Angle.y, Angle.z);
}

void RadialFundMatrix::optimizedLM(MatrixXd _MatF, ofVec2f _pp, ofVec2f _pc, double _dp, double _dc) {
	const int num = 15;
	int info;

	VectorXd p(num);
	p << _MatF(0, 0), _MatF(0, 1), _MatF(0, 2), _MatF(1, 0), _MatF(1, 1), _MatF(1, 2), _MatF(2, 0), _MatF(2, 1), _MatF(2, 2), _pp.x, _pp.y, _pc.x, _pc.y, _dp, _dc;

	vector<double> camX, camY, proX, proY;
	for (int i = 0; i < camPosSet.size(); i++) {
		camX.push_back(camPosSet[i].x);
		camY.push_back(camPosSet[i].y);
		proX.push_back(proPosSet[i].x);
		proY.push_back(proPosSet[i].y);
	}
}

ofVec2f RadialFundMatrix::calcInitialDistortion(MatrixXd _radialFundMat, float _ac, float _bc, float _ap, float _bp) {
	double ab_c = pow(_ac, 2) + pow(_bc, 2);
	double ab_p = pow(_ap, 2) + pow(_bp, 2);
	double _dc = _radialFundMat(0, 0) / ((ab_c*_radialFundMat(0, 0)) + (_ac*_radialFundMat(0, 1)) + (_bc*_radialFundMat(0, 2)) + _radialFundMat(0, 3));
	double _dp = _radialFundMat(0, 0) / ((ab_p*_radialFundMat(0, 0)) + (_ap*_radialFundMat(1, 0)) + (_bp*_radialFundMat(2, 0)) + _radialFundMat(3, 0));
	return ofVec2f(_dc, _dp);
}
MatrixXd RadialFundMatrix::calcInitialFundMat(MatrixXd _radialFundMat, float _ac, float _bc, float _ap, float _bp, float _dc, float _dp) {
	MatrixXd MatDc(3, 4);
	MatrixXd MatDp(3, 4);
	MatDc(0, 0) = _dc * _ac; MatDc(0, 1) = 1 - (2 * _dc * pow(_ac, 2)); MatDc(0, 2) = -2 * _dc * _ac*_bc;          MatDc(0, 3) = _dc * _ac*(pow(_ac, 2) + pow(_bc, 2));
	MatDc(1, 0) = _dc * _bc; MatDc(1, 1) = -2 * _dc * _ac;              MatDc(1, 2) = 1 - (2 * _dc * pow(_bc, 2)); MatDc(1, 3) = _dc * _bc*(pow(_ac, 2) + pow(_bc, 2));
	MatDc(2, 0) = _dc;       MatDc(2, 1) = -2 * _dc * _ac;              MatDc(2, 2) = -2 * _dc *_bc;               MatDc(2, 3) = 1 +  _dc *(pow(_ac, 2) + pow(_bc, 2));

	MatDp(0, 0) = _dp * _ap; MatDp(0, 1) = 1 - (2 * _dp * pow(_ap, 2)); MatDp(0, 2) = -2 * _dp * _ap*_bp;          MatDp(0, 3) = _dp * _ap*(pow(_ap, 2) + pow(_bp, 2));
	MatDp(1, 0) = _dp * _bp; MatDp(1, 1) = -2 * _dp * _ap;              MatDp(1, 2) = 1 - (2 * _dp * pow(_bp, 2)); MatDp(1, 3) = _dp * _bp*(pow(_ap, 2) + pow(_bp, 2));
	MatDp(2, 0) = _dp;       MatDp(2, 1) = -2 * _dp * _ap;              MatDp(2, 2) = -2 * _dp *_bp;               MatDp(2, 3) = 1 + _dp * (pow(_ap, 2) + pow(_bp, 2));
	std::cout << MatDc << std::endl;
	std::cout << MatDp << std::endl;
	MatrixXd MatF(3, 3);
	//MatrixXd mat = MatDp.transpose()*MatDp;
	MatrixXd MatDpT_inv = ((MatDp*MatDp.transpose()).inverse())*(MatDp);
	MatrixXd MatDcT_inv = ((MatDc*MatDc.transpose()).inverse())*(MatDc);
	MatrixXd mat1 = MatDc * MatDcT_inv.transpose();
	MatrixXd mat2 = MatDpT_inv * MatDp.transpose();
	printf("mat1, mat2\n");
	std::cout << mat1 << std::endl;
	std::cout << mat2 << std::endl;
	MatF = MatDpT_inv * _radialFundMat*MatDcT_inv.transpose();
	std::cout << MatF << std::endl;
	return MatF;
}

MatrixXd RadialFundMatrix::calcInitialEssentialMat(MatrixXd _fundMat, float _fc, float _fp, float _pc, float _qc, float _pp, float _qp) {
	MatKc = MatrixXd::Identity(3, 3);
	MatKp = MatrixXd::Identity(3, 3);
	MatKc(0, 0) = _fc; MatKc(1, 1) = _fc; MatKc(0, 2) = _pc; MatKc(1, 2) = _qc;
	MatKp(0, 0) = _fp; MatKp(1, 1) = _fp; MatKp(0, 2) = _pp; MatKp(1, 2) = _qp;
	
	MatrixXd MatE;
	MatE = MatKp.transpose().inverse()*_fundMat*MatKc.inverse();
	std::cout << MatE << std::endl;
	return MatE;
}


MatrixXd RadialFundMatrix::calcfundMat(vector<ofVec4f> _camPosSet, vector<ofVec4f> _proPosSet) {
	vector<ofVec4f> camPosSet = _camPosSet;
	vector<ofVec4f> proPosSet = _proPosSet;
	MatrixXd MatA(15, 16);
	for (int i = 0; i < camPosSet.size(); i++) {
		MatA(i, 0) = proPosSet[i].x*camPosSet[i].x;
		MatA(i, 1) = proPosSet[i].x*camPosSet[i].y;
		MatA(i, 2) = proPosSet[i].x*camPosSet[i].z;
		MatA(i, 3) = proPosSet[i].x*camPosSet[i].w;
		MatA(i, 4) = proPosSet[i].y*camPosSet[i].x;
		MatA(i, 5) = proPosSet[i].y*camPosSet[i].y;
		MatA(i, 6) = proPosSet[i].y*camPosSet[i].z;
		MatA(i, 7) = proPosSet[i].y*camPosSet[i].w;
		MatA(i, 8) = proPosSet[i].z*camPosSet[i].x;
		MatA(i, 9) = proPosSet[i].z*camPosSet[i].y;
		MatA(i, 10) = proPosSet[i].z*camPosSet[i].z;
		MatA(i, 11) = proPosSet[i].z*camPosSet[i].w;
		MatA(i, 12) = proPosSet[i].w*camPosSet[i].x;
		MatA(i, 13) = proPosSet[i].w*camPosSet[i].y;
		MatA(i, 14) = proPosSet[i].w*camPosSet[i].z;
		MatA(i, 15) = proPosSet[i].w*camPosSet[i].w;
	}
	//printf("matA");
	//std::cout << MatA << std::endl;
	JacobiSVD<MatrixXd> svd(MatA, ComputeThinU | ComputeThinV);
	MatrixXd matV = svd.matrixV();
	//printf("matV");
	//std::cout << matV << std::endl;
	VectorXd Vector_e(16);
	for (int i = 0; i < 16; i++) {
		Vector_e(i) = matV(i,14);
	}
	//printf("mate");
	//std::cout << Vector_e << std::endl;
	
	MatrixXd matE(4, 4);
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			matE(j, i) = matV(i * 4 + j, 14);
			//不確定事項：iとjは逆？
		}
	}

	JacobiSVD<MatrixXd> svdF(matE, ComputeThinU | ComputeThinV);
	MatrixXd matD = svdF.singularValues();
	//std::cout << matD << std::endl;
	//printf("\n");

	double d1, d2;
	d1 = matD(0, 0);
	d2 = matD(1, 0);
	Vector4d diag;
	diag << d1, d2, 0, 0;
	matD = diag.asDiagonal();
	//std::cout << matD << std::endl;


	MatrixXd matF = svdF.matrixU()*matD*(svdF.matrixV().transpose());
	//for (int i = 0; i < 4; i++) {
	//	for (int j = 0; j < 4; j++) {
	//		printf("%f,", matF(i, j));
	//	}
	//}
	//std::cout << matF(0) << std::endl;
	//std::cout << matF(1) << std::endl;
	//std::cout << matF(2) << std::endl;
	//std::cout << matF(3) << std::endl;
	printf("\n");
	return matE;
}

ofVec3f RadialFundMatrix::computeAnglesFromMatrix(MatrixXd _MatR) {
	double threshold = 0.001;
	double angle_x, angle_y, angle_z;

	if (abs(_MatR(2, 1) - 1.0) < threshold) { // R(2,1) = sin(x) = 1の時
		angle_x = PI / 2;
		angle_y = 0;
		angle_z = atan2(_MatR(1, 0), _MatR(0, 0));
	}
	else if (abs(_MatR(2, 1) + 1.0) < threshold) { // R(2,1) = sin(x) = -1の時
		angle_x = -PI / 2;
		angle_y = 0;
		angle_z = atan2(_MatR(1, 0), _MatR(0, 0));
	}
	else {
		angle_x = asin(_MatR(2, 1));
		angle_y = atan(-_MatR(2, 0)/ _MatR(2, 2) );
		angle_z = atan(-_MatR(0, 1)/ _MatR(1, 1) );
	}

	ofVec3f angle = ofVec3f(angle_x, angle_y, angle_z);
	return angle;
}