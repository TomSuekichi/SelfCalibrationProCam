#pragma once

#include "ofMain.h"
#include "ofxCVfunc.h"
#include "StructuredLight.h"
#include "RadialFundMatrix.h"
#include <Eigen/Core>
using Eigen::MatrixXd;

class ofApp : public ofBaseApp{

	public:
		void setup();
		void update();
		void draw();
		void setupProjectorWindow();
		void drawProjectorWindow(ofEventArgs & args);

		void keyPressed(int key);
		void keyReleased(int key);
		void mouseMoved(int x, int y );
		void mouseDragged(int x, int y, int button);
		void mousePressed(int x, int y, int button);
		void mouseReleased(int x, int y, int button);
		void mouseEntered(int x, int y);
		void mouseExited(int x, int y);
		void windowResized(int w, int h);
		void dragEvent(ofDragInfo dragInfo);
		void gotMessage(ofMessage msg);

		const int camW = 2048;
		const int camH = 1088;
		const int proW = 1280;
		const int proH = 720;

		ofVideoGrabber vidGrabber;
		ofImage colorImg;
		ofImage grayImg;

		StructuredLight structLight;
		int count;

		RadialFundMatrix radialFund;
		int num;
};
