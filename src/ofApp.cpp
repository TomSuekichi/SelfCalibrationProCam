#include "ofApp.h"

//--------------------------------------------------------------
void ofApp::setup(){
	vidGrabber.setPixelFormat(OF_PIXELS_NATIVE);
	vidGrabber.setDeviceID(0);
	vidGrabber.setVerbose(true);
	vidGrabber.setDesiredFrameRate(500);
	vidGrabber.initGrabber(camW, camH);
	ofSetVerticalSync(false);
	colorImg.allocate(camW, camH, OF_IMAGE_COLOR);
	grayImg.allocate(camW, camH, OF_IMAGE_GRAYSCALE);
	structLight.setup(camW, camH);
	count = 0;
	num = 0;
}

//--------------------------------------------------------------
void ofApp::update(){
	bool bNewFrame = false;
	vidGrabber.update();
	bNewFrame = vidGrabber.isFrameNew();
	if (bNewFrame) {
		colorImg.setFromPixels(vidGrabber.getPixels());
		structLight.update(count, colorImg);
		if (structLight.isCalculateRun == 1) {
			radialFund.setup(structLight.camPosSet, structLight.proPosSet);
			radialFund.calcLoop();
			structLight.isCalculateRun += 1;
		}
	}
}

//--------------------------------------------------------------
void ofApp::draw() {
	structLight.draw();
	if (structLight.isCalculateRun > 1) {
		radialFund.drawTest(num);
	}
}

void ofApp::setupProjectorWindow() {
	ofSetVerticalSync(false);

}

void ofApp::drawProjectorWindow(ofEventArgs & args) {
	structLight.drawGrayCode();
}

//--------------------------------------------------------------
void ofApp::keyPressed(int key){
	switch (key)
	{
	case 'g':
		count++;
		break;
	case ' ':
		num++;
		break;
	}
}

//--------------------------------------------------------------
void ofApp::keyReleased(int key){

}

//--------------------------------------------------------------
void ofApp::mouseMoved(int x, int y ){

}

//--------------------------------------------------------------
void ofApp::mouseDragged(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mousePressed(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mouseReleased(int x, int y, int button){

}

//--------------------------------------------------------------
void ofApp::mouseEntered(int x, int y){

}

//--------------------------------------------------------------
void ofApp::mouseExited(int x, int y){

}

//--------------------------------------------------------------
void ofApp::windowResized(int w, int h){

}

//--------------------------------------------------------------
void ofApp::gotMessage(ofMessage msg){

}

//--------------------------------------------------------------
void ofApp::dragEvent(ofDragInfo dragInfo){ 

}
