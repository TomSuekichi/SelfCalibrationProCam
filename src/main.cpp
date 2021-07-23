#include "ofMain.h"
#include "ofApp.h"
#include "ofAppGLFWWindow.h"
//========================================================================
int main() {
	ofGLFWWindowSettings settings;
	ofVec2f displaySize, projectorSize;
	ofVec2f displayPos, projectorPos;
	displaySize = ofVec2f(1280, 1440);
	projectorSize = ofVec2f(1280, 720);
	displayPos = ofVec2f(0, 0);
	//projectorPos = ofVec2f(1280, 0);
	projectorPos = ofVec2f(2560, 0);

	settings.setSize(displaySize.x, displaySize.y);
	settings.setPosition(displayPos);
	settings.resizable = true;
	settings.decorated = false;
	shared_ptr<ofAppBaseWindow> mainWindow = ofCreateWindow(settings);

	settings.setSize(projectorSize.x, projectorSize.y);
	settings.setPosition(projectorPos);
	settings.resizable = true;
	settings.shareContextWith = mainWindow;
	shared_ptr<ofAppBaseWindow> guiWindow = ofCreateWindow(settings);

	shared_ptr<ofApp> mainApp(new ofApp);
	mainApp->setupProjectorWindow();
	ofAddListener(guiWindow->events().draw, mainApp.get(), &ofApp::drawProjectorWindow);

	//ofRunApp(guiWindow, guiApp);
	ofRunApp(mainWindow, mainApp);
	ofRunMainLoop();
}
