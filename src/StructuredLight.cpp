#include "StructuredLight.h"

void StructuredLight::setup(int _camW, int _camH) {
	num_bit_x = int(log2f(proW)) + 1;
	num_bit_y = int(log2f(proH)) + 1;
	camW = _camW;
	camH = _camH;

	/*generate grayCode*/
	for (int i = 0; i < proW; i++) {
		vector<int> grayCode;
		calcGrayValue(&grayCode, i, num_bit_x);
		grayCodeSet_x.push_back(grayCode);
	}
	for (int i = 0; i < proH; i++) {
		vector<int> grayCode;
		calcGrayValue(&grayCode, i, num_bit_y);
		grayCodeSet_y.push_back(grayCode);
	}
	/*generate grayCode*/

	/*create grayCodeImage*/
	ofFbo fbo;
	for (int i = 0; i < num_bit_x; i++) {
		fbo.clear();
		fbo.allocate(proW, proH);
		fbo.begin();
		ofClear(0);
		ofFill();
		int patternColor;
		for (int j = 0; j < proW; j++) {
			patternColor = grayCodeSet_x[j][i];
			ofSetColor(patternColor * 255);
			ofDrawRectangle(j, 0, 1, proH);
		}
		fbo.end();
		grayCodeFboSet.push_back(fbo);
	}
	for (int i = 0; i < num_bit_y; i++) {
		fbo.clear();
		fbo.allocate(proW, proH);
		fbo.begin();
		ofClear(0);
		ofFill();
		int patternColor;
		for (int j = 0; j < proH; j++) {
			patternColor = grayCodeSet_y[j][i];
			ofSetColor(patternColor * 255);
			ofDrawRectangle(0, j, proW, 1);
		}
		fbo.end();
		grayCodeFboSet.push_back(fbo);
	}
	fbo.clear();
	fbo.allocate(proW, proH);
	fbo.begin();
	ofClear(0);
	ofFill();
	ofSetColor(255);
	ofDrawRectangle(0, 0, proW, proH);
	fbo.end();
	grayCodeFboSet.push_back(fbo);
	/*create grayCodeImage*/
	/*create grayCodeNegaImage*/
	for (int i = 0; i < num_bit_x; i++) {
		fbo.clear();
		fbo.allocate(proW, proH);
		fbo.begin();
		ofClear(0);
		ofFill();
		int patternColor;
		for (int j = 0; j < proW; j++) {
			patternColor = grayCodeSet_x[j][i];
			ofSetColor((1 - patternColor) * 255);
			ofDrawRectangle(j, 0, 1, proH);
		}
		fbo.end();
		grayCodeNegaFboSet.push_back(fbo);
	}
	for (int i = 0; i < num_bit_y; i++) {
		fbo.clear();
		fbo.allocate(proW, proH);
		fbo.begin();
		ofClear(0);
		ofFill();
		int patternColor;
		for (int j = 0; j < proH; j++) {
			patternColor = grayCodeSet_y[j][i];
			ofSetColor((1 - patternColor) * 255);
			ofDrawRectangle(0, j, proW, 1);
		}
		fbo.end();
		grayCodeNegaFboSet.push_back(fbo);
	}
	fbo.clear();
	fbo.allocate(proW, proH);
	fbo.begin();
	ofClear(0);
	ofFill();
	ofSetColor(0);
	ofDrawRectangle(0, 0, proW, proH);
	fbo.end();
	grayCodeNegaFboSet.push_back(fbo);
	/*create grayCodeNegaImage*/

	isSaved = false;
	count = 0;
	preCount = 0;
	isCalculateRun = 0;
	correspondFbo.allocate(camW, camH, GL_RGBA32F);
	isSetupFinished = true;
	baseImg.allocate(camW, camH, OF_IMAGE_COLOR);
}

void StructuredLight::update(int _count, ofImage& _baseImg) {
	count = _count;
	baseImg = _baseImg;
	if (count != preCount && count % 2 == 1) {
		//savedImgSet.push_back(baseImg);
		if (count / 2 < grayCodeFboSet.size()) {
			savedGrayImgSet.push_back(baseImg);
		}
		else {
			savedNegaImgSet.push_back(baseImg);
		}
	}
	preCount = count;

	if (savedNegaImgSet.size() == grayCodeFboSet.size()) {
		isCalculateRun += 1;
		printf("a");
	}
	if (isCalculateRun == 1) {
		printf("b");
		ofImage img;
		img.allocate(camW, camH, OF_IMAGE_COLOR);
		for (int i = 0; i < savedGrayImgSet.size(); i++) {
			subtract(savedGrayImgSet[i].getPixels(), savedNegaImgSet[i].getPixels(), &img.getPixels());
			savedImgSet.push_back(img);
		}

		correspondFbo.begin();
		ofClear(0);
		shader.load("", "shader.frag");
		ofFill();
		shader.begin();
		shader.setUniform2f("u_resolution", proW, proH);
		shader.setUniform2f("cam_resolution", camW, camH);
		shader.setUniformTexture("img0", savedImgSet[0].getTextureReference(), 0);
		shader.setUniformTexture("img1", savedImgSet[1].getTextureReference(), 1);
		shader.setUniformTexture("img2", savedImgSet[2].getTextureReference(), 2);
		shader.setUniformTexture("img3", savedImgSet[3].getTextureReference(), 3);
		shader.setUniformTexture("img4", savedImgSet[4].getTextureReference(), 4);
		shader.setUniformTexture("img5", savedImgSet[5].getTextureReference(), 5);
		shader.setUniformTexture("img6", savedImgSet[6].getTextureReference(), 6);
		shader.setUniformTexture("img7", savedImgSet[7].getTextureReference(), 7);
		shader.setUniformTexture("img8", savedImgSet[8].getTextureReference(), 8);
		shader.setUniformTexture("img9", savedImgSet[9].getTextureReference(), 9);
		shader.setUniformTexture("img10", savedImgSet[10].getTextureReference(), 10);
		shader.setUniformTexture("img11", savedImgSet[11].getTextureReference(), 11);
		shader.setUniformTexture("img12", savedImgSet[12].getTextureReference(), 12);
		shader.setUniformTexture("img13", savedImgSet[13].getTextureReference(), 13);
		shader.setUniformTexture("img14", savedImgSet[14].getTextureReference(), 14);
		shader.setUniformTexture("img15", savedImgSet[15].getTextureReference(), 15);
		shader.setUniformTexture("img16", savedImgSet[16].getTextureReference(), 16);
		shader.setUniformTexture("img17", savedImgSet[17].getTextureReference(), 17);
		shader.setUniformTexture("img18", savedImgSet[18].getTextureReference(), 18);
		shader.setUniformTexture("img19", savedImgSet[19].getTextureReference(), 19);
		shader.setUniformTexture("img20", savedImgSet[20].getTextureReference(), 20);
		shader.setUniformTexture("img21", savedImgSet[21].getTextureReference(), 21);
		ofRect(0, 0, camW, camH);
		shader.end();
		correspondFbo.end();

		correspondFbo.readToPixels(correspondPixels);
		for (int i = 0; i < camW; i++) {
			for (int j = 0; j < camH; j++) {
				ofFloatColor color = correspondPixels.getColor(i, j);
				if (color.b < 0.5) {
					ofVec2f p = ofVec2f(int(color.r * proW), int(color.g * proH));
					camPosSet.push_back(ofVec2f(i, j));
					proPosSet.push_back(ofVec2f(p.x, p.y));
				}
			}
		}
	}
}
void StructuredLight::draw() {
	for (int i = 0; i < savedImgSet.size(); i++) {
		savedImgSet[i].draw((i / 9) * (proW / 10), (i % 9) * (proH / 10), proW / 10, proH / 10);
	}
	correspondFbo.draw(0, 0, camW / 5, camH / 5);
	correspondFbo.draw(camW/5, 0, camW / 5, camH / 5);

}

void StructuredLight::drawGrayCode() {
	if (isSetupFinished) {
		if (count / 2 < grayCodeFboSet.size()) {
			grayCodeFboSet[int(count / 2)].draw(0, 0);
		}
		else {
			grayCodeNegaFboSet[int(count / 2 - grayCodeFboSet.size())].draw(0, 0);
		}
	}
}

void StructuredLight::generateGrayCode(vector<int>* grayValue, int value, int num_bit) {
	if (value == 0) {
		for (int i = 0; i < num_bit; i++) {
			grayValue->push_back(0);
		}
	}
	else {
		vector<int> tempBinary;
		while (value > 0) {
			tempBinary.push_back(value % 2);
			value = value / 2;
		}
		for (int i = 0; i < num_bit - tempBinary.size(); i++) {
			grayValue->push_back(0);
		}
		grayValue->push_back(tempBinary[tempBinary.size() - 1]);
		for (int j = tempBinary.size() - 1; j > 0; j--) {
			int num = tempBinary[j] + tempBinary[j - 1];
			num = num % 2;
			grayValue->push_back(num);
		}
	}
}


