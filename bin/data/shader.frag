
uniform sampler2DRect img0;
uniform sampler2DRect img1;
uniform sampler2DRect img2;
uniform sampler2DRect img3;
uniform sampler2DRect img4;
uniform sampler2DRect img5;
uniform sampler2DRect img6;
uniform sampler2DRect img7;
uniform sampler2DRect img8;
uniform sampler2DRect img9;
uniform sampler2DRect img10;
uniform sampler2DRect img11;
uniform sampler2DRect img12;
uniform sampler2DRect img13;
uniform sampler2DRect img14;
uniform sampler2DRect img15;
uniform sampler2DRect img16;
uniform sampler2DRect img17;
uniform sampler2DRect img18;
uniform sampler2DRect img19;
uniform sampler2DRect img20;
uniform sampler2DRect img21;

uniform vec2 u_resolution;
uniform vec2 cam_resolution;

void main() {
    vec4 color;
    vec2 st = gl_FragCoord.xy;

    float g0 = texture(img0, st).x;
    float g1 = texture(img1, st).x;
    float g2 = texture(img2, st).x;
    float g3 = texture(img3, st).x;
    float g4 = texture(img4, st).x;
    float g5 = texture(img5, st).x;
    float g6 = texture(img6, st).x;
    float g7 = texture(img7, st).x;
    float g8 = texture(img8, st).x;
    float g9 = texture(img9, st).x;
    float g10 = texture(img10, st).x;

    float g11 = texture(img11, st).x;
    float g12 = texture(img12, st).x;
    float g13 = texture(img13, st).x;
    float g14 = texture(img14, st).x;
    float g15 = texture(img15, st).x;
    float g16 = texture(img16, st).x;
    float g17 = texture(img17, st).x;
    float g18 = texture(img18, st).x;
    float g19 = texture(img19, st).x;
    float g20 = texture(img20, st).x;
    float g21 = texture(img21, st).x;


    /*x*/
    float[11] gx = { g0, g1, g2, g3, g4, g5, g6, g7, g8, g9, g10 };
    int grayCodeX[11];
    for(int i=0; i<11; i++)
    {
        int codeValue=1.0;
        if (gx[i] < 0.1) { codeValue = 0; }
        grayCodeX[i] = codeValue;
    }
    int binaryCodeX[11];
    binaryCodeX[0] = grayCodeX[0];
    for (int i=1; i<11; i++)
    {
        binaryCodeX[i] = (grayCodeX[i] + binaryCodeX[i - 1]) % 2;
    }

    int valueX = 0;
    for(int i=0; i<11; i++)
    {
        valueX += binaryCodeX[i]*pow(2, 10 - i);
    }
    /*x*/

    /*y*/
    float[10] gy = { g11, g12, g13, g14, g15, g16, g17, g18, g19, g20 };
    int grayCodeY[10];
    for (int i = 0; i < 10; i++)
    {
        int codeValue = 1.0;
        if (gy[i] < 0.1) { codeValue = 0; }
        grayCodeY[i] = codeValue;
    }
    int binaryCodeY[10];
    binaryCodeY[0] = grayCodeY[0];
    for (int i = 1; i < 10; i++)
    {
        binaryCodeY[i] = (grayCodeY[i] + binaryCodeY[i - 1]) % 2;
    }

    int valueY = 0;
    for (int i = 0; i < 10; i++)
    {
        valueY += binaryCodeY[i] * pow(2, 9 - i);
    }
    /*y*/

    color = vec4(valueX / u_resolution.x, valueY / u_resolution.y, 0.0, 1.0);

    if (g21 < 0.3)
    {
        color.z += 1.0;
    }

    gl_FragColor = color;
}

