
uniform sampler2DRect map;
uniform sampler2DRect img;

uniform vec2 u_resolution;
uniform vec2 cam_resolution;

void main() {
    vec4 color;
    vec2 st = gl_FragCoord.xy;
    st.y = u_resolution.y - st.y;
    vec2 pos = texture(map, st).xy;
    pos.x = pos.x* u_resolution.x;
    //pos.y = 1.0 - pos.y;
    pos.y = pos.y * u_resolution.y;
    float c = texture(img, pos).x;
    color = vec4(c, c, c, 1.0);
    gl_FragColor = color;
}

