#version 430

in vec3 vColor;
in float vAlpha;

out vec4 FragColor;

void main() {
    FragColor = vec4(vColor, clamp(vAlpha, 0.0, 1.0));
}
