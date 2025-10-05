#version 430
const int TRAIL_LENGTH = 64; // must match compute shader

struct Photon {
    vec2 s;        // (r, phi)
    float L;
    float h;
    int alive;
    float freqShift;
    vec2 trail[TRAIL_LENGTH];
    int trail_head; // index of next write position (newest is trail_head-1)
};

layout(std430, binding = 0) buffer PhotonBuffer {
    Photon photons[];
};

uniform mat4 proj;
uniform mat4 view;
uniform float zoom;

// Outputs to fragment shader
out vec3 vColor;
out float vAlpha;

void main() {
    uint photonID = uint(gl_InstanceID);
    uint pointID  = uint(gl_VertexID); // 0 = newest after remap

    Photon p = photons[photonID];

    // Cull dead photons by moving off-screen with zero alpha
    if (p.alive == 0) {
        gl_Position = vec4(2.0, 2.0, 0.0, 1.0);
        vColor = vec3(0.0);
        vAlpha = 0.0;
        return;
    }

    // Remap so that pointID 0 corresponds to the newest sample (trail_head-1)
    uint newest = uint((p.trail_head - 1 + TRAIL_LENGTH) % TRAIL_LENGTH);
    uint idx = (newest + pointID) % uint(TRAIL_LENGTH);
    vec2 pos = p.trail[idx]; // already Cartesian (x, y)
    gl_Position = proj * view * vec4(pos * zoom, 0.0, 1.0);
    gl_PointSize = 3.0;

    // Simple Doppler-based color from freqShift
    float shift = clamp(p.freqShift, 0.0, 2.0);
    vColor = vec3(shift, 0.2, 1.0 - shift * 0.5);

    // Fading alpha: newer = brighter; head is fully opaque
    vAlpha = clamp(1.0 - float(pointID) / float(TRAIL_LENGTH), 0.05, 1.0);
    if (pointID == 0u) {
        vAlpha = 1.0;
    }
}
