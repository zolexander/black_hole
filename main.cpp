#include "engine.hpp"

int main() {
    BlackholeSim::Engine engine;
    if (!engine.initGL("Photon Trails Demo", 1200, 800)) {
        return -1;
    }
    engine.run();
    return 0;
}