# black_hole
Black hole simulation project
Here is the black hole raw code, everything will be inside a src bin incase you want to copy the files
I'm writing this as I'm beginning this project (hopefully I complete it ;D) here is what I plan to do:

1. Ray-tracing : add ray tracing to the gravity simulation to simulate gravitational lensing
2. Accretion disk : simulate accreciate disk using the ray tracing + the halos
3. Spacetime curvature : demonstrate visually the "trapdoor in spacetime" that is black holes using spacetime grid
4. [optional] try to make it run realtime ;D

I hope it works :/



Edit: After completion of project - 

Thank you everyone for checking out the video, if you haven't it explains code in detail: https://www.youtube.com/watch?v=8-B6ryuBkCM

Quickly to run the file you should install opengl (GLFW) and OpenGL wrangler GLEW: https://glew.sourceforge.net/
I'm on windows and used msys2 for installation though.

How the code works:

for 2D: simple, just run 2D_lensing.cpp with the nessesary dependencies installed.

for 3D: black_hole.cpp and geodesic.comp work together to run the simuation faster using GPU, essentially it sends over a UBO and geodesic.comp runs heavy calculations using that data.
should work with nessesary dependencies installed, however I have only run it on windows with my GPU so am not sure!

LMK if you would like an in-depth explanation of how the code works aswell :)
