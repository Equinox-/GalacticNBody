export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:/opt/AMDAPP/lib/x86_64"
gcc -c -I/opt/AMDAPP/include/ -I/opt/AMDAPP/include/SDKUtil/ -I./ NBody.cpp -o NBody.o
g++ NBody.o -lOpenCL -lGLU -lGL -lglut -o NBody
