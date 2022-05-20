# g++ -std=c++17 main.cpp coords.cpp mesh.cpp parameterization.cpp -o main
# clang++ -std=c++17 -Xpreprocessor -fopenmp src/*.cpp -o main
clang++ -std=c++17 -fopenmp -lomp -I"$(brew --prefix libomp)/include" -L"$(brew --prefix libomp)/lib" src/*.cpp -o main
