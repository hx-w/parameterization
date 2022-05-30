if [ "$(uname)" = "Darwin" ]; then
    clang++ -std=c++17 -fopenmp -lomp -I"$(brew --prefix libomp)/include" -L"$(brew --prefix libomp)/lib" src/*.cpp -o main
elif [ "$(expr substr $(uname -s) 1 5)" = "Linux" ];then
    c++ -fopenmp -std=c++17 ./src/*.cpp -o main
elif [ "$(expr substr $(uname -s) 1 10)" = "MINGW32_NT" ];then    
    c++ -fopenmp -std=c++17 .\\src\\*.cpp -o main.exe
fi
