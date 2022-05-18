#include "parameterization.h"

using namespace RenderSpace;

int main() {
    Mesh ori;
    Mesh tar;

    ori.load_OBJ("models/model.obj");

    Parameterization param(&ori, &tar);
    param.parameterize();

    tar.save_OBJ("models/model_param.obj");
    return 0;
}