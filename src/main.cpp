#include "parameterization.h"

using namespace RenderSpace;

int main() {
    Mesh uns_mesh;
    Mesh param_mesh;
    Mesh str_mesh;

    uns_mesh.load_OBJ("models/uns.obj");

    Parameterization pmethod(&uns_mesh, &param_mesh, &str_mesh);
    pmethod.parameterize();
    pmethod.resample(100);

    param_mesh.save_OBJ("models/param.obj");
    str_mesh.save_OBJ("models/str.obj");
    return 0;
}