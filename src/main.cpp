#include <memory>
#include "parameterization.h"

using namespace RenderSpace;
using namespace std;

int main() {
    auto uns_mesh = make_shared<Mesh>();
    auto param_mesh = make_shared<Mesh>();
    auto str_mesh = make_shared<Mesh>();

    uns_mesh->load_OBJ("models/uns.obj");

    Parameterization pmethod(uns_mesh, param_mesh, str_mesh);
    float progress = 0.0f;
    pmethod.parameterize(ParamMethod::Laplace, progress, 500);

    param_mesh->save_OBJ("models/param.obj");
    str_mesh->save_OBJ("models/str.obj");
    return 0;
}
