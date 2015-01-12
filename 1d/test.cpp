#include "test.h"

void laplacian(std::array<const LocalVariable, 1> inputs,
               std::array<LocalVariable, 1> outputs,
               LocalMesh& mesh) {
    double u = inputs[0].val();
    double uL = inputs[0].nbr(0);
    double uR = inputs[0].nbr(1);
    outputs[0].val((uL + uR - 2 * u) / (mesh.dx * mesh.dx));
}

void derivative(std::array<const LocalVariable, 1> inputs,
                std::array<LocalVariable, 1> outputs,
                LocalMesh& mesh) {
    double uL = inputs[0].nbr(0);
    double uR = inputs[0].nbr(1);
    outputs[0].val((uR - uL) / (2 * mesh.dx));
}

Variable calcDudt(Variable u) {
    Variable uxx, diff, conv;

    Operator(laplacian) >> u << uxx;
    Operator(laplacian) >> (u + uxx) << diff;
    Operator(derivative) >> (0.5 * u * u) << conv;
    return -conv - diff;
}

int main()
{
    double dt = 0.01;
    Variable u;
    for (int iSteps = 0; iSteps < 100; ++ iSteps) {
        Variable u2 = u + 0.5* dt * calcDudt(u);
        u += dt * calcDudt(u2);
    }
}

