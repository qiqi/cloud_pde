#include "diamond.h"

const double DT = 0.005;
const int nStepsPerPixel = 1000, nPixel = 10;

void uxxStep0(const LocalInputs<1>& inputs,
              LocalOutputs<2>& outputs,
              const LocalMesh& mesh) {
    double u = inputs.val(0);
    outputs.val(0, u);
    double uL = inputs.nbr(0, 0), uR = inputs.nbr(0, 1);
    double uxx = (uL + uR - 2 * u) / (mesh.dx * mesh.dx);
    outputs.val(1, uxx);
}

void updateStep0(const LocalInputs<2>& inputs,
                 LocalOutputs<2>& outputs,
                 const LocalMesh& mesh) {
    double u = inputs.val(0);
    outputs.val(0, u);

    double uL = inputs.nbr(0, 0), uR = inputs.nbr(0, 1);
    double conv = (uR*uR - uL*uL) / (4 * mesh.dx);

    double uPlusUxx = u + inputs.val(1);
    double uPlusUxxL = uL + inputs.nbr(1, 0);
    double uPlusUxxR = uR + inputs.nbr(1, 1);
    double diff = (uPlusUxxL + uPlusUxxR - 2 * uPlusUxx)
                / (mesh.dx * mesh.dx);

    double dudt = -conv - diff;
    outputs.val(1, u + 0.5 * DT * dudt);
}

void uxxStep1(const LocalInputs<2>& inputs,
              LocalOutputs<3>& outputs,
              const LocalMesh& mesh) {
    double u0 = inputs.val(0);
    outputs.val(0, u0);

    double u = inputs.val(1);
    outputs.val(1, u);

    double uL = inputs.nbr(1, 0), uR = inputs.nbr(1, 1);
    double uxx = (uL + uR - 2 * u) / (mesh.dx * mesh.dx);
    outputs.val(2, uxx);
}

void updateStep1(const LocalInputs<3>& inputs,
                 LocalOutputs<1>& outputs,
                 const LocalMesh& mesh) {
    double u0 = inputs.val(0);

    double u = inputs.val(1);
    double uL = inputs.nbr(1, 0), uR = inputs.nbr(1, 1);
    double conv = (uR*uR - uL*uL) / (4 * mesh.dx);

    double uPlusUxx = u + inputs.val(2);
    double uPlusUxxL = uL + inputs.nbr(2, 0);
    double uPlusUxxR = uR + inputs.nbr(2, 1);
    double diff = (uPlusUxxL + uPlusUxxR - 2 * uPlusUxx)
                / (mesh.dx * mesh.dx);

    double dudt = -conv - diff;
    outputs.val(0, u0 + 0.5 * DT * dudt);
}
void init(LocalOutputs<1>& u, const LocalMesh& mesh) {
    const double PI = atan(1.0) * 4;
    u.val(0, cos(mesh.x / 128. * 19 * PI) * 2.);
}
