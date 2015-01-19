#include "diamond.h"

const double DT = 0.005;

void uxxStep0(const LocalVariables<1>& inputs,
              LocalVariables<2>& outputs,
              const LocalMesh& mesh) {
    double u = inputs.val(0);
    outputs.val(0, u);
    double uL = inputs.nbr(0, 0), uR = inputs.nbr(0, 1);
    double uxx = (uL + uR - 2 * u) / (mesh.dx * mesh.dx);
    outputs.val(1, uxx);
}

void updateStep0(const LocalVariables<2>& inputs,
                 LocalVariables<2>& outputs,
                 const LocalMesh& mesh) {
    double u = inputs.val(0);
    outputs.val(0, u);

    double uL = inputs.nbr(0, 0), uR = inputs.nbr(0, 1);
    double conv = (uR - uL) / (2 * mesh.dx);

    double uPlusUxx = u + inputs.val(1);
    double uPlusUxxL = uL + inputs.nbr(1, 0);
    double uPlusUxxR = uR + inputs.nbr(1, 1);
    double diff = (uPlusUxxL + uPlusUxxR - 2 * uPlusUxx)
                / (mesh.dx * mesh.dx);

    double dudt = -conv - diff;
    outputs.val(1, u + 0.5 * DT * dudt);
}

void uxxStep1(const LocalVariables<2>& inputs,
              LocalVariables<3>& outputs,
              const LocalMesh& mesh) {
    double u0 = inputs.val(0);
    outputs.val(0, u0);

    double u = inputs.val(1);
    outputs.val(1, u);

    double uL = inputs.nbr(1, 0), uR = inputs.nbr(1, 1);
    double uxx = (uL + uR - 2 * u) / (mesh.dx * mesh.dx);
    outputs.val(2, uxx);
}

void updateStep1(const LocalVariables<3>& inputs,
                 LocalVariables<1>& outputs,
                 const LocalMesh& mesh) {
    double u0 = inputs.val(0);

    double u = inputs.val(1);
    double uL = inputs.nbr(1, 0), uR = inputs.nbr(1, 1);
    double conv = (uR - uL) / (2 * mesh.dx);

    double uPlusUxx = u + inputs.val(2);
    double uPlusUxxL = uL + inputs.nbr(2, 0);
    double uPlusUxxR = uR + inputs.nbr(2, 1);
    double diff = (uPlusUxxL + uPlusUxxR - 2 * uPlusUxx)
                / (mesh.dx * mesh.dx);

    double dudt = -conv - diff;
    outputs.val(0, u0 + 0.5 * DT * dudt);
}

void init(LocalVariables<1>& u, const LocalMesh& mesh) {
    u.val(0, mesh.x);
}

int main()
{
    ClassicDiscretization1D disc(128 * 2, 0.5, init);
    disc.colorMap.red.set(0, 0., 128.);
    for (int iPng = 0; iPng < 200; ++iPng) {
        for (int iStep = 0; iStep < 500; ++iStep) {
            disc.applyOp(uxxStep0);
            disc.applyOp(updateStep0);
            disc.applyOp(uxxStep1);
            disc.applyOp(updateStep1);
        }
        disc.variablesToColor(iPng);
        disc.writePng("test.png");
    }
}

