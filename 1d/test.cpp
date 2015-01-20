#include <cmath>
#include <ctime>
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
    double conv = (uR*uR - uL*uL) / (4 * mesh.dx);

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
    double conv = (uR*uR - uL*uL) / (4 * mesh.dx);

    double uPlusUxx = u + inputs.val(2);
    double uPlusUxxL = uL + inputs.nbr(2, 0);
    double uPlusUxxR = uR + inputs.nbr(2, 1);
    double diff = (uPlusUxxL + uPlusUxxR - 2 * uPlusUxx)
                / (mesh.dx * mesh.dx);

    double dudt = -conv - diff;
    outputs.val(0, u0 + 0.5 * DT * dudt);
}
void init(LocalVariables<1>& u, const LocalMesh& mesh) {
    const double PI = atan(1.0) * 4;
    u.val(0, cos(mesh.x / 128. * 19 * PI) * 2.);
}

int main()
{
    const int nStepsPerPixel = 1000, nPixel = 100;

    ClassicDiscretization1D disc(16, 0.5, init);
    disc.colorMap.red.set(0, -2., 2.);

    std::clock_t startTime = std::clock();

    for (int iPixel = 0; iPixel < nPixel; ++iPixel) {
        for (int iStep = 0; iStep < nStepsPerPixel; ++iStep) {
            disc.applyOp(uxxStep0);
            disc.applyOp(updateStep0);
            disc.applyOp(uxxStep1);
            disc.applyOp(updateStep1);
        }
        disc.variablesToColor(iPixel);
        disc.writePng("test");
    }

    std::clock_t endTime = std::clock();
    double totalTime = (endTime - startTime) / (double)CLOCKS_PER_SEC;
    std::cout << totalTime * 1000000 / nStepsPerPixel / nPixel / 4
              << " microseconds per SubStep" << std::endl;
}

