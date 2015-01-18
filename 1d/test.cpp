#include "diamond.h"

const double DT = 0.01;

void uxxStep0(std::array<const LocalVariable*, 1> inputs,
               std::array<LocalVariable*, 2> outputs,
               LocalMesh& mesh) {
    const LocalVariable* u = inputs[0];
    outputs[0]->val(u->val());
    outputs[1]->val((u->nbr(0) + u->nbr(1) - 2 * u->val())
                  / (mesh.dx * mesh.dx));
}

void updateStep0(std::array<const LocalVariable*, 2> inputs,
                 std::array<LocalVariable*, 2> outputs,
                 LocalMesh& mesh) {
    const LocalVariable* u = inputs[0];
    double conv = (u->nbr(1) - u->nbr(0)) / (2 * mesh.dx);
    outputs[0]->val(u->val());

    const LocalVariable* uxx = inputs[1];
    double uPlusUxx = u->val() + uxx->val();
    double uPlusUxxL = u->nbr(0) + uxx->nbr(0);
    double uPlusUxxR = u->nbr(1) + uxx->nbr(1);
    double diff = (uPlusUxxL + uPlusUxxR - 2 * uPlusUxx)
                / (mesh.dx * mesh.dx);

    double dudt = -conv - diff;
    outputs[1]->val(u->val() + 0.5 * DT * dudt);
}

void uxxStep1(std::array<const LocalVariable*, 2> inputs,
               std::array<LocalVariable*, 3> outputs,
               LocalMesh& mesh) {
    const LocalVariable* u = inputs[0];
    const LocalVariable* uMid = inputs[1];
    outputs[0]->val(u->val());
    outputs[1]->val(uMid->val());
    outputs[2]->val((uMid->nbr(0) + uMid->nbr(1) - 2 * uMid->val())
                  / (mesh.dx * mesh.dx));
}

void updateStep1(std::array<const LocalVariable*, 3> inputs,
                 std::array<LocalVariable*, 1> outputs,
                 LocalMesh& mesh) {
    const LocalVariable* u = inputs[0];
    const LocalVariable* uMid = inputs[1];
    double conv = (uMid->nbr(1) - uMid->nbr(0)) / (2 * mesh.dx);

    const LocalVariable* uxx = inputs[2];
    double uPlusUxx = uMid->val() + uxx->val();
    double uPlusUxxL = uMid->nbr(0) + uxx->nbr(0);
    double uPlusUxxR = uMid->nbr(1) + uxx->nbr(1);
    double diff = (uPlusUxxL + uPlusUxxR - 2 * uPlusUxx)
                / (mesh.dx * mesh.dx);

    double dudt = -conv - diff;
    outputs[1]->val(u->val() + DT * dudt);
}

void init(std::array<LocalVariable*, 1> u, LocalMesh&) {
    u[0]->val(0.0);
}

int main()
{
    ClassicDiscretization1D disc(100, 0.01, init);
    for (int iStep = 0; iStep < 100; ++iStep) {
        disc.applyOp(uxxStep0);
        disc.applyOp(updateStep0);
        disc.applyOp(uxxStep1);
        disc.applyOp(updateStep1);
    }
}

