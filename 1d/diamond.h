#include<cassert>
#include<vector>
#include<array>

class LocalVariable {
    public:
    virtual double val() const = 0;
    virtual void val(double newVal) = 0;
    virtual double nbr(int i) const = 0;
};

class LocalMesh {
    public:
    double dx;
};

// sub-step

template<size_t numInput, size_t numOutput>
class SubStep {
    public:
    typedef void (&LocalOperator)(
         std::array<const LocalVariable*, numInput> inputs,
         LocalVariable* output, LocalMesh& mesh);

    SubStep (LocalOperator localOp);

    private:
    LocalOperator _localOp;
};

class Mesh {
    public:
    template<size_t numVar>
    Mesh(int numGrids, double dx,
         void (&LocalOperator)(
               std::array<LocalVariable*, numVar> inputs, LocalMesh& mesh));

    template<size_t numInput, size_t numOutput>
    void take(void (&LocalOperator)(
                std::array<const LocalVariable*, numInput> inputs,
                std::array<LocalVariable*, numOutput> outputs,
                LocalMesh& mesh));
};

