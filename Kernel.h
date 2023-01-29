#ifndef KERNEL_H
#define KERNEL_H

/// Abstract Kernel class
class Kernel {
public:
    /**
     * Return weight at OFFSET amino acids away from the boundary
     * @param offset Distance from the boundary (in amino acids).
     *               Offset coordinates are zero-based
     * @return kernel weight
     */
    virtual double getWeight(int offset) = 0;
    /**
     * Set scoring window width
     */
    virtual void setWidth(int width);
    /**
     * Sum of all kernel weights within a window, are under kernel
     */
    virtual double weightSum();
    virtual ~Kernel() {}
protected:
    int width;
};

/// Box Kernel
class BoxKernel : public Kernel {
public:
    double getWeight(int offset);
};

/// Triangular Kernel
class TriangularKernel : public Kernel {
public:
    double getWeight(int offset);
};

/// Parabolic Kernel
class ParabolicKernel : public Kernel {
public:
    double getWeight(int offset);
};

/// Triweight Kernel
class TriweightKernel : public Kernel {
public:
    double getWeight(int offset);
};

#endif /* KERNEL_H */
