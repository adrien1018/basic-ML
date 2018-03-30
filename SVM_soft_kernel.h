#ifndef SVM_SOFT_KERNEL_H_INCLUDED
#define SVM_SOFT_KERNEL_H_INCLUDED

#include "binary_class_kernel.h"

#include <cmath>
#include <CGAL/MP_Float.h>
#include <CGAL/QP_functions.h>

class SVMSoftKernel : public BinaryClassKernel<double> {
  using BinaryClassKernel<double>::dimension_;
  using BinaryClassKernel<double>::vectors_;
  using BinaryClassKernel<double>::kernel_;
  using BinaryClassKernel<double>::alpha_;
  using BinaryClassKernel<double>::bias_;
public:
  typedef double ValueType;

  SVMSoftKernel(size_t dim, KernelType ker) :
      BinaryClassKernel<double>(dim, ker) {}

  bool Solve(double C, double eps = 1. / 1048576) {
    CGAL::Quadratic_program<double> qp(CGAL::EQUAL, true, 0, true, C);
    for (size_t i = 0; i < vectors_.size(); i++) qp.set_c(i, -1);
    for (size_t i = 0; i < vectors_.size(); i++) {
      for (size_t j = 0; j <= i; j++)
        qp.set_d(i, j, kernel_(vectors_[i].first, vectors_[j].first) *
                        (vectors_[i].second * vectors_[j].second)),
      qp.set_a(i, 0, vectors_[i].second);
    }
    qp.set_b(0, 0);

    CGAL::Quadratic_program_solution<CGAL::MP_Float> sol =
        CGAL::solve_nonnegative_quadratic_program(qp, CGAL::MP_Float());
    if (sol.is_infeasible()) return false;

    alpha_.clear();
    size_t sv = 0; double mn = C;
    auto it = sol.variable_values_begin();
    for (size_t i = 0; i < vectors_.size(); i++, it++) {
      double val = CGAL::to_double(*it);
      if (fabs(val) > eps) alpha_.emplace_back(i, val);
      if (fabs(val - C / 2) < mn) mn = fabs(val - C / 2), sv = i;
    }
    bias_ = vectors_[sv].second;
    for (std::pair<size_t, double>& i : alpha_)
      bias_ -= kernel_(vectors_[sv].first, vectors_[i.first].first) *
          (vectors_[i.first].second * i.second);
    return true;
  }
};

#endif
