#ifndef SVM_SOFT_H_INCLUDED
#define SVM_SOFT_H_INCLUDED

#include "binary_class.h"

#include <algorithm>

#include <CGAL/MP_Float.h>
#include <CGAL/QP_functions.h>

class SVMSoft : public BinaryClass<double> {
  using BinaryClass<double>::dimension_;
  using BinaryClass<double>::vectors_;
  using BinaryClass<double>::line_;
protected:
  std::vector<double> xi_;
public:
  typedef double ValueType;

  SVMSoft(size_t dim) : BinaryClass<double>(dim) {}

  const std::vector<double>& GetXi() const { return xi_; }

  void Solve(double C) {
    if (C < 0) return;
    CGAL::Quadratic_program<double> qp(CGAL::LARGER, false, 0, false, 0);
    for (size_t i = 1; i <= dimension_; i++) qp.set_d(i, i, 1);

    for (size_t i = 0; i < vectors_.size(); i++) {
      qp.set_a(0, i, vectors_[i].second);
      for (size_t j = 1; j <= dimension_; j++)
        qp.set_a(j, i, vectors_[i].first[j] * vectors_[i].second);
      qp.set_b(i, 1);

      qp.set_a(dimension_ + 1 + i, i, 1);
      qp.set_c(dimension_ + 1 + i, C);
      qp.set_l(dimension_ + 1 + i, true, 0);
    }

    CGAL::Quadratic_program_solution<CGAL::MP_Float> sol =
        CGAL::solve_quadratic_program(qp, CGAL::MP_Float());

    auto it = sol.variable_values_begin();
    xi_.resize(vectors_.size());
    for (size_t i = 0; i <= dimension_; i++, it++)
      line_[i] = CGAL::to_double(*it);
    for (size_t i = 0; i < vectors_.size(); i++, it++)
      xi_[i] = CGAL::to_double(*it);
  }
};

#endif
