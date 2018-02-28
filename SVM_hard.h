#ifndef SVM_HARD_H_INCLUDED
#define SVM_HARD_H_INCLUDED

#include "binary_class.h"

#include <CGAL/MP_Float.h>
#include <CGAL/QP_functions.h>

class SVMHard : public BinaryClass<double> {
  using BinaryClass<double>::CheckPointIncorrect;
  using BinaryClass<double>::vectors_;
  using BinaryClass<double>::line_;
public:
  typedef double ValueType;

  SVMHard(size_t dim) : BinaryClass<double>(dim) {}

  bool Solve() {
    typedef CGAL::Quadratic_program<double> QProgram;
    typedef CGAL::Quadratic_program_solution<CGAL::MP_Float> Solution;
    QProgram qp(CGAL::LARGER, false, 0, false, 0);
    for (size_t i = 1; i <= dimension_; i++) qp.set_d(i, i, 1);
    for (size_t i = 0; i < vectors_.size(); i++) {
      qp.set_a(0, i, vectors_[i].second);
      for (size_t j = 1; j <= dimension_; j++)
        qp.set_a(j, i, vectors_[i].first[j] * vectors_[i].second);
      qp.set_b(i, 1);
    }

    Solution sol = CGAL::solve_quadratic_program(qp, CGAL::MP_Float());
    if (sol.is_infeasible()) return false;

    auto it = sol.variable_values_begin();
    for (size_t i = 0; i <= dimension_; i++, it++)
      line_[i] = CGAL::to_double(*it);
    return true;
  }
};

#endif
