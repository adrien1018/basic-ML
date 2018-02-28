#ifndef SVM_HARD_H_INCLUDED
#define SVM_HARD_H_INCLUDED

#include <valarray>
#include <vector>

#include <CGAL/MP_Float.h>
#include <CGAL/QP_functions.h>

class SVMHard {
public:
  typedef double ValueType;
private:
  size_t dimension_;
  std::vector<std::pair<std::valarray<ValueType>, int>> vectors_;
  std::valarray<ValueType> answer_;
  int sgn(const ValueType& a) const {
    return a > 0 ? 1 : -1;
  }
public:
  SVMHard(size_t dim) : dimension_(dim), answer_(0., dim + 1) {}

  const size_t size() const { return vectors_.size(); }

  const std::vector<std::pair<std::valarray<ValueType>, int>>&
  GetDatas() const {
    return vectors_;
  }

  const std::valarray<ValueType>& GetLine() const { return answer_; }

  void SetLine(const std::valarray<ValueType>& ln) {
    if (ln.size() != dimension_ + 1) return;
    answer_ = ln;
  }

  bool CheckPointIncorrect(size_t i) const {
    return sgn((vectors_[i].first * answer_).sum()) != vectors_[i].second;
  }

  size_t CountPointIncorrect() const {
    size_t ans = 0;
    for (auto& i : vectors_)
      if (sgn((i.first * answer_).sum()) != i.second) ans++;
    return ans;
  }

  void AddData(const std::valarray<ValueType>& data, bool cls) {
    if (data.size() != dimension_) return;
    std::valarray<ValueType> v(dimension_ + 1);
    v[0] = 1;
    v[std::slice(1, dimension_, 1)] = data;
    vectors_.emplace_back(std::move(v), cls ? 1 : -1);
  }

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
      answer_[i] = CGAL::to_double(*it);
    return true;
  }
};

#endif
