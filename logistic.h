#ifndef LOGISTIC_H_INCLUDED
#define LOGISTIC_H_INCLUDED

#include "binary_class.h"
#include <cmath>

template <class T> class Logistic : public BinaryClass<T> {
  using BinaryClass<T>::dimension_;
  using BinaryClass<T>::vectors_;
  using BinaryClass<T>::line_;
private:
  T theta(const T& a) const {
    return 1 / (1 + exp(-a));
  }
public:
  Logistic(size_t dim) : BinaryClass<T>(dim) {}

  std::valarray<T> GetNegPartialGradient(size_t i) const {
    return theta((vectors_[i].first * line_).sum() * (T)-vectors_[i].second) *
          (vectors_[i].first * (T)vectors_[i].second);
  }
  std::valarray<T> GetNegGradient() const {
    std::valarray<T> v(dimension_ + 1);
    for (size_t i = 0; i < vectors_.size(); i++) v += GetNegPartialGradient(i);
    return v / (T)vectors_.size();
  }

  void UpdateSGD(size_t i, T eta) {
    line_ += eta * GetNegPartialGradient(i);
  }
  void UpdateGD(T eta) {
    line_ += eta * GetNegGradient();
  }
};

#endif // LOGISTIC_H_INCLUDED
