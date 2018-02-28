#ifndef PLA_H_INCLUDED
#define PLA_H_INCLUDED

#include "binary_class.h"

template <class T> class PLA : public BinaryClass<T> {
  using BinaryClass<T>::CheckPointIncorrect;
  using BinaryClass<T>::vectors_;
  using BinaryClass<T>::line_;
public:
  PLA(size_t dim) : BinaryClass<T>(dim) {}
  bool CheckAndUpdate(size_t i) {
    if (!CheckPointIncorrect(i)) return false;
    line_ += vectors_[i].first * (T)vectors_[i].second;
    return true;
  }
};

#endif // PLA_H_INCLUDED
