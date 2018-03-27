#ifndef BINARY_CLASS_KERNEL_H_INCLUDED
#define BINARY_CLASS_KERNEL_H_INCLUDED

#include <valarray>
#include <vector>
#include <functional>

template <class T> class BinaryClassKernel {
public:
  typedef std::valarray<T> VectorType;
  typedef std::vector<std::pair<VectorType, int>> DataSetType;
  typedef std::function<T(const VectorType&, const VectorType&)> KernelType;
protected:
  size_t dimension_;
  DataSetType vectors_;
  KernelType kernel_;
  std::vector<std::pair<size_t, T>> alpha_;
  T bias_;
  int sgn(const T& a) const {
    return a > 0 ? 1 : -1;
  }
public:
  BinaryClassKernel(size_t dim, KernelType ker) :
      dimension_(dim), kernel_(ker) {}
  const size_t size() const { return vectors_.size(); }
  const DataSetType& GetDatas() const {
    return vectors_;
  }

  void SetKernel(KernelType ker) { kernel_ = ker; }
  KernelType GetKernel() const { return kernel_; }

  bool CheckPointIncorrect(size_t i) const {
    return Predict(vectors_[i].first) != vectors_[i].second;
  }
  size_t CountPointIncorrect() const {
    size_t ans = 0;
    for (auto& i : vectors_)
      if (Predict(i.first) != i.second) ans++;
    return ans;
  }

  const std::vector<std::pair<size_t, T>>& GetAlpha() const { return alpha_; }
  void SetAlpha(const std::vector<std::pair<size_t, T>>& a) { alpha_ = a; }
  const T& GetBias() const { return bias_; }
  void SetBias(const T& b) { bias_ = b; }

  void AddData(const VectorType& data, bool cls) {
    if (data.size() != dimension_) return;
    vectors_.emplace_back(data, cls ? 1 : -1);
  }
  int Predict(const VectorType& data) const {
    T res = bias_;
    for (const std::pair<size_t, T>& i : alpha_) {
      res += kernel_(vectors_[i.first].first, data) * vectors_[i.first].second
          * i.second;
    }
    return sgn(res);
  }
};

#endif
