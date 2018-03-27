#ifndef BINARY_CLASS_H_INCLUDED
#define BINARY_CLASS_H_INCLUDED

#include <valarray>
#include <vector>

template <class T> class BinaryClass {
public:
  typedef std::valarray<T> VectorType;
  typedef std::vector<std::pair<VectorType, int>> DataSetType;
protected:
  size_t dimension_;
  DataSetType vectors_;
  VectorType line_;
  int sgn(const T& a) const {
    return a > 0 ? 1 : -1;
  }
public:
  BinaryClass(size_t dim) : dimension_(dim), line_(T(0), dim + 1) {}
  const size_t size() const { return vectors_.size(); }
  const DataSetType& GetDatas() const { return vectors_; }
  const VectorType& GetLine() const { return line_; }
  void SetLine(const VectorType& ln) {
    if (ln.size() != dimension_ + 1) return;
    line_ = ln;
  }

  bool CheckPointIncorrect(size_t i) const {
    return sgn((vectors_[i].first * line_).sum()) != vectors_[i].second;
  }
  size_t CountPointIncorrect() const {
    size_t ans = 0;
    for (auto& i : vectors_)
      if (sgn((i.first * line_).sum()) != i.second) ans++;
    return ans;
  }

  void AddData(const VectorType& data, bool cls) {
    if (data.size() != dimension_) return;
    VectorType v(dimension_ + 1);
    v[0] = 1;
    v[std::slice(1, dimension_, 1)] = data;
    vectors_.emplace_back(std::move(v), cls ? 1 : -1);
  }
  bool Predict(const VectorType& data) const {
    return (data * line_).sum() > 0;
  }
};

#endif
