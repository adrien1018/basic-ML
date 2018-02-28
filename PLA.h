#ifndef PLA_H_INCLUDED
#define PLA_H_INCLUDED

#include <valarray>
#include <vector>

template <class T> class PLA2 {
private:
  size_t Dimension_;
  std::vector<std::pair<std::valarray<T>, int>> Vectors_;
  std::valarray<T> NowLine_;
  int sgn(const T& a) const {
    return a > 0 ? 1 : -1;
  }
public:
  PLA2(size_t dim) : Dimension_(dim), NowLine_(T(0), dim + 1) {}
  const size_t size() const { return Vectors_.size(); }
  const std::vector<std::pair<std::valarray<T>, int>>& GetDatas() const {
    return Vectors_;
  }
  const std::valarray<T>& GetLine() const { return NowLine_; }
  void SetLine(const std::valarray<T>& ln) {
    if (ln.size() != Dimension_ + 1) return;
    NowLine_ = ln;
  }
  bool CheckPointIncorrect(size_t i) const {
    return sgn((Vectors_[i].first * NowLine_).sum()) != Vectors_[i].second;
  }
  size_t CountPointIncorrect() const {
    size_t ans = 0;
    for (auto& i : Vectors_)
      if (sgn((i.first * NowLine_).sum()) != i.second) ans++;
    return ans;
  }
  bool CheckAndUpdate(size_t i) {
    if (!CheckPointIncorrect(i)) return false;
    NowLine_ += Vectors_[i].first * (T)Vectors_[i].second;
    return true;
  }
  void AddData(const std::valarray<T>& data, bool cls) {
    if (data.size() != Dimension_) return;
    std::valarray<T> v(Dimension_ + 1);
    v[0] = 1;
    v[std::slice(1, Dimension_, 1)] = data;
    Vectors_.emplace_back(std::move(v), cls ? 1 : -1);
  }
};

#endif // PLA_H_INCLUDED
