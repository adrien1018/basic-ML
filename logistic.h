#ifndef LOGISTIC_H_INCLUDED
#define LOGISTIC_H_INCLUDED

#include <valarray>
#include <vector>
#include <cmath>

template <class T> class Logistic {
private:
	size_t Dimension_;
	std::vector<std::pair<std::valarray<T>, int>> Vectors_;
	std::valarray<T> NowLine_;
	int sgn(const T& a) const {
		return a > 0 ? 1 : -1;
	}
	T theta(const T& a) const {
		return 1 / (1 + exp(-a));
	}
public:
	Logistic(size_t dim) : Dimension_(dim), NowLine_(T(0), dim + 1) {}
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

	std::valarray<T> GetNegPartialGradient(size_t i) const {
		return theta((Vectors_[i].first * NowLine_).sum() * (T)-Vectors_[i].second) *
					(Vectors_[i].first * (T)Vectors_[i].second);
	}
	std::valarray<T> GetNegGradient() const {
		std::valarray<T> v(Dimension_ + 1);
		for (size_t i = 0; i < Vectors_.size(); i++) v += GetNegPartialGradient(i);
		return v / (T)Vectors_.size();
	}

	void UpdateSGD(size_t i, T eta) {
		NowLine_ +=	eta * GetNegPartialGradient(i);
	}
	void UpdateGD(T eta) {
		NowLine_ +=	eta * GetNegGradient();
	}

	void AddData(const std::valarray<T>& data, bool cls) {
		if (data.size() != Dimension_) return;
		std::valarray<T> v(Dimension_ + 1);
		v[0] = 1;
		v[std::slice(1, Dimension_, 1)] = data;
		Vectors_.emplace_back(std::move(v), cls ? 1 : -1);
	}
};

#endif // LOGISTIC_H_INCLUDED
