#ifndef MultiindexSet_h
#define MultiindexSet_h

#include <numeric> // std::accumulate
#include <vector>  // std::vector

/**
 * Multiindex set base class
 */
class MultiindexSet {
  protected:

    int dim, order;

    std::vector<std::vector<int>> indices;

  public:

    MultiindexSet(const int dim, const int order) : dim(dim), order(order) { }

    virtual void        ConstructRecursively(std::vector<int>& index, const int level) = 0;

    inline unsigned int size() const
    {
      return indices.size();
    }

    inline std::vector<int> GetMultiindex(const int i) const
    {
      return indices[i];
    }
};

/**
 * Full tensor multiindex set
 */
class FullTensorMultiindexSet : public MultiindexSet {
  public:

    inline void ConstructRecursively(std::vector<int>& index, const int level) override
    {
      if (level == dim - 1) {
        for (int p = 0; p <= order; ++p) {
          index[level] = p;
          indices.push_back(index);
        }
      } else {
        for (int p = 0; p <= order; ++p) {
          index[level] = p;
          ConstructRecursively(index, level + 1);
        }
      }
    }

    FullTensorMultiindexSet(const int dim, const int order) : MultiindexSet(dim, order)
    {
      std::vector<int> index(dim, 0);
      ConstructRecursively(index, 0);
    }
};

/**
 * Total-order multiindex set, i.e., L1 norm of multiindex does not exceed maxOrder
 */
class TotalOrderMultiindexSet : public MultiindexSet {
  public:

    inline void ConstructRecursively(std::vector<int>& index, const int level) override
    {
      if (level == dim - 1) {
        for (int p = 0; p <= order; ++p) {
          index[level] = p;
          if (std::accumulate(index.begin(), index.end(), 0) > order) break;
          indices.push_back(index);
        }
      } else {
        for (int p = 0; p <= order; ++p) {
          index[level] = p;
          ConstructRecursively(index, level + 1);
        }
      }
    }

    TotalOrderMultiindexSet(const int dim, const int order) : MultiindexSet(dim, order)
    {
      std::vector<int> index(dim, 0);
      ConstructRecursively(index, 0);
    }
};

#endif // ifndef _MultiindexSet_h
