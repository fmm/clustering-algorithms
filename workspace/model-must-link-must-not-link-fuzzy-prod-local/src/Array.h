#ifndef ARRAY_H_
#define ARRAY_H_

#include "Includes.h"

template<class __Tp>
struct Array : public vector<__Tp> {
  Array(size_t __n = 0, __Tp value = __Tp()) : vector<__Tp>(__n,value) {
  }
  __Tp& operator[](size_t __n) {
    assert(__n < this->size());
    return vector<__Tp>::operator[](__n);
  }
  const __Tp& operator[](size_t __n) const {
    assert(__n < this->size());
    return vector<__Tp>::operator[](__n);
  }
};

#endif /* ARRAY_H_ */
