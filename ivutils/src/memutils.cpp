// функции по работе с памятью
#include "memutils.h"
#include <vector>

/// Пара: настоящий_адрес, выровненный_адрес.
struct aligned_ptr_t
{ void* aligned;
  void* original;
  aligned_ptr_t (void* _original, void* _aligned) : original(_original), aligned(_aligned) {};
};

/// Здесь хранится массив с парами: настоящий_адрес, выровненный_адрес.
std::vector<aligned_ptr_t> aligned_ptr_list;

/// Функция, которая записывает во внутренний массив aligned_ptr_list
/// пару: настоящий_адрес, выровненный_адрес.
/// (эта функция не привязана к типу данных)
void aligned_ptr_list_push (void* original, void* aligned)
{ aligned_ptr_list.push_back ( aligned_ptr_t (original, aligned) );
}

/// Функция, которая ищет пару : настоящий_адрес, выровненный_адрес
/// во внутреннем массиве aligned_ptr_list.
/// (эта функция не привязана к типу данных)
void* aligned_ptr_list_pop (void* Value)
{ void* ret = 0;
  // перебрать весь массив с парами: настоящий_адрес, выровненный_адрес
  for (std::vector<aligned_ptr_t>::iterator it = aligned_ptr_list.begin();
       it != aligned_ptr_list.end(); it++)
    // если нашлась пара, где выровненный указатель равен указанному
    if ( (*it).aligned == Value )
    { // вернуть настоящий адрес
      ret = (*it).original;
      // удалить из массива запись
      aligned_ptr_list.erase(it);
      // выйти из цикла
      break;
    }
  // если соответствующая запись не была найдена, будет возвращен 0
  return ret;
}


#include "stdio.h"

template<>
void CompareMem (double one, double two, size_t size, size_t i)
{ printf("CompareMem [%lu] [%lu] %g != %g\n", size, i, one, two); 
}

template<>
void CompareMem (float one, float two, size_t size, size_t i)
{ printf("CompareMem [%lu] [%lu] %g != %g\n", size, i, one, two); 
}

template<>
void CompareMem (int one, int two, size_t size, size_t i)
{ printf("CompareMem [%lu] [%lu] %i != %i\n", size, i, one, two); 
}

template<>
void CompareMem (media_t<double> one, media_t<double> two, size_t size, size_t i)
{ printf("CompareMem(): media size = %lu, i = %lu\r\n", size, i); 
}

template<>
void CompareMem (media_t<float> one, media_t<float> two, size_t size, size_t i)
{ printf("CompareMem(): media size = %lu, i = %lu\r\n", size, i); 
}

template<>
void CompareMem (media2_t<double> one, media2_t<double> two, size_t size, size_t i)
{ printf("CompareMem(): media2 size = %lu, i = %lu\r\n", size, i); 
}

template<>
void CompareMem (media2_t<float> one, media2_t<float> two, size_t size, size_t i)
{ printf("CompareMem(): media2 size = %lu, i = %lu\r\n", size, i); 
}

