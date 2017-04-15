#ifndef MEMUTILS_H
#define MEMUTILS_H

/// \ru @file memutils.h
/// \brief    Базовые функции по работе с памятью и некоторые базовые типы.
///           Если определено USE_CUDA, тогда #include "cuda_memutils.h"
///           сюда добавляются базовые функции для работы с CUDA памятью.

# include <stdlib.h>


// НЕКОТОРЫЕ ТИПЫ

// Структуры для хранения данных медиа коэффициентов
template <class Type>
struct media_t
{ Type c1, c2;
  media_t ()                   : c1(0),   c2(0)   {}
  media_t (Type _c1, Type _c2) : c1(_c1), c2(_c2) {}
  bool operator!= (const media_t& src) const
  { return ( (c1!=src.c1) || (c2!=src.c2) );
  }
  bool operator== (const media_t& src) const
  { return ( (c1==src.c1) && (c2==src.c2) );
  }
  /*bool operator< (const media_t& r) const
  { if (c1!=r.c1) // если не равно
    { if (less_than(c1,r.c1)) return true; // если c1<r.c1
      else                 return false;
    }
    if (c2!=r.c2) // если не равно
    { if (less_than(c2,r.c2)) return true;
      else                 return false;
    }
    return false;
  } */
};

template <class Type>
struct media2_t
{ Type c1, c2[2];
  media2_t (Type _c1=0, Type _c21=0, Type _c22=0) : c1(_c1)
  { c2[0]=_c21; c2[1]=_c22;
  }
  bool operator!= (const media2_t& src) const
  { return ( (c1!=src.c1) || (c2[0]!=src.c2[0]) || (c2[1]!=src.c2[1]) );
  }
  bool operator== (const media2_t& src) const
  { return ( (c1==src.c1) && (c2[0]==src.c2[0]) && (c2[1]==src.c2[1]) );
  }
  /*bool operator< (const media2& r) const
  { if (c1!=r.c1)
    { if (less_than(c1,r.c1)) return true;
      else                 return false;
    }
    for (int i=0; i<2; i++)
      if (c2[i] != r.c2[i])
      { if (less_than(c2[i],r.c2[i])) return true;
        else                       return false;
      }
    return false;
  }*/
};

template <class Type>
void CompareMem (Type one, Type two, size_t size, size_t i);

// РАБОТА С ПАМЯТЬЮ

/// Функция, которая записывает во внутренний массив aligned_ptr_list
/// пару: настоящий_адрес, выровненный_адрес.
/// (эта функция не привязана к типу данных (не шаблон)
///  и реализуется в файле memutils.cpp)
void aligned_ptr_list_push (void* original, void* aligned);
/// Функция, которая ищет пару : настоящий_адрес, выровненный_адрес
/// во внутреннем массиве aligned_ptr_list.
/// (эта функция не привязана к типу данных (не шаблон)
///  и реализуется в файле memutils.cpp)
void* aligned_ptr_list_pop (void* aligned);

/// выделить память, вернуть выровненный указатель на эту память
/// значение выровненного указателя кратно sizeOf(T)*mult
template <class T>
bool aligned_alloc (T **ptr, size_t size, size_t mult = 1)
{ // указатели
  T* original; T* aligned;
  // проверить что размер size кратен mult и, если не кратен, поправить.
  size_t strcount = size / mult;
  if (size % mult) strcount++;
  // выделить память (с запасом +1 под выравнивание)
  original = (T*) new char [sizeof(T) * mult * (strcount+1)];
  // вычислить сколько не хватает до выровненого адреса
  size_t strsize = mult*sizeof(T);
  size_t rest = ((size_t)original) % strsize;
  // выровнить адрес
  aligned = (T*) ( (size_t)original + ( (rest) ? (strsize-rest):0 ));
  // сохранить пару: настоящий_адрес, выровненный_адрес
  aligned_ptr_list_push ( (void*)original, (void*)aligned );
  // инициализировать объекты
  for (size_t i=0; i< size; i++)
    aligned[i] = T();
  // вернуть выровненный адрес
  *ptr = aligned;
  return true;
}

/// освободить память, на которую указывает выровненный указатель.
/// Если указанному указателю соответствует запись в aligned_ptr_list с оригинальным указателем, тогда
/// память освобождается и функция возварщает true, иначе функция возвращает false
template <class T>
bool aligned_free (T *ptr)
{ // получить настоящий адрес
  void* original = aligned_ptr_list_pop ((void*)ptr);
  // если адрес получен, тогда освободить память, и вернуть ОК
  if (original !=0)
  { delete[] (T*)original;
    return true;
  }
  // иначе вернуть ошибку
  else
    return false;
}

// если определен USE_CUDA включить сюда базовые функции
// для работы с CUDA памятью
#ifdef USE_CUDA
  #include "cuda_memutils.h"
#endif

#endif//MEMUTILS_H

