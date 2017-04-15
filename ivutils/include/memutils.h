#ifndef MEMUTILS_H
#define MEMUTILS_H

/// \ru @file memutils.h
/// \brief    ������� ������� �� ������ � ������� � ��������� ������� ����.
///           ���� ���������� USE_CUDA, ����� #include "cuda_memutils.h"
///           ���� ����������� ������� ������� ��� ������ � CUDA �������.

# include <stdlib.h>


// ��������� ����

// ��������� ��� �������� ������ ����� �������������
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
  { if (c1!=r.c1) // ���� �� �����
    { if (less_than(c1,r.c1)) return true; // ���� c1<r.c1
      else                 return false;
    }
    if (c2!=r.c2) // ���� �� �����
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

// ������ � �������

/// �������, ������� ���������� �� ���������� ������ aligned_ptr_list
/// ����: ���������_�����, �����������_�����.
/// (��� ������� �� ��������� � ���� ������ (�� ������)
///  � ����������� � ����� memutils.cpp)
void aligned_ptr_list_push (void* original, void* aligned);
/// �������, ������� ���� ���� : ���������_�����, �����������_�����
/// �� ���������� ������� aligned_ptr_list.
/// (��� ������� �� ��������� � ���� ������ (�� ������)
///  � ����������� � ����� memutils.cpp)
void* aligned_ptr_list_pop (void* aligned);

/// �������� ������, ������� ����������� ��������� �� ��� ������
/// �������� ������������ ��������� ������ sizeOf(T)*mult
template <class T>
bool aligned_alloc (T **ptr, size_t size, size_t mult = 1)
{ // ���������
  T* original; T* aligned;
  // ��������� ��� ������ size ������ mult �, ���� �� ������, ���������.
  size_t strcount = size / mult;
  if (size % mult) strcount++;
  // �������� ������ (� ������� +1 ��� ������������)
  original = (T*) new char [sizeof(T) * mult * (strcount+1)];
  // ��������� ������� �� ������� �� ����������� ������
  size_t strsize = mult*sizeof(T);
  size_t rest = ((size_t)original) % strsize;
  // ��������� �����
  aligned = (T*) ( (size_t)original + ( (rest) ? (strsize-rest):0 ));
  // ��������� ����: ���������_�����, �����������_�����
  aligned_ptr_list_push ( (void*)original, (void*)aligned );
  // ���������������� �������
  for (size_t i=0; i< size; i++)
    aligned[i] = T();
  // ������� ����������� �����
  *ptr = aligned;
  return true;
}

/// ���������� ������, �� ������� ��������� ����������� ���������.
/// ���� ���������� ��������� ������������� ������ � aligned_ptr_list � ������������ ����������, �����
/// ������ ������������� � ������� ���������� true, ����� ������� ���������� false
template <class T>
bool aligned_free (T *ptr)
{ // �������� ��������� �����
  void* original = aligned_ptr_list_pop ((void*)ptr);
  // ���� ����� �������, ����� ���������� ������, � ������� ��
  if (original !=0)
  { delete[] (T*)original;
    return true;
  }
  // ����� ������� ������
  else
    return false;
}

// ���� ��������� USE_CUDA �������� ���� ������� �������
// ��� ������ � CUDA �������
#ifdef USE_CUDA
  #include "cuda_memutils.h"
#endif

#endif//MEMUTILS_H

