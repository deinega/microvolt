// ������� �� ������ � �������
#include "memutils.h"
#include <vector>

/// ����: ���������_�����, �����������_�����.
struct aligned_ptr_t
{ void* aligned;
  void* original;
  aligned_ptr_t (void* _original, void* _aligned) : original(_original), aligned(_aligned) {};
};

/// ����� �������� ������ � ������: ���������_�����, �����������_�����.
std::vector<aligned_ptr_t> aligned_ptr_list;

/// �������, ������� ���������� �� ���������� ������ aligned_ptr_list
/// ����: ���������_�����, �����������_�����.
/// (��� ������� �� ��������� � ���� ������)
void aligned_ptr_list_push (void* original, void* aligned)
{ aligned_ptr_list.push_back ( aligned_ptr_t (original, aligned) );
}

/// �������, ������� ���� ���� : ���������_�����, �����������_�����
/// �� ���������� ������� aligned_ptr_list.
/// (��� ������� �� ��������� � ���� ������)
void* aligned_ptr_list_pop (void* Value)
{ void* ret = 0;
  // ��������� ���� ������ � ������: ���������_�����, �����������_�����
  for (std::vector<aligned_ptr_t>::iterator it = aligned_ptr_list.begin();
       it != aligned_ptr_list.end(); it++)
    // ���� ������� ����, ��� ����������� ��������� ����� ����������
    if ( (*it).aligned == Value )
    { // ������� ��������� �����
      ret = (*it).original;
      // ������� �� ������� ������
      aligned_ptr_list.erase(it);
      // ����� �� �����
      break;
    }
  // ���� ��������������� ������ �� ���� �������, ����� ��������� 0
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

