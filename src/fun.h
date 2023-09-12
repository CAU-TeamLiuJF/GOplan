#ifndef _FUN_H
#define _FUN_H

#include <errno.h>

extern int errno;

// search函数定义
int search(int key, int keys[], int length)
{
  int ret = -1;
  for (int i = 0; i < length; i++) {
    if (keys[i] == key) {
      ret = i;
      break;
    }
  }
  return ret;
}

// 生成n1-n2的数组(返回一维数组)
int* N_array(int* a, int n1, int n2)
{
  for (int i = n1; i < n2; i++) {
    a[i - n1] = i;
  }
  return a;
}
// 抽随机2个不相同的数(返回一维数组)
int* rand_array(int nb)
{
  static int randbox[2];
  int index = 0, check = 0, random = 0;
  while (index < 2) {
    random = rand() % nb;
    for (int j = 0; j < index; j++) {
      if (random != randbox[j]) {
        check = 0;
      }
      else {
        check = 1;
        break;
      }
    }
    if (check == 1) {
      continue;
    } else {
      randbox[index] = random;
      index++;
    }
  }
  return randbox;
}
#endif /* _FUN_H */
