#ifndef _FUN_H
#define _FUN_H

#include <vector>
#include <algorithm>
#include <cstdlib>

// search函数定义：用于在一个vector中查找key的索引
int search(int key, const std::vector<int>& keys) {
  auto it = std::find(keys.begin(), keys.end(), key);
  if (it != keys.end()) {
    return std::distance(keys.begin(), it);
  }
  return -1; // 如果未找到，返回-1
}

// 生成从n1到n2的数组(返回vector)
std::vector<int> N_array(int n1, int n2) {
  std::vector<int> a(n2 - n1);
  std::iota(a.begin(), a.end(), n1);
  return a;
}

// 抽取2个不相同的随机数(返回vector)
std::vector<int> rand_array(int nb) {
  std::vector<int> randbox(2);
  int index = 0, random = 0;
  while (index < 2) {
    random = rand() % nb;
    if (std::find(randbox.begin(), randbox.begin() + index, random) == randbox.begin() + index) {
      randbox[index] = random;
      index++;
    }
  }
  return randbox;
}

#endif /* _FUN_H */
