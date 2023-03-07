#include <iostream>
#include <functional>
#include <string>
#ifndef BLOOM_FILTER_H
#define BLOOM_FILTER_H

//Hash functions obtained from: https://www.partow.net/programming/hashfunctions/#AvailableHashFunctions

unsigned int APHash(const char* str, unsigned int length);

unsigned int DEKHash(const char* str, unsigned int length);

unsigned int DJBHash(const char* str, unsigned int length);

unsigned int SDBMHash(const char* str, unsigned int length);

unsigned int BKDRHash(const char* str, unsigned int length);

unsigned int ELFHash(const char* str, unsigned int length);

unsigned int PJWHash(const char* str, unsigned int length);

unsigned int JSHash(const char* str, unsigned int length);

unsigned int RSHash(const char* str, unsigned int length);

unsigned int stl_hash(std::string str);

class bloom_filter {
  std::string prob_array;
  int m;
  int k;
  int n;
  float p = 0.005;
  int num_of_elements = 0;

  public:
  bloom_filter(int m, int k, int n);

  void allocate_data();

  void insert(char str[]);

  int search(char str[]);
};

#endif