#include "./bloom_filter.h"

//Hash functions obtained from: https://www.partow.net/programming/hashfunctions/#AvailableHashFunctions

unsigned int APHash(const char* str, unsigned int length)
{
   unsigned int hash = 0xAAAAAAAA;
   unsigned int i    = 0;

   for (i = 0; i < length; ++str, ++i)
   {
      hash ^= ((i & 1) == 0) ? (  (hash <<  7) ^ (*str) * (hash >> 3)) :
                               (~((hash << 11) + ((*str) ^ (hash >> 5))));
   }

   return hash;
}

unsigned int DEKHash(const char* str, unsigned int length)
{
   unsigned int hash = length;
   unsigned int i    = 0;

   for (i = 0; i < length; ++str, ++i)
   {
      hash = ((hash << 5) ^ (hash >> 27)) ^ (*str);
   }

   return hash;
}

unsigned int DJBHash(const char* str, unsigned int length)
{
   unsigned int hash = 5381;
   unsigned int i    = 0;

   for (i = 0; i < length; ++str, ++i)
   {
      hash = ((hash << 5) + hash) + (*str);
   }

   return hash;
}

unsigned int SDBMHash(const char* str, unsigned int length)
{
   unsigned int hash = 0;
   unsigned int i    = 0;

   for (i = 0; i < length; ++str, ++i)
   {
      hash = (*str) + (hash << 6) + (hash << 16) - hash;
   }

   return hash;
}

unsigned int BKDRHash(const char* str, unsigned int length)
{
   unsigned int seed = 131; /* 31 131 1313 13131 131313 etc.. */
   unsigned int hash = 0;
   unsigned int i    = 0;

   for (i = 0; i < length; ++str, ++i)
   {
      hash = (hash * seed) + (*str);
   }

   return hash;
}

unsigned int ELFHash(const char* str, unsigned int length)
{
   unsigned int hash = 0;
   unsigned int x    = 0;
   unsigned int i    = 0;

   for (i = 0; i < length; ++str, ++i)
   {
      hash = (hash << 4) + (*str);

      if ((x = hash & 0xF0000000L) != 0)
      {
         hash ^= (x >> 24);
      }

      hash &= ~x;
   }

   return hash;
}

unsigned int PJWHash(const char* str, unsigned int length)
{
   const unsigned int BitsInUnsignedInt = (unsigned int)(sizeof(unsigned int) * 8);
   const unsigned int ThreeQuarters     = (unsigned int)((BitsInUnsignedInt  * 3) / 4);
   const unsigned int OneEighth         = (unsigned int)(BitsInUnsignedInt / 8);
   const unsigned int HighBits          =
                      (unsigned int)(0xFFFFFFFF) << (BitsInUnsignedInt - OneEighth);
   unsigned int hash = 0;
   unsigned int test = 0;
   unsigned int i    = 0;

   for (i = 0; i < length; ++str, ++i)
   {
      hash = (hash << OneEighth) + (*str);

      if ((test = hash & HighBits) != 0)
      {
         hash = (( hash ^ (test >> ThreeQuarters)) & (~HighBits));
      }
   }

   return hash;
}

unsigned int JSHash(const char* str, unsigned int length)
{
   unsigned int hash = 1315423911;
   unsigned int i    = 0;

   for (i = 0; i < length; ++str, ++i)
   {
      hash ^= ((hash << 5) + (*str) + (hash >> 2));
   }

   return hash;
}

unsigned int RSHash(const char* str, unsigned int length)
{
   unsigned int b    = 378551;
   unsigned int a    = 63689;
   unsigned int hash = 0;
   unsigned int i    = 0;

   for (i = 0; i < length; ++str, ++i)
   {
      hash = hash * a + (*str);
      a    = a * b;
   }

   return hash;
}

unsigned int stl_hash(std::string str) {
  std::hash<std::string> str_hash;

  return str_hash(str);
}


bloom_filter::bloom_filter(int m, int k, int n) {
  this->m = m;
  this->k = k;
  this->n = n;
}

void bloom_filter::allocate_data() {
  for(int i = 0; i < this->m; i++)
    this->prob_array.append("0");
}

void bloom_filter::insert(char str[]) {
  int hash[this->k];

  ///call hash functions here. if k is more than 8, this needs to be increased.
  hash[0] = APHash(str, strlen(str)) % this->m;
  hash[1] = DEKHash(str, strlen(str)) % this->m;
  hash[2] = DJBHash(str, strlen(str)) % this->m;
  hash[3] = SDBMHash(str, strlen(str)) % this->m;
  hash[4] = BKDRHash(str, strlen(str)) % this->m;
  hash[5] = JSHash(str, strlen(str)) % this->m;
  hash[6] = RSHash(str, strlen(str)) % this->m;
  hash[7] = stl_hash(str) % this->m;
  //hash[8] = PJWHash(str, strlen(str)) % this->m;
  //hash[9] = ELFHash(str, strlen(str)) % this->m;

  //std::cout<<hash[0]<<std::endl;

  for(int i = 0; i < this->k; i++)
    this->prob_array[hash[i]] = '1';

  return;
}

int bloom_filter::search(char str[]) {
  int hash[this->k];

  ///call hash functions here. if k is more than 8, this needs to be increased.
  hash[0] = APHash(str, strlen(str)) % this->m;
  hash[1] = DEKHash(str, strlen(str)) % this->m;
  hash[2] = DJBHash(str, strlen(str)) % this->m;
  hash[3] = SDBMHash(str, strlen(str)) % this->m;
  hash[4] = BKDRHash(str, strlen(str)) % this->m;
  hash[5] = JSHash(str, strlen(str)) % this->m;
  hash[6] = RSHash(str, strlen(str)) % this->m;
  hash[7] = stl_hash(str) % this->m;
  //hash[8] = PJWHash(str, strlen(str)) % this->m;
  //hash[9] = ELFHash(str, strlen(str)) % this->m;

  int not_found = 0;
  for(int i = 0; i < this->k; i++) {
    if(this->prob_array[hash[i]] == '0')
      not_found++;
  }

  if(not_found > 0)
    return 0;
  else
    return 1;
}
