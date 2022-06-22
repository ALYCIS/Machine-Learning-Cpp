#include <iostream>
#include <bitset>

using namespace std;


class Answer
{
public:

  //
  static int getBit(unsigned int value, int pos)
  {
    std::bitset<32> value_binaire(value);
    return value_binaire[pos];
  }

};


int main(){


  cout<< Answer::getBit(5, 0) << endl; // 1
  cout<< Answer::getBit(5, 1) << endl; // 0
  cout<< Answer::getBit(5, 2) << endl; // 1

  return 0;
}
