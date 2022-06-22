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

    // cout<<"Bitset de :"<<value<<endl;
    // for(int i = 0; i <32 ; i++)
    //   {
    //     cout<<value_binaire[31-i];
    //   }
    // cout<<endl;

    return value_binaire[pos];
  }

};






int main(){


  cout<< Answer::getBit(5, 0) << endl; // 1
  cout<< Answer::getBit(5, 1) << endl; // 0
  cout<< Answer::getBit(5, 2) << endl; // 1

  return 0;
}
