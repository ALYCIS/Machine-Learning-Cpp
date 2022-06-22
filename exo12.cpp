#include <string>
#include <iostream>

using namespace std;

class Answer
{
public:
  static bool isTwin(string a, string b){
    int n = a.size();


    // On converti le tout en miniscule ou en majuscule

    for(int i = 0; i < n; i++)
      {

        char caractere = a[i];
        auto cherche = b.find(caractere);

        if(cherche == string::npos) // On a pas retrouvé la lettre
          {
            // on converti le caractere
            if(65<= caractere && caractere <= 90) caractere = caractere + 32; // miniscule
            else caractere = caractere - 32; // Conversion majuscule

            auto cherche1 = b.find(caractere);
            if(cherche1 == string::npos) return false;
          }

      }

    return true;
  }

};



int main(){

  cout<< boolalpha << Answer::isTwin("Hello", "world") << endl; // false
  cout<< boolalpha << Answer::isTwin("abc", "bca") << endl; // true
  cout<< boolalpha << Answer::isTwin("Lookout", "Outlook") << endl; // false

  return 0;
}
