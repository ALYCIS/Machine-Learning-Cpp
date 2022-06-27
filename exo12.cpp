#include <string>
#include <iostream>

using namespace std;



class Answer
{
public:
  static bool isTwin(string a, string b){
    int n = a.size();
    int m = b.size();

    cout<<" a = "<<a<<" b= "<<b<<endl;

    if(n != m) return false;


    // On converti le tout en miniscule ou en majuscule
    int somme_valeur1 = 0;
    int somme_valeur2 = 0;

    for(int i = 0; i < n; i++)
      {

        char caractere1 = a[i];
        auto cherche1 = b.find(caractere1);


        // Pour le premier tableau
        if(cherche1 == string::npos) // On a pas retrouvé la lettre
          {
            // on converti le caractere
            if(65<= caractere1 && caractere1 <= 90) caractere1 = caractere1 + 32; // miniscule
            else if(97 <= caractere1 && caractere1 <= 122) caractere1 = caractere1 - 32; // Conversion majuscule

            auto cherche_test = b.find(caractere1);
            if(cherche_test == string::npos) return false;
          }

        // Pour le deuxieme tableau b
        char caractere2 = b[i];
        auto cherche2 = a.find(caractere2);

        if(cherche2 == string::npos) // On a pas retrouvé la lettre  
          {
            // on converti le caractere
            if(65<= caractere2 && caractere2 <= 90) caractere2 = caractere2 + 32; // miniscule
            else if(97 <= caractere2 && caractere2 <= 122) caractere2= caractere2 - 32; // Conversion majuscule

            auto cherche_test = a.find(caractere2);
            if(cherche_test == string::npos) return false;
          }

        somme_valeur1 += (65<= caractere1 && caractere1 <= 90)? (  int(caractere1 + 32) ): ( int (caractere1));

        somme_valeur2 += (65<= caractere2 && caractere2 <= 90)? (  int(caractere2 + 32) ): ( int (caractere2));
      }

    cout<<"somme1 =" << somme_valeur1<<" somme2 = "<<somme_valeur2<<endl;
    if(somme_valeur1 != somme_valeur2) return false;

    return true; 
  }
 
};
//kakouj     joujak


int main(){

  cout<< boolalpha << Answer::isTwin("Hello", "world") << endl; // false
  cout<< boolalpha << Answer::isTwin("aaba", "bcda") << endl; // true
  cout<< boolalpha << Answer::isTwin("Lookout", "Outlook") << endl; // true
  cout<< boolalpha << Answer::isTwin("kakouj", "AkoJJu") << endl; // false

  return 0;
}
