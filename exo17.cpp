public:
  static bool locateUniverseFormula(string & formulaPath)
  {

#include <string>
#include <iostream>
#include <cstring>
#include <filesystem>
using namespace std;

#define FORMULA_PATH "/tmp/documents/11985/universe-formula"

using namespace std;

class Answer
{
    string file_name = "universe-formula";
    string first_path {"/tmp/documents"};

    string absolute_path = first_path + "/" + file_name;

    

    std::filesystem::path fichier(absolute_path);

    //cout<<"Repertoire premier : "<< fichier.string() << std::endl;

    if(std::filesystem::exists(fichier))
      {
        formulaPath = absolute_path;
        //cout<<" Le chemin est : "<<formulaPath<<endl;
        return true;
      }

    for(auto const & f : std::filesystem::recursive_directory_iterator(first_path,std::filesystem::directory_options::follow_directory_symlink))
      {
        const std::string s = f.path();

        //cout<<"Repertoire : "<< s << std::endl;

        if( f.is_directory() )
          {
            absolute_path = f.path().string()+"/"+file_name;
            auto fichier_test = std::filesystem::path(absolute_path);

            bool  test = std::filesystem::exists(fichier_test);

            if(test)
              {
                // //cout<<file_name<<endl;
                // std::filesystem::directory_entry dossier_fichier(first_path);
                formulaPath = absolute_path;
            
                // cout<<"Le fichier existe ? = "<< test << endl;
                // cout<<" Le chemin est : "<<formulaPath<<endl;
                return true;
              }
          }
      }

    return false;
  }


};


int main()
{
  string pathResult;
  bool found = Answer::locateUniverseFormula(pathResult);

  cout<< "In this example, the formula path is :" <<endl;
  cout<< FORMULA_PATH << endl;
  cout << "and you found it at :"<<endl;
  cout << pathResult << endl;

  return 0;
}


