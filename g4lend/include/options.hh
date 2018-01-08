/*
 
 parse command line arguments
 
 options("blah",argc,argv)  return true if blah is in any argument that starts with '-'
 options(argc,argv)         return vector of arguments that do not start with '-'
 
 Doug Wright 2/28/15
 */


#include <vector>
#include <string>
using namespace std;

//....return true if arg is in any argument that begins with '-'
//    for example: -xyz
int options(const char* arg, int argc, char* argv[]){
    
    //....loop over all arguments, skip argv[0], since it is program name
    for( int i=1; i<argc; i++)
        //....return true if starts with '-' and contains arg
        if (argv[i][0]=='-' && (strstr(argv[i],arg) != NULL) ) return true;
    return false;
}

//....return all arguments that do NOT begin with '-'
vector<string> options(int argc, char* argv[]){
    vector<string> args;
    
    //....loop over all arguments, skip argv[0], since it is program name
    for( int i=1; i<argc; i++){
        if ( argv[i][0]=='-' ) continue;
        args.push_back(argv[i]);
    }
    return args;
}
