/** @file
 *  @author Georgios G. Vogiatzis (<gvog@chemeng.ntua.gr>)
 *  @version 1.0 (January 24, 2014)
 *  @brief The C++ source code file containing some auxiliary functions. 
 * 
 *  @section LICENSE
 *      
 *  Copyright (c) 2014 Georgios Vogiatzis. All rights reserved.
 *  
 *  This work is licensed under the terms of the MIT license.  
 *  For a copy, see <https://opensource.org/licenses/MIT>.
 */

#include<string>
#include<sstream>
#include<sys/time.h>
#include<vector>

#include<Auxiliary.h>

using namespace std;

std::vector<string> tokenize(std::string input_string) {

   string buf;
   stringstream ss(input_string);
   std::vector<string> tokens;

   while (ss >> buf)
      tokens.push_back(buf);

   return tokens;
}

void create_dir(string name_of_dir){
   struct stat st;
   if (stat(name_of_dir.c_str(), &st) != 0)
      mkdir(name_of_dir.c_str(), 0750);
}


double get_wall_time(){
    struct timeval time;
    if (gettimeofday(&time,NULL)){
        //  Handle error
        return 0;
    }
    return (double)time.tv_sec + (double)time.tv_usec * .000001;
}

double get_cpu_time(){
    return (double)clock() / CLOCKS_PER_SEC;
}


