/**
 * @author Roman Pahl
 */

#ifndef permory_helper_algorithm_hpp
#define permory_helper_algorithm_hpp

#include <vector>

#include "config.hpp"
#include "helper/tokenizer.hpp"

namespace Permory 
{
    // Copy tokens into vector implemented for various types
    template<class T> inline void copy_token(Tokenizer&, std::vector<T>&);
    template<> inline void copy_token<int>(Tokenizer& tok, std::vector<int>& v)
    {
        while(!tok.empty()){
            v.push_back(atoi(*tok));
            ++tok;
        }
    }
    template<> inline void copy_token<unsigned short int>(
            Tokenizer& tok, std::vector<unsigned short int>& v)
    {
        while(!tok.empty()){
            v.push_back(atoi(*tok));
            ++tok;
        }
    }
    template<> inline void copy_token<double>(
            Tokenizer& tok, std::vector<double>& v)
    {
        while(!tok.empty()){
            v.push_back(atof(*tok));
            ++tok;
        }
    }
    template<> inline void copy_token<char>(
            Tokenizer& tok, std::vector<char>& v)
    {
        while(!tok.empty()){
            v.push_back((*tok)[0]);
            ++tok;
        }
    }
    template<> inline void copy_token<std::string>(
            Tokenizer& tok, std::vector<std::string>& v)
    {
        while(!tok.empty()){
            v.push_back(std::string(*tok));
            ++tok;
        }
    }
} // namespace Permory

#endif // include guard

