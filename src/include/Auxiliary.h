/** @file 
 *  @author Georgios G. Vogiatzis (<gvog@chemeng.ntua.gr>)
 *  @version 1.0 (January 7, 2014)
 *  @brief  A C++ header file containing auxiliary type definitions and functions.
 * 
 *  @section LICENSE
 *      
 *  Copyright (c) 2014 Georgios Vogiatzis. All rights reserved.
 *  
 *  This work is licensed under the terms of the MIT license.  
 *  For a copy, see <https://opensource.org/licenses/MIT>.
 * 
 */

#ifndef AUXILIARY_H
#define	AUXILIARY_H

#include<string>
#include<sys/stat.h>
#include<sstream>
#include<string>
#include<vector>


using namespace std;

std::vector<string> tokenize(std::string input_string); ///< A function for tokenizing a string into
                                                        ///< substrings.
double get_wall_time();
double get_cpu_time();

template <typename T>
T StringToNumber ( const string &Text )//Text not by const reference so that the function can be used with a 
{                               //character array as argument
	stringstream ss(Text);
	T result;
	return ss >> result ? result : 0;
}

template <typename T>
string NumberToString ( T Number )
{
	stringstream ss;
	ss << Number;
	return ss.str();
}


void create_dir(string name_of_dir);



/** @param[in] x The floating-point number to be raised to an integer power.
 *  @param[in] n The power to which the point will be raised.
 */
static inline double powint(const double &x, const int n) {
   double yy, ww;

   if (x == 0.0) return 0.0;
   int nn = (n > 0) ? n : -n;
   ww = x;

   for (yy = 1.0; nn != 0; nn >>= 1, ww *= ww)
      if (nn & 1) yy *= ww;

   return (n > 0) ? yy : 1.0 / yy;
} ///< Optimized version of pow(x,n) with n being integer. It can be up to 10x faster than pow(x,y).



/** @class D3DVector
 *  @brief A class describing a vector with three components of type T.
 * 
 *  Adapted from: http://rosettacode.org/
 */
template< class T >
class D3DVector {

public :
   T x;     ///< The x-component of the D3DVector vector.
   T y;     ///< The y-component of the D3DVector vector.
   T z ;    ///< The z-component of the D3DVector vector.
   
   D3DVector( T a , T b , T c ) {
      x = a ;
      y = b ;
      z = c ;
   }///< The constructor of the D3DVector class, assigning value to each one of the three components.

   D3DVector(){
      x = 0.0;
      y = 0.0;
      z = 0.0;
   } ///< Constructor which initializes all of the three components of the vector to zero.

   void zero (){
      x = (T)0.0;
      y = (T)0.0;
      z = (T)0.0;
   }  ///< Set all components of the vector to zero.
   
   D3DVector operator=(const D3DVector & rhs){
      return (D3DVector(rhs.x, rhs.y, rhs.z));
   } ///< The assignment operator defined for the D3DVector vector.

   D3DVector operator-(const D3DVector & rhs){
      return (D3DVector(x-rhs.x, y-rhs.y, z-rhs.z));
   }///< The substraction operator defined for two D3DVector vectors.

   D3DVector operator+(const D3DVector & rhs){
      return (D3DVector(x+rhs.x, y+rhs.y, z+rhs.z));
   }///< The addition operator defined for two D3DVector vectors.

   D3DVector operator!(void){
      T inorm = 1.0 / sqrt(x*x + y*y + z*z);
      T nx = inorm * x;
      T ny = inorm * y;
      T nz = inorm * z;
      return (D3DVector(nx,ny,nz));
   }///< The operator ``!'' is used for calculating the norm of a D3DVector vector.

   T norm( ){
      return (sqrt(x*x + y*y + z*z));
   }///< A function calculating the norm of a D3DVector vector.

   T dotproduct (const D3DVector & rhs ) {
      T scalar = x * rhs.x + y * rhs.y + z * rhs.z ;
      return scalar ;
   }///< The dot (inner) product defined for two D3DVector vectors.


   /* gvog: DOT PRODUCT operator is defined: */
   T operator*(const D3DVector & rhs){
      T scalar = x * rhs.x + y * rhs.y + z * rhs.z ;
      return scalar ;
   }///< The operator ``*'' is used for the dot product between two D3DVector vectors.

   /* gvog: SCALAR PRODUCT operator is defined: */
   D3DVector operator*(T mult){
      return (D3DVector(x*mult, y*mult, z*mult));
   }///< The operator ``*'' defined for the scalar product between a single number and a D3DVector
    ///< vector.
   
   /* gvog: SCALAR DIVISION operator is defined: */
   D3DVector operator/(T mult){
         return (D3DVector(x/mult, y/mult, z/mult));
   }///< The operator ``/'' is used for the scalar division between a D3DVector vector and a number.

   D3DVector crossproduct ( const D3DVector & rhs ) {
      T a = y * rhs.z - z * rhs.y ;
      T b = z * rhs.x - x * rhs.z ;
      T c = x * rhs.y - y * rhs.x ;
      D3DVector product( a , b , c ) ;
      return product ;
   }///< The function caclculating the cross product between two vectors.

   D3DVector operator^ (const D3DVector & rhs){
      return crossproduct(rhs);
   }///< The operator ``^'' is used for the cross product between D3DVector vectors.

   D3DVector triplevec( D3DVector & a , D3DVector & b ) {
      return crossproduct ( a.crossproduct( b ) ) ;
   }///< The triple vector product between two D3DVector vectors.

   T triplescal( D3DVector & a, D3DVector & b ) {
      return dotproduct( a.crossproduct( b ) ) ;
   }///< The mixed product.

} ;

#endif	/* AUXILIARY_H */

