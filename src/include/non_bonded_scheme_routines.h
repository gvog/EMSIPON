/** @file
 *  @author Grigorios Megariotis (<gmegariotis@yahoo.gr>)
 *  @brief The header file containing the definitions of the functions of the smoothed non-bonded free energy estimation
 *         scheme.
 */

#ifndef NEW_SCHEME_ROUTINES_H
#define	NEW_SCHEME_ROUTINES_H

void phi_function(double x, double xj, double Rj, double delta, double &phi_value);
void integral_minus_infinity_to_x(double xasterisk, double xj, double Rj, double delta, double &integral_value);
void integral_x_to_plus_infinity(double xasterisk, double xj, double Rj, double delta, double &integral_value);

#endif	/* NEW_SCHEME_ROUTINES_H */

