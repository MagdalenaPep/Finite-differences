// Include all necessary libraries
#include <iostream>
#include <iomanip>
#include <valarray>
#include <cmath>
using namespace  std;

int main() {
    // Define the given variables: N, a, b, and dx
    int N = 31;               // number of grid points
    long double a = -1;       // left endpoint of the interval
    long double b = 1;        // right endpoint of the interval
    long double dx = (b - a) / N;   // grid spacing

    // Create empty valarrays for the xi values, f(xi) values, and derivative values
    valarray<long double> derivative_val(N + 1);
    valarray<long double> xi(N + 1);
    valarray<long double> f_xi(N + 1);

    // Fill the xi values and f(xi) values for all grid points
    for (int i = 0; i <= N; i++) {
        xi[i] = a + i * dx;    // xi = a + i*dx
        f_xi[i] = sin(3*xi[i]); // f(xi) = sin(3*xi)
    }

    // Find the values for the derivative_val valarray using the provided approximations for all values of f'(x)
    //The first instance of the derivative of the function is f'(0) (derivative_val(0)) not including the error
    derivative_val[0] = (-3 * f_xi[0] + 4 * f_xi[1] - f_xi[2]) / (2 * dx);

    //Using a for loop iterate through all values from i=1 to i=n-1 to find the derivative f'(i) not including error
    for (int i = 1; i < N; i++) {
        derivative_val[i] = (f_xi[i + 1] - f_xi[i - 1]) / (2 * dx);
    }
    //Finally find the derivative f(N) not including error
    derivative_val[N] = (f_xi[N - 2] - 4 * f_xi[N - 1] + 3 * f_xi[N]) / (2 * dx);

    // Calculate the absolute value of the error as e = f'(x)(numerical) - f'(x)(analytical)
    //Analytical derivative is 3cos(3*xi)
    valarray<long double> errors(N + 1);// Create a valarray error to store the error values for all grid points. The size of the valarray is N + 1.
    for (int i = 0; i <= N; i++) {
        errors[i] = abs(derivative_val[i] - 3*cos(3*xi[i]));
    }

    //Print the header of the output table with the column names.
    cout << "Error values (ei)" << endl;
    //Loop through all the grid points(from 0 to N) and print only the error values for each point
    for (int i = 0; i <= N; i++) {
        cout << setw(15) << fixed << setprecision(20) << errors[i] << " ";
    }
    cout << endl;

}