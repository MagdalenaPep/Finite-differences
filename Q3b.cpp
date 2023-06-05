// Include all necessary libraries
#include <iostream>
#include <iomanip>
#include <valarray>
#include <cmath>
#include <vector>
using namespace  std;

// Struct to calculate the weighted norm of a valarray copied from question 2C
struct WeightedNorm {
    int m; // This is the member variable, it will be used as the power in the calculation.
    WeightedNorm(int m) : m(m) {}
    // Function call operator to calculate the weighted norm of a valarray<double> using the specified power m.
    long double operator()(const valarray<long double> u) const {
        long double sum = 0.0; // This variable will store the sum for the weighted norm calculation.

        // For loop iterates through all elements of u, raises the absolute value of them to the power of m, and adds the result to sum.
        for (int i = 0; i < u.size(); ++i) {
            sum += pow(abs(u[i]), m);
        }
        return pow(sum, 1.0 / m); // Returns the m-th root of sum.
    }
};


int main() {
    // Define the given variables: N, a, b, and dx
    long double a = -1;       // left endpoint of the interval
    long double b = 1;        // right endpoint of the interval

    //Create a vector N to store the different grid-point numbers
    //Using pushback add elements to the back of the vector N which contains the n values through which a for lop will iterate
    vector<int> N;
    N.push_back(15);
    N.push_back(31);
    N.push_back(63);
    N.push_back(127);

    // Print the header for the output
    cout << " N+1 \t  N^2*Mean Error" << endl;

    // Iterate through every value in N to find the mean error
    for (size_t i = 0; i < N.size(); i++) {
        int N_val = N[i];
        long double d_x = (b - a) / N_val; // Calculate grid spacing
        //Define the vallerays
        valarray<long double> Derivative_val(N_val + 1);
        valarray<long double> Xi(N_val + 1);
        valarray<long double> F_xi(N_val + 1);

        //Code to fill in vallerays is the same  as from part 3a
        // Fill the xi values and f(xi) values for all grid points
        for (int i = 0; i <= N_val; i++) {
            Xi[i] = a + i * d_x; // Calculate xi
            F_xi[i] = sin(3 * Xi[i]); // Calculate f(xi)
        }
        // Find the values for the derivative_val valarray using the provided approximations for all values of f'(x)
        //The first instance of the derivative of the function is f'(0) (derivative_val(0)) not including the error
        Derivative_val[0] = (-3 * F_xi[0] + 4 * F_xi[1] - F_xi[2]) / (2 * d_x);

        //Using a for loop iterate through all values from i=1 to i=n-1 to find the derivative f'(i) not including error
        for (int i = 1; i < N_val; i++) {
            Derivative_val[i] = (F_xi[i + 1] - F_xi[i - 1]) / (2 * d_x);
        }
        //Finally find the derivative  f(N) not including error
        Derivative_val[N_val] = (F_xi[N_val - 2] - 4 * F_xi[N_val - 1] + 3 * F_xi[N_val]) / (2 * d_x);

        // Calculate the absolute value of the error as e = f'(x)(numerical) - f'(x)(analytical)
        //Analytical derivative is 3cos(3*xi)
        valarray<long double> error(N_val + 1);// Create a valarray error to store the error values for all grid points. The size of the valarray is N + 1.
        for (int i = 0; i <= N_val; i++) {
            error[i] = abs(Derivative_val[i] - 3*cos(3*Xi[i]));
        }


        // Calculate the mean error using the L1 norm
        WeightedNorm l1_norm(1);
        long double mean_error = l1_norm(error) / (N_val + 1);

        // Output N²⟨e⟩ to 20 points of accuracy
        long double result = (N_val*N_val) * mean_error;
        cout << N_val +1 << "   \t "     << setprecision(20) << result << endl;
    }

    return 0;
}