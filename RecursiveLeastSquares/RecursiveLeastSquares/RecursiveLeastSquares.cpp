// RecursiveLeastSquares.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include "rls.h"
#include <array>

double linf(double x) {
    return 0.5 * x + 1.5;
}

int main()
{
    auto r = new Rls<2>(100,0.999);

    std::array<double,2> tmp{ 0,1 };

    for (int i = -10; i < 10; i++) {
        tmp[0] = i;
        std::cout << tmp[0] << " " << tmp[1] <<  " " << linf(i) << "\n";
        r->update(tmp, linf(i));
        r->print();
    }
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
