#include <pybind11/pybind11.h>
#include <matplot/matplot.h>
#include <complex>
#include <pybind11/complex.h>
#include <vector>
#include <pybind11/stl.h>
#include <cmath>
using namespace std;

#define STRINGIFY(x) #x
#define MACRO_STRINGIFY(x) STRINGIFY(x)
#define M_PI 3.14159265358979323846

int add(int i, int j) {
    return i + j;
}

namespace py = pybind11;

PYBIND11_MODULE(_core, m) {
    m.doc() = R"pbdoc(
        Pybind11 example plugin
        -----------------------

        .. currentmodule:: scikit_build_example

        .. autosummary::
           :toctree: _generate

           add
           subtract
    )pbdoc";

    m.def("add", &add, R"pbdoc(
        Add two numbers

        Some other explanation about the add function.
    )pbdoc");

    m.def("subtract", [](int i, int j) { return i - j; }, R"pbdoc(
        Subtract two numbers

        Some other explanation about the subtract function.
    )pbdoc");

    m.def("sinus", [](double czestotliwosc, double poczatek, double koniec, double probek, double amplituda) { 
        vector<double> x;
        vector<double> y;
        x=matplot::linspace(poczatek,koniec,probek);
        for (int i=0; i<x.size();i++){
            y.push_back(sin(czestotliwosc*x[i])*amplituda);
        }
        matplot::plot(x,y);
        matplot::grid(matplot::on);
        matplot::show();
        return y;
    }, py::arg("czestotliwosc"),py::arg("poczatek"),py::arg("koniec"),py::arg("probek"),py::arg("amplituda"));

    m.def("cosinus", [](double czestotliwosc, double poczatek, double koniec, double probek, double amplituda) { 
        vector<double> x;
        vector<double> y;
        x=matplot::linspace(poczatek,koniec,probek);
        for (int i=0; i<x.size();i++){
            y.push_back(cos(czestotliwosc*x[i])*amplituda);
        }
        matplot::plot(x,y);
        matplot::grid(matplot::on);
        matplot::show();
    
    }, py::arg("czestotliwosc"),py::arg("poczatek"),py::arg("koniec"),py::arg("probek"),py::arg("amplituda"));

    m.def("prostokatny", [](double czestotliwosc, double poczatek, double koniec, double probek, double amplituda) { 
        vector<double> x;
        vector<double> y;
        x=matplot::linspace(poczatek,koniec,probek);
        for (int i=0; i<x.size();i++){
                if (cos(czestotliwosc * x[i]) * amplituda >= 0) {
                    y.push_back(amplituda);
                }
                else {
                    y.push_back(-amplituda);
    }
        }
        matplot::plot(x,y);
        matplot::grid(matplot::on);
        matplot::show();
    
    }, py::arg("czestotliwosc"),py::arg("poczatek"),py::arg("koniec"),py::arg("probek"),py::arg("amplituda"));

    m.def("piloksztautna", [](double czestotliwosc, double poczatek, double koniec, double probek, double amplituda) { 
        vector<double> x;
        vector<double> y;
        x=matplot::linspace(poczatek,koniec,probek);
        for (int i=0; i<x.size();i++){
            y.push_back((fmod((amplituda *x[i] * czestotliwosc),2))-1);
        }
        matplot::plot(x,y);
        matplot::grid(matplot::on);
        matplot::show();
    
    }, py::arg("czestotliwosc"),py::arg("poczatek"),py::arg("koniec"),py::arg("probek"),py::arg("amplituda"));

    m.def("DFT",[](vector<double> input, double poczatek, double koniec, int dl){
        vector<double> modul;
        complex<double> i(0.0,1.0);
        complex<double> zesp;
        vector<complex<double>> y;
        vector<double> x = matplot::linspace(poczatek,koniec,dl);
        double odwri;

        for(int t=0; t<dl; t++){
            zesp=0;
            for(int w=0; w<dl; w++){
                odwri = 2.0*M_PI*double(t)*double(w);
                zesp += input[w] * exp((-i*odwri)/double(dl));
            }
            y.push_back(zesp);
            modul.push_back(abs(zesp));
        }
        matplot::plot(x,modul);
        matplot::grid(matplot::on);
        matplot::show();
        return y;
    });

    m.def("odwr_transformata",[](vector<complex<double>> input, double poczatek, double koniec, int dl){
        vector<double> modul;
        complex<double> i(0.0,1.0);
        complex<double> zesp;
        vector<complex<double>> y;
        vector<double> x = matplot::linspace(poczatek,koniec,dl);
        double odwri;

        for(int t=0; t<dl; t++){
            zesp=0;
            for(int w=0; w<dl; w++){
                odwri = 2.0*M_PI*double(t)*double(w);
                zesp += input[w] * exp((i*odwri)/double(dl));
            }
            zesp=zesp/double(dl);
            y.push_back(zesp);
            modul.push_back(real(zesp));
        }
        matplot::plot(x,modul);
        matplot::grid(matplot::on);
        matplot::show();
    });




    m.attr("__version__") = "dev";
}
