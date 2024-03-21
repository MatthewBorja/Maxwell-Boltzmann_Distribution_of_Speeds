/*
 * File:   main.cpp
 * Author: Matthew Borja
 * Created on June 1, 2022, 5:27 PM
 * Purpose: Calculate the Maxwell-Boltzmann Distribution of Speeds
 *          and generate a probability versus velocity graph.
 */

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <omp.h>
#include <string>
#include <time.h>
#include <vector>
#include "pbPlots.hpp"
#include "supportLib.hpp"

using namespace std;

void Title();
double CtoK(double);
double mass(double, double);
double mbDos1(double, double, double, double);
double rmsSpeed(double, double, double);
double peakSpeed(double, double, double);
void scatterPlot(vector<double>&, vector<double>&, vector<double>&, vector<double>&, vector<double>&, vector<double>&);
void getYRMS(vector<double>, vector<double>&, vector<double>&);
void getYP(vector<double>&, vector<double>&, vector<double>&);

int main(int argc, char* argv[]) {
    //Define constants.
    double Kb = 1.38E-23;
    double Na = 6.02E+23;
    double R = Na * Kb;
    //Initialize variables.
    int upper_limit = strtol(argv[1], NULL, 10);
    int thread_count = strtol(argv[2], NULL, 10);
    double C, T, M, m, v1, v2 = 0.0;
    //Main Distribution
    vector<double> x, y;
    x.assign(upper_limit, 0.0);
    y.assign(upper_limit, 0.0);
    //RMS
    vector<double> xRMS, yRMS;
    //Peak
    vector<double> xP, yP;
    Title();
    //Get user input for parameters.
    cout << "Enter a molar mass (kg/mol): ";
    cin >> M;
    cout << "Enter a temperature (Celsius): ";
    cin >> C;
    /*cout << "Enter the lower limit of velocities: ";
    cin >> v1;
    cout << "Enter the upper limit of velocities: ";
    cin >> v2;*/
    clock_t begin = clock();
    T = CtoK(C);
    m = mass(M, Na);
    //Calculate RMS and Peak speeds.
    xRMS.push_back(rmsSpeed(Kb, T, m));
    xP.push_back(peakSpeed(R, T, M));
    //Determine the x,y coordinates over the provided.
#pragma omp parallel for num_threads(thread_count)
    for (int i = 0; i < upper_limit; i++) {
        //y.push_back(mbDos1(Kb, T, m, i));
        x[i] = i;
        y[i] = mbDos1(Kb, T, m, i);
    }
    getYRMS(xRMS, yRMS, y);
    getYP(xP, yP, y);
    scatterPlot(x, y, xRMS, yRMS, xP, yP);
    
    clock_t end = clock();
    double time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
    printf("The program took %f seconds to run.\n", time_spent);
    return 0;
}

void Title() {
    cout << "Maxwell-Boltzmann Distribution of Speeds Calculator\n"
        << "This program will calculate the distribution of speeds of particles in an ideal gas at any given temperature\n"
        << "and generate a probability versus velocity graph of the given parameters.\n";

}

//Convert Celsius to Kelvin.
double CtoK(double C) {
    return C + 273.15;
}

//Convert molecular mass to molar mass.
double mass(double M, double Na) {
    return M / Na;
}

//Calculate the Maxwell-Boltzmann Distribution of Speeds.
double mbDos1(double Kb, double T, double m, double v) {
    double result = ((4 / sqrt(3.14)) * pow((m / (2 * Kb * T)), (3 / 2)) * pow(v, 2) * exp((-m * pow(v, 2)) / (2 * Kb * T)));
    //cout << "MBDOS = " << result << endl;
    return result;
}

//Calculate the RMS speed.
double rmsSpeed(double Kb, double T, double m) {
    double result = sqrt((3 * Kb * T) / m);
    //cout << "RMS Speed = " << result << endl;
    return result;
}

//Calculate the peak speed.
double peakSpeed(double R, double T, double M) {
    double result = sqrt((2 * R * T) / M);
    //cout << "Peak speed = " << result << endl;
    return result;
}

//Find the probability for the RMS speed.
void getYRMS(vector<double> vRMS, vector<double>& yRMS, vector<double>& y) {
    int result;
    result = trunc(vRMS[0]);
    try {
        yRMS.push_back(y.at(result));
    }
    catch (out_of_range const& exc) {
        cout << "ERROR: " << exc.what() << ". "
             << "Probability for RMS speed not found. Try a larger velocity value.\n";
    }
    //yRMS.push_back(y[result]);
}

//Find the probability for the peak speed.
void getYP(vector<double>& xP, vector<double>& yP, vector<double>& y) {
    int result;
    result = trunc(xP[0]);
    result -= 0.015;
    try {
        yP.push_back(y.at(result));
    }
    catch (out_of_range const& exc) {
        cout << "ERROR: " << exc.what() << ". "
             << "Probability for peak speed not found. Try a larger velocity value.\n";
    }
    //yP.push_back(y[result]);
}

//Generate the graph.
void scatterPlot(vector<double>& x, vector<double>& y, vector<double>& x2, vector<double>& y2, vector<double>& x3, vector<double>& y3) {
    bool success;
    StringReference* errorMessage = new StringReference();
    RGBABitmapImageReference* imageReference = CreateRGBABitmapImageReference();

    ScatterPlotSeries* series = GetDefaultScatterPlotSeriesSettings();
    series->xs = &x;
    series->ys = &y;
    series->linearInterpolation = true;
    series->lineType = toVector(L"dashed");
    series->lineThickness = 2;
    series->color = CreateRGBColor(0, 0, 0);

    ScatterPlotSeries* series2 = GetDefaultScatterPlotSeriesSettings();
    series2->xs = &x2;
    series2->ys = &y2;
    series2->linearInterpolation = false;
    series2->pointType = toVector(L"dots");
    series2->color = CreateRGBColor(1, 0, 0);

    ScatterPlotSeries* series3 = GetDefaultScatterPlotSeriesSettings();
    series3->xs = &x3;
    series3->ys = &y3;
    series3->linearInterpolation = false;
    series3->pointType = toVector(L"dots");
    series3->color = CreateRGBColor(0, 0, 1);

    ScatterPlotSettings* settings = GetDefaultScatterPlotSettings();
    settings->width = 1200;
    settings->height = 800;
    settings->autoBoundaries = true;
    settings->autoPadding = true;
    settings->title = toVector(L"Maxwell-Boltzmann Distribution of Speeds");
    settings->xLabel = toVector(L"Velocity (m/s)");
    settings->yLabel = toVector(L"Probability");
    settings->scatterPlotSeries->push_back(series);
    settings->scatterPlotSeries->push_back(series2);
    settings->scatterPlotSeries->push_back(series3);

    success = DrawScatterPlotFromSettings(imageReference, settings, errorMessage);
    if (success) {
        vector<double>* pngdata = ConvertToPNG(imageReference->image);
        WriteToFile(pngdata, "Maxwell-Boltzmann_DistributionOfSpeeds.png");
        DeleteImage(imageReference->image);
        cout << "Graph 'Maxwell-Boltzmann_DistributionOfSpeeds.png' generated!\n"
            << "The curve is main distribution, the blue dot is the peak speed, and the red dot is the RMS speed.\n";
    }
    else {
        cerr << "Error: ";
        for (wchar_t c : *errorMessage->string) {
            wcerr << c;
        }
        cerr << endl;
    }
}