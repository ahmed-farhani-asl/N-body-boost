#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <boost/numeric/odeint.hpp>

using namespace std;
using namespace boost::numeric::odeint;

// Gravitational Constants
const double G = 6.67408e-11;
const double au = 1.495978707e11;
const double sm = 1.98847e30;
const double yr = 3.1471401e7;

double g = G * sm * pow(yr, 2) / pow(au, 3); // G in simulation units 
double dt = 0.005;                            // Time step
double end_time = 1000;
size_t steps = end_time / dt;

// Number of bodies (determined dynamically)
int N = 0;

// State vectors
vector<double> mass;
vector<double> x;

// Function to count lines in input file (to determine N)
int count_lines(const string &filename) {
    ifstream file(filename);

    int j = 0;
    string str;
    while(getline(file, str)) j++;

    return j;
    //return count(istreambuf_iterator<char>(file), istreambuf_iterator<char>(), '\n');
}

// Acceleration Calculation
void accelerate(const vector<double> &x, int i, double &ax, double &ay, double &az) {
    ax = ay = az = 0.0;

    for (int j = 0; j < N; j++) {
        if (j == i) continue;

        double dx = x[6 * j] - x[6 * i];
        double dy = x[6 * j + 1] - x[6 * i + 1];
        double dz = x[6 * j + 2] - x[6 * i + 2];

        double r2 = dx * dx + dy * dy + dz * dz;
        double r3_inv = 1.0 / (sqrt(r2) * r2);

        ax += mass[j] * dx * r3_inv;
        ay += mass[j] * dy * r3_inv;
        az += mass[j] * dz * r3_inv;
    }

    ax *= g;
    ay *= g;
    az *= g;
}

// ODEs for velocity & acceleration update
void ODEs(const vector<double> &x, vector<double> &dxdt, double /* t */) {
    
    for (int i = 0; i < N; i++) {
        dxdt[6 * i] = x[6 * i + 3];
        dxdt[6 * i + 1] = x[6 * i + 4];
        dxdt[6 * i + 2] = x[6 * i + 5];

        double ax, ay, az;
        accelerate(x, i, ax, ay, az);

        dxdt[6 * i + 3] = ax;
        dxdt[6 * i + 4] = ay;
        dxdt[6 * i + 5] = az;
    }
}

// Calculate Kinetic, Potential, and Total Energy
void compute_energy(const vector<double> &x, double &KE, double &PE) {
    KE = 0.0;
    PE = 0.0;

    // Kinetic Energy: Sum (1/2 * m * v^2)
    for (int i = 0; i < N; i++) {
        double v2 = x[6 * i + 3] * x[6 * i + 3] +
                    x[6 * i + 4] * x[6 * i + 4] +
                    x[6 * i + 5] * x[6 * i + 5];

        KE += 0.5 * mass[i] * v2;
    }

    // Potential Energy: Sum (-G * m_i * m_j / r_ij)
    for (int i = 0; i < N; i++) {
        for (int j = i + 1; j < N; j++) {
            double dx = x[6 * j] - x[6 * i];
            double dy = x[6 * j + 1] - x[6 * i + 1];
            double dz = x[6 * j + 2] - x[6 * i + 2];

            double r = sqrt(dx * dx + dy * dy + dz * dz);
            PE += -g * mass[i] * mass[j] / r;
        }
    }
}

// Read input data dynamically
void input() {
    ifstream read("input.txt");
    if (!read) {
        cerr << "Error: Could not open input.txt\n";
        exit(1);
    }

    N = count_lines("input.txt");  // Count bodies dynamically
    mass.resize(N);
    x.resize(6 * N);

    cout << "G = " << g << "\t dt = " << dt << "\t Simulation time = " << end_time << endl;
    cout << "Number of bodies: " << N << endl;

    for (int i = 0; i < N; i++) {
        read >> mass[i];
        for (int j = 0; j < 6; j++) read >> x[6 * i + j];
    }

    read.close();
}

// Output to CSV file
void output(int step, ofstream &csvFile) {
    double KE, PE;
    compute_energy(x, KE, PE);
    double E = KE + PE;

    csvFile << step * dt << "," << KE << "," << PE << "," << E; // Time & Energies
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < 6; j++) {
            csvFile << "," << x[6 * i + j];
        }
    }
    csvFile << "\n";
}

int main() {

    input();

    //runge_kutta4<vector<double>> stepper;
    runge_kutta_fehlberg78<vector<double>> stepper;
    stepper.adjust_size(x);

    ofstream csvFile("simulation_output.csv");
    if (!csvFile) {
        cerr << "Error: Could not open simulation_output.csv\n";
        exit(1);
    }

    // Write CSV header
    csvFile << "time,KE,PE,Total_Energy";
    for (int i = 0; i < N; i++) {
        csvFile << ",x" << i << ",y" << i << ",z" << i << ",vx" << i << ",vy" << i << ",vz" << i;
    }
    csvFile << "\n";

    double t = 0.0;
    for (size_t i = 0; i <= steps; i++, t += dt) {
        //if(i%5==0) 
        output(i, csvFile);
        stepper.do_step(ODEs, x, t, dt);
    }

    csvFile.close();
    cout << "Simulation completed. Output saved in simulation_output.csv" << endl;

    return 0;
}
