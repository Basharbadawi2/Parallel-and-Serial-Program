#include <iostream>
#include <vector>
#include <cmath>
#include <thread>
#include <mutex>

using namespace std;

class ChemicalSystem {
private:
    double H2O2, H2O, O2, K, Time;
    double dH2O2_dt, dH2O_dt, dO2_dt;
public:
    ChemicalSystem(double K, double H2O2, double H2O, double O2, double Time, double dH2O2_dt, double dO2_dt) {
        this->K = K;
        this->H2O2 = H2O2;
        this->H2O = H2O;
        this->O2 = O2;
        this->Time = Time;
        this->dH2O2_dt = dH2O2_dt;
        this->dO2_dt = dO2_dt;
    }

    vector<double> Rate_Function(const vector<double>& Concentrations) {
        H2O2 = Concentrations[0];
        H2O = Concentrations[1];
        O2 = Concentrations[2];

        dH2O2_dt = -K * pow(H2O2, 2);
        dH2O_dt = 2 * K * pow(H2O2, 2);
        dO2_dt = K * pow(H2O2, 2);

        vector<double> Rates = {dH2O2_dt, dH2O_dt, dO2_dt};
        return Rates;
    }

    void ParallelProgram(double Time) {
        this->Time = Time;
        vector<thread> Thread1;

        // Define the time step
        double dt = 0.01;

        // Define the initial concentrations
        vector<double> concentrations = {H2O2, H2O, O2};

        mutex mtx;
        auto simulate = [&](){
            for (double t = 0.0; t < Time; t += dt) {

                vector<double> rates = Rate_Function(concentrations);

                mtx.lock();
                concentrations[0] += rates[0] * dt;
                concentrations[1] += rates[1] * dt;
                concentrations[2] += rates[2] * dt;
                mtx.unlock();
            }
        };

        for (auto& t : Thread1) {
            t.join();
        }

        cout << "Parallel simulation: [H2O2] = " << concentrations[0] << ", [H2O] = " << concentrations[1] << ", [O2] = " << concentrations[2] << endl;
    }

    void Run_Simulation(double Start, double End, vector<double> concentrations) {

        double dt = 0.01;


        double t = Start;


        while (t < End) {

            vector<double> rates = Rate_Function(concentrations);
            concentrations[0] += rates[0] * dt;
            concentrations[1] += rates[1] * dt;
            concentrations[2] += rates[2] * dt;
            t += dt;
        }

        cout << "Serial simulation: [H2O2] = " << concentrations[0] << ", [H2O] = " << concentrations[1] << ", [O2] = " << concentrations[2] << endl;
    }
};

int main() {
    double K = 0.1;
    double H2O2 = 1.0;
    double H2O = 0.0;
    double O2 = 0.0;
    double Time = 10.0;
    double dH2O2_dt = 0.0;
    double dO2_dt = 0.0;

    ChemicalSystem system(K, H2O2, H2O, O2, Time, dH2O2_dt, dO2_dt);

    double Start = 0.0;
    double End = 10.0;
    vector<double> concentrations = {H2O2, H2O, O2};

    system.Run_Simulation(Start, End, concentrations);
    system.ParallelProgram(Time);

    return 0;
}
