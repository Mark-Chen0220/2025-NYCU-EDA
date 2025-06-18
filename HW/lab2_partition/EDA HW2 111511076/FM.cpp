#include "fmHeader.h"
#include <iostream>
#include <chrono>

using namespace std;

int main(int argc, char* argv[]) {
    Partitioner partitioner(argv[1]);

    using namespace std::chrono;
    auto start_time = high_resolution_clock::now();
    auto time_limit = seconds(25);
    int num_passes = 0;
    int prev_cut = INT_MAX;
    int patience = 100;
    int wait = 0;

    while (high_resolution_clock::now() - start_time < time_limit && wait < patience) {
        partitioner.run_fm_pass();
        int current_cut = partitioner.get_cut_count();
        if (current_cut < prev_cut) {
            prev_cut = current_cut;
            wait = 0;  // reset wait if we improve
        } else {
            wait++;    // increment wait if no improvement
        }
        ++num_passes;
    }
    // cout << "Total FM passes run: " << num_passes << "\n";
    // cout << "Final cut count: " << partitioner.get_cut_count() << "\n";
    // cout << "Balance Ratio: " <<  double(double(partitioner.get_group0count()) / partitioner.get_nodes().size()) << "\n";
    // auto end_time = high_resolution_clock::now();
    // duration<double> elapsed = end_time - start_time;
    // cout << "Elapsed time: " << elapsed.count() << " seconds\n";
    partitioner.print_status("output.txt");
    return 0;
}
