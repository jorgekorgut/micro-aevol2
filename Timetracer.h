#pragma once

#ifdef TRACES

#include <chrono>
#include <iostream>
#include <fstream>
#include <vector>

#ifdef USE_OMP

#include <omp.h>

#endif

#define CSV_HEADER "Gen,IndivID,Stamp,TimeStamp_Start,TimeStamp_End,Resource,Duration"
#define GET_TIME std::chrono::steady_clock::now().time_since_epoch().count()

#ifdef USE_OMP
#define GET_RESOURCE omp_get_thread_num()
#else
#define GET_RESOURCE 0
#endif

namespace time_tracer {
    std::vector<const char *> stamp_name;

    static std::ofstream trace_file;
    static long **starts = nullptr;
    static long **ends = nullptr;
    static int **resources = nullptr;

    static int traces_size = 0;

    static void init_tracer(const char *trace_file_name, int size, std::vector<const char *> stamps) {
        trace_file.open(trace_file_name, std::ofstream::trunc);
        trace_file << CSV_HEADER << std::endl;

        stamp_name = std::move(stamps);
        traces_size = size;

        starts = new long *[traces_size];
        ends = new long *[traces_size];
        resources = new int *[traces_size];

        for (int i = 0; i < traces_size; ++i) {
            starts[i] = new long[stamp_name.size()]{};
            ends[i] = new long[stamp_name.size()]{};
            resources[i] = new int[stamp_name.size()]{};
        }
    }

    static void stop_tracer() {
        for (int i = 0; i < traces_size; ++i) {
            delete[] starts[i];
            delete[] ends[i];
            delete[] resources[i];
        }

        delete[] starts;
        delete[] ends;
        delete[] resources;
        starts = nullptr;
        ends = nullptr;
        resources = nullptr;
        traces_size = 0;
        stamp_name.clear();

        trace_file.close();
    }

    static void timestamp_start(int indiv_id, int stamp_id) {
        starts[indiv_id][stamp_id] = GET_TIME;
    }

    static void timestamp_end(int indiv_id, int stamp_id) {
        ends[indiv_id][stamp_id] = GET_TIME;
        resources[indiv_id][stamp_id] = GET_RESOURCE;
    }

    static void set_traces() {
        for (int i = 0; i < traces_size; ++i) {
            for (int j = 0; j < stamp_name.size(); ++j) {
                starts[i][j] = 0;
            }
        }
    }

    static void write_traces(int generation) {
        for (int indiv_id = 0; indiv_id < traces_size; ++indiv_id) {
            for (int stamp = 0; stamp < stamp_name.size(); ++stamp) {
                if (starts[indiv_id][stamp])
                    trace_file << generation << "," << indiv_id
                               << "," << stamp_name[stamp]
                               << "," << starts[indiv_id][stamp]
                               << "," << ends[indiv_id][stamp]
                               << "," << resources[indiv_id][stamp]
                               << "," << ends[indiv_id][stamp] - starts[indiv_id][stamp]
                               << std::endl;
            }
        }
        trace_file.flush();
        set_traces();
    }
}

#define INIT_TRACER(file_name, nb_indivs, ...) time_tracer::init_tracer(file_name, nb_indivs, __VA_ARGS__);

#define TIMESTAMP(id, STAMP, BLOCK) { \
time_tracer::timestamp_start(id, STAMP); \
BLOCK \
time_tracer::timestamp_end(id, STAMP); \
}

#define FLUSH_TRACES(generation) time_tracer::write_traces(generation);
#define STOP_TRACER time_tracer::stop_tracer();

#else //#ifndef TREACES

#define INIT_TRACER(file_name, nb_indivs, ...)
#define TIMESTAMP(id, STAMP, BLOCK) BLOCK
#define FLUSH_TRACES(generation)
#define STOP_TRACER

#endif //TRACES